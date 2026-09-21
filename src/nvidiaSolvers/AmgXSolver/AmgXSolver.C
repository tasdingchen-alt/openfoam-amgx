/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | openfoam-amgx
   \\    /   O peration     |
    \\  /    A nd           |
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of openfoam-amgx.

    openfoam-amgx is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by the
    Free Software Foundation, either version 3 of the License, or (at your
    option) any later version.

    openfoam-amgx is distributed in the hope that it will be useful, but
    WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY
    or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License along
    with openfoam-amgx.  If not, see <https://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "AmgXSolver.H"
#include "processorLduInterface.H"
#include "Pstream.H"
#include "PstreamGlobals.H"

#include <cuda_runtime.h>

#include <sstream>
#include <vector>
#include <map>
#include <memory>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define AMGX_CHECK_CALL(call)                                                 \
    do                                                                        \
    {                                                                         \
        AMGX_RC amgxRc = (call);                                              \
        if (amgxRc != AMGX_RC_OK)                                             \
        {                                                                     \
            FatalErrorInFunction                                              \
                << "AmgX call failed: " #call                                 \
                << " (rc=" << int(amgxRc) << ")"                              \
                << exit(FatalError);                                          \
        }                                                                     \
    }                                                                         \
    while (false)

namespace Foam
{
    defineTypeNameAndDebug(AmgXSolver, 0);

    lduMatrix::solver::addsymMatrixConstructorToTable<AmgXSolver>
        addAmgXSolverSymMatrixConstructorToTable_;

    lduMatrix::solver::addasymMatrixConstructorToTable<AmgXSolver>
        addAmgXSolverAsymMatrixConstructorToTable_;
}


// * * * * * * * * * * * * * * * * *  State  * * * * * * * * * * * * * * * * //

Foam::AmgXSolver::State::State()
:
    cfg(nullptr),
    rsrc(nullptr),
    mat(nullptr),
    x(nullptr),
    b(nullptr),
    solver(nullptr),
    csr(),
    nCells(0),
    nGlobal(0),
    parallel(false),
    initialised(false),
    patternSet(false),
    needSetup(true)
{}


Foam::AmgXSolver::State::~State()
{
    // Intentionally do not destroy the AmgX handles here.
    //
    // The cache is a function-local static and is therefore destroyed at
    // process exit, by which point AmgX's own internal state (and possibly the
    // CUDA context) may already be gone. Calling AMGX_*_destroy then raises
    // "Invalid C wrapper" and corrupts the heap. The handles are released by
    // the operating system at exit; mid-run invalidation (mesh change) goes
    // through AmgXSolver::shutdown() which does destroy them.
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::AmgXSolver::AmgXSolver
(
    const word& fieldName,
    const lduMatrix& matrix,
    const FieldField<Field, scalar>& interfaceBouCoeffs,
    const FieldField<Field, scalar>& interfaceIntCoeffs,
    const lduInterfaceFieldPtrsList& interfaces,
    const dictionary& solverControls
)
:
    lduMatrix::solver
    (
        fieldName,
        matrix,
        interfaceBouCoeffs,
        interfaceIntCoeffs,
        interfaces,
        solverControls
    ),
    configFile_(),
    mode_("dDDI"),
    verbose_(false),
    setupEveryTime_(true),
    solverName_("FGMRES"),
    algorithm_("AGGREGATION"),
    smoother_("BLOCK_JACOBI"),
    coarseSolver_("NOSOLVER"),
    interpolator_("D2"),
    selector_("SIZE_2"),
    maxLevels_(50)
{
    // Solver options live in the "AmgX" sub-dictionary
    const dictionary& amgxDict = controlDict_.subOrEmptyDict("AmgX");

    configFile_ = amgxDict.lookupOrDefault<word>("configFile", word(""));
    mode_ = amgxDict.lookupOrDefault<word>("mode", "dDDI");
    verbose_ = amgxDict.lookupOrDefault<bool>("verbose", false);
    setupEveryTime_ = amgxDict.lookupOrDefault<bool>("setupEveryTime", true);

    solverName_ = amgxDict.lookupOrDefault<word>("solver", "FGMRES");
    algorithm_ = amgxDict.lookupOrDefault<word>("algorithm", "AGGREGATION");
    smoother_ = amgxDict.lookupOrDefault<word>("smoother", "BLOCK_JACOBI");
    coarseSolver_ = amgxDict.lookupOrDefault<word>("coarseSolver", "NOSOLVER");
    interpolator_ = amgxDict.lookupOrDefault<word>("interpolator", "D2");
    selector_ = amgxDict.lookupOrDefault<word>("selector", "SIZE_2");
    maxLevels_ = amgxDict.lookupOrDefault<label>("maxLevels", 50);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

Foam::AmgXSolver::~AmgXSolver()
{
    // The AmgX objects are cached process-wide and are not destroyed here.
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

Foam::AmgXSolver::State& Foam::AmgXSolver::cachedState(const std::string& key)
{
    static std::map<std::string, std::unique_ptr<State>> states;

    std::unique_ptr<State>& p = states[key];

    if (!p)
    {
        p.reset(new State());
    }

    return *p;
}


std::string Foam::AmgXSolver::defaultConfig() const
{
    std::ostringstream os;

    // AmgX convergence is relative to the initial residual. OpenFOAM's
    // normalised initial residual is O(1), so its absolute tolerance is a
    // good approximation of an AmgX relative tolerance when relTol is 0.
    const scalar amgxTol = (relTol_ > 0.0 ? relTol_ : tolerance_);

    os  << "{"
        << "\"config_version\": 2, "
        << "\"solver\": {"
        <<   "\"preconditioner\": {"
        <<     "\"error_scaling\": 0, "
        <<     "\"print_grid_stats\": 1, "
        <<     "\"algorithm\": \"" << algorithm_ << "\", "
        <<     "\"solver\": \"AMG\", "
        <<     "\"smoother\": \"" << smoother_ << "\", "
        <<     "\"presweeps\": 0, "
        <<     "\"selector\": \"" << selector_ << "\", "
        <<     "\"coarse_solver\": \"" << coarseSolver_ << "\", "
        <<     "\"max_iters\": 1, "
        <<     "\"min_coarse_rows\": 32, "
        <<     "\"relaxation_factor\": 0.75, "
        <<     "\"scope\": \"amg\", "
        <<     "\"max_levels\": " << maxLevels_ << ", "
        <<     "\"postsweeps\": 3, "
        <<     "\"cycle\": \"V\"";

    if (algorithm_ == "CLASSICAL")
    {
        // D1 is not supported in distributed settings
        os << ", \"interpolator\": \"" << interpolator_ << "\"";
    }

    os  <<   "}, "
        <<   "\"use_scalar_norm\": 1, "
        <<   "\"solver\": \"" << solverName_ << "\", "
        <<   "\"print_solve_stats\": " << (verbose_ ? 1 : 0) << ", "
        <<   "\"obtain_timings\": " << (verbose_ ? 1 : 0) << ", "
        <<   "\"max_iters\": " << maxIter_ << ", "
        <<   "\"monitor_residual\": 1, "
        <<   "\"gmres_n_restart\": 32, "
        <<   "\"convergence\": \"RELATIVE_INI\", "
        <<   "\"scope\": \"main\", "
        <<   "\"tolerance\": " << amgxTol << ", "
        <<   "\"norm\": \"L2\""
        << "}"
        << "}";

    return os.str();
}


AMGX_Mode Foam::AmgXSolver::readMode() const
{
    if (mode_ == "dDDI") return AMGX_mode_dDDI;
    if (mode_ == "dDFI") return AMGX_mode_dDFI;
    if (mode_ == "dFFI") return AMGX_mode_dFFI;
    if (mode_ == "hDDI") return AMGX_mode_hDDI;
    if (mode_ == "hDFI") return AMGX_mode_hDFI;
    if (mode_ == "hFFI") return AMGX_mode_hFFI;

    FatalErrorInFunction
        << "Unsupported AmgX mode " << mode_ << nl
        << "Supported: dDDI, dDFI, dFFI, hDDI, hDFI, hFFI"
        << exit(FatalError);

    return AMGX_mode_dDDI;
}


void Foam::AmgXSolver::initialise(State& st) const
{
    // AmgX requires a one-time process initialisation
    static bool amgxInitialised = false;

    if (!amgxInitialised)
    {
        AMGX_CHECK_CALL(AMGX_initialize());
        amgxInitialised = true;
    }

    if (!configFile_.empty())
    {
        AMGX_CHECK_CALL
        (
            AMGX_config_create_from_file(&st.cfg, configFile_.c_str())
        );
    }
    else
    {
        const std::string cfgStr = defaultConfig();

        if (verbose_)
        {
            Info<< "AmgX config: " << cfgStr << endl;
        }

        AMGX_CHECK_CALL(AMGX_config_create(&st.cfg, cfgStr.c_str()));
    }

    st.parallel = Pstream::parRun();

    if (st.parallel)
    {
        // Use the OpenFOAM world communicator (a persistent global) and
        // assign one device per rank
        int device = 0;
        int numDevices = 0;

        if (cudaGetDeviceCount(&numDevices) == cudaSuccess && numDevices > 0)
        {
            device = Pstream::myProcNo() % numDevices;
        }

        AMGX_CHECK_CALL
        (
            AMGX_resources_create
            (
                &st.rsrc,
                st.cfg,
                &PstreamGlobals::MPI_COMM_FOAM,
                1,
                &device
            )
        );
    }
    else
    {
        AMGX_CHECK_CALL(AMGX_resources_create_simple(&st.rsrc, st.cfg));
    }

    const AMGX_Mode mode = readMode();

    AMGX_CHECK_CALL(AMGX_matrix_create(&st.mat, st.rsrc, mode));
    AMGX_CHECK_CALL(AMGX_vector_create(&st.x, st.rsrc, mode));
    AMGX_CHECK_CALL(AMGX_vector_create(&st.b, st.rsrc, mode));
    AMGX_CHECK_CALL(AMGX_solver_create(&st.solver, st.rsrc, mode, st.cfg));

    st.initialised = true;
}


void Foam::AmgXSolver::shutdown(State& st) const
{
    if (st.initialised)
    {
        AMGX_solver_destroy(st.solver);
        AMGX_matrix_destroy(st.mat);
        AMGX_vector_destroy(st.x);
        AMGX_vector_destroy(st.b);
        AMGX_resources_destroy(st.rsrc);
        AMGX_config_destroy(st.cfg);

        st.cfg = nullptr;
        st.rsrc = nullptr;
        st.mat = nullptr;
        st.x = nullptr;
        st.b = nullptr;
        st.solver = nullptr;
        st.initialised = false;
    }
}


void Foam::AmgXSolver::buildHalos
(
    const label rowStart,
    List<ldu2csr::Entry>& halos
) const
{
    forAll(interfaces_, i)
    {
        if (!interfaces_.set(i))
        {
            continue;
        }

        const lduInterface& intf = interfaces_[i].interface();

        // Only processor interfaces are supported for now
        const processorLduInterface* pi =
            dynamic_cast<const processorLduInterface*>(&intf);

        if (!pi)
        {
            continue;
        }

        const labelUList& fc = intf.faceCells();
        const scalarField& bou = interfaceBouCoeffs_[i];

        // Global row index of the local cells on this interface
        labelField myG(fc.size());
        forAll(fc, j)
        {
            myG[j] = rowStart + fc[j];
        }

        // The matching global row indices on the neighbour process
        labelField nbrG(fc.size(), 0);

        if (pi->myProcNo() < pi->neighbProcNo())
        {
            pi->send(Pstream::commsTypes::blocking, myG);
            pi->receive(Pstream::commsTypes::blocking, nbrG);
        }
        else
        {
            pi->receive(Pstream::commsTypes::blocking, nbrG);
            pi->send(Pstream::commsTypes::blocking, myG);
        }

        // Halo entry: local row, global column, coupling coefficient
        forAll(fc, j)
        {
            ldu2csr::Entry e;
            e.row = fc[j];
            e.col = nbrG[j];
            e.value = -bou[j];
            halos.append(e);
        }
    }

    if (::getenv("AMGX_DEBUG_HALO"))
    {
        label minC = labelMax;
        label maxC = -labelMax;
        forAll(halos, k)
        {
            minC = min(minC, halos[k].col);
            maxC = max(maxC, halos[k].col);
        }
        Pout<< "rank " << Pstream::myProcNo()
            << " rowStart=" << rowStart
            << " nHalos=" << halos.size()
            << " colMin=" << minC
            << " colMax=" << maxC << endl;
    }
}


void Foam::AmgXSolver::globalRange
(
    const label nCells,
    label& rowStart,
    label& nGlobal
) const
{
    // Cells are numbered rank-contiguously: this rank owns the rows
    // [rowStart, rowStart + nCells) of the global matrix
    rowStart = 0;
    nGlobal = nCells;

    if (!Pstream::parRun())
    {
        return;
    }

    labelList nCellsPerProc(Pstream::nProcs(), 0);
    nCellsPerProc[Pstream::myProcNo()] = nCells;
    Pstream::gatherList(nCellsPerProc);
    Pstream::scatterList(nCellsPerProc);

    nGlobal = 0;
    forAll(nCellsPerProc, proc)
    {
        if (proc < Pstream::myProcNo())
        {
            rowStart += nCellsPerProc[proc];
        }
        nGlobal += nCellsPerProc[proc];
    }

    if (::getenv("AMGX_DEBUG_PAR"))
    {
        Pout<< "rank " << Pstream::myProcNo() << "/" << Pstream::nProcs()
            << " nCells=" << nCells << " nGlobal=" << nGlobal
            << " rowStart=" << rowStart
            << " sumLocal=" << sum(nCellsPerProc) << endl;
    }
}


std::string Foam::AmgXSolver::stateKey() const
{
    std::ostringstream key;

    key << fieldName_ << ':' << static_cast<const void*>(&matrix_.mesh())
        << ':' << relTol_ << ':' << tolerance_ << ':' << maxIter_
        << ':' << int(Pstream::parRun());

    return key.str();
}


void Foam::AmgXSolver::assembleCsr
(
    State& st,
    const label rowStart,
    const label nCells,
    const label nFaces
) const
{
    if (Pstream::parRun())
    {
        // IMPORTANT: the halo coefficients change between solves, so the
        // halo entries are rebuilt on every call. Caching them (e.g. for
        // performance) makes AmgX solve a stale matrix and stalls FGMRES.
        List<ldu2csr::Entry> halos;
        buildHalos(rowStart, halos);

        st.csr.reset(new ldu2csr(matrix_, rowStart, halos));
    }
    else if
    (
        !st.csr.valid()
     || st.csr->nRows() != nCells
     || st.csr->nNonZeros() != nCells + 2*nFaces
    )
    {
        st.csr.reset(new ldu2csr(matrix_));
        st.patternSet = false;
        st.needSetup = true;
    }
    else
    {
        st.csr->updateValues(matrix_);
    }
}


bool Foam::AmgXSolver::normaliseSign(State& st) const
{
    const scalar sign = (matrix_.diag()[0] < 0.0) ? -1.0 : 1.0;

    if (sign > 0.0)
    {
        return false;
    }

    scalar* vals = st.csr->values();

    for (label i = 0; i < st.csr->nNonZeros(); i++)
    {
        vals[i] = -vals[i];
    }

    return true;
}


Foam::scalar Foam::AmgXSolver::normalisedResidual
(
    const scalarField& psi,
    const scalarField& source,
    const direction cmpt
) const
{
    scalarField Apsi(psi.size());
    matrix_.Amul(Apsi, psi, interfaceBouCoeffs_, interfaces_, cmpt);

    const scalarField rA(source - Apsi);
    scalarField tmp(psi.size());

    const scalar nf = this->normFactor(psi, source, Apsi, tmp);

    return gSumMag(rA, matrix_.mesh().comm())/nf;
}


void Foam::AmgXSolver::uploadMatrix(State& st) const
{
    const label nRows = st.csr->nRows();
    const label nnz = st.csr->nNonZeros();

    if (!st.patternSet)
    {
        // AmgX dDDI uses int32 indices
        std::vector<int> rowPtr(nRows + 1);
        std::vector<int> colInd(nnz);

        for (label i = 0; i <= nRows; i++)
        {
            rowPtr[i] = int(st.csr->rowPtr()[i]);
        }

        for (label i = 0; i < nnz; i++)
        {
            colInd[i] = int(st.csr->colInd()[i]);
        }

        if (st.parallel)
        {
            // Distributed upload with rank-contiguous global numbering:
            // off-processor columns are halo entries with global indices
            int nrings = 1;
            AMGX_CHECK_CALL
            (
                AMGX_config_get_default_number_of_rings(st.cfg, &nrings)
            );

            AMGX_CHECK_CALL
            (
                AMGX_matrix_upload_all_global_32
                (
                    st.mat,
                    int(st.nGlobal),
                    int(st.nCells),
                    int(nnz),
                    1,
                    1,
                    rowPtr.data(),
                    colInd.data(),
                    st.csr->values(),
                    nullptr,
                    nrings,
                    nrings,
                    nullptr
                )
            );

            AMGX_CHECK_CALL(AMGX_vector_bind(st.x, st.mat));
            AMGX_CHECK_CALL(AMGX_vector_bind(st.b, st.mat));
        }
        else
        {
            AMGX_CHECK_CALL
            (
                AMGX_matrix_upload_all
                (
                    st.mat,
                    int(st.nCells),
                    int(nnz),
                    1,
                    1,
                    rowPtr.data(),
                    colInd.data(),
                    st.csr->values(),
                    nullptr
                )
            );
        }

        st.patternSet = true;
        st.needSetup = true;
    }
    else
    {
        AMGX_CHECK_CALL
        (
            AMGX_matrix_replace_coefficients
            (
                st.mat,
                int(st.nCells),
                int(nnz),
                st.csr->values(),
                nullptr
            )
        );
    }
}


Foam::label Foam::AmgXSolver::solveAmgX
(
    State& st,
    scalarField& psi,
    const scalarField& rhs
) const
{
    // x holds psi as the initial guess
    AMGX_CHECK_CALL(AMGX_vector_upload(st.b, st.nCells, 1, rhs.begin()));
    AMGX_CHECK_CALL(AMGX_vector_upload(st.x, st.nCells, 1, psi.begin()));

    AMGX_CHECK_CALL(AMGX_solver_solve(st.solver, st.b, st.x));
    AMGX_CHECK_CALL(AMGX_vector_download(st.x, psi.begin()));

    label iters = 0;
    AMGX_CHECK_CALL(AMGX_solver_get_iterations_number(st.solver, &iters));

    if (verbose_)
    {
        AMGX_SOLVE_STATUS status;
        AMGX_CHECK_CALL(AMGX_solver_get_status(st.solver, &status));

        Info<< "    AmgX status=" << int(status)
            << " iterations=" << iters << endl;
    }

    return iters;
}


Foam::solverPerformance Foam::AmgXSolver::solve
(
    scalarField& psi,
    const scalarField& source,
    const direction cmpt
) const
{
    solverPerformance solverPerf(typeName, fieldName_);

    const label nCells = psi.size();
    const label nFaces = matrix_.upper().size();

    // --- Rank-contiguous global cell numbering (parallel) ---
    label rowStart = 0;
    label nGlobal = nCells;
    globalRange(nCells, rowStart, nGlobal);

    // --- Process-wide cached AmgX state ---
    State& st = cachedState(stateKey());

    // Recreate the AmgX objects if the matrix structure changed
    if
    (
        st.initialised
     && (st.nCells != nCells || st.nGlobal != nGlobal
         || st.parallel != Pstream::parRun())
    )
    {
        shutdown(st);
        st.csr.clear();
        st.patternSet = false;
        st.needSetup = true;
    }

    if (!st.initialised)
    {
        initialise(st);
        st.nCells = nCells;
        st.nGlobal = nGlobal;
        st.parallel = Pstream::parRun();
    }

    // --- CSR assembly ---
    assembleCsr(st, rowStart, nCells, nFaces);

    // --- Sign normalisation: solve (-A) x = (-b) for a negative diagonal ---
    const bool signFlipped = normaliseSign(st);

    // --- Initial residual, OpenFOAM-normalised ---
    solverPerf.initialResidual() = normalisedResidual(psi, source, cmpt);
    solverPerf.finalResidual() = solverPerf.initialResidual();

    // --- Matrix upload / coefficient update ---
    uploadMatrix(st);

    // --- (Re)build the AMG hierarchy when required ---
    if (st.needSetup || setupEveryTime_)
    {
        AMGX_CHECK_CALL(AMGX_solver_setup(st.solver, st.mat));
        st.needSetup = false;
    }

    // --- RHS, negated to match the sign normalisation ---
    scalarField rhs(source);

    if (signFlipped)
    {
        rhs = -rhs;
    }

    // --- Upload vectors, solve, download ---
    solverPerf.nIterations() = solveAmgX(st, psi, rhs);

    AMGX_SOLVE_STATUS status;
    AMGX_CHECK_CALL(AMGX_solver_get_status(st.solver, &status));

    // Rebuild the hierarchy next time if this solve did not converge
    if (status != AMGX_SOLVE_SUCCESS)
    {
        st.needSetup = true;
    }

    // --- Final residual, OpenFOAM-normalised ---
    solverPerf.finalResidual() = normalisedResidual(psi, source, cmpt);

    solverPerf.checkConvergence(tolerance_, relTol_);

    return solverPerf;
}


// ************************************************************************* //
