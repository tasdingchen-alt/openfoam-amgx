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

Class
    Foam::AmgXSolver

Description
    GPU linear solver for OpenFOAM lduMatrices backed by NVIDIA AmgX.

    OpenFOAM creates a new solver object for every linear solve, so the AmgX
    objects (config, resources, matrix, vectors, solver) are kept in a
    process-wide cache keyed by field name and mesh. On subsequent solves only
    the matrix coefficients and the vectors are updated, which avoids the
    expensive AMG setup and repeated allocation.

    The lduMatrix is assembled into CSR and handed to AmgX. A JSON (or legacy)
    AmgX configuration can be supplied with the "configFile" entry; otherwise a
    sensible default is generated from the OpenFOAM controls.

    Example dictionary entry:

        solver      AmgX;
        AmgX
        {
            configFile      "amgx_config.json";   // optional
            mode            dDDI;                 // optional (default dDDI)
            verbose         false;                // optional
            setupEveryTime  false;                // optional
        }
        tolerance   1e-8;
        relTol      0.01;

SourceFiles
    AmgXSolver.C

\*---------------------------------------------------------------------------*/

#ifndef AmgXSolver_H
#define AmgXSolver_H

#include "lduMatrix.H"
#include "ldu2csr.H"

#include <amgx_c.h>

#include <string>
#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                         Class AmgXSolver Declaration
\*---------------------------------------------------------------------------*/

class AmgXSolver
:
    public lduMatrix::solver
{
public:

    //- Cached per-field AmgX state
    struct State
    {
        //- AmgX handles
        AMGX_config_handle cfg;
        AMGX_resources_handle rsrc;
        AMGX_matrix_handle mat;
        AMGX_vector_handle x;
        AMGX_vector_handle b;
        AMGX_solver_handle solver;

        //- CSR assembly (cached between solves)
        autoPtr<ldu2csr> csr;

        //- Structure of the cached matrix
        label nCells;
        label nGlobal;

        //- Parallel (MPI) build?
        bool parallel;

        //- Has AmgX been initialised?
        bool initialised;

        //- Has the matrix pattern been uploaded?
        bool patternSet;

        //- Does the AMG hierarchy need rebuilding before the next solve?
        bool needSetup;

        State();

        ~State();

        State(const State&) = delete;

        void operator=(const State&) = delete;
    };


private:

    // Private data

        //- Optional AmgX configuration file
        word configFile_;

        //- AmgX mode string (e.g. dDDI)
        word mode_;

        //- Verbose output
        bool verbose_;

        //- Rebuild the AMG hierarchy on every solve
        bool setupEveryTime_;

        //- Krylov solver (FGMRES, GMRES, PCG, ...)
        word solverName_;

        //- AMG algorithm (AGGREGATION, CLASSICAL)
        word algorithm_;

        //- Smoother (BLOCK_JACOBI, MULTICOLOR_DILU, ...)
        word smoother_;

        //- Coarse solver (NOSOLVER, DENSE_LU_SOLVER, ...)
        word coarseSolver_;

        //- Interpolator (CLASSICAL only; D2 is required in parallel)
        word interpolator_;

        //- Aggregation selector
        word selector_;

        //- Maximum number of AMG levels
        label maxLevels_;

        //- FGMRES restart length (Krylov workspace size)
        label gmresRestart_;


    // Private Member Functions

        //- Generate a default AmgX configuration from the OpenFOAM controls
        std::string defaultConfig() const;

        //- Map the mode string to the AmgX mode enum
        AMGX_Mode readMode() const;

        //- Create the AmgX objects for a state
        void initialise(State& st) const;

        //- Destroy the AmgX objects for a state
        void shutdown(State& st) const;

        //- Build the off-processor (halo) matrix entries
        void buildHalos
        (
            const label rowStart,
            List<ldu2csr::Entry>& halos
        ) const;

        //- Return the process-wide cached state for a key
        static State& cachedState(const std::string& key);

        //- Compute the rank-contiguous global row offset (rowStart) and the
        //  global row count (nGlobal) for the parallel numbering scheme
        void globalRange
        (
            const label nCells,
            label& rowStart,
            label& nGlobal
        ) const;

        //- Cache key: the AmgX configuration (tolerance, max_iters) is baked
        //  into the solver at creation, so the controls are part of the key
        std::string stateKey() const;

        //- Assemble the CSR matrix. In parallel the halo entries are rebuilt
        //  on every call so their coefficients stay up to date
        void assembleCsr
        (
            State& st,
            const label rowStart,
            const label nCells,
            const label nFaces
        ) const;

        //- Solve (-A, -b) instead of (A, b) when the diagonal is negative
        //  (OpenFOAM's pressure matrix is negative definite). AmgX's
        //  aggregation AMG assumes a positive diagonal. Returns true when
        //  the sign was flipped
        bool normaliseSign(State& st) const;

        //- OpenFOAM-normalised residual of (A psi - source)
        scalar normalisedResidual
        (
            const scalarField& psi,
            const scalarField& source,
            const direction cmpt
        ) const;

        //- Upload the matrix pattern or replace the coefficients
        void uploadMatrix(State& st) const;

        //- Upload the RHS and the initial guess, solve, download the result
        //  and return the iteration count
        label solveAmgX
        (
            State& st,
            scalarField& psi,
            const scalarField& rhs
        ) const;


public:

    //- Runtime type information
    TypeName("AmgX");


    // Constructors

        //- Construct from matrix components and solver controls
        AmgXSolver
        (
            const word& fieldName,
            const lduMatrix& matrix,
            const FieldField<Field, scalar>& interfaceBouCoeffs,
            const FieldField<Field, scalar>& interfaceIntCoeffs,
            const lduInterfaceFieldPtrsList& interfaces,
            const dictionary& solverControls
        );

        //- Disallow default bitwise copy construction
        AmgXSolver(const AmgXSolver&) = delete;


    //- Destructor
    virtual ~AmgXSolver();


    // Member Functions

        //- Solve the matrix with this solver
        virtual solverPerformance solve
        (
            scalarField& psi,
            const scalarField& source,
            const direction cmpt = 0
        ) const;


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const AmgXSolver&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
