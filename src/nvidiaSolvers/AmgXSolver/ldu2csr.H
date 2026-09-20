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
    Foam::ldu2csr

Description
    Assemble an OpenFOAM lduMatrix into CSR format.

    OpenFOAM stores the matrix as (diag, upper, lower) over the mesh faces.
    For an internal face f with owner = upperAddr()[f] and
    neighbour = lowerAddr()[f]:

        lower[f] = A(owner, neighbour)   -> row owner,     column neighbour
        upper[f] = A(neighbour, owner)   -> row neighbour, column owner
        diag[i]  = A(i, i)

    The diagonal already contains the implicit contribution of non-coupled
    boundary conditions (see fvMatrix::solveSegregated and addBoundaryDiag).

    In parallel, columns are global cell indices (local row offset = rowStart)
    and off-processor couplings are supplied as additional halo entries.

SourceFiles
    ldu2csr.C

\*---------------------------------------------------------------------------*/

#ifndef ldu2csr_H
#define ldu2csr_H

#include "lduMatrix.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

/*---------------------------------------------------------------------------*\
                           Class ldu2csr Declaration
\*---------------------------------------------------------------------------*/

class ldu2csr
{
public:

    //- Additional (off-processor) matrix entry
    struct Entry
    {
        label row;
        label col;
        scalar value;
    };


private:

    // Private data

        //- Number of rows (cells)
        label nRows_;

        //- Number of non-zeros
        label nNonZeros_;

        //- CSR row offsets (size nRows_ + 1)
        labelList rowPtr_;

        //- CSR column indices (size nNonZeros_), global in parallel
        labelList colInd_;

        //- CSR values (size nNonZeros_)
        scalarField values_;


    // Private Member Functions

        //- Assemble the CSR structure and values from the matrix
        void assemble
        (
            const lduMatrix& matrix,
            const label rowStart,
            const List<Entry>& extra
        );


public:

    // Constructors

        //- Construct from an lduMatrix (serial, local column indices)
        explicit ldu2csr(const lduMatrix& matrix);

        //- Construct from an lduMatrix with global columns and halo entries
        ldu2csr
        (
            const lduMatrix& matrix,
            const label rowStart,
            const List<Entry>& extra
        );

        //- Disallow default bitwise copy construction
        ldu2csr(const ldu2csr&) = delete;


    //- Destructor
    ~ldu2csr() = default;


    // Member Functions

        // Access

            label nRows() const
            {
                return nRows_;
            }

            label nNonZeros() const
            {
                return nNonZeros_;
            }

            const label* rowPtr() const
            {
                return rowPtr_.begin();
            }

            const label* colInd() const
            {
                return colInd_.begin();
            }

            const scalar* values() const
            {
                return values_.begin();
            }

            scalar* values()
            {
                return values_.begin();
            }


        //- Recompute the values (serial; structure regenerated from the
        //  matrix, expected to be identical for an unchanged mesh)
        void updateValues(const lduMatrix& matrix);


    // Member Operators

        //- Disallow default bitwise assignment
        void operator=(const ldu2csr&) = delete;
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
