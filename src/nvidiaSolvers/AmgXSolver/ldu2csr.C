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

#include "ldu2csr.H"

#include <vector>
#include <algorithm>
#include <utility>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::ldu2csr::assemble
(
    const lduMatrix& matrix,
    const label rowStart,
    const List<Entry>& extra
)
{
    const label nCells = matrix.diag().size();
    const label nFaces = matrix.upper().size();

    const labelUList& owner = matrix.lduAddr().upperAddr();
    const labelUList& neighbour = matrix.lduAddr().lowerAddr();

    const scalarField& diag = matrix.diag();
    const scalarField& upper = matrix.upper();
    const scalarField& lower = matrix.lower();

    // Collect (column, value) entries per row
    std::vector<std::vector<std::pair<label, scalar>>> rows(nCells);

    for (label cell = 0; cell < nCells; cell++)
    {
        rows[cell].push_back(std::make_pair(rowStart + cell, diag[cell]));
    }

    for (label f = 0; f < nFaces; f++)
    {
        const label o = owner[f];
        const label n = neighbour[f];

        // A(owner, neighbour) = lower
        rows[o].push_back(std::make_pair(rowStart + n, lower[f]));

        // A(neighbour, owner) = upper
        rows[n].push_back(std::make_pair(rowStart + o, upper[f]));
    }

    // Off-processor (halo) entries, already with global column indices
    forAll(extra, i)
    {
        const Entry& e = extra[i];
        rows[e.row].push_back(std::make_pair(e.col, e.value));
    }

    // Sort columns and merge duplicates, then count non-zeros
    nRows_ = nCells;
    nNonZeros_ = 0;

    for (label cell = 0; cell < nCells; cell++)
    {
        std::vector<std::pair<label, scalar>>& r = rows[cell];

        std::sort
        (
            r.begin(),
            r.end(),
            [](const std::pair<label, scalar>& a, const std::pair<label, scalar>& b)
            {
                return a.first < b.first;
            }
        );

        std::vector<std::pair<label, scalar>> merged;
        merged.reserve(r.size());

        for (const std::pair<label, scalar>& e : r)
        {
            if (!merged.empty() && merged.back().first == e.first)
            {
                merged.back().second += e.second;
            }
            else
            {
                merged.push_back(e);
            }
        }

        rows[cell].swap(merged);
        nNonZeros_ += label(rows[cell].size());
    }

    // Flatten
    rowPtr_.setSize(nRows_ + 1);
    colInd_.setSize(nNonZeros_);
    values_.setSize(nNonZeros_);

    label idx = 0;
    rowPtr_[0] = 0;

    for (label cell = 0; cell < nCells; cell++)
    {
        for (const std::pair<label, scalar>& e : rows[cell])
        {
            colInd_[idx] = e.first;
            values_[idx] = e.second;
            idx++;
        }

        rowPtr_[cell + 1] = idx;
    }
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

Foam::ldu2csr::ldu2csr(const lduMatrix& matrix)
:
    nRows_(0),
    nNonZeros_(0),
    rowPtr_(),
    colInd_(),
    values_()
{
    assemble(matrix, 0, List<Entry>());
}


Foam::ldu2csr::ldu2csr
(
    const lduMatrix& matrix,
    const label rowStart,
    const List<Entry>& extra
)
:
    nRows_(0),
    nNonZeros_(0),
    rowPtr_(),
    colInd_(),
    values_()
{
    assemble(matrix, rowStart, extra);
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void Foam::ldu2csr::updateValues(const lduMatrix& matrix)
{
    assemble(matrix, 0, List<Entry>());
}


// ************************************************************************* //
