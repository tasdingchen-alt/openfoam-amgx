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

#include <algorithm>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

void Foam::ldu2csr::gather
(
    const lduMatrix& matrix,
    const label rowStart,
    const List<Entry>& extra,
    std::vector<Row>& rows
) const
{
    const label nCells = matrix.diag().size();
    const label nFaces = matrix.upper().size();

    const labelUList& owner = matrix.lduAddr().upperAddr();
    const labelUList& neighbour = matrix.lduAddr().lowerAddr();

    const scalarField& diag = matrix.diag();
    const scalarField& upper = matrix.upper();
    const scalarField& lower = matrix.lower();

    rows.resize(nCells);

    // Diagonal entries
    for (label cell = 0; cell < nCells; cell++)
    {
        rows[cell].push_back(std::make_pair(rowStart + cell, diag[cell]));
    }

    // Off-diagonal entries of the internal faces
    for (label face = 0; face < nFaces; face++)
    {
        rows[owner[face]].push_back
        (
            std::make_pair(rowStart + neighbour[face], lower[face])
        );
        rows[neighbour[face]].push_back
        (
            std::make_pair(rowStart + owner[face], upper[face])
        );
    }

    // Extra (halo) entries already carry global column indices
    forAll(extra, i)
    {
        const Entry& e = extra[i];
        rows[e.row].push_back(std::make_pair(e.col, e.value));
    }
}


void Foam::ldu2csr::flatten(std::vector<Row>& rows)
{
    nRows_ = rows.size();
    nNonZeros_ = 0;

    for (label cell = 0; cell < nRows_; cell++)
    {
        Row& row = rows[cell];

        std::sort
        (
            row.begin(),
            row.end(),
            [](const std::pair<label, scalar>& a, const std::pair<label, scalar>& b)
            {
                return a.first < b.first;
            }
        );

        // Merge duplicates that point at the same column
        Row merged;
        merged.reserve(row.size());

        for (const std::pair<label, scalar>& entry : row)
        {
            if (!merged.empty() && merged.back().first == entry.first)
            {
                merged.back().second += entry.second;
            }
            else
            {
                merged.push_back(entry);
            }
        }

        row.swap(merged);
        nNonZeros_ += row.size();
    }

    // Pack the rows into the CSR arrays
    rowPtr_.setSize(nRows_ + 1);
    colInd_.setSize(nNonZeros_);
    values_.setSize(nNonZeros_);

    label idx = 0;
    rowPtr_[0] = 0;

    for (label cell = 0; cell < nRows_; cell++)
    {
        for (const std::pair<label, scalar>& entry : rows[cell])
        {
            colInd_[idx] = entry.first;
            values_[idx] = entry.second;
            idx++;
        }

        rowPtr_[cell + 1] = idx;
    }
}


// * * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * //

Foam::ldu2csr::ldu2csr(const lduMatrix& matrix)
:
    nRows_(0),
    nNonZeros_(0),
    rowPtr_(),
    colInd_(),
    values_()
{
    std::vector<Row> rows;
    gather(matrix, 0, List<Entry>(), rows);
    flatten(rows);
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
    std::vector<Row> rows;
    gather(matrix, rowStart, extra, rows);
    flatten(rows);
}


// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

void Foam::ldu2csr::updateValues(const lduMatrix& matrix)
{
    std::vector<Row> rows;
    gather(matrix, 0, List<Entry>(), rows);
    flatten(rows);
}


// ************************************************************************* //
