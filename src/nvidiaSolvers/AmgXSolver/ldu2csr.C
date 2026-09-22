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
    WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General
    Public License for more details.

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
        rows[cell].push_back(Item{rowStart + cell, diag[cell], 0, cell});
    }

    // Off-diagonal entries of the internal faces
    for (label face = 0; face < nFaces; face++)
    {
        rows[owner[face]].push_back
        (
            Item{rowStart + neighbour[face], lower[face], 1, face}
        );
        rows[neighbour[face]].push_back
        (
            Item{rowStart + owner[face], upper[face], 2, face}
        );
    }

    // Extra (halo) entries already carry global column indices
    forAll(extra, i)
    {
        const Entry& e = extra[i];
        rows[e.row].push_back(Item{e.col, e.value, 3, i});
    }
}


void Foam::ldu2csr::flatten(std::vector<Row>& rows)
{
    nRows_ = rows.size();
    nNonZeros_ = 0;

    bool duplicates = false;

    for (label cell = 0; cell < nRows_; cell++)
    {
        Row& row = rows[cell];

        std::sort
        (
            row.begin(),
            row.end(),
            [](const Item& a, const Item& b)
            {
                return a.col < b.col;
            }
        );

        // Merge duplicates that point at the same column
        Row merged;
        merged.reserve(row.size());

        for (const Item& entry : row)
        {
            if (!merged.empty() && merged.back().col == entry.col)
            {
                merged.back().val += entry.val;
                duplicates = true;
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

    cacheable_ = !duplicates;

    if (cacheable_)
    {
        srcType_.setSize(nNonZeros_);
        srcIndex_.setSize(nNonZeros_);
    }
    else
    {
        srcType_.clear();
        srcIndex_.clear();
    }

    label idx = 0;
    rowPtr_[0] = 0;

    for (label cell = 0; cell < nRows_; cell++)
    {
        for (const Item& entry : rows[cell])
        {
            colInd_[idx] = entry.col;
            values_[idx] = entry.val;

            if (cacheable_)
            {
                srcType_[idx] = entry.src;
                srcIndex_[idx] = entry.idx;
            }

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
    values_(),
    srcType_(),
    srcIndex_(),
    cacheable_(false)
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
    values_(),
    srcType_(),
    srcIndex_(),
    cacheable_(false)
{
    std::vector<Row> rows;
    gather(matrix, rowStart, extra, rows);
    flatten(rows);
}


// * * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * //

void Foam::ldu2csr::updateValues(const lduMatrix& matrix)
{
    updateValues(matrix, List<Entry>());
}


void Foam::ldu2csr::updateValues(const lduMatrix& matrix, const List<Entry>& extra)
{
    if
    (
        !cacheable_
     || srcType_.size() != nNonZeros_
     || rowPtr_.size() != nRows_ + 1
    )
    {
        // Structure is not reusable: rebuild everything
        std::vector<Row> rows;
        gather(matrix, 0, extra, rows);
        flatten(rows);
        return;
    }

    const scalarField& diag = matrix.diag();
    const scalarField& upper = matrix.upper();
    const scalarField& lower = matrix.lower();

    const label* srcType = srcType_.begin();
    const label* srcIndex = srcIndex_.begin();
    scalar* vals = values_.begin();
    const label nnz = nNonZeros_;

    #ifdef _OPENMP
    #pragma omp parallel for schedule(static)
    #endif
    for (label i = 0; i < nnz; i++)
    {
        switch (srcType[i])
        {
            case 0: vals[i] = diag[srcIndex[i]]; break;
            case 1: vals[i] = lower[srcIndex[i]]; break;
            case 2: vals[i] = upper[srcIndex[i]]; break;
            default: vals[i] = extra[srcIndex[i]].value; break;
        }
    }
}


// ************************************************************************* //
