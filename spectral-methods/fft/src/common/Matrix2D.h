#ifndef MATRIX2D_H
#define MATRIX2D_H

#include <iostream>
#include "config.h"

// ****************************************************************************
// Class:  Matrix2D
//
// Purpose:
//   Encapsulation of 2D matrices.
//
// Programmer:  Phil Roth
// Creation:    October 28, 2009
//
// ****************************************************************************
template<class T>
class Matrix2D
{
public:
    typedef T* restrict FlatDataPtr;
    typedef T* restrict* restrict DataPtr;
    typedef T* const restrict* const restrict ConstDataPtr;

private:
    size_t nRows;
    size_t nColumns;
    FlatDataPtr flatData;   // 1D array of data
    DataPtr data;           // data as 2D array (ptr to array of ptrs)

public:
    Matrix2D( size_t _nRows, size_t _nColumns )
      : nRows( _nRows ),
        nColumns( _nColumns ),
        flatData( new T[nRows * nColumns] ),
        data( new T*[nRows] )
    {
        for( size_t i = 0; i < nRows; i++ )
        {
            data[i] = &(flatData[i * nColumns]);
        }
    }

    ~Matrix2D( void )
    {
        delete[] data;
        data = NULL;

        delete[] flatData;
        flatData = NULL;
    }

    DataPtr GetData( void )
    {
        return data;
    }

    ConstDataPtr GetConstData( void ) const
    {
        return data;
    }
    
    FlatDataPtr GetFlatData( void )
    {
        return flatData;
    }

    size_t GetNumRows( void ) const { return nRows; }
    size_t GetNumColumns( void ) const { return nColumns; }

    size_t GetDataSize( void ) const { return nRows * nColumns * sizeof(T); }
};


template<class T>
std::ostream&
operator<<( std::ostream& s, const Matrix2D<T>& m )
{
    typename Matrix2D<T>::ConstDataPtr mdata = m.GetConstData();

    for( unsigned int i = 0; i < m.GetNumRows(); i++ )
    {
        for( unsigned int j = 0; j < m.GetNumColumns(); j++ )
        {
            if( j != 0 )
            {
                s << '\t';
            }
            s << mdata[i][j];
        }
        s << '\n';
    }
    return s;
}

#endif /* MATRIX2D_H */
