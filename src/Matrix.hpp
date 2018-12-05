//      Matrix.h
//      
//      Copyright 2013 Alessia Visconti <visconti@di.unito.it>
//      
//      This program is free software; you can redistribute it and/or modify
//      it under the terms of the GNU General Public License as published by
//      the Free Software Foundation; either version 3 of the License, or
//      (at your option) any later version.
//      
//      This program is distributed in the hope that it will be useful,
//      but WITHOUT ANY WARRANTY; without even the implied warranty of
//      MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
//      GNU General Public License for more details.
//      
//      You should have received a copy of the GNU General Public License
//      along with this program; if not, write to the Free Software
//      Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
//      MA 02110-1301, USA.


#ifndef MATRIX_H
#define MATRIX_H

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <functional>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string.h>
#include <vector>
#include <limits>
#include <float.h>

#include "utilities.h"

using namespace std;



/**
	\brief Matrix class. 
	
	It represent a gene expression matrix as a float matrix, where each row represent a gene, each column represent a condition and cells represent expression level values.
	 
 */

class Matrix {

private:

	floatmatrix m;
	floatvect readRow(string row); 
	floatmatrix readVector(istream &is); 
	
	

public:

/**
	\brief  Return an empty and uninitialized matrix.

	\return the matrix
*/

	Matrix() {};

/**
	\brief  Return an initialized matrix.
	
	\param fm the matrix entries
	\return the matrix
*/
	
	Matrix(floatmatrix& fm)
	{
		m = fm;
	};

/**
	\brief Destructor.
*/

	~Matrix() {	};
	
/**
	\brief Return a matrix saved in filename.

	\param filename filepath
	\return the matrix
*/	
	
	void loadFromFile(char* filename);
	
/**
	\brief Return a string representing the matrix
	
	\return the string representing the matrix
*/
	
	string to_string() ;
	
/**
	\brief Return the matrix row number
	
	\return number of rows
*/
	
	unsigned int getRowsNumber();
	
/**
	\brief Return the matrix column number
	
	\return number of column
*/

	unsigned int getColumnsNumber();

	
/**
	\brief Return the element in position (i,j)
	
	\param i row index
	\param j column index
	\return a matrix element
*/	

	float getElement(int i, int j);

/**
	\brief Set the element in position (i,j) to value
	
	\param i row index
	\param j column index
	\param value the value to set
*/	
		
	void setElement(int i, int j, float value);
	
/**
	\brief Return a copy of the matrix
	
	\return a copy of the matrix
*/

	Matrix copy();


/**
	\brief Return the matrix-vector product
	
	\param fv the vector
	\return the product 
*/		
	
	floatvect vector_product(floatvect fv);	

//====================================================================
//            			ISA biclustering							//
//====================================================================

/**
	\brief Return the matrix traspose
	
	\return the trasposed matrix
*/	
	
	Matrix traspose();

/**
	\brief Return the normalized matrix, i.e. a matrix having  0 mean and 1 standard deviation. 
	
	\return the normalized matrix
*/	
	
	void normalize();
	





} ;

#endif
