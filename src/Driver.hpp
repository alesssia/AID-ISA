//      Driver.h
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


#ifndef DRIVER_H
#define DRIVER_H

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
#include "Matrix.hpp"

using namespace std;




/**
	\brief Driver class. 
	
	Define a driver as a float matrix, where each row/column represent an object and cells represent distances.
	A driver contains the additional information used by AID algorithm [Visconti et al., Intelligent Data Analysis, 2013].
	 
 */

class Driver {

private:

	floatmatrix m;
	
	static floatmatrix readFloatVector(istream &is);
	static floatvect readFloatRow(string row);
	
	static stringvect readStringRow(string row);
	static stringmatrix readStringVector(istream &is);


/**
	\brief Return the normalized driver

*/	
	void normalize(); 	
	


public:

/**
	\brief  Return an empty and uninitialized driver.

	\return the driver
*/

	Driver() {};

/**
	\brief  Return an initialized driver.
	
	\param matrix the driver values
	\return the driver
*/	
	
	Driver(floatmatrix matrix)
	{
		m = matrix;
	}
	
/**
	\brief Destructor.
*/	
	
	~Driver() {};

/**
	\brief Return the driver row number
	
	\return number of rows
*/
	unsigned int getRowsNumber();
	
/**
	\brief Return the driver column number
	
	\return number of columns
*/	
	unsigned int getColumnsNumber();
	
/**
	\brief Return the value in position (i,j)
	
	\param i row index
	\param j column index
	\return a driver element
*/	
	
	
	float getElement(int i, int j);

/**
	\brief Return a driver saved in filename.

	\param filename filepath
	\return the driver
*/	
	void loadFromFile(char* filename); 
	
	
/**
	\brief Return a string representing the driver
	
	\return the string representing the drivers
*/		
	
	string to_string();



	


} ;

#endif
