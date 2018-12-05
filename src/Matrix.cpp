//      Matrix.cpp
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


#include "Matrix.hpp"

floatvect Matrix::readRow(string row) 
{
  floatvect retval;
  istringstream is(row);
  float num;
  while (is >> num) retval.push_back(num);
  return retval;
}

floatmatrix Matrix::readVector(istream &is) 
{
  string line;
  floatmatrix retval;
  while (getline(is, line))
    retval.push_back(readRow(line));
  return retval;
}


void Matrix::loadFromFile(char* filename) 
{
	filebuf fb;
	fb.open (filename, ios::in);
	istream is(&fb);
	if (!is) return;
	
	m = readVector(is);
		
	fb.close();
	return;
}


unsigned int Matrix::getRowsNumber() 
{
	return m.size();
}
	

unsigned int Matrix::getColumnsNumber() 
{
	return m[0].size();
}


float Matrix::getElement(int i, int j)
{
	return m[i][j];
}


void Matrix::setElement(int i, int j, float value)
{
	m[i][j] = value;
}


string Matrix::to_string()
{
	ostringstream output;
	for(unsigned int i=0; i<m.size(); i++)
		output << floatvectToString(m[i]) << endl;
	return output.str();
}


Matrix Matrix::copy()
{
	floatmatrix fm;
	for (unsigned int i=0; i<this->m.size(); i++)
	{
		floatvect fv;
		for (unsigned int j=0; j<this->m[0].size(); j++)
			fv.push_back(m[i][j]);
		fm.push_back(fv);
	}
	
	return Matrix(fm);
}



//====================================================================
//            			ISA biclustering							//
//====================================================================



Matrix Matrix::traspose() 
{
	floatmatrix n;
	for(unsigned int j=0; j<m[0].size(); j++) {
		floatvect v;
		for(unsigned int i=0; i<m.size(); i++) 
			v.push_back(m[i][j]);
		n.push_back(v);
	}	
	
	return Matrix(n);
}



void Matrix::normalize()
{
	floatvect tmp;
	for(unsigned int i=0; i<m.size(); i++) 
		for(unsigned int j=0; j<m[i].size(); j++)
			tmp.push_back(m[i][j]);
	
	float avg = vect_mean(tmp);
	float var = vect_variance(tmp);
		
	for (unsigned int i=0; i<m.size(); i++) 
		for(unsigned int j=0; j<m[i].size(); j++)	
			m[i][j] = (m[i][j] - avg)/var;
}

floatvect Matrix::vector_product(floatvect fv)
{
	floatmatrix fm = this->m;
	floatvect rv;
	for (unsigned int i=0; i<fm.size(); i++) 
		rv.push_back(0.0);
	
	for (unsigned int i=0; i<fm.size(); i++)
       for (unsigned int j=0; j<fm[i].size(); j++)
			rv[i] += fm[i][j]*fv[j];
			
    return rv;       
}

