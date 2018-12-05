//      Driver.cpp
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


#include "Driver.hpp"

floatvect Driver::readFloatRow(string row) 
{
  floatvect retval;
  istringstream is(row);
  float num;
  while (is >> num) retval.push_back(num);
  return retval;
}

floatmatrix Driver::readFloatVector(istream &is) 
{
  string line;
  floatmatrix retval;
  while (getline(is, line))
    retval.push_back(readFloatRow(line));
  return retval;
}


void Driver::loadFromFile(char* filename) 
{
	filebuf fb;
	fb.open (filename, ios::in);
	istream is(&fb);
	if (!is) return;
	
	m = readFloatVector(is);
		
	fb.close();
	return;
}


unsigned int Driver::getRowsNumber() 
{
	return m.size();
}
	
unsigned int Driver::getColumnsNumber() 
{
	return m[0].size();
}


float Driver::getElement(int i, int j)
{
	return m[i][j];
}



string Driver::to_string()
{
	ostringstream output;
	for(unsigned int i=0; i<m.size(); i++)
	{
		for(unsigned int j=0; j<m[i].size(); j++) 
		{	
			output << m[i][j];
			output << "\t";
		}
		output << "\n";
	}
	return output.str();
}



void Driver::normalize()
{
	float max = FLT_MIN;
	for(unsigned int i=0; i<m.size(); i++)
		for(unsigned int j=0; j<m[i].size(); j++)
			if (m[i][j] > max) max = m[i][j];
	
	for(unsigned int i=0; i<m.size(); i++)
		for(unsigned int j=0; j<m[i].size(); j++)
			m[i][j] = m[i][j]/max;
}
