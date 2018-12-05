//      utilities.h
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


#ifndef FUNCTION_H
#define FUNCTION_H

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
#include <map>


using namespace std;

//==================
//    Typedef
//==================

typedef vector< vector <float > > floatmatrix;
typedef vector< vector <int > > intmatrix;
typedef vector<float> floatvect;
typedef vector<float>::iterator floatvect_it;
typedef vector<int> intvect;
typedef vector<int>::iterator intvect_it;
typedef vector<string> stringvect;
typedef vector< vector <string > > stringmatrix;
typedef vector< map<string, int> > hashvector;
typedef map<string, int> hash;



//==================
//    Constants
//==================

// --- ISA ---
static const int seed_ratio = 10; //percentage of gene belonging to initial random seed w.r.t. gene pools
static const float uniform_score = 1.0; //starting value of genes belonging to the initial random seed
static const int max_isa_runs = 100;  //number beyond which AID-ISA diverges

//Following values have been choosen according to [Ihmels et al, 2002] and [Ihmels et al, 2004] 
//Condition thresholds vary according to the R package "eisa" by G. Csárdi. In the original paper it was fixed to 2.0

static const float min_gene_threshold = 1.8; //the min gene threshold
static const float max_gene_threshold = 4.0; //the max gene threshold
static const float gene_threshold_step = 0.1; //the step for gene threshold computation
static const float min_condition_threshold = 2.0; //the min condition threshold
static const float max_condition_threshold = 2.0; //the max condition threshold
static const float condition_threshold_step = 0.5; //the step for condition threshold computation



//==================
//    Functions Signature
//==================


/**
	\brief Return the mean of a vector
	 
	\param v the vector
	\return the mean 
*/


	static float vect_mean(floatvect v)
	{
		return accumulate(v.begin(), v.end(), 0.0)/v.size();
	}

/**
	\brief Return the variance of a vector
	 
	\param v the vector
	\return the variance 
*/


	static float vect_variance(floatvect v)
	{
		float mean = vect_mean(v);
		float sum = 0.0;
		for(unsigned int i=0; i<v.size(); i++)
			sum += pow((v[i] - mean), 2);
		
		return sum/(v.size()-1);
	}

/**
	\brief Return the standartd deviation of a vector
	 
	\param v the vector
	\return the standard deviation
*/

	static float vect_std(floatvect v)
	{
			return sqrt(vect_variance(v));
	}

/**
	\brief Return a random value included in [0, max]
	 
	\param max_value the max value
	\return a random value
*/
	static int random_value(int max_value)
	{
		srand (time(NULL));
		return rand()%max_value;
	}
	

/**
	\brief Return a list of object saved in filename.

	\param filename filepath
	\return the list of objects
*/
	static stringvect loadListFromFile(char* filename)
	{
		stringvect s;
		filebuf fb;
		fb.open (filename, ios::in);
		istream is(&fb);
		if (!is) return s;
		
		string line;
	    while (getline(is, line))
			s.push_back(line);
			
		fb.close();
		return s;
	}


/**
	\brief Save to a string in a filepath
	
	\param filepath the filepath
	\param msg the string
*/
	static void printToFile(const char* filepath, string msg)
	{
		ofstream file;
		file.open (filepath);
		file << msg;
		file.close();
	}


/**
	\brief Returns a list of label where the prefix string is followed by an increasing integer number up to max_index
	
	\param prefix the label prefix
	\param max_index the max label value
	\return the list of labels
*/

	static vector<string> get_mock_elements(const string &prefix, unsigned int max_index)
	{
		vector<string> vect;
		for(unsigned int i=0; i<max_index; i++)
		{
			stringstream el;
			el << prefix << i;
			vect.push_back(el.str());
		}
		return vect;
	}

/**
	\brief Returns whether a value is contained in vector
	
	\param vector the vector
	\param value the value
	\return true idf value is contained in vector, false otherwise
*/

	static bool include(intvect vector, int value)
	{
		intvect_it it = find(vector.begin(), vector.end(), value);
		if (it == vector.end()) return false;
		return true;
	}



//---------------------------- helpers ---------------------------------

/**
	\brief Return a vector that contains objects (labels) in a list that are marked as 'present' in a index vector
	
	It is a helper for the bicluster to_human_string() method.
	
	\param index the index vector
	\param list the list of objects
	\return a vector of objects
*/
	static stringvect intvectToStringvect(intvect index, stringvect list)
	{
		stringvect s;
		for(unsigned int j=0; j<index.size(); j++) 
			s.push_back(list[index[j]]);
		return s;
	}
/**
	\brief Return a string codifying a vector of objects
	
	It is a helper for the bicluster to_human_string() method.
	
	\param v the vector
	\return the string 
*/	
	static string stringvectToString(stringvect v)
	{
		ostringstream output;
		for (unsigned int j=0; j<v.size(); j++) output << v[j] << "\t";
		return output.str();
	}

/**
	\brief Return a string codifying a vector of floats
	
	It is a helper for the Matrix to_string method.
	
	\param v the vector
	\return the string 
*/	

	static string floatvectToString(floatvect v)
	{
		ostringstream output;
		for (unsigned int j=0; j<v.size(); j++) output << v[j] << "\t";
		return output.str();
	}

	


	


#endif
