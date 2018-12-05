//      Cluster.h
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


#ifndef CLUSTER_H
#define CLUSTER_H

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

#include "utilities.h"
#include "Matrix.hpp"
#include "Driver.hpp"



using namespace std;



/**
	\brief Cluster class. 
	
	Define a cluster as index: if the index entry is set to zero the corespondant object does not belog to the cluster.
	Entries set to other values represent objects belonging to the cluster. 
	
 */

class Cluster {

private:

	floatvect values;

/**
	\brief Return objects (genes/conditions) that pass a statistical test.
	Objects having a value far from the mean more than threshold times the standar deviation are retained.
	
	If eveluated standard deviation is zero, standard deviation expected for random fluctuation is used.
	
	\param threshold distance in standard deviations
	\param n cluster signature size
*/

	void filter(float threshold, unsigned int n);	

/** 
	\brief  Return the mean of objects (genes/conditions) in a cluster signature
	
	\param n cluster signature size
*/

	void average(unsigned int n);


/**
	\brief Return the average object distance
	
	\param index cluster
	\param driver distance matrix
*/
	static float compute_averange_distance(intvect index, Driver driver);

/**
	\brief Return the cluster centroid.
	
	A cluster centroid is the object more similar to all the other objects belonging to the cluster, i.e., the one closest (minimun distance), to all the other objects belonging to the cluster
	
	\param index cluster
	\param driver distance matrix
*/	
	static int selectCentroid(intvect index, Driver driver); 


/**
	\brief Return the reduced cluster.
	
	It removes objects having a distance wrt the cluster centroid smaller or equal than reduce_threshold times the averange distance  within all the cluster objects.
	
	\param driver distance matrix
	\param reduce_coefficient reduction threshold
*/	
	void reduce(Driver driver, float reduce_coefficient); 
	
/**
	\brief Return the expanded cluster
	
	It adds objects having a distance wrt the cluster centroid smaller or equal than expand_threshold times the averange distance within all the cluster objects.
	
	\param driver distance matrix
	\param expand_coefficient expansion threshold
*/	
	void expand(Driver driver, float expand_coefficient); 
	
/**
	\brief Return a cluster where the values of indices in iv are set to 1.0, whilist other values are set to 0.0. 
	
	\param iv the index to set to 1.0
*/	

	void resetValues(intvect iv);

public:

/**
	\brief  Return an empty and uninitialized cluster.
	
	\return the cluster
*/

	Cluster() 	{};

/**
	\brief  Return a cluster with n elements set to an intial value.
	
	\param n cluster size
	\param initial_value the default value for cluster entries
	\return the cluster
*/
	
	Cluster(unsigned int n, float initial_value)
	{
		for(unsigned int i=0; i<n; i++)
			values.push_back(initial_value);
	}
	
/**
	\brief  Return an initialized cluster.
	
	\param fv the cluster entries
	\return the cluster
*/
	
	Cluster(floatvect& fv)
	{
		values = fv;
	};

/**
	\brief Destructor.
*/

	~Cluster() {	};

/**
	\brief Return the cluster.
	
	Each element index refer to an index in the original matrix. An entry set to zero indicates that the object does not belog to the cluster. Other values represent correlation and indicate objects that belong to the cluster.
	
	\return the cluster
*/

	
	floatvect getCluster();

	
/**
	\brief Set the cluster object i to the value v
	
	\param i object index
	\param v value
*/		
	void setValue(int i, float v);
	

	
/**
	\brief Return the value of the object i
	
	\param i object index
	\return object value
*/		
	float getValue(int i);

/**
	\brief Return cluster objects (indices)
	
	\return cluster objects
*/
	intvect getElements();

/**
	\brief Return the cluster size.
	Only objects with value not set to zero are counted.
	
	\return cluster size
*/	
	unsigned int size();


/**
	\brief Return a string representing the cluster
	The output is composed by two lines: 
	 - the first line contains the indices of objects belonging to the cluser;
	 - the second line contains their correlation values

	\return the string representing the cluster
*/	
	string to_string();

/**
   \brief Return whether two clusters are equal. 
   Two clusters are equal iff they include the same objects, despite object values
   
   \param c cluster to compare
   \return true, if the clusters are equal, false otherwise
*/ 

  	
	bool equal(Cluster c);
	
/**
	\brief Return a copy of the cluster
	
	\return a copy of the cluster
*/
	Cluster copy();
	
	

//====================================================================
//            			ISA biclustering							//
//====================================================================


/**
	\brief Return the cluster signature according to the SA algorithm [Ihmels et al., Nat Genet, 2002].
	
	\param E data matrix
	\param threshold objects threshold
	\return cluster signature
*/

	Cluster calculate(Matrix E, float threshold); 

/**
	\brief Return the initial seed according to the SA algorithm [Ihmels et al., Nat Genet, 2002].
	
	\param n number of objects belonging to the cluster
*/

	void setRandomSeed(int n);

	

//====================================================================
//            			Data Driven Activities						//
//====================================================================

/**
	\brief Return the results of a run of the AID algorithm [Visconti et al., Intelligent Data Analysis, 2013].
	
	\param driver distance matrix
	\param reduce_coefficient reduction threshold
	\param expand_coefficient expansion threshold
*/

	void drive(Driver driver, float reduce_coefficient, float expand_coefficient);

} ;

#endif
