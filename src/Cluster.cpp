//      Cluster.cpp
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


#include "Cluster.hpp"

	
floatvect Cluster::getCluster()
{
	return values;
}

	
void Cluster::setValue(int i, float v)
{
	values[i] = v;
}	



float Cluster::getValue(int i)
{
	return values[i];
}
	
string Cluster::to_string()
{
	ostringstream output_index;
	ostringstream output_values;
	for(unsigned int i=0; i<values.size(); i++)
	{
		if (!(values[i] == 0.0))
		{
			output_index << i << "\t";
			output_values << values[i] << "\t";
		}
	}
	return (output_index.str() + "\n" + output_values.str()); 
}
	
bool Cluster::equal(Cluster c)
{
	floatvect a = this->values;
	floatvect b = c.values;
	
	bool result = true;
	unsigned int i=0;
	while (result && i<a.size())
	{
		if ((a[i] == 0.0 && b[i] != 0.0) || (a[i] != 0.0 && b[i] == 0.0)) result = false;
		i++;
	}
		
	return result;
}	
	
Cluster Cluster::copy()
{
	floatvect fv;
	for (unsigned int i=0; i<this->values.size(); i++)
		fv.push_back(this->values[i]);
	
	return Cluster(fv);
}	

intvect Cluster::getElements()
{
	intvect iv;
	for (unsigned int i=0; i<this->values.size(); i++)
		if (values[i] != 0.0) iv.push_back(i);	//element with 0 value does not belong to the cluster
		
	return iv;
}
	

//====================================================================
//            			ISA biclustering							//
//====================================================================


void Cluster::filter(float threshold, unsigned int n)
{
	float avg = vect_mean(this->values);
	float sigma = vect_std(this->values);
	
	if (sigma < 0.0001 && sigma > -0.0001) sigma = 1.0/sqrt(n); //random fluctuation
	
	threshold *= sigma;
	
	for (unsigned int i=0; i<this->values.size(); i++)
		if (fabs(this->values[i] - avg) < threshold) this->values[i] = 0.0;
		
}


void Cluster::average(unsigned int n)
{	
	for(unsigned int i=0; i<this->values.size(); i++)
		this->values[i] /= n; 
}

unsigned int Cluster::size()
{
	int n = 0;
	for(unsigned int i=0; i<this->values.size(); i++)
		if (this->values[i] != 0.0) n++;	
		
	return n;
}


Cluster Cluster::calculate(Matrix E, float threshold)
{
	unsigned int n = this->size();
	Cluster cluster;
	cluster.values = E.vector_product(this->values);
	cluster.average(n);
	cluster.filter(threshold, n);
	return cluster;
}

void Cluster::setRandomSeed(int n) 
{
	intvect v;
	intvect_it v_it;
	
	int max = this->values.size();
	
	int i = 0;
	while (i<n) 
	{
		int index = random_value(max);
		v_it = find(v.begin(), v.end(), index); //check for n *different* values
		if (v_it == v.end()) 
		{
			v.push_back(index);
			i++;
		}
	}
	
	//the score of all the elements of the random seed are inizialized to a uniform value
	for(unsigned int j=0; j<v.size(); j++)
		values[v[j]] = uniform_score;
}


void Cluster::resetValues(intvect iv)
{
	unsigned int n = this->values.size();
	this->values.resize(0);
	for(unsigned int i=0; i<n; i++)
		this->values.push_back(0.0);
	for(unsigned int i=0; i<iv.size(); i++)
		this->values[iv[i]] = 1.0;
}


//====================================================================
//            			Data Driven Activities						//
//====================================================================


float Cluster::compute_averange_distance(intvect index, Driver driver)
{
	unsigned int n = index.size();
	int count = 0; 
	
	float sum = 0.0;
	for(unsigned int i=0; i<n-1; i++)
		for(unsigned int j=i+1; j<n; j++) //it computes the distance between each object excluded itself and the driver is triangular
		{
			float value = driver.getElement(index[i], index[j]);
			if (value != -1)	//if no information was available, this distance is not considered
			{
				count ++;
				sum += value;
			}
		}	
		
	return sum/count;	
}


int Cluster::selectCentroid(intvect index, Driver driver)
{
	int min_distance = numeric_limits<int>::max();
	int centroid = -1;
	
	for (unsigned int i=0; i<index.size(); i++)
	{	
		float sum = 0.0;
		for(unsigned int j=0; j<index.size(); j++)
		{
			float value = driver.getElement(index[i], index[j]);
			if (value != -1)
				sum += value; //if no information was available, this distance is not considered
		}
		
		if (sum < min_distance) 
		{
			min_distance = sum;
			centroid = index[i];
		}
	}
	
	return centroid;
}


void Cluster::reduce(Driver driver, float reduce_coefficient)
{
	//it retains only the objects with a distance wrt the cluster centroid 
	//smaller or equal than the averange distance within all the cluster objects
	
	intvect index = this->getElements();
	
	if (index.size() < 2) return; //reduction is useless
	
	float thresold = compute_averange_distance(index, driver);
	int centroid = selectCentroid(index, driver);
	if (centroid == -1) return;
	
	intvect to_retain;
	
	to_retain.push_back(centroid);
	for(unsigned int i=0; i<index.size(); i++) //it checks only the objects belonging to the cluster
		if (driver.getElement(index[i], centroid) <= (thresold*reduce_coefficient)) to_retain.push_back(index[i]); //if no information is available, the distance is set to be -1, then the object is retained 'by default'
	
	this->resetValues(to_retain);
}

void Cluster::expand(Driver driver, float expand_coefficient)
{
	//it joins only the objects with a distance wrt the cluster centroid 
	//smaller or equal than the averange distance within all the cluster objects
	
	intvect index = this->getElements();
	
	if (index.size() == 0) return;
	float thresold = compute_averange_distance(index, driver);
	if (thresold < 0.0001 && thresold > -0.0001) return; //no information available for this cluster
	
	int centroid = selectCentroid(index, driver);
	if (centroid == -1) return;
	
	intvect to_retain;
	for(unsigned int d=0; d<driver.getRowsNumber(); d++) //it checks all the objects belonging to the data set
	{
		if (include(index, d)) to_retain.push_back(d); //it belongs already to the cluster
		else if (driver.getElement(d, centroid) <= (thresold*expand_coefficient) && driver.getElement(d, centroid) >= 0) to_retain.push_back(d); //if no information is available, the distance is set to be -1, and the object MUST NOT be added (if the second check is not performed, it'll added 'by default')
	}
	
	this->resetValues(to_retain);
}

void Cluster::drive(Driver driver, float reduce_coefficient, float expand_coefficient)
{
	this->reduce(driver, reduce_coefficient);
	this->expand(driver, expand_coefficient);
}

