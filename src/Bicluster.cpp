//      Bicluster.cpp
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

#include "Bicluster.hpp"

Cluster Bicluster::getGeneCluster()
{
	return gene;
}

Cluster Bicluster::getConditionCluster()
{
	return condition;
}

void Bicluster::setGeneCluster(Cluster& g)
{
	gene = g;
}

void Bicluster::setConditionCluster(Cluster& c)
{
	condition = c;
}

//TODO : ameliorate me, please!
string Bicluster::to_humanString(stringvect geneList, stringvect conditionList)
{
	intvect g_index;
	for (unsigned int i=0; i<this->gene.getCluster().size(); i++)
		if (this->gene.getValue(i) != 0.0) g_index.push_back(i);
	
	stringvect geneHumanCluster = intvectToStringvect(g_index, geneList);
	
	intvect c_index;
	for (unsigned int i=0; i<this->condition.getCluster().size(); i++)
		if (this->condition.getValue(i) != 0.0) c_index.push_back(i);
	
	stringvect conditionHumanCluster  = intvectToStringvect(c_index, conditionList);
	
	ostringstream output;
	output << stringvectToString(geneHumanCluster) << endl; 
	output << stringvectToString(conditionHumanCluster) << endl; 
	return output.str();
}

string Bicluster::to_string()
{
	ostringstream output;
	output << gene.to_string() << endl;
	output << condition.to_string() << endl;
	return output.str();
}

Bicluster Bicluster::copy()
{
	Cluster g = this->gene.copy();
	Cluster c = this->condition.copy();
	
	return Bicluster(g, c);
}

bool Bicluster::equal(Bicluster b)
{
	return (this->gene.equal(b.gene)) && (this->condition.equal(b.condition));
}




bool Bicluster::include(vector<Bicluster> cv)
{
	bool found = false;
	unsigned int i = 0;
	while (!found && i < cv.size())
	{
		found = this->equal(cv[i]);
		i++;
	}
	
	return found;
}




//====================================================================
//            			AID-ISA biclustering						//
//====================================================================


void Bicluster::signatureAlgorithm(Matrix E_R, Matrix E_C, float r_threshold, float c_threshold, Driver row_driver, Driver col_driver,  float reduce_coefficient, float expand_coefficient, unsigned int dd_row, unsigned int dd_col)
{
	Cluster row = this->gene; //reference gene set
	Cluster col = row.calculate(E_R, c_threshold); //is the condition signature!
	if (dd_col) col.drive(col_driver, reduce_coefficient, expand_coefficient);
	row = col.calculate(E_C, r_threshold); //is the gene signature!
	if (dd_row) row.drive(row_driver, reduce_coefficient, expand_coefficient);
		
	this->gene = row;
	this->condition = col;
}


void Bicluster::iterativeSignatureAlgorithm(Matrix E_R, Matrix E_C, float r_threshold, float c_threshold, Driver row_driver, Driver col_driver,  float reduce_coefficient, float expand_coefficient, unsigned int dd_row, unsigned int dd_col)
{
	if (dd_row) this->gene.drive(row_driver, reduce_coefficient, expand_coefficient);
	bool loop = true;
	int i = 0;
	while(loop)
	{
		Cluster g_seed = this->gene;
		Cluster c_seed = this->condition;
		this->signatureAlgorithm(E_R, E_C, r_threshold, c_threshold, row_driver, col_driver, reduce_coefficient, expand_coefficient, dd_row, dd_col);
		i++;
		
		Cluster new_g_seed = this->gene;
		Cluster new_c_seed = this->condition;
		
		//solution found
		if (g_seed.equal(new_g_seed) && c_seed.equal(new_c_seed) ) loop = false; 
		if (i > max_isa_runs) //it diverges
		{
			//the algorithm returns a void bicluster
			Cluster void_cluster (this->gene.getCluster().size(), 0.0);
			this->gene = void_cluster; //if gene cluster if void, also condition cluster will be void
			loop = false; 
		}
	}
}



void Bicluster::initializeSignature(unsigned int num_genes)
{
	Cluster g(num_genes, 0.0); //it starts by using a void signature (codify by the value 0)
	this->gene = g;
	
	int seed_number = (num_genes/100.0)*seed_ratio;
	this->gene.setRandomSeed(seed_number);
}
