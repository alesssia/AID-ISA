//      Bicluster.hpp
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


#ifndef BICLUSTER_H
#define BICLUSTER_H

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
#include <limits>

#include "utilities.h"
#include "Matrix.hpp"
#include "Cluster.hpp"



using namespace std;


/**
	\brief Bicluster class. 
	
	Define biclusters as a pair of Clusters: 
	gene represents the row cluster; 
	condition the column cluster.
	
	\see Cluster.hpp
 */


class Bicluster {


private:

	Cluster gene;
	Cluster condition;

/**
	\brief Return the results of a run of the AID-SA algorithm [Visconti et al., Intelligent Data Analysis, 2013].
	
	It identifies a transcriptional module given a normalized gene expression matrix and its traspose. 
	AID-SA starts from a gene cluster and it evaluates the correspondant condition cluster. Then, the obtained 
	condition cluster is thresholded and it is exploited to evaluate the new gene cluster. The gene cluster
	is also thresholded.

	\see Cluster.calculate, for a detailed description of these steps.
	
		
	\param E_R the transposed and normalized gene expression matrix
	\param E_C the normalized gene expression matrix
	\param r_threshold gene threshols (SA parameter)
	\param c_threshold condition threshold (SA parameter)
	\param row_driver distance matrix for the gene dimension
	\param col_driver distance matrix for the condition dimension
	\param reduce_coefficient reduction threshold (AID parameter)
	\param expand_coefficient expansion threshold (AID parameter)
	\param dd_row set if AID is performed on gene dimension
	\param dd_col set if AID is performed on condition dimension
*/	

	void signatureAlgorithm(Matrix E_R, Matrix E_C, float r_threshold, float c_threshold, Driver row_driver, Driver col_driver, float reduce_coefficient, float expand_coefficient, unsigned int dd_row, unsigned int dd_col);
	


public:


/**
	\brief  Return an empty and uninitialized bicluster.

	\return the bicluster
*/

	Bicluster() {};

/**
	\brief  Returns an initialized bicluster.
	
	\param g row dimension
	\param c column dimension
	\return the bicluster
*/
	
	Bicluster(Cluster& g, Cluster& c)
	{
		gene = g;
		condition = c;
	};
	

/**
	\brief Destructor.
*/
	
	~Bicluster() {	};

/**
	\brief Return cluster on the gene dimension 
	
	\return gene cluster
*/
	Cluster getGeneCluster();

/**
	\brief  Return cluster on the condition dimension 
	
	\return condition cluster
*/
	Cluster getConditionCluster();
	
/**
	\brief Set the cluster on the gene dimension to g
	
	\param g gene cluster to set
*/	
	void setGeneCluster(Cluster& g);
	
/**
	\brief Set the cluster on the condition dimension to c
	
	\param c condition cluster to set
*/
	void setConditionCluster(Cluster& c);
	
/**
	\brief Return a copy of the bicluster
	
	\return a copy of the bicluster
*/
	Bicluster copy();

/**
   \brief Return whether two biclusters are equal. 
   Two biclusters are equal iff their cluster are equals, i.e. their index in the dta matrix are the same.
   
   \param b bicluster to compare
   \return true, if the biclusters are equal, false otherwise
*/   
	bool equal(Bicluster b);
	

/**
	\brief Return whether a vector of Biclusters contains a given Bicluster
	 
	\param cv vector of biclusters
	\return true if the vectors contains the bicluster, false otherwise
*/

	bool include(vector<Bicluster> cv);
	
/**
	\brief Return a string representing the bicluster
	The string is composed by four lines:
	the first (third) line represents the indices of the genes (conditions) that belong to the bicluster;
	the second (forth) line represents the correlation scores.

	\return the string representing the bicluster
*/
	string to_string(); 

/**
	\brief Return a string representing the bicluster using the label geneList and conditionList
	The output is composed by four lines, two for each cluster.
	The string is composed by two lines: the first one lists the gene belonging to the bicluster; the second line lists the conditions belonging to the biclusters.

	\param geneList gene labels
	\param conditionList condition labels
	\return the string representing the bicluster
*/	
	string to_humanString(stringvect geneList, stringvect conditionList);




	
//====================================================================
//            			AID-ISA biclustering						//
//====================================================================

/**
	\brief Return the initial gene seed for ISA. 
	
	First, it sets the gene signature to have an entry for each gene in the matrix. 
	Note that it is not necessary to set the dimension of condition signature: it will be evaluated in each SA run. Second, it randomly select a number of gene indices and set them to an uniform value. 
	
	Both the number of genes and the uniform value are global parameter.
	
	\param num_genes number of genes in the data set
*/
	void initializeSignature(unsigned int num_genes);
	

/**
	\brief Return the results of the AID-ISA algorithm [Visconti et al., Intelligent Data Analysis, 2013].
	
	It evaluates the AID-SA algorithm until the convegence criteria is reached (i.e., the element in both the gene and the condition does not change) or until the initial seed is proved to be divergent (i.e. the number of
	iteration exceded a global parameter). 
	In the latter case a void bicluser is returned.
	
	\param E_R the transposed and normalized gene expression matrix
	\param E_C the normalized gene expression matrix
	\param r_threshold gene threshols (SA parameter)
	\param c_threshold condition threshold (SA parameter)
	\param row_driver distance matrix for the gene dimension
	\param col_driver distance matrix for the condition dimension
	\param reduce_coefficient reduction threshold (AID parameter)
	\param expand_coefficient expansion threshold (AID parameter)
	\param dd_row set if AID is performed on gene dimension
	\param dd_col set if AID is performed on condition dimension
*/	

	void iterativeSignatureAlgorithm(Matrix E_R, Matrix E_C, float r_threshold, float c_threshold, Driver row_driver, Driver col_driver, float reduce_coefficient, float expand_coefficient, unsigned int dd_row, unsigned int dd_col);



 



	

} ;

#endif
