//      AID-ISA.cpp
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

#include <boost/program_options.hpp>
#include "sysexits.h"

#include <cstdlib>
#include <iostream>
#include <vector>
#include <fstream>
#include <string.h>
#include <sstream>

#include "utilities.h"
#include "Matrix.hpp"
#include "Bicluster.hpp"
#include "Driver.hpp"

using namespace std;
using namespace boost::program_options;

typedef vector<Bicluster> Biclustervect;



int main(int argc, char** argv)
{
	
	/*
	 * List of command line parameters, and object used throughtout 
	 * the code
	 */ 

	
	string input_filename;
	string output_filename;
	bool if_row_driver;
	bool if_col_driver;
	unsigned int runs_number;
	float delta_expand;
	float delta_reduce;
	string gene_filename;
	string condition_filename;
	string gene_driver_filename;
	string condition_driver_filename;
	
	Matrix E;
	Driver gene_driver;
	Driver condition_driver;		    
	stringvect geneList;
	stringvect conditionList;

	/*
	 * Define and parse the command line parameters
	 * Data are loaded, including additional information sources.
	 *  
	 * If no list of row/col name are load, mock names
	 * are created.
	 */
	
	try
    {	
		options_description generic("Generic options");
		generic.add_options()
			("help,h", "produce help message and exit")
			("version,v", "display version information and exit");
			
		options_description mandatory("Parameters");
		mandatory.add_options()
			("input,i", value<string>(&input_filename),  "gene expression data input filepath");
			
		options_description parameter("Optional parameters");
		parameter.add_options()
			("gene_ida?,G", value<bool>(&if_row_driver)->default_value(false), "whether additional information are provided for gene dimension")
			("condition_ida?,C", value<bool>(&if_col_driver)->default_value(false), "whether additional information are provided for condition dimension")
			("gene_information,g", value<string>(&gene_driver_filename), "additional information for gene dimension")
			("condition_information,c", value<string>(&condition_driver_filename), "additional information for condition dimension")
			("output,o", value<string>(&output_filename),  "AID-ISA output filepath")
			("runs,n", value<unsigned int>(&runs_number)->default_value(10), "number of random initial seeds to use")
			("d_reduction,r", value<float>(&delta_reduce)->default_value(2.0), "delta for AID reduction step")
			("d_expansion,e", value<float>(&delta_expand)->default_value(0.5), "delta for AID expansion step")
			("gene_labels,x", value<string>(&gene_filename ),  "gene labels")
			("condition_label,y", value<string>(&condition_filename),  "condition labels");
        
        options_description cmdline_options;
        cmdline_options.add(generic).add(mandatory).add(parameter);
            	
		positional_options_description pos;
		pos.add("input", -1);
        
        variables_map vm;
        store(command_line_parser(argc, argv).options(cmdline_options).positional(pos).run(), vm);
		notify(vm);    
		
		//--------------------------------------------------------------
		
		if (vm.count("help")) 
		{
			cout << "Usage: AID-ISA input gene_ida? condition_ida? [gene_information, condition_information, output, runs, d_reduction, d_expansion, gene_labels, condition_labels]" << endl << cmdline_options << endl;
			cout << endl << "If gene_ida? is true gene_information MUST be supplied" << endl;
			cout << "If condition_isa? is true  condition_information MUST be supplied" << endl;
			
			return EX_OK;
		}
		 
		if (vm.count("version"))
		{
			cout << "AID-IDA version 1.0 -- June 2013" << endl;
			cout << "Please report bugs to <visconti@di.unito.it>" << endl;
			return EX_OK;
		}
		
		//--------------------------------------------------------------
		
		cout << endl << "###################   AID-ISA   ###################" << endl << endl;
		
		//read the mandatory
		if (!vm.count("input"))
		{
			cout << "Input file missed." << endl;
			cout << endl << "###################################################" << endl << endl;
			cout << "Usage: AID-ISA input [options]"  << endl << cmdline_options << endl;
			return EX_USAGE;
		}
		
		//Shall I use gene information?
		if(if_row_driver && !vm.count("gene_information"))
		{
			cerr << "ERROR: If gene_ida? is true gene_information MUST be supplied" << endl;
			cout << endl << "###################################################" << endl << endl;
			cout << "Usage: " << endl << cmdline_options << endl;
			return EX_USAGE;
		}
		else if(if_row_driver && vm.count("gene_information"))
		{
			cout << "Loading additional information (gene)..." << endl;	
			gene_driver.loadFromFile(const_cast<char *>(gene_driver_filename.c_str()));
			if (gene_driver.getRowsNumber() == 0)
			{
				cout << "ERROR: no additional information provided in file '" << gene_driver_filename << "'" << endl;
				return EX_DATAERR;
			}
			cout << "\t done." << endl;
		}

		//Shall I use condition information?
		if(if_col_driver && !vm.count("condition_information"))
		{
			cerr << "ERROR: If condition_ida? is true condition_information MUST be supplied" << endl; 
			cout << endl << "###################################################" << endl << endl;
			cout << "Usage: " << endl << cmdline_options << endl;
			return EX_USAGE;
		}
		else if(if_col_driver && vm.count("condition_information"))
		{
			cout << "Loading additional information (conditions)..." << endl;	
			condition_driver.loadFromFile(const_cast<char *>(condition_driver_filename.c_str()));
			if (condition_driver.getRowsNumber() == 0)
			{
				cout << "ERROR: no additional information provided in file '" << gene_driver_filename << "'" << endl;
				return EX_DATAERR;
			}
			cout << "\t done." << endl;			
		}
	
		//--------------------------------------------------------------
		
		if(!if_row_driver && !if_col_driver)
			cout << "The usage of additional information is not required: ISA algorithm will be evaluated" << endl << endl;
		
		
		//reading data (parameter already checked)
		cout << "Loading data..." << endl;	
	    E.loadFromFile(const_cast<char *>(input_filename.c_str()));
		
		if (E.getRowsNumber() == 0) 
		{
			cout << "ERROR: I cannot open the file '" << input_filename << "'" << endl;
			return EX_DATAERR;
		}
		cout << "\t done." << endl;


		
	
		//Gene labels read -if not present mock symbols are used 
		if (!vm.count("gene_labels"))
		{
			geneList = get_mock_elements("R", E.getRowsNumber());
		}
		else 
		{
			cout << "Loading labels (gene)..." << endl;	
			geneList = loadListFromFile(const_cast<char *>(gene_filename.c_str()));
			if (geneList.size() != E.getRowsNumber()) 
			{
				cout << "WARNING: I cannot open the file '" << gene_filename << "' or there are some problem with the gene list size" << endl;
				cout << "\t I will return the numeric signature" << endl;
				geneList = get_mock_elements("R", E.getRowsNumber());
			}
			cout << "\t done." << endl;
		}

		//Gene labels read -if not present mock symbol are used 
		if (!vm.count("condition_labels"))
		{
			conditionList = get_mock_elements("C", E.getColumnsNumber());
		}
		else 
		{
			cout << "Loading labels (condition)..." << endl;	
			conditionList = loadListFromFile(const_cast<char *>(condition_filename.c_str()));
			if (conditionList.size()!= E.getColumnsNumber()) 
			{
				cout << "WARNING: I cannot open the file '" << condition_filename << "' or there are some problem with the condition list size" << endl;
				cout << "\t I will return the numeric signature" << endl;
				conditionList = get_mock_elements("C", E.getColumnsNumber());
			}
			cout << "\t done." << endl;			
		}
		
		if (!vm.count("output"))
		{
			output_filename = input_filename + ".out";
			cout << "Data will be saved in " << output_filename << endl;
		}		
								
    }
	catch (unknown_option e)
    {
		cerr << e.what() << endl;
		return EX_USAGE;
    }
	catch (exception e)
	{
		cerr << e.what() << endl;
		return EX_SOFTWARE;
	}
	
	
	/*
	 * Data pre-processing: ISA works shows best performance if exploit 
	 * normalized matrices (one for the gene signature evaluation and 
	 * one for the condition signature evaluation)
	 */
	
	
	cout << endl << "Data pre-processing..." << endl;
	Matrix E_g = E.traspose();
	E_g.normalize();
	Matrix E_c = E.copy();
	E_c.normalize();
	cout << "\t done" << endl;
	
	/*
	 * Running AID-ISA
	 */
	
	Biclustervect results;

		//AID-ISA starts from a random sparse seed. 
		//In this way it is completly stochastic, and at each run it may give
		//different outputs for the same thresholds.
	
	cout << endl << "AID-ISA starts" << endl;
	for(unsigned int r=0; r<runs_number; r++)
	{
		cout << "\tRun: " << r << "/" << runs_number << endl;

		Bicluster initial_signature;
		initial_signature.initializeSignature(E.getRowsNumber());
			
		//the gene_threshold determine the resolution of the modular decomposition. 
		//By varying it, it is possible to discover multiple biclusters.
		float condition_threshold = min_condition_threshold;
		while(condition_threshold <= max_condition_threshold)
		{
			float gene_threshold = min_gene_threshold;
			while(gene_threshold <= max_gene_threshold)
			{
				//each seed will be evaluated on all the possible gene_threshold
				Bicluster signature = initial_signature.copy();
				signature.iterativeSignatureAlgorithm(E_g, E_c, gene_threshold, condition_threshold, gene_driver, condition_driver, delta_reduce, delta_expand, if_row_driver, if_col_driver);
				
				//void bicluster are discarded
				if (signature.getGeneCluster().size() != 0)
				{
					if (!signature.include(results)) //check if the evaluated bicluster is already known
					{
						Bicluster b = signature.copy();
						results.push_back(b);
					}
				}
	
				gene_threshold += gene_threshold_step;
			}
			condition_threshold += condition_threshold_step;
		}
	}
	
	/*
	 * Results are saved
	 */
	
	cout << endl << "Saving results..." << endl;

	for(unsigned int i=0; i<results.size(); i++)
	{
		ostringstream output_str;
		//bicluster size [row, col]
		output_str << "[" << results[i].getGeneCluster().size() << ", " << results[i].getConditionCluster().size() << "]" << endl;
		//map the signature (that are index!) in the name of genes/experiments	
		output_str << results[i].to_humanString(geneList, conditionList) << endl;
		printToFile(const_cast<char *>(output_filename.c_str()), output_str.str());
	}
	
	cout << "\t done" << endl;
	cout << endl << "###################################################" << endl << endl;
	return EX_OK;

}



