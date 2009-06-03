//
// File: bppConsense.cpp
// Created by: Julien Dutheil
// Created on: Jun Wed 06 11:17 2007
//

/*
Copyright or © or Copr. CNRS

This software is a computer program whose purpose is to estimate
phylogenies and evolutionary parameters from a dataset according to
the maximum likelihood principle.

This software is governed by the CeCILL  license under French law and
abiding by the rules of distribution of free software.  You can  use, 
modify and/ or redistribute the software under the terms of the CeCILL
license as circulated by CEA, CNRS and INRIA at the following URL
"http://www.cecill.info". 

As a counterpart to the access to the source code and  rights to copy,
modify and redistribute granted by the license, users are provided only
with a limited warranty  and the software's author,  the holder of the
economic rights,  and the successive licensors  have only  limited
liability. 

In this respect, the user's attention is drawn to the risks associated
with loading,  using,  modifying and/or developing or reproducing the
software by the user in light of its specific status of free software,
that may mean  that it is complicated to manipulate,  and  that  also
therefore means  that it is reserved for developers  and  experienced
professionals having in-depth computer knowledge. Users are therefore
encouraged to load and test the software's suitability as regards their
requirements in conditions enabling the security of their systems and/or 
data to be ensured and,  more generally, to use and operate it in the 
same conditions as regards security. 

The fact that you are presently reading this means that you have had
knowledge of the CeCILL license and that you accept its terms.
*/

// From the STL:
#include <iostream>
#include <iomanip>

using namespace std;

// From PhylLib:
#include <Phyl/Tree.h>
#include <Phyl/Newick.h>
#include <Phyl/PhylogeneticsApplicationTools.h>

// From Utils:
#include <Utils/AttributesTools.h>
#include <Utils/ApplicationTools.h>
#include <Utils/FileTools.h>
#include <Utils/TextTools.h>

using namespace bpp;

void help()
{
  *ApplicationTools::message << "__________________________________________________________________________" << endl;
  *ApplicationTools::message << "bppconsense parameter1_name=parameter1_value" << endl;
  *ApplicationTools::message << "      parameter2_name=parameter2_value ... param=option_file" << endl;
  *ApplicationTools::message << endl;
  *ApplicationTools::message << "  Refer to the Bio++ Program Suite Manual for a list of available options." << endl;
  *ApplicationTools::message << "__________________________________________________________________________" << endl;
}

int main(int args, char ** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*       Bio++ Consensus and Bootstrap Methods, version 0.3.0     *" << endl;
  cout << "* Authors: J. Dutheil                       Created     06/06/07 *" << endl;
  cout << "*          N. Galtier                       Last Modif. 02/06/09 *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  if(args == 1)
  {
    help();
    exit(0);
  }
  
  try {
  
  ApplicationTools::startTimer();

  cout << "Parsing options:" << endl;
  
  map<string, string> params = AttributesTools::parseOptions(args, argv);
 
  vector<Tree*> list = PhylogeneticsApplicationTools::getTrees(params);

  Tree* tree = 0;
  string treeMethod = ApplicationTools::getStringParameter("tree", params, "consensus");
  if(treeMethod == "input")
  {
    tree = PhylogeneticsApplicationTools::getTree(params);
    ApplicationTools::displayResult("Number of leaves", tree->getNumberOfLeaves());
  }
  else if(treeMethod == "consensus")
  {
    double threshold = ApplicationTools::getDoubleParameter("tree_consensus.threshold", params, 0);
    ApplicationTools::displayResult("Consensus threshold", TextTools::toString(threshold));
    ApplicationTools::displayTask("Computing consensus tree");
    tree = TreeTools::thresholdConsensus(list, threshold, true);
    ApplicationTools::displayTaskDone();
  }
  else throw Exception("Unknown input tree method: " + treeMethod);
  
  ApplicationTools::displayTask("Compute bootstrap values");
  TreeTools::computeBootstrapValues(*tree, list);
  ApplicationTools::displayTaskDone();

  //Write resulting tree:
  PhylogeneticsApplicationTools::writeTree(*tree, params);
  for(unsigned int i = 0; i < list.size(); i++) delete list[i];
  delete tree;

  cout << "BppConsense's done. Bye." << endl;
  ApplicationTools::displayTime("Total execution time:");

	}
  catch(exception & e)
  {
		cout << e.what() << endl;
		exit(-1);
	}

	return (0);
}

