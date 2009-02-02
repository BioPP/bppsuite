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
  *ApplicationTools::message << "input.list.file            | [path] toward multi-trees file (newick)      " << endl;
  *ApplicationTools::message << "tree                       | [input|consensus] reference tree to use      " << endl;
  *ApplicationTools::message << "  tree_input.file          | [path] toward reference tree file            " << endl;
  *ApplicationTools::message << "  tree_consensus.threshold | [[0,1]] threshold to use in consensus        " << endl;
  *ApplicationTools::message << "                           | 0   = fully resolved                         " << endl;
  *ApplicationTools::message << "                           | 0.5 = majority rule                          " << endl;
  *ApplicationTools::message << "                           | 1   = strict                                 " << endl;
  *ApplicationTools::message << "                           | or any other value                           " << endl;
  *ApplicationTools::message << "___________________________|______________________________________________" << endl;
  PhylogeneticsApplicationTools::printOutputTreeHelp();
}

int main(int args, char ** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*       Bio++ Consensus and Bootstrap Methods, version 0.2.0     *" << endl;
  cout << "* Authors: J. Dutheil                       Created     06/06/07 *" << endl;
  cout << "*          N. Galtier                       Last Modif. 02/02/09 *" << endl;
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

  Newick newick;
  string listPath = ApplicationTools::getAFilePath("input.list.file", params);
  ApplicationTools::displayResult("Input list file", listPath);
  if(listPath == "none") throw Exception("You must provide an input tree list file.");
  
  vector<Tree *> list;
  newick.read(listPath, list);
  ApplicationTools::displayResult("Number of trees found", TextTools::toString(list.size()));

  Tree * tree = NULL;
  string treeMethod = ApplicationTools::getStringParameter("tree", params, "consensus");
  if(treeMethod == "input")
  {
    string treePath = ApplicationTools::getAFilePath("tree_input.file", params);
    if(treePath == "none") throw Exception("You must provide a file name for 'tree_input.file'.");
    ApplicationTools::displayResult("Input tree file", treePath);
    tree = newick.read(treePath);
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

