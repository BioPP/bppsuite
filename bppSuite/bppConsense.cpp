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

#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/KeyvalTools.h>

// From PhylLib:
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>

using namespace bpp;

void help()
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
  (*ApplicationTools::message << "bppconsense parameter1_name=parameter1_value").endLine();
  (*ApplicationTools::message << "      parameter2_name=parameter2_value ... param=option_file").endLine();
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "  Refer to the Bio++ Program Suite Manual for a list of available options.").endLine();
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
}

int main(int args, char ** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*       Bio++ Consensus and Bootstrap Methods, version 0.3.0     *" << endl;
  cout << "* Authors: J. Dutheil                       Created     06/06/07 *" << endl;
  cout << "*          N. Galtier                       Last Modif. 08/08/09 *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  if(args == 1)
  {
    help();
    return 0;
  }
  
  try {
  
  BppApplication bppconsense(args, argv, "BppConsense");
  bppconsense.startTimer();

  vector<Tree*> list = PhylogeneticsApplicationTools::getTrees(bppconsense.getParams());

  Tree* tree = 0;
  string treeMethod = ApplicationTools::getStringParameter("tree", bppconsense.getParams(), "consensus");
  string cmdName;
  map<string, string> cmdArgs;
  KeyvalTools::parseProcedure(treeMethod, cmdName, cmdArgs);
  if(cmdName == "Input")
  {
    tree = PhylogeneticsApplicationTools::getTree(bppconsense.getParams());
    ApplicationTools::displayResult("Number of leaves", tree->getNumberOfLeaves());
  }
  else if(cmdName == "Consensus")
  {
    double threshold = ApplicationTools::getDoubleParameter("threshold", cmdArgs, 0);
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
  PhylogeneticsApplicationTools::writeTree(*tree, bppconsense.getParams());
  for (unsigned int i = 0; i < list.size(); i++)
    delete list[i];
  delete tree;

  bppconsense.done();

	}
  catch(exception & e)
  {
		cout << e.what() << endl;
		return 1;
	}

	return (0);
}

