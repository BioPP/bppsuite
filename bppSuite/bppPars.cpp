//
// File: bppPars.cpp
// Created by: Julien Dutheil
// Created on: May Sat 05 15:09 2007
// From file bppML.cpp
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

// From SeqLib:
#include <Seq/Alphabet.h>
#include <Seq/VectorSiteContainer.h>
#include <Seq/SiteTools.h>
#include <Seq/SequenceApplicationTools.h>

// From PhylLib:
#include <Phyl/Tree.h>
#include <Phyl/PatternTools.h>
#include <Phyl/PhylogeneticsApplicationTools.h>
#include <Phyl/DRTreeParsimonyScore.h>
#include <Phyl/OptimizationTools.h>
#include <Phyl/Newick.h>

// From Utils:
#include <Utils/BppApplication.h>
#include <Utils/ApplicationTools.h>
#include <Utils/FileTools.h>
#include <Utils/TextTools.h>

using namespace bpp;

void help()
{
  *ApplicationTools::message << "__________________________________________________________________________" << endl;
  *ApplicationTools::message << "bpppars parameter1_name=parameter1_value parameter2_name=parameter2_value"  << endl;
  *ApplicationTools::message << "      ... param=option_file" << endl;
  *ApplicationTools::message << endl;
  *ApplicationTools::message << "    Refer to the Bio++ Program Suite Manual for list of available options." << endl;
  *ApplicationTools::message << "__________________________________________________________________________" << endl;
}

int main(int args, char ** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*             Bio++ Parsimony Methods, version 0.1.0             *" << endl;
  cout << "* Author: J. Dutheil                        Created     05/05/07 *" << endl;
  cout << "*                                           Last Modif. 08/08/09 *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  if(args == 1)
  {
    help();
    return 0;
  }
  
  try {
 
  BppApplication bpppars(args, argv, "BppPars");
  bpppars.startTimer();

	Alphabet * alphabet = SequenceApplicationTools::getAlphabet(bpppars.getParams(), "", false);

	VectorSiteContainer * allSites = SequenceApplicationTools::getSiteContainer(alphabet, bpppars.getParams());
	
	VectorSiteContainer * sites = SequenceApplicationTools::getSitesToAnalyse(* allSites, bpppars.getParams());
	delete allSites;

  ApplicationTools::displayResult("Number of sequences", TextTools::toString(sites->getNumberOfSequences()));
  ApplicationTools::displayResult("Number of sites", TextTools::toString(sites->getNumberOfSites()));
	
  // Get the initial tree
  Tree* tree = 0;
  string initTree = ApplicationTools::getStringParameter("init.tree", bpppars.getParams(), "user", "", false, false);
  ApplicationTools::displayResult("Input tree", initTree);
  if(initTree == "user")
  {
    tree = PhylogeneticsApplicationTools::getTree(bpppars.getParams());
    ApplicationTools::displayResult("Number of leaves", TextTools::toString(tree->getNumberOfLeaves()));
  }
  else if(initTree == "random")
  {
    vector<string> names = sites->getSequencesNames();
    tree = TreeTemplateTools::getRandomTree(names);
    tree->setBranchLengths(1.);
  }
  else throw Exception("Unknown init tree method.");
	
  ApplicationTools::displayTask("Initializing parsimony");
  DRTreeParsimonyScore * tp = new DRTreeParsimonyScore(*tree, *sites, false);
  delete tree;
  ApplicationTools::displayTaskDone();
  double score = tp->getScore();
  ApplicationTools::displayResult("Initial parsimony score", TextTools::toString(score, 15));
  bool optTopo = ApplicationTools::getBooleanParameter("optimization.topology", bpppars.getParams(), false);
  ApplicationTools::displayResult("Optimize topology", optTopo ? "yes" : "no");
  if(optTopo)
  {
    tp = OptimizationTools::optimizeTreeNNI(tp, 1);
    score = tp->getScore();
    ApplicationTools::displayResult("Final parsimony score", TextTools::toString(score, 15));
  }
  tree = new TreeTemplate<Node>(*tp->getTree());
  
	PhylogeneticsApplicationTools::writeTree(*tree, bpppars.getParams());
  
  //Bootstrap:
  unsigned int nbBS = ApplicationTools::getParameter<unsigned int>("bootstrap.number", bpppars.getParams(), 0);
  if(nbBS > 0)
  {
    ApplicationTools::displayResult("Number of bootstrap samples", TextTools::toString(nbBS));
    const Tree * initTree = tree;
    if(!optTopo)
    {
      tp = OptimizationTools::optimizeTreeNNI(tp, 1);
      initTree = tp->getTree();
    }

    
    string bsTreesPath = ApplicationTools::getAFilePath("bootstrap.output.file", bpppars.getParams(), false, false);
    ofstream *out = NULL;
    if(bsTreesPath != "none")
    {
      ApplicationTools::displayResult("Bootstrap trees stored in file", bsTreesPath);
      out = new ofstream(bsTreesPath.c_str(), ios::out);
    }
    Newick newick;

    ApplicationTools::displayTask("Bootstrapping", true);
    vector<Tree *> bsTrees(nbBS);
    for(unsigned int i = 0; i < nbBS; i++)
    {
      ApplicationTools::displayGauge(i, nbBS-1, '=');
      VectorSiteContainer * sample = SiteContainerTools::bootstrapSites(*sites);
      DRTreeParsimonyScore * tp = new DRTreeParsimonyScore(*initTree, *sample, false);
      tp = OptimizationTools::optimizeTreeNNI(tp, 0);
      bsTrees[i] = new TreeTemplate<Node>(*tp->getTree());
      if(out && i==0) newick.write(*bsTrees[i], bsTreesPath, true);
      if(out && i>0) newick.write(*bsTrees[i], bsTreesPath, false);
      delete tp;
      delete sample;
    }
    if(out) out->close();
    if(out) delete out;
    ApplicationTools::displayTaskDone();
    

    ApplicationTools::displayTask("Compute bootstrap values");
    TreeTools::computeBootstrapValues(*tree, bsTrees);
    ApplicationTools::displayTaskDone();
    for(unsigned int i = 0; i < nbBS; i++) delete bsTrees[i];

    //Write resulting tree:
    PhylogeneticsApplicationTools::writeTree(*tree, bpppars.getParams());
  }

	delete sites;
  delete tp;
	delete alphabet;

  bpppars.done();
	}
  catch(exception & e)
  {
		cout << e.what() << endl;
		return 1;
	}

	return 0;
}

