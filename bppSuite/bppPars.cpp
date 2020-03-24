//
// File: bppPars.cpp
// Created by: Julien Dutheil
// Created on: May Sat 05 15:09 2007
// From file bppML.cpp
//

/*
Copyright or © or Copr. Bio++ Development Team

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

// From bpp-core:
#include <Bpp/Version.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

// From bpp-phyl:
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/PatternTools.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Io/Newick.h>

using namespace bpp;

void help()
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
  (*ApplicationTools::message << "bpppars parameter1_name=parameter1_value parameter2_name=parameter2_value").endLine();
  (*ApplicationTools::message << "      ... param=option_file").endLine();
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "    Refer to the Bio++ Program Suite Manual for list of available options.").endLine();
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
}

int main(int args, char ** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*             Bio++ Parsimony Methods, version " << BPP_VERSION << "             *" << endl;
  cout << "* Author: J. Dutheil                        Created     05/05/07 *" << endl;
  cout << "*                                           Last Modif. " << BPP_REL_DATE << " *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  if (args == 1)
  {
    help();
    return 0;
  }
  
  try {
 
  BppApplication bpppars(args, argv, "BppPars");
  bpppars.startTimer();

	Alphabet* alphabet = SequenceApplicationTools::getAlphabet(bpppars.getParams(), "", false);
  
  bool includeGaps = ApplicationTools::getBooleanParameter("use.gaps", bpppars.getParams(), false, "", false, false);
  ApplicationTools::displayBooleanResult("Use gaps", includeGaps);

	VectorSiteContainer* allSites = SequenceApplicationTools::getSiteContainer(alphabet, bpppars.getParams());
	
	VectorSiteContainer* sites = SequenceApplicationTools::getSitesToAnalyse(* allSites, bpppars.getParams(), "", true, !includeGaps, true);
	delete allSites;
  
  ApplicationTools::displayResult("Number of sequences", TextTools::toString(sites->getNumberOfSequences()));
  ApplicationTools::displayResult("Number of sites", TextTools::toString(sites->getNumberOfSites()));
	
  // Get the initial tree
  Tree* tree = 0;
  string initTreeOpt = ApplicationTools::getStringParameter("init.tree", bpppars.getParams(), "user", "", false, false);
  ApplicationTools::displayResult("Input tree", initTreeOpt);
  if (initTreeOpt == "user")
  {
    tree = PhylogeneticsApplicationTools::getTree(bpppars.getParams());
    ApplicationTools::displayResult("Number of leaves", TextTools::toString(tree->getNumberOfLeaves()));
  }
  else if (initTreeOpt == "random")
  {
    vector<string> names = sites->getSequencesNames();
    tree = TreeTemplateTools::getRandomTree(names, false);
    tree->setBranchLengths(1.);
  }
  else throw Exception("Unknown init tree method.");
	
  ApplicationTools::displayTask("Initializing parsimony");
  DRTreeParsimonyScore* tp = new DRTreeParsimonyScore(*tree, *sites, false, includeGaps);
  delete tree;
  ApplicationTools::displayTaskDone();
  double score = tp->getScore();
  ApplicationTools::displayResult("Initial parsimony score", TextTools::toString(score, 15));
  bool optTopo = ApplicationTools::getBooleanParameter("optimization.topology", bpppars.getParams(), false);
  ApplicationTools::displayResult("Optimize topology", optTopo ? "yes" : "no");
  if (optTopo)
  {
    tp = OptimizationTools::optimizeTreeNNI(tp, 1);
    score = tp->getScore();
    ApplicationTools::displayResult("Final parsimony score", TextTools::toString(score, 15));
  }
  tree = new TreeTemplate<Node>(tp->getTree());
  
	PhylogeneticsApplicationTools::writeTree(*tree, bpppars.getParams());
  
  //Bootstrap:
  unsigned int nbBS = ApplicationTools::getParameter<unsigned int>("bootstrap.number", bpppars.getParams(), 0);
  if (nbBS > 0)
  {
    ApplicationTools::displayResult("Number of bootstrap samples", TextTools::toString(nbBS));
    const Tree* initTree = tree;
    if (!optTopo)
    {
      tp = OptimizationTools::optimizeTreeNNI(tp, 1);
      initTree = &tp->getTree();
    }

    
    string bsTreesPath = ApplicationTools::getAFilePath("bootstrap.output.file", bpppars.getParams(), false, false);
    ofstream *out = 0;
    if (bsTreesPath != "none")
    {
      ApplicationTools::displayResult("Bootstrap trees stored in file", bsTreesPath);
      out = new ofstream(bsTreesPath.c_str(), ios::out);
    }
    Newick newick;

    ApplicationTools::displayTask("Bootstrapping", true);
    vector<Tree*> bsTrees(nbBS);
    for (unsigned int i = 0; i < nbBS; i++)
    {
      ApplicationTools::displayGauge(i, nbBS - 1, '=');
      VectorSiteContainer* sample = SiteContainerTools::bootstrapSites(*sites);
      DRTreeParsimonyScore* tpRep = new DRTreeParsimonyScore(*initTree, *sample, false);
      tpRep = OptimizationTools::optimizeTreeNNI(tpRep, 0);
      bsTrees[i] = new TreeTemplate<Node>(tpRep->getTree());
      if (out && i==0) newick.writeTree(*bsTrees[i], bsTreesPath, true);
      if (out && i>0) newick.writeTree(*bsTrees[i], bsTreesPath, false);
      delete tpRep;
      delete sample;
    }
    if(out) out->close();
    if(out) delete out;
    ApplicationTools::displayTaskDone();
    

    ApplicationTools::displayTask("Compute bootstrap values", true);
    TreeTools::computeBootstrapValues(*tree, bsTrees);
    ApplicationTools::displayTaskDone();
    for (unsigned int i = 0; i < nbBS; i++)
      delete bsTrees[i];

    //Write resulting tree:
    PhylogeneticsApplicationTools::writeTree(*tree, bpppars.getParams());
  }

	delete sites;
  delete tp;
	delete alphabet;

  bpppars.done();
	}
  catch (exception & e)
  {
		cout << e.what() << endl;
		return 1;
	}

	return 0;
}

