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
#include <Bpp/Text/TextTools.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

// From PhylLib:
#include <Bpp/Phyl/Tree/Tree.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Parsimony/DRTreeParsimonyScore.h>

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

    std::shared_ptr<const bpp::Alphabet> alphabet = SequenceApplicationTools::getAlphabet(bpppars.getParams(), "", false);
  
    bool includeGaps = ApplicationTools::getBooleanParameter("use.gaps", bpppars.getParams(), false, "", false, false);
    ApplicationTools::displayBooleanResult("Use gaps", includeGaps);

    auto allSites = SequenceApplicationTools::getSiteContainer(alphabet, bpppars.getParams());
	
    shared_ptr<VectorSiteContainer> sites = SequenceApplicationTools::getSitesToAnalyse(* allSites, bpppars.getParams(), "", true, !includeGaps, true);

    ApplicationTools::displayResult("Number of sequences", TextTools::toString(sites->getNumberOfSequences()));
    ApplicationTools::displayResult("Number of sites", TextTools::toString(sites->getNumberOfSites()));
	
    if (sites->getNumberOfSequences()==0 || sites->getNumberOfSites()==0)
      throw Exception("Empty data.");

    // Get the initial tree
    std::shared_ptr<Tree> tree = nullptr;
    string initTreeOpt = ApplicationTools::getStringParameter("init.tree", bpppars.getParams(), "user", "", false, false);
    ApplicationTools::displayResult("Input tree", initTreeOpt);
    if (initTreeOpt == "user")
    {
      tree = PhylogeneticsApplicationTools::getTree(bpppars.getParams());
      ApplicationTools::displayResult("Number of leaves", TextTools::toString(tree->getNumberOfLeaves()));
    }
    else if (initTreeOpt == "random")
    {
      vector<string> names = sites->getSequenceNames();
      tree = TreeTemplateTools::getRandomTree(names, false);
      tree->setBranchLengths(1.);
    }
    else throw Exception("Unknown init tree method.");
	
    ApplicationTools::displayTask("Initializing parsimony");
    auto treen = make_shared<TreeTemplate<Node>>(*tree);
    auto tp = std::make_unique<DRTreeParsimonyScore>(treen, sites, false, includeGaps);

    ApplicationTools::displayTaskDone();
    double score = tp->getScore();
    ApplicationTools::displayResult("Initial parsimony score", TextTools::toString(score, 15));
    tree = make_shared<TreeTemplate<Node>>(tp->tree());
  
    PhylogeneticsApplicationTools::writeTree(*tree, bpppars.getParams());
  
    bpppars.done();
  }
  catch (exception & e)
  {
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}

