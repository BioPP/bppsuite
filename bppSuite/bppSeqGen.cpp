//
// File: bppSeqGen.cpp
// Created by: Julien Dutheil
// Created on: Oct Mon 24 18:50 2005
//

/*
Copyright or © or Copr. CNRS

This software is a computer program whose purpose is to simulate sequence
data according to a phylogenetic tree and an evolutionary model.

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
#include <fstream>
#include <iomanip>

// From SeqLib:
#include <Seq/Alphabet.h>
#include <Seq/VectorSiteContainer.h>
#include <Seq/SequenceApplicationTools.h>

// From PhylLib:
#include <Phyl/TreeTemplate.h>
#include <Phyl/PhylogeneticsApplicationTools.h>
#include <Phyl/NonHomogeneousSequenceSimulator.h>
#include <Phyl/SequenceSimulationTools.h>
#include <Phyl/SubstitutionModelSetTools.h>
#include <Phyl/Newick.h>

// From NumCalc:
#include <NumCalc/DiscreteDistribution.h>
#include <NumCalc/ConstantDistribution.h>
#include <NumCalc/DataTable.h>

// From Utils:
#include <Utils/AttributesTools.h>
#include <Utils/FileTools.h>
#include <Utils/ApplicationTools.h>
#include <Utils/Number.h>

using namespace std;

void help()
{
  *ApplicationTools::message << "__________________________________________________________________________" << endl;
  *ApplicationTools::message << "param                         | a parameter file to parse" << endl;
  *ApplicationTools::message << "tree.file                     | tree file path (Newick format)" << endl;
  *ApplicationTools::message << "alphabet                      | the alphabet to use [DNA|RNA|Proteins]" << endl;
  *ApplicationTools::message << "number_of_sites               | number of site to simulate" << endl;
  *ApplicationTools::message << "______________________________|___________________________________________" << endl;
  PhylogeneticsApplicationTools::printSubstitutionModelHelp();
  PhylogeneticsApplicationTools::printRateDistributionHelp();
  SequenceApplicationTools::printOutputSequenceHelp();
}

int main(int args, char ** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*            Bio++ Sequence Generator, version 0.2               *" << endl;
  cout << "* Author: J. Dutheil                                             *" << endl;
  cout << "*         B. Boussau                        Last Modif. 06/10/07 *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;
  
  if(args == 1)
  {
    help();
    exit(0);
  }
  
  try {

  cout << "Parsing options:" << endl;
  
  // Get the parameters from command line:
  map<string, string> cmdParams = AttributesTools::getAttributesMap(
    AttributesTools::getVector(args, argv), "=");

  // Look for a specified file with parameters:
  map<string, string> params;
  if(cmdParams.find("param") != cmdParams.end())
  {
    string file = cmdParams["param"];
    if(!FileTools::fileExists(file))
    {
      cerr << "Parameter file not found." << endl;
      exit(-1);
    }
    else
    {
      params = AttributesTools::getAttributesMapFromFile(file, "=");
      // Actualize attributes with ones passed to command line:
      AttributesTools::actualizeAttributesMap(params, cmdParams);
    }
  }
  else
  {
    params = cmdParams;
  }

  Alphabet * alphabet = SequenceApplicationTools::getAlphabet(params, "", false);

  TreeTemplate<Node> * tree = PhylogeneticsApplicationTools::getTree(params);
  ApplicationTools::displayResult("Number of leaves", TextTools::toString(tree->getNumberOfLeaves()));
  ApplicationTools::displayResult("Number of sons at root", TextTools::toString(tree->getRootNode()->getNumberOfSons()));

  string treeWIdPath = ApplicationTools::getAFilePath("output.tree.path", params, false, false);
  if(treeWIdPath != "none")
  {
    vector<Node *> nodes = tree->getNodes();
    for(unsigned int i = 0; i < nodes.size(); i++)
    {
      if(nodes[i]->isLeaf())
        nodes[i]->setName(TextTools::toString(nodes[i]->getId()) + "_" + nodes[i]->getName());
      else
        nodes[i]->setBranchProperty("NodeId", String(TextTools::toString(nodes[i]->getId())));
    }
    Newick treeWriter;
    treeWriter.enableExtendedBootstrapProperty("NodeId");
    ApplicationTools::displayResult("Writing tagged tree to", treeWIdPath);
    treeWriter.write(*tree, treeWIdPath);
    delete tree;
    cout << "BppSegGen's done." << endl;
    exit(0);
  }

  string infosFile = ApplicationTools::getAFilePath("input.infos", params, false, true);
  ApplicationTools::displayResult("Site information", infosFile);

  string nhOpt = ApplicationTools::getStringParameter("nonhomogeneous", params, "no", "", true, false);
  ApplicationTools::displayResult("Heterogeneous model", nhOpt);

  SubstitutionModelSet * modelSet = NULL;

  //Homogeneous case:
  if(nhOpt == "no")
  {
    SubstitutionModel * model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet,NULL, params);
    modelSet = SubstitutionModelSetTools::createHomogeneousModelSet(model, tree);
  }
  //Galtier-Gouy case:
  else if(nhOpt == "one_per_branch")
  {
    SubstitutionModel * model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet,NULL, params);
    vector<string> globalParameters = ApplicationTools::getVectorParameter<string>("nonhomogeneous_one_per_branch.shared_parameters", params, ',', "");
    modelSet = SubstitutionModelSetTools::createNonHomogeneousModelSet(model, tree, globalParameters); 
  }
  //General case;
  else if(nhOpt == "general")
  {
    modelSet = PhylogeneticsApplicationTools::getSubstitutionModelSet(alphabet,NULL, params);
  }
  else throw Exception("Unknown non-homogeneous option: " + nhOpt);

	DiscreteDistribution * rDist = NULL;
  NonHomogeneousSequenceSimulator * seqsim = NULL;
  SiteContainer * sites = NULL;
  if(infosFile != "none")
  {
    ifstream in(infosFile.c_str());
    DataTable * infos = DataTable::read(in, "\t");
    rDist = new ConstantDistribution(1.);
    seqsim = new NonHomogeneousSequenceSimulator(modelSet, rDist, tree);
    unsigned int nbSites = infos->getNumberOfRows();
    ApplicationTools::displayResult("Number of sites", TextTools::toString(nbSites));
    vector<double> rates(nbSites);
    vector<string> ratesStrings = infos->getColumn(string("pr"));
    for(unsigned int i = 0; i < nbSites; i++)
    {
      rates[i] = TextTools::toDouble(ratesStrings[i]);
    }
    ApplicationTools::displayTask("Perform simulations");
    sites = SequenceSimulationTools::simulateSites(*seqsim, rates);
    ApplicationTools::displayTaskDone();
  }
  else
  {
    if(modelSet->getNumberOfStates() > modelSet->getAlphabet()->getSize())
    {
      //Markov-modulated Markov model!
      rDist = new ConstantDistribution(1.);
    }
    else
    {
	    rDist = PhylogeneticsApplicationTools::getRateDistribution(params);
    }
    seqsim = new NonHomogeneousSequenceSimulator(modelSet, rDist, tree);
    unsigned int nbSites = ApplicationTools::getParameter<unsigned int>("number_of_sites", params, 100);
    ApplicationTools::displayResult("Number of sites", TextTools::toString(nbSites));
    ApplicationTools::displayTask("Perform simulations");
    sites = seqsim->simulate(nbSites);
    ApplicationTools::displayTaskDone();
  }
  
  // Write to file:
  SequenceApplicationTools::writeSequenceFile(*sites, params);

  delete alphabet;
  delete tree;
  delete rDist;
  delete seqsim;

  } catch(exception & e) {
    cout << e.what() << endl;
    exit(-1);
  }

  cout << "BppSegGen's done." << endl;

  return (0);
}

