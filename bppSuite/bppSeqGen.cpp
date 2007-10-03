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
#include <Phyl/MutationProcess.h>
#include <Phyl/NonHomogeneousSequenceSimulator.h>
#include <Phyl/SequenceSimulationTools.h>

// From NumCalc:
#include <NumCalc/DiscreteDistribution.h>
#include <NumCalc/ConstantDistribution.h>
#include <NumCalc/DataTable.h>

// From Utils:
#include <Utils/AttributesTools.h>
#include <Utils/FileTools.h>
#include <Utils/ApplicationTools.h>

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
  cout << "*            Bio++ Sequence Generator, version 0.1               *" << endl;
  cout << "* Author: J. Dutheil                        Last Modif. 24/10/05 *" << endl;
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

  string infosFile = ApplicationTools::getAFilePath("input.infos", params, false, true);
  ApplicationTools::displayResult("Site information", infosFile);

  SubstitutionModel * model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet,NULL, params);
  
  DiscreteDistribution * rDist = NULL;
  HomogeneousSequenceSimulator * seqsim = NULL;
  SiteContainer * sites = NULL;
  if(infosFile != "none")
  {
    ifstream in(infosFile.c_str());
    DataTable * infos = DataTable::read(in, "\t");
    rDist = new ConstantDistribution(1.);
    seqsim = new HomogeneousSequenceSimulator(model, rDist, tree);
    unsigned int nbSites = infos->getNumberOfRows();
    vector<double> rates(nbSites);
    vector<string> ratesStrings = infos->getColumn(string("pr"));
    for(unsigned int i = 0; i < nbSites; i++)
    {
      rates[i] = TextTools::toDouble(ratesStrings[i]);
    }
    sites = SequenceSimulationTools::simulateSites(*seqsim, rates);
  }
  else
  {
    rDist = PhylogeneticsApplicationTools::getRateDistribution(params);
    seqsim = new HomogeneousSequenceSimulator(model, rDist, tree);
    unsigned int nbSites = ApplicationTools::getParameter<unsigned int>("number_of_sites", params, 100);
    sites = seqsim->simulate(nbSites);
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

  return (0);
}

