//
// File: bppDist.cpp
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
#include <Phyl/DistanceEstimation.h>
#include <Phyl/IODistanceMatrixFactory.h>
#include <Phyl/PGMA.h>
#include <Phyl/NeighborJoining.h>
#include <Phyl/BioNJ.h>
#include <Phyl/OptimizationTools.h>
#include <Phyl/MarkovModulatedSubstitutionModel.h>
#include <Phyl/Newick.h>

// From NumCalc:
#include <NumCalc/DiscreteDistribution.h>
#include <NumCalc/ConstantDistribution.h>

// From Utils:
#include <Utils/AttributesTools.h>
#include <Utils/ApplicationTools.h>
#include <Utils/FileTools.h>
#include <Utils/TextTools.h>

using namespace bpp;

void help()
{
  *ApplicationTools::message << "__________________________________________________________________________" << endl;
  SequenceApplicationTools::printInputAlignmentHelp();
  PhylogeneticsApplicationTools::printInputTreeHelp();
  PhylogeneticsApplicationTools::printOutputTreeHelp();
}

int main(int args, char ** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*              Bio++ Distance Methods, version 0.3.0             *" << endl;
  cout << "* Author: J. Dutheil                        Created     05/05/07 *" << endl;
  cout << "*                                           Last Modif. 08/05/09 *" << endl;
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

  Alphabet * alphabet = SequenceApplicationTools::getAlphabet(params, "", false);

  VectorSiteContainer * allSites = SequenceApplicationTools::getSiteContainer(alphabet, params);
  
  VectorSiteContainer * sites = SequenceApplicationTools::getSitesToAnalyse(* allSites, params);
  delete allSites;

  ApplicationTools::displayResult("Number of sequences", TextTools::toString(sites->getNumberOfSequences()));
  ApplicationTools::displayResult("Number of sites", TextTools::toString(sites->getNumberOfSites()));
  
  SubstitutionModel * model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, sites, params);
  
	DiscreteDistribution * rDist = NULL;
  if(model->getNumberOfStates() > model->getAlphabet()->getSize())
  {
    //Markov-modulated Markov model!
    rDist = new ConstantDistribution(1.);
  }
  else
  {
	  rDist = PhylogeneticsApplicationTools::getRateDistribution(params);
  }
   
  DistanceEstimation distEstimation(model, rDist, sites, 1, false);
 
  string method = ApplicationTools::getStringParameter("method", params, "nj");
  ApplicationTools::displayResult("Tree reconstruction method", method);
  TreeTemplate<Node> * tree;
  AgglomerativeDistanceMethod * distMethod = NULL;
  if(method == "wpgma")
  {
    PGMA * wpgma = new PGMA(true);
    distMethod = wpgma;
  }
  else if(method == "upgma")
  {
    PGMA * upgma = new PGMA(false);
    distMethod = upgma;
  }
  else if(method == "nj")
  {
    NeighborJoining * nj = new NeighborJoining();
    nj->outputPositiveLengths(true);
    distMethod = nj;
  }
  else if(method == "bionj")
  {
    BioNJ * bionj = new BioNJ();
    bionj->outputPositiveLengths(true);
    distMethod = bionj;
  }
  else throw Exception("Unknown tree reconstruction method.");
  
  string type = ApplicationTools::getStringParameter("optimization.method", params, "init");
  ApplicationTools::displayResult("Model parameters estimation method", type);
  if(type == "init") type = OptimizationTools::DISTANCEMETHOD_INIT;
  else if(type == "pairwise") type = OptimizationTools::DISTANCEMETHOD_PAIRWISE;
  else if(type == "iterations") type = OptimizationTools::DISTANCEMETHOD_ITERATIONS;
  else throw Exception("Unknown parameter estimation procedure '" + type + "'.");
  
	unsigned int optVerbose = ApplicationTools::getParameter<unsigned int>("optimization.verbose", params, 2);
	
	string mhPath = ApplicationTools::getAFilePath("optimization.message_handler", params, false, false);
	ostream * messenger = 
		(mhPath == "none") ? NULL :
			(mhPath == "std") ? &cout :
				new ofstream(mhPath.c_str(), ios::out);
	ApplicationTools::displayResult("Message handler", mhPath);

	string prPath = ApplicationTools::getAFilePath("optimization.profiler", params, false, false);
	ostream * profiler = 
		(prPath == "none") ? NULL :
			(prPath == "std") ? &cout :
				new ofstream(prPath.c_str(), ios::out);
	if(profiler != NULL) (*profiler) << setprecision(20);
	ApplicationTools::displayResult("Profiler", prPath);

	// Should I ignore some parameters?
  ParameterList allParameters = model->getParameters();
  allParameters.addParameters(rDist->getParameters());
	ParameterList parametersToIgnore;
  string paramListDesc = ApplicationTools::getStringParameter("optimization.ignore_parameter", params, "", "", true, false);
	bool ignoreBrLen = false;
  StringTokenizer st(paramListDesc, ",");
	while(st.hasMoreToken())
  {
		try
    {
      string param = st.nextToken();
      if(param == "BrLen")
        ignoreBrLen = true;
      else
      {
        if (allParameters.hasParameter(param))
        {
          Parameter* p = &allParameters.getParameter(param);
          parametersToIgnore.addParameter(*p);
        }
        else ApplicationTools::displayWarning("Parameter '" + param + "' not found."); 
      }
		} 
    catch(ParameterNotFoundException pnfe)
    {
			ApplicationTools::displayError("Parameter '" + pnfe.getParameter() + "' not found, and so can't be ignored!");
		}
	}
	
	unsigned int nbEvalMax = ApplicationTools::getParameter<unsigned int>("optimization.max_number_f_eval", params, 1000000);
	ApplicationTools::displayResult("Max # ML evaluations", TextTools::toString(nbEvalMax));
	
	double tolerance = ApplicationTools::getDoubleParameter("optimization.tolerance", params, .000001);
	ApplicationTools::displayResult("Tolerance", TextTools::toString(tolerance));
	
  //Here it is:
  ofstream out("warnings", ios::out);
  ApplicationTools::warning = &out;
  tree = OptimizationTools::buildDistanceTree(distEstimation, *distMethod, parametersToIgnore, !ignoreBrLen, false, type, tolerance, nbEvalMax, profiler, messenger, optVerbose);
  ApplicationTools::warning = &cout;

  string matrixPath = ApplicationTools::getAFilePath("output.matrix.file", params, false, false, "", false);
  if(matrixPath != "none")
  {
    ApplicationTools::displayResult("Output matrix file", matrixPath);
    ODistanceMatrix * odm = IODistanceMatrixFactory().createWriter(IODistanceMatrixFactory::PHYLIP_FORMAT);
    odm->write(*distEstimation.getMatrix(), matrixPath, true);
    delete odm;
  }

  PhylogeneticsApplicationTools::writeTree(*tree, params);
  
  //Output some parameters:
  if(type == OptimizationTools::DISTANCEMETHOD_ITERATIONS)
  {
    // Write parameters to screen:
    ParameterList parameters = model->getParameters();
    for(unsigned int i = 0; i < parameters.size(); i++)
    {
		  ApplicationTools::displayResult(parameters[i]->getName(), TextTools::toString(parameters[i]->getValue()));
    }
    parameters = rDist->getParameters();
    for(unsigned int i = 0; i < parameters.size(); i++)
    {
		  ApplicationTools::displayResult(parameters[i]->getName(), TextTools::toString(parameters[i]->getValue()));
    }
    // Write parameters to file:
	  string parametersFile = ApplicationTools::getAFilePath("output.estimates", params, false, false);
    if(parametersFile != "none")
    {
		  ofstream out(parametersFile.c_str(), ios::out);
      parameters = model->getParameters();
      for(unsigned int i = 0; i < parameters.size(); i++)
      {
        out << parameters[i]->getName() << " = " << parameters[i]->getValue() << endl;
      }
      parameters = rDist->getParameters();
      for(unsigned int i = 0; i < parameters.size(); i++)
      {
        out << parameters[i]->getName() << " = " << parameters[i]->getValue() << endl;
      }
      out.close();
    }
  }
 
  //Bootstrap:
  unsigned int nbBS = ApplicationTools::getParameter<unsigned int>("bootstrap.number", params, 0);
  if(nbBS > 0)
  {
    ApplicationTools::displayResult("Number of bootstrap samples", TextTools::toString(nbBS));
    bool approx = ApplicationTools::getBooleanParameter("bootstrap.approximate", params, true);
    ApplicationTools::displayResult("Use approximate bootstrap", TextTools::toString(approx ? "yes" : "no"));
    if(approx)
    {
      type = OptimizationTools::DISTANCEMETHOD_INIT;
      parametersToIgnore = allParameters;
      ignoreBrLen = true;
    }
    bool bootstrapVerbose = ApplicationTools::getBooleanParameter("bootstrap.verbose", params, false, "", true, false);
 
    string bsTreesPath = ApplicationTools::getAFilePath("bootstrap.output.file", params, false, false);
    ofstream *out = NULL;
    if(bsTreesPath != "none")
    {
      ApplicationTools::displayResult("Bootstrap trees stored in file", bsTreesPath);
      out = new ofstream(bsTreesPath.c_str(), ios::out);
    }
    Newick newick;
    
    vector<Tree *> bsTrees(nbBS);
    ApplicationTools::displayTask("Bootstrapping", true);
    for(unsigned int i = 0; i < nbBS; i++)
    {
      ApplicationTools::displayGauge(i, nbBS-1, '=');
      VectorSiteContainer * sample = SiteContainerTools::bootstrapSites(*sites);
      if(approx) model->setFreqFromData(*sample);
      distEstimation.setData(sample);
      bsTrees[i] = OptimizationTools::buildDistanceTree(
          distEstimation,
          *distMethod,
          parametersToIgnore,
          ignoreBrLen,
          false,
          type,
          tolerance,
          nbEvalMax,
          NULL,
          NULL,
          (bootstrapVerbose ? 1 : 0)
        );
      if(out && i == 0) newick.write(*bsTrees[i], bsTreesPath, true);
      if(out && i >  0) newick.write(*bsTrees[i], bsTreesPath, false);
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
    PhylogeneticsApplicationTools::writeTree(*tree, params);
  }
    
  delete alphabet;
  delete sites;
  delete model;
  delete rDist;
  delete distMethod;
  delete tree;
  cout << "BppDist's done. Bye." << endl;
  ApplicationTools::displayTime("Total execution time:");

  }
  catch(exception & e)
  {
    cout << e.what() << endl;
    exit(-1);
  }

  return (0);
}

