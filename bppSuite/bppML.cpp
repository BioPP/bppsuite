//
// File: bppML.cpp
// Created by: Julien Dutheil
// Created on: Dec Sat 03 16:41 2005
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

// From SeqLib:
#include <Seq/Alphabet.h>
#include <Seq/VectorSiteContainer.h>
#include <Seq/SiteTools.h>
#include <Seq/SequenceApplicationTools.h>

// From PhylLib:
#include <Phyl/Tree.h>
#include <Phyl/DiscreteRatesAcrossSitesTreeLikelihood.h>
#include <Phyl/HomogeneousTreeLikelihood.h>
#include <Phyl/DRHomogeneousTreeLikelihood.h>
#include <Phyl/NNIHomogeneousTreeLikelihood.h>
#include <Phyl/ClockTreeLikelihood.h>
#include <Phyl/PatternTools.h>
#include <Phyl/PhylogeneticsApplicationTools.h>
#include <Phyl/MarginalAncestralStateReconstruction.h>
#include <Phyl/OptimizationTools.h>
#include <Phyl/RASTools.h>
#include <Phyl/Newick.h>

// From NumCalc:
#include <NumCalc/DiscreteDistribution.h>
#include <NumCalc/ConstantDistribution.h>
#include <NumCalc/DataTable.h>
#include <NumCalc/MatrixTools.h>
#include <NumCalc/VectorTools.h>
#include <NumCalc/AutoParameter.h>

// From Utils:
#include <Utils/AttributesTools.h>
#include <Utils/ApplicationTools.h>
#include <Utils/FileTools.h>
#include <Utils/TextTools.h>

using namespace std;

/******************************************************************************/

void help()
{
  *ApplicationTools::message << "__________________________________________________________________________" << endl;
  SequenceApplicationTools::printInputAlignmentHelp();
  PhylogeneticsApplicationTools::printInputTreeHelp();
  PhylogeneticsApplicationTools::printSubstitutionModelHelp();
  PhylogeneticsApplicationTools::printRateDistributionHelp();
  PhylogeneticsApplicationTools::printCovarionModelHelp();
  PhylogeneticsApplicationTools::printOptimizationHelp(true, false);
  PhylogeneticsApplicationTools::printOutputTreeHelp();
  *ApplicationTools::message << "output.infos                  | file where to write site infos" << endl;
  *ApplicationTools::message << "output.estimates              | file where to write estimated parameter values" << endl;
  *ApplicationTools::message << "______________________________|___________________________________________" << endl;
}

int main(int args, char ** argv)
{
	cout << "******************************************************************" << endl;
	cout << "*       Bio++ Maximum Likelihood Computation, version 1.1.0      *" << endl;
	cout << "* Author: J. Dutheil                        Last Modif. 19/03/07 *" << endl;
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

	VectorSiteContainer * allSites = SequenceApplicationTools::getSiteContainer(alphabet, params);
	
	VectorSiteContainer * sites = SequenceApplicationTools::getSitesToAnalyse(* allSites, params);
	delete allSites;

  ApplicationTools::displayResult("Number of sequences", TextTools::toString(sites->getNumberOfSequences()));
  ApplicationTools::displayResult("Number of sites", TextTools::toString(sites->getNumberOfSites()));
	
	SubstitutionModel * model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, sites, params);
	
	DiscreteDistribution * rDist = PhylogeneticsApplicationTools::getRateDistribution(params);

  string covarion = ApplicationTools::getStringParameter("covarion", params, "none", "", false, false);
  SubstitutionModel * modelCov = NULL;
  DiscreteDistribution * rDistCov = NULL;
  if(covarion != "none")
  {
    modelCov = model; 
    rDistCov = rDist;
    rDist = new ConstantDistribution(1.);
    model = PhylogeneticsApplicationTools::getCovarionProcess(modelCov, rDistCov, params, "", false, true);
  }
 
  // Get the initial tree
  TreeTemplate<Node> * tree = NULL;
  string initTree = ApplicationTools::getStringParameter("init.tree", params, "user", "", false, false);
  ApplicationTools::displayResult("Input tree", initTree);
  if(initTree == "user")
  {
    tree = PhylogeneticsApplicationTools::getTree(params);
    ApplicationTools::displayResult("Number of leaves", TextTools::toString(tree->getNumberOfLeaves()));
  }
  else if(initTree == "random")
  {
    vector<string> names = sites->getSequencesNames();
    tree = TreeTemplateTools::getRandomTree(names);
    tree->setBranchLengths(1.);
  }
  else throw Exception("Unknown init tree method.");
	
  // Try to write the current tree to file. This will be overwritten by the optimized tree,
  // but allow to check file existence before running optimization!
	PhylogeneticsApplicationTools::writeTree(*tree, params);
	
	bool computeLikelihood = ApplicationTools::getBooleanParameter("compute.likelihood", params, true, "", false, false);
  if(!computeLikelihood)
  {
    delete alphabet;
	  delete sites;
	  delete tree;
    cout << "BppML's done. Bye." << endl;
    return(0);
  }
  
  // Setting branch lengths?
  string initBrLenMethod = ApplicationTools::getStringParameter("init.brlen.method", params, "input", "", true, false);
  ApplicationTools::displayResult("Branch lengths", TextTools::toString(initBrLenMethod));
  if(initBrLenMethod == "input")
  {
    //Do nothing!
  }
  else if(initBrLenMethod == "equal")
  {
    double value = ApplicationTools::getDoubleParameter("init.brlen.method_equal.value", params, 0.1, "", true, false);
    if(value <= 0) throw Exception("Value for branch length must be superior to 0");
    ApplicationTools::displayResult("Branch lengths set to", TextTools::toString(value));
    tree->setBranchLengths(value);
  }
  else if(initBrLenMethod == "clock")
  {
    TreeTools::convertToClockTree(*tree, tree->getRootId(), true);
  }
  else if(initBrLenMethod == "grafen")
  {
    string grafenHeight = ApplicationTools::getStringParameter("init.brlen.method_grafen.height", params, "input", "", true, false);
    double h;
    if(grafenHeight == "input")
    {
      h = TreeTools::getHeight(*tree, tree->getRootId());
    }
    else
    {
      h = TextTools::toDouble(grafenHeight);
      if(h <= 0) throw Exception("Height must be positive in Grafen's method.");
    }
    ApplicationTools::displayResult("Total height", TextTools::toString(h));
    
    double rho = ApplicationTools::getDoubleParameter("init.brlen.method_grafen.rho", params, 1., "", true, false);
    ApplicationTools::displayResult("Grafen's rho", TextTools::toString(rho));
    TreeTools::computeBranchLengthsGrafen(*tree, rho);
    double nh = TreeTools::getHeight(*tree, tree->getRootId());
    tree->scaleTree(h/nh);
  }
  else throw Exception("Method '" + initBrLenMethod + "' unknown for computing branch lengths.");

	AbstractHomogeneousTreeLikelihood *tl;
  string optimizeClock = ApplicationTools::getStringParameter("optimization.clock", params, "no", "", true, false);
  ApplicationTools::displayResult("Clock", optimizeClock);
  bool optimizeTopo = ApplicationTools::getBooleanParameter("optimization.topology", params, false, "", true, false);
  ApplicationTools::displayTask("Initializing likelihood");
  if(optimizeClock == "global") tl = new ClockTreeLikelihood(*tree, *sites, model, rDist, true, true);
  else if(optimizeClock == "no")
  {
    if(optimizeTopo)
      tl = new NNIHomogeneousTreeLikelihood(*tree, *sites, model, rDist, true, true);
    else
      tl = new DRHomogeneousTreeLikelihood(*tree, *sites, model, rDist, true, true);
  }
  else throw Exception("Unknown option for optimization.clock: " + optimizeClock);
  tl->initialize();
 
  ApplicationTools::displayTaskDone();
  delete tree;
		
  double logL = tl->getValue();
  if(isinf(logL))
  {
    // This may be due to null branch lengths, leading to null likelihood!
    ApplicationTools::displayWarning("!!! Warning!!! Initial likelihood is zero.");
    ApplicationTools::displayWarning("!!! This may be due to branch length == 0.");
    ApplicationTools::displayWarning("!!! All null branch lengths will be set to 0.000001.");
    vector<Node*> nodes = tree->getNodes();
    for(unsigned int i = 0; i < nodes.size(); i++)
    {
      if(nodes[i]->hasDistanceToFather() && nodes[i]->getDistanceToFather() < 0.000001) nodes[i]->setDistanceToFather(0.000001);
    }
    tl->initParameters();
	  logL = tl->f(tl->getParameters());
  }
  ApplicationTools::displayResult("Initial likelihood", TextTools::toString(logL, 15));
  if(isinf(logL))
  {
    ApplicationTools::displayError("!!! Unexpected initial likelihood == 0.");
    ApplicationTools::displayError("!!! Looking at each site:");
    for(unsigned int i = 0; i < sites->getNumberOfSites(); i++)
    {
      *ApplicationTools::error << "Site " << sites->getSite(i)->getPosition() << "\tlog likelihood = " << tl->getLogLikelihoodForASite(i) << endl;
    }
    ApplicationTools::displayError("!!! 0 values (inf in log) may be due to computer overflow, particularily if datasets are big (>~500 sequences).");
    exit(-1);
  }

  if(optimizeClock == "global")
  {
    PhylogeneticsApplicationTools::optimizeParameters(dynamic_cast<ClockTreeLikelihood *>(tl), params);
  }
  else
  {
	  tl = dynamic_cast<AbstractHomogeneousTreeLikelihood *>(
        PhylogeneticsApplicationTools::optimizeParameters(tl, params));
  }
	
	const TreeTemplate<Node> * treeOpt = dynamic_cast<const TreeTemplate<Node> *>(tl->getTree());
	PhylogeneticsApplicationTools::writeTree(* treeOpt, params);
	
  // Write parameters to screen:
	ApplicationTools::displayResult("Log likelihood", TextTools::toString(tl->getValue(), 15));
  ParameterList parameters = tl->getSubstitutionModelParameters();
  for(unsigned int i = 0; i < parameters.size(); i++)
  {
		ApplicationTools::displayResult(parameters[i]->getName(), TextTools::toString(parameters[i]->getValue()));
  }
  parameters = tl->getRateDistributionParameters();
  for(unsigned int i = 0; i < parameters.size(); i++)
  {
		ApplicationTools::displayResult(parameters[i]->getName(), TextTools::toString(parameters[i]->getValue()));
  }
  // Write parameters to file:
	string parametersFile = ApplicationTools::getAFilePath("output.estimates", params, false, false);
  if(parametersFile != "none")
  {
		ofstream out(parametersFile.c_str(), ios::out);
		out << "Log likelihood = " << tl->getValue() << endl;
    parameters = tl->getSubstitutionModelParameters();
    for(unsigned int i = 0; i < parameters.size(); i++)
    {
      out << parameters[i]->getName() << " = " << parameters[i]->getValue() << endl;
    }
    parameters = tl->getRateDistributionParameters();
    for(unsigned int i = 0; i < parameters.size(); i++)
    {
      out << parameters[i]->getName() << " = " << parameters[i]->getValue() << endl;
    }
    out.close();
  }

  // Getting posterior rate class distribution:
	DiscreteDistribution * prDist = RASTools::getPosteriorRateDistribution(* tl);
	ApplicationTools::displayMessage("\nPosterior rate distribution for dataset:\n");
	if(ApplicationTools::message) prDist->print(*ApplicationTools::message);
	ApplicationTools::displayMessage("\n");
  delete prDist;
 
	// Write infos to file:
	string infosFile = ApplicationTools::getAFilePath("output.infos", params, false, false);
	if(infosFile != "none")
  {
		ApplicationTools::displayResult("Alignment information logfile", infosFile);
		ofstream out(infosFile.c_str(), ios::out);
		
		// Get the rate class with maximum posterior probability:
		vector<unsigned int> classes = tl->getRateClassWithMaxPostProbOfEachSite();
		// Get the posterior rate, i.e. rate averaged over all posterior probabilities:
		Vdouble rates = tl->getPosteriorRateOfEachSite();

		vector<string> colNames;
		colNames.push_back("is.complete");
		colNames.push_back("is.constant");
		colNames.push_back("lnL");
		colNames.push_back("rc");
		colNames.push_back("pr");
		vector<string> row(5);
		DataTable * infos = new DataTable(colNames);
		
		for(unsigned int i = 0; i < sites->getNumberOfSites(); i++)
    {
			double lnL = tl->getLogLikelihoodForASite(i);
			const Site * currentSite = sites->getSite(i);
			int currentSitePosition = currentSite->getPosition();
			int isCompl = (SiteTools::isComplete(* currentSite) ? 1 : 0);
			int isConst = (SiteTools::isConstant(* currentSite) ? 1 : 0);
			row[0] = TextTools::toString(isCompl);
			row[1] = TextTools::toString(isConst);
			row[2] = TextTools::toString(lnL);
			row[3] = TextTools::toString(classes[i]);
			row[4] = TextTools::toString(rates[i]);
			infos->addRow(string("Site"+TextTools::toString(currentSitePosition)), row);
		}

		// Not compatible with covarion model...
    //MarginalAncestralStateReconstruction masr(* dynamic_cast<DRHomogeneousTreeLikelihood *>(tl));
		//Sequence * ancestors = masr.getAncestralSequenceForNode(tree->getRootNode());
		//vector<string> ancestorColumn(sites->getNumberOfSites());
		//for(unsigned int i = 0; i < ancestorColumn.size(); i++) {
		//	ancestorColumn[i] = alphabet->intToChar((*ancestors)[i]);
		//}
		//infos->addColumn("ancestral.state", ancestorColumn);
		//delete ancestors;

		DataTable::write(*infos, out, "\t");

		delete infos;
	}
	
  
  
  //Bootstrap:
  unsigned int nbBS = ApplicationTools::getParameter<unsigned int>("bootstrap.number", params, 0);
  if(nbBS > 0 && optimizeClock != "no")
  {
    ApplicationTools::displayError("Bootstrap is not supported with clock trees.");
  }
  if(nbBS > 0 && optimizeClock == "no")
  {
    ApplicationTools::displayResult("Number of bootstrap samples", TextTools::toString(nbBS));
    bool approx = ApplicationTools::getBooleanParameter("bootstrap.approximate", params, true);
    ApplicationTools::displayResult("Use approximate bootstrap", TextTools::toString(approx ? "yes" : "no"));
    
    const Tree * initTree = tree;
    if(!optimizeTopo)
    {
      tl = OptimizationTools::optimizeTreeNNI(dynamic_cast<NNIHomogeneousTreeLikelihood *>(tl), 1);
      initTree = tl->getTree();
    }
    
    string bsTreesPath = ApplicationTools::getAFilePath("bootstrap.output.file", params, false, false);
    ofstream *out = NULL;
    if(bsTreesPath != "none")
    {
      ApplicationTools::displayResult("Bootstrap trees stored in file", bsTreesPath);
      out = new ofstream(bsTreesPath.c_str(), ios::out);
    }
    Newick newick;
    ParameterList paramsToIgnore = tl->getSubstitutionModelParameters();
    paramsToIgnore.addParameters(tl->getRateDistributionParameters());

    ApplicationTools::displayTask("Bootstrapping", true);
    vector<Tree *> bsTrees(nbBS);
    for(unsigned int i = 0; i < nbBS; i++)
    {
      ApplicationTools::displayGauge(i, nbBS-1, '=');
      VectorSiteContainer * sample = SiteContainerTools::bootstrapSites(*sites);
      NNIHomogeneousTreeLikelihood * tl = new NNIHomogeneousTreeLikelihood(*initTree, *sample, model, rDist, true, false);
      tl->initialize();
      if(approx)
      {
        for(unsigned int j = 0; j < paramsToIgnore.size(); j++)
          tl->ignoreParameter(paramsToIgnore[j]->getName());
      }
      else
      {
        model->setFreqFromData(*sample);
      }

	    tl = dynamic_cast<NNIHomogeneousTreeLikelihood *>(
        PhylogeneticsApplicationTools::optimizeParameters(tl, params));
      //tl = OptimizationTools::optimizeTreeNNI(tl, false, 100, 100, 1000000, 1, 0);
      bsTrees[i] = new TreeTemplate<Node>(*tl->getTree());
      if(out && i==0) newick.write(*bsTrees[i], bsTreesPath, true);
      if(out && i>0) newick.write(*bsTrees[i], bsTreesPath, false);
      delete tl;
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
  if(modelCov != NULL) delete modelCov;
  if(rDistCov != NULL) delete rDistCov;
	delete tl;
  cout << "BppML's done. Bye." << endl;
  ApplicationTools::displayTime("Total execution time:");
 
  }
  catch(exception & e)
  {
		cout << e.what() << endl;
		exit(-1);
	}

	return (0);
}

