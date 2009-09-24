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

using namespace std;

// From SeqLib:
#include <Seq/Alphabet.h>
#include <Seq/VectorSiteContainer.h>
#include <Seq/SiteTools.h>
#include <Seq/SequenceApplicationTools.h>

// From PhylLib:
#include <Phyl/Tree.h>
#include <Phyl/DiscreteRatesAcrossSitesTreeLikelihood.h>
#include <Phyl/RHomogeneousTreeLikelihood.h>
#include <Phyl/DRHomogeneousTreeLikelihood.h>
#include <Phyl/NNIHomogeneousTreeLikelihood.h>
#include <Phyl/RHomogeneousClockTreeLikelihood.h>
#include <Phyl/RNonHomogeneousTreeLikelihood.h>
#include <Phyl/DRNonHomogeneousTreeLikelihood.h>
#include <Phyl/PatternTools.h>
#include <Phyl/PhylogeneticsApplicationTools.h>
#include <Phyl/MarginalAncestralStateReconstruction.h>
#include <Phyl/OptimizationTools.h>
#include <Phyl/RASTools.h>
#include <Phyl/Newick.h>
#include <Phyl/MarkovModulatedSubstitutionModel.h>
#include <Phyl/SubstitutionModelSet.h>
#include <Phyl/SubstitutionModelSetTools.h>
#include <Phyl/UserProteinSubstitutionModel.h>
#include <Phyl/UserProteinSubstitutionModelF.h>

// From NumCalc:
#include <NumCalc/DiscreteDistribution.h>
#include <NumCalc/ConstantDistribution.h>
#include <NumCalc/DataTable.h>
#include <NumCalc/MatrixTools.h>
#include <NumCalc/VectorTools.h>
#include <NumCalc/AutoParameter.h>

// From Utils:
#include <Utils/BppApplication.h>
#include <Utils/ApplicationTools.h>
#include <Utils/FileTools.h>
#include <Utils/TextTools.h>
#include <Utils/KeyvalTools.h>

using namespace bpp;

/******************************************************************************/

void help()
{
  *ApplicationTools::message << "__________________________________________________________________________" << endl;
  *ApplicationTools::message << "bppml parameter1_name=parameter1_value parameter2_name=parameter2_value"    << endl;
  *ApplicationTools::message << "      ... param=option_file" << endl;
  *ApplicationTools::message << endl;
  *ApplicationTools::message << "  Refer to the Bio++ Program Suite Manual for a list of available options." << endl;
  *ApplicationTools::message << "__________________________________________________________________________" << endl;
}

int main(int args, char ** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*       Bio++ Maximum Likelihood Computation, version 1.3.0      *" << endl;
  cout << "*                                                                *" << endl; 
  cout << "* Authors: J. Dutheil                       Last Modif. 08/08/09 *" << endl;
  cout << "*          B. Boussau                                            *" << endl;
  cout << "*          L. Gueguen                                            *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  if (args == 1)
  {
    help();
    return 0;
  }
  
  try {

  BppApplication bppml(args, argv, "BppML");
  bppml.startTimer();

  Alphabet* alphabet = SequenceApplicationTools::getAlphabet(bppml.getParams(), "", false);

  VectorSiteContainer* allSites = SequenceApplicationTools::getSiteContainer(alphabet, bppml.getParams());
  
  VectorSiteContainer* sites = SequenceApplicationTools::getSitesToAnalyse(* allSites, bppml.getParams(), "", true, false);
  delete allSites;

  ApplicationTools::displayResult("Number of sequences", TextTools::toString(sites->getNumberOfSequences()));
  ApplicationTools::displayResult("Number of sites", TextTools::toString(sites->getNumberOfSites()));
  
  // Get the initial tree
  Tree* tree = 0;
  string initTree = ApplicationTools::getStringParameter("init.tree", bppml.getParams(), "user", "", false, false);
  ApplicationTools::displayResult("Input tree", initTree);
  if (initTree == "user")
  {
    tree = PhylogeneticsApplicationTools::getTree(bppml.getParams());
    ApplicationTools::displayResult("Number of leaves", TextTools::toString(tree->getNumberOfLeaves()));
  }
  else if (initTree == "random")
  {
    vector<string> names = sites->getSequencesNames();
    tree = TreeTemplateTools::getRandomTree(names);
    tree->setBranchLengths(1.);
  }
  else throw Exception("Unknown init tree method.");
  
  // Try to write the current tree to file. This will be overwritten by the optimized tree,
  // but allow to check file existence before running optimization!
  PhylogeneticsApplicationTools::writeTree(*tree, bppml.getParams());
  
  bool computeLikelihood = ApplicationTools::getBooleanParameter("compute.likelihood", bppml.getParams(), true, "", false, false);
  if (!computeLikelihood)
  {
    delete alphabet;
    delete sites;
    delete tree;
    cout << "BppML's done. Bye." << endl;
    return(0);
  }
  
  // Setting branch lengths?
  string initBrLenMethod = ApplicationTools::getStringParameter("init.brlen.method", bppml.getParams(), "Input", "", true, false);
  string cmdName;
  map<string, string> cmdArgs;
  KeyvalTools::parseProcedure(initBrLenMethod, cmdName, cmdArgs);
  if (cmdName == "Input")
  {
    //Do nothing!
  }
  else if (cmdName == "Equal")
  {
    double value = ApplicationTools::getDoubleParameter("value", cmdArgs, 0.1, "", true, false);
    if (value <= 0)
      throw Exception("Value for branch length must be superior to 0");
    ApplicationTools::displayResult("Branch lengths set to", value);
    tree->setBranchLengths(value);
  }
  else if (cmdName == "Clock")
  {
    TreeTools::convertToClockTree(*tree, tree->getRootId(), true);
  }
  else if (cmdName == "Grafen")
  {
    string grafenHeight = ApplicationTools::getStringParameter("height", cmdArgs, "Input", "", true, false);
    double h;
    if (grafenHeight == "input")
    {
      h = TreeTools::getHeight(*tree, tree->getRootId());
    }
    else
    {
      h = TextTools::toDouble(grafenHeight);
      if (h <= 0) throw Exception("Height must be positive in Grafen's method.");
    }
    ApplicationTools::displayResult("Total height", TextTools::toString(h));
    
    double rho = ApplicationTools::getDoubleParameter("rho", cmdArgs, 1., "", true, false);
    ApplicationTools::displayResult("Grafen's rho", rho);
    TreeTools::computeBranchLengthsGrafen(*tree, rho);
    double nh = TreeTools::getHeight(*tree, tree->getRootId());
    tree->scaleTree(h/nh);
  }
  else throw Exception("Method '" + initBrLenMethod + "' unknown for computing branch lengths.");
  ApplicationTools::displayResult("Branch lengths", cmdName);

  string treeWIdPath = ApplicationTools::getAFilePath("output.tree_ids.file", bppml.getParams(), false, false);
  if (treeWIdPath != "none")
  {
    TreeTemplate<Node> ttree(*tree);
    vector<Node *> nodes = ttree.getNodes();
    for (unsigned int i = 0; i < nodes.size(); i++)
    {
      if (nodes[i]->isLeaf())
        nodes[i]->setName(TextTools::toString(nodes[i]->getId()) + "_" + nodes[i]->getName());
      else
        nodes[i]->setBranchProperty("NodeId", String(TextTools::toString(nodes[i]->getId())));
    }
    Newick treeWriter;
    treeWriter.enableExtendedBootstrapProperty("NodeId");
    ApplicationTools::displayResult("Writing tagged tree to", treeWIdPath);
    treeWriter.write(ttree, treeWIdPath);
    delete tree;
    cout << "BppML's done." << endl;
    exit(0);
  }

  DiscreteRatesAcrossSitesTreeLikelihood *tl;
  string optimizeClock = ApplicationTools::getStringParameter("optimization.clock", bppml.getParams(), "no", "", true, false);
  ApplicationTools::displayResult("Clock", optimizeClock);
  string nhOpt = ApplicationTools::getStringParameter("nonhomogeneous", bppml.getParams(), "no", "", true, false);
  ApplicationTools::displayResult("Heterogeneous model", nhOpt);

  bool optimizeTopo = ApplicationTools::getBooleanParameter("optimization.topology", bppml.getParams(), false, "", true, false);
  unsigned int nbBS = ApplicationTools::getParameter<unsigned int>("bootstrap.number", bppml.getParams(), 0, "", true, false);
  
  SubstitutionModel*    model    = 0;
  SubstitutionModelSet* modelSet = 0;
  DiscreteDistribution* rDist    = 0;

  if (optimizeClock == "global")
  {
    model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, sites, bppml.getParams());
    if (model->getName() != "RE08") SiteContainerTools::changeGapsToUnknownCharacters(*sites);
    if (model->getNumberOfStates() >= 2*model->getAlphabet()->getSize())
    {
      //Markov-modulated Markov model!
      rDist = new ConstantDistribution(1.);
    }
    else
    {
      rDist = PhylogeneticsApplicationTools::getRateDistribution(bppml.getParams());
    }
    tl = new RHomogeneousClockTreeLikelihood(*tree, *sites, model, rDist, true, true);
  }
  else if (optimizeClock == "no")
  {
    if (optimizeTopo || nbBS > 0)
    {
      if (nhOpt != "no")
        throw Exception("Topology estimation with NH model not supported yet, sorry :(");
      model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, sites, bppml.getParams());
      if (model->getName() != "RE08") SiteContainerTools::changeGapsToUnknownCharacters(*sites);
      if (model->getNumberOfStates() >= 2*model->getAlphabet()->getSize())
      {
        //Markov-modulated Markov model!
        rDist = new ConstantDistribution(1.);
      }
      else
      {
        rDist = PhylogeneticsApplicationTools::getRateDistribution(bppml.getParams());
      }
      tl = new NNIHomogeneousTreeLikelihood(*tree, *sites, model, rDist, true, true);
    }
    else if (nhOpt == "no")
    { 
      model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, sites, bppml.getParams());
      if (model->getName() != "RE08") SiteContainerTools::changeGapsToUnknownCharacters(*sites);
      if (model->getNumberOfStates() >= 2*model->getAlphabet()->getSize())
      {
        //Markov-modulated Markov model!
        rDist = new ConstantDistribution(1.);
      }
      else
      {
        rDist = PhylogeneticsApplicationTools::getRateDistribution(bppml.getParams());
      }
      string recursion = ApplicationTools::getStringParameter("likelihood.recursion", bppml.getParams(), "simple", "", true, false);
      ApplicationTools::displayResult("Likelihood recursion", recursion);
      if (recursion == "simple")
      {
        string compression = ApplicationTools::getStringParameter("likelihood.recursion_simple.compression", bppml.getParams(), "recursive", "", true, false);
        ApplicationTools::displayResult("Likelihood data compression", compression);
        if (compression == "simple")
          tl = new RHomogeneousTreeLikelihood(*tree, *sites, model, rDist, true, true, false);
        else if (compression == "recursive")
          tl = new RHomogeneousTreeLikelihood(*tree, *sites, model, rDist, true, true, true);
        else throw Exception("Unknown likelihood data compression method: " + compression);
      }
      else if (recursion == "double")
      {
        tl = new DRHomogeneousTreeLikelihood(*tree, *sites, model, rDist, true);
      }
      else throw Exception("Unknown recursion option: " + recursion);
    }
    else if (nhOpt == "one_per_branch")
    {
      model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, sites, bppml.getParams());
      if (model->getName() != "RE08") SiteContainerTools::changeGapsToUnknownCharacters(*sites);
      if (model->getNumberOfStates() >= 2*model->getAlphabet()->getSize())
      {
        //Markov-modulated Markov model!
        rDist = new ConstantDistribution(1.);
      }
      else
      {
        rDist = PhylogeneticsApplicationTools::getRateDistribution(bppml.getParams());
      }
      vector<double> rateFreqs;
      if (model->getNumberOfStates() != alphabet->getSize())
      {
        //Markov-Modulated Markov Model...
        unsigned int n =(unsigned int)(model->getNumberOfStates() / alphabet->getSize());
        rateFreqs = vector<double>(n, 1./(double)n); // Equal rates assumed for now, may be changed later (actually, in the most general case,
                                                     // we should assume a rate distribution for the root also!!!  
      }
      FrequenciesSet* rootFreqs = PhylogeneticsApplicationTools::getRootFrequenciesSet(alphabet, sites, bppml.getParams(), rateFreqs);
      vector<string> globalParameters = ApplicationTools::getVectorParameter<string>("nonhomogeneous_one_per_branch.shared_parameters", bppml.getParams(), ',', "");
      modelSet = SubstitutionModelSetTools::createNonHomogeneousModelSet(model, rootFreqs, tree, globalParameters); 
      model = 0;
      
      string recursion = ApplicationTools::getStringParameter("likelihood.recursion", bppml.getParams(), "simple", "", true, false);
      ApplicationTools::displayResult("Likelihood recursion", recursion);
      if (recursion == "simple")
        tl = new RNonHomogeneousTreeLikelihood(*tree, *sites, modelSet, rDist, true, true);
      else if (recursion == "double")
        tl = new DRNonHomogeneousTreeLikelihood(*tree, *sites, modelSet, rDist, true);
      else throw Exception("Unknown recursion option: " + recursion);
    }
    else if (nhOpt == "general")
    {
      modelSet = PhylogeneticsApplicationTools::getSubstitutionModelSet(alphabet, sites, bppml.getParams());
      if (modelSet->getModel(0)->getName() != "RE08") SiteContainerTools::changeGapsToUnknownCharacters(*sites);
      if (modelSet->getNumberOfStates() >= 2*modelSet->getAlphabet()->getSize())
      {
        //Markov-modulated Markov model!
        rDist = new ConstantDistribution(1.);
      }
      else
      {
        rDist = PhylogeneticsApplicationTools::getRateDistribution(bppml.getParams());
      }

      string recursion = ApplicationTools::getStringParameter("likelihood.recursion", bppml.getParams(), "simple", "", true, false);
      ApplicationTools::displayResult("Likelihood recursion", recursion);
      if (recursion == "simple")
        tl = new RNonHomogeneousTreeLikelihood(*tree, *sites, modelSet, rDist, true, true);
      else if (recursion == "double")
        tl = new DRNonHomogeneousTreeLikelihood(*tree, *sites, modelSet, rDist, true);
      else throw Exception("Unknown recursion option: " + recursion);
    }
    else throw Exception("Unknown option for nonhomogeneous: " + nhOpt);
  }
  else throw Exception("Unknown option for optimization.clock: " + optimizeClock);
  tl->initialize();
 
  delete tree;
    
  double logL = tl->getValue();
  if (isinf(logL))
  {
    // This may be due to null branch lengths, leading to null likelihood!
    ApplicationTools::displayWarning("!!! Warning!!! Initial likelihood is zero.");
    ApplicationTools::displayWarning("!!! This may be due to branch length == 0.");
    ApplicationTools::displayWarning("!!! All null branch lengths will be set to 0.000001.");
    ParameterList pl = tl->getBranchLengthsParameters();
    for (unsigned int i = 0; i < pl.size(); i++)
    {
      if (pl[i].getValue() < 0.000001) pl[i].setValue(0.000001);
    }
    tl->matchParametersValues(pl);
    logL = tl->getValue();
  }
  ApplicationTools::displayResult("Initial log likelihood", TextTools::toString(-logL, 15));
  if (isinf(logL))
  {
    ApplicationTools::displayError("!!! Unexpected initial likelihood == 0.");
    ApplicationTools::displayError("!!! Looking at each site:");
    for (unsigned int i = 0; i < sites->getNumberOfSites(); i++)
    {
      *ApplicationTools::error << "Site " << sites->getSite(i).getPosition() << "\tlog likelihood = " << tl->getLogLikelihoodForASite(i) << endl;
    }
    ApplicationTools::displayError("!!! 0 values (inf in log) may be due to computer overflow, particularily if datasets are big (>~500 sequences).");
    exit(-1);
  }

  if (optimizeClock == "global")
  {
    PhylogeneticsApplicationTools::optimizeParameters(
        dynamic_cast<DiscreteRatesAcrossSitesClockTreeLikelihood *>(tl), tl->getParameters(), bppml.getParams());
  }
  else
  {
    tl = dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood *>(
        PhylogeneticsApplicationTools::optimizeParameters(tl, tl->getParameters(), bppml.getParams()));
  }
  
  tree = new TreeTemplate<Node>(tl->getTree());
  PhylogeneticsApplicationTools::writeTree(*tree, bppml.getParams());
  
  // Write parameters to screen:
  ApplicationTools::displayResult("Log likelihood", TextTools::toString(- tl->getValue(), 15));
  ParameterList parameters = tl->getSubstitutionModelParameters();
  for (unsigned int i = 0; i < parameters.size(); i++)
  {
    ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
  }
  parameters = tl->getRateDistributionParameters();
  for (unsigned int i = 0; i < parameters.size(); i++)
  {
    ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
  }
  // Write parameters to file:
  string parametersFile = ApplicationTools::getAFilePath("output.estimates", bppml.getParams(), false, false);
  ApplicationTools::displayResult("Output estimates to file", parametersFile);
  if (parametersFile != "none")
  {
    ofstream out(parametersFile.c_str(), ios::out);
    out << "# Log likelihood = " << (- tl->getValue()) << endl;
    out << endl;
    out << "# Substitution model parameters:" << endl;
    if (modelSet)
    {
      modelSet->matchParametersValues(tl->getParameters());
      PhylogeneticsApplicationTools::printParameters(modelSet, out);
    }
    else
    {
      model->matchParametersValues(tl->getParameters());
      PhylogeneticsApplicationTools::printParameters(model, out);
    }
    out << endl;
    out << "# Rate distribution parameters:" << endl;
    rDist->matchParametersValues(tl->getParameters());
    PhylogeneticsApplicationTools::printParameters(rDist, out);
    out.close();
  }

  // Getting posterior rate class distribution:
  DiscreteDistribution * prDist = RASTools::getPosteriorRateDistribution(* tl);
  ApplicationTools::displayMessage("\nPosterior rate distribution for dataset:\n");
  if (ApplicationTools::message) prDist->print(*ApplicationTools::message);
  ApplicationTools::displayMessage("\n");
  delete prDist;
 
  // Write infos to file:
  string infosFile = ApplicationTools::getAFilePath("output.infos", bppml.getParams(), false, false);
  if (infosFile != "none")
  {
    ApplicationTools::displayResult("Alignment information logfile", infosFile);
    ofstream out(infosFile.c_str(), ios::out);
    
    // Get the rate class with maximum posterior probability:
    vector<unsigned int> classes = tl->getRateClassWithMaxPostProbOfEachSite();
    // Get the posterior rate, i.e. rate averaged over all posterior probabilities:
    Vdouble rates = tl->getPosteriorRateOfEachSite();

    vector<string> colNames;
    colNames.push_back("Sites");
    colNames.push_back("is.complete");
    colNames.push_back("is.constant");
    colNames.push_back("lnL");
    colNames.push_back("rc");
    colNames.push_back("pr");
    vector<string> row(6);
    DataTable * infos = new DataTable(colNames);
    
    for (unsigned int i = 0; i < sites->getNumberOfSites(); i++)
    {
      double lnL = tl->getLogLikelihoodForASite(i);
      const Site* currentSite = &sites->getSite(i);
      int currentSitePosition = currentSite->getPosition();
      int isCompl = (SiteTools::isComplete(* currentSite) ? 1 : 0);
      int isConst = (SiteTools::isConstant(* currentSite) ? 1 : 0);
      row[0] = (string("[" + TextTools::toString(currentSitePosition) + "]"));
      row[1] = TextTools::toString(isCompl);
      row[2] = TextTools::toString(isConst);
      row[3] = TextTools::toString(lnL);
      row[4] = TextTools::toString(classes[i]);
      row[5] = TextTools::toString(rates[i]);
      infos->addRow(row);
    }

    DataTable::write(*infos, out, "\t");

    delete infos;
  }
  
  
  
  //Bootstrap:
  if (nbBS > 0 && optimizeClock != "no")
  {
    ApplicationTools::displayError("Bootstrap is not supported with clock trees.");
  }
  if (nbBS > 0 && optimizeClock == "no")
  {
    ApplicationTools::displayResult("Number of bootstrap samples", TextTools::toString(nbBS));
    bool approx = ApplicationTools::getBooleanParameter("bootstrap.approximate", bppml.getParams(), true);
    ApplicationTools::displayResult("Use approximate bootstrap", TextTools::toString(approx ? "yes" : "no"));
    bool bootstrapVerbose = ApplicationTools::getBooleanParameter("bootstrap.verbose", bppml.getParams(), false, "", true, false);
    
    const Tree* initTree = tree;
    if (!bootstrapVerbose) bppml.getParam("optimization.verbose") = "0";
    bppml.getParam("optimization.profiler") = "none";
    bppml.getParam("optimization.messageHandler") = "none";
    if (!optimizeTopo)
    {
      bppml.getParam("optimization.topology") = "yes";
      tl = dynamic_cast<NNIHomogeneousTreeLikelihood *>(
          PhylogeneticsApplicationTools::optimizeParameters(tl, tl->getParameters(), bppml.getParams(), "", true, false));
      initTree = &tl->getTree();
    }
    
    string bsTreesPath = ApplicationTools::getAFilePath("bootstrap.output.file", bppml.getParams(), false, false);
    ofstream *out = 0;
    if (bsTreesPath != "none")
    {
      ApplicationTools::displayResult("Bootstrap trees stored in file", bsTreesPath);
      out = new ofstream(bsTreesPath.c_str(), ios::out);
    }
    Newick newick;
    ParameterList paramsToIgnore = tl->getSubstitutionModelParameters();
    paramsToIgnore.addParameters(tl->getRateDistributionParameters());

    ApplicationTools::displayTask("Bootstrapping", true);
    vector<Tree *> bsTrees(nbBS);
    for (unsigned int i = 0; i < nbBS; i++)
    {
      ApplicationTools::displayGauge(i, nbBS-1, '=');
      VectorSiteContainer * sample = SiteContainerTools::bootstrapSites(*sites);
      if (!approx)
      {
        model->setFreqFromData(*sample);
      }
      NNIHomogeneousTreeLikelihood *tlrep = new NNIHomogeneousTreeLikelihood(*initTree, *sample, model, rDist, true, false);
      tlrep->initialize();
      ParameterList parameters = tlrep->getParameters();
      if (approx)
      {
        parameters.deleteParameters(paramsToIgnore.getParameterNames());
      }
      tlrep = dynamic_cast<NNIHomogeneousTreeLikelihood *>(
          PhylogeneticsApplicationTools::optimizeParameters(tlrep, parameters, bppml.getParams(), "", true, false));
      bsTrees[i] = new TreeTemplate<Node>(tlrep->getTree());
      if (out && i == 0) newick.write(*bsTrees[i], bsTreesPath, true);
      if (out && i >  0) newick.write(*bsTrees[i], bsTreesPath, false);
      delete tlrep;
      delete sample;
    }
    if (out) out->close();
    if (out) delete out;
    ApplicationTools::displayTaskDone();
    

    ApplicationTools::displayTask("Compute bootstrap values");
    TreeTools::computeBootstrapValues(*tree, bsTrees);
    ApplicationTools::displayTaskDone();
    for (unsigned int i = 0; i < nbBS; i++) delete bsTrees[i];

    //Write resulting tree:
    PhylogeneticsApplicationTools::writeTree(*tree, bppml.getParams());
  }


  delete alphabet;
  delete sites;
  if (model)    delete model;
  if (modelSet) delete modelSet;
  delete rDist;
  delete tl;
  delete tree;
  bppml.done();
 
  }
  catch(exception & e)
  {
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}

