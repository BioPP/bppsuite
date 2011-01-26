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
#include <limits>

using namespace std;

#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/KeyvalTools.h>

// From SeqLib:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

// From PhylLib:
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/Likelihood.all>
#include <Bpp/Phyl/PatternTools.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Model.all>
#include <Bpp/Phyl/Io/Newick.h>

using namespace bpp;

/******************************************************************************/

void help()
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
  (*ApplicationTools::message << "bppml parameter1_name=parameter1_value parameter2_name=parameter2_value").endLine();
  (*ApplicationTools::message << "      ... param=option_file").endLine();
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "  Refer to the Bio++ Program Suite Manual for a list of available options.").endLine();
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
}

int main(int args, char** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*       Bio++ Maximum Likelihood Computation, version 1.4.0      *" << endl;
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

  try
  {
    BppApplication bppml(args, argv, "BppML");
    bppml.startTimer();

    Alphabet* alphabet = SequenceApplicationTools::getAlphabet(bppml.getParams(), "", false);

    VectorSiteContainer* allSites = SequenceApplicationTools::getSiteContainer(alphabet, bppml.getParams());

    VectorSiteContainer* sites = SequenceApplicationTools::getSitesToAnalyse(*allSites, bppml.getParams(), "", true, false);
    delete allSites;

    ApplicationTools::displayResult("Number of sequences", TextTools::toString(sites->getNumberOfSequences()));
    ApplicationTools::displayResult("Number of sites", TextTools::toString(sites->getNumberOfSites()));

    // Get the initial tree
    Tree* tree = 0;
    string initTreeOpt = ApplicationTools::getStringParameter("init.tree", bppml.getParams(), "user", "", false, false);
    ApplicationTools::displayResult("Input tree", initTreeOpt);
    if (initTreeOpt == "user")
    {
      tree = PhylogeneticsApplicationTools::getTree(bppml.getParams());
      ApplicationTools::displayResult("Number of leaves", TextTools::toString(tree->getNumberOfLeaves()));
    }
    else if (initTreeOpt == "random")
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
      return 0;
    }

    // Setting branch lengths?
    string initBrLenMethod = ApplicationTools::getStringParameter("init.brlen.method", bppml.getParams(), "Input", "", true, false);
    string cmdName;
    map<string, string> cmdArgs;
    KeyvalTools::parseProcedure(initBrLenMethod, cmdName, cmdArgs);
    if (cmdName == "Input")
    {
      // Do nothing!
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
      tree->scaleTree(h / nh);
    }
    else throw Exception("Method '" + initBrLenMethod + "' unknown for computing branch lengths.");
    ApplicationTools::displayResult("Branch lengths", cmdName);

    string treeWIdPath = ApplicationTools::getAFilePath("output.tree_ids.file", bppml.getParams(), false, false);
    if (treeWIdPath != "none")
    {
      TreeTemplate<Node> ttree(*tree);
      vector<Node*> nodes = ttree.getNodes();
      for (unsigned int i = 0; i < nodes.size(); i++)
      {
        if (nodes[i]->isLeaf())
          nodes[i]->setName(TextTools::toString(nodes[i]->getId()) + "_" + nodes[i]->getName());
        else
          nodes[i]->setBranchProperty("NodeId", BppString(TextTools::toString(nodes[i]->getId())));
      }
      Newick treeWriter;
      treeWriter.enableExtendedBootstrapProperty("NodeId");
      ApplicationTools::displayResult("Writing tagged tree to", treeWIdPath);
      treeWriter.write(ttree, treeWIdPath);
      delete tree;
      cout << "BppML's done." << endl;
      exit(0);
    }

    DiscreteRatesAcrossSitesTreeLikelihood* tl;
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
      if (model->getNumberOfStates() >= 2 * model->getAlphabet()->getSize())
      {
        // Markov-modulated Markov model!
        rDist = new ConstantDistribution(1., true);
      }
      else
      {
        rDist = PhylogeneticsApplicationTools::getRateDistribution(bppml.getParams());
      }
      if (dynamic_cast<MixedSubstitutionModel*>(model) == 0)
        tl = new RHomogeneousClockTreeLikelihood(*tree, *sites, model, rDist, true, true);
      else
        throw Exception("Molecular clock with Mixed model not supported yet, sorry :(");
    }
    else if (optimizeClock == "no")
    {
      if (optimizeTopo || nbBS > 0)
      {
        if (nhOpt != "no")
          throw Exception("Topology estimation with NH model not supported yet, sorry :(");
        model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, sites, bppml.getParams());
        if (model->getName() != "RE08") SiteContainerTools::changeGapsToUnknownCharacters(*sites);
        if (model->getNumberOfStates() >= 2 * model->getAlphabet()->getSize())
        {
          // Markov-modulated Markov model!
          rDist = new ConstantDistribution(1., true);
        }
        else
        {
          rDist = PhylogeneticsApplicationTools::getRateDistribution(bppml.getParams());
        }
        if (dynamic_cast<MixedSubstitutionModel*>(model) == 0)
          tl = new NNIHomogeneousTreeLikelihood(*tree, *sites, model, rDist, true, true);
        else
          throw Exception("Topology estimation with Mixed model not supported yet, sorry :(");
      }
      else if (nhOpt == "no")
      {
        model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, sites, bppml.getParams());
        if (model->getName() != "RE08") SiteContainerTools::changeGapsToUnknownCharacters(*sites);
        if (model->getNumberOfStates() >= 2 * model->getAlphabet()->getSize())
        {
          // Markov-modulated Markov model!
          rDist = new ConstantDistribution(1., true);
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
            if (dynamic_cast<MixedSubstitutionModel*>(model) == 0)
              tl = new RHomogeneousTreeLikelihood(*tree, *sites, model, rDist, true, true, false);
            else
              tl = new RHomogeneousMixedTreeLikelihood(*tree, *sites, model, rDist, true, true, false);

          else if (compression == "recursive")
            if (dynamic_cast<MixedSubstitutionModel*>(model) == 0)
              tl = new RHomogeneousTreeLikelihood(*tree, *sites, model, rDist, true, true, true);
            else
              tl = new RHomogeneousMixedTreeLikelihood(*tree, *sites, model, rDist, true, true, true);

          else throw Exception("Unknown likelihood data compression method: " + compression);
        }
        else if (recursion == "double")
        {
          if (dynamic_cast<MixedSubstitutionModel*>(model) == 0)
            tl = new DRHomogeneousTreeLikelihood(*tree, *sites, model, rDist, true);
          else
            tl = new DRHomogeneousMixedTreeLikelihood(*tree, *sites, model, rDist, true);
        }
        else throw Exception("Unknown recursion option: " + recursion);
      }
      else if (nhOpt == "one_per_branch")
      {
        model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, sites, bppml.getParams());
        if (model->getName() != "RE08") SiteContainerTools::changeGapsToUnknownCharacters(*sites);
        if (model->getNumberOfStates() >= 2 * model->getAlphabet()->getSize())
        {
          // Markov-modulated Markov model!
          rDist = new ConstantDistribution(1., true);
        }
        else
        {
          rDist = PhylogeneticsApplicationTools::getRateDistribution(bppml.getParams());
        }
        vector<double> rateFreqs;
        if (model->getNumberOfStates() != alphabet->getSize())
        {
          // Markov-Modulated Markov Model...
          unsigned int n = (unsigned int)(model->getNumberOfStates() / alphabet->getSize());
          rateFreqs = vector<double>(n, 1. / (double)n); // Equal rates assumed for now, may be changed later (actually, in the most general case,
                                                         // we should assume a rate distribution for the root also!!!
        }

        bool stationarity = ApplicationTools::getBooleanParameter("nonhomogeneous.stationarity", bppml.getParams(), false, "", false, false);
        FrequenciesSet* rootFreqs = 0;
        if (!stationarity)
        {
          rootFreqs = PhylogeneticsApplicationTools::getRootFrequenciesSet(alphabet, sites, bppml.getParams(), rateFreqs);
          stationarity = !rootFreqs;
        }
        ApplicationTools::displayBooleanResult("Stationarity assumed", stationarity);
   
        vector<string> globalParameters = ApplicationTools::getVectorParameter<string>("nonhomogeneous_one_per_branch.shared_parameters", bppml.getParams(), ',', "");
        for (unsigned int i = 0; i < globalParameters.size(); i++)
          ApplicationTools::displayResult("Global parameter", globalParameters[i]);
        modelSet = SubstitutionModelSetTools::createNonHomogeneousModelSet(model, rootFreqs, tree, globalParameters);
        model = 0;

        string recursion = ApplicationTools::getStringParameter("likelihood.recursion", bppml.getParams(), "simple", "", true, false);
        ApplicationTools::displayResult("Likelihood recursion", recursion);
        if (recursion == "simple")
        {
          if (modelSet->hasMixedSubstitutionModel())
            tl = new RNonHomogeneousMixedTreeLikelihood(*tree, *sites, modelSet, rDist, true, true);
          else
            tl = new RNonHomogeneousTreeLikelihood(*tree, *sites, modelSet, rDist, true, true);
        }
        else if (recursion == "double")
        {
          if (modelSet->hasMixedSubstitutionModel())
            throw Exception("Double recursion with non homogeneous mixed models is not implemented yet.");
              //            tl = new DRNonHomogeneousMixedTreeLikelihood(*tree, *sites, modelSet, rDist, true);
          else
            tl = new DRNonHomogeneousTreeLikelihood(*tree, *sites, modelSet, rDist, true);
        }
        else throw Exception("Unknown recursion option: " + recursion);
      }
      else if (nhOpt == "general")
      {
        modelSet = PhylogeneticsApplicationTools::getSubstitutionModelSet(alphabet, sites, bppml.getParams());
        if (modelSet->getModel(0)->getName() != "RE08") SiteContainerTools::changeGapsToUnknownCharacters(*sites);
        if (modelSet->getNumberOfStates() >= 2 * modelSet->getAlphabet()->getSize())
        {
          // Markov-modulated Markov model!
          rDist = new ConstantDistribution(1., true);
        }
        else
        {
          rDist = PhylogeneticsApplicationTools::getRateDistribution(bppml.getParams());
        }

        string recursion = ApplicationTools::getStringParameter("likelihood.recursion", bppml.getParams(), "simple", "", true, false);
        ApplicationTools::displayResult("Likelihood recursion", recursion);
        if (recursion == "simple")
        {
          if (modelSet->hasMixedSubstitutionModel())
            tl = new RNonHomogeneousMixedTreeLikelihood(*tree, *sites, modelSet, rDist, true, true);
          else
            tl = new RNonHomogeneousTreeLikelihood(*tree, *sites, modelSet, rDist, true, true);
        }
        else if (recursion == "double")
          if (modelSet->hasMixedSubstitutionModel())
            throw Exception("Double recursion with non homogeneous mixed models is not implemented yet.");
              //            tl = new DRNonHomogeneousMixedTreeLikelihood(*tree, *sites, modelSet, rDist, true);
          else
            tl = new DRNonHomogeneousTreeLikelihood(*tree, *sites, modelSet, rDist, true);
        else throw Exception("Unknown recursion option: " + recursion);
      }
      else throw Exception("Unknown option for nonhomogeneous: " + nhOpt);
    }
    else throw Exception("Unknown option for optimization.clock: " + optimizeClock);

    tl->initialize();

    delete tree;

    string paramNameFile = ApplicationTools::getAFilePath("output.parameter_names.file", bppml.getParams(), false, false);
    if (paramNameFile != "none") {
      ofstream pnfile(paramNameFile.c_str(), ios::out);
      ParameterList pl = tl->getParameters();
      for (unsigned int i = 0; i < pl.size(); ++i) {
        pnfile << pl[i].getName() << endl;
      }
      pnfile.close();
      cout << "BppML's done." << endl;
      exit(0);
    }

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
      CodonAlphabet *pca=dynamic_cast<CodonAlphabet*>(alphabet);
      if (pca) {
        bool f = false;
        unsigned int  s;
        for (unsigned int i = 0; i < sites->getNumberOfSites(); i++) {
          if (isinf(tl->getLogLikelihoodForASite(i))) {
            const Site& site=sites->getSite(i);
            s = site.size();
            for (unsigned int j = 0; j < s; j++) {
              if (pca->isStop(site.getValue(j))) {
                (*ApplicationTools::error << "Stop Codon at site " << site.getPosition() << " in sequence " << sites->getSequence(j).getName()).endLine();
                f = true;
              }
            }
          }
        }
        if (f)
          exit(-1);
      }
      bool removeSaturated = ApplicationTools::getBooleanParameter("input.sequence.remove_saturated_sites", bppml.getParams(), false, "", true, false);
      if (!removeSaturated) {
        ApplicationTools::displayError("!!! Looking at each site:");
        for (unsigned int i = 0; i < sites->getNumberOfSites(); i++) {
          (*ApplicationTools::error << "Site " << sites->getSite(i).getPosition() << "\tlog likelihood = " << tl->getLogLikelihoodForASite(i)).endLine();
        }
        ApplicationTools::displayError("!!! 0 values (inf in log) may be due to computer overflow, particularily if datasets are big (>~500 sequences).");
        ApplicationTools::displayError("!!! You may want to try input.sequence.remove_saturated_sites = yes to ignore positions with likelihood 0.");
        exit(1);
      } else {
        ApplicationTools::displayBooleanResult("Saturated site removal enabled", true);
        for (unsigned int i = sites->getNumberOfSites(); i > 0; --i) {
          if (isinf(tl->getLogLikelihoodForASite(i - 1))) {
            ApplicationTools::displayResult("Ignore saturated site", sites->getSite(i - 1).getPosition());
            sites->deleteSite(i - 1);
          }
        }
        ApplicationTools::displayResult("Number of sites retained", sites->getNumberOfSites());
        tl->setData(*sites);
        tl->initialize();
        logL = tl->getValue();
        if (isinf(logL)) {
          ApplicationTools::displayError("This should not happen. Exiting now.");
          exit(1);
        }
        ApplicationTools::displayResult("Initial log likelihood", TextTools::toString(-logL, 15));
      }
    }

    if (optimizeClock == "global")
    {
      PhylogeneticsApplicationTools::optimizeParameters(
        dynamic_cast<DiscreteRatesAcrossSitesClockTreeLikelihood*>(tl), tl->getParameters(), bppml.getParams());
    }
    else
    {
      tl = dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood*>(
        PhylogeneticsApplicationTools::optimizeParameters(tl, tl->getParameters(), bppml.getParams()));
    }

    tree = new TreeTemplate<Node>(tl->getTree());
    PhylogeneticsApplicationTools::writeTree(*tree, bppml.getParams());

    // Write parameters to screen:
    ApplicationTools::displayResult("Log likelihood", TextTools::toString(-tl->getValue(), 15));
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

    // Checking convergence:
    PhylogeneticsApplicationTools::checkEstimatedParameters(tl->getParameters());

    // Write parameters to file:
    string parametersFile = ApplicationTools::getAFilePath("output.estimates", bppml.getParams(), false, false);
    ApplicationTools::displayResult("Output estimates to file", parametersFile);
    if (parametersFile != "none")
    {
      StlOutputStream out(auto_ptr<ostream>(new ofstream(parametersFile.c_str(), ios::out)));
      out << "# Log likelihood = ";
      out.setPrecision(20) << (-tl->getValue());
      out.endLine();
      out.endLine();
      out << "# Substitution model parameters:";
      out.endLine();
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
      out.endLine();
      (out << "# Rate distribution parameters:").endLine();
      rDist->matchParametersValues(tl->getParameters());
      PhylogeneticsApplicationTools::printParameters(rDist, out);
    }

    // Getting posterior rate class distribution:
    DiscreteDistribution* prDist = RASTools::getPosteriorRateDistribution(*tl);
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
      DataTable* infos = new DataTable(colNames);

      for (unsigned int i = 0; i < sites->getNumberOfSites(); i++)
      {
        double lnL = tl->getLogLikelihoodForASite(i);
        const Site* currentSite = &sites->getSite(i);
        int currentSitePosition = currentSite->getPosition();
        string isCompl = "NA";
        string isConst = "NA";
        try { isCompl = (SiteTools::isComplete(*currentSite) ? "1" : "0"); }
        catch(EmptySiteException& ex) {}
        try { isConst = (SiteTools::isConstant(*currentSite) ? "1" : "0"); }
        catch(EmptySiteException& ex) {}
        row[0] = (string("[" + TextTools::toString(currentSitePosition) + "]"));
        row[1] = isCompl;
        row[2] = isConst;
        row[3] = TextTools::toString(lnL);
        row[4] = TextTools::toString(classes[i]);
        row[5] = TextTools::toString(rates[i]);
        infos->addRow(row);
      }

      DataTable::write(*infos, out, "\t");

      delete infos;
    }


    // Bootstrap:
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
        tl = dynamic_cast<NNIHomogeneousTreeLikelihood*>(
            PhylogeneticsApplicationTools::optimizeParameters(tl, tl->getParameters(), bppml.getParams(), "", true, false));
        initTree = &tl->getTree();
      }

      string bsTreesPath = ApplicationTools::getAFilePath("bootstrap.output.file", bppml.getParams(), false, false);
      ofstream* out = 0;
      if (bsTreesPath != "none")
      {
        ApplicationTools::displayResult("Bootstrap trees stored in file", bsTreesPath);
        out = new ofstream(bsTreesPath.c_str(), ios::out);
      }
      Newick newick;
      ParameterList paramsToIgnore = tl->getSubstitutionModelParameters();
      paramsToIgnore.addParameters(tl->getRateDistributionParameters());

      ApplicationTools::displayTask("Bootstrapping", true);
      vector<Tree*> bsTrees(nbBS);
      for (unsigned int i = 0; i < nbBS; i++)
      {
        ApplicationTools::displayGauge(i, nbBS - 1, '=');
        VectorSiteContainer* sample = SiteContainerTools::bootstrapSites(*sites);
        if (!approx)
        {
          model->setFreqFromData(*sample);
        }

        if (dynamic_cast<MixedSubstitutionModel*>(model) != NULL)
          throw Exception("Bootstrap estimation with Mixed model not supported yet, sorry :(");

        NNIHomogeneousTreeLikelihood* tlRep = new NNIHomogeneousTreeLikelihood(*initTree, *sample, model, rDist, true, false);
        tlRep->initialize();
        ParameterList parametersRep = tlRep->getParameters();
        if (approx)
        {
          parametersRep.deleteParameters(paramsToIgnore.getParameterNames());
        }
        tlRep = dynamic_cast<NNIHomogeneousTreeLikelihood*>(
          PhylogeneticsApplicationTools::optimizeParameters(tlRep, parametersRep, bppml.getParams(), "", true, false));
        bsTrees[i] = new TreeTemplate<Node>(tlRep->getTree());
        if (out && i == 0) newick.write(*bsTrees[i], bsTreesPath, true);
        if (out && i >  0) newick.write(*bsTrees[i], bsTreesPath, false);
        delete tlRep;
        delete sample;
      }
      if (out) out->close();
      if (out) delete out;
      ApplicationTools::displayTaskDone();


      ApplicationTools::displayTask("Compute bootstrap values");
      TreeTools::computeBootstrapValues(*tree, bsTrees);
      ApplicationTools::displayTaskDone();
      for (unsigned int i = 0; i < nbBS; i++)
      {
        delete bsTrees[i];
      }

      // Write resulting tree:
      PhylogeneticsApplicationTools::writeTree(*tree, bppml.getParams());
    }


    delete alphabet;
    delete sites;
    if (model) delete model;
    if (modelSet) delete modelSet;
    delete rDist;
    delete tl;
    delete tree;
    bppml.done();
  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}

