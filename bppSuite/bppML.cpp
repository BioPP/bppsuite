//
// File: bppML.cpp
// Created by: Julien Dutheil
// Created on: Dec Sat 03 16:41 2005
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
#include <limits>

using namespace std;

#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/Hmm/FullHmmTransitionMatrix.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/KeyvalTools.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

// From bpp-phyl:
#include <Bpp/Phyl/Tree/Tree.h>
#include <Bpp/Phyl/Likelihood.all>
#include <Bpp/Phyl/PatternTools.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Model.all>
#include <Bpp/Phyl/Model/Protein/CoalaCore.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Model/FrequenciesSet/MvaFrequenciesSet.h>
#include <Bpp/Phyl/Io/Newick.h>

#include <Bpp/Phyl/NewLikelihood/SingleRecursiveTreeLikelihoodCalculation.h>
#include <Bpp/Phyl/NewLikelihood/DoubleRecursiveTreeLikelihoodCalculation.h>
#include <Bpp/Phyl/NewLikelihood/SingleProcessPhyloLikelihood.h>
#include <Bpp/Phyl/NewLikelihood/MixturePhyloLikelihood.h>
#include <Bpp/Phyl/NewLikelihood/HmmPhyloLikelihood.h>
#include <Bpp/Phyl/NewLikelihood/SingleDataPhyloLikelihood.h>
#include <Bpp/Phyl/NewLikelihood/SumOfDataPhyloLikelihood.h>
#include <Bpp/Phyl/NewLikelihood/AutoCorrelationPhyloLikelihood.h>
#include <Bpp/Phyl/NewLikelihood/SubstitutionProcessCollection.h>

using namespace bpp;

using namespace newlik;

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
  cout << "*       Bio++ Maximum Likelihood Computation, version 2          *" << endl;
  cout << "*                                                                *" << endl;
  cout << "* Authors: J. Dutheil                       Last Modif. 10/02/14 *" << endl;
  cout << "*          B. Boussau                                            *" << endl;
  cout << "*          L. Guéguen                                            *" << endl;
  cout << "*          M. Groussin                                           *" << endl;
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

    ///// Alphabet
    
    Alphabet* alphabet = SequenceApplicationTools::getAlphabet(bppml.getParams(), "", false);
    auto_ptr<GeneticCode> gCode;
    CodonAlphabet* codonAlphabet = dynamic_cast<CodonAlphabet*>(alphabet);
    if (codonAlphabet) {
      string codeDesc = ApplicationTools::getStringParameter("genetic_code", bppml.getParams(), "Standard", "", true, true);
      ApplicationTools::displayResult("Genetic Code", codeDesc);
      
      gCode.reset(SequenceApplicationTools::getGeneticCode(codonAlphabet->getNucleicAlphabet(), codeDesc));
    }


    ////// Data
    
    VectorSiteContainer* allSites=0;
    vector<SiteContainer* > vSites;
    VectorSiteContainer* sites;

    try {
      allSites = SequenceApplicationTools::getSiteContainer(alphabet, bppml.getParams(), "", false);
      
      sites = SequenceApplicationTools::getSitesToAnalyse(*allSites, bppml.getParams(), "", true, false);
      
      ApplicationTools::displayResult("Number of sequences", TextTools::toString(sites->getNumberOfSequences()));
      ApplicationTools::displayResult("Number of sites", TextTools::toString(sites->getNumberOfSites()));

      vSites.push_back(sites);
    }
    catch(Exception& ex) {
      size_t i=1;

      while (true)
      {
        try {
          allSites = SequenceApplicationTools::getSiteContainer(alphabet, bppml.getParams(), TextTools::toString(i), false);
          sites = SequenceApplicationTools::getSitesToAnalyse(*allSites, bppml.getParams(), TextTools::toString(i), true, false);
      
          ApplicationTools::displayResult("Number of sequences", TextTools::toString(sites->getNumberOfSequences()));
          ApplicationTools::displayResult("Number of sites", TextTools::toString(sites->getNumberOfSites()));

          vSites.push_back(sites);
          i++;
        }
        catch(Exception& ex2) {
          break;
        }
      }
    }
    
    if (vSites.size() == 0)
      throw Exception("Missing data input.sequence.file option");

    /////// Get the initial tree  
    
    vector<Tree*> vTree;
    
    string initTreeOpt = ApplicationTools::getStringParameter("init.tree", bppml.getParams(), "user", "", false, 1);

    ApplicationTools::displayResult("Input tree", initTreeOpt);

    if (initTreeOpt == "user")
    {
      vTree = PhylogeneticsApplicationTools::getTrees(bppml.getParams());

      if (vTree.size()==0){
        ApplicationTools::displayError("!!! Empty file tree.");
        exit(-1);
      }
          
      for (size_t i=0; i<vTree.size(); i++)
        ApplicationTools::displayResult("Number of leaves", TextTools::toString(vTree[i]->getNumberOfLeaves()));
    }
    else if (initTreeOpt == "random")
    {
      vector<string> names = vSites[0]->getSequencesNames();
      Tree* tree = TreeTemplateTools::getRandomTree(names);
      tree->setBranchLengths(1.);
      vTree.push_back(tree);
    }
    else throw Exception("Unknown init tree method.");

    // Try to write the current tree to file. This will be overwritten by the optimized tree,
    // but allow to check file existence before running optimization!

    vector<const Tree*> vcTree;
    for (size_t i=0; i<vTree.size(); i++)
      vcTree.push_back(vTree[i]);
      
    PhylogeneticsApplicationTools::writeTrees(vcTree, bppml.getParams());
    
    bool computeLikelihood = ApplicationTools::getBooleanParameter("compute.likelihood", bppml.getParams(), true, "", false, 1);
    if (!computeLikelihood)
    {
      delete alphabet;
      for (size_t i=0;i<vSites.size(); i++)
        delete vSites[i];
      
      for (size_t i=0; i<vTree.size(); i++)
        delete vTree[i];
      cout << "BppML's done. Bye." << endl;
      return 0;
    }

    ////////////
    // Setting branch lengths?
    string initBrLenMethod = ApplicationTools::getStringParameter("init.brlen.method", bppml.getParams(), "Input", "", true, 1);
    string cmdName;
    map<string, string> cmdArgs;
    KeyvalTools::parseProcedure(initBrLenMethod, cmdName, cmdArgs);
    if (cmdName == "Input")
    {
      // Is the root has to be moved to the midpoint position along the branch that contains it ? If no, do nothing!
      string midPointRootBrLengths = ApplicationTools::getStringParameter("midPointRootBrLengths", cmdArgs, "no", "", true, 2);
      if(midPointRootBrLengths == "yes")
        for (size_t i=0; i< vTree.size(); i++)
          TreeTools::constrainedMidPointRooting(*vTree[i]);
    }
    else if (cmdName == "Equal")
    {
      double value = ApplicationTools::getDoubleParameter("value", cmdArgs, 0.1, "", true, 2);
      if (value <= 0)
        throw Exception("Value for branch length must be superior to 0");
      ApplicationTools::displayResult("Branch lengths set to", value);
      for (size_t i=0; i< vTree.size(); i++)
        vTree[i]->setBranchLengths(value);
    }
    else if (cmdName == "Clock")
    {
      for (size_t i=0; i< vTree.size(); i++)
        TreeTools::convertToClockTree(*vTree[i], vTree[i]->getRootId(), true);
    }
    else if (cmdName == "Grafen")
    {
      string grafenHeight = ApplicationTools::getStringParameter("height", cmdArgs, "input", "", true, 2);
      double h;
      for (size_t i=0; i< vTree.size(); i++)
      {
        Tree* tree=vTree[i];
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

        double rho = ApplicationTools::getDoubleParameter("rho", cmdArgs, 1., "", true, 2);
        ApplicationTools::displayResult("Grafen's rho", rho);
        TreeTools::computeBranchLengthsGrafen(*tree, rho);
        double nh = TreeTools::getHeight(*tree, tree->getRootId());
        tree->scaleTree(h / nh);
      }
    }
    else
      throw Exception("Method '" + initBrLenMethod + "' unknown for computing branch lengths.");
    ApplicationTools::displayResult("Branch lengths", cmdName);
  
    string treeWIdPath = ApplicationTools::getAFilePath("output.tree_ids.file", bppml.getParams(), false, false, "", true, "none", 1);

    if (treeWIdPath != "none")
    {
      Newick treeWriter;
      treeWriter.enableExtendedBootstrapProperty("NodeId");
      ApplicationTools::displayResult("Writing tagged tree to", treeWIdPath);

      for (size_t iv=0; iv< vTree.size(); iv++)
      {
        TreeTemplate<Node> ttree(*vTree[iv]);
        vector<Node*> nodes = ttree.getNodes();
        for (size_t i = 0; i < nodes.size(); i++)
        {
          if (nodes[i]->isLeaf())
            nodes[i]->setName(TextTools::toString(nodes[i]->getId()) + "_" + nodes[i]->getName());
          else
            nodes[i]->setBranchProperty("NodeId", BppString(TextTools::toString(nodes[i]->getId())));
        }
        treeWriter.write(ttree, treeWIdPath, iv==0);
        delete vTree[iv];
      }
      cout << "BppML's done." << endl;
      exit(0);
    }

    /////////////////
    // Computing stuff
    
    DiscreteRatesAcrossSitesTreeLikelihood* tl_old = 0;
    PhyloLikelihood* tl_new = 0;
    SubstitutionProcessCollection* SPC = 0;
    
    bool checkTree    = ApplicationTools::getBooleanParameter("input.tree.check_root", bppml.getParams(), true, "", true, 2);
    bool optimizeTopo = ApplicationTools::getBooleanParameter("optimization.topology", bppml.getParams(), false, "", true, 1);
    unsigned int nbBS = ApplicationTools::getParameter<unsigned int>("bootstrap.number", bppml.getParams(), 0, "", true, 1);
    string collection = ApplicationTools::getStringParameter("collection", bppml.getParams(), "", "", true, 1);
    
    
    SubstitutionModel*    model    = 0;
    SubstitutionModelSet* modelSet = 0;
    DiscreteDistribution* rDist    = 0;

    map<string, string> unparsedparams;
    /// Topology estimation
    
    if (optimizeTopo || nbBS > 0)
    {
      if (collection != "")
        throw Exception("Topology estimation in collections not supported yet, sorry.");
      
      string nhOpt = ApplicationTools::getStringParameter("nonhomogeneous", bppml.getParams(), "no", "", true, false);
      ApplicationTools::displayResult("Heterogeneous model", nhOpt);

      if (nhOpt != "no")
        throw Exception("Topology estimation with NH model not supported yet, sorry :(");
      
      model = PhylogeneticsApplicationTools::getSubstitutionModels(alphabet, gCode.get(), vSites, bppml.getParams(), unparsedparams)[0];
      
      if (model->getName() != "RE08") SiteContainerTools::changeGapsToUnknownCharacters(*vSites[0]);
      if (model->getNumberOfStates() >= 2 * model->getAlphabet()->getSize())
      {
        // Markov-modulated Markov model!
        rDist = new ConstantRateDistribution();
      }
      else
      {
        rDist = PhylogeneticsApplicationTools::getRateDistributions(bppml.getParams())[0];
      }
      if (dynamic_cast<MixedSubstitutionModel*>(model) == 0)
        tl_old = new NNIHomogeneousTreeLikelihood(*vTree[0], *vSites[0], model, rDist, checkTree, true);
      else
        throw Exception("Topology estimation with Mixed model not supported yet, sorry :(");
    }
      
    /// Constant Topology
      
    else
    {
      map<size_t, DiscreteDistribution*> mDist = PhylogeneticsApplicationTools::getRateDistributions(bppml.getParams());

      map<size_t, SubstitutionModel*> mMod = PhylogeneticsApplicationTools::getSubstitutionModels(alphabet, gCode.get(), vSites, bppml.getParams(), unparsedparams);

      map<size_t, FrequenciesSet*> mRootFreq = PhylogeneticsApplicationTools::getRootFrequenciesSets(alphabet, gCode.get(), vSites, bppml.getParams());

      SPC=PhylogeneticsApplicationTools::getSubstitutionProcessCollection(alphabet, gCode.get(), vSites, vTree, mMod, mRootFreq, mDist, bppml.getParams(), unparsedparams);
      
      map<size_t, PhyloLikelihood*> mPhyl=PhylogeneticsApplicationTools::getPhyloLikelihoods(*SPC, vSites, bppml.getParams(), unparsedparams);

      std::vector<SingleDataPhyloLikelihood*> vPhyl;

      map<size_t, PhyloLikelihood*>::iterator it;

      for (it=mPhyl.begin(); it != mPhyl.end(); it++)
        if (dynamic_cast<SingleDataPhyloLikelihood*>(it->second))
          vPhyl.push_back(dynamic_cast<SingleDataPhyloLikelihood*>(it->second));

      if (vPhyl.size()==1)
        tl_new=vPhyl[0];
      else
        tl_new=new SumOfDataPhyloLikelihood(vPhyl);      
    }
    
    tl_new->getParameters().printParameters(cerr);
    
      
    //Listing parameters
    string paramNameFile = ApplicationTools::getAFilePath("output.parameter_names.file", bppml.getParams(), false, false, "", true, "none", 1);
    if (paramNameFile != "none") {
      ApplicationTools::displayResult("List parameters to", paramNameFile);
      ofstream pnfile(paramNameFile.c_str(), ios::out);
        
      ParameterList pl;
      if (tl_old)
        pl=tl_old->getParameters();
      else
        pl=tl_new->getParameters();
        
      for (unsigned int i = 0; i < pl.size(); ++i) {
        pnfile << pl[i].getName() << endl;
      }
      pnfile.close();
      cout << "BppML's done." << endl;
      exit(0);
    }


    // Old optimization
    
    if (tl_old){
      //Check initial likelihood:
      double logL = tl_old->getValue();
      if (isinf(logL))
      {
        // This may be due to null branch lengths, leading to null likelihood!
        ApplicationTools::displayWarning("!!! Warning!!! Initial likelihood is zero.");
        ApplicationTools::displayWarning("!!! This may be due to branch length == 0.");
        ApplicationTools::displayWarning("!!! All null branch lengths will be set to 0.000001.");
        ParameterList pl = tl_old->getBranchLengthsParameters();
        for (unsigned int i = 0; i < pl.size(); i++)
        {
          if (pl[i].getValue() < 0.000001) pl[i].setValue(0.000001);
        }
        tl_old->matchParametersValues(pl);
        logL = tl_old->getValue();
      }
      ApplicationTools::displayResult("Initial log likelihood", TextTools::toString(-logL, 15));
      if (isinf(logL))
      {
        ApplicationTools::displayError("!!! Unexpected initial likelihood == 0.");
        if (codonAlphabet)
        {
          bool f = false;
          size_t s;
          for (size_t i = 0; i < vSites[0]->getNumberOfSites(); i++) {
            if (isinf(tl_old->getLogLikelihoodForASite(i))) {
              const Site& site = vSites[0]->getSite(i);
              s = site.size();
              for (size_t j = 0; j < s; j++) {
                if (gCode->isStop(site.getValue(j))) {
                  (*ApplicationTools::error << "Stop Codon at site " << site.getPosition() << " in sequence " << vSites[0]->getSequence(j).getName()).endLine();
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
          for (unsigned int i = 0; i < vSites[0]->getNumberOfSites(); i++) {
            (*ApplicationTools::error << "Site " << vSites[0]->getSite(i).getPosition() << "\tlog likelihood = " << tl_old->getLogLikelihoodForASite(i)).endLine();
          }
          ApplicationTools::displayError("!!! 0 values (inf in log) may be due to computer overflow, particularily if datasets are big (>~500 sequences).");
          ApplicationTools::displayError("!!! You may want to try input.sequence.remove_saturated_sites = yes to ignore positions with likelihood 0.");
          exit(1);
        } else {
          ApplicationTools::displayBooleanResult("Saturated site removal enabled", true);
          for (size_t i = vSites[0]->getNumberOfSites(); i > 0; --i) {
            if (isinf(tl_old->getLogLikelihoodForASite(i - 1))) {
              ApplicationTools::displayResult("Ignore saturated site", vSites[0]->getSite(i - 1).getPosition());
              vSites[0]->deleteSite(i - 1);
            }
          }
          ApplicationTools::displayResult("Number of sites retained", vSites[0]->getNumberOfSites());
          tl_old->setData(*vSites[0]);
          tl_old->initialize();
          logL = tl_old->getValue();
          if (isinf(logL)) {
            ApplicationTools::displayError("This should not happen. Exiting now.");
            exit(1);
          }
          ApplicationTools::displayResult("Initial log likelihood", TextTools::toString(-logL, 15));
        }
      }
        
      tl_old = dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood*>(
        PhylogeneticsApplicationTools::optimizeParameters(tl_old, tl_old->getParameters(), bppml.getParams()));
        
      Tree* tree = new TreeTemplate<Node>(tl_old->getTree());
      PhylogeneticsApplicationTools::writeTree(*tree, bppml.getParams());
        
      // Write parameters to screen:
      ApplicationTools::displayResult("Log likelihood", TextTools::toString(-tl_old->getValue(), 15));
      ParameterList parameters = tl_old->getSubstitutionModelParameters();
      for (unsigned int i = 0; i < parameters.size(); i++)
      {
        ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
      }
      parameters = tl_old->getRateDistributionParameters();
      for (unsigned int i = 0; i < parameters.size(); i++)
      {
        ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
      }
        
      // Checking convergence:
      PhylogeneticsApplicationTools::checkEstimatedParameters(tl_old->getParameters());
      
      // Write parameters to file:
      string parametersFile = ApplicationTools::getAFilePath("output.estimates", bppml.getParams(), false, false);
      ApplicationTools::displayResult("Output estimates to file", parametersFile);
      if (parametersFile != "none")
      {
        StlOutputStream out(new ofstream(parametersFile.c_str(), ios::out));
        out << "# Log likelihood = ";
        out.setPrecision(20) << (-tl_old->getValue());
        out.endLine();
        out << "# Number of sites = ";
        out.setPrecision(20) << vSites[0]->getNumberOfSites();
        out.endLine();
        out.endLine();
        out << "# Substitution model parameters:";
        out.endLine();
        if (modelSet)
        {
          modelSet->matchParametersValues(tl_old->getParameters());
          PhylogeneticsApplicationTools::printParameters(modelSet, out);
        }
        else
        {
          model->matchParametersValues(tl_old->getParameters());
          PhylogeneticsApplicationTools::printParameters(model, out);
        }
        out.endLine();
        (out << "# Rate distribution parameters:").endLine();
        rDist->matchParametersValues(tl_old->getParameters());
        PhylogeneticsApplicationTools::printParameters(rDist, out);
      }
        
      // Getting posterior rate class distribution:
      DiscreteDistribution* prDist = RASTools::getPosteriorRateDistribution(*tl_old);
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
        vector<size_t> classes = tl_old->getRateClassWithMaxPostProbOfEachSite();
          
        // Get the posterior rate, i.e. rate averaged over all posterior probabilities:
        Vdouble rates = tl_old->getPosteriorRateOfEachSite();
          
        vector<string> colNames;
        colNames.push_back("Sites");
        colNames.push_back("is.complete");
        colNames.push_back("is.constant");
        colNames.push_back("lnL");
        colNames.push_back("rc");
        colNames.push_back("pr");
        vector<string> row(6);
        DataTable* infos = new DataTable(colNames);
          
        for (unsigned int i = 0; i < vSites[0]->getNumberOfSites(); i++)
        {
          double lnL = tl_old->getLogLikelihoodForASite(i);
          const Site* currentSite = &vSites[0]->getSite(i);
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
    }
    // New optimization
    else // tl_new!=0 
    {
      //Check initial likelihood:
      double logL = tl_new->getValue();
      cerr << logL << endl;
      
      if (isinf(logL))
      {
        // This may be due to null branch lengths, leading to null likelihood!
        ApplicationTools::displayWarning("!!! Warning!!! Initial likelihood is zero.");
        ApplicationTools::displayWarning("!!! This may be due to branch length == 0.");
        ApplicationTools::displayWarning("!!! All null branch lengths will be set to 0.000001.");
        ParameterList pl = tl_new->getBranchLengthsParameters();
        for (unsigned int i = 0; i < pl.size(); i++)
        {
          if (pl[i].getValue() < 0.000001) pl[i].setValue(0.000001);
        }
        tl_new->matchParametersValues(pl);
        logL = tl_new->getValue();
      }
      ApplicationTools::displayResult("Initial log likelihood", TextTools::toString(-logL, 15));
      if (isinf(logL))
      {
        ApplicationTools::displayError("!!! Unexpected initial likelihood == 0.");

        vector<SingleDataPhyloLikelihood*> vSD;
        if (dynamic_cast<SingleDataPhyloLikelihood*>(tl_new)!=NULL)
          vSD.push_back(dynamic_cast<SingleDataPhyloLikelihood*>(tl_new));
        else{
          MultiDataPhyloLikelihood* mDP=dynamic_cast<MultiDataPhyloLikelihood*>(tl_new);
          
          for (size_t nSD=0; nSD< mDP->getNumberOfSingleDataPhyloLikelihoods(); nSD++)
            vSD.push_back(mDP->getSingleDataPhylolikelihood(nSD));
        }

        for (size_t nSD=0; nSD != vSD.size(); nSD++)
        {
          ApplicationTools::displayWarning("Checking for phylolikelihood " + TextTools::toString(nSD));
          
          if (isinf(vSD[nSD]->getValue()))
          {
            SingleDataPhyloLikelihood* sDP=dynamic_cast<SingleDataPhyloLikelihood*>(vSD[nSD]);
            /// !!! Not economic
            SiteContainer* vData=sDP->getData()->clone();
            
            if (codonAlphabet)
            {
              bool f = false;
              size_t s;
              for (size_t i = 0; i < vData->getNumberOfSites(); i++) {
                if (isinf(sDP->getLogLikelihoodForASite(i))) {
                  const Site& site = vData->getSite(i);
                  s = site.size();
                  for (size_t j = 0; j < s; j++) {
                    if (gCode->isStop(site.getValue(j))) {
                      (*ApplicationTools::error << "Stop Codon at site " << site.getPosition() << " in sequence " << vData->getSequence(j).getName()).endLine();
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
              for (unsigned int i = 0; i < vData->getNumberOfSites(); i++) {
                (*ApplicationTools::error << "Site " << vData->getSite(i).getPosition() << "\tlog likelihood = " << sDP->getLogLikelihoodForASite(i)).endLine();
              }
              ApplicationTools::displayError("!!! 0 values (inf in log) may be due to computer overflow, particularily if datasets are big (>~500 sequences).");
              ApplicationTools::displayError("!!! You may want to try input.sequence.remove_saturated_sites = yes to ignore positions with likelihood 0.");
              exit(1);
            } else {
              ApplicationTools::displayBooleanResult("Saturated site removal enabled", true);
              for (size_t i = vData->getNumberOfSites(); i > 0; --i) {
                if (isinf(sDP->getLogLikelihoodForASite(i - 1))) {
                  ApplicationTools::displayResult("Ignore saturated site", vData->getSite(i - 1).getPosition());
                  vData->deleteSite(i - 1);
                }
              }
              ApplicationTools::displayResult("Number of sites retained", vSites[0]->getNumberOfSites());

              sDP->setData(*vData);
              logL = sDP->getValue();
              if (isinf(logL)) {
                ApplicationTools::displayError("This should not happen. Exiting now.");
                exit(1);
              }
              ApplicationTools::displayResult("Initial log likelihood", TextTools::toString(-logL, 15));
            }
          }
        }
      }
      
      tl_new = PhylogeneticsApplicationTools::optimizeParameters(tl_new, tl_new->getParameters(), bppml.getParams());

      std::vector<const TreeTemplate<Node>* > vTNree = SPC->getTrees();
      PhylogeneticsApplicationTools::writeTrees(vTNree, bppml.getParams());
      
      // Write parameters to screen:
      ApplicationTools::displayResult("Log likelihood", TextTools::toString(-tl_new->getValue(), 15));
      ParameterList parameters = tl_new->getParameters();
      parameters.deleteParameters(tl_new->getBranchLengthsParameters().getParameterNames(),false);
              
      for (unsigned int i = 0; i < parameters.size(); i++)
        ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
      
      // Checking convergence:
      PhylogeneticsApplicationTools::checkEstimatedParameters(tl_new->getParameters());
      
      // Write parameters to file:
      string parametersFile = ApplicationTools::getAFilePath("output.estimates", bppml.getParams(), false, false);
      ApplicationTools::displayResult("Process estimates to file", parametersFile);
      if (parametersFile != "none")
      {
        StlOutputStream out(new ofstream(parametersFile.c_str(), ios::out));
        PhylogeneticsApplicationTools::printParameters(tl_new, out);
      }

      // Write infos to file:
      //     probabilities of rate discrete distributions
      //     site infos : lnL, class (or process in case of collection) posterior probability distribution

      string infosFile = ApplicationTools::getAFilePath("output.infos", bppml.getParams(), false, false);
      if (infosFile != "none")
      {
        ApplicationTools::displayResult("Alignment information logfile", infosFile);
        StlOutputStream out(new ofstream(infosFile.c_str(), ios::out));
                  
        PhylogeneticsApplicationTools::printAnalysisInformation(tl_new, out);
      }
              
      ////////////////////////////////////////////
      // Bootstrap:
    
      string optimizeClock = ApplicationTools::getStringParameter("optimization.clock", bppml.getParams(), "None", "", true, false);
      if (nbBS > 0 && optimizeClock != "None")
      {
        ApplicationTools::displayError("Bootstrap is not supported with clock trees.");
      }
      if (nbBS > 0 && optimizeClock == "None")
      {
        ApplicationTools::displayResult("Number of bootstrap samples", TextTools::toString(nbBS));
        bool approx = ApplicationTools::getBooleanParameter("bootstrap.approximate", bppml.getParams(), true);
        ApplicationTools::displayResult("Use approximate bootstrap", TextTools::toString(approx ? "yes" : "no"));
        bool bootstrapVerbose = ApplicationTools::getBooleanParameter("bootstrap.verbose", bppml.getParams(), false, "", true, false);

        const Tree* initTree = vTree[0];
        if (!bootstrapVerbose) bppml.getParam("optimization.verbose") = "0";
        bppml.getParam("optimization.profiler") = "none";
        bppml.getParam("optimization.messageHandler") = "none";
        if (!optimizeTopo)
        {
          bppml.getParam("optimization.topology") = "yes";
          tl_old = dynamic_cast<NNIHomogeneousTreeLikelihood*>(
            PhylogeneticsApplicationTools::optimizeParameters(tl_old, tl_old->getParameters(), bppml.getParams(), "", true, false));
          initTree = &tl_old->getTree();
        }

        string bsTreesPath = ApplicationTools::getAFilePath("bootstrap.output.file", bppml.getParams(), false, false);
        ofstream* out = 0;
        if (bsTreesPath != "none")
        {
          ApplicationTools::displayResult("Bootstrap trees stored in file", bsTreesPath);
          out = new ofstream(bsTreesPath.c_str(), ios::out);
        }
        Newick newick;
        ParameterList paramsToIgnore = tl_old->getSubstitutionModelParameters();
        paramsToIgnore.addParameters(tl_old->getRateDistributionParameters());

        ApplicationTools::displayTask("Bootstrapping", true);
        vector<Tree*> bsTrees(nbBS);
        for (unsigned int i = 0; i < nbBS; i++)
        {
          ApplicationTools::displayGauge(i, nbBS - 1, '=');
          VectorSiteContainer* sample = SiteContainerTools::bootstrapSites(*vSites[0]);
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
        TreeTools::computeBootstrapValues(*vTree[0], bsTrees);
        ApplicationTools::displayTaskDone();
        for (unsigned int i = 0; i < nbBS; i++)
        {
          delete bsTrees[i];
        }

        // Write resulting tree:
        vcTree.clear();
        for (size_t i=0; i<vTree.size(); i++)
          vcTree.push_back(vTree[i]);
      
        PhylogeneticsApplicationTools::writeTrees(vcTree, bppml.getParams());
      }


      delete alphabet;
      for (size_t i=0; i< vSites.size(); i++)
        delete vSites[i];
      
      if (model) delete model;
      if (modelSet) delete modelSet;
      delete rDist;
      if (tl_old)
        delete tl_old;
      if (tl_new)
        delete tl_new;
      bppml.done();
    }
  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    return 1;
  }
  
return 0;
}

