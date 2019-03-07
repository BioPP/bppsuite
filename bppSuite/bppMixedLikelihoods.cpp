//
// File: bppMixedLikelihoods.cpp
// Created by: Laurent Guéguen
// Created on: lundi 12 novembre 2012, à 07h 02
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
#include <Bpp/App/NumCalcApplicationTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/App/NumCalcApplicationTools.h>

// From bpp-seq:
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SequenceContainerTools.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

// From PhylLib:
#include <Bpp/Phyl/Tree/Tree.h>
#include <Bpp/Phyl/Tree/TreeTemplate.h>
#include <Bpp/Phyl/PatternTools.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
#include <Bpp/Phyl/Model/AbstractBiblioMixedTransitionModel.h>
#include <Bpp/Phyl/Model/MixedTransitionModel.h>
#include <Bpp/Phyl/Model/MixtureOfATransitionModel.h>
#include <Bpp/Phyl/Model/MixtureOfTransitionModels.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/RNonHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Phyl/Io/Newick.h>

using namespace bpp;

/******************************************************************************/

void help()
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
  (*ApplicationTools::message << "bppmixedlikelihoods parameter1_name=parameter1_value ").endLine();
  (*ApplicationTools::message << "      parameter2_name=parameter2_value ... param=option_file").endLine();
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "  Refer to the Bio++ Program Suite Manual for a list of available options.").endLine();
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
}

int main(int args, char** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*     Bio++ Computation of site likelihoods inside mixed models  *" << endl;
  cout << "*                        Version " << BPP_VERSION << ".                          *" << endl;
  cout << "* Author: L. Guéguen                       Last Modif.: " << BPP_REL_DATE << " *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  if (args == 1)
  {
    help();
    return 0;
  }

  try
  {
    BppApplication bppmixedlikelihoods(args, argv, "BppMixedLikelihoods");
    bppmixedlikelihoods.startTimer();

    Alphabet* alphabet = SequenceApplicationTools::getAlphabet(bppmixedlikelihoods.getParams(), "", false);
    unique_ptr<GeneticCode> gCode;
    CodonAlphabet* codonAlphabet = dynamic_cast<CodonAlphabet*>(alphabet);
    if (codonAlphabet) {
      string codeDesc = ApplicationTools::getStringParameter("genetic_code", bppmixedlikelihoods.getParams(), "Standard", "", true, true);
      ApplicationTools::displayResult("Genetic Code", codeDesc);
      
      gCode.reset(SequenceApplicationTools::getGeneticCode(codonAlphabet->getNucleicAlphabet(), codeDesc));
    }

    // get the data

    AlignedValuesContainer* allSites = SequenceApplicationTools::getAlignedContainer(alphabet, bppmixedlikelihoods.getParams());

    AlignedValuesContainer* sites = SequenceApplicationTools::getSitesToAnalyse(*allSites, bppmixedlikelihoods.getParams(), "", true, false);
    delete allSites;

    ApplicationTools::displayResult("Number of sequences", TextTools::toString(sites->getNumberOfSequences()));
    ApplicationTools::displayResult("Number of sites", TextTools::toString(sites->getNumberOfSites()));

    // Get the tree
    Tree* tree = PhylogeneticsApplicationTools::getTree(bppmixedlikelihoods.getParams());
    ApplicationTools::displayResult("Number of leaves", TextTools::toString(tree->getNumberOfLeaves()));


    AbstractDiscreteRatesAcrossSitesTreeLikelihood* tl;
    string nhOpt = ApplicationTools::getStringParameter("nonhomogeneous", bppmixedlikelihoods.getParams(), "no", "", true, false);
    ApplicationTools::displayResult("Heterogeneous model", nhOpt);

    MixedTransitionModel* model       = 0;
    MixedSubstitutionModelSet* modelSet = 0;
    DiscreteDistribution* rDist         = 0;

    map<string, string> unparsedparams;

    if (nhOpt == "no")
    {
      model = dynamic_cast<MixedTransitionModel*>(PhylogeneticsApplicationTools::getTransitionModel(alphabet, gCode.get(), sites, bppmixedlikelihoods.getParams(), unparsedparams));
      if (model == 0)
      {
        cout << "Model is not a Mixed model" << endl;
        exit(0);
      }

      SiteContainerTools::changeGapsToUnknownCharacters(*sites);
      if (model->getNumberOfStates() > model->getAlphabet()->getSize())
      {
        // Markov-modulated Markov model!
        rDist = new ConstantRateDistribution();
      }
      else
      {
        rDist = PhylogeneticsApplicationTools::getRateDistribution(bppmixedlikelihoods.getParams());
      }
      tl = new RHomogeneousMixedTreeLikelihood(*tree, *sites, model, rDist, true);
    }
    else if (nhOpt == "one_per_branch")
    {
      model = dynamic_cast<MixedTransitionModel*>(PhylogeneticsApplicationTools::getTransitionModel(alphabet, gCode.get(), sites, bppmixedlikelihoods.getParams(), unparsedparams));
      if (model == 0)
      {
        cout << "Model is not a Mixed model" << endl;
        exit(0);
      }

      SiteContainerTools::changeGapsToUnknownCharacters(*sites);
      if (model->getNumberOfStates() > model->getAlphabet()->getSize())
      {
        // Markov-modulated Markov model!
        rDist = new ConstantRateDistribution();
      }
      else
      {
        rDist = PhylogeneticsApplicationTools::getRateDistribution(bppmixedlikelihoods.getParams());
      }
      vector<double> rateFreqs;
      if (model->getNumberOfStates() != alphabet->getSize())
      {
        // Markov-Modulated Markov Model...
        unsigned int n = (unsigned int)(model->getNumberOfStates() / alphabet->getSize());
        rateFreqs = vector<double>(n, 1. / (double)n); // Equal rates assumed for now, may be changed later (actually, in the most general case,
        // we should assume a rate distribution for the root also!!!
      }
      
      std::map<std::string, std::string> aliasFreqNames;

      FrequenciesSet* rootFreqs = PhylogeneticsApplicationTools::getRootFrequenciesSet(alphabet, gCode.get(), sites, bppmixedlikelihoods.getParams(), aliasFreqNames, rateFreqs);

      string descGlobal = ApplicationTools::getStringParameter("nonhomogeneous_one_per_branch.shared_parameters", bppmixedlikelihoods.getParams(), "", "", true, 1);

      NestedStringTokenizer nst(descGlobal,"[","]",",");
      const deque<string>& descGlobalParameters=nst.getTokens();

      map<string, vector<Vint> > globalParameters;
      for (const auto& desc:descGlobalParameters)
      {
        size_t post=desc.rfind("_");
        if (post==std::string::npos || post==desc.size()-1 || desc[post+1]!='[')
          globalParameters[desc]={};
        else
        {
          string key=desc.substr(0,post);
          Vint sint=NumCalcApplicationTools::seqFromString(desc.substr(post+2, desc.size()-post-3));
          if (globalParameters.find(key)==globalParameters.end())
            globalParameters[key]=vector<Vint>(1, sint);
          else
            globalParameters[key].push_back(sint);
        }
      }

      for (const auto& globpar:globalParameters)
      {
        ApplicationTools::displayResult("Global parameter", globpar.first);
        if (globpar.second.size()==0)
        {
          string all="All nodes";
          ApplicationTools::displayResult(" shared between nodes", all);
        }
        else
          for (const auto& vint:globpar.second)
            ApplicationTools::displayResult(" shared between nodes", VectorTools::paste(vint,","));
      }
      
      modelSet = dynamic_cast<MixedSubstitutionModelSet*>(SubstitutionModelSetTools::createNonHomogeneousModelSet(model, rootFreqs, tree, aliasFreqNames, globalParameters));
      model = 0;
      tl = new RNonHomogeneousMixedTreeLikelihood(*tree, *sites, modelSet, rDist, true);
    }
    else if (nhOpt == "general")
    {
      modelSet = dynamic_cast<MixedSubstitutionModelSet*>(PhylogeneticsApplicationTools::getSubstitutionModelSet(alphabet, gCode.get(), sites, bppmixedlikelihoods.getParams()));
      if (modelSet == 0)
      {
        cout << "Missing a Mixed model" << endl;
        exit(0);
      }

      SiteContainerTools::changeGapsToUnknownCharacters(*sites);
      if (modelSet->getNumberOfStates() > modelSet->getAlphabet()->getSize())
      {
        // Markov-modulated Markov model!
        rDist = new ConstantDistribution(1.);
      }
      else
      {
        rDist = PhylogeneticsApplicationTools::getRateDistribution(bppmixedlikelihoods.getParams());
      }
      tl = new RNonHomogeneousMixedTreeLikelihood(*tree, *sites, modelSet, rDist, true);
    }
    else
      throw Exception("Unknown option for nonhomogeneous: " + nhOpt);

    tl->initialize();

    double logL = tl->getValue();
    if (std::isinf(logL))
    {
      // This may be due to null branch lengths, leading to null likelihood!
      ApplicationTools::displayWarning("!!! Warning!!! Likelihood is zero.");
      ApplicationTools::displayWarning("!!! This may be due to branch length == 0.");
      ApplicationTools::displayWarning("!!! All null branch lengths will be set to 0.000001.");
      ParameterList pl = tl->getBranchLengthsParameters();
      for (unsigned int i = 0; i < pl.size(); i++)
      {
        if (pl[i].getValue() < 0.000001)
          pl[i].setValue(0.000001);
      }
      tl->matchParametersValues(pl);
      logL = tl->getValue();
    }
    if (std::isinf(logL))
    {
      ApplicationTools::displayError("!!! Unexpected likelihood == 0.");
      ApplicationTools::displayError("!!! Looking at each site:");
      for (unsigned int i = 0; i < sites->getNumberOfSites(); i++)
      {
        (*ApplicationTools::error << "Site " << sites->getSymbolListSite(i).getPosition() << "\tlog likelihood = " << tl->getLogLikelihoodForASite(i)).endLine();
      }
      ApplicationTools::displayError("!!! 0 values (inf in log) may be due to computer overflow, particularily if datasets are big (>~500 sequences).");
      exit(-1);
    }


    // Write parameters to screen:
    ApplicationTools::displayResult("Log likelihood", TextTools::toString(tl->getValue(), 15));
    ParameterList parameters = tl->getSubstitutionModelParameters();
    for (size_t i = 0; i < parameters.size(); i++)
    {
      ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
    }
    parameters = tl->getRateDistributionParameters();
    for (size_t i = 0; i < parameters.size(); i++)
    {
      ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
    }


    // /////////////////////////////////////////////
    // Getting likelihoods per submodel

    string outputFile;
    outputFile = ApplicationTools::getAFilePath("output.likelihoods.file", bppmixedlikelihoods.getParams(), true, false);
    ApplicationTools::displayResult("Output file for likelihoods", outputFile);
    ofstream out(outputFile.c_str(), ios::out);

    size_t nSites = sites->getNumberOfSites();

    size_t nummodel = ApplicationTools::getParameter<size_t>("likelihoods.model_number", bppmixedlikelihoods.getParams(), 1, "", true, true);

    string parname = ApplicationTools::getStringParameter("likelihoods.parameter_name", bppmixedlikelihoods.getParams(), "", "", true, false);

    if (modelSet && ((nummodel <= 0) || (nummodel > modelSet->getNumberOfModels())))
    {
      ApplicationTools::displayError("Bad number of model " + TextTools::toString(nummodel) + ".");
      exit(-1);
    }

    MixedTransitionModel* p0 = dynamic_cast<MixedTransitionModel*>(model ? model : modelSet->getModel(nummodel - 1));

    if (!p0)
    {
      ApplicationTools::displayError("Model " + TextTools::toString(nummodel) + " is not a Mixed Model.");
      exit(-1);
    }

    const AbstractBiblioMixedTransitionModel* ptmp = dynamic_cast<const AbstractBiblioMixedTransitionModel*>(p0);
    if (ptmp) {
      p0 = ptmp->getMixedModel().clone();

      if (nhOpt == "no")
        model = p0;
      else {
        modelSet->replaceModel(nummodel-1, p0);
        modelSet->isFullySetUpFor(*tree);
      }
    }
    
    //////////////////////////////////////////////////
    // Case of a MixtureOfSubstitutionModels

    MixtureOfTransitionModels* pMSM = dynamic_cast<MixtureOfTransitionModels*>(p0);

    if (pMSM)
    {
      vector<string> colNames;
      colNames.push_back("Sites");

      size_t nummod = pMSM->getNumberOfModels();
      for (unsigned int i = 0; i < nummod; i++)
      {
        colNames.push_back(pMSM->getNModel(i)->getName());
      }

      DataTable* rates = new DataTable(nSites, colNames.size());
      rates->setColumnNames(colNames);

      for (unsigned int i = 0; i < nSites; i++)
      {
        const CruxSymbolListSite& currentSite = sites->getSymbolListSite(i);
        int currentSitePosition = currentSite.getPosition();
        (*rates)(i, "Sites") = string("[" + TextTools::toString(currentSitePosition) + "]");
      }

      Vdouble vprob = pMSM->getProbabilities();
      for (unsigned int i = 0; i < nummod; i++)
      {
        string modname = pMSM->getNModel(i)->getName();

        for (unsigned int j = 0; j < nummod; j++)
        {
          pMSM->setNProbability(j, (j == i) ? 1 : 0);
        }

        if (tl)
          delete tl;

        if (nhOpt == "no")
          tl = new RHomogeneousMixedTreeLikelihood(*tree, *sites, model, rDist, true, false, true);
        else
          tl = new RNonHomogeneousMixedTreeLikelihood(*tree, *sites, modelSet, rDist, false, true);

        tl->initialize();
        logL = tl->getValue();
        Vdouble Vd = tl->getLogLikelihoodPerSite();
        for (unsigned int j = 0; j < nSites; j++)
        {
          (*rates)(j, modname) = TextTools::toString(Vd[j]);
        }

        ApplicationTools::displayMessage("\n");
        ApplicationTools::displayMessage("Model " + modname + ":");
        ApplicationTools::displayResult("Log likelihood", TextTools::toString(tl->getValue(), 15));
        ApplicationTools::displayResult("Probability", TextTools::toString(vprob[i], 15));
      }

      DataTable::write(*rates, out, "\t");
    }

    //////////////////////////////////////////////////
    // Case of a MixtureOfASubstitutionModel

    else
    {
      MixtureOfATransitionModel* pMSM2 = dynamic_cast<MixtureOfATransitionModel*>(p0);
      if (pMSM2 != NULL)
      {
        size_t nummod = pMSM2->getNumberOfModels();
        if (parname == "")
        {
          ParameterList pl=pMSM2->getParameters();

          for (size_t i2 = 0; i2 < pl.size(); i2++)
          {
            string pl2n = pl[i2].getName();

            if (dynamic_cast<const ConstantDistribution*>(pMSM2->getDistribution(pl2n))==NULL)
            {
              parname=pl2n;

              while (parname.size()>0 && pMSM2->getDistribution(parname)==NULL)
                parname=pl2n.substr(0,pl2n.rfind("_"));

              if (parname.size()>0){
                ApplicationTools::displayResult("likelihoods.parameter_name", parname);
                break;
              }
            }
          }
        }

        if (parname == "")
        {
          ApplicationTools::displayError("Argument likelihoods.parameter_name is required.");
          exit(-1);
        }

        vector< Vint > vvnmod;
        size_t i2 = 0;
        while (i2 < nummod)
        {
          string par2 = parname + "_" + TextTools::toString(i2 + 1);
          Vint vnmod = pMSM2->getSubmodelNumbers(par2);
          if (vnmod.size() == 0)
            break;
          vvnmod.push_back(vnmod);
          i2++;
        }

        size_t nbcl = vvnmod.size();
        if (nbcl==0)
          throw Exception("Parameter " + parname + " is not mixed.");
        
        Vdouble vprob = pMSM2->getProbabilities();

        vector<vector<double> > vvprob;
        vector<double> vsprob;
        
        for (size_t i = 0; i < nbcl; i++)
        {
          vector<double> vprob2;
          for (size_t j = 0; j < vvnmod[i].size(); j++)
          {
            vprob2.push_back(vprob[static_cast<size_t>(vvnmod[i][j])]);
          }

          vvprob.push_back(vprob2);
          vsprob.push_back(VectorTools::sum(vvprob[i]));
        }

        vector<string> colNames;
        colNames.push_back("Sites");

        Vdouble dval;
        for (size_t i = 0; i < nbcl; i++)
        {
          const TransitionModel* pSM = pMSM2->getNModel(static_cast<size_t>(vvnmod[i][0]));
          double valPar = pSM->getParameterValue(pSM->getParameterNameWithoutNamespace(parname));
          dval.push_back(valPar);
          colNames.push_back("Ll_" + parname + "=" + TextTools::toString(valPar));
        }
        for (size_t i = 0; i < nbcl; i++)
          colNames.push_back("Pr_" + parname + "=" + TextTools::toString(dval[i]));

        colNames.push_back("mean");

        DataTable* rates = new DataTable(nSites, colNames.size());
        rates->setColumnNames(colNames);

        for (size_t i = 0; i < nSites; i++)
        {
          const CruxSymbolListSite& currentSite = sites->getSymbolListSite(i);
          int currentSitePosition = currentSite.getPosition();
          (*rates)(i,"Sites")=TextTools::toString(currentSitePosition);
        }

        VVdouble vvd;

          
        vector<double> vRates = pMSM2->getVRates();

        for (size_t i = 0; i < nbcl; ++i)
        {
          string par2 = parname + "_" + TextTools::toString(i + 1);
          
          for (unsigned int j = 0; j < nummod; ++j)
            pMSM2->setNProbability(j, 0);

          for (size_t j = 0; j < vvprob[i].size(); ++j)
            pMSM2->setNProbability(static_cast<size_t>(vvnmod[i][j]), vvprob[i][j] / vsprob[i]);

          if (tl)
            delete tl;

          if (nhOpt == "no")
            tl = new RHomogeneousMixedTreeLikelihood(*tree, *sites, model, rDist, true, false, true);
          else
            tl = new RNonHomogeneousMixedTreeLikelihood(*tree, *sites, modelSet, rDist, false, true);

          tl->initialize();
          logL = tl->getValue();
          Vdouble vd = tl->getLogLikelihoodPerSite();

          for (unsigned int j = 0; j < nSites; j++)
            (*rates)(j, i + 1) = TextTools::toString(vd[j]);

          vvd.push_back(vd);

          ApplicationTools::displayMessage("\n");
          ApplicationTools::displayMessage("Parameter " + par2 + "=" + TextTools::toString(dval[i]) + " with rate=" + TextTools::toString(vRates[i]));

          ApplicationTools::displayResult("Log likelihood", TextTools::toString(tl->getValue(), 15));
          ApplicationTools::displayResult("Probability", TextTools::toString(vsprob[i], 15));
        }

        for (size_t j = 0; j < nSites; j++)
        {
          Vdouble vd;
          for (size_t i = 0; i < nbcl; i++)
            vd.push_back(std::log(vsprob[i])+vvd[i][j]);
          
          VectorTools::logNorm(vd);
          for (size_t i = 0; i < nbcl; i++)
            (*rates)(j,nbcl + i + 1) = TextTools::toString(std::exp(vd[i]));
          (*rates)(j, 2 * nbcl + 1) = TextTools::toString(VectorTools::sumExp(vd, dval));
        }

        DataTable::write(*rates, out, "\t");
      }
    }

    delete alphabet;
    delete sites;
    if (model)
      delete model;
    if (modelSet)
      delete modelSet;
    delete rDist;
    delete tl;
    delete tree;
    ApplicationTools::displayMessage("\n");
    bppmixedlikelihoods.done();
  }

  catch (exception& e)
  {
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}

