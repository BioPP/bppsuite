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

// From bpp-core:
#include <Bpp/Version.h>
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
#include <Bpp/Phyl/Tree/TreeTemplate.h>
#include <Bpp/Phyl/Likelihood/RHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/DRHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/RNonHomogeneousMixedTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/DRNonHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/RASTools.h>

#include <Bpp/Phyl/PatternTools.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
#include <Bpp/Phyl/Model/MixedSubstitutionModel.h>
#include <Bpp/Phyl/Model/Protein/CoalaCore.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Model/FrequenciesSet/MvaFrequenciesSet.h>
#include <Bpp/Phyl/Io/Newick.h>

#include <Bpp/Phyl/NewLikelihood/PhyloLikelihoods/PhyloLikelihoodContainer.h>
#include <Bpp/Phyl/NewLikelihood/PhyloLikelihoods/SetOfAbstractPhyloLikelihood.h>

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
  cout << "*       Bio++ Maximum Likelihood Computation, version " << BPP_VERSION << "      *" << endl;
  cout << "*                                                                *" << endl;
  cout << "* Authors: J. Dutheil                       Last Modif. " << BPP_REL_DATE << " *" << endl;
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

    map<string, string> unparsedparams;

    ///// Alphabet

    Alphabet* alphabet = SequenceApplicationTools::getAlphabet(bppml.getParams(), "", false);
    unique_ptr<GeneticCode> gCode;
    CodonAlphabet* codonAlphabet = dynamic_cast<CodonAlphabet*>(alphabet);
    if (codonAlphabet) {
      string codeDesc = ApplicationTools::getStringParameter("genetic_code", bppml.getParams(), "Standard", "", true, true);
      ApplicationTools::displayResult("Genetic Code", codeDesc);

      gCode.reset(SequenceApplicationTools::getGeneticCode(codonAlphabet->getNucleicAlphabet(), codeDesc));
    }


    ////// Get the map of the sequences

    map<size_t, AlignedValuesContainer*> mSites = SequenceApplicationTools::getAlignedContainers(alphabet, bppml.getParams());


    if (mSites.size() == 0)
      throw Exception("Missing data input.sequence.file option");

    /////// Get the map of initial trees

    map<size_t, PhyloTree*> mpTree=PhylogeneticsApplicationTools::getPhyloTrees(bppml.getParams(), mSites, unparsedparams);

    // Try to write the current tree to file. This will be overwritten
    // by the optimized tree, but allow to check file existence before
    // running optimization!

    map<size_t, PhyloTree*>::const_iterator itp;
    map<size_t, Tree*>::const_iterator it;
    vector<const PhyloTree*> vcpTree;
    
    for (itp = mpTree.begin(); itp != mpTree.end(); itp++)
      vcpTree.push_back(itp->second);
    
    PhylogeneticsApplicationTools::writeTrees(vcpTree, bppml.getParams());

    map<size_t, Tree*> mTree=PhylogeneticsApplicationTools::getTrees(bppml.getParams(), mSites, unparsedparams);
    
    string treeWIdPath = ApplicationTools::getAFilePath("output.tree_ids.file", bppml.getParams(), false, false, "", true, "none", 1);

    if (treeWIdPath != "none")
    {
      Newick treeWriter;
      // treeWriter.enableExtendedBootstrapProperty("NodeId");
      // ApplicationTools::displayResult("Writing tagged tree to", treeWIdPath);

      for (itp = mpTree.begin(); itp != mpTree.end(); itp++)
      {
        PhyloTree* ttree=itp->second;
        vector<std::shared_ptr<PhyloNode> > nodes = ttree->getAllNodes();
        
        for (size_t i = 0; i < nodes.size(); i++)
        {
          if (ttree->isLeaf(nodes[i]))
            nodes[i]->setName(TextTools::toString(ttree->getNodeIndex(nodes[i]) + "_" + nodes[i]->getName()));
          else
          {
            shared_ptr<PhyloBranch> branch=ttree->hasFather(nodes[i])?ttree->getEdgeToFather(nodes[i]):0;
            if (branch)
              branch->setProperty("BranchId", BppString(TextTools::toString(ttree->getNodeIndex(nodes[i]))));
          }
        }

        treeWriter.write(*itp->second, treeWIdPath, itp==mpTree.begin());
        delete itp->second;
      }
      cout << "BppML's done." << endl;
      exit(0);
    }

    /////////////////
    // Computing stuff

    DiscreteRatesAcrossSitesTreeLikelihood* tl_old = 0;
    PhyloLikelihood* tl_new = 0;
    SubstitutionProcessCollection* SPC = 0;
    map<size_t, SequenceEvolution*> mSeqEvol;
    PhyloLikelihoodContainer* mPhyl=0;
    
    bool checkTree    = ApplicationTools::getBooleanParameter("input.tree.check_root", bppml.getParams(), true, "", true, 2);
    bool optimizeTopo = ApplicationTools::getBooleanParameter("optimization.topology", bppml.getParams(), false, "", true, 1);
    unsigned int nbBS = ApplicationTools::getParameter<unsigned int>("bootstrap.number", bppml.getParams(), 0, "", true, 1);
    string collection = ApplicationTools::getStringParameter("collection", bppml.getParams(), "", "", true, 1);


    TransitionModel*    model    = 0;
    SubstitutionModelSet* modelSet = 0;
    DiscreteDistribution* rDist    = 0;

    Tree* firstTree = mTree.begin()->second;

    /// Topology estimation

    if (optimizeTopo || nbBS > 0)
    {
      if (collection != "")
        throw Exception("Topology estimation in collections not supported yet, sorry.");

      string nhOpt = ApplicationTools::getStringParameter("nonhomogeneous", bppml.getParams(), "no", "", true, false);
      ApplicationTools::displayResult("Heterogeneous model", nhOpt);

      if (nhOpt != "no")
        throw Exception("Topology estimation with NH model not supported yet, sorry :(");

      model = PhylogeneticsApplicationTools::getTransitionModels(alphabet, gCode.get(), mSites, bppml.getParams(), unparsedparams)[0];

      if (model->getName() != "RE08")
        for (auto  itc : mSites)
          SiteContainerTools::changeGapsToUnknownCharacters(*itc.second);

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
        tl_old = new NNIHomogeneousTreeLikelihood(*firstTree, *mSites.begin()->second, model, rDist, checkTree, true);
      else
        throw Exception("Topology estimation with Mixed model not supported yet, sorry :(");
    }

    /// Constant Topology

    else
    {
      map<size_t, DiscreteDistribution*> mDist = PhylogeneticsApplicationTools::getRateDistributions(bppml.getParams());

      map<size_t, TransitionModel*> mMod = PhylogeneticsApplicationTools::getTransitionModels(alphabet, gCode.get(), mSites, bppml.getParams(), unparsedparams);

      map<size_t, FrequenciesSet*> mRootFreq = PhylogeneticsApplicationTools::getRootFrequenciesSets(alphabet, gCode.get(), mSites, bppml.getParams(), unparsedparams);

      SPC=PhylogeneticsApplicationTools::getSubstitutionProcessCollection(alphabet, gCode.get(), mpTree, mMod, mRootFreq, mDist, bppml.getParams(), unparsedparams);

      mSeqEvol = PhylogeneticsApplicationTools::getSequenceEvolutions(*SPC, bppml.getParams(), unparsedparams);

      for (auto  itc : mSites)
        SiteContainerTools::changeGapsToUnknownCharacters(*itc.second);

      mPhyl=PhylogeneticsApplicationTools::getPhyloLikelihoodContainer(*SPC, mSeqEvol, mSites, bppml.getParams());

      // filter to Single Data PhyloLikelihoods

      if (!mPhyl->hasPhyloLikelihood(0))
        throw Exception("Missing phyloLikelihoods.");

      tl_new=(*mPhyl)[0];
      
    }
    
    ApplicationTools::displayMessage("");
    
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
        
      for (size_t i = 0; i < pl.size(); ++i) {
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
      if (std::isinf(logL))
      {
        // This may be due to null branch lengths, leading to null likelihood!
        ApplicationTools::displayWarning("!!! Warning!!! Initial likelihood is zero.");
        ApplicationTools::displayWarning("!!! This may be due to branch length == 0.");
        ApplicationTools::displayWarning("!!! All null branch lengths will be set to 0.001.");
        ParameterList pl = tl_old->getBranchLengthsParameters();
        for (unsigned int i = 0; i < pl.size(); i++)
        {
          if (pl[i].getValue() < 0.001)
            pl[i].setValue(0.001);
        }
        tl_old->matchParametersValues(pl);
        logL = tl_old->getValue();
      }
      ApplicationTools::displayResult("Initial log likelihood", TextTools::toString(-logL, 15));
      if (std::isinf(logL))
      {
        ApplicationTools::displayError("!!! Unexpected initial likelihood == 0.");
        if (codonAlphabet)
        {
          SiteContainer* vSC=dynamic_cast<SiteContainer*>(mSites.begin()->second);
          ProbabilisticSiteContainer* pSC=dynamic_cast<ProbabilisticSiteContainer*>(mSites.begin()->second);
          
          bool f = false;
          size_t s;
          for (size_t i = 0; i < mSites.begin()->second->getNumberOfSites(); i++) {
            if (std::isinf(tl_old->getLogLikelihoodForASite(i))) {
              if (vSC)
              {
                const Site& site = vSC->getSite(i);
                s = site.size();
                for (size_t j = 0; j < s; j++) {
                  if (gCode->isStop(site.getValue(j))) {
                    (*ApplicationTools::error << "Stop Codon at site " << site.getPosition() << " in sequence " << mSites.begin()->second->getName(j)).endLine();
                    f = true;
                  }
                }
              }
              else
              {
                const std::shared_ptr<ProbabilisticSite> site = pSC->getSite(i);
                s = site->size();
                for (size_t j = 0; j < s; j++)
                {
                  bool g=false;
                  for (int st = 0; !g && st < (int)alphabet->getSize(); st++)
                    g = (site->getStateValueAt(j,st)!=0 && !gCode->isStop(st));
                  
                  if (!g)
                  {
                    (*ApplicationTools::error << "Only stop Codons at site " << site->getPosition() << " in sequence " << mSites.begin()->second->getName(j)).endLine();
                    f = true;
                  }
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
          for (unsigned int i = 0; i < mSites.begin()->second->getNumberOfSites(); i++) {
            (*ApplicationTools::error << "Site " << mSites.begin()->second->getSymbolListSite(i).getPosition() << "\tlog likelihood = " << tl_old->getLogLikelihoodForASite(i)).endLine();
          }
          ApplicationTools::displayError("!!! 0 values (inf in log) may be due to computer overflow, particularily if datasets are big (>~500 sequences).");
          ApplicationTools::displayError("!!! You may want to try input.sequence.remove_saturated_sites = yes to ignore positions with likelihood 0.");
          exit(1);
        } else {
          ApplicationTools::displayBooleanResult("Saturated site removal enabled", true);
          for (size_t i = mSites.begin()->second->getNumberOfSites(); i > 0; --i) {
            if (std::isinf(tl_old->getLogLikelihoodForASite(i - 1))) {
              ApplicationTools::displayResult("Ignore saturated site", mSites.begin()->second->getSymbolListSite(i - 1).getPosition());
              mSites.begin()->second->deleteSites(i - 1, 1);
            }
          }
          ApplicationTools::displayResult("Number of sites retained", mSites.begin()->second->getNumberOfSites());
          tl_old->setData(*mSites.begin()->second);
          tl_old->initialize();
          logL = tl_old->getValue();
          if (std::isinf(logL)) {
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
      bool withAlias = ApplicationTools::getBooleanParameter("output.estimates.withalias", bppml.getParams(), true, "", false, 1);

      if (parametersFile != "none")
      {
        StlOutputStream out(new ofstream(parametersFile.c_str(), ios::out));
        out << "# Log likelihood = ";
        out.setPrecision(20) << (-tl_old->getValue());
        out.endLine();
        out << "# Number of sites = ";
        out.setPrecision(20) << mSites.begin()->second->getNumberOfSites();
        out.endLine();
        out.endLine();
        out << "# Substitution model parameters:";
        out.endLine();
        if (modelSet)
        {
          modelSet->matchParametersValues(tl_old->getParameters());
          PhylogeneticsApplicationTools::printParameters(modelSet, out, 1, withAlias);
        }
        else
        {
          model->matchParametersValues(tl_old->getParameters());
          PhylogeneticsApplicationTools::printParameters(model, out, 1);
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

        for (unsigned int i = 0; i < mSites.begin()->second->getNumberOfSites(); i++)
        {
          double lnL = tl_old->getLogLikelihoodForASite(i);
          const CruxSymbolListSite& currentSite = mSites.begin()->second->getSymbolListSite(i);
          int currentSitePosition = currentSite.getPosition();
          string isCompl = "NA";
          string isConst = "NA";
          try { isCompl = (SiteTools::isComplete(currentSite) ? "1" : "0"); }
          catch(EmptySiteException& ex) {}
          try { isConst = (SiteTools::isConstant(currentSite) ? "1" : "0"); }
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

      if (std::isinf(logL))
      {
        // This may be due to null branch lengths, leading to null likelihood!
        ApplicationTools::displayWarning("!!! Warning!!! Initial likelihood is zero.");
        ApplicationTools::displayWarning("!!! This may be due to branch length == 0.");
        ApplicationTools::displayWarning("!!! All null branch lengths will be set to 0.000001.");
        ParameterList pl = tl_new->getBranchLengthParameters();
        for (unsigned int i = 0; i < pl.size(); i++)
        {
          if (pl[i].getValue() < 0.000001) pl[i].setValue(0.001);
        }
        tl_new->matchParametersValues(pl);
        logL = tl_new->getValue();
      }
      ApplicationTools::displayResult("Initial log likelihood", TextTools::toString(-logL, 15));
      if (std::isinf(logL))
      {
        ApplicationTools::displayError("!!! Unexpected initial likelihood == 0.");

        map<size_t, AbstractSingleDataPhyloLikelihood*> mSD;

        if (dynamic_cast<AbstractSingleDataPhyloLikelihood*>(tl_new)!=NULL)
          mSD[1]=dynamic_cast<AbstractSingleDataPhyloLikelihood*>(tl_new);
        else{
          SetOfAbstractPhyloLikelihood* sOAP=dynamic_cast<SetOfAbstractPhyloLikelihood*>(tl_new);
          if (sOAP!=NULL)
          {
            const vector<size_t>& nSD=sOAP->getNumbersOfPhyloLikelihoods();
          
            for (size_t iSD=0; iSD< nSD.size(); iSD++)
            {
              AbstractSingleDataPhyloLikelihood* pASDP=dynamic_cast<AbstractSingleDataPhyloLikelihood*>(sOAP->getAbstractPhyloLikelihood(nSD[iSD]));
            
              if (pASDP!=NULL)
                mSD[nSD[iSD]]=pASDP;
            }
          }
        }

        for (map<size_t, AbstractSingleDataPhyloLikelihood*>::iterator itm=mSD.begin();itm!=mSD.end(); itm++)
        {
          ApplicationTools::displayWarning("Checking for phyloLikelihood " + TextTools::toString(itm->first));
          
          if (std::isinf(itm->second->getValue()))
          {
            AbstractSingleDataPhyloLikelihood* sDP=itm->second;
            /// !!! Not economic
            AlignedValuesContainer* vData=sDP->getData()->clone();

            SiteContainer* vSC=dynamic_cast<SiteContainer*>(vData);
            ProbabilisticSiteContainer* pSC=dynamic_cast<ProbabilisticSiteContainer*>(vData);

            if (codonAlphabet)
            {
              bool f = false;
              size_t s;
              for (size_t i = 0; i < vData->getNumberOfSites(); i++) {
                if (std::isinf(sDP->getLogLikelihoodForASite(i))) {
                  if (vSC)
                  {
                    const Site& site = vSC->getSite(i);
                    s = site.size();
                    for (size_t j = 0; j < s; j++) {
                      if (gCode->isStop(site.getValue(j))) {
                        (*ApplicationTools::error << "Stop Codon at site " << site.getPosition() << " in sequence " << vData->getName(j)).endLine();
                        f = true;
                      }
                    }
                  }
                  else
                  {
                    const std::shared_ptr<ProbabilisticSite> site = pSC->getSite(i);
                    s = site->size();
                    for (size_t j = 0; j < s; j++)
                    {
                      bool g=false;
                      for (int st = 0; !g && st < (int)alphabet->getSize(); st++)
                        g = (site->getStateValueAt(j,st)!=0 && !gCode->isStop(st));
                    
                      if (!g)
                      {
                        (*ApplicationTools::error << "Only stop Codons at site " << site->getPosition() << " in sequence " << vData->getName(j)).endLine();
                        f = true;
                      }
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
                (*ApplicationTools::error << "Site " << vData->getSymbolListSite(i).getPosition() << "\tlog likelihood = " << sDP->getLogLikelihoodForASite(i)).endLine();
              }
              ApplicationTools::displayError("!!! 0 values (inf in log) may be due to computer overflow, particularily if datasets are big (>~500 sequences).");
              ApplicationTools::displayError("!!! You may want to try input.sequence.remove_saturated_sites = yes to ignore positions with likelihood 0.");
              exit(1);
            } else {
              ApplicationTools::displayBooleanResult("Saturated site removal enabled", true);
              for (size_t i = vData->getNumberOfSites(); i > 0; --i) {
                if (std::isinf(sDP->getLogLikelihoodForASite(i - 1))) {
                  ApplicationTools::displayResult("Ignore saturated site", vData->getSymbolListSite(i - 1).getPosition());
                  vData->deleteSites(i - 1, i);
                }
              }
              ApplicationTools::displayResult("Number of sites retained", mSites.begin()->second->getNumberOfSites());

              sDP->setData(*vData);
              logL = sDP->getValue();
              if (std::isinf(logL)) {
                ApplicationTools::displayError("This should not happen. Exiting now.");
                exit(1);
              }
              ApplicationTools::displayResult("Initial log likelihood", TextTools::toString(-logL, 15));
            }
          }
        }
      }

      // First `true` means that default is to optimize model parameters.
      if(ApplicationTools::getBooleanParameter("optimization.model_parameters", bppml.getParams(), true, "", true, 1))
        tl_new = PhylogeneticsApplicationTools::optimizeParameters(tl_new, tl_new->getParameters(), bppml.getParams());
      else
        tl_new = PhylogeneticsApplicationTools::optimizeParameters(tl_new, tl_new->getBranchLengthParameters(), bppml.getParams());

      PhylogeneticsApplicationTools::writeTrees(*SPC, bppml.getParams());

      // Write parameters to screen:
      ApplicationTools::displayResult("Log likelihood", TextTools::toString(-tl_new->getValue(), 15));
      ParameterList parameters = tl_new->getParameters();
      parameters.deleteParameters(tl_new->getBranchLengthParameters().getParameterNames(),false);
              
      for (size_t i = 0; i < parameters.size(); i++)
        ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));

      // Checking convergence:
      PhylogeneticsApplicationTools::checkEstimatedParameters(tl_new->getParameters());

      // Write parameters to file:
      string parametersFile = ApplicationTools::getAFilePath("output.estimates", bppml.getParams(), false, false);
      bool withAlias = ApplicationTools::getBooleanParameter("output.estimates.withalias", bppml.getParams(), true, "", false, 1);

      ApplicationTools::displayResult("Process estimates to file", parametersFile);
      
      if (parametersFile != "none")
      {
        StlOutputStream out(new ofstream(parametersFile.c_str(), ios::out));

        PhylogeneticsApplicationTools::printParameters(*mPhyl, out);

        PhylogeneticsApplicationTools::printParameters(SPC, out, 1, withAlias);

        for (map<size_t, SequenceEvolution*>::const_iterator it2=mSeqEvol.begin(); it2!=mSeqEvol.end(); it2++)
        {
          PhylogeneticsApplicationTools::printParameters(it2->second, out, it2->first);
          out.endLine();
        }

      }

      // Write infos to file:
      //     probabilities of rate discrete distributions
      //     site infos : lnL, class (or process in case of collection) posterior probability distribution

      string infosFile = ApplicationTools::getAFilePath("output.infos", bppml.getParams(), false, false);
      if (infosFile != "none")
      {
        ApplicationTools::displayResult("Alignment information logfile", infosFile);
        PhylogeneticsApplicationTools::printAnalysisInformation(*mPhyl, infosFile);
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

        const Tree* initTree = firstTree;
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
          AlignedValuesContainer* sample = SiteContainerTools::bootstrapSites(*mSites.begin()->second);
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
        TreeTools::computeBootstrapValues(*firstTree, bsTrees);
        ApplicationTools::displayTaskDone();
        for (unsigned int i = 0; i < nbBS; i++)
        {
          delete bsTrees[i];
        }

        // Write resulting tree:
        vector<const Tree*> vcTree;
    
        for (it = mTree.begin(); it != mTree.end(); it++)
          vcTree.push_back(it->second);

        PhylogeneticsApplicationTools::writeTrees(vcTree, bppml.getParams());
      }

      delete alphabet;

      for (auto itc : mSites)
        delete itc.second;

      for (it = mTree.begin(); it != mTree.end(); it++)
        delete it->second;

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

