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
#include <Bpp/Numeric/DataTable.h>

// // From bpp-seq:
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/SiteTools.h>

// // From bpp-phyl:
#include <Bpp/Phyl/App/BppPhylogeneticsApplication.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Legacy/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Legacy/Likelihood/RASTools.h>
#include <Bpp/Phyl/Legacy/Likelihood/NNIHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Model/MixedTransitionModel.h>
#include <Bpp/Phyl/Io/Newick.h>

using namespace bpp;

/******************************************************************************/

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

  try
  {
    BppPhylogeneticsApplication bppml(args, argv, "bppml");

    if (args == 1)
    {
      bppml.help("bppml");
      return 0;
    }

    bppml.startTimer();

    map<string, string> unparsedParams;

    Context context;
    
    ///// Alphabet

    unique_ptr<Alphabet> alphabet(bppml.getAlphabet());

    /// GeneticCode
    
    unique_ptr<GeneticCode> gCode(bppml.getGeneticCode(alphabet.get()));

    ////// Get the map of the sequences

    map<size_t, AlignedValuesContainer*> mSites = bppml.getAlignmentsMap(alphabet.get());

    /////// Get the map of initial trees

    auto mpTree = bppml.getPhyloTreesMap(mSites, unparsedParams);

    // Try to write the current tree to file. This will be overwritten
    // by the optimized tree, but allow to check file existence before
    // running optimization!

    vector<const PhyloTree*> vcpTree;
    
    for (const auto& pTree : mpTree)
      vcpTree.push_back(pTree.second.get());
    
    PhylogeneticsApplicationTools::writeTrees(vcpTree, bppml.getParams(),"output.","",true,false,true);

    ApplicationTools::displayWarning("Reading trees for oldlik version : to be removed when not needed.");

    map<size_t, Tree*> mTree=PhylogeneticsApplicationToolsOld::getTrees(bppml.getParams(), mSites, unparsedParams);
    

    /////////////////
    // Computing stuff

    unique_ptr<DiscreteRatesAcrossSitesTreeLikelihood> tl_old;

    PhyloLikelihood* tl_new = 0;
    
    unique_ptr<SubstitutionProcessCollection> SPC;
    map<size_t, SequenceEvolution*> mSeqEvol;

    shared_ptr<PhyloLikelihoodContainer> mPhyl=0;
    
    bool checkTree    = ApplicationTools::getBooleanParameter("input.tree.check_root", bppml.getParams(), true, "", true, 2);
    bool optimizeTopo = ApplicationTools::getBooleanParameter("optimization.topology", bppml.getParams(), false, "", true, 1);
    unsigned int nbBS = ApplicationTools::getParameter<unsigned int>("bootstrap.number", bppml.getParams(), 0, "", true, 1);

    shared_ptr<BranchModel>    model;
    shared_ptr<TransitionModel>    tmodel; // for legacy 
    unique_ptr<SubstitutionModelSet> modelSet;
    shared_ptr<DiscreteDistribution> rDist;

    Tree* firstTree = mTree.begin()->second;

    /// Topology estimation

    ////////////
    // If optimize topology
    
    if (optimizeTopo || nbBS > 0)
    {
      string nhOpt = ApplicationTools::getStringParameter("nonhomogeneous", bppml.getParams(), "no", "", true, false);
      ApplicationTools::displayResult("Heterogeneous model", nhOpt);

      if (nhOpt != "no")
        throw Exception("Topology estimation with NH model not supported yet, sorry :(");

      tmodel=dynamic_pointer_cast<TransitionModel>(PhylogeneticsApplicationTools::getBranchModels(alphabet.get(), gCode.get(), mSites, bppml.getParams(), unparsedParams).begin()->second);

      if (tmodel->getName() != "RE08")
        for (auto  itc : mSites)
          SiteContainerTools::changeGapsToUnknownCharacters(*itc.second);

      if (tmodel->getNumberOfStates() >= 2 * tmodel->getAlphabet()->getSize())
      {
        // Markov-modulated Markov model!
        rDist=std::make_shared<ConstantRateDistribution>();
      }
      else
      {
        rDist=PhylogeneticsApplicationTools::getRateDistributions(bppml.getParams())[0];
      }
      if (dynamic_cast<MixedTransitionModel*>(tmodel.get()) == 0)
        tl_old.reset(new NNIHomogeneousTreeLikelihood(*firstTree, *mSites.begin()->second, tmodel.get(), rDist.get(), checkTree, true));
      else
        throw Exception("Topology estimation with Mixed model not supported yet, sorry :(");
    }

    /// Constant Topology

    else
    {
      // No RE08 model (yet)
      for (auto  itc : mSites)
        SiteContainerTools::changeGapsToUnknownCharacters(*itc.second);

      SPC.reset(bppml.getCollection(alphabet.get(), gCode.get(), mSites, mpTree, unparsedParams));

      mSeqEvol = bppml.getProcesses(*SPC, unparsedParams);

      mPhyl=bppml.getPhyloLikelihoods(context, mSeqEvol, *SPC, mSites);
      
      // retrieve Phylo 0, aka result phylolikelihood
      
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

      for (auto& it : mSeqEvol)
        delete it.second;
      
      for (size_t i = 0; i < pl.size(); ++i) {
        pnfile << pl[i].getName() << endl;
      }
      pnfile.close();
      cout << "BppML's done." << endl;
      exit(0);
    }

    //Output trees
    string treeWIdPath = ApplicationTools::getAFilePath("output.tree_ids.file", bppml.getParams(), false, false, "", true, "none", 1);
    if (treeWIdPath != "none")
    {
      bppml.getParams()["output_ids.tree.file"]=treeWIdPath;
      
      PhylogeneticsApplicationTools::writeTrees(*SPC, bppml.getParams(), "output_ids.", "", true, true, false, true);

      ApplicationTools::displayResult("Writing tagged tree to", treeWIdPath + "_...");
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
        if (AlphabetTools::isCodonAlphabet(alphabet.get()))
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

      tl_old.reset(dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood*>(
                     PhylogeneticsApplicationToolsOld::optimizeParameters(tl_old.get(), tl_old->getParameters(), bppml.getParams())));

      Tree* tree = new TreeTemplate<Node>(tl_old->getTree());
      PhylogeneticsApplicationToolsOld::writeTree(*tree, bppml.getParams());

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
          PhylogeneticsApplicationToolsOld::printParameters(modelSet.get(), out, 1, withAlias);
        }
        else
        {
          tmodel->matchParametersValues(tl_old->getParameters());
          PhylogeneticsApplicationTools::printParameters(tmodel.get(), out, 1);
        }
        out.endLine();
        (out << "# Rate distribution parameters:").endLine();
        rDist->matchParametersValues(tl_old->getParameters());
        PhylogeneticsApplicationTools::printParameters(rDist.get(), out);
      }

      // Getting posterior rate class distribution:
      unique_ptr<DiscreteDistribution> prDist(RASTools::getPosteriorRateDistribution(*tl_old));
      ApplicationTools::displayMessage("\nPosterior rate distribution for dataset:\n");
      if (ApplicationTools::message)
        prDist->print(*ApplicationTools::message);
      ApplicationTools::displayMessage("\n");

      // Write infos to file:
      string infosFile = ApplicationTools::getAFilePath("output.infos", bppml.getParams(), false, false);
      if (infosFile != "none")
      {
        ApplicationTools::displayResult("Alignment information logfile", infosFile);
        ofstream out(infosFile.c_str(), ios::out);

        // Get the rate class with maximum posterior probability:
        vector<size_t> classes = tl_old->getRateClassWithMaxPostProbPerSite();

        // Get the posterior rate, i.e. rate averaged over all posterior probabilities:
        Vdouble rates = tl_old->getPosteriorRatePerSite();

        vector<string> colNames;
        colNames.push_back("Sites");
        colNames.push_back("is.complete");
        colNames.push_back("is.constant");
        colNames.push_back("lnL");
        colNames.push_back("rc");
        colNames.push_back("pr");
        vector<string> row(6);
        unique_ptr<DataTable> infos(new DataTable(colNames));

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
      }
    }
    // New optimization
    else // tl_new!=0
    {
      //Check initial likelihood:
      
      bppml.fixLikelihood(alphabet.get(), gCode.get(), tl_new);

      // First `true` means that default is to optimize model parameters.
      if(ApplicationTools::getBooleanParameter("optimization.model_parameters", bppml.getParams(), true, "", true, 1))
        tl_new=PhylogeneticsApplicationTools::optimizeParameters(tl_new, tl_new->getParameters(), bppml.getParams());
      else
        tl_new=PhylogeneticsApplicationTools::optimizeParameters(tl_new, tl_new->getBranchLengthParameters(), bppml.getParams());

      SPC->matchParametersValues(tl_new->getParameters());
      
      PhylogeneticsApplicationTools::writeTrees(*SPC, bppml.getParams(), "output.", "", true, true, true);

      
      // Write parameters to screen:
      bppml.displayParameters(*tl_new);

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

        PhylogeneticsApplicationTools::printParameters(SPC.get(), out, 1, withAlias);

        for (map<size_t, SequenceEvolution*>::const_iterator it2=mSeqEvol.begin(); it2!=mSeqEvol.end(); it2++)
        {
          PhylogeneticsApplicationTools::printParameters(it2->second, out, it2->first);
          out.endLine();
        }

        PhylogeneticsApplicationTools::writeTrees(*SPC, bppml.getParams(), "output.", "",true,true,false,false);

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
          tl_old.reset(dynamic_cast<NNIHomogeneousTreeLikelihood*>(
                         PhylogeneticsApplicationToolsOld::optimizeParameters(tl_old.get(), tl_old->getParameters(), bppml.getParams(), "", true, false)));
          initTree = &tl_old->getTree();
        }

        string bsTreesPath = ApplicationTools::getAFilePath("bootstrap.output.file", bppml.getParams(), false, false);
        unique_ptr<ofstream> out;
        if (bsTreesPath != "none")
        {
          ApplicationTools::displayResult("Bootstrap trees stored in file", bsTreesPath);
          out.reset(new ofstream(bsTreesPath.c_str(), ios::out));
        }
        Newick newick;
        ParameterList paramsToIgnore = tl_old->getSubstitutionModelParameters();
        paramsToIgnore.addParameters(tl_old->getRateDistributionParameters());

        ApplicationTools::displayTask("Bootstrapping", true);
        vector<Tree*> bsTrees(nbBS);
        for (unsigned int i = 0; i < nbBS; i++)
        {
          ApplicationTools::displayGauge(i, nbBS - 1, '=');
          unique_ptr<AlignedValuesContainer> sample(SiteContainerTools::bootstrapSites(*mSites.begin()->second));
          if (!approx)
          {
            tmodel->setFreqFromData(*sample);
          }

          if (dynamic_cast<MixedTransitionModel*>(tmodel.get()) != NULL)
            throw Exception("Bootstrap estimation with Mixed model not supported yet, sorry :(");

          unique_ptr<NNIHomogeneousTreeLikelihood> tlRep(new NNIHomogeneousTreeLikelihood(*initTree, *sample, tmodel.get(), rDist.get(), true, false));
          tlRep->initialize();
          ParameterList parametersRep = tlRep->getParameters();
          if (approx)
          {
            parametersRep.deleteParameters(paramsToIgnore.getParameterNames());
          }
          tlRep.reset(dynamic_cast<NNIHomogeneousTreeLikelihood*>(
                        PhylogeneticsApplicationToolsOld::optimizeParameters(tlRep.get(), parametersRep, bppml.getParams(), "", true, false)));
          bsTrees[i] = new TreeTemplate<Node>(tlRep->getTree());
          if (out && i == 0) newick.writeTree(*bsTrees[i], bsTreesPath, true);
          if (out && i >  0) newick.writeTree(*bsTrees[i], bsTreesPath, false);
        }
        if (out)
          out->close();
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
    
        for (const auto& tree : mTree)
          vcTree.push_back(tree.second);

        PhylogeneticsApplicationToolsOld::writeTrees(vcTree, bppml.getParams());
      }

      for (auto& it : mSeqEvol)
        delete it.second;

      for (auto& itc : mSites)
        delete itc.second;

      for (const auto& tree : mTree)
        delete tree.second;

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

