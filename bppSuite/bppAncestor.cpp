//
// File: bppAncestor.cpp
// Created by: Julien Dutheil
// Created on: Sep Wed 10 14:14 2008
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

#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/KeyvalTools.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/AutoParameter.h>

// From SeqLib:
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SequenceContainerTools.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/Io/BppOAlignmentWriterFormat.h>

// From PhylLib:
#include <Bpp/Phyl/Tree/Tree.h>
#include <Bpp/Phyl/Likelihood.all>
#include <Bpp/Phyl/PatternTools.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Model/MarkovModulatedSubstitutionModel.h>
#include <Bpp/Phyl/Model/SubstitutionModelSet.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Io/Newick.h>

// From Newlik:
#include <Bpp/Phyl/NewLikelihood/DoubleRecursiveTreeLikelihoodCalculation.h>
#include <Bpp/Phyl/NewLikelihood/SingleProcessPhyloLikelihood.h>
#include <Bpp/Phyl/NewLikelihood/SingleDataPhyloLikelihood.h>
#include <Bpp/Phyl/NewLikelihood/SumOfDataPhyloLikelihood.h>
#include <Bpp/Phyl/NewLikelihood/SubstitutionProcessCollection.h>
#include <Bpp/Phyl/NewLikelihood/MarginalAncestralReconstruction.h>

using namespace bpp;
using namespace newlik;

/******************************************************************************/

void help()
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
  (*ApplicationTools::message << "bppancestor parameter1_name=parameter1_value ").endLine();
  (*ApplicationTools::message << "      parameter2_name=parameter2_value ... param=option_file").endLine();
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "  Refer to the Bio++ Program Suite Manual for a list of available options.").endLine();
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
}

int main(int args, char ** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*     Bio++ Ancestral Sequence Reconstruction, version 2.2.0     *" << endl;
  cout << "* Authors: J. Dutheil                       Created on: 10/09/08 *" << endl;
  cout << "*          B. Boussau                       Last Modif: 25/09/14 *" << endl;
  cout << "*          L. Guéguen                       Last Modif: 22/12/14 *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  if (args == 1)
  {
    help();
    return 0;
  }
  
  try {

    
    BppApplication bppancestor(args, argv, "BppAncestor");
    bppancestor.startTimer();

    map<string, string> allParams=bppancestor.getParams();
    
    Alphabet* alphabet = SequenceApplicationTools::getAlphabet(allParams, "", false);
    auto_ptr<GeneticCode> gCode;
    CodonAlphabet* codonAlphabet = dynamic_cast<CodonAlphabet*>(alphabet);
    if (codonAlphabet) {
      string codeDesc = ApplicationTools::getStringParameter("genetic_code", allParams, "Standard", "", true, true);
      ApplicationTools::displayResult("Genetic Code", codeDesc);
      
      gCode.reset(SequenceApplicationTools::getGeneticCode(codonAlphabet->getNucleicAlphabet(), codeDesc));
    }

  
    ////// Get the map of the sequences 

    map<size_t, SiteContainer*> mSites = SequenceApplicationTools::getSiteContainers(alphabet, allParams);

    if (mSites.size() == 0)
      throw Exception("Missing data input.sequence.file option");

    /////// Get the map of initial trees
    
    map<size_t, Tree*> mTree=PhylogeneticsApplicationTools::getTrees(allParams, mSites);

    // Try to write the current tree to file. This will be overwritten
    // by the optimized tree, but allow to check file existence before
    // running optimization!

    vector<const Tree*> vcTree;
    
    for (map<size_t, Tree*>::const_iterator it = mTree.begin(); it != mTree.end(); it++)
      vcTree.push_back(it->second);
    
    PhylogeneticsApplicationTools::writeTrees(vcTree, allParams);

    
    bool computeLikelihood = ApplicationTools::getBooleanParameter("compute.likelihood", allParams, true, "", false, 1);
    if (!computeLikelihood)
    {
      delete alphabet;
      
      for (map<size_t, SiteContainer*>::iterator itc=mSites.begin(); itc != mSites.end(); itc++)
        delete itc->second;
      
      for (map<size_t, Tree*>::const_iterator it = mTree.begin(); it != mTree.end(); it++)
        delete it->second;
      cout << "Bppancestor's done. Bye." << endl;
      return 0;
    }

    /////////////////
    // Computing stuff
    
    PhyloLikelihood* tl = 0;
    SubstitutionProcessCollection* SPC = 0;
  
    string collection = ApplicationTools::getStringParameter("collection", allParams, "", "", true, 1);

    map<string, string> unparsedparams;

    map<size_t, DiscreteDistribution*> mDist = PhylogeneticsApplicationTools::getRateDistributions(allParams);

    map<size_t, SubstitutionModel*> mMod = PhylogeneticsApplicationTools::getSubstitutionModels(alphabet, gCode.get(), mSites, allParams, unparsedparams);

    map<size_t, FrequenciesSet*> mRootFreq = PhylogeneticsApplicationTools::getRootFrequenciesSets(alphabet, gCode.get(), mSites, allParams, unparsedparams);

    SPC=PhylogeneticsApplicationTools::getSubstitutionProcessCollection(alphabet, gCode.get(), mTree, mMod, mRootFreq, mDist, allParams, unparsedparams);

  
// Missing check
//  if (model->getName() != "RE08") SiteContainerTools::changeGapsToUnknownCharacters(*sites);

    for (map<size_t, SiteContainer*>::iterator itc=mSites.begin(); itc != mSites.end(); itc++)
      SiteContainerTools::changeGapsToUnknownCharacters(*itc->second);

    // Change all phylolikelihoods recursion parameters to "Double"

    vector<string> phylosName=ApplicationTools::matchingParameters("phylo*", allParams);
    vector<size_t> phylosNum;
    for (size_t i=0; i< phylosName.size(); i++)
    {
      size_t poseq=phylosName[i].find("=");
      phylosNum.push_back((size_t)TextTools::toInt(phylosName[i].substr(5,poseq-5)));
    }

    map<size_t, PhyloLikelihood*> mPhylo;
  
    map<string, string> mrec;
    mrec["recursion"]="double";

    for (size_t mPi=0; mPi< phylosNum.size(); mPi++)
    {
      string phyloName = "";

      string phyloDesc = ApplicationTools::getStringParameter("phylo", allParams, "Single", TextTools::toString(phylosNum[mPi]), 0);

      string newPhyloDesc=KeyvalTools::changeKeyvals(phyloDesc, mrec);
      
      allParams["phylo"+TextTools::toString(phylosNum[mPi])]=newPhyloDesc;
    }
    
    // get the Double recursive phylo likelihoods
  
    map<size_t, PhyloLikelihood*> mPhyl=PhylogeneticsApplicationTools::getPhyloLikelihoods(*SPC, mSites, allParams, unparsedparams);

    map<size_t, SingleDataPhyloLikelihood*> mPhyl2;
      
    for (map<size_t, PhyloLikelihood*>::const_iterator itm=mPhyl.begin(); itm != mPhyl.end(); itm++)
      if (dynamic_cast<SingleDataPhyloLikelihood*>(itm->second))
        mPhyl2[itm->first]=dynamic_cast<SingleDataPhyloLikelihood*>(itm->second);

    if (mPhyl2.size()==1)
      tl=mPhyl2.begin()->second;
    else
      tl=new SumOfDataPhyloLikelihood(mPhyl2);      

    //////////////////////////////////////////////
    /// Infinite likelihood
  
    double logL = tl->getValue();
    if (isinf(logL))
    {
      // This may be due to null branch lengths, leading to null likelihood!
      ApplicationTools::displayWarning("!!! Warning!!! Initial likelihood is zero.");
      ApplicationTools::displayWarning("!!! This may be due to branch length == 0.");
      ApplicationTools::displayWarning("!!! All null branch lengths will be set to 0.001.");
      ParameterList pl = tl->getBranchLengthParameters();
      for (unsigned int i = 0; i < pl.size(); i++)
      {
        if (pl[i].getValue() < 0.001)
          pl[i].setValue(0.001);
      }
      tl->matchParametersValues(pl);
      logL = tl->getValue();
    }
  
    ApplicationTools::displayResult("Initial log likelihood", TextTools::toString(-logL, 15));
    if (isinf(logL))
    {
      ApplicationTools::displayError("!!! Unexpected initial likelihood == 0.");

      map<size_t, SingleDataPhyloLikelihood*> mSD;
        
      if (dynamic_cast<SingleDataPhyloLikelihood*>(tl)!=NULL)
        mSD[1]=dynamic_cast<SingleDataPhyloLikelihood*>(tl);
      else{
        MultiDataPhyloLikelihood* mDP=dynamic_cast<MultiDataPhyloLikelihood*>(tl);
        vector<size_t> nSD=mDP->getNumbersOfSingleDataPhyloLikelihoods();
          
        for (size_t iSD=0; iSD< nSD.size(); iSD++)
          mSD[nSD[iSD]]=mDP->getSingleDataPhylolikelihood(nSD[iSD]);
      }


      for (map<size_t, SingleDataPhyloLikelihood*>::iterator itm=mSD.begin();itm!=mSD.end(); itm++)
      {
        ApplicationTools::displayWarning("Checking for phylolikelihood " + TextTools::toString(itm->first));
          
        if (isinf(itm->second->getValue()))
        {
          SingleDataPhyloLikelihood* sDP=itm->second;
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
            
          bool removeSaturated = ApplicationTools::getBooleanParameter("input.sequence.remove_saturated_sites", allParams, false, "", true, false);
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
            ApplicationTools::displayResult("Number of sites retained", mSites.begin()->second->getNumberOfSites());

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

    // Write parameters to screen:
    //    ApplicationTools::displayResult("Log likelihood", TextTools::toString(-tl->getValue(), 15));
    ParameterList parameters = tl->getParameters();
    parameters.deleteParameters(tl->getBranchLengthParameters().getParameterNames(),false);
              
    for (unsigned int i = 0; i < parameters.size(); i++)
      ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));


    // Reconstruct ancestral sequences:
    string reconstruction = ApplicationTools::getStringParameter("asr.method", bppancestor.getParams(), "marginal", "", true, false);
    ApplicationTools::displayResult("Ancestral state reconstruction method", reconstruction);
    bool probs = false;

  
    bool probMethod = false;
    if (reconstruction == "none")
    {
      //do nothing
    } else if (reconstruction == "marginal")
    {
      probMethod = true;
    } else
      throw Exception("Unknown ancestral state reconstruction method: " + reconstruction);

    string sequenceFilePath = ApplicationTools::getAFilePath("output.sequence.file", allParams, true, false, "", false, "none", 1);

    string sequenceFormat="";
  
    if (sequenceFilePath!="none")
      sequenceFormat   = ApplicationTools::getStringParameter("output.sequence.format", allParams, "Fasta", "", false, 1);

    string outputSitesFile, outputNodesFile;
  
    /// map of the Single Process 
    map<size_t, SingleDataPhyloLikelihood*> mSD;
        
    if (dynamic_cast<SingleDataPhyloLikelihood*>(tl)!=NULL)
      mSD[1]=dynamic_cast<SingleDataPhyloLikelihood*>(tl);
    else{
      MultiDataPhyloLikelihood* mDP=dynamic_cast<MultiDataPhyloLikelihood*>(tl);
      vector<size_t> nSD=mDP->getNumbersOfSingleDataPhyloLikelihoods();
      
      for (size_t iSD=0; iSD< nSD.size(); iSD++)
        mSD[nSD[iSD]]=mDP->getSingleDataPhylolikelihood(nSD[iSD]);
    }

    bool sample=false;
    unsigned int nbSamples=0;
    bool addSitesExtant=false;
    bool addNodesExtant = false;

    // ASR
  
    if (probMethod)
    {
      probs = ApplicationTools::getBooleanParameter("asr.probabilities", allParams, false, "", true, false);
      ApplicationTools::displayResult("Output probabilities", probs ? "yes" : "no");

      outputSitesFile = ApplicationTools::getAFilePath("output.sites.file", allParams, false, false);

      sample = ApplicationTools::getBooleanParameter("asr.sample", allParams, false, "", true, false);

      ApplicationTools::displayResult("Sample from posterior distribution", sample ? "yes" : "no");

      if (sample)
        nbSamples = ApplicationTools::getParameter<unsigned int>("asr.sample.number", allParams, 1, "", true, false);

      addSitesExtant = ApplicationTools::getBooleanParameter("asr.add_extant", allParams, false, "", true, false);
    }

    // Nodes Frequencies
  
    outputNodesFile = ApplicationTools::getAFilePath("output.nodes.file", allParams, false, false);

    if (outputNodesFile!="none")
      addNodesExtant = ApplicationTools::getBooleanParameter("output.nodes.add_extant", allParams, false, "", true, false);

  
    /// per Single Process 
    
    for (map<size_t, SingleDataPhyloLikelihood*>::iterator itm=mSD.begin();itm!=mSD.end(); itm++)
    {
      
      SingleProcessPhyloLikelihood* sPP=dynamic_cast<SingleProcessPhyloLikelihood*>(itm->second);

      if (sPP==NULL){
        ApplicationTools::displayWarning("Multi Process ancestral reconstruction not implemented");
        ApplicationTools::displayWarning("for phylolikelihood " + TextTools::toString(itm->first));
        continue;
      }

      const SiteContainer* sites = sPP->getData();

      DoubleRecursiveTreeLikelihoodCalculation* pDR=dynamic_cast<DoubleRecursiveTreeLikelihoodCalculation*>(sPP->getLikelihoodCalculation());
    
      if (!pDR)
      {
        ApplicationTools::displayWarning("Not double recursive calculation for phylo likelihood " + TextTools::toString(itm->first));
        continue;
      }

      if (probMethod)
      {    
        AncestralStateReconstruction *asr = new MarginalAncestralReconstruction(pDR);

        size_t nbStates=sPP->getNumberOfStates();
      
        // Write infos to file:
        if (outputSitesFile != "none")
        {
          ApplicationTools::displayResult("Phylo " + TextTools::toString(itm->first) + " : Output file for sites", outputSitesFile + "_" + TextTools::toString(itm->first));
          string outF=outputSitesFile + "_" + TextTools::toString(itm->first);
          ofstream out(outF.c_str(), ios::out);
          TreeTemplate<Node> ttree(sPP->getTree());
          vector<Node *> nodes = ttree.getInnerNodes();
          size_t nbNodes = nodes.size();
    
          // Get the class with maximum posterior probability:
          vector<size_t> classes = sPP->getClassWithMaxPostProbOfEachSite();
          // Get the posterior rate, i.e. rate averaged over all posterior probabilities:
          Vdouble rates = sPP->getPosteriorRateOfEachSite();
          // Get the ancestral sequences:
          vector<Sequence*> sequences(nbNodes);
          vector<VVdouble*> probabilities(nbNodes);
        
          vector<string> colNames;
          colNames.push_back("Sites");
          colNames.push_back("is.complete");
          colNames.push_back("is.constant");
          colNames.push_back("lnL");
          colNames.push_back("rc");
          colNames.push_back("pr");
          for (size_t i = 0; i < nbNodes; i++) {
            Node *node = nodes[i];
            colNames.push_back("max." + TextTools::toString(node->getId()));
            if (probs) {
              probabilities[i] = new VVdouble();
              //The cast will have to be updated when more probabilistic method will be available:
              sequences[i] = dynamic_cast<MarginalAncestralReconstruction *>(asr)->getAncestralSequenceForNode(node->getId(), probabilities[i], false);
            
              for (unsigned int j = 0; j < nbStates; j++) {
                colNames.push_back("prob." + TextTools::toString(node->getId()) + "." + alphabet->intToChar((int)j));
              }
            }
            else
            {
              if (node->isLeaf()) {
              
              } else {
                sequences[i] = asr->getAncestralSequenceForNode(node->getId());
              }
            }
          }
        
          //Now fill the table:
          vector<string> row(colNames.size());
          DataTable* infos = new DataTable(colNames);

          for (size_t i = 0; i < sites->getNumberOfSites(); i++)
          {
            double lnL = sPP->getLogLikelihoodForASite(i);
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

            unsigned int k = 6;
            for (unsigned int j = 0; j < nbNodes; j++) {
              row[k] = sequences[j]->getChar(i);
              k++;
              if (probs) {
                for (unsigned int l = 0; l < nbStates; l++) {
                  row[k] = TextTools::toString((*probabilities[j])[i][l]);
                  k++;
                }
              }
            }
          
            infos->addRow(row);
          }
        
          DataTable::write(*infos, out, "\t");
        
          delete infos;
        }

      
        /// Ancestral sequences
      
        SiteContainer* asSites = 0;
      
        if (sample)
        {
          asSites = new AlignedSequenceContainer(alphabet);

          for (unsigned int i = 0; i < nbSamples; i++)
          {
            ApplicationTools::displayGauge(i, nbSamples-1, '=');
            SequenceContainer *sampleSites = dynamic_cast<MarginalAncestralReconstruction *>(asr)->getAncestralSequences(true);
            vector<string> names = sampleSites->getSequencesNames();
            for (unsigned int j = 0; j < names.size(); j++)
              names[j] += "_" + TextTools::toString(i+1);
            sampleSites->setSequencesNames(names, false);
            SequenceContainerTools::append(*asSites, *sampleSites);
            delete sampleSites;
          }
          ApplicationTools::message->endLine();
        }
        else
          asSites = asr->getAncestralSequences();
      
        //Add existing sequence to output?
        if (addSitesExtant) 
          SequenceContainerTools::append(*asSites, *sites);


        //Write output:
        BppOAlignmentWriterFormat bppoWriter(1);
        auto_ptr<OAlignment> oAln(bppoWriter.read(sequenceFormat));
        ApplicationTools::displayResult("Output alignment file ", sequenceFilePath + "_" + TextTools::toString(itm->first));
        ApplicationTools::displayResult("Output alignment format ", oAln->getFormatName());

        // Write sequences:
        oAln->writeAlignment(sequenceFilePath + "_" + TextTools::toString(itm->first), *sites, true);
    
        delete asSites;
        delete asr;
      }

      /// end of ancestral reconstruction
    
      if (outputNodesFile != "none")
      {
        ApplicationTools::displayResult("Output file for nodes", outputNodesFile + "_" + TextTools::toString(itm->first));
        string outF=outputNodesFile  + "_" + TextTools::toString(itm->first);
        
        ofstream out(outF.c_str(), ios::out);
        map<int, vector<double> > frequencies;
        pDR->getAncestralFrequencies(frequencies, addNodesExtant);
      
        vector<string> colNames;
        colNames.push_back("Nodes");
        for (size_t i = 0; i < sPP->getNumberOfStates(); i++)
          colNames.push_back("exp" + sPP->getAlphabet()->intToChar(i));
        for (size_t i = 0; i < sPP->getNumberOfStates(); i++)
          colNames.push_back("eb" + sPP->getAlphabet()->intToChar(i));
      
        //Now fill the table:
        vector<string> row(colNames.size());
        DataTable* infos = new DataTable(colNames);
      
        for (map<int, vector<double> >::iterator itf = frequencies.begin(); itf != frequencies.end(); itf++)
        {
          row[0] = TextTools::toString(itf->first);
          Vdouble ebFreqs = pDR->getPosteriorStateFrequencies(itf->first);
          for (size_t i = 0; i < sPP->getNumberOfStates(); i++)
          {
            row[i + 1] = TextTools::toString(itf->second[i]);
          }
          for (size_t i = 0; i < sPP->getNumberOfStates(); i++)
          {
            row[i + sPP->getNumberOfStates() + 1] = TextTools::toString(ebFreqs[i]);
          }
          infos->addRow(row);
        }
       
        DataTable::write(*infos, out, "\t");
       
        delete infos;
      }
    }
  
    delete alphabet;
    bppancestor.done();
  
  }
  catch (exception & e)
  {
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}

