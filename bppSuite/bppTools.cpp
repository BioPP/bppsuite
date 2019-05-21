//
// File: bppTools.cpp
// Created by: Laurent Guéguen
// Created on: mardi 28 novembre 2017, à 09h 05
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

#include "bppTools.h"

#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>

using namespace std;
using namespace bpp;

/******************************************************************************/

void bppTools::help(const string& program)
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
  (*ApplicationTools::message << program << " parameter1_name=parameter1_value parameter2_name=parameter2_value").endLine();
  (*ApplicationTools::message << "      ... param=option_file").endLine();
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "  Refer to the Bio++ Program Suite Manual for a list of available options.").endLine();
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
}

Alphabet* bppTools::getAlphabet(const map<string, string>& params)
{
  return SequenceApplicationTools::getAlphabet(params, "", false);
}
  
GeneticCode* bppTools::getGeneticCode(const map<string, string>& params,
                                      const Alphabet* alphabet)
{
  const CodonAlphabet* codonAlphabet = dynamic_cast<const CodonAlphabet*>(alphabet);
  if (codonAlphabet) {
    string codeDesc = ApplicationTools::getStringParameter("genetic_code", params, "Standard", "", true, true);
    ApplicationTools::displayResult("Genetic Code", codeDesc);
    
    return SequenceApplicationTools::getGeneticCode(codonAlphabet->getNucleicAlphabet()->clone(), codeDesc);
  }
  else
    return 0;
}
  
map<size_t, AlignedValuesContainer*> bppTools::getAlignmentsMap(const map<string, string>& params,
                                                                const Alphabet* alphabet,
                                                                bool changeGapsToUnknownCharacters,
                                                                bool optionalData)
{
  map<size_t, AlignedValuesContainer*> mSites = SequenceApplicationTools::getAlignedContainers(alphabet, params);
  
  if (!optionalData && mSites.size() == 0)
    throw Exception("Missing data input.sequence.file option");

  if (changeGapsToUnknownCharacters)
    for (auto itc : mSites)
      SiteContainerTools::changeGapsToUnknownCharacters(*itc.second);

  return mSites;
}
  
map<size_t, PhyloTree*> bppTools::getPhyloTreesMap(const map<string, string>& params,
                                                   const map<size_t, AlignedValuesContainer*>& mSites,
                                                   map<string, string>& unparsedparams)
{
  map<size_t, PhyloTree*> mpTree=PhylogeneticsApplicationTools::getPhyloTrees(params, mSites, unparsedparams);

  // Scaling of trees:
  double scale = ApplicationTools::getDoubleParameter("input.tree.scale", params, 1, "", false, false);

  if (scale != 1) {
    ApplicationTools::displayResult("Trees are scaled by", scale);
    
    for (auto it : mpTree) {
      it.second -> scaleTree(scale);
    }
  }

  return mpTree;
}

SubstitutionProcessCollection* bppTools::getCollection(const map<string, string>& params,
                                                       const Alphabet* alphabet,
                                                       const GeneticCode* gCode,
                                                       const map<size_t, AlignedValuesContainer*>& mSites,
                                                       map<string, string>& unparsedparams)
{
  map<size_t, PhyloTree*> mpTree = getPhyloTreesMap(params, mSites, unparsedparams);

  SubstitutionProcessCollection* SPC= getCollection(params, alphabet, gCode, mSites, mpTree, unparsedparams);
      
  return SPC;
}

SubstitutionProcessCollection* bppTools::getCollection(const map<string, string>& params,
                                                       const Alphabet* alphabet,
                                                       const GeneticCode* gCode,
                                                       const map<size_t, AlignedValuesContainer*>& mSites,
                                                       const map<size_t, PhyloTree*>& mpTree,
                                                       map<string, string>& unparsedparams)
{
  map<size_t, DiscreteDistribution*> mDist = PhylogeneticsApplicationTools::getRateDistributions(params);
  
  map<size_t, TransitionModel*> mMod = PhylogeneticsApplicationTools::getTransitionModels(alphabet, gCode, mSites, params, unparsedparams);
  
  map<size_t, FrequenciesSet*> mRootFreq = PhylogeneticsApplicationTools::getRootFrequenciesSets(alphabet, gCode, mSites, params, unparsedparams);
      
  SubstitutionProcessCollection* SPC=PhylogeneticsApplicationTools::getSubstitutionProcessCollection(alphabet, gCode, mpTree, mMod, mRootFreq, mDist, params, unparsedparams);
      
  return SPC;
}


map<size_t, SequenceEvolution*> bppTools::getProcesses(const map<string, string>& params,
                                                       SubstitutionProcessCollection& collection,
                                                       map<string, string>& unparsedparams)
{
  return PhylogeneticsApplicationTools::getSequenceEvolutions(collection, params, unparsedparams);
}


PhyloLikelihoodContainer* bppTools::getPhyloLikelihoods(const map<string, string>& params,
                                                        dataflow::Context& context,
                                                        map<size_t, SequenceEvolution*> mSeqEvol, 
                                                        SubstitutionProcessCollection& collection,
                                                        const map<size_t, AlignedValuesContainer*>& mSites)
{
  return PhylogeneticsApplicationTools::getPhyloLikelihoodContainer(context, collection, mSeqEvol, mSites, params);
}


PhyloLikelihood* bppTools::getResultPhyloLikelihood(const std::map<std::string, std::string>& params,
                                                    dataflow::Context& context,
                                                    const Alphabet* alphabet,
                                                    const GeneticCode* gCode,
                                                    std::map<std::string, std::string>& unparsedparams)
{
  map<size_t, AlignedValuesContainer*> mSites = getAlignmentsMap(params, alphabet, true);

  map<size_t, PhyloTree*> mTrees = bppTools::getPhyloTreesMap(params, mSites, unparsedparams);

  SubstitutionProcessCollection* SPC = bppTools::getCollection(params, alphabet, gCode, mSites, mTrees, unparsedparams);

  map<size_t, SequenceEvolution*> mSeqEvol = PhylogeneticsApplicationTools::getSequenceEvolutions(*SPC, params, unparsedparams);

  PhyloLikelihoodContainer* mPhyl=PhylogeneticsApplicationTools::getPhyloLikelihoodContainer(context, *SPC, mSeqEvol, mSites, params);

  if (!mPhyl->hasPhyloLikelihood(0))
    throw Exception("Missing phyloLikelihoods.");

  return (*mPhyl)[0];
}



void bppTools::fixLikelihood(const map<string, string>& params,
                             const Alphabet* alphabet,
                             const GeneticCode* gCode,
                             PhyloLikelihood* phylolik)
{
  double logL = phylolik->getValue();

  if (std::isinf(logL))
  {
    // This may be due to null branch lengths, leading to null likelihood!
    ApplicationTools::displayWarning("!!! Warning!!! Initial likelihood is zero.");
    ApplicationTools::displayWarning("!!! This may be due to branch length == 0.");
    ApplicationTools::displayWarning("!!! All null branch lengths will be set to 0.000001.");
    ParameterList pl = phylolik->getBranchLengthParameters();
    for (unsigned int i = 0; i < pl.size(); i++)
    {
      if (pl[i].getValue() < 0.000001) pl[i].setValue(0.001);
    }
    phylolik->matchParametersValues(pl);
    logL = phylolik->getValue();
  }
  
  ApplicationTools::displayMessage("");
  ApplicationTools::displayResult("Initial log likelihood", TextTools::toString(-logL, 15));
  if (std::isinf(logL))
  {
    ApplicationTools::displayError("!!! Unexpected initial likelihood == 0.");

    map<size_t, AbstractSingleDataPhyloLikelihood*> mSD;

    if (dynamic_cast<AbstractSingleDataPhyloLikelihood*>(phylolik)!=NULL)
      mSD[1]=dynamic_cast<AbstractSingleDataPhyloLikelihood*>(phylolik);
    else{
      SetOfAbstractPhyloLikelihood* sOAP=dynamic_cast<SetOfAbstractPhyloLikelihood*>(phylolik);
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

        if (AlphabetTools::isCodonAlphabet(alphabet))
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

        // Then remove saturated positions
            
        bool removeSaturated = ApplicationTools::getBooleanParameter("input.sequence.remove_saturated_sites", params, false, "", true, false);

        if (!removeSaturated) {
          ApplicationTools::displayError("!!! Looking at each site:");
          for (unsigned int i = 0; i < vData->getNumberOfSites(); i++) {
            (*ApplicationTools::error << "Site " << vData->getSymbolListSite(i).getPosition() << "\tlog likelihood = " << sDP->getLogLikelihoodForASite(i)).endLine();
          }
          ApplicationTools::displayError("!!! 0 values (inf in log) may be due to computer overflow, particularily if datasets are big (>~500 sequences).");
          ApplicationTools::displayError("!!! You may want to try input.sequence.remove_saturated_sites = yes to ignore positions with likelihood 0.");
          exit(1);
        }
        else
        {
          ApplicationTools::displayBooleanResult("Saturated site removal enabled", true);
          for (size_t i = vData->getNumberOfSites(); i > 0; --i) {
            if (std::isinf(sDP->getLogLikelihoodForASite(i - 1))) {
              ApplicationTools::displayResult("Ignore saturated site", vData->getSymbolListSite(i - 1).getPosition());
              vData->deleteSites(i - 1, i);
            }
          }
          ApplicationTools::displayResult("Number of sites retained", vData->getNumberOfSites());

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
}


void bppTools::displayParameters(const PhyloLikelihood& tl, bool displaylL)
{
  // Write parameters to screen:
  if (displaylL)
    ApplicationTools::displayResult("Log likelihood", TextTools::toString(-tl.getValue(), 15));

  if (tl.getNumberOfParameters()-tl.getBranchLengthParameters().size()>=30)
    ApplicationTools::displayMessage("Too many parameters for screen output!");
  else
  {
    ParameterList parameters = tl.getParameters();
    parameters.deleteParameters(tl.getBranchLengthParameters().getParameterNames(),false);
    for (unsigned int i = 0; i < parameters.size(); i++)
      ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
  }
}

