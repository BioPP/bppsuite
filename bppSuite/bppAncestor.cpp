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

// From bpp-core:
#include <Bpp/Version.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/Numeric/DataTable.h>

// From bpp-seq:
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/Container/SequenceContainerTools.h>
#include <Bpp/Seq/Io/BppOAlignmentWriterFormat.h>

// From PhylLib:

// From Newlik:
#include <Bpp/Phyl/NewLikelihood/MarginalAncestralReconstruction.h>

#include <Bpp/Phyl/NewLikelihood/PhyloLikelihoods/OneProcessSequencePhyloLikelihood.h>
#include <Bpp/Phyl/NewLikelihood/PhyloLikelihoods/SetOfAbstractPhyloLikelihood.h>
#include <Bpp/Phyl/NewLikelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>

#include "bppTools.h"

using namespace bpp;

/******************************************************************************/

int main(int args, char ** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*     Bio++ Ancestral Sequence Reconstruction, version " << BPP_VERSION << "     *" << endl;
  cout << "* Authors: J. Dutheil                       Created on: 10/09/08 *" << endl;
  cout << "*          B. Boussau                       Last Modif: 25/09/14 *" << endl;
  cout << "*          L. Guéguen                       Last Modif: 22/12/14 *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  if (args == 1)
  {
    bppTools::help("bppAncestor");
    return 0;
  }
  
  try {

    
    BppApplication bppancestor(args, argv, "bppancestor");
    bppancestor.startTimer();

    map<string, string> allParams=bppancestor.getParams();
    map<string, string> unparsedparams;

    unique_ptr<Alphabet> alphabet(bppTools::getAlphabet(allParams));
    unique_ptr<GeneticCode> gCode(bppTools::getGeneticCode(allParams, alphabet.get()));

// Missing check
//  if (model->getName() != "RE08") SiteContainerTools::changeGapsToUnknownCharacters(*sites);

    // get the result phylo likelihood

    PhyloLikelihood* tl=bppTools::getResultPhyloLikelihood(allParams, alphabet.get(), gCode.get(), unparsedparams);
    
    bppTools::fixLikelihood(allParams, alphabet.get(), gCode.get(), tl);
    
    bppTools::displayParameters(*tl);

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

    string sequenceFilePath = ApplicationTools::getAFilePath("output.sequence.file", allParams, false, false, "", false, "none", 1);

    string sequenceFormat="";
  
    if (sequenceFilePath!="none")
      sequenceFormat   = ApplicationTools::getStringParameter("output.sequence.format", allParams, "Fasta", "", false, 1);

    string outputSitesFile, outputNodesFile;
  
    /// map of the Single Process 
    map<size_t, AbstractSingleDataPhyloLikelihood*> mSD;
        
    if (dynamic_cast<AbstractSingleDataPhyloLikelihood*>(tl)!=NULL)
      mSD[1]=dynamic_cast<AbstractSingleDataPhyloLikelihood*>(tl);
    else{
      SetOfAbstractPhyloLikelihood* sOAP=dynamic_cast<SetOfAbstractPhyloLikelihood*>(tl);
      if (sOAP)
      {
        const vector<size_t>& nSD=sOAP->getNumbersOfPhyloLikelihoods();
        
        for (size_t iSD=0; iSD< nSD.size(); iSD++){
          AbstractSingleDataPhyloLikelihood* pASDP=dynamic_cast<AbstractSingleDataPhyloLikelihood*>(sOAP->getAbstractPhyloLikelihood(nSD[iSD]));
          if (pASDP!=NULL)
            mSD[nSD[iSD]]=pASDP;
        }
      }
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
    
    for (map<size_t, AbstractSingleDataPhyloLikelihood*>::iterator itm=mSD.begin();itm!=mSD.end(); itm++)
    {
      SingleProcessPhyloLikelihood* sPP=dynamic_cast<SingleProcessPhyloLikelihood*>(itm->second);
      OneProcessSequencePhyloLikelihood* oPSP=dynamic_cast<OneProcessSequencePhyloLikelihood*>(itm->second);

      if ((sPP==NULL) && (oPSP==NULL))
      {
        ApplicationTools::displayWarning("Multi Process ancestral reconstruction not implemented");
        ApplicationTools::displayWarning("for phyloLikelihood " + TextTools::toString(itm->first));
        continue;
      }

      const AlignedValuesContainer* sites;
      if (sPP!=NULL)
        sites=sPP->getData();
      else
        sites=oPSP->getData();
        
      AbstractLikelihoodTreeCalculation* pDR;
      
      if (sPP!=NULL)
        pDR=dynamic_cast<AbstractLikelihoodTreeCalculation*>(sPP->getLikelihoodCalculation());
      else
        pDR=dynamic_cast<AbstractLikelihoodTreeCalculation*>(oPSP->getLikelihoodCalculation());

      if (probMethod)
      {    
        AncestralStateReconstruction *asr = new MarginalAncestralReconstruction(pDR);

        size_t nbStates=sPP->getNumberOfStates();

        /////////////////////////////////////
        // Write infos to file:
        
        if (outputSitesFile != "none")
        {
          ApplicationTools::displayResult("Phylo " + TextTools::toString(itm->first) + " : Output file for sites", outputSitesFile + "_" + TextTools::toString(itm->first));
          string outF=outputSitesFile + "_" + TextTools::toString(itm->first);
          ofstream out(outF.c_str(), ios::out);
          PhyloTree ttree(sPP->getTree());
          vector<shared_ptr<PhyloNode> > nodes = ttree.getAllNodes();
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
            shared_ptr<PhyloNode> node = nodes[i];
            colNames.push_back("max." + TextTools::toString(ttree.getNodeIndex(node)));
            if (probs) {
              probabilities[i] = new VVdouble();
              
              //The cast will have to be updated when more probabilistic method will be available:
              sequences[i] = dynamic_cast<MarginalAncestralReconstruction *>(asr)->getAncestralSequenceForNode(ttree.getNodeIndex(node), probabilities[i], false);

              for (unsigned int j = 0; j < nbStates; j++) {
                colNames.push_back("prob." + TextTools::toString(ttree.getNodeIndex(node)) + "." + alphabet->intToChar((int)j));
              }
            }
            else
              sequences[i] = asr->getAncestralSequenceForNode(ttree.getNodeIndex(node));
          }

          //Now fill the table:
          vector<string> row(colNames.size());
          DataTable* infos = new DataTable(colNames);

          sPP->computeLikelihood();

          for (size_t i = 0; i < sites->getNumberOfSites(); i++)
          {
            double lnL = sPP->getLogLikelihoodForASite(i);
            const CruxSymbolListSite& currentSite = sites->getSymbolListSite(i);
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
        
        ////////////////////////////////////////////////
        /// Ancestral sequences
        
        if (sequenceFilePath!="none")
        {
          SiteContainer* asSites = 0;
          
          if (sample)
          {
            asSites = new VectorSiteContainer(alphabet.get());
            
            for (unsigned int i = 0; i < nbSamples; i++)
            {
              ApplicationTools::displayGauge(i, nbSamples-1, '=');
              SequenceContainer *sampleSites = dynamic_cast<MarginalAncestralReconstruction *>(asr)->getAncestralSequences(true);
              vector<string> names = sampleSites->getSequencesNames();
              for (unsigned int j = 0; j < names.size(); j++){
                names[j] += "_" + TextTools::toString(i+1);
              }
              
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
          {
            const SiteContainer* vSC=dynamic_cast<const SiteContainer*>(sites);
            if (!vSC)
              ApplicationTools::displayMessage("Output sequences not possible with probabilistic sequences.");
            else
              SequenceContainerTools::append(*asSites, *vSC);
          }
          
          //Write output:
          BppOAlignmentWriterFormat bppoWriter(1);
          unique_ptr<OAlignment> oAln(bppoWriter.read(sequenceFormat));
          ApplicationTools::displayResult("Output alignment file ", sequenceFilePath + "_" + TextTools::toString(itm->first));
          ApplicationTools::displayResult("Output alignment format ", oAln->getFormatName());
          
          // Write sequences:
          oAln->writeAlignment(sequenceFilePath + "_" + TextTools::toString(itm->first), *asSites, true);
          
          delete asSites;
          delete asr;
        }
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
          colNames.push_back("exp" + sPP->getData()->getAlphabet()->intToChar((int)i));
        for (size_t i = 0; i < sPP->getNumberOfStates(); i++)
          colNames.push_back("eb" + sPP->getData()->getAlphabet()->intToChar((int)i));
      
        //Now fill the table:
        vector<string> row(colNames.size());
        DataTable* infos = new DataTable(colNames);
      
        for (map<int, vector<double> >::iterator itf = frequencies.begin(); itf != frequencies.end(); itf++)
        {
          row[0] = TextTools::toString(itf->first);
          Vdouble ebFreqs = pDR->getLikelihoodData().getPosteriorStateFrequencies(itf->first);
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
  
    bppancestor.done();
  
  }
  catch (exception & e)
  {
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}

