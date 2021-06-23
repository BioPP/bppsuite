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
#include <Bpp/Numeric/DataTable.h>

// From bpp-seq:
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/Container/SequenceContainerTools.h>
#include <Bpp/Seq/Io/BppOAlignmentWriterFormat.h>

// From bpp-phyl:
#include <Bpp/Phyl/App/BppPhylogeneticsApplication.h>
#include <Bpp/Phyl/NewLikelihood/MarginalAncestralReconstruction.h>
#include <Bpp/Phyl/NewLikelihood/PhyloLikelihoods/OneProcessSequencePhyloLikelihood.h>
#include <Bpp/Phyl/NewLikelihood/PhyloLikelihoods/SetOfAbstractPhyloLikelihood.h>
#include <Bpp/Phyl/NewLikelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>

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

  try {

    BppPhylogeneticsApplication bppancestor(args, argv, "bppancestor");

    if (args == 1)
    {
      bppancestor.help("bppAncestor");
      return 0;
    }
  
    bppancestor.startTimer();
    
    Context context;
    map<string, string> allParams=bppancestor.getParams();
    map<string, string> unparsedParams;

    unique_ptr<Alphabet> alphabet(bppancestor.getAlphabet());
    unique_ptr<GeneticCode> gCode(bppancestor.getGeneticCode(alphabet.get()));

    // Missing check
    //  if (model->getName() != "RE08") SiteContainerTools::changeGapsToUnknownCharacters(*sites);

    // get the result phylo likelihood
    map<size_t, AlignedValuesContainer*> mSites = bppancestor.getAlignmentsMap(alphabet.get());
    auto mpTree = bppancestor.getPhyloTreesMap(mSites, unparsedParams);
    auto SPC=bppancestor.getCollection(alphabet.get(), gCode.get(), mSites, mpTree, unparsedParams);
    auto mSeqEvol = bppancestor.getProcesses(*SPC, unparsedParams);
    auto mPhyl=bppancestor.getPhyloLikelihoods(context, mSeqEvol, *SPC, mSites);
      
    // retrieve Phylo 0, aka result phylolikelihood
      
    if (!mPhyl->hasPhyloLikelihood(0))
      throw Exception("Missing phyloLikelihoods.");

    PhyloLikelihood* tl=(*mPhyl)[0];
    
    bppancestor.fixLikelihood(alphabet.get(), gCode.get(), tl);
    
    bppancestor.displayParameters(*tl, false);

    ApplicationTools::displayMessage("");

    //////////////////////////////////////
    // Reconstruct ancestral sequences:

    
    /// map of the Single Data Process 
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

    ////////////////////////
    /// Options

    // Ancestral information
    // Sites
    
    string outputSitesFile = ApplicationTools::getAFilePath("output.sites.file", allParams, false, false);

    bool probs(false);

    if (outputSitesFile!="none")
    {
      probs = ApplicationTools::getBooleanParameter("output.sites.probabilities", allParams, false, "", true, 1);
      if (!probs)
        probs = ApplicationTools::getBooleanParameter("asr.probabilities", allParams, false, "", true, 1);
        
      ApplicationTools::displayResult("Output site probabilities", probs ? "yes" : "no");
    }
    
    // Nodes
    
    bool addNodesExtant = false;
    string outputNodesFile = ApplicationTools::getAFilePath("output.nodes.file", allParams, false, false, "", false, "none", 1);
    if (outputNodesFile!="none")
    {
      addNodesExtant = ApplicationTools::getBooleanParameter("output.nodes.add_extant", allParams, false, "", true, 1);    
      ApplicationTools::displayResult("Output extant nodes", addNodesExtant ? "yes" : "no");
    }
    
    
    // ASR
    
    string sequenceFilePath = ApplicationTools::getAFilePath("asr.sequence.file", allParams, false, false, "", false, "none", 1);
    if (sequenceFilePath=="none")
      sequenceFilePath = ApplicationTools::getAFilePath("output.sequence.file", allParams, false, false, "", false, "none", 1);
    
    string sequenceFormat="";
    bool sample=false;
    unsigned int nbSamples=0;
    bool addSitesExtant=false;

    if (sequenceFilePath!="none")
    {
      sequenceFormat   = ApplicationTools::getStringParameter("asr.sequence.format", allParams, "none", "", false, 1);
      if (sequenceFormat=="none")
        sequenceFormat   = ApplicationTools::getStringParameter("output.sequence.format", allParams, "Fasta", "", false, 1);

      sample = ApplicationTools::getBooleanParameter("asr.sample", allParams, false, "", true, 1);
    
      ApplicationTools::displayResult("Sample from posterior distribution", sample ? "yes" : "no");
    
      if (sample)
        nbSamples = ApplicationTools::getParameter<unsigned int>("asr.sample.number", allParams, 1, "", true, false);
      
      addSitesExtant = ApplicationTools::getBooleanParameter("asr.add_extant", allParams, false, "", true, 1);
      ApplicationTools::displayResult("ASR extant", addSitesExtant ? "yes" : "no");
    }
    

    /////////////////////////////////  
    /// per Single Data Process 
    
    for (auto& itm:mSD)
    {
      SingleProcessPhyloLikelihood* sPP=dynamic_cast<SingleProcessPhyloLikelihood*>(itm.second);
      OneProcessSequencePhyloLikelihood* oPSP=dynamic_cast<OneProcessSequencePhyloLikelihood*>(itm.second);

      if ((sPP==NULL) && (oPSP==NULL))
      {
        ApplicationTools::displayWarning("Multi Process ancestral reconstruction not implemented");
        ApplicationTools::displayWarning("for phyloLikelihood " + TextTools::toString(itm.first));
        continue;
      }

      const AlignedValuesContainer* sites=sPP?sPP->getData():oPSP->getData();
            
      auto pDR=sPP?sPP->getLikelihoodCalculationSingleProcess(): oPSP->getLikelihoodCalculationSingleProcess();

      // Only Marginal reconstruction method
      AncestralStateReconstruction *asr = new MarginalAncestralReconstruction(pDR);

      size_t nbStates=sPP?sPP->getNumberOfStates():oPSP->getNumberOfStates();

      ApplicationTools::displayMessage("\nPhylo " + TextTools::toString(itm.first));
          
      /////////////////////////////////////
      // Write sites infos to file:

      if (outputSitesFile != "none")
      {
        string outF=outputSitesFile + "_" + TextTools::toString(itm.first);
        ApplicationTools::displayResult(" Output file for sites", outF);
        ofstream out(outF.c_str(), ios::out);
        PhyloTree ttree(sPP?sPP->getTree():oPSP->getTree());
        vector<shared_ptr<PhyloNode> > nodes = ttree.getAllNodes();
        size_t nbNodes = nodes.size();

        // Get the class with maximum posterior probability:
        vector<size_t> classes = sPP?sPP->getClassWithMaxPostProbPerSite():oPSP->getClassWithMaxPostProbPerSite();
        // Get the posterior rate, i.e. rate averaged over all posterior probabilities:

        Vdouble rates = sPP?sPP->getPosteriorRatePerSite():oPSP->getPosteriorRatePerSite();

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

        for (size_t i = 0; i < sites->getNumberOfSites(); i++)
        {
          double lnL = sPP?sPP->getLogLikelihoodForASite(i):oPSP->getLogLikelihoodForASite(i);
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

      
      /////////////////////////////////////
      // Write nodes infos to file:

      if (outputNodesFile != "none")
      {        
        string outF=outputNodesFile  + "_" + TextTools::toString(itm.first);
        ApplicationTools::displayResult(" Output file for nodes", outF);

        ofstream out(outF.c_str(), ios::out);
        // map<int, vector<double> > frequencies;

        const auto& tree = pDR->getSubstitutionProcess().getParametrizablePhyloTree();
        
        auto allIndex = addNodesExtant? tree.getAllNodesIndexes(): tree.getAllInnerNodesIndexes();

        // Former output of bppAncestor
        
        // if (oPSP)
        //   oPSP->getAncestralFrequencies(frequencies, addNodesExtant);
        // else
        //   sPP->getAncestralFrequencies(frequencies, addNodesExtant);

        vector<string> colNames;
        colNames.push_back("Nodes");
        // for (size_t i = 0; i < nbStates; i++)
        //   colNames.push_back("exp" + (sPP?sPP->getData()->getAlphabet()->intToChar((int)i):oPSP->getData()->getAlphabet()->intToChar((int)i)));
        for (size_t i = 0; i < nbStates; i++)
          colNames.push_back("eb" + (sPP?sPP->getData()->getAlphabet()->intToChar((int)i):oPSP->getData()->getAlphabet()->intToChar((int)i)));
      
        //Now fill the table:
        vector<string> row(colNames.size());
        DataTable* infos = new DataTable(colNames);

        for (const auto& index:allIndex)
        {
          row[0] = TextTools::toString(index);

          Vdouble ebFreqs =sPP?sPP->getPosteriorStateFrequencies(index): oPSP->getPosteriorStateFrequencies(index);

          // for (size_t i = 0; i < nbStates; i++)
          // {
          //   row[i + 1] = TextTools::toString(itf.second[i]);
          // }
          for (size_t i = 0; i < nbStates; i++)
          {
            row[i + 1] = TextTools::toString(ebFreqs[i]);
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
        
        //Write output:
        BppOAlignmentWriterFormat bppoWriter(1);
        unique_ptr<OAlignment> oAln(bppoWriter.read(sequenceFormat));
        ApplicationTools::displayResult(" Output alignment file ", sequenceFilePath + "_" + TextTools::toString(itm.first));
        ApplicationTools::displayResult(" Output alignment format ", oAln->getFormatName());

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
              
            sampleSites->setSequencesNames(names, true);

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
            ApplicationTools::displayWarning("Output extant sequences not possible with probabilistic sequences.");
          else
            SequenceContainerTools::append(*asSites, *vSC);
        }
          
          
        // Write sequences:
        oAln->writeAlignment(sequenceFilePath + "_" + TextTools::toString(itm.first), *asSites, true);
          
        delete asSites;
        delete asr;
      }      

      /// end of ancestral reconstruction

    }
  
    ApplicationTools::displayMessage("");
    bppancestor.done();
  
  }
  catch (exception & e)
  {
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}

