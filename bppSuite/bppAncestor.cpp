// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

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
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Likelihood/MarginalAncestralReconstruction.h>
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/OneProcessSequencePhyloLikelihood.h>
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/AlignedPhyloLikelihoodSet.h>
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>

using namespace bpp;

/******************************************************************************/

int main(int args, char** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*     Bio++ Ancestral Sequence Reconstruction, version " << BPP_VERSION << "     *" << endl;
  cout << "* Authors: J. Dutheil                       Created on: 10/09/08 *" << endl;
  cout << "*          B. Boussau                       Last Modif: 25/09/14 *" << endl;
  cout << "*          L. Guéguen                       Last Modif: 22/12/14 *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  try
  {
    BppPhylogeneticsApplication bppancestor(args, argv, "bppancestor");

    if (args == 1)
    {
      bppancestor.help("bppAncestor");
      return 0;
    }

    bppancestor.startTimer();

    Context context;
    map<string, string> allParams = bppancestor.getParams();
    map<string, string> unparsedParams;

    shared_ptr<const Alphabet> alphabet(bppancestor.getAlphabet());
    shared_ptr<const GeneticCode> gCode(bppancestor.getGeneticCode(alphabet));

    // Missing check
    //  if (model->getName() != "RE08") SiteContainerTools::changeGapsToUnknownCharacters(*sites);

    // get the result phylo likelihood
    auto mSitesuniq = bppancestor.getConstAlignmentsMap(alphabet, true);

    const std::map<size_t, std::shared_ptr<const AlignmentDataInterface >> mSites = PhylogeneticsApplicationTools::uniqueToSharedMap<const TemplateAlignmentDataInterface<string>>(mSitesuniq);

    auto mpTree = bppancestor.getPhyloTreesMap(mSites, unparsedParams);
    shared_ptr<SubstitutionProcessCollection> SPC = bppancestor.getCollection(alphabet, gCode, mSites, mpTree, unparsedParams);
    auto mSeqEvoltmp = bppancestor.getProcesses(SPC, unparsedParams);
    auto mSeqEvol = PhylogeneticsApplicationTools::uniqueToSharedMap<SequenceEvolution>(mSeqEvoltmp);

    auto mPhyl = bppancestor.getPhyloLikelihoods(context, mSeqEvol, SPC, mSites);

    // retrieve Phylo 0, aka result phylolikelihood

    if (!mPhyl->hasPhyloLikelihood(0))
      throw Exception("Missing phyloLikelihoods.");

    auto tl = (*mPhyl)[0];

    bppancestor.fixLikelihood(alphabet, gCode, tl);

    bppancestor.displayParameters(*tl, false);

    ApplicationTools::displayMessage("");

    //////////////////////////////////////
    // Reconstruct ancestral sequences:


    /// map of the Single Data Process
    map<size_t, shared_ptr<AbstractSingleDataPhyloLikelihood>> mSD;

    if (dynamic_pointer_cast<AbstractSingleDataPhyloLikelihood>(tl) != NULL)
      mSD[1] = dynamic_pointer_cast<AbstractSingleDataPhyloLikelihood>(tl);
    else
    {
      auto sOAP = dynamic_pointer_cast<AlignedPhyloLikelihoodSetInterface>(tl);
      if (sOAP)
      {
        const vector<size_t>& nSD = sOAP->getNumbersOfPhyloLikelihoods();

        for (size_t iSD = 0; iSD < nSD.size(); iSD++)
        {
          auto pASDP = dynamic_pointer_cast<AbstractSingleDataPhyloLikelihood>(sOAP->getPhyloLikelihood(nSD[iSD]));
          if (pASDP != NULL)
            mSD[nSD[iSD]] = pASDP;
        }
      }
    }

    ////////////////////////
    /// Options

    // Ancestral information
    // Sites

    string outputSitesFile = ApplicationTools::getAFilePath("output.sites.file", allParams, false, false);

    bool probs(false);

    if (outputSitesFile != "none")
    {
      probs = ApplicationTools::getBooleanParameter("output.sites.probabilities", allParams, false, "", true, 1);
      if (!probs)
        probs = ApplicationTools::getBooleanParameter("asr.probabilities", allParams, false, "", true, 1);

      ApplicationTools::displayResult("Output site probabilities", probs ? "yes" : "no");
    }

    // Nodes

    bool addNodesExtant = false;
    string outputNodesFile = ApplicationTools::getAFilePath("output.nodes.file", allParams, false, false, "", false, "none", 1);

    if (outputNodesFile != "none")
    {
      addNodesExtant = ApplicationTools::getBooleanParameter("output.nodes.add_extant", allParams, false, "", true, 1);
      ApplicationTools::displayResult("Output extant nodes", addNodesExtant ? "yes" : "no");
    }

    // ASR

    string sequenceFilePath = ApplicationTools::getAFilePath("asr.sequence.file", allParams, false, false, "", false, "none", 1);

    if (sequenceFilePath == "none")
      sequenceFilePath = ApplicationTools::getAFilePath("output.sequence.file", allParams, false, false, "", false, "none", 1);

    string sequenceFormat = "";
    bool sample = false;
    unsigned int nbSamples = 0;
    bool addSitesExtant = false;

    if (sequenceFilePath != "none")
    {
      sequenceFormat   = ApplicationTools::getStringParameter("asr.sequence.format", allParams, "none", "", false, 1);
      if (sequenceFormat == "none")
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
      auto sPP = dynamic_pointer_cast<SingleProcessPhyloLikelihood>(itm.second);
      auto oPSP = dynamic_pointer_cast<OneProcessSequencePhyloLikelihood>(itm.second);

      if ((sPP == NULL) && (oPSP == NULL))
      {
        ApplicationTools::displayWarning("Multi Process ancestral reconstruction not implemented");
        ApplicationTools::displayWarning("for phyloLikelihood " + TextTools::toString(itm.first));
        continue;
      }

      auto pDR = sPP ? sPP->getLikelihoodCalculationSingleProcess() : oPSP->getLikelihoodCalculationSingleProcess();
      auto sites = std::shared_ptr<const AlignmentDataInterface>(sPP ? sPP->getData() : oPSP->getData());

      // Only Marginal reconstruction method

      VectorSequenceContainer vSC(alphabet);

      if (addSitesExtant)
      {
        auto vSC0 = dynamic_pointer_cast<const SiteContainerInterface>(sites);
        if (!vSC0)
          ApplicationTools::displayWarning("Output extant sequences not possible with probabilistic sequences.");
        // Keep only leaves
        const auto tree = sPP ? sPP->tree() : oPSP->tree();
        const auto leaves = tree->getAllLeavesNames();

        SequenceContainerTools::getSelectedSequences(*vSC0, leaves, vSC);
      }

      AncestralStateReconstruction* asr = new MarginalAncestralReconstruction(pDR);

      size_t nbStates = sPP ? sPP->getNumberOfStates() : oPSP->getNumberOfStates();

      const auto& stateMap = pDR->stateMap();

      ApplicationTools::displayMessage("\nPhylo " + TextTools::toString(itm.first));

      /////////////////////////////////////
      // Write sites infos to file:

      if (outputSitesFile != "none")
      {
        string outF = outputSitesFile + "_" + TextTools::toString(itm.first);
        ApplicationTools::displayResult(" Output file for sites", outF);
        ofstream out(outF.c_str(), ios::out);
        auto ttree(sPP ? sPP->tree() : oPSP->tree());
        vector<shared_ptr<PhyloNode>> nodes = ttree->getAllNodes();
        size_t nbNodes = nodes.size();

        // Get the class with maximum posterior probability:
        vector<size_t> classes = sPP ? sPP->getClassWithMaxPostProbPerSite() : oPSP->getClassWithMaxPostProbPerSite();
        // Get the posterior rate, i.e. rate averaged over all posterior probabilities:

        Vdouble rates = sPP ? sPP->getPosteriorRatePerSite() : oPSP->getPosteriorRatePerSite();

        // Get the ancestral sequences:
        vector<unique_ptr<Sequence>> sequences(nbNodes);
        vector<VVdouble*> probabilities(nbNodes);

        vector<string> colNames;
        colNames.push_back("Sites");
        colNames.push_back("is.complete");
        colNames.push_back("is.constant");
        colNames.push_back("lnL");
        colNames.push_back("rc");
        colNames.push_back("pr");

        for (size_t i = 0; i < nbNodes; i++)
        {
          shared_ptr<PhyloNode> node = nodes[i];
          auto nodeindex = ttree->getNodeIndex(node);
          colNames.push_back("max." + TextTools::toString(nodeindex));
          if (probs)
          {
            probabilities[i] = new VVdouble();

            // The cast will have to be updated when more probabilistic method will be available:
            sequences[i] = dynamic_cast<MarginalAncestralReconstruction*>(asr)->getAncestralSequenceForNode(nodeindex, probabilities[i], false);

            for (unsigned int j = 0; j < nbStates; j++)
            {
              colNames.push_back("prob." + TextTools::toString(nodeindex) + "." + stateMap.getStateDescription(j));
            }
          }
          else
            sequences[i] = asr->getAncestralSequenceForNode(nodeindex);
        }

        // Now fill the table:
        vector<string> row(colNames.size());
        auto infos = make_unique<DataTable>(colNames);

        for (size_t i = 0; i < sites->getNumberOfSites(); i++)
        {
          double lnL = sPP ? sPP->getLogLikelihoodForASite(i) : oPSP->getLogLikelihoodForASite(i);

          auto& currentSite = sites->site(i);
          int currentSitePosition = currentSite.getCoordinate();
          string isCompl = "NA";
          string isConst = "NA";
          try
          {
            isCompl = (SiteTools::isComplete(currentSite) ? "1" : "0");
          }
          catch (EmptySiteException& ex)
          {}
          try
          {
            isConst = (SiteTools::isConstant(currentSite) ? "1" : "0");
          }
          catch (EmptySiteException& ex)
          {}
          row[0] = (string("[" + TextTools::toString(currentSitePosition) + "]"));
          row[1] = isCompl;
          row[2] = isConst;
          row[3] = TextTools::toString(lnL);
          row[4] = TextTools::toString(classes[i]);
          row[5] = TextTools::toString(rates[i]);

          unsigned int k = 6;
          for (unsigned int j = 0; j < nbNodes; j++)
          {
            row[k] = sequences[j]->getChar(i);
            k++;
            if (probs)
            {
              for (unsigned int l = 0; l < nbStates; l++)
              {
                row[k] = TextTools::toString((*probabilities[j])[i][l]);
                k++;
              }
            }
          }

          infos->addRow(row);
        }

        DataTable::write(*infos, out, "\t");
      }


      /////////////////////////////////////
      // Write nodes infos to file:

      if (outputNodesFile != "none")
      {
        string outF = outputNodesFile  + "_" + TextTools::toString(itm.first);
        ApplicationTools::displayResult(" Output file for nodes", outF);

        ofstream out(outF.c_str(), ios::out);
        // map<int, vector<double> > frequencies;

        auto tree = pDR->getSubstitutionProcess()->getParametrizablePhyloTree();

        auto allIndex = addNodesExtant ? tree->getAllNodesIndexes() : tree->getAllInnerNodesIndexes();

        // Former output of bppAncestor


        vector<string> colNames;
        colNames.push_back("Nodes");
        // for (size_t i = 0; i < nbStates; i++)
        //   colNames.push_back("exp" + (sPP?sPP->getData()->getAlphabet()->intToChar((int)i):oPSP->getData()->getAlphabet()->intToChar((int)i)));
        for (size_t i = 0; i < nbStates; i++)
        {
          colNames.push_back("eb" + stateMap.getStateDescription(i));
        }

        // Now fill the table:
        vector<string> row(colNames.size());
        auto infos = make_unique<DataTable>(colNames);

        for (const auto& index:allIndex)
        {
          row[0] = TextTools::toString(index);

          Vdouble ebFreqs = sPP ? sPP->getPosteriorStateFrequencies(index) : oPSP->getPosteriorStateFrequencies(index);

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
      }

      ////////////////////////////////////////////////
      /// Ancestral sequences

      if (sequenceFilePath != "none")
      {
        shared_ptr<AlignedSequenceContainer> asSites = nullptr;

        // Write output:
        BppOAlignmentWriterFormat bppoWriter(1);
        unique_ptr<OAlignment> oAln(bppoWriter.read(sequenceFormat));
        ApplicationTools::displayResult(" Output alignment file ", sequenceFilePath + "_" + TextTools::toString(itm.first));
        ApplicationTools::displayResult(" Output alignment format ", oAln->getFormatName());

        if (sample)
        {
          asSites = make_shared<AlignedSequenceContainer>(alphabet);

          for (unsigned int i = 0; i < nbSamples; i++)
          {
            ApplicationTools::displayGauge(i, nbSamples - 1, '=');
            shared_ptr<AlignmentDataInterface> sampleSites = dynamic_cast<MarginalAncestralReconstruction*>(asr)->getAncestralSequences(true);
            vector<string> names = sampleSites->getSequenceNames();

            for (unsigned int j = 0; j < names.size(); j++)
            {
              names[j] += "_" + TextTools::toString(i + 1);
            }

            sampleSites->setSequenceNames(names, true);

            SequenceContainerTools::append(*asSites, *dynamic_pointer_cast<const SequenceContainerInterface>(sampleSites));
          }
          ApplicationTools::message->endLine();
        }
        else
          asSites = asr->getAncestralSequences();

        // Add existing sequence to output?
        if (addSitesExtant)
          SequenceContainerTools::append(*asSites, vSC);

        // Write sequences:
        oAln->writeAlignment(sequenceFilePath + "_" + TextTools::toString(itm.first), *asSites, true);
      }

      /// end of ancestral reconstruction
    }

    ApplicationTools::displayMessage("");
    bppancestor.done();
  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}
