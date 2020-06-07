//
// File: bppPopStats.cpp
// Created by: Julien Dutheil
// Created on: Jun Wed 24 12:04 2015
//

/*
   Copyright or © or Copr. Bio++ Development Team

   This software is a computer program whose purpose is to simulate sequence
   data according to a phylogenetic tree and an evolutionary model.

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
#include <fstream>
#include <iomanip>
#include <memory>
#include <cmath>

using namespace std;

// From bpp-core:
#include <Bpp/Version.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/KeyvalTools.h>

// From bpp-seq:
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/CodonSiteTools.h>
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

// From bpp-phyl:
#include <Bpp/Phyl/Tree/Tree.h>
#include <Bpp/Phyl/Model/Nucleotide/K80.h>
#include <Bpp/Phyl/Model/Codon/YN98.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Distance/DistanceEstimation.h>
#include <Bpp/Phyl/Distance/BioNJ.h>
#include <Bpp/Phyl/Likelihood/DRHomogeneousTreeLikelihood.h>
#include <Bpp/Phyl/Likelihood/MarginalAncestralStateReconstruction.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>

// From bpp-popgen
#include <Bpp/PopGen/PolymorphismSequenceContainer.h>
#include <Bpp/PopGen/PolymorphismSequenceContainerTools.h>
#include <Bpp/PopGen/SequenceStatistics.h>

using namespace bpp;

void help()
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
  (*ApplicationTools::message << "bpppopstats parameter1_name=parameter1_value").endLine();
  (*ApplicationTools::message << "      parameter2_name=parameter2_value ... param=option_file").endLine();
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "  Refer to the Bio++ Program Suite Manual for a list of available options.").endLine();
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
}

int main(int args, char** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*              Bio++ Population Statistics, version " << BPP_VERSION << "        *" << endl;
  cout << "* Author: J. Dutheil                        Last Modif. " << BPP_REL_DATE << " *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  if (args == 1)
  {
    help();
    return 0;
  }

  BppApplication bpppopstats(args, argv, "BppPopStats");
  bpppopstats.startTimer();

  string logFile = ApplicationTools::getAFilePath("logfile", bpppopstats.getParams(), false, false);
  unique_ptr<ofstream> cLog;
  if (logFile != "none")
    cLog.reset(new ofstream(logFile.c_str(), ios::out));

  //This counts instances of each tool, in case one is used several times, for instance with different options:
  map<string, unsigned int> toolCounter;

  try
  {
    // Get alphabet
    Alphabet* alphabet = SequenceApplicationTools::getAlphabet(bpppopstats.getParams(), "", false, true, true);

    // Get the genetic code, if codon alphabet
    unique_ptr<GeneticCode> gCode;
    CodonAlphabet* codonAlphabet = dynamic_cast<CodonAlphabet*>(alphabet);
    if (codonAlphabet) {
      string codeDesc = ApplicationTools::getStringParameter("genetic_code", bpppopstats.getParams(), "Standard", "", true, true);
      ApplicationTools::displayResult("Genetic Code", codeDesc);
      gCode.reset(SequenceApplicationTools::getGeneticCode(codonAlphabet->shareNucleicAlphabet(), codeDesc));
    }

    unique_ptr<PolymorphismSequenceContainer> psc;
    if (ApplicationTools::parameterExists("input.sequence.file.ingroup", bpppopstats.getParams())) {
      // Get the ingroup alignment:
      unique_ptr<SiteContainer> sitesIn(SequenceApplicationTools::getSiteContainer(alphabet, bpppopstats.getParams(), ".ingroup", false, true));
      psc.reset(new PolymorphismSequenceContainer(*sitesIn));
      if (ApplicationTools::parameterExists("input.sequence.file.outgroup", bpppopstats.getParams())) {
        // Get the outgroup alignment:
        unique_ptr<SiteContainer> sitesOut(SequenceApplicationTools::getSiteContainer(alphabet, bpppopstats.getParams(), ".outgroup", false, true));
        SequenceContainerTools::append(*psc, *sitesOut);
        for (size_t i = sitesIn->getNumberOfSequences(); i < psc->getNumberOfSequences(); ++i) {
          psc->setAsOutgroupMember(i);
        }
      }
    } else {
      //Everything in one file
      unique_ptr<SiteContainer> sites(SequenceApplicationTools::getSiteContainer(alphabet, bpppopstats.getParams(), "", false, true));
      psc.reset(new PolymorphismSequenceContainer(*sites));
      if (ApplicationTools::parameterExists("input.sequence.outgroup.index", bpppopstats.getParams())) {
        vector<size_t> outgroups = ApplicationTools::getVectorParameter<size_t>("input.sequence.outgroup.index", bpppopstats.getParams(), ',', "");
        for (auto g : outgroups) {
          psc->setAsOutgroupMember(g-1);
        }
      }
      if (ApplicationTools::parameterExists("input.sequence.outgroup.name", bpppopstats.getParams())) {
        vector<string> outgroups = ApplicationTools::getVectorParameter<string>("input.sequence.outgroup.name", bpppopstats.getParams(), ',', "");
        for (auto g : outgroups) {
          ApplicationTools::displayResult("Sequence from outgroup", g);
          psc->setAsOutgroupMember(g);
        }
      }
    }

    // Take care of stop codons:
    string stopCodonOpt = ApplicationTools::getStringParameter("input.sequence.stop_codons_policy", bpppopstats.getParams(), "Keep", "", true, true);
    ApplicationTools::displayResult("Stop codons policy", stopCodonOpt);

    if (stopCodonOpt == "Keep") {
      //do nothing
    } else if (stopCodonOpt == "RemoveIfLast") {
      if (CodonSiteTools::hasStop(psc->getSite(psc->getNumberOfSites() - 1), *gCode)) {
        psc->deleteSite(psc->getNumberOfSites() - 1);
        ApplicationTools::displayMessage("Info: last site contained a stop codon and was discarded.");
        if (logFile != "none")
          *cLog << "# Info: last site contained a stop codon and was discarded." << endl;
      }
    } else if (stopCodonOpt == "RemoveAll") {
      size_t l1 = psc->getNumberOfSites();
      SiteContainerTools::removeStopCodonSites(*psc, *gCode);
      size_t l2 = psc->getNumberOfSites();
      if (l2 != l1) {
        ApplicationTools::displayMessage("Info: discarded " + TextTools::toString(l1 - l2) + " sites with stop codons.");
        if (logFile != "none")
          *cLog << "# Info: discarded " << (l1 - l2) << " sites with stop codons." << endl;
      }
    } else {
      throw Exception("Unrecognized option for input.sequence.stop_codons_policy: " + stopCodonOpt);
    }

    shared_ptr<PolymorphismSequenceContainer> pscIn;
    shared_ptr<PolymorphismSequenceContainer> pscOut;

    if (psc->hasOutgroup()) {
      pscIn.reset(PolymorphismSequenceContainerTools::extractIngroup(*psc));
      pscOut.reset(PolymorphismSequenceContainerTools::extractOutgroup(*psc));
    } else {
      pscIn = std::move(psc);
    }
    ApplicationTools::displayResult("Number of sequences in ingroup", pscIn->getNumberOfSequences());
    ApplicationTools::displayResult("Number of sequences in outgroup", pscOut.get() ? pscOut->getNumberOfSequences() : 0);
   
    // Shall we estimate some parameters first?
    
    bool estimateTsTv = ApplicationTools::getBooleanParameter("estimate.kappa", bpppopstats.getParams(), false, "", false, 1);
    double kappa = 1;
    double omega = -1;
    
    bool estimateAncestor = ApplicationTools::getBooleanParameter("estimate.ancestor", bpppopstats.getParams(), false, "", false, 1);
    if (estimateAncestor & ! pscOut)
      throw Exception("Error: an outgroup sequence is needed for estimating ancestral states.");

    bool fitModel = estimateTsTv || estimateAncestor;

    // Fit a model for later use:
    unique_ptr<Tree> tree;
    unique_ptr<SubstitutionModel> model;
    unique_ptr<DiscreteDistribution> rDist;
    DRTreeLikelihood* treeLik = nullptr;
    unique_ptr<Sequence> ancestralSequence;
    if (fitModel) {
      // Get the alignment:
      bool sampleIngroup = ApplicationTools::getBooleanParameter("estimate.sample_ingroup", bpppopstats.getParams(), true);
      size_t sampleIngroupSize = 0;
      if (sampleIngroup) {
        sampleIngroupSize = ApplicationTools::getParameter<size_t>("estimate.sample_ingroup.size", bpppopstats.getParams(), 10);
        if (sampleIngroupSize > pscIn->getNumberOfSequences()) {
          ApplicationTools::displayWarning("Sample size higher than number of sequence. No sampling performed.");
          sampleIngroup = false;
        }
      }
      unique_ptr<AlignedSequenceContainer> aln;
      if (sampleIngroup) {
        ApplicationTools::displayResult("Nb of ingroup sequences for model fitting", sampleIngroupSize);
        aln.reset(new AlignedSequenceContainer(pscIn->getAlphabet()));
        vector<string> selection(sampleIngroupSize);
        RandomTools::getSample(pscIn->getSequencesNames(), selection, false); 
        SequenceContainerTools::getSelectedSequences(*pscIn, selection ,*aln);
      } else {
        aln.reset(new AlignedSequenceContainer(*pscIn));
      }
      if (pscOut) {
        aln->addSequence(pscOut->getSequence(0)); //As for now, we only consider one sequence as outgroup, the first one.
      }

      // Get a tree:
      string treeOpt = ApplicationTools::getStringParameter("input.tree.method", bpppopstats.getParams(), "bionj", "");
      if (codonAlphabet) {
        unique_ptr<FrequenciesSet> freqSet(new FixedCodonFrequenciesSet(gCode.get()));
        model.reset(new YN98(gCode.get(), freqSet.release()));
      } else {
        model.reset(new K80(&AlphabetTools::DNA_ALPHABET));
      } //Note: proteins not supported!
      rDist.reset(new ConstantRateDistribution()); 
      if (treeOpt == "user") {
        tree.reset(PhylogeneticsApplicationTools::getTree(bpppopstats.getParams()));
      } else if (treeOpt == "bionj") {
        ApplicationTools::displayTask("Estimating distance matrix", true);
        //DistanceEstimation distEstimation(model->clone(), rDist->clone(), aln.get(), 1, false);
        //distEstimation.computeMatrix();
        //unique_ptr<DistanceMatrix> matrix(distEstimation.getMatrix());
        unique_ptr<DistanceMatrix> matrix(SiteContainerTools::computeSimilarityMatrix(*aln, true, SiteContainerTools::SIMILARITY_NOGAP, true));
        
        ApplicationTools::displayTaskDone();
        ApplicationTools::displayTask("Computing BioNJ tree", true);
        BioNJ bionj(false, true);
        bionj.setDistanceMatrix(*matrix);
        bionj.computeTree();
        ApplicationTools::displayTaskDone();
        tree.reset(bionj.getTree());
      } else {
        throw Exception("Invalid input.tree.method. Should be either 'user' or 'bionj'.");
      }

      // Create a likelihood object:
      treeLik = new DRHomogeneousTreeLikelihood(*tree, *aln, model.get(), rDist.get());
      treeLik->initialize();
      if (std::isinf(treeLik->getValue()))
        throw Exception("Error: null likelihood. Possible cause: stop codon or numerical underflow (too many sequences).");
      // Optimize parameters:
      treeLik = dynamic_cast<DRTreeLikelihood*>(PhylogeneticsApplicationTools::optimizeParameters(treeLik, treeLik->getParameters(), bpppopstats.getParams(), "", true, true, 2));

      // Get kappa:              
      if (estimateTsTv) {
        kappa = model->getParameter("kappa").getValue();
        ApplicationTools::displayResult("Transition / transversions ratio", kappa);
      }
      if (estimateAncestor) {
        MarginalAncestralStateReconstruction asr(treeLik);
        int outgroupId = tree->getLeafId(pscOut->getSequence(0).getName());
        ancestralSequence.reset(asr.getAncestralSequenceForNode(tree->getFatherId(outgroupId)));
      }
      if (codonAlphabet) {
        omega = model->getParameter("omega").getValue();
      }
    }
    if (treeLik)
      delete treeLik; //Not needed anymore.

    // +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    
    // Compute statistics
    vector<string> actions = ApplicationTools::getVectorParameter<string>("pop.stats", bpppopstats.getParams(), ',', "", "", false, 1);

    for (size_t a = 0; a < actions.size(); a++)
    {
      string cmdName;
      map<string, string> cmdArgs;
      KeyvalTools::parseProcedure(actions[a], cmdName, cmdArgs);
      toolCounter[cmdName]++;

      // +-------------------+
      // | Frequencies       |
      // +-------------------+
      if (cmdName == "SiteFrequencies")
      {
        unsigned int s = SequenceStatistics::numberOfPolymorphicSites(*pscIn);
        ApplicationTools::displayResult("Number of segregating sites:", s);
        unsigned int nsg = SequenceStatistics::numberOfSingletons(*pscIn);
        ApplicationTools::displayResult("Number of singletons:", nsg);
        //Print to logfile:
        if (logFile != "none") {
          *cLog << "# Site frequencies" << endl;
          *cLog << "NbSegSites" << (toolCounter[cmdName] > 1 ? TextTools::toString(toolCounter[cmdName]) : "") << " = " << s << endl;
          *cLog << "NbSingl" << (toolCounter[cmdName] > 1 ? TextTools::toString(toolCounter[cmdName]) : "") << " = " << nsg << endl;
        }
      }

      // +-------------------+
      // | Watterson's theta |
      // +-------------------+
      else if (cmdName == "Watterson75")
      {
        double thetaW75 = SequenceStatistics::watterson75(*pscIn, true, true, true);
        ApplicationTools::displayResult("Watterson's (1975) theta:", thetaW75);
        //Print to logfile:
        if (logFile != "none") {
          *cLog << "# Watterson's (1975) theta" << endl;
          *cLog << "thetaW75" << (toolCounter[cmdName] > 1 ? TextTools::toString(toolCounter[cmdName]) : "") << " = " << thetaW75 << endl;
        }
      }

      // +-------------+
      // | Tajima's pi |
      // +-------------+
      else if (cmdName == "Tajima83")
      {
        double piT83 = SequenceStatistics::tajima83(*pscIn, true, true, true);
        ApplicationTools::displayResult("Tajima's (1983) pi:", piT83);
        //Print to logfile:
        if (logFile != "none") {
          *cLog << "# Tajima's (1983) pi" << endl;
          *cLog << "piT83" << (toolCounter[cmdName] > 1 ? TextTools::toString(toolCounter[cmdName]) : "") << " = " << piT83 << endl;
        }
      }

      // +------------+
      // | Tajima's D |
      // +------------+
      else if (cmdName == "TajimaD")
      {
        string positions = ApplicationTools::getStringParameter("positions", cmdArgs, "all", "", false, 1);
        shared_ptr<PolymorphismSequenceContainer> pscTmp;
        if ((positions == "synonymous" || positions == "non-synonymous") && !codonAlphabet)
          throw Exception("Error: synonymous and non-synonymous positions can only be defined with a codon alphabet.");
        if (positions == "synonymous") {
          pscTmp.reset(PolymorphismSequenceContainerTools::getSynonymousSites(*pscIn, *gCode));
        } else if (positions == "non-synonymous") {
          pscTmp.reset(PolymorphismSequenceContainerTools::getNonSynonymousSites(*pscIn, *gCode));
        } else if (positions == "all") {
          pscTmp = pscIn;
        } else throw Exception("Unrecognized option for argument 'positions': " + positions);

        if (SequenceStatistics::numberOfPolymorphicSites(*pscTmp) > 0) {
          double tajimaD = SequenceStatistics::tajimaDss(*pscTmp, true, true);
          ApplicationTools::displayResult("Tajima's (1989) D:", tajimaD);
          //Print to logfile:
          if (logFile != "none") {
            *cLog << "# Tajima's (1989) D (" << positions << " sites)" << endl;
            *cLog << "tajD" << (toolCounter[cmdName] > 1 ? TextTools::toString(toolCounter[cmdName]) : "") << " = " << tajimaD << endl;
          }
        } else {
          ApplicationTools::displayResult<string>("Tajima's (1989) D:", "NA (0 polymorphic sites)");
          if (logFile != "none") {
            *cLog << "# Tajima's (1989) D (" << positions << " sites)" << endl;
            *cLog << "tajD" << (toolCounter[cmdName] > 1 ? TextTools::toString(toolCounter[cmdName]) : "") << " = NA" << endl;
          }
          
        }
      }

      // +-----------+
      // | FuAndLiD* |
      // +-----------+
      else if (cmdName == "FuAndLiDStar")
      {
        bool useTotMut = ApplicationTools::getBooleanParameter("tot_mut", cmdArgs, true, "", false, 1);
        double flDstar = SequenceStatistics::fuLiDStar(*pscIn, !useTotMut);
        ApplicationTools::displayResult("Fu and Li's (1993) D*:", flDstar);
        ApplicationTools::displayResult("  computed using", (useTotMut ? "total number of mutations" : "number of segregating sites"));
        //Print to logfile:
        if (logFile != "none") {
          *cLog << "# Fu and Li's (1993) D*" << endl;
          if (useTotMut)
            *cLog << "fuLiDstarTotMut" << (toolCounter[cmdName] > 1 ? TextTools::toString(toolCounter[cmdName]) : "") << " = " << flDstar << endl;
          else
            *cLog << "fuLiDstarSegSit" << (toolCounter[cmdName] > 1 ? TextTools::toString(toolCounter[cmdName]) : "") << " = " << flDstar << endl;
        }
      }

      // +-----------+
      // | FuAndLiF* |
      // +-----------+
      else if (cmdName == "FuAndLiFStar")
      {
        bool useTotMut = ApplicationTools::getBooleanParameter("tot_mut", cmdArgs, true, "", false, 1);
        double flFstar = SequenceStatistics::fuLiFStar(*pscIn, !useTotMut);
        ApplicationTools::displayResult("Fu and Li (1993)'s F*:", flFstar);
        ApplicationTools::displayResult("  computed using", (useTotMut ? "total number of mutations" : "number of segregating sites"));
        //Print to logfile:
        if (logFile != "none") {
          *cLog << "# Fu and Li's (1993) F*" << endl;
          if (useTotMut)
            *cLog << "fuLiFstarTotMut" << (toolCounter[cmdName] > 1 ? TextTools::toString(toolCounter[cmdName]) : "") << " = " << flFstar << endl;
          else
            *cLog << "fuLiFstarSegSit" << (toolCounter[cmdName] > 1 ? TextTools::toString(toolCounter[cmdName]) : "") << " = " << flFstar << endl;
        }
      }
      // +-----------+
      // | PiN / PiS |
      // +-----------+
      else if (cmdName == "PiN_PiS")
      {
        if (!codonAlphabet) {
          throw Exception("PiN_PiS can only be used with a codon alignment. Check the input alphabet!");
        }
        double piS = SequenceStatistics::piSynonymous(*pscIn, *gCode);
        double piN = SequenceStatistics::piNonSynonymous(*pscIn, *gCode);
        double nbS = SequenceStatistics::meanNumberOfSynonymousSites(*pscIn, *gCode, kappa);
        double nbN = SequenceStatistics::meanNumberOfNonSynonymousSites(*pscIn, *gCode, kappa);
        double r = (piN / nbN) / (piS / nbS);
        ApplicationTools::displayResult("PiN:", piN);
        ApplicationTools::displayResult("PiS:", piS);
        ApplicationTools::displayResult("#N:", nbN);
        ApplicationTools::displayResult("#S:", nbS);
        ApplicationTools::displayResult("PiN / PiS (corrected for #N and #S):", r);
        if (fitModel) {
          ApplicationTools::displayResult("Omega (YN98 model):", omega);
        }

        if (logFile != "none") {
          *cLog << "# PiN and PiS" << endl;
          *cLog << "PiN" << (toolCounter[cmdName] > 1 ? TextTools::toString(toolCounter[cmdName]) : "") << " = " << piN << endl;
          *cLog << "PiS" << (toolCounter[cmdName] > 1 ? TextTools::toString(toolCounter[cmdName]) : "") << " = " << piS << endl;
          *cLog << "NbN" << (toolCounter[cmdName] > 1 ? TextTools::toString(toolCounter[cmdName]) : "") << " = " << nbN << endl;
          *cLog << "NbS" << (toolCounter[cmdName] > 1 ? TextTools::toString(toolCounter[cmdName]) : "") << " = " << nbS << endl;
          if (fitModel) {
            *cLog << "Omega" << (toolCounter[cmdName] > 1 ? TextTools::toString(toolCounter[cmdName]) : "") << " = " << omega << endl;
          }
        }
      }

      // +---------+
      // | dN / dS |
      // +---------+
      else if (cmdName == "dN_dS")
      {
        if (!codonAlphabet) {
          throw Exception("dN_dS can only be used with a codon alignment. Check the input alphabet!");
        }
        //Get consensus sequences:
        unique_ptr<SiteContainer> alnIn(pscIn->toSiteContainer());
        unique_ptr<SiteContainer> alnOut(pscOut->toSiteContainer());
        unique_ptr<Sequence> consensusIn(SiteContainerTools::getConsensus(*alnIn, "consIn", true, false));
        unique_ptr<Sequence> consensusOut(SiteContainerTools::getConsensus(*alnOut, "consOut", true, false));
        unique_ptr<AlignedSequenceContainer> alnCons(new AlignedSequenceContainer(codonAlphabet));
        alnCons->addSequence(*consensusIn);
        alnCons->addSequence(*consensusOut);
        unique_ptr<FrequenciesSet> freqSetDiv(new FixedCodonFrequenciesSet(gCode.get()));
        YN98* modelDiv = new YN98(gCode.get(), freqSetDiv.release());
        DiscreteDistribution* rDistDiv = new ConstantRateDistribution(); 
        DistanceEstimation distEstimation(modelDiv, rDistDiv, alnCons.get(), 0, false);
        distEstimation.setAdditionalParameters(modelDiv->getIndependentParameters());
        distEstimation.computeMatrix();
        unique_ptr<DistanceMatrix> matrix(distEstimation.getMatrix());
        ApplicationTools::displayResult("Yang and Nielsen's Omega (dN/dS):", modelDiv->getParameter("omega").getValue());
        ApplicationTools::displayResult("Yang and Nielsen's Kappa:", modelDiv->getParameter("kappa").getValue());
        ApplicationTools::displayResult("Yang and Nielsen's Distance:", (*matrix)(1,0));
        if (logFile != "none") {
          *cLog << "# dN and dS (Yang and Nielsen's 1998 substitution model)" << endl;
          *cLog << "OmegaDiv" << (toolCounter[cmdName] > 1 ? TextTools::toString(toolCounter[cmdName]) : "") << " = " << modelDiv->getParameter("omega").getValue() << endl;
          *cLog << "KappaDiv" << (toolCounter[cmdName] > 1 ? TextTools::toString(toolCounter[cmdName]) : "") << " = " << modelDiv->getParameter("kappa").getValue() << endl;
          *cLog << "DistanceDiv" << (toolCounter[cmdName] > 1 ? TextTools::toString(toolCounter[cmdName]) : "") << " = " << (*matrix)(1,0) << endl;
        }
      }

      // +---------+
      // | MK test |
      // +---------+
      else if (cmdName == "MKT")
      {
        if (!codonAlphabet) {
          throw Exception("MacDonald-Kreitman test can only be performed on a codon alignment. Check the input alphabet!");
        }
        if (!pscOut) {
          throw Exception("MacDonald-Kreitman test requires at least one outgroup sequence.");
        }
        vector<unsigned int> mktable = SequenceStatistics::mkTable(*pscIn, *pscOut, *gCode);
        ApplicationTools::displayResult("MK table, Pa:", mktable[0]);
        ApplicationTools::displayResult("MK table, Ps:", mktable[1]);
        ApplicationTools::displayResult("MK table, Da:", mktable[2]);
        ApplicationTools::displayResult("MK table, Ds:", mktable[3]);
        if (logFile != "none") {
          *cLog << "# MK table" << endl;
          *cLog << "# Pa Ps Da Ds" << endl;
          *cLog << "MKtable" << (toolCounter[cmdName] > 1 ? TextTools::toString(toolCounter[cmdName]) : "") << " = " << mktable[0] << " " << mktable[1] << " " << mktable[2] << " " << mktable[3] << endl;
        }
      }

      // +-----------------------+
      // | Codon site statistics |
      // +-----------------------+
      else if (cmdName == "CodonSiteStatistics")
      {
        if (!codonAlphabet) {
          throw Exception("CodonSiteStatstics can only be used with a codon alignment. Check the input alphabet!");
        }
        string path = ApplicationTools::getAFilePath("output.file", cmdArgs, true, false);
        if (path == "none") throw Exception("You must specify an ouptut file for CodonSiteStatistics"); 
        ApplicationTools::displayResult("Site statistics output to:", path);
        ofstream out(path.c_str(), ios::out);
        out << "Site\tMissingDataFrequency\tNbAlleles\tMinorAlleleFrequency\tMajorAlleleFrequency\tMinorAllele\tMajorAllele";
        out << "\tMeanNumberSynPos\tIsSynPoly\tIs4Degenerated\tPiN\tPiS";
        bool outgroup = (pscOut && pscOut->getNumberOfSequences() == 1);
        if (outgroup) {
          out << "\tOutgroupAllele";
        }
        if (estimateAncestor) {
          out << "\tAncestralAllele";
        }
        if (outgroup) {
          out << "\tMeanNumberSynPosDiv\tdN\tdS";
        }
        out << endl;

        unique_ptr<SiteContainer> sites(pscIn->toSiteContainer());
        for (size_t i = 0; i < sites->getNumberOfSites(); ++i) {
          const Site& site = sites->getSite(i);
          map<int, size_t> counts;
          SymbolListTools::getCounts(site, counts);
          size_t minFreq = site.size() + 1;
          size_t maxFreq = 0;
          int minState = -1;
          int maxState = -1;
          size_t nbAlleles = 0;
          size_t nbMissing = 0;
          for (map<int, size_t>::iterator it = counts.begin(); it != counts.end(); it++)
          {
            if (!alphabet->isUnresolved(it->first)
              && !alphabet->isGap(it->first)) {
              nbAlleles++; 
              if (it->second != 0) {
                if (it->second < minFreq) {
                  minFreq = it->second;
                  minState = it->first;        
                }
                if (it->second > maxFreq) {
                  maxFreq = it->second;
                  maxState = it->first;
                }
              }
            } else {
              nbMissing += it->second;            
            }
          }
          
          out << site.getPosition() << "\t";
          out << nbMissing << "\t";
          out << nbAlleles << "\t";
          out << minFreq << "\t";
          out << maxFreq << "\t";
          out << alphabet->intToChar(minState) << "\t";
          out << alphabet->intToChar(maxState) << "\t";
          if (estimateAncestor) {
            out << CodonSiteTools::numberOfSynonymousPositions(ancestralSequence->getValue(i), *gCode, kappa) << "\t";
          } else {
            out << CodonSiteTools::meanNumberOfSynonymousPositions(site, *gCode, kappa) << "\t";
          }
          out << CodonSiteTools::isSynonymousPolymorphic(site, *gCode) << "\t";
          out << CodonSiteTools::isFourFoldDegenerated(site, *gCode) << "\t";
          out << CodonSiteTools::piNonSynonymous(site, *gCode) << "\t";
          out << CodonSiteTools::piSynonymous(site, *gCode);
          if (outgroup) {
           out << "\t" << pscOut->getSequence(0).getChar(i); 
          }
          if (estimateAncestor) {
           out << "\t" << ancestralSequence->getChar(i); 
          }
          if (outgroup) {
            //Add divergence
            int outgroupState = pscOut->getSequence(0)[i];
            if (codonAlphabet->isUnresolved(outgroupState) || codonAlphabet->isGap(outgroupState)) {
              out << "\tNA\tNA\tNA";
            } else {
              //Average over outgroup (Note: minState and maxState are identical in this case)
              out << "\t" << (CodonSiteTools::numberOfSynonymousPositions(outgroupState, *gCode, kappa) +
                                     CodonSiteTools::numberOfSynonymousPositions(minState, *gCode, kappa)) / 2.;
              if (nbAlleles == 1) {
                //Compare with outgroup:
                if (site[0] == outgroupState) {
                  out << "\t0\t0";
                } else {
                  //This is a real substitution:
                  double nt = static_cast<double>(CodonSiteTools::numberOfDifferences(outgroupState, minState, *codonAlphabet));
                  double ns = CodonSiteTools::numberOfSynonymousDifferences(outgroupState, minState, *gCode); 
                  out << "\t" << (nt - ns) << "\t" << ns;
                }
              } else {
                //Site is polymorphic, this is not a substitution    
                out << "\t0\t0";
              }
            }
          }
          out << endl;
        }
      }
      
      else throw Exception("Unknown operation " + cmdName + ".");
    }
 
    // We're done!
    bpppopstats.done();
  }
  catch (exception& e)
  {
    if (logFile != "none")
      *cLog << "# Error: " << e.what() << endl;
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}

