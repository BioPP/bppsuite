// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

// From the STL:
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

// From bpp-core:
#include <Bpp/Version.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Text/KeyvalTools.h>
#include <Bpp/App/NumCalcApplicationTools.h>

// From bpp-seq:
#include <Bpp/Seq/SequenceTools.h>
#include <Bpp/Seq/Io/BppOAlignmentWriterFormat.h>

// From bpp-phy:
#include <Bpp/Phyl/App/BppPhylogeneticsApplication.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Simulation/SequenceSimulationTools.h>

#include <Bpp/Phyl/Simulation/EvolutionSequenceSimulator.h>
#include <Bpp/Phyl/Simulation/SimpleSubstitutionProcessSequenceSimulator.h>
#include <Bpp/Phyl/Simulation/GivenDataSubstitutionProcessSequenceSimulator.h>
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/PhyloLikelihoodContainer.h>
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/OneProcessSequencePhyloLikelihood.h>
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>

using namespace bpp;

int main(int args, char** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*            Bio++ Sequence Generator, version " << BPP_VERSION << "             *" << endl;
  cout << "*                                                                *" << endl;
  cout << "* Authors: J. Dutheil                                            *" << endl;
  cout << "*          B. Boussau                       Last Modif. " << BPP_REL_DATE << " *" << endl;
  cout << "*          L. Gueguen                                            *" << endl;
  cout << "*          M. Groussin                                           *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  try
  {
    BppPhylogeneticsApplication bppseqgen(args, argv, "bppseqgen");

    if (args == 1)
    {
      bppseqgen.help("bppSeqGen");
      return 0;
    }

    bppseqgen.startTimer();
    map<string, string> unparsedParams;

    Context context;

    std::shared_ptr<const Alphabet> alphabet = bppseqgen.getAlphabet();

    auto gCode = bppseqgen.getGeneticCode(alphabet);

    // Write to file:
    BppOAlignmentWriterFormat bppoWriter(1);

    ////// Get the optional map of the sequences

    auto mSitesuniq = bppseqgen.getConstAlignmentsMap(alphabet, true, true);

    auto mSites = PhylogeneticsApplicationTools::uniqueToSharedMap<const AlignmentDataInterface>(mSitesuniq);

    /// collection

    std::shared_ptr<SubstitutionProcessCollection> spc = bppseqgen.getCollection(alphabet, gCode, mSites, unparsedParams);

    auto mSeqEvoluniq = bppseqgen.getProcesses(spc, unparsedParams);

    auto mSeqEvol = PhylogeneticsApplicationTools::uniqueToSharedMap<SequenceEvolution>(mSeqEvoluniq);

    /// Get optional phylolikelihoods (in case of posterior simulation)


    bppseqgen.setWarningLevel(0);
    std::shared_ptr<PhyloLikelihoodContainer> phyloCont =  bppseqgen.getPhyloLikelihoods(context, mSeqEvol, spc, mSites, "");

    bppseqgen.setWarningLevel(1);

    /// ///////////////////
    // Simulators
    vector<string> vSimulName = ApplicationTools::matchingParameters("simul*", bppseqgen.getParams());

    if (vSimulName.size() == 0)
    {
      ApplicationTools::displayWarning("Did not find any descriptor matching `simul*`, so no simulation performed.");
    }

    for (size_t nS = 0; nS < vSimulName.size(); nS++)
    {
      size_t poseq = vSimulName[nS].find("=");
      string suff = vSimulName[nS].substr(5, poseq - 5);

      size_t num = static_cast<size_t>(TextTools::toInt(suff));

      string simulDesc = ApplicationTools::getStringParameter(vSimulName[nS], bppseqgen.getParams(), "", "", true, true);

      map<string, string> argsim;
      string simulName;

      KeyvalTools::parseProcedure(simulDesc, simulName, argsim);

      ApplicationTools::displayMessage("");
      ApplicationTools::displayMessage("Simulation " + TextTools::toString(num));

      //////////////////////

      if (argsim.find("process") == argsim.end() && (argsim.find("phylo") == argsim.end()))
        throw BadIntegerException("bppseqgen. Missing process or phylo argument for simul:", (int)num);


      // Root states
      vector<size_t> states;
      vector<double> rates;
      bool withStates = false;
      bool withRates = false;
      size_t nbSites = 0;

      auto sm = spc->getModel(1)->getStateMap();

      // Data or info at root
      if (argsim.find("data") != argsim.end())
      {
        size_t indData = (size_t)ApplicationTools::getIntParameter("data", argsim, 1, "", true, 0);

        if (mSites.find(indData) == mSites.end())
          throw BadIntegerException("bppseqgen : Unknown data number:", (int)indData);

        auto data = mSites.at(indData);

        auto nseq = RandomTools::giveIntRandomNumberBetweenZeroAndEntry(data->getNumberOfSequences());

        ApplicationTools::displayResult("Using root sequence", data->sequence(nseq).getName());

        string siteSet = ApplicationTools::getStringParameter("input.site.selection", argsim, "none", "", true, 1);

        nbSites = data->getNumberOfSites();

        if (siteSet != "none")
        {
          vector<size_t> vSite;
          try
          {
            vector<int> vSite1 = NumCalcApplicationTools::seqFromString(siteSet, ",", ":");
            for (size_t i = 0; i < vSite1.size(); ++i)
            {
              int x = (vSite1[i] >= 0 ? vSite1[i] : static_cast<int>(nbSites) + vSite1[i]);
              if (x >= 0)
                vSite.push_back(static_cast<size_t>(x));
              else
                throw Exception("bppseqgen. Incorrect negative index for site selection: " + TextTools::toString(x));
            }
          }
          catch (Exception& e)
          {
            string seln;
            map<string, string> selArgs;
            KeyvalTools::parseProcedure(siteSet, seln, selArgs);
            if (seln == "Sample")
            {
              size_t n = ApplicationTools::getParameter<size_t>("n", selArgs, nbSites, "", true, 1);
              bool replace = ApplicationTools::getBooleanParameter("replace", selArgs, false, "", true, 1);

              vSite.resize(n);
              vector<size_t> vPos;
              for (size_t p = 0; p < nbSites; ++p)
              {
                vPos.push_back(p);
              }

              RandomTools::getSample(vPos, vSite, replace);
            }
          }

          nbSites = vSite.size();
        }

        ApplicationTools::displayResult("Number of sites", nbSites);
        states.resize(nbSites);

        withStates = true;

        size_t nbStates = alphabet->getSize();
        std::vector<double> probstate(nbStates);

        auto resChar = alphabet->getResolvedChars();

        for (size_t i = 0; i < nbSites; ++i)
        {
          for (size_t j = 0; j < nbStates; j++)
          {
            probstate[j] = data->getStateValueAt(i, nseq, alphabet->getIntCodeAt(j + 1));
          }

          auto pchar = RandomTools::pickOne<string>(resChar, probstate, true);
          states[i] = RandomTools::pickOne<size_t>(sm->getModelStates(pchar));
        }
      }
      else
      {
        string infosFile = ApplicationTools::getAFilePath("input.infos", argsim, false, true);

        if (infosFile != "none")
        {
          ApplicationTools::displayResult("Site information", infosFile);
          ifstream in(infosFile.c_str());
          auto infos = DataTable::read(in, "\t");
          nbSites = infos->getNumberOfRows();
//      ApplicationTools::displayResult("Number of sites", TextTools::toString(nbSites));
          string rateCol = ApplicationTools::getStringParameter("input.infos.rates", argsim, "pr", "", true, true);
          string stateCol = ApplicationTools::getStringParameter("input.infos.states", argsim, "none", "", true, true);
          withRates = rateCol != "none";
          withStates = stateCol != "none";


          // Specific input data
          if (withRates)
          {
            rates.resize(nbSites);
            vector<string> ratesStrings = infos->getColumn(rateCol);
            for (size_t i = 0; i < nbSites; i++)
            {
              rates[i] = TextTools::toDouble(ratesStrings[i]);
            }
          }
          if (withStates)
          {
            vector<string> ancestralStates = infos->getColumn(stateCol);

            states.resize(nbSites);
            for (size_t i = 0; i < nbSites; i++)
            {
              int alphabetState = alphabet->charToInt(ancestralStates[i]);
              // If a generic character is provided, we pick one state randomly from the possible ones:
              if (alphabet->isUnresolved(alphabetState))
                alphabetState = RandomTools::pickOne<int>(alphabet->getAlias(alphabetState));
              states[i] = RandomTools::pickOne<size_t>(sm->getModelStates(alphabetState));
            }

            string siteSet = ApplicationTools::getStringParameter("input.site.selection", argsim, "none", "", true, 1);

            if (siteSet != "none")
            {
              vector<size_t> vSite;
              try
              {
                vector<int> vSite1 = NumCalcApplicationTools::seqFromString(siteSet, ",", ":");
                for (size_t i = 0; i < vSite1.size(); ++i)
                {
                  int x = (vSite1[i] >= 0 ? vSite1[i] : static_cast<int>(nbSites) + vSite1[i]);
                  if (x >= 0)
                    vSite.push_back(static_cast<size_t>(x));
                  else
                    throw Exception("bppseqgen : incorrect negative index: " + TextTools::toString(x));
                }
              }
              catch (Exception& e)
              {
                string seln;
                map<string, string> selArgs;
                KeyvalTools::parseProcedure(siteSet, seln, selArgs);
                if (seln == "Sample")
                {
                  size_t n = ApplicationTools::getParameter<size_t>("n", selArgs, nbSites, "", true, 1);
                  bool replace = ApplicationTools::getBooleanParameter("replace", selArgs, false, "", true, 1);

                  vSite.resize(n);
                  vector<size_t> vPos;
                  for (size_t p = 0; p < nbSites; ++p)
                  {
                    vPos.push_back(p);
                  }

                  RandomTools::getSample(vPos, vSite, replace);
                }
              }

              nbSites = vSite.size();

              vector<size_t> newStates(nbSites);
              vector<double> newRates(nbSites);

              for (size_t ni = 0; ni < nbSites; ++ni)
              {
                newStates[ni] = states[vSite[ni]];
                newRates[ni]  = rates[vSite[ni]];
              }

              states = newStates;
              rates = newRates;
            }
          }
        }
      }

      ////////////////////////////////////
      /////// Process

      unique_ptr<SequenceSimulatorInterface> ss;

      // for null nodes in posterior phylo simulations
      vector<unsigned int> nodesId(0);


      if (argsim.find("process") != argsim.end())
      {
        size_t indProcess = (size_t)ApplicationTools::getIntParameter("process", argsim, 1, "", true, 0);

        ApplicationTools::displayResult(" Process", TextTools::toString(indProcess));

        if (!spc->hasSubstitutionProcessNumber(indProcess))
        {
          if (mSeqEvol.find(indProcess) == mSeqEvol.end())

            throw BadIntegerException("bppseqgen. Unknown process number:", (int)indProcess);

          ss = make_unique<EvolutionSequenceSimulator>(*mSeqEvol.find(indProcess)->second);
        }
        else
        {
          ss = make_unique<SimpleSubstitutionProcessSequenceSimulator>(spc->getSubstitutionProcess(indProcess));
        }
      }
      else
      {
        size_t indPhylo = (size_t)ApplicationTools::getIntParameter("phylo", argsim, 1, "", true, 0);

        if (!phyloCont)
          throw BadIntegerException("bppseqgen. Empty phylocontainer for simul:", (int)num);

        auto phylo = phyloCont->getPhyloLikelihood(indPhylo);

        ApplicationTools::displayResult(" Phylolikelihood", TextTools::toString(indPhylo));

        if (!phylo)
          throw BadIntegerException("bppseqgen. Unknown phylo number:", (int)indPhylo);

        auto spph = dynamic_pointer_cast<SingleProcessPhyloLikelihood>(phylo);
        auto opsp = dynamic_pointer_cast<OneProcessSequencePhyloLikelihood>(phylo);

        if (!spph && !opsp)
          throw BadIntegerException("bppseqgen. Posterior simulation not implemented for this kind of phylolikelihood. Ask developers.", (int)num);

        std::shared_ptr<LikelihoodCalculationSingleProcess> lcsp = spph ? spph->getLikelihoodCalculationSingleProcess() :
          opsp->getLikelihoodCalculationSingleProcess();

        // Get nodes with no phyloLik
        auto descnodes = ApplicationTools::getStringParameter("nullnodes", argsim, "", "", true, 0);

        if (descnodes!="")
        {
          auto tree= lcsp->substitutionProcess().getParametrizablePhyloTree();
        
          if (descnodes == "All")
          {
            nodesId = tree->getEdgeIndexes(tree->getSubtreeEdges(tree->getRoot()));
          }
          else if (descnodes == "Leaves")
          {
            nodesId = tree->getNodeIndexes(tree->getLeavesUnderNode(tree->getRoot()));
          }
          else if (descnodes == "NoLeaves")
          {
            auto allIds = tree->getEdgeIndexes(tree->getSubtreeEdges(tree->getRoot()));
            auto leavesId = tree->getNodeIndexes(tree->getLeavesUnderNode(tree->getRoot()));
            VectorTools::diff(allIds, leavesId, nodesId);
          }
          else
            nodesId = ApplicationTools::getVectorParameter<unsigned int>("nullnodes", argsim, ',', ':', "", "", true, 1);
        }

        
        ///////////////: Specific sites? 
        if (argsim.find("pos") == argsim.end())// Sequence simulation similar to the data, number_of_sites will not be used
          ss = make_unique<GivenDataSubstitutionProcessSequenceSimulator>(lcsp, nodesId);
        else
        {
          size_t pos = (size_t)ApplicationTools::getIntParameter("pos", argsim, 1, "", true, 0);
          ApplicationTools::displayResult(" Position", TextTools::toString(pos));
          
          ss = make_unique<SimpleSubstitutionProcessSequenceSimulator>(lcsp, pos, false, nodesId);
        }
      }

      // output

      string mfnames = ApplicationTools::getAFilePath("output.sequence.file", argsim, true, false, "", false, "none", true);
      ApplicationTools::displayResult(" Output file", mfnames);

      string mformats = ApplicationTools::getStringParameter("output.sequence.format", argsim, "Fasta", "", false, true);
      ApplicationTools::displayResult(" Output format", mformats);

      auto mintern = ApplicationTools::getBooleanParameter("output.internal.sequences", argsim, false, "", true, 1);
      ApplicationTools::displayBooleanResult(" Output internal", mintern);

      if (nodesId.size()!=0)
        ApplicationTools::displayResult(" No phylo nodes", VectorTools::paste(nodesId, ", "));
        
      auto gds = dynamic_cast<GivenDataSubstitutionProcessSequenceSimulator*>(ss.get());
      auto pps = dynamic_cast<SubstitutionProcessSequenceSimulator*>(ss.get());

      size_t nbmin = gds ? gds->getNumberOfSites() : pps ? pps->getNumberOfSites() : 100;

      nbSites = (size_t)ApplicationTools::getIntParameter("number_of_sites", argsim, (int)nbmin, "", false, 0);
      ApplicationTools::displayResult(" Number of sites", TextTools::toString(nbSites));


      ss->outputInternalSequences(mintern);
      std::shared_ptr<SiteContainerInterface> sites = 0;

      if (withStates || withRates)
      {
        auto pss = dynamic_cast<SubstitutionProcessSequenceSimulator*>(ss.get());

        if (withStates)
          if (withRates)
            sites = pss ? pss->simulate(rates, states) : SequenceSimulationTools::simulateSites(*ss, rates, states);
          else
            sites = pss ? pss->simulate(states) : SequenceSimulationTools::simulateSites(*ss, states);
        else
          sites = pss ? pss->simulate(rates) : SequenceSimulationTools::simulateSites(*ss, rates);

        ApplicationTools::displayTaskDone();
      }
      else
      {
        ApplicationTools::displayMessage("");
        ApplicationTools::displayTask("Perform simulations");

        sites = ss->simulate(nbSites);
        ApplicationTools::displayTaskDone();
        ApplicationTools::displayMessage("");
      }

      unique_ptr<OAlignment> oAln(bppoWriter.read(mformats));
      // ApplicationTools::displayResult("Output alignment file ", filenames[it.first]);
      // ApplicationTools::displayResult("Output alignment format ", oAln->getFormatName());

      oAln->writeAlignment(mfnames, *sites, true);
    }

    bppseqgen.done();
  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}
