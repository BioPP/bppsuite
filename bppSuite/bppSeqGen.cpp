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

// Classes to postpone generation in simulation part, and not in reading part


template<class T>
class SampleClass 
{
private:
  std::vector<T> states_;
  bool replace_;
  bool ordered_; // fix the sampling
  
public:
  SampleClass() : states_(), replace_(false), ordered_(false){};

  SampleClass(const SampleClass<T>& po) : states_(po.states_), replace_(po.replace_), ordered_(po.ordered_) {};

  SampleClass<T>& operator=(const SampleClass<T>& po)
  {
    states_=po.states_;
    replace_=po.replace_;
    ordered_=po.ordered_;
    return *this;
  }

  void setStates(const std::vector<T>& states)
  {
    states_=states;
  }

  void setStates(const T& state)
  {
    states_=std::vector<T>(1,state);
  }

  const std::vector<T> getStates() const
  {
    return states_;
  }

  void setOrdered(bool ordered)
  {
    ordered_=ordered;
  }

  bool getOrdered() const
  {
    return ordered_;
  }

  void setReplace(bool replace)
  {
    replace_=replace;
  }

  bool getReplace() const
  {
    return replace_;
  }

  T pickOne() const
  {
    if (states_.size()==1)
      return states_[0];
    else
      return RandomTools::pickOne<T>(states_);
  }

  /*
   *@brief extract a vector of a given length
   *
   * If ordered_, returns size first elements, otherwise samples from
   * it.
   */

  std::vector<T> pickSeveral(size_t size) const
  {
    std::vector<T> vS(size);

    if (!ordered_)
      RandomTools::getSample(states_, vS, replace_);
    else
    {
      if (size>states_.size())
        throw BadSizeException("pickSeveral: too many sites to extract.",size,states_.size());
      
      vS.insert(vS.begin(), states_.begin(), states_.begin() + static_cast<bpp::ModelPath::PathNode::difference_type>(size));
    }
    
    return vS;
  }

};


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

    size_t nbStates = alphabet->getSize();
    auto resChar = alphabet->getResolvedChars();

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

      // output

      string mfnames = ApplicationTools::getAFilePath("output.sequence.file", argsim, true, false, "", false, "none", true);

      string mformats = ApplicationTools::getStringParameter("output.sequence.format", argsim, "Fasta", "", false, true);
      ApplicationTools::displayResult(" Output format", mformats);

      auto mintern = ApplicationTools::getBooleanParameter("output.internal.sequences", argsim, false, "", true, 1);
      ApplicationTools::displayBooleanResult(" Output internal", mintern);

      
      //////////////////////////////////////////////////////
      /////// Process
      size_t nbSites = 0;

      unique_ptr<SequenceSimulatorInterface> ss;

      // for null nodes in posterior phylo simulations
      vector<unsigned int> nodesId(0);


      if (argsim.find("process") != argsim.end())
      {
        size_t indProcess = (size_t)ApplicationTools::getIntParameter("process", argsim, 1, "", true, 0);

        ApplicationTools::displayResult(" Process", TextTools::toString(indProcess));

        if (!spc->hasSubstitutionProcessNumber(indProcess))  // Sequence process
        {
          if (mSeqEvol.find(indProcess) == mSeqEvol.end())

            throw BadIntegerException("bppseqgen. Unknown process number:", (int)indProcess);

          ss = make_unique<EvolutionSequenceSimulator>(*mSeqEvol.find(indProcess)->second);
        } 
        else // Site Process
        {
          ss = make_unique<SimpleSubstitutionProcessSequenceSimulator>(spc->getSubstitutionProcess(indProcess));
        }
      }
      else // Phylo process
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

        
        ///////////////: Simulation of a specific site
        if (argsim.find("pos") == argsim.end())// Sequence simulation similar to the data, number_of_sites will not be used
          ss = make_unique<GivenDataSubstitutionProcessSequenceSimulator>(lcsp, nodesId);
        else
        {
          size_t pos = (size_t)ApplicationTools::getIntParameter("pos", argsim, 1, "", true, 0);
          ApplicationTools::displayResult(" Position", TextTools::toString(pos));
          
          ss = make_unique<SimpleSubstitutionProcessSequenceSimulator>(lcsp, pos, false, nodesId);
        }
      }

      ss->outputInternalSequences(mintern);

      if (nodesId.size()!=0)
        ApplicationTools::displayResult(" No phylo nodes", VectorTools::paste(nodesId, ", "));

      
      auto gds = dynamic_cast<GivenDataSubstitutionProcessSequenceSimulator*>(ss.get());  // Constrained by data
      auto pps = dynamic_cast<SubstitutionProcessSequenceSimulator*>(ss.get());// Constrained by sequence structure

      // Number of sites from process
      nbSites = gds ? gds->getNumberOfSites() : pps ? pps->getNumberOfSites() : 0;

      // Specified number of sites
      size_t nbmin = (size_t)ApplicationTools::getIntParameter("number_of_sites", argsim, 0, "", false, 2);
      if (nbSites==0 || nbmin<nbSites)
        nbSites=nbmin;
      
      ////////////////////////////////////////////////
      // Root states
      vector<SampleClass<size_t>> states;
      vector<double> rates;
      bool withStates = false; // using states at root
      bool withRates = false;  // using rates

      auto sm = spc->getModel(1)->getStateMap();

      // Data or info at root
      auto rootData = ApplicationTools::getStringParameter("root.data", argsim, "", "", true, 0);
      auto withData = (rootData!="");

      std::shared_ptr<const AlignmentDataInterface> data;
      bool sampleseq = false;
      SampleClass<uint> nseq;

      string infosFile = ApplicationTools::getAFilePath("input.infos", argsim, false, false,"",1,"",1);

      // withData
      if (withData)
      {
        withStates=true;
        
        auto paro = rootData.find("(");
        if (paro==std::string::npos)
          throw Exception("Bad syntax for root.data: need SequenceFrom(<int>) or SitesFrom(<int>)");
        auto pref= rootData.substr(0,paro);
        auto parc = rootData.find(")",paro);
        unsigned int indData = TextTools::to<unsigned int>(rootData.substr(paro+1, parc-paro-1));
        
        if (mSites.find(indData) == mSites.end())
          throw BadIntegerException("bppseqgen : Unknown data number:", (int)indData);
        
        data = mSites.at(indData);
        if (nbSites==0 || data->getNumberOfSites()<nbSites)
          nbSites=data->getNumberOfSites();
        
        sampleseq=(pref=="SequenceFrom");
        Vuint vseq(data->getNumberOfSequences());
        std::iota(vseq.begin(), vseq.end(), 0);
          
        nseq.setStates(vseq);
      }

      // end withData
      
      /// Info file      
      if (infosFile != "none")
      {
        ApplicationTools::displayResult("Site information", infosFile);
        ifstream in(infosFile.c_str());
        auto infos = DataTable::read(in, "\t");

        if (nbSites==0 || (infos->getNumberOfRows()< nbSites))
          nbSites = infos->getNumberOfRows();

        string rateCol = ApplicationTools::getStringParameter("input.infos.rates", argsim, "none", "", true, true);
        string stateCol = ApplicationTools::getStringParameter("input.infos.states", argsim, "none", "", true, true);

        withRates = rateCol != "none";

        // Specific input files
        if (withRates)
        {
          rates.resize(nbSites);
          vector<string> ratesStrings = infos->getColumn(rateCol);
          for (size_t i = 0; i < nbSites; i++)
          {
            rates[i] = TextTools::toDouble(ratesStrings[i]);
          }
        }
        
        if (!withData) // if no Data already seen
        {
          withStates = (stateCol != "none");
          if (withStates)
          {
            vector<string> ancestralStates = infos->getColumn(stateCol);
            
            states.resize(nbSites);
            for (size_t i = 0; i < nbSites; i++)
            {
              vector<size_t> vstates;
              int alphabetState = alphabet->charToInt(ancestralStates[i]);
              // If a generic character is provided, we pick one state randomly from the possible ones:
              if (alphabet->isUnresolved(alphabetState))
                for (const auto& vs : alphabet->getAlias(alphabetState))
                  for (const auto& j : sm->getModelStates(vs))
                    vstates.push_back(j);
              else
                vstates=sm->getModelStates(alphabetState);
              
              states[i].setStates(vstates);
            }
          }
        }
      }

      // Select sites
      string siteSet = ApplicationTools::getStringParameter("input.site.selection", argsim, "none", "", true, 1);
      if (siteSet[0] == '(')
        siteSet = siteSet.substr(1, siteSet.size() - 2);

      SampleClass<size_t> vSite;
      
      if (siteSet != "none")
      {
        try
        {
          vector<int> vSite1 = NumCalcApplicationTools::seqFromString(siteSet, ",", ":");
          vector<size_t> vS;
          for (size_t i = 0; i < vSite1.size(); ++i)
          {
            int x = (vSite1[i] >= 0 ? vSite1[i] : static_cast<int>(nbSites) + vSite1[i]);
            if (x >= 0)
              if ((nbSites!=0) && (x>=(int)nbSites)) // because up to now sites are contiguous
                throw BadSizeException("bppseqgen. Incorrect too large index for site selection: " , nbSites, static_cast<size_t>(x));
              else                
                vS.push_back((size_t)x);
            else
              throw Exception("bppseqgen. Incorrect negative index for site selection: " + TextTools::toString(x));
          }
          vSite.setOrdered(true);
          vSite.setStates(vS);
          nbSites = vS.size();
        }
        catch (Exception& e)
        {
          if (nbSites==0)
            nbSites=100;
          string seln;
          map<string, string> selArgs;
          KeyvalTools::parseProcedure(siteSet, seln, selArgs);
          if (seln == "Sample")
          {
            size_t n = ApplicationTools::getParameter<size_t>("n", selArgs, nbSites, "", true, 1);
            if (n==0)
              throw Exception("Missing simul length for Sample in  input.site.selection " + siteSet);
            bool replace = ApplicationTools::getBooleanParameter("replace", selArgs, false, "", true, 1);
            
            vector<size_t> vPos(nbSites);
            std::iota(vPos.begin(), vPos.end(), 0);
            
            vSite.setReplace(replace);
            vSite.setStates(vPos);
            nbSites = n;
          }
          else
            throw Exception("Unknown description of input.site.selection " + siteSet);
        }
      }
      else
      {
        if (nbSites==0)
          nbSites=100;

        vector<size_t> vPos(nbSites);
        std::iota(vPos.begin(), vPos.end(), 0);
          
        vSite.setOrdered(true);
        vSite.setStates(vPos);
      }
      
      ApplicationTools::displayResult(" Number of sites", nbSites);
      
      // Type of simulation

      size_t nbsimul=1;
      if (simulName=="Multiple")
        nbsimul = (size_t)ApplicationTools::getIntParameter("n", argsim, 1, "", true, 0);

      unique_ptr<OAlignment> oAln(bppoWriter.read(mformats));

      for (size_t isimul=0; isimul<nbsimul; isimul++)
      {
        vector<size_t> newStates(nbSites);
        vector<double> newRates(nbSites);

        auto vsitesR = vSite.pickSeveral(nbSites); // Take sites 
        
        if (withData)
        {
          size_t iseq=0;
          
          if (sampleseq)
          {
            iseq = nseq.pickOne();
            ApplicationTools::displayResult("Using root sequence", data->sequence(iseq).getName());
          }
          else
            ApplicationTools::displayMessage(" Using root sequence from site-sampling in alignment ");

          std::vector<double> probstate(nbStates);
          for (size_t i = 0; i < nbSites; ++i)
          {
            if (!sampleseq)
              iseq = nseq.pickOne();
            
            for (size_t j = 0; j < nbStates; j++)
              probstate[j] = data->getStateValueAt(vsitesR[i], iseq, alphabet->getIntCodeAt(j + 1));
            
            string pchar;
            if (VectorTools::sum(probstate)>0.001)
              pchar = RandomTools::pickOne<string>(resChar, probstate, true);
            else
              pchar = alphabet->intToChar(alphabet->getGapCharacterCode());
            
            newStates[i] = RandomTools::pickOne<size_t>(sm->getModelStates(pchar));
            newRates[i] = withRates?rates[vsitesR[i]]:1;
          }
        }
        
        std::shared_ptr<SiteContainerInterface> sites = 0;

        ApplicationTools::displayTask(" Perform simulations " + TextTools::toString(nS+1) + (nbsimul>1?"_"+TextTools::toString(isimul+1):""));
        ApplicationTools::displayTaskDone();

        if (withRates || withStates)
        {
          if (withStates)
            if (withRates)
              sites = pps ? pps->simulate(newRates, newStates) : SequenceSimulationTools::simulateSites(*ss, newRates, newStates);
            else
              sites = pps ? pps->simulate(newStates) : SequenceSimulationTools::simulateSites(*ss, newStates);
          else
            sites = pps ? pps->simulate(newRates) : SequenceSimulationTools::simulateSites(*ss, newRates);
        }
        else
          sites = ss->simulate(nbSites);
        
        string outnf = (nbsimul > 1 ? mfnames + "_" + TextTools::toString(isimul+1):mfnames);
        
        ApplicationTools::displayResult(" Output file", outnf);
        ApplicationTools::displayMessage("");

        oAln->writeAlignment(outnf, *sites, true);
      }
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
