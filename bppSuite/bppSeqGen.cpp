//
// File: bppSeqGen.cpp
// Created by: Julien Dutheil
// Created on: Oct Mon 24 18:50 2005
//

/*
  Copyright or � or Copr. Bio++ Development Team

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
#include <Bpp/Phyl/Simulation/SequenceSimulationTools.h>
#include <Bpp/Phyl/Simulation/EvolutionSequenceSimulator.h>
#include <Bpp/Phyl/Simulation/SimpleSubstitutionProcessSequenceSimulator.h>
#include <Bpp/Phyl/Simulation/GivenDataSubstitutionProcessSequenceSimulator.h>
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/PhyloLikelihoodContainer.h>
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/OneProcessSequencePhyloLikelihood.h>
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>


using namespace bpp;


map<size_t, SequenceSimulator*> readSimul(const SubstitutionProcessCollection& spc,
                                          const map<size_t, SequenceEvolution*>& mSeqEvol,
                                          std::shared_ptr<PhyloLikelihoodContainer> phyloCont,
                                          const map<string, string>& params,
                                          map<size_t, string>& mfnames,
                                          map<size_t, string>& mformats,
                                          map<size_t, bool>& mintern,
                                          map<size_t, size_t>& mlength)
{

  map<size_t, SequenceSimulator*> mSim;

  vector<string> vSimulName=ApplicationTools::matchingParameters("simul*", params);

  if (vSimulName.size() == 0) {
    ApplicationTools::displayWarning("Did not find any descriptor matching `simul*`, so no simulation performed.");
  }
  
  SequenceSimulator* ss;
  
  for (size_t nS=0; nS< vSimulName.size(); nS++)
  {
    size_t poseq=vSimulName[nS].find("=");
    string suff = vSimulName[nS].substr(5,poseq-5);

    size_t num=static_cast<size_t>(TextTools::toInt(suff));

    string simulDesc=ApplicationTools::getStringParameter(vSimulName[nS], params, "", "", true, true);

    map<string, string> args;
    string simulName;

    KeyvalTools::parseProcedure(simulDesc, simulName, args);

    // Process
    ApplicationTools::displayMessage("");
    ApplicationTools::displayMessage("Simulation " + TextTools::toString(num));

    if (args.find("process")==args.end() && (args.find("phylo")==args.end()))
      throw BadIntegerException("bppseqgen::readSimul. Missing process and phylo argument for simul:",(int)num);

    if (args.find("process")!=args.end())
    {
      size_t indProcess=(size_t)ApplicationTools::getIntParameter("process", args, 1, "", true, 0);

      ApplicationTools::displayResult(" Process", TextTools::toString(indProcess));

      if (! spc.hasSubstitutionProcessNumber(indProcess))
      {
        if (mSeqEvol.find(indProcess)==mSeqEvol.end())
          
          throw BadIntegerException("bppseqgen::readSimul. Unknown process number:",(int)indProcess);
        
        ss= new EvolutionSequenceSimulator(*mSeqEvol.find(indProcess)->second);
      }
      else
      {
        ss=new SimpleSubstitutionProcessSequenceSimulator(spc.getSubstitutionProcess(indProcess));
      }
    }
    else
    {
      size_t indPhylo=(size_t)ApplicationTools::getIntParameter("phylo", args, 1, "", true, 0);

      if (!phyloCont)
        throw BadIntegerException("bppseqgen::readSimul. Empty phylocontainer for simul:",(int)num);

      auto phylo = phyloCont->getPhyloLikelihood(indPhylo);

      ApplicationTools::displayResult(" Phylolikelihood", TextTools::toString(indPhylo));

      if (!phylo)
        throw BadIntegerException("bppseqgen::readSimul. Unknown phylo number:",(int)indPhylo);

      auto spph = dynamic_pointer_cast<SingleProcessPhyloLikelihood>(phylo);
      auto opsp = dynamic_pointer_cast<OneProcessSequencePhyloLikelihood>(phylo);

      if (!spph && !opsp)
        throw BadIntegerException("bppseqgen::readSimul. Posterior simulation not implemented for this kind of phylolikelihood. Ask developpers.",(int)num);

      std::shared_ptr<LikelihoodCalculationSingleProcess> lcsp = spph?spph->getLikelihoodCalculationSingleProcess():
        opsp->getLikelihoodCalculationSingleProcess();

      if (args.find("pos")==args.end()) // Sequence simulation similar to the data, number_of_sites will not be used
        ss=new GivenDataSubstitutionProcessSequenceSimulator(lcsp);
      else
      {        
        size_t pos=(size_t)ApplicationTools::getIntParameter("pos", args, 1, "", true, 0);
        ApplicationTools::displayResult(" Position", TextTools::toString(pos));

        ss=new SimpleSubstitutionProcessSequenceSimulator(lcsp, pos);
      }
    }
    
    // output

    mfnames[num] = ApplicationTools::getAFilePath("output.sequence.file", args, true, false, "", false, "none", true);
    ApplicationTools::displayResult(" Output file", mfnames[num]);
    
    mformats[num] = ApplicationTools::getStringParameter("output.sequence.format", args, "Fasta", "", false, true);
    ApplicationTools::displayResult(" Output format", mformats[num]);

    mintern[num] = ApplicationTools::getBooleanParameter("output.internal.sequences", args, false, "", true, 1);
    ApplicationTools::displayBooleanResult(" Output internal", mintern[num]);


    if (dynamic_cast<GivenDataSubstitutionProcessSequenceSimulator*>(ss)==0)
    {
      mlength[num] = (size_t)ApplicationTools::getIntParameter("number_of_sites", args, 100, "", false, 0);
      ApplicationTools::displayResult(" Number of sites", TextTools::toString(mlength[num]));
    }

    mSim[num]=ss;
  }

  
  return mSim;
}


int main(int args, char ** argv)
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

  try {

    BppPhylogeneticsApplication bppseqgen(args, argv, "bppseqgen");

    if(args == 1)
    {
      bppseqgen.help("bppSeqGen");
      return 0;
    }

    bppseqgen.startTimer();
    map<string, string> unparsedParams;

    Context context;

    Alphabet* alphabet = bppseqgen.getAlphabet();

    unique_ptr<GeneticCode> gCode(bppseqgen.getGeneticCode(alphabet));

    ////// Get the optional map of the sequences
    
    map<size_t, AlignedValuesContainer*> mSites = bppseqgen.getAlignmentsMap(alphabet, true, true);

    /// collection
    
    SubstitutionProcessCollection* SPC = bppseqgen.getCollection(alphabet, gCode.get(), mSites, unparsedParams);

    map<size_t, SequenceEvolution*> mSeqEvol = bppseqgen.getProcesses(*SPC, unparsedParams);

    /// Get optional phylolikelihoods (in case of posterior simulation)

    auto phyloCont =  bppseqgen.getPhyloLikelihoods(context, mSeqEvol, *SPC, mSites, "", 0);

    /*******************************************/
    /*     Starting sequence                   */
    /*******************************************/

    std::shared_ptr<SiteContainer> sites = 0;
    size_t nbSites = 0;

    bool withStates = false;
    bool withRates = false;
    vector<size_t> states;
    vector<double> rates;

    const StateMap& sm=SPC->getModel(1)->getStateMap();
    
    string infosFile = ApplicationTools::getAFilePath("input.infos", bppseqgen.getParams(), false, true);

    if (infosFile != "none")
    {
      ApplicationTools::displayResult("Site information", infosFile);
      ifstream in(infosFile.c_str());
      DataTable* infos = DataTable::read(in, "\t");
      nbSites = infos->getNumberOfRows();
//      ApplicationTools::displayResult("Number of sites", TextTools::toString(nbSites));
      string rateCol = ApplicationTools::getStringParameter("input.infos.rates", bppseqgen.getParams(), "pr", "", true, true);
      string stateCol = ApplicationTools::getStringParameter("input.infos.states", bppseqgen.getParams(), "none", "", true, true);
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
          //If a generic character is provided, we pick one state randomly from the possible ones:
          if (alphabet->isUnresolved(alphabetState))
            alphabetState = RandomTools::pickOne<int>(alphabet->getAlias(alphabetState));
          states[i] = RandomTools::pickOne<size_t>(sm.getModelStates(alphabetState));
        }

        string siteSet = ApplicationTools::getStringParameter("input.site.selection", bppseqgen.getParams(), "none", "", true, 1);

        if (siteSet != "none")
        {
          vector<size_t> vSite;
          try {
            vector<int> vSite1 = NumCalcApplicationTools::seqFromString(siteSet);
            for (size_t i = 0; i < vSite1.size(); ++i){
              int x = (vSite1[i] >= 0 ? vSite1[i] : static_cast<int>(nbSites) + vSite1[i]);
              if (x >= 0)
                vSite.push_back(static_cast<size_t>(x));
              else
                throw Exception("SequenceApplicationTools::getSiteContainer(). Incorrect negative index: " + TextTools::toString(x));
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
                vPos.push_back(p);

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
    else
    {
      try {
        VectorSiteContainer* allSeq = 0;
        
        allSeq = SequenceApplicationTools::getSiteContainer(alphabet, bppseqgen.getParams());
        
        if (allSeq->getNumberOfSequences() > 0)
        {
          Sequence* pseq = SequenceTools::getSequenceWithCompleteSites(allSeq->getSequence(0));

          nbSites = pseq->size();
          states.resize(nbSites);
          withStates = true;

          for (size_t i = 0; i < nbSites; ++i) {
            states[i] = RandomTools::pickOne<size_t>(sm.getModelStates((*pseq)[i]));
          }
//          ApplicationTools::displayResult("Number of sites", TextTools::toString(nbSites));

          delete pseq;
        }
      }
      catch (Exception& e)
      {
        // cout << e.what() << endl;
        // return 1;
      }
    }

  
    /*******************/
    /* Simulations     */
    /*******************/

    map<size_t, string> filenames;
    map<size_t, string> formats;
    map<size_t, bool> internal;
    map<size_t, size_t> lengths;

    auto mSim=readSimul(*SPC, mSeqEvol, phyloCont, bppseqgen.getParams(),filenames, formats, internal, lengths);

    for (auto& it : mSim)
    {
      SequenceSimulator& seqsim = *it.second;
      seqsim.outputInternalSequences(internal[it.first]);
      
      if (withStates || withRates)
      {
        SubstitutionProcessSequenceSimulator* pss=dynamic_cast<SubstitutionProcessSequenceSimulator*>(&seqsim);

        if (withStates)
          if (withRates)
            sites = pss?pss->simulate(rates, states):SequenceSimulationTools::simulateSites(seqsim, rates, states);
          else
            sites = pss?pss->simulate(states):SequenceSimulationTools::simulateSites(seqsim, states);
        else
          sites = pss?pss->simulate(rates):SequenceSimulationTools::simulateSites(seqsim, rates);
        
        ApplicationTools::displayTaskDone();
      }
      else
      {
        size_t nSites=(lengths[it.first]==0?nbSites:lengths[it.first]);
      
        ApplicationTools::displayMessage("");
        ApplicationTools::displayTask("Perform simulations");

        sites = seqsim.simulate(nSites);
        ApplicationTools::displayTaskDone();
        ApplicationTools::displayMessage("");
      }

      // Write to file:
      BppOAlignmentWriterFormat bppoWriter(1);
      unique_ptr<OAlignment> oAln(bppoWriter.read(formats[it.first]));
      // ApplicationTools::displayResult("Output alignment file ", filenames[it.first]);
      // ApplicationTools::displayResult("Output alignment format ", oAln->getFormatName());

      oAln->writeAlignment(filenames[it.first], *sites, true);

    }

    for (auto& it : mSim)
      delete it.second;

    delete alphabet;

    bppseqgen.done();
  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}
