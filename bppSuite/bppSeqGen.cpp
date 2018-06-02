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
#include <Bpp/App/BppApplication.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Text/KeyvalTools.h>
#include <Bpp/App/NumCalcApplicationTools.h>

// From bpp-seq:
#include <Bpp/Seq/SequenceTools.h>
#include <Bpp/Seq/Io/BppOAlignmentWriterFormat.h>

// From PhylLib:
#include <Bpp/Phyl/Simulation/SequenceSimulationTools.h>
#include <Bpp/Phyl/Simulation/EvolutionSequenceSimulator.h>

#include "bppTools.h"


using namespace bpp;


map<size_t, SequenceSimulator*> readSimul(const SubstitutionProcessCollection& spc,
                                          const map<size_t, SequenceEvolution*>& mSeqEvol,
                                          map<string, string>& params,
                                          map<size_t, string>& mfnames,
                                          map<size_t, string>& mformats,
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

    while (args.find("process")==args.end())
      throw BadIntegerException("bppseqgen::readSimul. Missing process argument for simul:",(int)num);

    size_t indProcess=(size_t)ApplicationTools::getIntParameter("process", args, 1, "", true, 0);

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
    
    
    // output

    mfnames[num] = ApplicationTools::getAFilePath("output.sequence.file", args, true, false, "", false, "none", true);

    mformats[num] = ApplicationTools::getStringParameter("output.sequence.format", args, "Fasta", "", false, true);

    mlength[num] = (size_t)ApplicationTools::getIntParameter("number_of_sites", args, 100, "", false, 0);

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

  if(args == 1)
  {
    bppTools::help("bppSeqGen");
    return 0;
  }

  try {

    BppApplication bppseqgen(args, argv, "bppseqgen");
    bppseqgen.startTimer();
    map<string, string> unparsedparams;

    Alphabet* alphabet = bppTools::getAlphabet(bppseqgen.getParams());

    unique_ptr<GeneticCode> gCode(bppTools::getGeneticCode(bppseqgen.getParams(), alphabet));

    ////// Get the optional map of the sequences
    
    map<size_t, AlignedValuesContainer*> mSites = bppTools::getAlignmentsMap(bppseqgen.getParams(), alphabet, true, true);


    /// collection
    
    SubstitutionProcessCollection* SPC = bppTools::getCollection(bppseqgen.getParams(), alphabet, gCode.get(), mSites, unparsedparams);

    map<size_t, SequenceEvolution*> mSeqEvol = bppTools::getProcesses(bppseqgen.getParams(), *SPC, unparsedparams);

    /*******************************************/
    /*     Starting sequence                   */
    /*******************************************/

    SiteContainer* sites = 0;
    size_t nbSites = 0;

    bool outputInternalSequences = ApplicationTools::getBooleanParameter("output.internal.sequences", bppseqgen.getParams(), false, "", true, 1);

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
      ApplicationTools::displayResult("Number of sites", TextTools::toString(nbSites));
      string rateCol = ApplicationTools::getStringParameter("input.infos.rates", bppseqgen.getParams(), "pr", "", true, true);
      string stateCol = ApplicationTools::getStringParameter("input.infos.states", bppseqgen.getParams(), "none", "", true, true);
      withRates = rateCol != "none";
      withStates = stateCol != "none";

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
          ApplicationTools::displayResult("Number of sites", TextTools::toString(nbSites));

          delete pseq;
        }
      }
      catch (Exception& e)
      {
        // cout << e.what() << endl;
        // return 1;
      }
    }
  
    // if (nbSites == 0)
    //   nbSites = ApplicationTools::getParameter<size_t>("number_of_sites", bppseqgen.getParams(), 100);

    /*******************/
    /* Simulations     */
    /*******************/

    map<size_t, string> filenames;
    map<size_t, string> formats;
    map<size_t, size_t> lengths;

    map<size_t, SequenceSimulator*> mSim=readSimul(*SPC, mSeqEvol, bppseqgen.getParams(),filenames, formats, lengths);

    for (auto& it : mSim)
    {
      SequenceSimulator& seqsim = *it.second;

      if (withStates || withRates)
      {
        SimpleSubstitutionProcessSequenceSimulator* ps=dynamic_cast<SimpleSubstitutionProcessSequenceSimulator*>(&seqsim);

        if (ps)
        {
          ps->outputInternalSequences(outputInternalSequences);
          if (withStates)
            if (withRates)
              sites = SequenceSimulationTools::simulateSites(*ps, rates, states);
            else
              sites = SequenceSimulationTools::simulateSites(*ps, states);
          else
            sites = SequenceSimulationTools::simulateSites(*ps, rates);
        }
        else
        {
          SubstitutionProcessSequenceSimulator* pss=dynamic_cast<SubstitutionProcessSequenceSimulator*>(&seqsim);

          if (pss)
          {
            pss->outputInternalSequences(outputInternalSequences);
            if (withStates)
              if (withRates)
                sites = pss->simulate(rates, states);
              else
                sites = pss->simulate(states);
            else
              sites = pss->simulate(rates);
          }
        }
        ApplicationTools::displayTaskDone();
      }

      else
      {
        size_t nSites=(lengths[it.first]==0?nbSites:lengths[it.first]);
      
        ApplicationTools::displayResult("Number of sites", TextTools::toString(nSites));
        ApplicationTools::displayTask("Perform simulations");

        sites = seqsim.simulate(nSites);
        ApplicationTools::displayTaskDone();
      }

      // Write to file:
      BppOAlignmentWriterFormat bppoWriter(1);
      unique_ptr<OAlignment> oAln(bppoWriter.read(formats[it.first]));
      ApplicationTools::displayResult("Output alignment file ", filenames[it.first]);
      ApplicationTools::displayResult("Output alignment format ", oAln->getFormatName());

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
