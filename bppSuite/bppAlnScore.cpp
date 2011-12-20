//
// File: bppAlnScore.cpp
// Created by: Julien Dutheil
// Created on: Dec Thu 15 16:16 2011
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

using namespace std;

#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Text/TextTools.h>

// From SeqLib:
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/Io.all>
#include <Bpp/Seq/Container.all>
#include <Bpp/Seq/SequenceTools.h>

using namespace bpp;

void help()
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
  (*ApplicationTools::message << "bppalnscore parameter1_name=parameter1_value").endLine();
  (*ApplicationTools::message << "      parameter2_name=parameter2_value ... param=option_file").endLine();
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "  Refer to the Bio++ Program Suite Manual for a list of available options.").endLine();
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
}

int main(int args, char** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*              Bio++ Alignment Score, version 0.1                *" << endl;
  cout << "* Author: J. Dutheil                        Last Modif. 15/12/11 *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;
  
  if (args == 1)
  {
    help();
    return 0;
  }
  
  try {

  BppApplication bppalnscore(args, argv, "BppAlnScore");
  bppalnscore.startTimer();
  
  // Get alphabet
  Alphabet* alphabet = SequenceApplicationTools::getAlphabet(bppalnscore.getParams(), "", false, true, true);

  // Get the test alignment:
  auto_ptr<SiteContainer> sitesTest(SequenceApplicationTools::getSiteContainer(alphabet, bppalnscore.getParams(), ".test", false, true));

  // Get the reference alignment:
  auto_ptr<SiteContainer> sitesRef(SequenceApplicationTools::getSiteContainer(alphabet, bppalnscore.getParams(), ".ref", false, true));
  
  //We check if the two alignments are compatible:
  vector<string> namesTest = sitesTest->getSequencesNames();
  vector<string> namesRef  = sitesRef->getSequencesNames();
  if (namesTest != namesRef) {
    ApplicationTools::displayTask("Reorder sequences in ref. alignment", true);
    auto_ptr<AlignedSequenceContainer> tmp(new AlignedSequenceContainer(sitesRef->getAlphabet()));
    for (size_t i = 0; i < namesTest.size(); ++i) {
      ApplicationTools::displayGauge(i, namesTest.size() - 1);
      try {
        tmp->addSequence(sitesRef->getSequence(namesTest[i]));
      } catch(SequenceNotFoundException& ex) {
        throw Exception("ERROR!!! Reference alignment should contain the same sequences as the test alignment!");
      }
    }
    ApplicationTools::displayTaskDone();
    sitesRef = tmp;
  }

  //Build alignment indexes:
  RowMatrix<unsigned int> indexTest, indexRef;
  SiteContainerTools::getSequencePositions(*sitesTest, indexTest);
  SiteContainerTools::getSequencePositions(*sitesRef,  indexRef);

  //Now build scores:
  int na = ApplicationTools::getIntParameter("score.na", bppalnscore.getParams(), 0);
  ApplicationTools::displayResult("NA value to used", na);
  vector<int> cs = SiteContainerTools::getColumnScores(indexTest, indexRef, na);
  vector<double> sps = SiteContainerTools::getSumOfPairsScores(indexTest, indexRef, static_cast<double>(na));

  //Should scores be averaged for words?
  size_t wsize = ApplicationTools::getParameter<size_t>("score.word_size", bppalnscore.getParams(), 1);
  size_t phase = 0;
  if (wsize > 1) {
    ApplicationTools::displayResult("Scores uniformized for words of size", wsize);
    string phaseOpt = ApplicationTools::getStringParameter("score.phase", bppalnscore.getParams(), "1");
    if (TextTools::isDecimalInteger(phaseOpt)) {
      phase = TextTools::toInt(phaseOpt);
      if (phase == 0)
        throw Exception("ERROR: positions are 1-based.");
      phase--;
    } else {
      //We look for the first occurrence of the given motif:
      try {
        BasicSequence motif("motif", phaseOpt, sitesTest->getAlphabet());
        ApplicationTools::displayResult("Phase based on 1st occurence of", motif.toString());
        unsigned int pos = sitesTest->getNumberOfSites();
        for (unsigned int i = 0; i < sitesTest->getNumberOfSequences(); ++i) {
          unsigned int p = SequenceTools::findFirstOf(sitesTest->getSequence(i), motif);
          if (p < pos)
            pos = p;
        }
        phase = pos;
      } catch(Exception& ex) {
        throw Exception("Error, unvalid motif specified for phase option.");
      }
    }
    ApplicationTools::displayResult("First word starts at", phase + 1);

    //Now perform the smoothing:
    size_t i;
    for (i = 0; i < phase; ++i) {
      cs[i] = 0;
      sps[i] = 0;
    }
    for (; i + wsize < cs.size(); i += wsize) {
      //First compute minimum criterion:
      int csmin = 1;
      double spsmin = 1;
      for (size_t j = i; j < i + wsize; ++j) {
        if (cs[j] < csmin) csmin = cs[j];
        if (sps[j] < spsmin) spsmin = sps[j];
      }
      //Assign min to all positions in word:
      for (size_t j = i; j < i + wsize; ++j) {
        cs[j] = csmin;
        sps[j] = spsmin;
      }
    }
    for (; i < cs.size(); ++i) {
      cs[i] = 0;
      sps[i] = 0;
    }
  }

  //Output scores to file:
  string outputScores = ApplicationTools::getAFilePath("output.scores", bppalnscore.getParams(), true, false);
  if (outputScores != "none") {
    ApplicationTools::displayResult("Output scores to", outputScores);
    ofstream output(outputScores.c_str(), ios::out);
    output << "Site\tColumnScore\tSumOfPairsScore" << endl;
    for (size_t i = 0; i < cs.size(); ++i) {
      output << sitesTest->getSite(i).getPosition() << "\t" << cs[i] << "\t" << sps[i] << endl;
    }
    output.close();
  }

  //Create a sequence filter:
  string outputFilter = ApplicationTools::getAFilePath("output.mase", bppalnscore.getParams(), true, false);
  if (outputFilter != "none") {
    ApplicationTools::displayResult("Output mase with site filter to", outputFilter);
    double spsThreshold = ApplicationTools::getDoubleParameter("output.sps_thresholds", bppalnscore.getParams(), 0.8);
    ApplicationTools::displayResult("Threshold for SPS", spsThreshold);
    
    MultiRange<unsigned int> csRanges;
    MultiRange<unsigned int> spsRanges;
    for (size_t i = 0; i < cs.size(); ++i) {
      if (cs[i] == 1) csRanges.addRange(Range<unsigned int>(i, i + 1));
      if (sps[i] >= spsThreshold) spsRanges.addRange(Range<unsigned int>(i, i + 1));
    }
    MaseHeader header;
    header.setSiteSelection("CS", csRanges);
    header.setSiteSelection("SPS", spsRanges);
    Mase writer;
    writer.write(outputFilter, *sitesTest, header);
  }

  //We're done!
  bppalnscore.done();

  } catch(exception & e) {
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}

