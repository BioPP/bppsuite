// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

// From the STL:
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

#include <Bpp/Version.h>
#include <Bpp/Numeric/Range.h>

// // From bpp-seq:
#include <Bpp/Seq/Io/Mase.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/SequenceTools.h>
#include <Bpp/Seq/App/BppSequenceApplication.h>

using namespace bpp;

int main(int args, char** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*              Bio++ Alignment Score, version " << BPP_VERSION << "              *" << endl;
  cout << "* Author: J. Dutheil                        Last Modif. " << BPP_REL_DATE << " *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  try
  {
    BppSequenceApplication bppalnscore(args, argv, "BppAlnScore");

    if (args == 1)
    {
      bppalnscore.help("bppalnscore");
      return 0;
    }

    bppalnscore.startTimer();

    // Get alphabet
    shared_ptr<Alphabet> alphabet(bppalnscore.getAlphabet());

    // Get the test alignment:
    auto sitesTest = SequenceApplicationTools::getSiteContainer(alphabet, bppalnscore.getParams(), ".test", false, true);

    // Get the reference alignment:
    shared_ptr<SiteContainerInterface> sitesRef;

    sitesRef = SequenceApplicationTools::getSiteContainer(alphabet, bppalnscore.getParams(), ".ref", false, true);

    // We check if the two alignments are compatible:
    vector<string> namesTest = sitesTest->getSequenceNames();
    vector<string> namesRef  = sitesRef->getSequenceNames();
    if (namesTest != namesRef)
    {
      ApplicationTools::displayTask("Reorder sequences in ref. alignment", true);
      unique_ptr<AlignedSequenceContainer> tmp(new AlignedSequenceContainer(sitesRef->getAlphabet()));
      for (size_t i = 0; i < namesTest.size(); ++i)
      {
        ApplicationTools::displayGauge(i, namesTest.size() - 1);
        try
        {
          auto seq = unique_ptr<Sequence>(sitesRef->sequence(namesTest[i]).clone());
          tmp->addSequence(namesRef[i], seq);
        }
        catch (SequenceNotFoundException& ex)
        {
          throw Exception("ERROR!!! Reference alignment should contain the same sequences as the test alignment!");
        }
      }
      ApplicationTools::displayTaskDone();

      sitesRef = shared_ptr<SiteContainerInterface>(tmp.release());
    }

    // Build alignment indexes:
    RowMatrix<size_t> indexTest, indexRef;
    SiteContainerTools::getSequencePositions(*sitesTest, indexTest);
    SiteContainerTools::getSequencePositions(*sitesRef,  indexRef);

    // Now build scores:
    int na = ApplicationTools::getIntParameter("score.na", bppalnscore.getParams(), 0);
    ApplicationTools::displayResult("NA value to used", na);
    vector<int> cs = SiteContainerTools::getColumnScores(indexTest, indexRef, na);
    vector<double> sps = SiteContainerTools::getSumOfPairsScores(indexTest, indexRef, static_cast<double>(na));

    // Should scores be averaged for words?
    size_t wsize = ApplicationTools::getParameter<size_t>("score.word_size", bppalnscore.getParams(), 1);
    size_t phase = 0;
    if (wsize > 1)
    {
      ApplicationTools::displayResult("Scores uniformized for words of size", wsize);
      string phaseOpt = ApplicationTools::getStringParameter("score.phase", bppalnscore.getParams(), "1");
      if (TextTools::isDecimalInteger(phaseOpt))
      {
        phase = TextTools::to<size_t>(phaseOpt);
        if (phase == 0)
          throw Exception("ERROR: positions are 1-based.");
        phase--;
      }
      else
      {
        // We look for the first occurrence of the given motif:
        try
        {
          auto alpha = sitesTest->getAlphabet();
          Sequence motif("motif", phaseOpt, alpha);
          ApplicationTools::displayResult("Phase based on 1st occurrence of", motif.toString());
          size_t pos = sitesTest->getNumberOfSites();
          for (size_t i = 0; i < sitesTest->getNumberOfSequences(); ++i)
          {
            size_t p = SequenceTools::findFirstOf(sitesTest->sequence(i), motif);
            if (p < pos)
              pos = p;
          }
          phase = pos;
        }
        catch (Exception& ex)
        {
          throw Exception("Error, invalid motif specified for phase option.");
        }
      }
      ApplicationTools::displayResult("First word starts at", phase + 1);

      // Now perform the smoothing:
      size_t i;
      for (i = 0; i < phase; ++i)
      {
        cs[i] = 0;
        sps[i] = 0;
      }
      for ( ; i + wsize <= cs.size(); i += wsize)
      {
        // First compute minimum criterion:
        int csmin = 1;
        double spsmin = 1;
        for (size_t j = i; j < i + wsize; ++j)
        {
          if (cs[j] < csmin)
            csmin = cs[j];
          if (sps[j] < spsmin)
            spsmin = sps[j];
        }
        // Assign min to all positions in word:
        for (size_t j = i; j < i + wsize; ++j)
        {
          cs[j] = csmin;
          sps[j] = spsmin;
        }
      }
      for ( ; i < cs.size(); ++i)
      {
        cs[i] = 0;
        sps[i] = 0;
      }
    }

    // Output scores to SGED file:
    string outputScores = ApplicationTools::getAFilePath("output.scores", bppalnscore.getParams(), false, false);
    if (outputScores != "none")
    {
      ApplicationTools::displayResult("Output scores to", outputScores);
      ofstream output(outputScores.c_str(), ios::out);
      output << "Site\tColumnScore\tSumOfPairsScore" << endl;
      for (size_t i = 0; i < cs.size(); ++i)
      {
        output << "[" << sitesTest->site(i).getCoordinate() << "]\t" << cs[i] << "\t" << sps[i] << endl;
      }
      output.close();
    }

    // Create a sequence filter:
    string outputFilter = ApplicationTools::getAFilePath("output.mase", bppalnscore.getParams(), false, false);
    if (outputFilter != "none")
    {
      ApplicationTools::displayResult("Output mase with site filter to", outputFilter);
      double spsThreshold = ApplicationTools::getDoubleParameter("output.sps_thresholds", bppalnscore.getParams(), 0.8);
      ApplicationTools::displayResult("Threshold for SPS", spsThreshold);

      MultiRange<size_t> csRanges;
      MultiRange<size_t> spsRanges;
      size_t csBeg = 0, spsBeg = 0, csEnd = 0, spsEnd = 0;
      size_t s = alphabet->getStateCodingSize();
      for (size_t i = 0; i < cs.size(); ++i)
      {
        if (cs[i] == 1 && i > 0 && cs[i - 1] != 1)
          csBeg = i;
        if (cs[i] != 1 && i > 0 && cs[i - 1] == 1)
        {
          csEnd = i;
          csRanges.addRange(Range<size_t>(csBeg * s, csEnd * s));
        }

        if (sps[i] >= spsThreshold && i > 0 && sps[i - 1] < spsThreshold)
          spsBeg = i;
        if (sps[i] < spsThreshold && i > 0 && sps[i - 1] >= spsThreshold)
        {
          spsEnd = i;
          spsRanges.addRange(Range<size_t>(spsBeg * s, spsEnd * s));
        }
      }
      // Add the last range if any:
      if (cs.back() == 1)
        csRanges.addRange(Range<size_t>(csBeg * s, cs.size() * s));
      if (sps.back() >= spsThreshold)
        spsRanges.addRange(Range<size_t>(spsBeg * s, sps.size() * s));

      MaseHeader header;
      header.setSiteSelection("CS", csRanges);
      header.setSiteSelection("SPS", spsRanges);
      Mase writer;
      writer.writeMeta(outputFilter, *sitesTest, header);
    }

    // Create an SGED index file(s):
    string outputIndexCs = ApplicationTools::getAFilePath("output.index.cs", bppalnscore.getParams(), false, false);
    ofstream indexOutCs(outputIndexCs, ios::out);
    if (outputIndexCs != "none")
    {
      ApplicationTools::displayResult("Output CS index to", outputIndexCs);

      indexOutCs << "# SGED index file version 1.00" << endl;
      indexOutCs << "# SGED index start" << endl;
      indexOutCs << "AlnPos,OrigPos" << endl;

      size_t pos = 0;
      for (size_t i = 0; i < cs.size(); ++i)
      {
        if (cs[i] == 1)
        {
          indexOutCs << ++pos << "," << sitesTest->site(i).getCoordinate() << endl;
        }
      }
    }

    string outputIndexSps = ApplicationTools::getAFilePath("output.index.sps", bppalnscore.getParams(), false, false);
    ofstream indexOutSps(outputIndexSps, ios::out);
    if (outputIndexSps != "none")
    {
      ApplicationTools::displayResult("Output SPS index to", outputIndexSps);
      double spsThreshold = ApplicationTools::getDoubleParameter("output.sps_thresholds", bppalnscore.getParams(), 0.8);

      indexOutSps << "# SGED index file version 1.00" << endl;
      indexOutSps << "# SGED index start" << endl;
      indexOutSps << "AlnPos,OrigPos" << endl;

      size_t pos = 0;
      for (size_t i = 0; i < sps.size(); ++i)
      {
        if (sps[i] >= spsThreshold)
        {
          indexOutSps << ++pos << "," << sitesTest->site(i).getCoordinate() << endl;
        }
      }
    }

    // We're done!
    bppalnscore.done();
  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}
