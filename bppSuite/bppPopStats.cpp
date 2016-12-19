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

using namespace std;

#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/KeyvalTools.h>

// From bpp-seq:
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/CodonSiteTools.h>
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

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
  cout << "*              Bio++ Population Statistics, version 2.3.0        *" << endl;
  cout << "* Author: J. Dutheil                        Last Modif. 24/06/15 *" << endl;
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
  unique_ptr<ofstream> clog;
  if (logFile != "none")
    clog.reset(new ofstream(logFile.c_str(), ios::out));

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
      gCode.reset(SequenceApplicationTools::getGeneticCode(codonAlphabet->getNucleicAlphabet(), codeDesc));
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
          *clog << "# Info: last site contained a stop codon and was discarded." << endl;
      }
    } else if (stopCodonOpt == "RemoveAll") {
      size_t l1 = psc->getNumberOfSites();
      SiteContainerTools::removeStopCodonSites(*psc, *gCode);
      size_t l2 = psc->getNumberOfSites();
      if (l2 != l1) {
        ApplicationTools::displayMessage("Info: discarded " + TextTools::toString(l2 - l1) + " sites with stop codons.");
        if (logFile != "none")
          *clog << "# Info: discarded " << (l2 - l1) << " sites with stop codons." << endl;
      }
    } else {
      throw Exception("Unrecognized option for input.sequence.stop_codons_policy: " + stopCodonOpt);
    }

    unique_ptr<PolymorphismSequenceContainer> pscIn;
    unique_ptr<PolymorphismSequenceContainer> pscOut;

    if (psc->hasOutgroup()) {
      pscIn.reset(PolymorphismSequenceContainerTools::extractIngroup(*psc));
      pscOut.reset(PolymorphismSequenceContainerTools::extractOutgroup(*psc));
    } else {
      pscIn = std::move(psc);
    }
    ApplicationTools::displayResult("Number of sequences in ingroup", pscIn->getNumberOfSequences());
    ApplicationTools::displayResult("Number of sequences in outgroup", pscOut.get() ? pscOut->getNumberOfSequences() : 0);
    
    // +-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+-+
    
    // Compute statistics
    vector<string> actions = ApplicationTools::getVectorParameter<string>("pop.stats", bpppopstats.getParams(), ',', "", "", false, 1);

    for (size_t a = 0; a < actions.size(); a++)
    {
      string cmdName;
      map<string, string> cmdArgs;
      KeyvalTools::parseProcedure(actions[a], cmdName, cmdArgs);

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
          *clog << "# Site frequencies" << endl;
          *clog << "NbSegSites = " << s << endl;
          *clog << "NbSingl = " << nsg << endl;
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
          *clog << "# Watterson's (1975) theta" << endl;
          *clog << "thetaW75 = " << thetaW75 << endl;
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
          *clog << "# Tajima's (1983) pi" << endl;
          *clog << "piT83 = " << piT83 << endl;
        }
      }

      // +------------+
      // | Tajima's D |
      // +------------+
      else if (cmdName == "TajimaD")
      {
        double tajimaD = SequenceStatistics::tajimaDss(*pscIn, true, true);
        ApplicationTools::displayResult("Tajima's (1989) D:", tajimaD);
        //Print to logfile:
        if (logFile != "none") {
          *clog << "# Tajima's (1989) D" << endl;
          *clog << "tajD = " << tajimaD << endl;
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
          *clog << "# Fu and Li's (1993) D*" << endl;
          if (useTotMut)
            *clog << "fuLiDstarTotMut = " << flDstar << endl;
          else
            *clog << "fuLiDstarSegSit = " << flDstar << endl;
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
          *clog << "# Fu and Li's (1993) F*" << endl;
          if (useTotMut)
            *clog << "fuLiFstarTotMut = " << flFstar << endl;
          else
            *clog << "fuLiFstarSegSit = " << flFstar << endl;
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
        double nbS = SequenceStatistics::meanNumberOfSynonymousSites(*pscIn, *gCode);
        double nbN = SequenceStatistics::meanNumberOfNonSynonymousSites(*pscIn, *gCode);
        double r = (piN / nbN) / (piS / nbS);
        ApplicationTools::displayResult("PiN:", piN);
        ApplicationTools::displayResult("PiS:", piS);
        ApplicationTools::displayResult("#N:", nbN);
        ApplicationTools::displayResult("#S:", nbS);
        ApplicationTools::displayResult("PiN / PiS (corrected for #N and #S):", r);
        if (logFile != "none") {
          *clog << "# PiN and PiS" << endl;
          *clog << "PiN = " << piN << endl;
          *clog << "PiS = " << piS << endl;
          *clog << "NbN = " << nbN << endl;
          *clog << "NbS = " << nbS << endl;
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
        if (!pscOut.get()) {
          throw Exception("MacDonald-Kreitman test requires at least one outgroup sequence.");
        }
        vector<unsigned int> mktable = SequenceStatistics::mkTable(*pscIn, *pscOut, *gCode);
        ApplicationTools::displayResult("MK table, Pa:", mktable[0]);
        ApplicationTools::displayResult("MK table, Ps:", mktable[1]);
        ApplicationTools::displayResult("MK table, Da:", mktable[2]);
        ApplicationTools::displayResult("MK table, Ds:", mktable[3]);
        if (logFile != "none") {
          *clog << "# MK table" << endl;
          *clog << "# Pa Ps Da Ds" << endl;
          *clog << "MKtable = " << mktable[0] << " " << mktable[1] << " " << mktable[2] << " " << mktable[3] << endl;
        }
      }

      // +---------------------+
      // | Codon site statistics |
      // +---------------------+
      else if (cmdName == "CodonSiteStatistics")
      {
        if (!codonAlphabet) {
          throw Exception("CodonSiteStatstics can only be used with a codon alignment. Check the input alphabet!");
        }
        string path = ApplicationTools::getAFilePath("output.file", cmdArgs, true, false);
        if (path == "none") throw Exception("You must specify an ouptut file for CodonSiteStatistics"); 
        ApplicationTools::displayResult("Site statistics output to:", path);
        ofstream out(path.c_str(), ios::out);
        out << "Site\tIsComplete\tNbAlleles\tMinorAlleleFrequency\tMajorAlleleFrequency\tMinorAllele\tMajorAllele\tMeanNumberSynPos\tIsSynPoly\tIs4Degenerated\tPiN\tPiS" << endl;
        unique_ptr<SiteContainer> sites(pscIn->toSiteContainer());
        for (size_t i = 0; i < sites->getNumberOfSites(); ++i) {
          const Site& site = sites->getSite(i);
          out << site.getPosition() << "\t";
          out << SiteTools::isComplete(site) << "\t";
          out << SiteTools::getNumberOfDistinctCharacters(site) << "\t";
          out << SiteTools::getMinorAlleleFrequency(site) << "\t";
          out << SiteTools::getMajorAlleleFrequency(site) << "\t";
          out << alphabet->intToChar(SiteTools::getMinorAllele(site)) << "\t";
          out << alphabet->intToChar(SiteTools::getMajorAllele(site)) << "\t";
          out << CodonSiteTools::meanNumberOfSynonymousPositions(site, *gCode) << "\t";
          out << CodonSiteTools::isSynonymousPolymorphic(site, *gCode) << "\t";
          out << CodonSiteTools::isFourFoldDegenerated(site, *gCode) << "\t";
          out << CodonSiteTools::piNonSynonymous(site, *gCode) << "\t";
          out << CodonSiteTools::piSynonymous(site, *gCode) << endl;
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
      *clog << "# Error: " << e.what() << endl;
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}

