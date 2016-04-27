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
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

// From bpp-popgen
#include <Bpp/PopGen/PolymorphismSequenceContainer.h>
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

  try
  {
    BppApplication bpppopstats(args, argv, "BppPopStats");
    bpppopstats.startTimer();

    // Get alphabet
    Alphabet* alphabet = SequenceApplicationTools::getAlphabet(bpppopstats.getParams(), "", false, true, true);

    // Get the ingroup alignment:
    unique_ptr<SiteContainer> sites(SequenceApplicationTools::getSiteContainer(alphabet, bpppopstats.getParams(), ".ingroup", false, true));
    unique_ptr<PolymorphismSequenceContainer> psc(new PolymorphismSequenceContainer(*sites));

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
        unsigned int s = SequenceStatistics::numberOfPolymorphicSites(*psc);
        ApplicationTools::displayResult("Number of segregating sites:", s);
        unsigned int nsg = SequenceStatistics::numberOfSingletons(*psc);
        ApplicationTools::displayResult("Number of singletons:", nsg);
      }

      // +-------------------+
      // | Watterson's theta |
      // +-------------------+
      if (cmdName == "Watterson75")
      {
        double thetaW75 = SequenceStatistics::watterson75(*psc, true, true, true);
        ApplicationTools::displayResult("Watterson (1975)'s theta:", thetaW75);
      }

      // +----------------+
      // | Tajima's theta |
      // +----------------+
      if (cmdName == "Tajima83")
      {
        double thetaT83 = SequenceStatistics::tajima83(*psc, true, true, true);
        ApplicationTools::displayResult("Tajima (1983)'s theta:", thetaT83);
      }

      // +------------+
      // | Tajima's D |
      // +------------+
      if (cmdName == "TajimaD")
      {
        double tajimaD = SequenceStatistics::tajimaDss(*psc, true, true);
        ApplicationTools::displayResult("Tajima (1989)'s D:", tajimaD);
      }

      // +-----------+
      // | FuAndLiD* |
      // +-----------+
      if (cmdName == "FuAndLiDStar")
      {
        bool useTotMut = ApplicationTools::getBooleanParameter("tot_mut", cmdArgs, true, "", false, 1);
        double flDstar = SequenceStatistics::fuLiDStar(*psc, !useTotMut);
        ApplicationTools::displayResult("Fu and Li (1993)'s D*:", flDstar);
        ApplicationTools::displayResult("  computed using", (useTotMut ? "total number of mutations" : "number of segregating sites"));
      }

      // +-----------+
      // | FuAndLiF* |
      // +-----------+
      if (cmdName == "FuAndLiFStar")
      {
        bool useTotMut = ApplicationTools::getBooleanParameter("tot_mut", cmdArgs, true, "", false, 1);
        double flFstar = SequenceStatistics::fuLiFStar(*psc, !useTotMut);
        ApplicationTools::displayResult("Fu and Li (1993)'s F*:", flFstar);
        ApplicationTools::displayResult("  computed using", (useTotMut ? "total number of mutations" : "number of segregating sites"));
      }


    }
 
    // We're done!
    bpppopstats.done();
  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}

