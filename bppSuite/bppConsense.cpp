// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

// From the STL:
#include <iostream>
#include <iomanip>

using namespace std;

// From bpp-core:
#include <Bpp/Version.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/KeyvalTools.h>

// From PhylLib:
#include <Bpp/Phyl/Tree/Tree.h>
#include <Bpp/Phyl/Tree/TreeTemplate.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>

using namespace bpp;

void help()
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
  (*ApplicationTools::message << "bppconsense parameter1_name=parameter1_value").endLine();
  (*ApplicationTools::message << "      parameter2_name=parameter2_value ... param=option_file").endLine();
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "  Refer to the Bio++ Program Suite Manual for a list of available options.").endLine();
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
}

int main(int args, char** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*       Bio++ Consensus and Bootstrap Methods, version " << BPP_VERSION << "     *" << endl;
  cout << "* Authors: J. Dutheil                       Created     06/06/07 *" << endl;
  cout << "*          N. Galtier                       Last Modif. " << BPP_REL_DATE << " *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  if (args == 1)
  {
    help();
    return 0;
  }

  try
  {
    BppApplication bppconsense(args, argv, "BppConsense");
    bppconsense.startTimer();

    auto list = PhylogeneticsApplicationTools::getTrees(bppconsense.getParams());

    unique_ptr<Tree> tree = nullptr;
    string treeMethod = ApplicationTools::getStringParameter("tree", bppconsense.getParams(), "Consensus", "", false, 1);
    string cmdName;
    map<string, string> cmdArgs;
    KeyvalTools::parseProcedure(treeMethod, cmdName, cmdArgs);
    if (cmdName == "Input")
    {
      tree = PhylogeneticsApplicationTools::getTree(bppconsense.getParams());
      ApplicationTools::displayResult("Number of leaves", tree->getNumberOfLeaves());
    }
    else if (cmdName == "Consensus")
    {
      double threshold = ApplicationTools::getDoubleParameter("threshold", cmdArgs, 0, "", false, 1);
      ApplicationTools::displayResult("Consensus threshold", TextTools::toString(threshold));
      ApplicationTools::displayTask("Computing consensus tree");
      tree = TreeTools::thresholdConsensus(list, threshold, true);
      ApplicationTools::displayTaskDone();
    }
    else
      throw Exception("Unknown input tree method: " + treeMethod);

    ApplicationTools::displayTask("Compute bootstrap values");

    int bsformat = ApplicationTools::getIntParameter("bootstrap.format", bppconsense.getParams(), 0, "", false, 1);
    TreeTools::computeBootstrapValues(*tree, list, true, bsformat);
    ApplicationTools::displayTaskDone();

    // Write resulting tree:
    PhylogeneticsApplicationTools::writeTree(*tree, bppconsense.getParams());

    bppconsense.done();
  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}
