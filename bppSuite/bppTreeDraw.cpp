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
#include <Bpp/Text/KeyvalTools.h>
#include <Bpp/Graphics/Svg/SvgGraphicDevice.h>
#include <Bpp/Graphics/Latex/PgfGraphicDevice.h>
#include <Bpp/Graphics/Fig/XFigGraphicDevice.h>

// From bpp-phyl:
#include <Bpp/Phyl/Tree/Tree.h>

#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Graphics/PhylogramPlot.h>
#include <Bpp/Phyl/Graphics/CladogramPlot.h>
#include <Bpp/Phyl/Graphics/TreeDrawingDisplayControler.h>

using namespace bpp;

/******************************************************************************/

void help()
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
  (*ApplicationTools::message << "bpptreedraw parameter1_name=parameter1_value parameter2_name=parameter2_value").endLine();
  (*ApplicationTools::message << "      ... param=option_file").endLine();
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "  Refer to the Bio++ Program Suite Manual for a list of available options.").endLine();
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
}

int main(int args, char** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*       Bio++ Tree Drawing program, version " << BPP_VERSION << ".               *" << endl;
  cout << "*                                                                *" << endl;
  cout << "* Authors: J. Dutheil                       Last Modif. " << BPP_REL_DATE << " *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  if (args == 1)
  {
    help();
    return 0;
  }

  try
  {
    BppApplication bpptreedraw(args, argv, "BppTreeDraw");
    bpptreedraw.startTimer();

    // Get the tree to plot:
    auto tree = PhylogeneticsApplicationTools::getTree(bpptreedraw.getParams());
    ApplicationTools::displayResult("Number of leaves", TextTools::toString(tree->getNumberOfLeaves()));

    // Get the graphic device:
    unique_ptr<GraphicDevice> gd = 0;
    string outputPath = ApplicationTools::getAFilePath("output.drawing.file", bpptreedraw.getParams(), true, false, "", false);
    ofstream file(outputPath.c_str(), ios::out);
    string graphicTypeCmd = ApplicationTools::getStringParameter("output.drawing.format", bpptreedraw.getParams(), "Svg");
    string graphicType;
    map<string, string> graphicTypeArgs;
    KeyvalTools::parseProcedure(graphicTypeCmd, graphicType, graphicTypeArgs);
    if (graphicType == "Svg")
    {
      gd = make_unique<SvgGraphicDevice>(file);
    }
    else if (graphicType == "Inkscape")
    {
      gd = make_unique<SvgGraphicDevice>(file, true);
    }
    else if (graphicType == "Xfig")
    {
      gd = make_unique<XFigGraphicDevice>(file);
      dynamic_cast<XFigGraphicDevice*>(gd.get())->setFontFlag(XFigGraphicDevice::FONTFLAG_POSTSCRIPT);
    }
    else if (graphicType == "Pgf")
    {
      gd = make_unique<PgfGraphicDevice>(file, 0.045);
    }
    else
      throw Exception("Unknown output format: " + graphicType);

    // Get the tree plotter:
    unique_ptr<TreeDrawing> td = 0;
    string plotTypeCmd = ApplicationTools::getStringParameter("output.drawing.plot", bpptreedraw.getParams(), "Cladogram");
    string plotType;
    map<string, string> plotTypeArgs;
    KeyvalTools::parseProcedure(plotTypeCmd, plotType, plotTypeArgs);
    if (plotType == "Cladogram")
    {
      td = make_unique<CladogramPlot>();
    }
    else if (plotType == "Phylogram")
    {
      td = make_unique<PhylogramPlot>();
    }
    else
      throw Exception("Unknown output format: " + plotType);
    td->setTree(tree.release());
    ApplicationTools::displayResult("Plot type", plotType);
    double xunit = ApplicationTools::getDoubleParameter("xu", plotTypeArgs, 10);
    double yunit = ApplicationTools::getDoubleParameter("yu", plotTypeArgs, 10);
    td->setXUnit(xunit);
    td->setYUnit(yunit);
    string hOrientation = ApplicationTools::getStringParameter("direction.h", plotTypeArgs, "left2right");
    if (hOrientation == "left2right")
    {
      dynamic_cast<AbstractDendrogramPlot*>(td.get())->setHorizontalOrientation(AbstractDendrogramPlot::ORIENTATION_LEFT_TO_RIGHT);
    }
    else if (hOrientation == "right2left")
    {
      dynamic_cast<AbstractDendrogramPlot*>(td.get())->setHorizontalOrientation(AbstractDendrogramPlot::ORIENTATION_RIGHT_TO_LEFT);
    }
    else
      throw Exception("Unknown orientation option: " + hOrientation);
    string vOrientation = ApplicationTools::getStringParameter("direction.v", plotTypeArgs, "top2bottom");
    if (vOrientation == "top2bottom")
    {
      dynamic_cast<AbstractDendrogramPlot*>(td.get())->setVerticalOrientation(AbstractDendrogramPlot::ORIENTATION_TOP_TO_BOTTOM);
    }
    else if (vOrientation == "bottom2top")
    {
      dynamic_cast<AbstractDendrogramPlot*>(td.get())->setVerticalOrientation(AbstractDendrogramPlot::ORIENTATION_BOTTOM_TO_TOP);
    }
    else
      throw Exception("Unknown orientation option: " + vOrientation);

    // Plotting option:
    TreeDrawingSettings tds;
    auto controler = make_unique<BasicTreeDrawingDisplayControler>(&tds);
    controler->registerTreeDrawing(td.get());

    bool drawLeafNames       = ApplicationTools::getBooleanParameter("draw.leaves", plotTypeArgs, true);
    bool drawNodesId         = ApplicationTools::getBooleanParameter("draw.ids", plotTypeArgs, false);
    bool drawBranchLengths   = ApplicationTools::getBooleanParameter("draw.brlen", plotTypeArgs, false);
    bool drawBootstrapValues = ApplicationTools::getBooleanParameter("draw.bs", plotTypeArgs, false);

    controler->enableListener(controler->PROPERTY_LEAF_NAMES, drawLeafNames);
    controler->enableListener(controler->PROPERTY_NODE_IDS, drawNodesId);
    controler->enableListener(controler->PROPERTY_BRANCH_LENGTHS, drawBranchLengths);
    controler->enableListener(controler->PROPERTY_BOOTSTRAP_VALUES, drawBootstrapValues);

    ApplicationTools::displayBooleanResult("Draw leaf names", drawLeafNames);
    ApplicationTools::displayBooleanResult("Draw node ids", drawNodesId);
    ApplicationTools::displayBooleanResult("Draw branch lengths", drawBranchLengths);
    ApplicationTools::displayBooleanResult("Draw bootstrap values", drawBootstrapValues);

    // Now draw the tree:
    gd->begin();
    td->plot(*gd);
    gd->end();

    // Finishing things:
    file.close();

    bpptreedraw.done();
  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}
