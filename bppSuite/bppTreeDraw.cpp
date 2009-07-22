//
// File: bppTreeDraw.cpp
// Created by: Julien Dutheil
// Created on: Jul Tue 21 13:40 2009
//

/*
Copyright or © or Copr. CNRS

This software is a computer program whose purpose is to draw phylogenies.

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
#include <iomanip>

using namespace std;

// From Utils:
#include <Utils/AttributesTools.h>
#include <Utils/ApplicationTools.h>
#include <Utils/KeyvalTools.h>
#include <Utils/graphics>

// From PhylLib:
#include <Phyl/Tree.h>
#include <Phyl/PhylogeneticsApplicationTools.h>
#include <Phyl/CladogramPlot.h>
#include <Phyl/PhylogramPlot.h>
using namespace bpp;

/******************************************************************************/

void help()
{
  *ApplicationTools::message << "__________________________________________________________________________" << endl;
  *ApplicationTools::message << "bpptreedraw parameter1_name=parameter1_value parameter2_name=parameter2_value"    << endl;
  *ApplicationTools::message << "      ... param=option_file" << endl;
  *ApplicationTools::message << endl;
  *ApplicationTools::message << "  Refer to the Bio++ Program Suite Manual for a list of available options." << endl;
  *ApplicationTools::message << "__________________________________________________________________________" << endl;
}

int main(int args, char ** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*       Bio++ Tree Drawing program, version 0.1.0                *" << endl;
  cout << "*                                                                *" << endl; 
  cout << "* Authors: J. Dutheil                       Last Modif. 21/07/09 *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  if (args == 1)
  {
    help();
    return 0;
  }
  
  try {

  ApplicationTools::startTimer();

  cout << "Parsing options:" << endl;
  
  map<string, string> params = AttributesTools::parseOptions(args, argv);

  // Get the tree to plot:
  Tree* tree = PhylogeneticsApplicationTools::getTree(params);
  ApplicationTools::displayResult("Number of leaves", TextTools::toString(tree->getNumberOfLeaves()));
  
  // Get the graphic device:
  GraphicDevice* gd = 0;
	string outputPath = ApplicationTools::getAFilePath("output.drawing.file", params, true, false, "", false);
  ofstream file(outputPath.c_str(), ios::out);
  string graphicType = ApplicationTools::getStringParameter("output.drawing.format", params, "svg");
  if (graphicType == "svg")
  {
    gd = new SVGGraphicDevice(file);
  }
  else if (graphicType == "inkscape")
  {
    gd = new SVGGraphicDevice(file, true);
  }
  else if (graphicType == "xfig")
  {
    gd = new XFigGraphicDevice(file);
    dynamic_cast<XFigGraphicDevice *>(gd)->setFontFlag(XFigGraphicDevice::FONTFLAG_POSTSCRIPT);
  }
  else if (graphicType == "pgf")
  {
    gd = new PGFGraphicDevice(file, 0.045);
  }
  else throw Exception("Unknown output format: " + graphicType);

  // Get the tree plotter:
  TreeDrawing* td = 0;
  string plotType = ApplicationTools::getStringParameter("output.drawing.plot", params, "cladogram");
  string name;
  map<string, string> args;
  KeyvalTools::parseProcedure(plotType, name, args);
  if(name == "cladogram")
  {
    td = new CladogramPlot(tree);
  }
  else if(name == "phylogram")
  {
    td = new PhylogramPlot(tree);
  }
  else throw Exception("Unknown output format: " + plotType);
  ApplicationTools::displayResult("Plot type", name);
  double xunit = ApplicationTools::getDoubleParameter("xu", args, 10);
  double yunit = ApplicationTools::getDoubleParameter("yu", args, 10);
  td->setXUnit(xunit);
  td->setYUnit(yunit);

  //Now draw the tree:
  gd->begin();
  td->plot(*gd);
  gd->end();

  //Finishing things:
  file.close();
  delete tree;
  delete td;
  delete gd;

  cout << "BppTreeDraw's done. Bye." << endl;
  ApplicationTools::displayTime("Total execution time:");
 
  }
  catch(exception & e)
  {
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}

