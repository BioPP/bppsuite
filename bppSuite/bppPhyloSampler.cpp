//
// File: PhyloSample.cpp
// Created by: Julien Dutheil
// Created on: Sunday, December 2nd 2007 16:48
//

/*
Copyright or © or Copr. CNRS

This software is a computer program whose purpose is to estimate
phylogenies and evolutionary parameters from a dataset according to
the maximum likelihood principle.

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

// From NumCalc:
#include <NumCalc/DataTable.h>

// From SeqLib:
#include <Seq/alphabets>
#include <Seq/containers>
#include <Seq/ioseq>
#include <Seq/SiteTools.h>
#include <Seq/SequenceTools.h>
#include <Seq/SequenceApplicationTools.h>

// From Utils:
#include <Utils/AttributesTools.h>
#include <Utils/ApplicationTools.h>
#include <Utils/FileTools.h>
#include <Utils/TextTools.h>

// From PhylLib:
#include <Phyl/trees>
#include <Phyl/PhylogeneticsApplicationTools.h>
#include <Phyl/PhylipDistanceMatrixFormat.h>

using namespace bpp;

void help()
{
  *ApplicationTools::message << "__________________________________________________________________________" << endl;
  *ApplicationTools::message << "bppphysamp arg1=value1 arg2=value2 etc" << endl;
  *ApplicationTools::message << "or" << endl;
  *ApplicationTools::message << "bppphylosamp param=optionfile" << endl;
  *ApplicationTools::message << "__________________________________________________________________________" << endl;
}

int main(int args, char ** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*           Bio++ Phylogenetic Sampler, version 0.1              *" << endl;
  cout << "* Author: J. Dutheil                        Last Modif. 02/03/09 *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;
  
  if(args == 1)
  {
    help();
    exit(0);
  }
  
  try {

  ApplicationTools::startTimer();

  cout << "Parsing options:" << endl;
  
  map<string, string> params = AttributesTools::parseOptions(args, argv);

  //Get sequences:
  Alphabet* alphabet      = SequenceApplicationTools::getAlphabet(params);
  SequenceContainer* seqs = SequenceApplicationTools::getSequenceContainer(alphabet, params);

  string inputMethod = ApplicationTools::getStringParameter("input.method", params, "tree");
  ApplicationTools::displayResult("Input method", inputMethod);

  DistanceMatrix* dist = NULL;
  if(inputMethod == "tree")
  {
    Tree* tree = PhylogeneticsApplicationTools::getTree(params);
    dist = TreeTemplateTools::getDistanceMatrix(*tree);
  }
  else if(inputMethod == "matrix")
  {
    string distPath = ApplicationTools::getAFilePath("input.matrix", params, true, true);
    PhylipDistanceMatrixFormat matIO;
    dist = matIO.read(distPath);
  }
  else throw Exception("Unknown input method: " + inputMethod);

  string critMeth = ApplicationTools::getStringParameter("choice_criterion", params, "length");
  ApplicationTools::displayResult("Sequence choice criterion", critMeth);

  double threshold = ApplicationTools::getDoubleParameter("threshold", params, 0.01);
  ApplicationTools::displayResult("Distance threshold", threshold);

  //Compute lengths:
  vector<string> seqNames;
  vector<unsigned int> seqLen(dist->size());
  string name;
  for(unsigned int i = 0; i < dist->size(); i++)
  {
    name = dist->getName(i);
    if(critMeth == "length.complete")
      seqLen[i] = SequenceTools::getNumberOfCompleteSites(*seqs->getSequence(name));
    else
      seqLen[i] = SequenceTools::getNumberOfSites(*seqs->getSequence(name));
    seqNames.push_back(name);
  }

  unsigned int rm = 0;
  for(unsigned int i = 0; i < dist->size()-1; i++)
  {
    for(unsigned int j = i+1; j < dist->size(); j++)
    {
      //if((*dist)(i,j) < -log(0.)) cout << (*dist)(i, j) << endl;
      if((*dist)(i, j) <= threshold)
      {
        //We need to chose between the two sequences:
        if(critMeth == "length" || critMeth == "length.complete")
        {
          if(seqLen[i] > seqLen[j]) rm = j; else rm = i;
        }
        else if(critMeth == "random")
          if(RandomTools::flipCoin()) rm = j; else rm = i;
        else throw Exception("Unknown cirterion: " + critMeth);

        //Remove sequence in list:
        unsigned int pos = VectorTools::which(seqNames, dist->getName(rm));
        ApplicationTools::displayResult("Remove sequence", seqNames[pos]);
        seqNames.erase(seqNames.begin() + pos); 
        
        //Ignore all distances from this sequence:
        for(unsigned int k = 0; k < dist->size(); k++)
        {
          (*dist)(rm, k) = (*dist)(k, rm) = -log(0.); 
        }
      }
    }
  }
  ApplicationTools::displayResult("Number of sequences kept:", seqNames.size());

  //Write sequences to file:
  VectorSequenceContainer vsc(alphabet);
  for(unsigned int i = 0; i < seqNames.size(); i++)
    vsc.addSequence(* seqs->getSequence(seqNames[i]));
   
  SequenceApplicationTools::writeSequenceFile(vsc, params);
  cout << "Bio++ PhyloSampler's done. Bye." << endl;
  ApplicationTools::displayTime("Total execution time:");

  }
  catch(exception & e)
  {
    cout << endl;
    cout << "_____________________________________________________" << endl;
    cout << "ERROR!!!" << endl;
    cout << e.what() << endl;
    exit(-1);
  }

  return (0);
}

