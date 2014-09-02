//
// File: PhyloSample.cpp
// Created by: Julien Dutheil
// Created on: Sunday, December 2nd 2007 16:48
//

/*
Copyright or © or Copr. Bio++ Development Team

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

#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Numeric/Random/RandomTools.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet.all>
#include <Bpp/Seq/Container.all>
#include <Bpp/Seq/Io.all>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/SequenceTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

// From bpp-phyl:
#include <Bpp/Phyl/Tree/Tree.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Io/PhylipDistanceMatrixFormat.h>

using namespace bpp;

void help()
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
  (*ApplicationTools::message << "bppphysamp parameter1_name=parameter1_value").endLine();
  (*ApplicationTools::message << "      parameter2_name=parameter2_value ... param=option_file").endLine();
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "  Refer to the Bio++ Program Suite Manual for a list of available options.").endLine();
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
}

class Index {
  public:
    double distance;
    unsigned int i1, i2;

  public:
    Index(double dist, unsigned int i, unsigned int j) : distance(dist), i1(i), i2(j) {}

  public:
    bool operator==(const Index& index) const { return distance == index.distance; }
    bool operator<(const Index& index) const { return distance < index.distance; }
};

class Test {
  private:
    unsigned int pos_;

  public:
    Test(unsigned int pos) : pos_(pos) {}
    
  public:
    bool operator()(const Index& index) { return index.i1 == pos_ || index.i2 == pos_; }
};

int main(int args, char ** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*           Bio++ Phylogenetic Sampler, version 0.3              *" << endl;
  cout << "* Author: J. Dutheil                        Last Modif. 06/06/14 *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;
  
  if(args == 1)
  {
    help();
    return 0;
  }
  
  try {

  BppApplication bppphysamp(args, argv, "BppPhySamp");
  bppphysamp.startTimer();

  //Get sequences:
  Alphabet* alphabet      = SequenceApplicationTools::getAlphabet(bppphysamp.getParams());
  SequenceContainer* seqs = SequenceApplicationTools::getSequenceContainer(alphabet, bppphysamp.getParams());

  string inputMethod = ApplicationTools::getStringParameter("input.method", bppphysamp.getParams(), "tree");
  ApplicationTools::displayResult("Input method", inputMethod);

  auto_ptr< DistanceMatrix > dist;
  auto_ptr< TreeTemplate<Node> > tree;
  if(inputMethod == "tree")
  {
    tree.reset(dynamic_cast<TreeTemplate<Node> *>(PhylogeneticsApplicationTools::getTree(bppphysamp.getParams())));
    dist.reset(TreeTemplateTools::getDistanceMatrix(*tree));
  }
  else if(inputMethod == "matrix")
  {
    string distPath = ApplicationTools::getAFilePath("input.matrix", bppphysamp.getParams(), true, true);
    PhylipDistanceMatrixFormat matIO;
    dist.reset(matIO.read(distPath));
  }
  else throw Exception("Unknown input method: " + inputMethod);

  string deleteMeth = ApplicationTools::getStringParameter("deletion_method", bppphysamp.getParams(), "threshold");
  ApplicationTools::displayResult("Deletion method", deleteMeth);

  string critMeth = ApplicationTools::getStringParameter("choice_criterion", bppphysamp.getParams(), "length");
  ApplicationTools::displayResult("Sequence choice criterion", critMeth);

  //Compute lengths:
  vector<string> seqNames;
  vector<size_t> seqLen(dist->size());
  string name;
  for(size_t i = 0; i < dist->size(); i++)
  {
    name = dist->getName(i);
    if (critMeth == "length.complete")
      seqLen[i] = SequenceTools::getNumberOfCompleteSites(seqs->getSequence(name));
    else
      seqLen[i] = SequenceTools::getNumberOfSites(seqs->getSequence(name));
    seqNames.push_back(name);
  }

  //Sort matrix entries:
  vector<Index> distances;
  for (unsigned int i = 0; i < dist->size()-1; i++)
    for (unsigned int j = i+1; j < dist->size(); j++)
      distances.push_back(Index((*dist)(i, j), i , j));
  sort(distances.begin(), distances.end());

  if (deleteMeth == "random")
  {
    unsigned int sampleSize = ApplicationTools::getParameter<unsigned int>("sample_size", bppphysamp.getParams(), 10);
    ApplicationTools::displayResult("Sample size", sampleSize);
    vector<string> sample(sampleSize);
    RandomTools::getSample(seqNames, sample, false);
    seqNames = sample;
    
    double mini = -log(0.);
    for (unsigned int i =  0; i < seqNames.size() - 1; ++i)
      for (unsigned int j = i + 1; j < seqNames.size(); ++j)
      {
        double d = (*dist)(seqNames[i], seqNames[j]);
        if (d < mini) mini = d;
      }
    ApplicationTools::displayResult("Minimal distance in final data set:", mini);
  }
  else if (deleteMeth == "threshold")
  {
    double threshold = ApplicationTools::getDoubleParameter("threshold", bppphysamp.getParams(), 0.01);
    ApplicationTools::displayResult("Distance threshold", threshold);

    unsigned int rm = 0;
    while (distances[0].distance <= threshold)
    {
      //We need to chose between the two sequences:
      if (critMeth == "length" || critMeth == "length.complete")
      {
        if (seqLen[distances[0].i1] > seqLen[distances[0].i2]) rm = distances[0].i2;
        else rm = distances[0].i1;
      }
      else if (critMeth == "random")
      {
        if (RandomTools::flipCoin()) rm = distances[0].i2;
        else rm = distances[0].i1;
      }
      else throw Exception("Unknown criterion: " + critMeth);

      //Remove sequence in list:
      size_t pos = VectorTools::which(seqNames, dist->getName(rm));
      ApplicationTools::displayResult("Remove sequence", seqNames[pos]);
      seqNames.erase(seqNames.begin() + static_cast<ptrdiff_t>(pos)); 
        
      //Ignore all distances from this sequence:
      remove_if(distances.begin(), distances.end(), Test(rm));
      if (distances.size() == 0)
        throw Exception("Error, all sequences have been removed with this criterion!");
    }
    ApplicationTools::displayResult("Number of sequences kept:", seqNames.size());
  }
  else if (deleteMeth == "sample")
  {
    unsigned int sampleSize = ApplicationTools::getParameter<unsigned int>("sample_size", bppphysamp.getParams(), 10);
    ApplicationTools::displayResult("Sample size", sampleSize);
    
    unsigned int rm = 0;
    while (seqNames.size() > sampleSize)
    {
      //We need to chose between the two sequences:
      if (critMeth == "length" || critMeth == "length.complete")
      {
        if (seqLen[distances[0].i1] > seqLen[distances[0].i2]) rm = distances[0].i2;
        else rm = distances[0].i1;
      }
      else if (critMeth == "random")
      {
        if (RandomTools::flipCoin()) rm = distances[0].i2;
        else rm = distances[0].i1;
      }
      else throw Exception("Unknown criterion: " + critMeth);

      //Remove sequence in list:
      size_t pos = VectorTools::which(seqNames, dist->getName(rm));
      ApplicationTools::displayResult("Remove sequence", seqNames[pos]);
      seqNames.erase(seqNames.begin() + static_cast<ptrdiff_t>(pos)); 
        
      //Ignore all distances from this sequence:
      remove_if(distances.begin(), distances.end(), Test(rm));
    }
    ApplicationTools::displayResult("Minimal distance in final data set:", distances[0].distance);
  }
  else throw Exception("Unknown deletion method: " + deleteMeth + ".");

  //Write sequences to file:
  AlignedSequenceContainer asc(alphabet);
  for (size_t i = 0; i < seqNames.size(); i++)
    asc.addSequence(seqs->getSequence(seqNames[i]));
   
  SequenceApplicationTools::writeAlignmentFile(asc, bppphysamp.getParams());

  //Write tree file:
  if (ApplicationTools::getStringParameter("output.tree.file", bppphysamp.getParams(), "None") != "None") {
    for (size_t i = 0; i < seqNames.size(); ++i) {
      TreeTemplateTools::dropLeaf(*tree, seqNames[i]);
    }
    PhylogeneticsApplicationTools::writeTree(*tree, bppphysamp.getParams(), "output.", "", true, true, false);
  }

  bppphysamp.done();
  }
  catch (exception& e)
  {
    cout << endl;
    cout << "_____________________________________________________" << endl;
    cout << "ERROR!!!" << endl;
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}

