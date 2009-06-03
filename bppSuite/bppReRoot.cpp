//
// File: bppReRoot.cpp
// Created by: Celine Scornavacca
// Created on: Jan Tue 15 18:15 2008
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

using namespace std;

// From PhylLib:
#include <Phyl/IOTree.h>
#include <Phyl/Newick.h>
#include <Phyl/Tree.h>
#include <Phyl/Node.h>
#include <Phyl/TreeExceptions.h>
#include <Phyl/TreeTemplateTools.h>
#include <Phyl/PhylogeneticsApplicationTools.h>

// From NumCalc:
#include <NumCalc/VectorTools.h>

// From Utils:
#include <Utils/FileTools.h>
#include <Utils/TextTools.h>
#include <Utils/StringTokenizer.h>
#include <Utils/ApplicationTools.h>
#include <Utils/AttributesTools.h>

using namespace bpp;

typedef TreeTemplate<Node> MyTree;


void help()
{
  *ApplicationTools::message << "bppreroot parameter1_name=parameter1_value"    << endl;
  *ApplicationTools::message << "      parameter2_name=parameter2_value ... param=option_file" << endl;
  *ApplicationTools::message << endl;
  *ApplicationTools::message << "  Refer to the Bio++ Program Suite Manual for a list of available options." << endl;
}




int main(int args, char ** argv)
{
  
  cout << "******************************************************************" << endl;
  cout << "*                  Bio++ ReRoot, version 0.1.3                   *" << endl;
  cout << "* Author: C. Scornavacca                    Created     15/01/08 *" << endl;
  cout << "*                                           Last Modif. 03/06/09 *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  if(args == 1)
  {
    help();
    return 0;
  }
  
  try {
  
  ApplicationTools::startTimer();

  cout << "Parsing options:" << endl;
  
  map<string, string> params = AttributesTools::parseOptions(args, argv);

  Newick newick;
  string listPath = ApplicationTools::getAFilePath("input.list.file", params);
  ApplicationTools::displayResult("Input list file", listPath);
  if(listPath == "none") throw Exception("You must provide an input tree list file.");
     
  string outgroupsPath = ApplicationTools::getAFilePath("outgroups.file", params);
  ApplicationTools::displayResult("Outgroups file", outgroupsPath);
  if(outgroupsPath == "none") throw Exception("You must provide an outgroup list file.");
       
  string outputPath = ApplicationTools::getAFilePath("output.trees.file", params, true, false);
  ApplicationTools::displayResult("Output file", outputPath);
  if(outputPath == "none") throw Exception("You must provide an output file.");
     
  bool printOption = ApplicationTools::getBooleanParameter("print.option", params, false);
  bool tryAgain = ApplicationTools::getBooleanParameter("tryAgain.option", params, true);
          
  vector<Tree*> tempTrees;
  vector<MyTree*> trees;

  //ApplicationTools::displayResult("Number of trees found", TextTools::toString(trees.size()));
            
  const string path = outgroupsPath;  
  ifstream file(path.c_str(), ios::in);
  string temp, description, taxon;

  vector < vector<string> > levelOutgroup;
    
  //Reading outgroup levels  
  while (!file.eof()) 
  {
    vector <string> tempTaxa;
    getline(file, temp, '\n');  
    StringTokenizer line = StringTokenizer::StringTokenizer(temp, " ,"); 
    while (line.hasMoreToken())
    {
      tempTaxa.push_back(line.nextToken());  
    }
    levelOutgroup.push_back(tempTaxa);
  }
  file.close();  
  
  const string path2 = listPath;  
  ifstream treePath(path2.c_str(), ios::in);  

  if(! treePath) { throw IOException ("Newick::read: failed to read from stream"); }

  string temp2, description2;// Initialization
  string::size_type index;  
  
  int k = 0;
  
  while(!treePath.eof())
  {
    k++;
    bool printOrNot =true;
    Tree * tempTree = NULL;

    getline(treePath, temp2, '\n');  // Copy current line in temporary string
   
    index = temp2.find(";");
    if(temp2 != "")
    {
      if(index != string::npos)
      {
        description2 += temp2.substr(0, index + 1);
        tempTree = TreeTemplateTools::parenthesisToTree(description2);    
        description2 = temp2.substr(index + 1);       
      }
      else description2 += temp;

      MyTree * tree = dynamic_cast <MyTree* >(tempTree);
      //ApplicationTools::displayGauge(tr, trees.size() - 1, '=');

      vector<string> leavesTree;      
      leavesTree = (* tree).getLeavesNames();  
  
      unsigned int numNodes = tree->getNumberOfNodes() - 1;
      unsigned int numNodeWithBranchLength = 0;
      vector<Node *>  nodes = tree->getNodes();
      for(unsigned int i = 0; i < nodes.size(); i++)
      {
        if(nodes[i]->hasDistanceToFather())
          numNodeWithBranchLength++;
      }
      if((numNodes != numNodeWithBranchLength) && (numNodeWithBranchLength != 0))\
      {
        cout << "Could not execute due to a source tree with missing branch lengths \n(reminder: a source tree must either have no branch length, either length for all branches\n";
        exit(-1);
      }
      vector<string> outGroup;
      bool found = false;
      bool analyseOutgroupLevel = true;
      for (unsigned int t = 0; t < levelOutgroup.size() && analyseOutgroupLevel; t++)
      {      
        outGroup.clear();
        vector<string>::iterator Iterator;  
        for(Iterator = levelOutgroup[t].begin(); Iterator != levelOutgroup[t].end(); Iterator++ )
        {
          if(VectorTools::contains(leavesTree, *Iterator))
          {
            outGroup.push_back(*Iterator);
          }
        }
        if(outGroup.size() > 0)
        {
          vector<string> remainingTaxa;
          VectorTools::diff(leavesTree, outGroup, remainingTaxa);
          if(remainingTaxa.size() > 0)
          {
            tree->newOutGroup(* tree->getNode(remainingTaxa[0]));
            Node * newRoot = tree->getNode(outGroup[0]);
            vector<string>  tempLeaves = TreeTemplateTools::getLeavesNames(* newRoot);
           
            while(newRoot->hasFather() && !(VectorTools::containsAll(tempLeaves, outGroup)))
            {   
              newRoot = newRoot->getFather();
              tempLeaves = TreeTemplateTools::getLeavesNames(* newRoot);
            }
          
            tempLeaves = TreeTemplateTools::getLeavesNames(* newRoot);
            std::sort(tempLeaves.begin(), tempLeaves.end());
      
            if(tempLeaves.size() == outGroup.size())
            {
              tree->newOutGroup(* newRoot);
              found = true;
              analyseOutgroupLevel = false;
            }
            else
            {
              bool monophylOk = true;

              for(unsigned f = 0; f < newRoot->getNumberOfSons() && monophylOk; f++)
              {
                vector<string>  tempLeaves = TreeTemplateTools::getLeavesNames(* newRoot->getSon(f));
                vector<string> diff;
                VectorTools::diff(outGroup, tempLeaves, diff);

                unsigned int difference = diff.size();
                if(!( (difference == 0) || (difference == tempLeaves.size()) ) )
                {
                  //The proposed outgroup is not monophyletic. The analysis for this tree is interrupted
                  //No more outgroup are analysed
                  monophylOk = false;
                }
              }
              if(monophylOk)
              {
                vector<string>  tempLeaves = TreeTemplateTools::getLeavesNames(* newRoot);

                std::sort(tempLeaves.begin(), tempLeaves.end());      
                if(tempLeaves.size() != leavesTree.size())
                {
                  MyTree * low = new MyTree(* TreeTemplateTools::cloneSubtree<Node>(* newRoot));
                  tree->newOutGroup(* newRoot);
                  Node * sonUpper;
                  vector<string>  tempLeaves2 = TreeTemplateTools::getLeavesNames(* (tree->getRootNode())->getSon(0));
                  std::sort(tempLeaves2.begin(), tempLeaves2.end());
                  if((VectorTools::vectorIntersection(tempLeaves2,outGroup).size()) !=0)
                  {
                    sonUpper = (tree->getRootNode())->getSon(1);
                  }
                  else
                  {
                    sonUpper = (tree->getRootNode())->getSon(0);
                  }
                  int ident = TreeTools::getMaxId(* low, low->getRootId());
                  vector <Node *> nodesTemp= TreeTemplateTools::getNodes( * sonUpper);
                  for(unsigned int F = 0; F < nodesTemp.size(); F++)
                    ( * nodesTemp[F]).setId(ident + F + 1);
                  low->getRootNode()->addSon(* sonUpper);
                  tree = low;
                }
                //A good outgroup was found

                found = true;
                analyseOutgroupLevel = false;
              }
            }                  
          }
          if(!tryAgain)
            analyseOutgroupLevel = false;
        }    
      }
      if(!found)
      {  
        if(!printOption)
          printOrNot = false;
        else
          printOrNot = true;
        cout << "Sorry but I can't root your tree " << k << " ; or none of the taxa in your list is present in the tree or the outgroup is not monophyletic!\n";
      }
      else
      {
        printOrNot = (true);
        tree->resetNodesId();
      }
      if(printOrNot)
      {
        if(k == 1)
          newick.write(* tree, outputPath, true);
        else
          newick.write(* tree, outputPath, false);
      }  
    

      delete tree;
      //delete(tempTree);
    }
  }
  ApplicationTools::displayTaskDone();
     
  //Write rooted trees:  

     
  for(unsigned int i = 0; i < trees.size(); i++) delete trees[i];
    
  cout << "Bio++ ReRoot's done. Bye." << endl;
      
  ApplicationTools::displayTime("Total execution time:");

  }
  catch(exception & e)
  {
    cout << e.what() << endl;
    return 1;
  }

  return 0;
};







