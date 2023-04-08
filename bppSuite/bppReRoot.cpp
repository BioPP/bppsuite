//
// File: bppReRoot.cpp
// Created by: Celine Scornavacca
// Created on: Jan Tue 15 18:15 2008
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

using namespace std;

// From bpp-core:
#include <Bpp/Version.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/StringTokenizer.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>

// From PhylLib:
#include <Bpp/Phyl/Tree/Tree.h>
#include <Bpp/Phyl/Tree/Node.h>
#include <Bpp/Phyl/Tree/TreeExceptions.h>
#include <Bpp/Phyl/Tree/TreeTemplateTools.h>
// From bpp-phyl:
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>

using namespace bpp;

typedef TreeTemplate<Node> MyTree;


void help()
{
  (*ApplicationTools::message << "bppreroot parameter1_name=parameter1_value").endLine();
  (*ApplicationTools::message << "      parameter2_name=parameter2_value ... param=option_file").endLine();
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "  Refer to the Bio++ Program Suite Manual for a list of available options.").endLine();
}




int main(int args, char ** argv)
{
  
  cout << "******************************************************************" << endl;
  cout << "*                  Bio++ ReRoot, version " << BPP_VERSION << "                   *" << endl;
  cout << "* Author: C. Scornavacca                    Created     15/01/08 *" << endl;
  cout << "*                                           Last Modif. " << BPP_REL_DATE << " *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  if(args == 1)
  {
    help();
    return 0;
  }
  
  try {
  
  BppApplication bppreroot(args, argv, "BppReRoot");
  bppreroot.startTimer();

  Newick newick;
  string listPath = ApplicationTools::getAFilePath("input.list.file", bppreroot.getParams());
  ApplicationTools::displayResult("Input list file", listPath);
  if(listPath == "none") throw Exception("You must provide an input tree list file.");
     
  string outgroupsPath = ApplicationTools::getAFilePath("outgroups.file", bppreroot.getParams());
  ApplicationTools::displayResult("Outgroups file", outgroupsPath);
  if(outgroupsPath == "none") throw Exception("You must provide an outgroup list file.");
       
  string outputPath = ApplicationTools::getAFilePath("output.trees.file", bppreroot.getParams(), true, false);
  ApplicationTools::displayResult("Output file", outputPath);
  if(outputPath == "none") throw Exception("You must provide an output file.");
     
  bool printOption = ApplicationTools::getBooleanParameter("print.option", bppreroot.getParams(), false);
  bool tryAgain = ApplicationTools::getBooleanParameter("tryAgain.option", bppreroot.getParams(), true);
          
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
    StringTokenizer line = StringTokenizer(temp, " ,"); 
    while (line.hasMoreToken())
    {
      tempTaxa.push_back(line.nextToken());  
    }
    levelOutgroup.push_back(tempTaxa);
  }
  file.close();  
  
  const string path2 = listPath;  
  ifstream treePath(path2.c_str(), ios::in);  

  if (!treePath) { throw IOException ("Newick::read: failed to read from stream"); }

  string temp2, description2;// Initialization
  string::size_type index;  
  
  int k = 0;
  
  while (!treePath.eof())
  {
    k++;
    bool printOrNot =true;
    unique_ptr<Tree> tempTree;

    getline(treePath, temp2, '\n');  // Copy current line in temporary string
   
    index = temp2.find(";");
    if (temp2 != "")
    {
      if (index != string::npos)
      {
        description2 += temp2.substr(0, index + 1);
        tempTree = TreeTemplateTools::parenthesisToTree(description2);    
        description2 = temp2.substr(index + 1);       
      }
      else description2 += temp;

      MyTree* tree = dynamic_cast <MyTree* >(tempTree.get());
      //ApplicationTools::displayGauge(tr, trees.size() - 1, '=');

      vector<string> leavesTree;      
      leavesTree = (*tree).getLeavesNames();  
  
      size_t numNodes = tree->getNumberOfNodes() - 1;
      size_t numNodeWithBranchLength = 0;
      vector<Node *>  nodes = tree->getNodes();
      for (size_t i = 0; i < nodes.size(); i++)
      {
        if(nodes[i]->hasDistanceToFather())
          numNodeWithBranchLength++;
      }
      if ((numNodes != numNodeWithBranchLength) && (numNodeWithBranchLength != 0))\
      {
        cout << "Could not execute due to a source tree with missing branch lengths \n(reminder: a source tree must either have no branch length, either length for all branches\n";
        exit(-1);
      }
      vector<string> outGroup;
      bool found = false;
      bool analyseOutgroupLevel = true;
      for (size_t t = 0; t < levelOutgroup.size() && analyseOutgroupLevel; t++)
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
            tree->newOutGroup(tree->getNode(remainingTaxa[0]));
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
              tree->newOutGroup(newRoot);
              found = true;
              analyseOutgroupLevel = false;
            }
            else
            {
              bool monophylOk = true;

              for (size_t f = 0; f < newRoot->getNumberOfSons() && monophylOk; f++)
              {
                tempLeaves = TreeTemplateTools::getLeavesNames(*newRoot->getSon(f));
                vector<string> diff;
                VectorTools::diff(outGroup, tempLeaves, diff);

                size_t difference = diff.size();
                if (!( (difference == 0) || (difference == tempLeaves.size()) ) )
                {
                  //The proposed outgroup is not monophyletic. The analysis for this tree is interrupted
                  //No more outgroup are analysed
                  monophylOk = false;
                }
              }
              if (monophylOk)
              {
                tempLeaves = TreeTemplateTools::getLeavesNames(* newRoot);

                std::sort(tempLeaves.begin(), tempLeaves.end());      
                if (tempLeaves.size() != leavesTree.size())
                {
                  MyTree* low = new MyTree(TreeTemplateTools::cloneSubtree<Node>(* newRoot));
                  tree->newOutGroup(newRoot);
                  Node* sonUpper;
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
                  int ident = TreeTools::getMaxId(*low, low->getRootId());
                  vector <Node *> nodesTemp= TreeTemplateTools::getNodes( * sonUpper);
                  for(size_t F = 0; F < nodesTemp.size(); F++)
                    nodesTemp[F]->setId(ident + static_cast<int>(F + 1));
                  low->getRootNode()->addSon(sonUpper);
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
      if (!found)
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
      if (printOrNot)
      {
        if(k == 1)
          newick.writeTree(* tree, outputPath, true);
        else
          newick.writeTree(* tree, outputPath, false);
      }  

      delete tree;
    }
  }
  ApplicationTools::displayTaskDone();
     
  //Write rooted trees:  
  for (size_t i = 0; i < trees.size(); i++) delete trees[i];
    
  bppreroot.done();
  }
  catch(exception & e)
  {
    cout << e.what() << endl;
    return 1;
  }

  return 0;
};







