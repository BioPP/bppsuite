//
// File: bppSeqGen.cpp
// Created by: Julien Dutheil
// Created on: Oct Mon 24 18:50 2005
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

using namespace std;

// From bpp-core:
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Numeric/Number.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Numeric/Random/RandomTools.h>
#include <Bpp/Text/KeyvalTools.h>
#include <Bpp/App/NumCalcApplicationTools.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SequenceContainerTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/SequenceTools.h>

// From bpp-phyl:
#include <Bpp/Phyl/TreeTemplate.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Simulation.all>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Model/FrequenciesSet/MvaFrequenciesSet.h>
#include <Bpp/Phyl/Io/Newick.h>

using namespace bpp;

/**
 * @brief Read trees from an input file, with segment annotations.
 */
void readTrees(ifstream& file, vector<Tree*>& trees, vector<double>& pos) throw (Exception)
{
  string line = "";
  double begin, end;
  string::size_type index1, index2, index3;
  double previousPos = 0;
  pos.push_back(0);
  string newickStr;
  while (!file.eof())
  {
    string tmp = TextTools::removeSurroundingWhiteSpaces(FileTools::getNextLine(file));
    if (tmp.size() == 0 || tmp.substr(0, 1) == "#") continue;
    line += tmp;
        
    index1 = line.find_first_of(" \t");
    if (index1 == string::npos) throw Exception("Error when parsing tree file: now begining position.");
    index2 = line.find_first_of(" \t", index1 + 1);
    if (index2 == string::npos) throw Exception("Error when parsing tree file: now ending position.");
    begin  = TextTools::toDouble(line.substr(0, index1));
    end    = TextTools::toDouble(line.substr(index1 + 1, index2 - index1 - 1));
    index3 = line.find_first_of(";", index2 + 1);
    while (index3 == string::npos)
    {
      if (file.eof()) throw Exception("Error when parsing tree file: incomplete tree.");
      line += FileTools::getNextLine(file);
      index3 = line.find_first_of(";", index3);
    }
    newickStr = line.substr(index2 + 1, index3 - index2);
    TreeTemplate<Node>* t = TreeTemplateTools::parenthesisToTree(newickStr);
    if (trees.size() > 0)
    {
      //Check leave names:
      if (!VectorTools::haveSameElements(t->getLeavesNames(), trees[trees.size()-1]->getLeavesNames()))
        throw Exception("Error: all trees must have the same leaf names.");
    }
    trees.push_back(t);
    if(begin != previousPos) throw Exception("Error when parsing tree file: segments do not match: " + TextTools::toString(begin) + " against " + TextTools::toString(previousPos) + ".");
    pos.push_back(end);
    previousPos = end;

    line = line.substr(index3 + 1);
  }
}

void help()
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
  (*ApplicationTools::message << "bppseqgen parameter1_name=parameter1_value").endLine();
  (*ApplicationTools::message << "      parameter2_name=parameter2_value ... param=option_file").endLine();
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "  Refer to the Bio++ Program Suite Manual for a list of available options.").endLine();
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
}

int main(int args, char ** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*            Bio++ Sequence Generator, version 2.2.0             *" << endl;
  cout << "*                                                                *" << endl;
  cout << "* Authors: J. Dutheil                                            *" << endl;
  cout << "*          B. Boussau                       Last Modif. 25/09/14 *" << endl;
  cout << "*          L. Gueguen                                            *" << endl;
  cout << "*          M. Groussin                                           *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;
  
  if(args == 1)
  {
    help();
    return 0;
  }
  
  try {

  BppApplication bppseqgen(args, argv, "BppSeqGen");
  bppseqgen.startTimer();

  Alphabet* alphabet = SequenceApplicationTools::getAlphabet(bppseqgen.getParams(), "", false);
  auto_ptr<GeneticCode> gCode;
  CodonAlphabet* codonAlphabet = dynamic_cast<CodonAlphabet*>(alphabet);
  if (codonAlphabet) {
    string codeDesc = ApplicationTools::getStringParameter("genetic_code", bppseqgen.getParams(), "Standard", "", true, true);
    ApplicationTools::displayResult("Genetic Code", codeDesc);
      
    gCode.reset(SequenceApplicationTools::getGeneticCode(codonAlphabet->getNucleicAlphabet(), codeDesc));
  }


  /**************************/
  /* Trees                  */
  /**************************/

  
  vector<Tree*> trees;
  vector<double> positions;
  string inputTrees = ApplicationTools::getStringParameter("input.tree.method", bppseqgen.getParams(), "single", "", true, false);
  if (inputTrees == "single")
  {
    trees.push_back(PhylogeneticsApplicationTools::getTree(bppseqgen.getParams()));
    positions.push_back(0);
    positions.push_back(1);
    ApplicationTools::displayResult("Number of leaves", TextTools::toString(trees[0]->getNumberOfLeaves()));
    string treeWIdPath = ApplicationTools::getAFilePath("output.tree_ids.file", bppseqgen.getParams(), false, false);
    if (treeWIdPath != "none")
    {
      TreeTemplate<Node> ttree(*trees[0]);
      vector<Node*> nodes = ttree.getNodes();
      for (size_t i = 0; i < nodes.size(); i++)
      {
        if (nodes[i]->isLeaf())
          nodes[i]->setName(TextTools::toString(nodes[i]->getId()) + "_" + nodes[i]->getName());
        else
          nodes[i]->setBranchProperty("NodeId", BppString(TextTools::toString(nodes[i]->getId())));
      }
      Newick treeWriter;
      treeWriter.enableExtendedBootstrapProperty("NodeId");
      ApplicationTools::displayResult("Writing tagged tree to", treeWIdPath);
      treeWriter.write(ttree, treeWIdPath);
      delete trees[0];
      cout << "BppSegGen's done." << endl;
      exit(0);
    }
  }
  else if (inputTrees == "multiple")
  {
    string treesPath = ApplicationTools::getAFilePath("input.tree.file", bppseqgen.getParams(), false, true);
    ApplicationTools::displayResult("Trees file", treesPath);
    ifstream treesFile(treesPath.c_str(), ios::in);
    readTrees(treesFile, trees, positions);
  }
  else throw Exception("Unknown input.tree.method option: " + inputTrees);


  /**********************************/
  /*  Models                        */
  /**********************************/

  
  string nhOpt = ApplicationTools::getStringParameter("nonhomogeneous", bppseqgen.getParams(), "no", "", true, false);
  ApplicationTools::displayResult("Heterogeneous model", nhOpt);

  SubstitutionModelSet* modelSet = 0;

  //Homogeneous case:
  if (nhOpt == "no")
  {
    SubstitutionModel* model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, gCode.get(), 0, bppseqgen.getParams());
    FrequenciesSet* fSet = new FixedFrequenciesSet(model->getAlphabet(), model->getFrequencies());
    modelSet = SubstitutionModelSetTools::createHomogeneousModelSet(model, fSet, trees[0]);
  }
  //Galtier-Gouy case:
  else if (nhOpt == "one_per_branch")
  {
    if(inputTrees == "multiple")
      throw Exception("Multiple input trees cannot be used with non-homogeneous simulations.");
    SubstitutionModel* model = 0;
    string modelName = ApplicationTools::getStringParameter("model", bppseqgen.getParams(), "");
    if (!TextTools::hasSubstring(modelName,"COaLA"))
      model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, gCode.get(), 0, bppseqgen.getParams());
    else
    {
      //COaLA model
      VectorSiteContainer* allSitesAln = 0;
      allSitesAln = SequenceApplicationTools::getSiteContainer(alphabet, bppseqgen.getParams());
      model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, gCode.get(), allSitesAln, bppseqgen.getParams());
    }

    vector<string> globalParameters = ApplicationTools::getVectorParameter<string>("nonhomogeneous_one_per_branch.shared_parameters", bppseqgen.getParams(), ',', "");
    vector<double> rateFreqs;
    if (model->getNumberOfStates() != alphabet->getSize())
    {
      //Markov-Modulated Markov Model...
      unsigned int n = static_cast<unsigned int>(model->getNumberOfStates() / alphabet->getSize());
      rateFreqs = vector<double>(n, 1./static_cast<double>(n)); // Equal rates assumed for now, may be changed later (actually, in the most general case,
                                                   // we should assume a rate distribution for the root also!!!  
    }
    FrequenciesSet* rootFreqs = PhylogeneticsApplicationTools::getRootFrequenciesSet(alphabet, gCode.get(), 0, bppseqgen.getParams(), rateFreqs);
    string freqDescription = ApplicationTools::getStringParameter("nonhomogeneous.root_freq", bppseqgen.getParams(), "Full(init=observed)");
    if (freqDescription.substr(0,10) == "MVAprotein")
    {
      dynamic_cast<MvaFrequenciesSet*>(rootFreqs)->setModelName("MVAprotein");   
      dynamic_cast<MvaFrequenciesSet*>(rootFreqs)->initSet(dynamic_cast<CoalaCore*>(model));      
    }
    modelSet = SubstitutionModelSetTools::createNonHomogeneousModelSet(model, rootFreqs, trees[0], globalParameters); 
  }
  //General case:
  else if (nhOpt == "general")
  {
    if (inputTrees == "multiple")
      throw Exception("Multiple input trees cannot be used with non-homogeneous simulations.");
    string modelName = ApplicationTools::getStringParameter("model1",bppseqgen.getParams(),"");
    if (!TextTools::hasSubstring(modelName,"COaLA"))
     modelSet = PhylogeneticsApplicationTools::getSubstitutionModelSet(alphabet, gCode.get(), 0, bppseqgen.getParams());
    else
    {
      //COaLA model
      VectorSiteContainer* allSitesAln = 0;
      allSitesAln = SequenceApplicationTools::getSiteContainer(alphabet, bppseqgen.getParams());
      modelSet = PhylogeneticsApplicationTools::getSubstitutionModelSet(alphabet, gCode.get(), allSitesAln, bppseqgen.getParams());
    } 
  }
  else throw Exception("Unknown non-homogeneous option: " + nhOpt);

  if (dynamic_cast<MixedSubstitutionModelSet*>(modelSet))
    throw Exception("Non-homogeneous mixed substitution sequence generation not implemented, sorry!");

  /*******************************************/
  /*     Starting sequence                   */
  /*******************************************/

  DiscreteDistribution* rDist = 0;
  NonHomogeneousSequenceSimulator* seqsim = 0;
  SiteContainer* sites = 0;
  size_t nbSites = 0;

  string infosFile = ApplicationTools::getAFilePath("input.infos", bppseqgen.getParams(), false, true);

  bool withStates = false;
  bool withRates = false;
  vector<size_t> states;
  vector<double> rates;
  
  if (infosFile != "none")
  {
    ApplicationTools::displayResult("Site information", infosFile);
    ifstream in(infosFile.c_str());
    DataTable* infos = DataTable::read(in, "\t");
    nbSites = infos->getNumberOfRows();
    ApplicationTools::displayResult("Number of sites", TextTools::toString(nbSites));
    string rateCol = ApplicationTools::getStringParameter("input.infos.rates", bppseqgen.getParams(), "pr", "", true, true);
    string stateCol = ApplicationTools::getStringParameter("input.infos.states", bppseqgen.getParams(), "none", "", true, true);
    withRates = rateCol != "none";
    withStates = stateCol != "none";
  
    if (withRates)
    {
      rDist = new ConstantRateDistribution();
      rates.resize(nbSites);
      vector<string> ratesStrings = infos->getColumn(rateCol);
      for (size_t i = 0; i < nbSites; i++)
      {
        rates[i] = TextTools::toDouble(ratesStrings[i]);
      }
    }
    if (withStates)
    {
      vector<string> ancestralStates = infos->getColumn(stateCol);
      states.resize(nbSites);
      for (size_t i = 0; i < nbSites; i++)
      {
        int alphabetState = alphabet->charToInt(ancestralStates[i]);
        //If a generic character is provided, we pick one state randomly from the possible ones:
        if (alphabet->isUnresolved(alphabetState))
          alphabetState = RandomTools::pickOne<int>(alphabet->getAlias(alphabetState));
        states[i] = RandomTools::pickOne<size_t>(modelSet->getModelStates(alphabetState));
      }

      string siteSet = ApplicationTools::getStringParameter("input.site.selection", bppseqgen.getParams(), "none", "", true, 1);

      if (siteSet != "none")
      {
        vector<size_t> vSite;
        try {
          vector<int> vSite1 = NumCalcApplicationTools::seqFromString(siteSet);
          for (size_t i = 0; i < vSite1.size(); ++i){
            int x = (vSite1[i] >= 0 ? vSite1[i] : static_cast<int>(nbSites) + vSite1[i]);
            if (x >= 0)
              vSite.push_back(static_cast<size_t>(x));
            else
              throw Exception("SequenceApplicationTools::getSiteContainer(). Incorrect negative index: " + TextTools::toString(x));
          }
        }
        catch (Exception& e)
        {
          string seln;
          map<string, string> selArgs;
          KeyvalTools::parseProcedure(siteSet, seln, selArgs);
          if (seln == "Sample")
          {
            size_t n = ApplicationTools::getParameter<size_t>("n", selArgs, nbSites, "", true, 1);
            bool replace = ApplicationTools::getBooleanParameter("replace", selArgs, false, "", true, 1);

            vSite.resize(n);
            vector<size_t> vPos;
            for (size_t p = 0; p < nbSites; ++p)
              vPos.push_back(p);
          
            RandomTools::getSample(vPos, vSite, replace);
          }
        }
        
        nbSites = vSite.size();

        vector<size_t> newStates(nbSites);
        vector<double> newRates(nbSites);

        for (size_t ni = 0; ni < nbSites; ++ni)
        {
          newStates[ni] = states[vSite[ni]];
          newRates[ni]  = rates[vSite[ni]];
        }

        states = newStates;
        rates = newRates;
      }
    }
  }
  else
  {
    try {
      VectorSiteContainer* allSeq = 0;
      allSeq = SequenceApplicationTools::getSiteContainer(alphabet, bppseqgen.getParams());
      
      if (allSeq->getNumberOfSequences() > 0)
      {  
        Sequence* pseq = SequenceTools::getSequenceWithCompleteSites(allSeq->getSequence(0));
        
        nbSites = pseq->size();
        states.resize(nbSites);
        withStates = true;
        
	for (size_t i = 0; i < nbSites; ++i) {
          states[i] = RandomTools::pickOne(modelSet->getModelStates((*pseq)[i]));
        }
        ApplicationTools::displayResult("Number of sites", TextTools::toString(nbSites));
        
        delete pseq;
      }
    }
    catch (Exception& e)
    {
    }
      
  }

  if (rDist == 0)
  {
    if (modelSet->getNumberOfStates() > modelSet->getAlphabet()->getSize())
    {
      //Markov-modulated Markov model!
      rDist = new ConstantRateDistribution();
    }
    else
    {
      rDist = PhylogeneticsApplicationTools::getRateDistribution(bppseqgen.getParams());
    }
  }


  if (nbSites == 0)
    nbSites = ApplicationTools::getParameter<size_t>("number_of_sites", bppseqgen.getParams(), 100);
  
  /*******************/
  /* Simulations     */
  /*******************/

  if (withStates)
  {
    if (trees.size() == 1)
    {
      seqsim = new NonHomogeneousSequenceSimulator(modelSet, rDist, trees[0]);
      ApplicationTools::displayTask("Perform simulations");
      if (withRates)
        if (withStates)
          sites = SequenceSimulationTools::simulateSites(*seqsim, rates, states);
        else
          sites = SequenceSimulationTools::simulateSites(*seqsim, rates);
      else
        if (withStates){
          sites = SequenceSimulationTools::simulateSites(*seqsim, states);
        }      
        else
          throw Exception("Error! Info file should contain either site specific rates of ancestral states or both.");
      
      delete seqsim;    
    }
    else
    {
      ApplicationTools::displayTask("Perform simulations", true);
      ApplicationTools::displayGauge(0, trees.size() - 1, '=');
      seqsim = new NonHomogeneousSequenceSimulator(modelSet, rDist, trees[0]);
      ptrdiff_t previousPos = 0;
      ptrdiff_t currentPos = static_cast<ptrdiff_t>(round(positions[1]*static_cast<double>(nbSites)));
      vector<double> tmpRates;
      if (withRates)
        tmpRates = vector<double>(rates.begin() + previousPos, rates.begin() + currentPos);
      vector<size_t> tmpStates;
      if (withStates)
        tmpStates = vector<size_t>(states.begin() + previousPos, states.begin() + currentPos);
      SequenceContainer* tmpCont1 = 0;
      if (withRates)
        if (withStates)
          tmpCont1 = SequenceSimulationTools::simulateSites(*seqsim, tmpRates, tmpStates);
        else
          tmpCont1 = SequenceSimulationTools::simulateSites(*seqsim, tmpRates);
      else
        if (withStates)
          tmpCont1 = SequenceSimulationTools::simulateSites(*seqsim, tmpStates);
        else
          throw Exception("Error! Info file should contain either site specific rates of ancestral states or both.");
      previousPos = currentPos;
      delete seqsim;
      
      for(size_t i = 1; i < trees.size(); i++)
      {
        ApplicationTools::displayGauge(i, trees.size() - 1, '=');
        seqsim = new NonHomogeneousSequenceSimulator(modelSet, rDist, trees[i]);
        currentPos = static_cast<ptrdiff_t>(round(positions[i+1]) * static_cast<double>(nbSites));
        if (withRates)
          tmpRates = vector<double>(rates.begin() + previousPos + 1, rates.begin() + currentPos);
        if (withStates)
          tmpStates = vector<size_t>(states.begin() + previousPos + 1, states.begin() + currentPos);
        SequenceContainer* tmpCont2 = 0;
        if (withRates)
          if (withStates)
            tmpCont2 = SequenceSimulationTools::simulateSites(*seqsim, tmpRates, tmpStates);
          else     
            tmpCont2 = SequenceSimulationTools::simulateSites(*seqsim, tmpRates);
        else
          if (withStates)
            tmpCont2 = SequenceSimulationTools::simulateSites(*seqsim, tmpStates);
          else
            throw Exception("Error! Info file should contain either site specific rates of ancestral states or both.");
        previousPos = currentPos;
        delete seqsim;
        VectorSequenceContainer* mergedCont = new VectorSequenceContainer(alphabet);
        SequenceContainerTools::merge(*tmpCont1, *tmpCont2, *mergedCont);
        delete tmpCont1;
        delete tmpCont2;
        tmpCont1 = mergedCont;
      }
      sites = new VectorSiteContainer(*tmpCont1);
      delete tmpCont1;
    }
    ApplicationTools::displayTaskDone();
  }
  else
  {
    if (modelSet->getNumberOfStates() > modelSet->getAlphabet()->getSize())
    {
      //Markov-modulated Markov model!
      rDist = new ConstantRateDistribution();
    }
    else
    {
      rDist = PhylogeneticsApplicationTools::getRateDistribution(bppseqgen.getParams());
    }
    
    if (trees.size() == 1)
    {
      seqsim = new NonHomogeneousSequenceSimulator(modelSet, rDist, trees[0]);
      ApplicationTools::displayResult("Number of sites", TextTools::toString(nbSites));
      ApplicationTools::displayTask("Perform simulations");
      sites = seqsim->simulate(nbSites);
      ApplicationTools::displayTaskDone();
    }
    else
    {
      ApplicationTools::displayTask("Perform simulations", true);
      ApplicationTools::displayGauge(0, trees.size() - 1, '=');
      seqsim = new NonHomogeneousSequenceSimulator(modelSet, rDist, trees[0]);
      size_t previousPos = 0;
      size_t currentPos = static_cast<unsigned int>(round(positions[1] * static_cast<double>(nbSites)));
      SequenceContainer* tmpCont1 = seqsim->simulate(currentPos - previousPos);
      previousPos = currentPos;
      delete seqsim;
 
      for (size_t i = 1; i < trees.size(); i++)
      {
        ApplicationTools::displayGauge(i, trees.size() - 1, '=');
        seqsim = new NonHomogeneousSequenceSimulator(modelSet, rDist, trees[i]);
        currentPos = static_cast<unsigned int>(round(positions[i+1] * static_cast<double>(nbSites)));
        SequenceContainer* tmpCont2 = seqsim->simulate(currentPos - previousPos);
        previousPos = currentPos;
        delete seqsim;
        VectorSequenceContainer* mergedCont = new VectorSequenceContainer(alphabet);
        SequenceContainerTools::merge(*tmpCont1, *tmpCont2, *mergedCont);
        delete tmpCont1;
        delete tmpCont2;
        tmpCont1 = mergedCont;
      }
      sites = new VectorSiteContainer(*tmpCont1);
      ApplicationTools::displayTaskDone();
      delete tmpCont1;
    }
  }
  
  // Write to file:
  SequenceApplicationTools::writeAlignmentFile(*sites, bppseqgen.getParams());
  
  delete alphabet;
  for (size_t i = 0; i < trees.size(); i++)
    delete trees[i];
  delete rDist;

  bppseqgen.done();
  
}
  catch (exception& e)
  {
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}

