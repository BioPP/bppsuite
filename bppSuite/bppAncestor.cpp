//
// File: bppAncestor.cpp
// Created by: Julien Dutheil
// Created on: Sep Wed 10 14:14 2008
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

#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/AutoParameter.h>

// From SeqLib:
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SequenceContainerTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

// From PhylLib:
#include <Bpp/Phyl/Tree.h>
#include <Bpp/Phyl/Likelihood.all>
#include <Bpp/Phyl/PatternTools.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Model/MarkovModulatedSubstitutionModel.h>
#include <Bpp/Phyl/Model/SubstitutionModelSet.h>
#include <Bpp/Phyl/Model/SubstitutionModelSetTools.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Io/Newick.h>

using namespace bpp;

/******************************************************************************/

void help()
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
  (*ApplicationTools::message << "bppancestor parameter1_name=parameter1_value ").endLine();
  (*ApplicationTools::message << "      parameter2_name=parameter2_value ... param=option_file").endLine();
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "  Refer to the Bio++ Program Suite Manual for a list of available options.").endLine();
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
}

int main(int args, char ** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*     Bio++ Ancestral Sequence Reconstruction, version 0.5.0     *" << endl;
  cout << "* Authors: J. Dutheil                       Created on: 10/09/08 *" << endl;
  cout << "*          B. Boussau                       Last Modif: 17/06/11 *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  if (args == 1)
  {
    help();
    return 0;
  }
  
  try {

  BppApplication bppancestor(args, argv, "BppAncestor");
  bppancestor.startTimer();

  Alphabet* alphabet = SequenceApplicationTools::getAlphabet(bppancestor.getParams(), "", false);

  VectorSiteContainer* allSites = SequenceApplicationTools::getSiteContainer(alphabet, bppancestor.getParams());
  
  VectorSiteContainer* sites = SequenceApplicationTools::getSitesToAnalyse(* allSites, bppancestor.getParams(), "", true, false);
  delete allSites;

  ApplicationTools::displayResult("Number of sequences", TextTools::toString(sites->getNumberOfSequences()));
  ApplicationTools::displayResult("Number of sites", TextTools::toString(sites->getNumberOfSites()));
  
  // Get the initial tree
  Tree* tree = PhylogeneticsApplicationTools::getTree(bppancestor.getParams());
  ApplicationTools::displayResult("Number of leaves", TextTools::toString(tree->getNumberOfLeaves()));
  
  string treeWIdPath = ApplicationTools::getAFilePath("output.tree_ids.file", bppancestor.getParams(), false, false);
  if (treeWIdPath != "none")
  {
    TreeTemplate<Node> ttree(*tree);
    vector<Node *> nodes = ttree.getNodes();
    for(unsigned int i = 0; i < nodes.size(); i++)
    {
      if(nodes[i]->isLeaf())
        nodes[i]->setName(TextTools::toString(nodes[i]->getId()) + "_" + nodes[i]->getName());
      else
        nodes[i]->setBranchProperty("NodeId", BppString(TextTools::toString(nodes[i]->getId())));
    }
    Newick treeWriter;
    treeWriter.enableExtendedBootstrapProperty("NodeId");
    ApplicationTools::displayResult("Writing tagged tree to", treeWIdPath);
    treeWriter.write(ttree, treeWIdPath);
    delete tree;
    cout << "BppAncestor's done." << endl;
    exit(0);
  }

  DRTreeLikelihood *tl;
  string nhOpt = ApplicationTools::getStringParameter("nonhomogeneous", bppancestor.getParams(), "no", "", true, false);
  ApplicationTools::displayResult("Heterogeneous model", nhOpt);

  SubstitutionModel    *model    = 0;
  SubstitutionModelSet *modelSet = 0;
  DiscreteDistribution *rDist    = 0;
  unsigned int nbStates;

  if (nhOpt == "no")
  {  
    model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, sites, bppancestor.getParams());
    if (model->getName() != "RE08") SiteContainerTools::changeGapsToUnknownCharacters(*sites);
    if (model->getNumberOfStates() > model->getAlphabet()->getSize())
    {
      //Markov-modulated Markov model!
      rDist = new ConstantRateDistribution();
    }
    else
    {
      rDist = PhylogeneticsApplicationTools::getRateDistribution(bppancestor.getParams());
    }
    tl = new DRHomogeneousTreeLikelihood(*tree, *sites, model, rDist, true);
    nbStates = model->getNumberOfStates();
  }
  else if (nhOpt == "one_per_branch")
  {
    model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, sites, bppancestor.getParams());
    if (model->getName() != "RE08") SiteContainerTools::changeGapsToUnknownCharacters(*sites);
    if (model->getNumberOfStates() > model->getAlphabet()->getSize())
    {
      //Markov-modulated Markov model!
      rDist = new ConstantRateDistribution();
    }
    else
    {
      rDist = PhylogeneticsApplicationTools::getRateDistribution(bppancestor.getParams());
    }
    vector<double> rateFreqs;
    if (model->getNumberOfStates() != alphabet->getSize())
    {
      //Markov-Modulated Markov Model...
      unsigned int n =(unsigned int)(model->getNumberOfStates() / alphabet->getSize());
      rateFreqs = vector<double>(n, 1./(double)n); // Equal rates assumed for now, may be changed later (actually, in the most general case,
                                                   // we should assume a rate distribution for the root also!!!  
    }
    FrequenciesSet * rootFreqs = PhylogeneticsApplicationTools::getRootFrequenciesSet(alphabet, sites, bppancestor.getParams(), rateFreqs);
    vector<string> globalParameters = ApplicationTools::getVectorParameter<string>("nonhomogeneous_one_per_branch.shared_parameters", bppancestor.getParams(), ',', "");
    modelSet = SubstitutionModelSetTools::createNonHomogeneousModelSet(model, rootFreqs, tree, globalParameters); 
    model = 0;
    tl = new DRNonHomogeneousTreeLikelihood(*tree, *sites, modelSet, rDist, true);
    nbStates = modelSet->getNumberOfStates();
  }
  else if (nhOpt == "general")
  {
    modelSet = PhylogeneticsApplicationTools::getSubstitutionModelSet(alphabet, sites, bppancestor.getParams());
    if (modelSet->getModel(0)->getName() != "RE08") SiteContainerTools::changeGapsToUnknownCharacters(*sites);
    if (modelSet->getNumberOfStates() > modelSet->getAlphabet()->getSize())
    {
      //Markov-modulated Markov model!
      rDist = new ConstantDistribution(1.);
    }
    else
    {
      rDist = PhylogeneticsApplicationTools::getRateDistribution(bppancestor.getParams());
    }
    tl = new DRNonHomogeneousTreeLikelihood(*tree, *sites, modelSet, rDist, true);
    nbStates = modelSet->getNumberOfStates();
  }
  else throw Exception("Unknown option for nonhomogeneous: " + nhOpt);
  tl->initialize();
 
  delete tree;
    
  double logL = tl->getValue();
  if (isinf(logL))
  {
    // This may be due to null branch lengths, leading to null likelihood!
    ApplicationTools::displayWarning("!!! Warning!!! Likelihood is zero.");
    ApplicationTools::displayWarning("!!! This may be due to branch length == 0.");
    ApplicationTools::displayWarning("!!! All null branch lengths will be set to 0.000001.");
    ParameterList pl = tl->getBranchLengthsParameters();
    for(unsigned int i = 0; i < pl.size(); i++)
    {
      if(pl[i].getValue() < 0.000001) pl[i].setValue(0.000001);
    }
    tl->matchParametersValues(pl);
    logL = tl->getValue();
  }
  if (isinf(logL))
  {
    ApplicationTools::displayError("!!! Unexpected likelihood == 0.");
    ApplicationTools::displayError("!!! Looking at each site:");
    for(unsigned int i = 0; i < sites->getNumberOfSites(); i++)
    {
      (*ApplicationTools::error << "Site " << sites->getSite(i).getPosition() << "\tlog likelihood = " << tl->getLogLikelihoodForASite(i)).endLine();
    }
    ApplicationTools::displayError("!!! 0 values (inf in log) may be due to computer overflow, particularily if datasets are big (>~500 sequences).");
    exit(-1);
  }
  tree = new TreeTemplate<Node>(tl->getTree());

  // Write parameters to screen:
  ApplicationTools::displayResult("Log likelihood", TextTools::toString(tl->getValue(), 15));
  ParameterList parameters = tl->getSubstitutionModelParameters();
  for (unsigned int i = 0; i < parameters.size(); i++)
  {
    ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
  }
  parameters = tl->getRateDistributionParameters();
  for (unsigned int i = 0; i < parameters.size(); i++)
  {
    ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
  }

  // Getting posterior rate class distribution:
  DiscreteDistribution* prDist = RASTools::getPosteriorRateDistribution(*tl);
  ApplicationTools::displayMessage("\nPosterior rate distribution for dataset:\n");
  if (ApplicationTools::message) prDist->print(*ApplicationTools::message);
  ApplicationTools::displayMessage("\n");
  delete prDist;

  // Reconstruct ancestral sequences:
  string reconstruction = ApplicationTools::getStringParameter("asr.method", bppancestor.getParams(), "marginal", "", true, false);
  ApplicationTools::displayResult("Ancestral state reconstruction method", reconstruction);
  bool probs = false;

  AncestralStateReconstruction *asr = 0;
  bool probMethod = false;
  if (reconstruction == "none")
  {
    //do nothing
  } else if (reconstruction == "marginal") {
    asr = new MarginalAncestralStateReconstruction(tl);
    probMethod = true;
  } else
    throw Exception("Unknown ancestral state reconstruction method: " + reconstruction);

  string outputFile;
  if (asr) {
    if (probMethod)
    {
      probs = ApplicationTools::getBooleanParameter("asr.probabilities", bppancestor.getParams(), false, "", true, false);
      ApplicationTools::displayResult("Output probabilities", probs ? "yes" : "no");
    }

    // Write infos to file:
    outputFile = ApplicationTools::getAFilePath("output.sites.file", bppancestor.getParams(), false, false);
    if (outputFile != "none")
    {
      ApplicationTools::displayResult("Output file for sites", outputFile);
      ofstream out(outputFile.c_str(), ios::out);
      TreeTemplate<Node> ttree(*tree);
      vector<Node *> nodes = ttree.getInnerNodes();
      size_t nbNodes = nodes.size();
    
      // Get the rate class with maximum posterior probability:
      vector<size_t> classes = tl->getRateClassWithMaxPostProbOfEachSite();
      // Get the posterior rate, i.e. rate averaged over all posterior probabilities:
      Vdouble rates = tl->getPosteriorRateOfEachSite();
      // Get the ancestral sequences:
      vector<Sequence*> sequences(nbNodes);
      vector<VVdouble*> probabilities(nbNodes);

      vector<string> colNames;
      colNames.push_back("Sites");
      colNames.push_back("is.complete");
      colNames.push_back("is.constant");
      colNames.push_back("lnL");
      colNames.push_back("rc");
      colNames.push_back("pr");
      for (size_t i = 0; i < nbNodes; i++) {
        Node *node = nodes[i];
        colNames.push_back("max." + TextTools::toString(node->getId()));
        if (probs) {
          probabilities[i] = new VVdouble();
          //The cast will have to be updated when more probabilistic method will be available:
          sequences[i] = dynamic_cast<MarginalAncestralStateReconstruction *>(asr)->getAncestralSequenceForNode(node->getId(), probabilities[i], false);

          for (unsigned int j = 0; j < nbStates; j++) {
            colNames.push_back("prob." + TextTools::toString(node->getId()) + "." + alphabet->intToChar((int)j));
          }
        }
        else
        {
          if (node->isLeaf()) {

          } else {
            sequences[i] = asr->getAncestralSequenceForNode(node->getId());
          }
        }
      }

      //Now fill the table:
      vector<string> row(colNames.size());
      DataTable* infos = new DataTable(colNames);
    
      for (size_t i = 0; i < sites->getNumberOfSites(); i++)
      {
        double lnL = tl->getLogLikelihoodForASite(i);
        const Site* currentSite = &sites->getSite(i);
        int currentSitePosition = currentSite->getPosition();
        string isCompl = "NA";
        string isConst = "NA";
        try { isCompl = (SiteTools::isComplete(*currentSite) ? "1" : "0"); }
        catch(EmptySiteException& ex) {}
        try { isConst = (SiteTools::isConstant(*currentSite) ? "1" : "0"); }
        catch(EmptySiteException& ex) {}
        row[0] = (string("[" + TextTools::toString(currentSitePosition) + "]"));
        row[1] = isCompl;
        row[2] = isConst;
        row[3] = TextTools::toString(lnL);
        row[4] = TextTools::toString(classes[i]);
        row[5] = TextTools::toString(rates[i]);

        unsigned int k = 6;
        for (unsigned int j = 0; j < nbNodes; j++) {
          row[k] = sequences[j]->getChar(i);
          k++;
          if (probs) {
            for (unsigned int l = 0; l < nbStates; l++) {
              row[k] = TextTools::toString((*probabilities[j])[i][l]);
              k++;
            }
          }
        }

        infos->addRow(row);
      }

      DataTable::write(*infos, out, "\t");

      delete infos;
    }

    SiteContainer* asSites = 0;
    if (probMethod)
    {
      bool sample = ApplicationTools::getBooleanParameter("asr.sample", bppancestor.getParams(), false, "", true, false);
      ApplicationTools::displayResult("Sample from posterior distribution", sample ? "yes" : "no");
      if (sample)
      {
        unsigned int nbSamples = ApplicationTools::getParameter<unsigned int>("asr.sample.number", bppancestor.getParams(), 1, "", true, false);
        asSites = new AlignedSequenceContainer(alphabet);
        for (unsigned int i = 0; i < nbSamples; i++)
        {
          ApplicationTools::displayGauge(i, nbSamples-1, '=');
          SequenceContainer *sampleSites = dynamic_cast<MarginalAncestralStateReconstruction *>(asr)->getAncestralSequences(true);
          vector<string> names = sampleSites->getSequencesNames();
          for (unsigned int j = 0; j < names.size(); j++)
            names[j] += "_" + TextTools::toString(i+1);
          sampleSites->setSequencesNames(names, false);
          SequenceContainerTools::append(*asSites, *sampleSites);
          delete sampleSites;
        }
        ApplicationTools::message->endLine();
      }
      else
      {
        asSites = asr->getAncestralSequences();
      }
    }
    else
    {
      asSites = asr->getAncestralSequences();
    }
  
    //Add existing sequence to output?
    bool addExtant = ApplicationTools::getBooleanParameter("asr.add_extant", bppancestor.getParams(), false, "", true, false);
    if (addExtant) {
      SequenceContainerTools::append(*asSites, *sites);
    }

    //Write output:
    if (ApplicationTools::getStringParameter("output.sequence.file", bppancestor.getParams(), "none") != "none") {
      SequenceApplicationTools::writeAlignmentFile(*asSites, bppancestor.getParams());
    }
    delete asSites;

    delete asr;
  }

  outputFile = ApplicationTools::getAFilePath("output.nodes.file", bppancestor.getParams(), false, false);
  if (outputFile != "none")
  {
    ApplicationTools::displayResult("Output file for nodes", outputFile);
    ofstream out(outputFile.c_str(), ios::out);

    //Add existing sequence to output?
    bool addExtant = ApplicationTools::getBooleanParameter("output.nodes.add_extant", bppancestor.getParams(), false, "", true, false);
    
    map<int, vector<double> > frequencies;
    TreeLikelihoodTools::getAncestralFrequencies(*tl, frequencies, addExtant);
    
    vector<string> colNames;
    colNames.push_back("Nodes");
    for (unsigned int i = 0; i < tl->getNumberOfStates(); i++)
      colNames.push_back("exp" + alphabet->intToChar(i));
    for (unsigned int i = 0; i < tl->getNumberOfStates(); i++)
      colNames.push_back("eb" + alphabet->intToChar(i));

    //Now fill the table:
    vector<string> row(colNames.size());
    DataTable* infos = new DataTable(colNames);
    
    for (map<int, vector<double> >::iterator it = frequencies.begin(); it != frequencies.end(); it++)
    {
      row[0] = TextTools::toString(it->first);
      Vdouble ebFreqs = DRTreeLikelihoodTools::getPosteriorStateFrequencies(*tl, it->first);
      for (unsigned int i = 0; i < tl->getNumberOfStates(); i++)
      {
        row[i + 1] = TextTools::toString(it->second[i]);
      }
      for (unsigned int i = 0; i < tl->getNumberOfStates(); i++)
      {
        row[i + tl->getNumberOfStates() + 1] = TextTools::toString(ebFreqs[i]);
      }
      infos->addRow(row);
    }
    
    DataTable::write(*infos, out, "\t");

    delete infos;
  }


  delete alphabet;
  delete sites;
  if(model)    delete model;
  if(modelSet) delete modelSet;
  delete rDist;
  delete tl;
  delete tree;
  bppancestor.done();

  }
  catch (exception & e)
  {
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}

