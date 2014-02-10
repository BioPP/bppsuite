//
// File: bppML.cpp
// Created by: Julien Dutheil
// Created on: Dec Sat 03 16:41 2005
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
#include <limits>

using namespace std;

#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Numeric/Prob/ConstantDistribution.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Numeric/Matrix/MatrixTools.h>
#include <Bpp/Numeric/Hmm/FullHmmTransitionMatrix.h>
#include <Bpp/Numeric/VectorTools.h>
#include <Bpp/Numeric/AutoParameter.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/KeyvalTools.h>

// From SeqLib:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>

// From PhylLib:
#include <Bpp/Phyl/Tree/Tree.h>
#include <Bpp/Phyl/Likelihood.all>
#include <Bpp/Phyl/PatternTools.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Model.all>
#include <Bpp/Phyl/Model/Protein/CoalaCore.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>
#include <Bpp/Phyl/Model/FrequenciesSet/MvaFrequenciesSet.h>
#include <Bpp/Phyl/Io/Newick.h>

#include <Bpp/Phyl/NewLikelihood/SingleRecursiveTreeLikelihoodCalculation.h>
#include <Bpp/Phyl/NewLikelihood/MixturePhyloLikelihood.h>
#include <Bpp/Phyl/NewLikelihood/HmmPhyloLikelihood.h>
#include <Bpp/Phyl/NewLikelihood/AutoCorrelationPhyloLikelihood.h>
#include <Bpp/Phyl/NewLikelihood/SubstitutionProcessCollection.h>

using namespace bpp;

using namespace newlik;

/******************************************************************************/

void help()
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
  (*ApplicationTools::message << "bppml parameter1_name=parameter1_value parameter2_name=parameter2_value").endLine();
  (*ApplicationTools::message << "      ... param=option_file").endLine();
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "  Refer to the Bio++ Program Suite Manual for a list of available options.").endLine();
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
}

int main(int args, char** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*       Bio++ Maximum Likelihood Computation, version 2          *" << endl;
  cout << "*                                                                *" << endl;
  cout << "* Authors: J. Dutheil                       Last Modif. 10/02/14 *" << endl;
  cout << "*          B. Boussau                                            *" << endl;
  cout << "*          L. Guéguen                                            *" << endl;
  cout << "*          M. Groussin                                           *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  if (args == 1)
    {
      help();
      return 0;
    }

  try
    {
      BppApplication bppml(args, argv, "BppML");
      bppml.startTimer();

      ///// Alphabet
    
      Alphabet* alphabet = SequenceApplicationTools::getAlphabet(bppml.getParams(), "", false);
      auto_ptr<GeneticCode> gCode;
      CodonAlphabet* codonAlphabet = dynamic_cast<CodonAlphabet*>(alphabet);
      if (codonAlphabet) {
        string codeDesc = ApplicationTools::getStringParameter("genetic_code", bppml.getParams(), "Standard", "", true, true);
        gCode.reset(SequenceApplicationTools::getGeneticCode(codonAlphabet->getNucleicAlphabet(), codeDesc));
      }


      ////// Data
    
      VectorSiteContainer* allSites = SequenceApplicationTools::getSiteContainer(alphabet, bppml.getParams());

      VectorSiteContainer* sites = SequenceApplicationTools::getSitesToAnalyse(*allSites, bppml.getParams(), "", true, false);
      delete allSites;

      ApplicationTools::displayResult("Number of sequences", TextTools::toString(sites->getNumberOfSequences()));
      ApplicationTools::displayResult("Number of sites", TextTools::toString(sites->getNumberOfSites()));


      /////// Get the initial tree  
    
      vector<Tree*> vTree;
    
      string initTreeOpt = ApplicationTools::getStringParameter("init.tree", bppml.getParams(), "user", "", false, false);
      
      ApplicationTools::displayResult("Input tree", initTreeOpt);

      if (initTreeOpt == "user")
        {
          vTree = PhylogeneticsApplicationTools::getTrees(bppml.getParams());

          if (vTree.size()==0){
            ApplicationTools::displayError("!!! Empty file tree.");
            exit(-1);
          }
          
          for (size_t i=0; i<vTree.size(); i++)
            ApplicationTools::displayResult("Number of leaves", TextTools::toString(vTree[i]->getNumberOfLeaves()));
        }
      else if (initTreeOpt == "random")
        {
          vector<string> names = sites->getSequencesNames();
          Tree* tree = TreeTemplateTools::getRandomTree(names);
          tree->setBranchLengths(1.);
          vTree.push_back(tree);
        }
      else throw Exception("Unknown init tree method.");

      // Try to write the current tree to file. This will be overwritten by the optimized tree,
      // but allow to check file existence before running optimization!

      vector<const Tree*> vcTree;
      for (size_t i=0; i<vTree.size(); i++)
        vcTree.push_back(vTree[i]);
      
      PhylogeneticsApplicationTools::writeTrees(vcTree, bppml.getParams());
    
      bool computeLikelihood = ApplicationTools::getBooleanParameter("compute.likelihood", bppml.getParams(), true, "", false, false);
      if (!computeLikelihood)
        {
          delete alphabet;
          delete sites;
          for (size_t i=0; i<vTree.size(); i++)
            delete vTree[i];
          cout << "BppML's done. Bye." << endl;
          return 0;
        }

      ////////////
      // Setting branch lengths?
    
      string initBrLenMethod = ApplicationTools::getStringParameter("init.brlen.method", bppml.getParams(), "Input", "", true, false);
      string cmdName;
      map<string, string> cmdArgs;
      KeyvalTools::parseProcedure(initBrLenMethod, cmdName, cmdArgs);
      if (cmdName == "Input")
        {
          // Is the root has to be moved to the midpoint position along the branch that contains it ? If no, do nothing!
          string midPointRootBrLengths = ApplicationTools::getStringParameter("midPointRootBrLengths", cmdArgs, "no", "", true, false);
          if(midPointRootBrLengths == "yes")
            for (size_t i=0; i< vTree.size(); i++)
              TreeTools::constrainedMidPointRooting(*vTree[i]);
        }
      else if (cmdName == "Equal")
        {
          double value = ApplicationTools::getDoubleParameter("value", cmdArgs, 0.1, "", true, false);
          if (value <= 0)
            throw Exception("Value for branch length must be superior to 0");
          ApplicationTools::displayResult("Branch lengths set to", value);
          for (size_t i=0; i< vTree.size(); i++)
            vTree[i]->setBranchLengths(value);
        }
      else if (cmdName == "Clock")
        {
          for (size_t i=0; i< vTree.size(); i++)
            TreeTools::convertToClockTree(*vTree[i], vTree[i]->getRootId(), true);
        }
      else if (cmdName == "Grafen")
        {
          string grafenHeight = ApplicationTools::getStringParameter("height", cmdArgs, "input", "", true, false);
          double h;
          for (size_t i=0; i< vTree.size(); i++)
            {
              Tree* tree=vTree[i];
              if (grafenHeight == "input")
                {
                  h = TreeTools::getHeight(*tree, tree->getRootId());
                }
              else
                {
                  h = TextTools::toDouble(grafenHeight);
                  if (h <= 0) throw Exception("Height must be positive in Grafen's method.");
                }
              ApplicationTools::displayResult("Total height", TextTools::toString(h));

              double rho = ApplicationTools::getDoubleParameter("rho", cmdArgs, 1., "", true, false);
              ApplicationTools::displayResult("Grafen's rho", rho);
              TreeTools::computeBranchLengthsGrafen(*tree, rho);
              double nh = TreeTools::getHeight(*tree, tree->getRootId());
              tree->scaleTree(h / nh);
            }
        }
      else
        throw Exception("Method '" + initBrLenMethod + "' unknown for computing branch lengths.");
      ApplicationTools::displayResult("Branch lengths", cmdName);
  
      string treeWIdPath = ApplicationTools::getAFilePath("output.tree_ids.file", bppml.getParams(), false, false);
      if (treeWIdPath != "none")
        {
          Newick treeWriter;
          treeWriter.enableExtendedBootstrapProperty("NodeId");
          ApplicationTools::displayResult("Writing tagged tree to", treeWIdPath);

          for (size_t iv=0; iv< vTree.size(); iv++)
            {
              TreeTemplate<Node> ttree(*vTree[iv]);
              vector<Node*> nodes = ttree.getNodes();
              for (size_t i = 0; i < nodes.size(); i++)
                {
                  if (nodes[i]->isLeaf())
                    nodes[i]->setName(TextTools::toString(nodes[i]->getId()) + "_" + nodes[i]->getName());
                  else
                    nodes[i]->setBranchProperty("NodeId", BppString(TextTools::toString(nodes[i]->getId())));
                }
              treeWriter.write(ttree, treeWIdPath, iv==0);
              delete vTree[iv];
            }
          cout << "BppML's done." << endl;
          exit(0);
        }


      /////////////////
      // Computing stuff
    
      DiscreteRatesAcrossSitesTreeLikelihood* tl_old = 0;
      newlik::PhyloLikelihood* tl_new = 0;

      bool checkTree    = ApplicationTools::getBooleanParameter("input.tree.check_root", bppml.getParams(), true, "", true, false);
      bool optimizeTopo = ApplicationTools::getBooleanParameter("optimization.topology", bppml.getParams(), false, "", true, false);
      unsigned int nbBS = ApplicationTools::getParameter<unsigned int>("bootstrap.number", bppml.getParams(), 0, "", true, false);
      string collection = ApplicationTools::getStringParameter("collection", bppml.getParams(), "", "", true, false);
    
    
      SubstitutionModel*    model    = 0;
      SubstitutionModelSet* modelSet = 0;
      DiscreteDistribution* rDist    = 0;

      map<string, string> unparsedparams;
      /// Topology estimation
    
      if (optimizeTopo || nbBS > 0)
        {
          if (collection != "")
            throw Exception("Topology estimation in collections not supported yet, sorry.");
      
          string nhOpt = ApplicationTools::getStringParameter("nonhomogeneous", bppml.getParams(), "no", "", true, false);
          ApplicationTools::displayResult("Heterogeneous model", nhOpt);

          if (nhOpt != "no")
            throw Exception("Topology estimation with NH model not supported yet, sorry :(");
          model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, gCode.get(), sites, bppml.getParams(), unparsedparams);
          if (model->getName() != "RE08") SiteContainerTools::changeGapsToUnknownCharacters(*sites);
          if (model->getNumberOfStates() >= 2 * model->getAlphabet()->getSize())
            {
              // Markov-modulated Markov model!
              rDist = new ConstantRateDistribution();
            }
          else
            {
              rDist = PhylogeneticsApplicationTools::getRateDistribution(bppml.getParams());
            }
          if (dynamic_cast<MixedSubstitutionModel*>(model) == 0)
            tl_old = new NNIHomogeneousTreeLikelihood(*vTree[0], *sites, model, rDist, checkTree, true);
          else
            throw Exception("Topology estimation with Mixed model not supported yet, sorry :(");
        }

      /// Constant Topology

    
      else
        {
          string recursion = ApplicationTools::getStringParameter("likelihood.recursion", bppml.getParams(), "simple", "", true, false);
          ApplicationTools::displayResult("Likelihood recursion", recursion);

          // Double recursion
          
          if (recursion == "double")
            {
              if (collection != "")
                throw Exception("Double recursion in collections not supported yet, sorry.");
      
              string nhOpt = ApplicationTools::getStringParameter("nonhomogeneous", bppml.getParams(), "no", "", true, false);
              ApplicationTools::displayResult("Heterogeneous model", nhOpt);


              if (nhOpt == "no")
                {
                  model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, gCode.get(), sites, bppml.getParams(), unparsedparams);
                  if (model->getName() != "RE08") SiteContainerTools::changeGapsToUnknownCharacters(*sites);
                  if (model->getNumberOfStates() >= 2 * model->getAlphabet()->getSize())
                    {
                      // Markov-modulated Markov model!
                      rDist = new ConstantRateDistribution();
                    }
                  else
                    rDist = PhylogeneticsApplicationTools::getRateDistribution(bppml.getParams());

                  if (dynamic_cast<MixedSubstitutionModel*>(model))
                    tl_old = new DRHomogeneousMixedTreeLikelihood(*vTree[0], *sites, model, rDist, checkTree);
                  else
                    tl_old = new DRHomogeneousTreeLikelihood(*vTree[0], *sites, model, rDist, checkTree);
                }
        
              else if (nhOpt == "one_per_branch")
                {
                  model = PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, gCode.get(), sites, bppml.getParams(), unparsedparams);
                  if (model->getName() != "RE08") SiteContainerTools::changeGapsToUnknownCharacters(*sites);
                  if (model->getNumberOfStates() >= 2 * model->getAlphabet()->getSize())
                    {
                      // Markov-modulated Markov model!
                      rDist = new ConstantRateDistribution();
                    }
                  else
                    rDist = PhylogeneticsApplicationTools::getRateDistribution(bppml.getParams());

                  vector<double> rateFreqs;
                  if (model->getNumberOfStates() != alphabet->getSize())
                    {
                      // Markov-Modulated Markov Model...
                      unsigned int n = (unsigned int)(model->getNumberOfStates() / alphabet->getSize());
                      rateFreqs = vector<double>(n, 1. / static_cast<double>(n)); // Equal rates assumed for now, may be changed later (actually, in the most general case,
                      // we should assume a rate distribution for the root also!!!
                    }

                  bool stationarity = ApplicationTools::getBooleanParameter("nonhomogeneous.stationarity", bppml.getParams(), false, "", false, false);
                  FrequenciesSet* rootFreqs = 0;
                  if (!stationarity)
                    {
                      rootFreqs = PhylogeneticsApplicationTools::getRootFrequenciesSet(alphabet, gCode.get(), sites, bppml.getParams(), rateFreqs);
                      stationarity = !rootFreqs;
                      string freqDescription = ApplicationTools::getStringParameter("nonhomogeneous.root_freq", bppml.getParams(), "");
                      if (freqDescription == "MVAprotein")
                        {
                          if (dynamic_cast<CoalaCore*>(model))
                            {
                              dynamic_cast<MvaFrequenciesSet*>(rootFreqs)->setModelName("MVAprotein");
                              dynamic_cast<MvaFrequenciesSet*>(rootFreqs)->initSet(dynamic_cast<CoalaCore*>(model)); 
                            }
                          else
                            throw Exception("The MVAprotein frequencies set at the root can only be used if a COaLA model is used on branches.");
                        }
                    }
                  ApplicationTools::displayBooleanResult("Stationarity assumed", stationarity);
   
                  vector<string> globalParameters = ApplicationTools::getVectorParameter<string>("nonhomogeneous_one_per_branch.shared_parameters", bppml.getParams(), ',', "");
                  for (unsigned int i = 0; i < globalParameters.size(); i++)
                    ApplicationTools::displayResult("Global parameter", globalParameters[i]);
                  modelSet = SubstitutionModelSetTools::createNonHomogeneousModelSet(model, rootFreqs, vTree[0], globalParameters);
                  model = 0;

                  if (dynamic_cast<MixedSubstitutionModelSet*>(modelSet))
                    throw Exception("Double recursion with non homogeneous mixed models is not implemented yet.");
                  else
                    tl_old = new DRNonHomogeneousTreeLikelihood(*vTree[0], *sites, modelSet, rDist, true);
                }
              else if (nhOpt == "general")
                {
                  modelSet = PhylogeneticsApplicationTools::getSubstitutionModelSet(alphabet, gCode.get(), sites, bppml.getParams());
                  if (modelSet->getModel(0)->getName() != "RE08") SiteContainerTools::changeGapsToUnknownCharacters(*sites);
                  if (modelSet->getNumberOfStates() >= 2 * modelSet->getAlphabet()->getSize())
                    {
                      // Markov-modulated Markov model!
                      rDist = new ConstantRateDistribution();
                    }
                  else
                    rDist = PhylogeneticsApplicationTools::getRateDistribution(bppml.getParams());

                  if (dynamic_cast<MixedSubstitutionModelSet*>(modelSet))
                    throw Exception("Double recursion with non homogeneous mixed models is not implemented yet.");
                  else
                    tl_old = new DRNonHomogeneousTreeLikelihood(*vTree[0], *sites, modelSet, rDist, true);
                }
              else throw Exception("Unknown option for nonhomogeneous: " + nhOpt);

              tl_old->initialize();

              for (size_t i=0; i< vTree.size(); i++)
                delete vTree[i];
    
            }

          // Simple recursion
      
          else if (recursion=="simple")
            {
              SiteContainerTools::changeGapsToUnknownCharacters(*sites);

              ///  Collection
              if (collection!="")
                {
                  string collName;
                  map<string, string> collArgs;
                  KeyvalTools::parseProcedure(collection, collName, collArgs);

                  if (recursion!="simple")
                    throw Exception("Double recursion is not implemented in collections yet, sorry.");

                  SubstitutionProcessCollection* SPC=PhylogeneticsApplicationTools::getSubstitutionProcessCollection(alphabet, gCode.get(), sites, vTree, bppml.getParams());

                  
                  ApplicationTools::displayResult("Collection type", collName);

                  if (collName=="Mixture")
                  {
                    MixturePhyloLikelihood* pMP = new MixturePhyloLikelihood(*sites, SPC);
                    tl_new = pMP;

                    size_t nbP = pMP->getCollection()->getNumberOfSubstitutionProcess();
                    
                    std::vector<double> vprob=ApplicationTools::getVectorParameter<double>("probas", collArgs, ',', "("+VectorTools::paste(std::vector<double>(nbP,1./(double)nbP))+")");
                    if (vprob.size()!=1)
                    {
                      if (vprob.size()!=pMP->getCollection()->getNumberOfSubstitutionProcess())
                        throw BadSizeException("Wrong size of probas description in Mixture", vprob.size(), pMP->getCollection()->getNumberOfSubstitutionProcess());
                      Simplex si(vprob);
                      pMP->setSubProcessProb(si);
                    }
                  }
                  else if (collName=="HMM")
                  {
                    HmmPhyloLikelihood* pMP = new HmmPhyloLikelihood(*sites, SPC);
                    
                    tl_new = pMP;
                    size_t nbP = pMP->getCollection()->getNumberOfSubstitutionProcess();

                    string vs="("+VectorTools::paste(std::vector<double>(nbP,1./(double)nbP),",")+")";
                    string vvs="(";
                    for (size_t i=0;i<nbP;i++)
                      vvs+=(i==0?"":",")+vs;
                    vvs+=")";

                    RowMatrix<double> mat=ApplicationTools::getMatrixParameter<double>("probas", collArgs, ',', vvs);
                    
                    FullHmmTransitionMatrix fhtm(pMP->getHmmTransitionMatrix().getHmmStateAlphabet(), pMP->getNamespace());
                    fhtm.setTransitionProbabilities(mat);
                    
                    pMP->matchParametersValues(fhtm.getParameters());
                    
                  }
                  else if (collName=="AutoCorr")
                  {
                    AutoCorrelationPhyloLikelihood* pMP = new AutoCorrelationPhyloLikelihood(*sites, SPC);
                    
                    tl_new = pMP;
                    size_t nbP = pMP->getCollection()->getNumberOfSubstitutionProcess();

                    string vs="("+VectorTools::paste(std::vector<double>(nbP,1./(double)nbP),",")+")";

                    vector<double> v=ApplicationTools::getVectorParameter<double>("probas", collArgs, ',', vs);

                    ParameterList pl;

                    for (size_t i=0;i<v.size();i++)
                      pl.addParameter(Parameter("lambda"+TextTools::toString(i+1),v[i]));
                    
                        
                    pMP->matchParametersValues(pl);
                    
                  }
                  else
                    throw Exception("Unknown Multiple Phylogeny description : "+ collName);
                  
                }
              else
                {
                  string nhOpt = ApplicationTools::getStringParameter("nonhomogeneous", bppml.getParams(), "no", "", true, false);
                  ApplicationTools::displayResult("Heterogeneous model", nhOpt);

                  auto_ptr<SubstitutionProcess> process(PhylogeneticsApplicationTools::getSubstitutionProcess(alphabet, gCode.get(), sites, vTree, bppml.getParams()));
        
                  if (process->getSubstitutionModel(0,0).getName() != "RE08")
                    SiteContainerTools::changeGapsToUnknownCharacters(*sites);

                  string compression = ApplicationTools::getStringParameter("likelihood.recursion_simple.compression", bppml.getParams(), "recursive", "", true, false);
                  ApplicationTools::displayResult("Likelihood data compression", compression);

                  SingleRecursiveTreeLikelihoodCalculation* tlcomp;
                  
                  if (compression == "simple")
                  {
                    if (dynamic_cast<const MixedSubstitutionModel*>(&process->getSubstitutionModel(0,0)))
                      throw Exception("Simple recursion process with mixed models is not implemented yet.");
                    else 
                      tlcomp = new SingleRecursiveTreeLikelihoodCalculation(*sites, process.get(), true, compression=="simple");
                  }
                  else
                  {
                    if (dynamic_cast<const MixedSubstitutionModel*>(&process->getSubstitutionModel(0,0)))
                      throw Exception("Recursive recursion process with mixed models is not implemented yet.");
                    else 
                      tlcomp = new SingleRecursiveTreeLikelihoodCalculation(*sites, process.get(), true, compression=="recursive");
                  }
                  tl_new = new SinglePhyloLikelihood(process.release(), tlcomp);
                }
            }
          else
            throw Exception("Unknown recursion option: " + recursion);

          
          //Listing parameters
          string paramNameFile = ApplicationTools::getAFilePath("output.parameter_names.file", bppml.getParams(), false, false);
          if (paramNameFile != "none") {
            ApplicationTools::displayResult("List parameters to", paramNameFile);
            ofstream pnfile(paramNameFile.c_str(), ios::out);
      
            ParameterList pl;
            if (tl_old)
              pl=tl_old->getParameters();
            else
              pl=tl_new->getParameters();

            for (unsigned int i = 0; i < pl.size(); ++i) {
              pnfile << pl[i].getName() << endl;
            }
            pnfile.close();
            cout << "BppML's done." << endl;
            exit(0);
          }


          // Old optimization
    
          if (tl_old){
            //Check initial likelihood:
            double logL = tl_old->getValue();
            if (isinf(logL))
              {
                // This may be due to null branch lengths, leading to null likelihood!
                ApplicationTools::displayWarning("!!! Warning!!! Initial likelihood is zero.");
                ApplicationTools::displayWarning("!!! This may be due to branch length == 0.");
                ApplicationTools::displayWarning("!!! All null branch lengths will be set to 0.000001.");
                ParameterList pl = tl_old->getBranchLengthsParameters();
                for (unsigned int i = 0; i < pl.size(); i++)
                  {
                    if (pl[i].getValue() < 0.000001) pl[i].setValue(0.000001);
                  }
                tl_old->matchParametersValues(pl);
                logL = tl_old->getValue();
              }
            ApplicationTools::displayResult("Initial log likelihood", TextTools::toString(-logL, 15));
            if (isinf(logL))
              {
                ApplicationTools::displayError("!!! Unexpected initial likelihood == 0.");
                if (codonAlphabet)
                  {
                    bool f = false;
                    size_t s;
                    for (size_t i = 0; i < sites->getNumberOfSites(); i++) {
                      if (isinf(tl_old->getLogLikelihoodForASite(i))) {
                        const Site& site = sites->getSite(i);
                        s = site.size();
                        for (size_t j = 0; j < s; j++) {
                          if (gCode->isStop(site.getValue(j))) {
                            (*ApplicationTools::error << "Stop Codon at site " << site.getPosition() << " in sequence " << sites->getSequence(j).getName()).endLine();
                            f = true;
                          }
                        }
                      }
                    }
                    if (f)
                      exit(-1);
                  }
                bool removeSaturated = ApplicationTools::getBooleanParameter("input.sequence.remove_saturated_sites", bppml.getParams(), false, "", true, false);
                if (!removeSaturated) {
                  ApplicationTools::displayError("!!! Looking at each site:");
                  for (unsigned int i = 0; i < sites->getNumberOfSites(); i++) {
                    (*ApplicationTools::error << "Site " << sites->getSite(i).getPosition() << "\tlog likelihood = " << tl_old->getLogLikelihoodForASite(i)).endLine();
                  }
                  ApplicationTools::displayError("!!! 0 values (inf in log) may be due to computer overflow, particularily if datasets are big (>~500 sequences).");
                  ApplicationTools::displayError("!!! You may want to try input.sequence.remove_saturated_sites = yes to ignore positions with likelihood 0.");
                  exit(1);
                } else {
                  ApplicationTools::displayBooleanResult("Saturated site removal enabled", true);
                  for (size_t i = sites->getNumberOfSites(); i > 0; --i) {
                    if (isinf(tl_old->getLogLikelihoodForASite(i - 1))) {
                      ApplicationTools::displayResult("Ignore saturated site", sites->getSite(i - 1).getPosition());
                      sites->deleteSite(i - 1);
                    }
                  }
                  ApplicationTools::displayResult("Number of sites retained", sites->getNumberOfSites());
                  tl_old->setData(*sites);
                  tl_old->initialize();
                  logL = tl_old->getValue();
                  if (isinf(logL)) {
                    ApplicationTools::displayError("This should not happen. Exiting now.");
                    exit(1);
                  }
                  ApplicationTools::displayResult("Initial log likelihood", TextTools::toString(-logL, 15));
                }
              }
      
            tl_old = dynamic_cast<DiscreteRatesAcrossSitesTreeLikelihood*>(
                                                                           PhylogeneticsApplicationTools::optimizeParameters(tl_old, tl_old->getParameters(), bppml.getParams()));

            Tree* tree = new TreeTemplate<Node>(tl_old->getTree());
            PhylogeneticsApplicationTools::writeTree(*tree, bppml.getParams());

            // Write parameters to screen:
            ApplicationTools::displayResult("Log likelihood", TextTools::toString(-tl_old->getValue(), 15));
            ParameterList parameters = tl_old->getSubstitutionModelParameters();
            for (unsigned int i = 0; i < parameters.size(); i++)
              {
                ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
              }
            parameters = tl_old->getRateDistributionParameters();
            for (unsigned int i = 0; i < parameters.size(); i++)
              {
                ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
              }
      
            // Checking convergence:
            PhylogeneticsApplicationTools::checkEstimatedParameters(tl_old->getParameters());
      
            // Write parameters to file:
            string parametersFile = ApplicationTools::getAFilePath("output.estimates", bppml.getParams(), false, false);
            ApplicationTools::displayResult("Output estimates to file", parametersFile);
            if (parametersFile != "none")
              {
                StlOutputStream out(new ofstream(parametersFile.c_str(), ios::out));
                out << "# Log likelihood = ";
                out.setPrecision(20) << (-tl_old->getValue());
                out.endLine();
                out << "# Number of sites = ";
                out.setPrecision(20) << sites->getNumberOfSites();
                out.endLine();
                out.endLine();
                out << "# Substitution model parameters:";
                out.endLine();
                if (modelSet)
                  {
                    modelSet->matchParametersValues(tl_old->getParameters());
                    PhylogeneticsApplicationTools::printParameters(modelSet, out);
                  }
                else
                  {
                    model->matchParametersValues(tl_old->getParameters());
                    PhylogeneticsApplicationTools::printParameters(model, out);
                  }
                out.endLine();
                (out << "# Rate distribution parameters:").endLine();
                rDist->matchParametersValues(tl_old->getParameters());
                PhylogeneticsApplicationTools::printParameters(rDist, out);
              }
      
            // Getting posterior rate class distribution:
            DiscreteDistribution* prDist = RASTools::getPosteriorRateDistribution(*tl_old);
            ApplicationTools::displayMessage("\nPosterior rate distribution for dataset:\n");
            if (ApplicationTools::message) prDist->print(*ApplicationTools::message);
            ApplicationTools::displayMessage("\n");
            delete prDist;
      
            // Write infos to file:
            string infosFile = ApplicationTools::getAFilePath("output.infos", bppml.getParams(), false, false);
            if (infosFile != "none")
              {
                ApplicationTools::displayResult("Alignment information logfile", infosFile);
                ofstream out(infosFile.c_str(), ios::out);
          
                // Get the rate class with maximum posterior probability:
                vector<size_t> classes = tl_old->getRateClassWithMaxPostProbOfEachSite();

                // Get the posterior rate, i.e. rate averaged over all posterior probabilities:
                Vdouble rates = tl_old->getPosteriorRateOfEachSite();
          
                vector<string> colNames;
                colNames.push_back("Sites");
                colNames.push_back("is.complete");
                colNames.push_back("is.constant");
                colNames.push_back("lnL");
                colNames.push_back("rc");
                colNames.push_back("pr");
                vector<string> row(6);
                DataTable* infos = new DataTable(colNames);
          
                for (unsigned int i = 0; i < sites->getNumberOfSites(); i++)
                  {
                    double lnL = tl_old->getLogLikelihoodForASite(i);
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
                    infos->addRow(row);
                  }
          
                DataTable::write(*infos, out, "\t");
          
                delete infos;
              }
          }
          // New optimization
          else // tl_new!=0 
            {
              //Check initial likelihood:
              double logL = tl_new->getValue();
              if (isinf(logL))
                {
                  // This may be due to null branch lengths, leading to null likelihood!
                  ApplicationTools::displayWarning("!!! Warning!!! Initial likelihood is zero.");
                  ApplicationTools::displayWarning("!!! This may be due to branch length == 0.");
                  ApplicationTools::displayWarning("!!! All null branch lengths will be set to 0.000001.");
                  ParameterList pl = tl_new->getBranchLengthsParameters();
                  for (unsigned int i = 0; i < pl.size(); i++)
                    {
                      if (pl[i].getValue() < 0.000001) pl[i].setValue(0.000001);
                    }
                  tl_new->matchParametersValues(pl);
                  logL = tl_new->getValue();
                }
              ApplicationTools::displayResult("Initial log likelihood", TextTools::toString(-logL, 15));
              if (isinf(logL))
                {
                  ApplicationTools::displayError("!!! Unexpected initial likelihood == 0.");
                  if (codonAlphabet)
                    {
                      bool f = false;
                      size_t s;
                      for (size_t i = 0; i < sites->getNumberOfSites(); i++) {
                        if (isinf(tl_new->getLogLikelihoodForASite(i))) {
                          const Site& site = sites->getSite(i);
                          s = site.size();
                          for (size_t j = 0; j < s; j++) {
                            if (gCode->isStop(site.getValue(j))) {
                              (*ApplicationTools::error << "Stop Codon at site " << site.getPosition() << " in sequence " << sites->getSequence(j).getName()).endLine();
                              f = true;
                            }
                          }
                        }
                      }
                      if (f)
                        exit(-1);
                    }
                  bool removeSaturated = ApplicationTools::getBooleanParameter("input.sequence.remove_saturated_sites", bppml.getParams(), false, "", true, false);
                  if (!removeSaturated) {
                    ApplicationTools::displayError("!!! Looking at each site:");
                    for (unsigned int i = 0; i < sites->getNumberOfSites(); i++) {
                      (*ApplicationTools::error << "Site " << sites->getSite(i).getPosition() << "\tlog likelihood = " << tl_new->getLogLikelihoodForASite(i)).endLine();
                    }
                    ApplicationTools::displayError("!!! 0 values (inf in log) may be due to computer overflow, particularily if datasets are big (>~500 sequences).");
                    ApplicationTools::displayError("!!! You may want to try input.sequence.remove_saturated_sites = yes to ignore positions with likelihood 0.");
                    exit(1);
                  } else {
                    ApplicationTools::displayBooleanResult("Saturated site removal enabled", true);
                    for (size_t i = sites->getNumberOfSites(); i > 0; --i) {
                      if (isinf(tl_new->getLogLikelihoodForASite(i - 1))) {
                        ApplicationTools::displayResult("Ignore saturated site", sites->getSite(i - 1).getPosition());
                        sites->deleteSite(i - 1);
                      }
                    }
                    ApplicationTools::displayResult("Number of sites retained", sites->getNumberOfSites());

                    tl_new->setData(*sites);
                    logL = tl_new->getValue();
                    if (isinf(logL)) {
                      ApplicationTools::displayError("This should not happen. Exiting now.");
                      exit(1);
                    }
                    ApplicationTools::displayResult("Initial log likelihood", TextTools::toString(-logL, 15));
                  }
                }
      
              tl_new = PhylogeneticsApplicationTools::optimizeParameters(tl_new, tl_new->getParameters(), bppml.getParams());

              if (dynamic_cast<SinglePhyloLikelihood*>(tl_new)!=NULL){
                Tree* tree = new TreeTemplate<Node>((dynamic_cast<SinglePhyloLikelihood*>(tl_new))->getTree());
                PhylogeneticsApplicationTools::writeTree(*tree, bppml.getParams());
              }
              else {
                std::vector<const TreeTemplate<Node>* > vTNree = dynamic_cast<MultiPhyloLikelihood*>(tl_new)->getTrees();
                PhylogeneticsApplicationTools::writeTrees(vTNree, bppml.getParams());
              }
              
              // Write parameters to screen:
              ApplicationTools::displayResult("Log likelihood", TextTools::toString(-tl_new->getValue(), 15));
              ParameterList parameters = tl_new->getParameters();
              parameters.deleteParameters(tl_new->getBranchLengthsParameters().getParameterNames(),false);
              
              for (unsigned int i = 0; i < parameters.size(); i++)
                ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
      
              // Checking convergence:
              PhylogeneticsApplicationTools::checkEstimatedParameters(tl_new->getParameters());
      
              // Write parameters to file:
              string parametersFile = ApplicationTools::getAFilePath("output.estimates", bppml.getParams(), false, false);
              ApplicationTools::displayResult("Process estimates to file", parametersFile);
              if (parametersFile != "none")
                {
                  StlOutputStream out(new ofstream(parametersFile.c_str(), ios::out));
                  PhylogeneticsApplicationTools::printParameters(tl_new, out);
                }

              // Write infos to file:
              //     probabilities of rate discrete distributions
              //     site infos : lnL, class (or process in case of collection) posterior probability distribution

              string infosFile = ApplicationTools::getAFilePath("output.infos", bppml.getParams(), false, false);
              if (infosFile != "none")
              {
                ApplicationTools::displayResult("Alignment information logfile", infosFile);
                StlOutputStream out(new ofstream(infosFile.c_str(), ios::out));
                  
                PhylogeneticsApplicationTools::printAnalysisInformation(tl_new, out);
              }
              
              ////////////////////////////////////////////
              // Bootstrap:
    
              string optimizeClock = ApplicationTools::getStringParameter("optimization.clock", bppml.getParams(), "None", "", true, false);
              if (nbBS > 0 && optimizeClock != "None")
                {
                  ApplicationTools::displayError("Bootstrap is not supported with clock trees.");
                }
              if (nbBS > 0 && optimizeClock == "None")
                {
                  ApplicationTools::displayResult("Number of bootstrap samples", TextTools::toString(nbBS));
                  bool approx = ApplicationTools::getBooleanParameter("bootstrap.approximate", bppml.getParams(), true);
                  ApplicationTools::displayResult("Use approximate bootstrap", TextTools::toString(approx ? "yes" : "no"));
                  bool bootstrapVerbose = ApplicationTools::getBooleanParameter("bootstrap.verbose", bppml.getParams(), false, "", true, false);

                  const Tree* initTree = vTree[0];
                  if (!bootstrapVerbose) bppml.getParam("optimization.verbose") = "0";
                  bppml.getParam("optimization.profiler") = "none";
                  bppml.getParam("optimization.messageHandler") = "none";
                  if (!optimizeTopo)
                    {
                      bppml.getParam("optimization.topology") = "yes";
                      tl_old = dynamic_cast<NNIHomogeneousTreeLikelihood*>(
                                                                           PhylogeneticsApplicationTools::optimizeParameters(tl_old, tl_old->getParameters(), bppml.getParams(), "", true, false));
                      initTree = &tl_old->getTree();
                    }

                  string bsTreesPath = ApplicationTools::getAFilePath("bootstrap.output.file", bppml.getParams(), false, false);
                  ofstream* out = 0;
                  if (bsTreesPath != "none")
                    {
                      ApplicationTools::displayResult("Bootstrap trees stored in file", bsTreesPath);
                      out = new ofstream(bsTreesPath.c_str(), ios::out);
                    }
                  Newick newick;
                  ParameterList paramsToIgnore = tl_old->getSubstitutionModelParameters();
                  paramsToIgnore.addParameters(tl_old->getRateDistributionParameters());

                  ApplicationTools::displayTask("Bootstrapping", true);
                  vector<Tree*> bsTrees(nbBS);
                  for (unsigned int i = 0; i < nbBS; i++)
                    {
                      ApplicationTools::displayGauge(i, nbBS - 1, '=');
                      VectorSiteContainer* sample = SiteContainerTools::bootstrapSites(*sites);
                      if (!approx)
                        {
                          model->setFreqFromData(*sample);
                        }

                      if (dynamic_cast<MixedSubstitutionModel*>(model) != NULL)
                        throw Exception("Bootstrap estimation with Mixed model not supported yet, sorry :(");

                      NNIHomogeneousTreeLikelihood* tlRep = new NNIHomogeneousTreeLikelihood(*initTree, *sample, model, rDist, true, false);
                      tlRep->initialize();
                      ParameterList parametersRep = tlRep->getParameters();
                      if (approx)
                        {
                          parametersRep.deleteParameters(paramsToIgnore.getParameterNames());
                        }
                      tlRep = dynamic_cast<NNIHomogeneousTreeLikelihood*>(
                                                                          PhylogeneticsApplicationTools::optimizeParameters(tlRep, parametersRep, bppml.getParams(), "", true, false));
                      bsTrees[i] = new TreeTemplate<Node>(tlRep->getTree());
                      if (out && i == 0) newick.write(*bsTrees[i], bsTreesPath, true);
                      if (out && i >  0) newick.write(*bsTrees[i], bsTreesPath, false);
                      delete tlRep;
                      delete sample;
                    }
                  if (out) out->close();
                  if (out) delete out;
                  ApplicationTools::displayTaskDone();


                  ApplicationTools::displayTask("Compute bootstrap values");
                  TreeTools::computeBootstrapValues(*vTree[0], bsTrees);
                  ApplicationTools::displayTaskDone();
                  for (unsigned int i = 0; i < nbBS; i++)
                    {
                      delete bsTrees[i];
                    }

                  // Write resulting tree:
                  vcTree.clear();
                  for (size_t i=0; i<vTree.size(); i++)
                    vcTree.push_back(vTree[i]);
      
                  PhylogeneticsApplicationTools::writeTrees(vcTree, bppml.getParams());
                }


              delete alphabet;
              delete sites;
              if (model) delete model;
              if (modelSet) delete modelSet;
              delete rDist;
              delete tl_old;
              for (size_t i=0; i< vTree.size(); i++)
                delete vTree[i];
              bppml.done();
            }
        }
    }
  catch (exception& e)
    {
      cout << e.what() << endl;
      return 1;
    }
  
  return 0;
}

