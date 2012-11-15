//
// File: bppMixedLikelihoods.cpp
// Created by: Laurent Guéguen
// Created on: lundi 12 novembre 2012, à 07h 02
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
#include <Bpp/Phyl/Model.all>
#include <Bpp/Phyl/Io/Newick.h>

using namespace bpp;

/******************************************************************************/

void help()
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
  (*ApplicationTools::message << "bppmixedlikelihoods parameter1_name=parameter1_value ").endLine();
  (*ApplicationTools::message << "      parameter2_name=parameter2_value ... param=option_file").endLine();
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "  Refer to the Bio++ Program Suite Manual for a list of available options.").endLine();
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
}

int main(int args, char ** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*     Bio++ Computation of site likelihoods inside mixed models  *" << endl;
  cout << "* Author: L. Guéguen                        Created on: 12/11/12 *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  if (args == 1)
  {
    help();
    return 0;
  }
  
  try {
    
    BppApplication bppmixedlikelihoods(args, argv, "BppMixedLikelihoods");
    bppmixedlikelihoods.startTimer();

    Alphabet* alphabet = SequenceApplicationTools::getAlphabet(bppmixedlikelihoods.getParams(), "", false);

    // get the data
    
    VectorSiteContainer* allSites = SequenceApplicationTools::getSiteContainer(alphabet, bppmixedlikelihoods.getParams());
    
    VectorSiteContainer* sites = SequenceApplicationTools::getSitesToAnalyse(* allSites, bppmixedlikelihoods.getParams(), "", true, false);
    delete allSites;
    
    ApplicationTools::displayResult("Number of sequences", TextTools::toString(sites->getNumberOfSequences()));
    ApplicationTools::displayResult("Number of sites", TextTools::toString(sites->getNumberOfSites()));
    
    // Get the tree
    Tree* tree = PhylogeneticsApplicationTools::getTree(bppmixedlikelihoods.getParams());
    ApplicationTools::displayResult("Number of leaves", TextTools::toString(tree->getNumberOfLeaves()));
    

    
    AbstractDiscreteRatesAcrossSitesTreeLikelihood *tl;
    string nhOpt = ApplicationTools::getStringParameter("nonhomogeneous", bppmixedlikelihoods.getParams(), "no", "", true, false);
    ApplicationTools::displayResult("Heterogeneous model", nhOpt);
    
    MixedSubstitutionModel    *model    = 0;
    MixedSubstitutionModelSet *modelSet = 0;
    DiscreteDistribution *rDist    = 0;
    unsigned int nbStates;

    if (nhOpt == "no")
      {  
        model = dynamic_cast<MixedSubstitutionModel*>(PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, sites, bppmixedlikelihoods.getParams()));
        if (model==0){
          cout << "Model is not a Mixed model" << endl;
          exit(0);
        }
        
        SiteContainerTools::changeGapsToUnknownCharacters(*sites);
        if(model->getNumberOfStates() > model->getAlphabet()->getSize())
          {
            //Markov-modulated Markov model!
            rDist = new ConstantDistribution(1., true);
          }
        else
          {
            rDist = PhylogeneticsApplicationTools::getRateDistribution(bppmixedlikelihoods.getParams());
          }
        tl = new RHomogeneousMixedTreeLikelihood(*tree, *sites, model, rDist, true);
        nbStates = model->getNumberOfStates();
      }
    else if (nhOpt == "one_per_branch")
      {
        model = dynamic_cast<MixedSubstitutionModel*>(PhylogeneticsApplicationTools::getSubstitutionModel(alphabet, sites, bppmixedlikelihoods.getParams()));
        if (model==0){
          cout << "Model is not a Mixed model" << endl;
          exit(0);
        }
        
        SiteContainerTools::changeGapsToUnknownCharacters(*sites);
        if (model->getNumberOfStates() > model->getAlphabet()->getSize())
          {
            //Markov-modulated Markov model!
            rDist = new ConstantDistribution(1., true);
          }
        else
          {
            rDist = PhylogeneticsApplicationTools::getRateDistribution(bppmixedlikelihoods.getParams());
          }
        vector<double> rateFreqs;
        if (model->getNumberOfStates() != alphabet->getSize())
          {
            //Markov-Modulated Markov Model...
            unsigned int n =(unsigned int)(model->getNumberOfStates() / alphabet->getSize());
            rateFreqs = vector<double>(n, 1./(double)n); // Equal rates assumed for now, may be changed later (actually, in the most general case,
            // we should assume a rate distribution for the root also!!!  
          }
        FrequenciesSet * rootFreqs = PhylogeneticsApplicationTools::getRootFrequenciesSet(alphabet, sites, bppmixedlikelihoods.getParams(), rateFreqs);
        vector<string> globalParameters = ApplicationTools::getVectorParameter<string>("nonhomogeneous_one_per_branch.shared_parameters", bppmixedlikelihoods.getParams(), ',', "");
        modelSet = dynamic_cast<MixedSubstitutionModelSet*>(SubstitutionModelSetTools::createNonHomogeneousModelSet(model, rootFreqs, tree, globalParameters));
        model = 0;
        tl = new RNonHomogeneousMixedTreeLikelihood(*tree, *sites, modelSet, rDist, true);
        nbStates = modelSet->getNumberOfStates();
      }
    else if (nhOpt == "general")
      {
        modelSet = dynamic_cast<MixedSubstitutionModelSet*>(PhylogeneticsApplicationTools::getSubstitutionModelSet(alphabet, sites, bppmixedlikelihoods.getParams()));
        if (modelSet==0){
          cout << "Missing a Mixed model" << endl;
          exit(0);
        }
        
        SiteContainerTools::changeGapsToUnknownCharacters(*sites);
        if (modelSet->getNumberOfStates() > modelSet->getAlphabet()->getSize())
          {
            //Markov-modulated Markov model!
            rDist = new ConstantDistribution(1.);
          }
        else
          {
            rDist = PhylogeneticsApplicationTools::getRateDistribution(bppmixedlikelihoods.getParams());
          }
        tl = new RNonHomogeneousMixedTreeLikelihood(*tree, *sites, modelSet, rDist, true);
        nbStates = modelSet->getNumberOfStates();
      }
    else throw Exception("Unknown option for nonhomogeneous: " + nhOpt);
    tl->initialize();
    
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

    
    ///////////////////////////////////////////////
    // Getting likelihoods per submodel

    string outputFile;
    outputFile = ApplicationTools::getAFilePath("output.likelihoods.file", bppmixedlikelihoods.getParams(), true, false);
    ApplicationTools::displayResult("Output file for likelihoods", outputFile);
    ofstream out(outputFile.c_str(), ios::out);

    unsigned int nSites=sites->getNumberOfSites();
    
    unsigned int nummodel=(unsigned int)ApplicationTools::getIntParameter("likelihoods.model_number", bppmixedlikelihoods.getParams(), 1, "", true, true);

    string parname=ApplicationTools::getStringParameter("likelihoods.parameter_name",bppmixedlikelihoods.getParams(), "", "", true, false);

    if (modelSet && ((nummodel<=0) || (nummodel>modelSet->getNumberOfModels()))){
      ApplicationTools::displayError("Bad number of model " + TextTools::toString(nummodel) + ".");
      exit(-1);
    }
    
    MixedSubstitutionModel* p0;      
    p0=dynamic_cast<MixedSubstitutionModel*>(model?model:modelSet->getModel(nummodel-1));

    if (p0==NULL){
      ApplicationTools::displayError("Model " + TextTools::toString(nummodel) + " is not a Mixed Model.");
      exit(-1);
    }
    
    if (dynamic_cast<AbstractBiblioMixedSubstitutionModel*>(p0)!=NULL)
      p0=dynamic_cast<AbstractBiblioMixedSubstitutionModel*>(p0)->getMixedModel();

    
    // Case of a MixtureOfSubstitutionModels
    
    MixtureOfSubstitutionModels* pMSM=dynamic_cast<MixtureOfSubstitutionModels*>(p0);
    if (pMSM!=NULL){
      vector<string> colNames;
      colNames.push_back("Sites");
      
      unsigned int nummod=pMSM->getNumberOfModels();
      for (unsigned int i=0;i<nummod;i++)
        colNames.push_back(pMSM->getNModel(i)->getName());
      
      DataTable* rates=new DataTable(nSites, colNames.size());
      rates->setColumnNames(colNames);
      
      for (unsigned int i=0;i<nSites;i++){
        const Site* currentSite = &sites->getSite(i);
        int currentSitePosition = currentSite->getPosition();
        (*rates)(i,"Sites")=string("[" + TextTools::toString(currentSitePosition) + "]");
      }
      
      Vdouble vprob=pMSM->getProbabilities();
      for (unsigned int i=0;i<nummod;i++){
        string modname=pMSM->getNModel(i)->getName();
        
        for (unsigned int j=0;j<nummod;j++)
          pMSM->setNProbability(j,(j==i)?1:0);
        
        if (tl)
          delete tl;
        
        if (nhOpt=="no")
          tl = new RHomogeneousMixedTreeLikelihood(*tree, *sites, model, rDist, true, false, true);
        else 
          tl = new RNonHomogeneousMixedTreeLikelihood(*tree, *sites, modelSet, rDist, false, true);
        
        tl->initialize();
        logL = tl->getValue();
        Vdouble Vd=tl->getLogLikelihoodForEachSite();
        for (unsigned int j=0;j<nSites;j++)
          (*rates)(j,modname)=TextTools::toString(Vd[j]);
        
        ApplicationTools::displayMessage("\n");
        ApplicationTools::displayMessage("Model " + modname + ":");
        ApplicationTools::displayResult("Log likelihood", TextTools::toString(tl->getValue(), 15)); 
        ApplicationTools::displayResult("Probability", TextTools::toString(vprob[i], 15)); 
      }
      
      DataTable::write(*rates, out, "\t");
    }

    // Case of a MixtureOfASubstitutionModel

    else {
      MixtureOfASubstitutionModel* pMSM2=dynamic_cast<MixtureOfASubstitutionModel*>(p0);
      if (pMSM2!=NULL){
        if (parname==""){
          ApplicationTools::displayError("Argument likelihoods.parameter_name is required.");
          exit(-1);
        }

        unsigned int nummod=pMSM2->getNumberOfModels();

        vector<vector<int> > vvnmod;
        unsigned int i2=0;
        while (i2<nummod){
          string par2=parname+"_"+TextTools::toString(i2+1);
          Vint vnmod=pMSM2->getSubmodelNumbers(par2);
          if (vnmod.size()==0)
            break;
          vvnmod.push_back(vnmod);
          i2++;
        }

        unsigned int nbcl=vvnmod.size();
        
        Vdouble vprob=pMSM2->getProbabilities();
        
        vector<vector<double> > vvprob;
        for (unsigned int i=0; i<nbcl; i++){
          vector<double> vprob2;
          for (unsigned int j=0;j<vvnmod[i].size();j++)
            vprob2.push_back(vprob[vvnmod[i][j]]);
          
          vvprob.push_back(vprob2);
        }
        
        vector<string> colNames;
        colNames.push_back("Sites");

        Vdouble dval;
        for (unsigned int i=0;i<nbcl;i++){
          SubstitutionModel* pSM=pMSM2->getNModel(vvnmod[i][0]);
          double valPar=pSM->getParameterValue(pSM->getParameterNameWithoutNamespace(parname));
          dval.push_back(valPar);
          colNames.push_back(parname+"="+TextTools::toString(valPar));
        }
        colNames.push_back("mean");
        
        DataTable* rates=new DataTable(nSites, colNames.size());
        rates->setColumnNames(colNames);
      
        for (unsigned int i=0;i<nSites;i++){
          const Site* currentSite = &sites->getSite(i);
          int currentSitePosition = currentSite->getPosition();
          (*rates)(i,"Sites")=string("[" + TextTools::toString(currentSitePosition) + "]");
        }

        VVdouble vvd;
        
        for (unsigned int i=0;i<nbcl;i++){
          string par2=parname+"_"+TextTools::toString(i+1);
          for (unsigned int j=0;j<nummod;j++)
            pMSM2->setNProbability(j,0);


          double s=VectorTools::sum(vvprob[i]);
          for (unsigned int j=0;j<vvprob[i].size();j++)
            pMSM2->setNProbability(vvnmod[i][j],vvprob[i][j]/s);
          
          if (tl)
            delete tl;
          
          if (nhOpt=="no")
            tl = new RHomogeneousMixedTreeLikelihood(*tree, *sites, model, rDist, true, false, true);
          else 
            tl = new RNonHomogeneousMixedTreeLikelihood(*tree, *sites, modelSet, rDist, false, true);
        
          tl->initialize();
          logL = tl->getValue();
          Vdouble vd=tl->getLogLikelihoodForEachSite();

          for (unsigned int j=0;j<nSites;j++)
            (*rates)(j,i+1)=TextTools::toString(vd[j]);

          vvd.push_back(vd);
          
          ApplicationTools::displayMessage("\n");
          ApplicationTools::displayMessage("Parameter " + par2 + ":");
          
          ApplicationTools::displayResult("Log likelihood", TextTools::toString(tl->getValue(), 15));
          ApplicationTools::displayResult("Conditional probability", TextTools::toString(s, 15)); 
        }

        for (unsigned int j=0;j<nSites;j++){
          Vdouble vd;
          for (unsigned int i=0;i<nbcl;i++)
            vd.push_back(vvd[i][j]);

          VectorTools::logNorm(vd);
          (*rates)(j,nbcl+1)=TextTools::toString(VectorTools::sumExp(vd,dval));
        }
        
        DataTable::write(*rates, out, "\t");
      }
    }

    delete alphabet;
    delete sites;
    if(model)    delete model;
    if(modelSet) delete modelSet;
    delete rDist;
    delete tl;
    delete tree;
    ApplicationTools::displayMessage("\n");
    bppmixedlikelihoods.done();
  }    
    
  catch (exception & e)
    {
      cout << e.what() << endl;
      return 1;
    }
  
  return 0;
}

