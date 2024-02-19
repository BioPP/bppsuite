//
// File: bppDist.cpp
// Created by: Julien Dutheil
// Created on: May Sat 05 15:09 2007
// From file bppML.cpp
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

// From bpp-core:
#include <Bpp/Version.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/App/BppSequenceApplication.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Text/KeyvalTools.h>


// From PhylLib:
#include <Bpp/Phyl/Tree/Tree.h>
#include <Bpp/Phyl/PatternTools.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Io/Newick.h>
#include <Bpp/Phyl/Io/IoDistanceMatrixFactory.h>
#include <Bpp/Phyl/Distance/DistanceEstimation.h>
#include <Bpp/Phyl/Distance/PGMA.h>
#include <Bpp/Phyl/Distance/BioNJ.h>
#include <Bpp/Phyl/OptimizationTools.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>

using namespace bpp;

void help()
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
  (*ApplicationTools::message << "bppdist parameter1_name=parameter1_value parameter2_name=parameter2_value").endLine();
  (*ApplicationTools::message << "      ... param=option_file").endLine();
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "  Refer to the Bio++ Program Suite Manual for a list of available options.").endLine();
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
}

int main(int args, char ** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*              Bio++ Distance Methods, version " << BPP_VERSION << "             *" << endl;
  cout << "* Author: J. Dutheil                        Created     05/05/07 *" << endl;
  cout << "*                                           Last Modif. " << BPP_REL_DATE << " *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  if(args == 1)
  {
    help();
    return 0;
  }
  
  try {

    BppSequenceApplication bppdist(args, argv, "BppDist");
    bppdist.startTimer();

    std::shared_ptr<const Alphabet> alphabet = bppdist.getAlphabet();
    
    /// GeneticCode
    
    auto gCode(bppdist.getGeneticCode(alphabet));

    //sites
    
    auto allSites = SequenceApplicationTools::getSiteContainer(alphabet, bppdist.getParams());
  
    shared_ptr<VectorSiteContainer> sites = SequenceApplicationTools::getSitesToAnalyse(* allSites, bppdist.getParams());
    
    // model
    
    map<string, string> unparsedparams;

    std::shared_ptr<BranchModelInterface> model = PhylogeneticsApplicationTools::getBranchModel(alphabet, gCode, sites, bppdist.getParams(), unparsedparams);

    if (!model)
      ApplicationTools::displayMessage("bppDist available only for regular transition models, not the given one.");

    
    std::shared_ptr<DiscreteDistributionInterface> rDist;
    if (model->getNumberOfStates() > model->getAlphabet()->getSize())
    {
      //Markov-modulated Markov model!
      rDist = std::shared_ptr<ConstantRateDistribution>();
    }
    else
    {
      rDist = std::shared_ptr<DiscreteDistributionInterface>(PhylogeneticsApplicationTools::getRateDistribution(bppdist.getParams()));
    }
   
    DistanceEstimation distEstimation(model, rDist, sites, 1, false);
 
    string method = ApplicationTools::getStringParameter("method", bppdist.getParams(), "nj");
    ApplicationTools::displayResult("Tree reconstruction method", method);

    unique_ptr<TreeTemplate<Node>> tree;
    unique_ptr<AgglomerativeDistanceMethodInterface> distMethod;
    if(method == "wpgma")
    {
      distMethod = make_unique<PGMA>(true);
    }
    else if(method == "upgma")
    {
      distMethod = make_unique<PGMA>(false);
    }
    else if(method == "nj")
    {
      auto nj = make_unique<NeighborJoining>();
      nj->outputPositiveLengths(true);
      distMethod = move(nj);
    }
    else if(method == "bionj")
    {
      auto bionj = make_unique<BioNJ>();
      bionj->outputPositiveLengths(true);
      distMethod = move(bionj);
    }
    else throw Exception("Unknown tree reconstruction method.");
  
    string type = ApplicationTools::getStringParameter("optimization.method", bppdist.getParams(), "init");
    ApplicationTools::displayResult("Model parameters estimation method", type);
    if (type == "init") type = OptimizationTools::DISTANCEMETHOD_INIT;
    else if (type == "pairwise") type = OptimizationTools::DISTANCEMETHOD_PAIRWISE;
    else if (type == "iterations") type = OptimizationTools::DISTANCEMETHOD_ITERATIONS;
    else throw Exception("Unknown parameter estimation procedure '" + type + "'.");
  
    unsigned int optVerbose = ApplicationTools::getParameter<unsigned int>("optimization.verbose", bppdist.getParams(), 2);
	
    string mhPath = ApplicationTools::getAFilePath("optimization.message_handler", bppdist.getParams(), false, false);
    auto messenger = shared_ptr<OutputStream>(
      (mhPath == "none") ? 0 :
      (mhPath == "std") ? ApplicationTools::message.get() :
      new StlOutputStream(make_unique<ofstream>(mhPath.c_str(), ios::out)));
    
    ApplicationTools::displayResult("Message handler", mhPath);

    string prPath = ApplicationTools::getAFilePath("optimization.profiler", bppdist.getParams(), false, false);
    auto profiler = shared_ptr<OutputStream>(
      (prPath == "none") ? 0 :
      (prPath == "std") ? ApplicationTools::message.get() :
      new StlOutputStream(make_unique<ofstream>(prPath.c_str(), ios::out)));
    
    if(profiler) profiler->setPrecision(20);
    ApplicationTools::displayResult("Profiler", prPath);

    // Should I ignore some parameters?
    ParameterList allParameters = model->getParameters();
    allParameters.addParameters(rDist->getParameters());
    ParameterList parametersToIgnore;
    string paramListDesc = ApplicationTools::getStringParameter("optimization.ignore_parameter", bppdist.getParams(), "", "", true, false);
    bool ignoreBrLen = false;
    StringTokenizer st(paramListDesc, ",");
    while (st.hasMoreToken())
    {
      try
      {
        string param = st.nextToken();
        if (param == "BrLen")
          ignoreBrLen = true;
        else
        {
          if (allParameters.hasParameter(param))
          {
            Parameter* p = &allParameters.parameter(param);
            parametersToIgnore.addParameter(*p);
          }
          else ApplicationTools::displayWarning("Parameter '" + param + "' not found."); 
        }
      } 
      catch (ParameterNotFoundException& pnfe)
      {
        ApplicationTools::displayError("Parameter '" + pnfe.parameter() + "' not found, and so can't be ignored!");
      }
    }
	
    unsigned int nbEvalMax = ApplicationTools::getParameter<unsigned int>("optimization.max_number_f_eval", bppdist.getParams(), 1000000);
    ApplicationTools::displayResult("Max # ML evaluations", TextTools::toString(nbEvalMax));
	
    double tolerance = ApplicationTools::getDoubleParameter("optimization.tolerance", bppdist.getParams(), .000001);
    ApplicationTools::displayResult("Tolerance", TextTools::toString(tolerance));
	
    //Here it is:
    ofstream warn("warnings", ios::out);
    ApplicationTools::warning=std::shared_ptr<OutputStream>(dynamic_cast<OutputStream*>(new StlOutputStreamWrapper(&warn))); 
   tree = OptimizationTools::buildDistanceTree(distEstimation, *distMethod, parametersToIgnore, !ignoreBrLen, type, tolerance, nbEvalMax, profiler, messenger, optVerbose);
    warn.close();
//    delete ApplicationTools::warning;
    ApplicationTools::warning = ApplicationTools::message;

    string matrixPath = ApplicationTools::getAFilePath("output.matrix.file", bppdist.getParams(), false, false, "", false);
    if (matrixPath != "none")
    {
      ApplicationTools::displayResult("Output matrix file", matrixPath);
      string matrixFormat = ApplicationTools::getAFilePath("output.matrix.format", bppdist.getParams(), false, false, "", false);
      string format = "";
      bool extended = false;
      std::map<std::string, std::string> unparsedArguments_;
      KeyvalTools::parseProcedure(matrixFormat, format, unparsedArguments_);
      if (unparsedArguments_.find("type") != unparsedArguments_.end())
      {
        if (unparsedArguments_["type"] == "extended")
        {
          extended = true;
        }     
        else if (unparsedArguments_["type"] == "classic")
          extended = false;
        else
          ApplicationTools::displayWarning("Argument '" +
                                           unparsedArguments_["type"] + "' for parameter 'Phylip#type' is unknown. " +
                                           "Default used instead: not extended.");
      }    
      else
        ApplicationTools::displayWarning("Argument 'Phylip#type' not found. Default used instead: not extended.");
    
      ODistanceMatrix* odm = IODistanceMatrixFactory().createWriter(IODistanceMatrixFactory::PHYLIP_FORMAT, extended);
      odm->writeDistanceMatrix(*distEstimation.getMatrix(), matrixPath, true);
      delete odm;
    }
    PhylogeneticsApplicationTools::writeTree(*tree, bppdist.getParams());
  
    //Output some parameters:
    if (type == OptimizationTools::DISTANCEMETHOD_ITERATIONS)
    {
      // Write parameters to screen:
      ParameterList parameters = model->getParameters();
      for (unsigned int i = 0; i < parameters.size(); i++)
      {
        ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
      }
      parameters = rDist->getParameters();
      for (unsigned int i = 0; i < parameters.size(); i++)
      {
        ApplicationTools::displayResult(parameters[i].getName(), TextTools::toString(parameters[i].getValue()));
      }
      // Write parameters to file:
      string parametersFile = ApplicationTools::getAFilePath("output.estimates", bppdist.getParams(), false, false);
      if (parametersFile != "none")
      {
        ofstream out(parametersFile.c_str(), ios::out);
        parameters = model->getParameters();
        for (unsigned int i = 0; i < parameters.size(); i++)
        {
          out << parameters[i].getName() << " = " << parameters[i].getValue() << endl;
        }
        parameters = rDist->getParameters();
        for (unsigned int i = 0; i < parameters.size(); i++)
        {
          out << parameters[i].getName() << " = " << parameters[i].getValue() << endl;
        }
        out.close();
      }
    }
 
    //Bootstrap:
    unsigned int nbBS = ApplicationTools::getParameter<unsigned int>("bootstrap.number", bppdist.getParams(), 0);
    if(nbBS > 0)
    {
      ApplicationTools::displayResult("Number of bootstrap samples", TextTools::toString(nbBS));
      bool approx = ApplicationTools::getBooleanParameter("bootstrap.approximate", bppdist.getParams(), true);
      ApplicationTools::displayResult("Use approximate bootstrap", TextTools::toString(approx ? "yes" : "no"));
      if(approx)
      {
        type = OptimizationTools::DISTANCEMETHOD_INIT;
        parametersToIgnore = allParameters;
        ignoreBrLen = true;
      }
      bool bootstrapVerbose = ApplicationTools::getBooleanParameter("bootstrap.verbose", bppdist.getParams(), false, "", true, false);
 
      string bsTreesPath = ApplicationTools::getAFilePath("bootstrap.output.file", bppdist.getParams(), false, false);
      ofstream *out = NULL;
      if(bsTreesPath != "none")
      {
        ApplicationTools::displayResult("Bootstrap trees stored in file", bsTreesPath);
        out = new ofstream(bsTreesPath.c_str(), ios::out);
      }
      Newick newick;
    
      vector<std::unique_ptr<Tree> > bsTrees(nbBS);
      ApplicationTools::displayTask("Bootstrapping", true);
      for(unsigned int i = 0; i < nbBS; i++)
      {
        ApplicationTools::displayGauge(i, nbBS-1, '=');
        shared_ptr<VectorSiteContainer> sample = SiteContainerTools::bootstrapSites(*sites);
        auto tm=dynamic_pointer_cast<TransitionModelInterface>(model);
        if (approx && tm)
          tm->setFreqFromData(*sample);
        distEstimation.setData(sample);
        bsTrees[i] = OptimizationTools::buildDistanceTree(
          distEstimation,
          *distMethod,
          parametersToIgnore,
          ignoreBrLen,
          type,
          tolerance,
          nbEvalMax,
          NULL,
          NULL,
          (bootstrapVerbose ? 1 : 0)
          );
        if(out && i == 0) newick.writeTree(*bsTrees[i], bsTreesPath, true);
        if(out && i >  0) newick.writeTree(*bsTrees[i], bsTreesPath, false);
      }
      if(out) out->close();
      if(out) delete out;
      ApplicationTools::displayTaskDone();
      ApplicationTools::displayTask("Compute bootstrap values");
      TreeTools::computeBootstrapValues(*tree, bsTrees);
      ApplicationTools::displayTaskDone();

      //Write resulting tree:
      PhylogeneticsApplicationTools::writeTree(*tree, bppdist.getParams());
    }
    
    bppdist.done();
  }
      
  catch(exception & e)
  {
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}

