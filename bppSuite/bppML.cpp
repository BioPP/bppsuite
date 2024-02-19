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

// From bpp-core:
#include <Bpp/Version.h>
#include <Bpp/Numeric/DataTable.h>

// // From bpp-seq:
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/SiteTools.h>

// // From bpp-phyl:
#include <Bpp/Phyl/App/BppPhylogeneticsApplication.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Model/MixedTransitionModel.h>

using namespace bpp;

/******************************************************************************/

int main(int args, char** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*       Bio++ Maximum Likelihood Computation, version " << BPP_VERSION << "      *" << endl;
  cout << "*                                                                *" << endl;
  cout << "* Authors: J. Dutheil                       Last Modif. " << BPP_REL_DATE << " *" << endl;
  cout << "*          B. Boussau                                            *" << endl;
  cout << "*          L. Guéguen                                            *" << endl;
  cout << "*          M. Groussin                                           *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  try
  {
    BppPhylogeneticsApplication bppml(args, argv, "bppml");

    if (args == 1)
    {
      bppml.help("bppml");
      return 0;
    }

    bppml.startTimer();

    map<string, string> unparsedParams;

    Context context;
    
    ///// Alphabet

    std::shared_ptr<const Alphabet> alphabet(bppml.getAlphabet());

    /// GeneticCode
    
    auto gCode(bppml.getGeneticCode(alphabet));

    ////// Get the map of the sequences

    auto mSitesuniq = bppml.getConstAlignmentsMap(alphabet, true);

    const std::map<size_t, std::shared_ptr<const AlignmentDataInterface > > mSites = PhylogeneticsApplicationTools::uniqueToSharedMap<const TemplateAlignmentDataInterface<string>>(mSitesuniq);

    /////// Get the map of initial trees

    auto mpTree = bppml.getPhyloTreesMap(mSites, unparsedParams);

    // Try to write the current tree to file. This will be overwritten
    // by the optimized tree, but allow to check file existence before
    // running optimization!

    vector<const PhyloTree*> vcpTree;
    
    for (const auto& pTree : mpTree)
      vcpTree.push_back(pTree.second.get());
    
    PhylogeneticsApplicationTools::writePhyloTrees(vcpTree, bppml.getParams(),"output.","",true,false,true);


    /////////////////
    // Computing stuff

    std::shared_ptr<PhyloLikelihoodInterface> tl_new = 0;
    
    shared_ptr<SubstitutionProcessCollection> SPC;

    shared_ptr<PhyloLikelihoodContainer> mPhyl=0;
    
    shared_ptr<BranchModelInterface>    model;
    shared_ptr<TransitionModelInterface>    tmodel; // for legacy 
    shared_ptr<DiscreteDistributionInterface> rDist;

    SPC = bppml.getCollection(alphabet, gCode, mSites, mpTree, unparsedParams);

    auto mSeqEvoltmp = bppml.getProcesses(SPC, unparsedParams);
    
    auto mSeqEvol = PhylogeneticsApplicationTools::uniqueToSharedMap<SequenceEvolution>(mSeqEvoltmp);
    
    mPhyl=bppml.getPhyloLikelihoods(context, mSeqEvol, SPC, mSites);
    
    // retrieve Phylo 0, aka result phylolikelihood
    
    if (!mPhyl->hasPhyloLikelihood(0))
      throw Exception("Missing phyloLikelihoods.");
    
    tl_new=(*mPhyl)[0];

    
    ApplicationTools::displayMessage("");
    
    //Listing parameters
    string paramNameFile = ApplicationTools::getAFilePath("output.parameter_names.file", bppml.getParams(), false, false, "", true, "none", 1);

    if (paramNameFile != "none") {
      ApplicationTools::displayResult("List parameters to", paramNameFile);
      ofstream pnfile(paramNameFile.c_str(), ios::out);

      ParameterList pl=tl_new->getParameters();

      for (size_t i = 0; i < pl.size(); ++i) {
        pnfile << pl[i].getName() << endl;
      }
      pnfile.close();
      cout << "BppML's done." << endl;
      exit(0);
    }

    //Output trees
    string treeWIdPath = ApplicationTools::getAFilePath("output.tree_ids.file", bppml.getParams(), false, false, "", true, "none", 1);
    if (treeWIdPath != "none")
    {
      bppml.getParams()["output_ids.tree.file"]=treeWIdPath;
      
      PhylogeneticsApplicationTools::writePhyloTrees(*SPC, bppml.getParams(), "output_ids.", "", true, true, false, true);

      ApplicationTools::displayResult("Writing tagged tree to", treeWIdPath + "_...");
      cout << "BppML's done." << endl;
      exit(0);
    }


    //Check initial likelihood:
      
    bppml.fixLikelihood(alphabet, gCode, tl_new);
    
    // First `true` means that default is to optimize model parameters.
    if(ApplicationTools::getBooleanParameter("optimization.model_parameters", bppml.getParams(), true, "", true, 1))
      tl_new=PhylogeneticsApplicationTools::optimizeParameters(tl_new, tl_new->getParameters(), bppml.getParams());
    else
      tl_new=PhylogeneticsApplicationTools::optimizeParameters(tl_new, tl_new->getBranchLengthParameters(), bppml.getParams());
    
    SPC->matchParametersValues(tl_new->getParameters());
    
    PhylogeneticsApplicationTools::writePhyloTrees(*SPC, bppml.getParams(), "output.", "", true, true, true);
    
    
    // Write parameters to screen:
    bppml.displayParameters(*tl_new);
    
    // Checking convergence:
    PhylogeneticsApplicationTools::checkEstimatedParameters(tl_new->getParameters());
    
    // Write parameters to file:
    string parametersFile = ApplicationTools::getAFilePath("output.estimates", bppml.getParams(), false, false);
    bool withAlias = ApplicationTools::getBooleanParameter("output.estimates.withalias", bppml.getParams(), true, "", false, 1);
    
    ApplicationTools::displayResult("output.estimates", parametersFile);
    
    if (parametersFile != "none")
    {
      StlOutputStream out(make_unique<ofstream>(parametersFile.c_str(), ios::out));
      
      PhylogeneticsApplicationTools::printParameters(*mPhyl, out);

      PhylogeneticsApplicationTools::printParameters(*SPC, out, 1, withAlias);
      
      for (const auto& it2:mSeqEvol)
      {
        PhylogeneticsApplicationTools::printParameters(*it2.second, out, it2.first);
        out.endLine();
      }
      
      PhylogeneticsApplicationTools::writePhyloTrees(*SPC, bppml.getParams(), "output.", "",true,true,false,false);
        
    }

    // Write infos to file:
    //     probabilities of rate discrete distributions
    //     site infos : lnL, class (or process in case of collection) posterior probability distribution
    
    string infosFile = ApplicationTools::getAFilePath("output.infos", bppml.getParams(), false, false);
    if (infosFile != "none")
    {
      ApplicationTools::displayResult("Alignment information logfile", infosFile);
      PhylogeneticsApplicationTools::printAnalysisInformation(*mPhyl, infosFile);
    }

    bppml.done();
  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    return 1;
  }

return 0;
}

