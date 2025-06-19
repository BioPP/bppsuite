// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

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
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>


using namespace bpp;

void calculateExpectedRatePerSite(const PhyloLikelihoodInterface& tl);
void normalizeVector(vector<double>& data);
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

    const std::map<size_t, std::shared_ptr<const AlignmentDataInterface >> mSites = PhylogeneticsApplicationTools::uniqueToSharedMap<const TemplateAlignmentDataInterface<string>>(mSitesuniq);

    /////// Get the map of initial trees

    auto mpTree = bppml.getPhyloTreesMap(mSites, unparsedParams);

    // Try to write the current tree to file. This will be overwritten
    // by the optimized tree, but allow to check file existence before
    // running optimization!

    vector<const PhyloTree*> vcpTree;

    for (const auto& pTree : mpTree)
    {
      vcpTree.push_back(pTree.second.get());
    }

    PhylogeneticsApplicationTools::writePhyloTrees(vcpTree, bppml.getParams(), "output.", "", true, false, true);


    /////////////////
    // Computing stuff

    std::shared_ptr<PhyloLikelihoodInterface> tl_new = 0;

    shared_ptr<SubstitutionProcessCollection> SPC;

    shared_ptr<PhyloLikelihoodContainer> mPhyl = 0;

    shared_ptr<BranchModelInterface>    model;
    shared_ptr<TransitionModelInterface>    tmodel; // for legacy
    shared_ptr<DiscreteDistributionInterface> rDist;

    SPC = bppml.getCollection(alphabet, gCode, mSites, mpTree, unparsedParams);

    auto mSeqEvoltmp = bppml.getProcesses(SPC, unparsedParams);

    auto mSeqEvol = PhylogeneticsApplicationTools::uniqueToSharedMap<SequenceEvolution>(mSeqEvoltmp);

    mPhyl = bppml.getPhyloLikelihoods(context, mSeqEvol, SPC, mSites);

    // retrieve Phylo 0, aka result phylolikelihood

    if (!mPhyl->hasPhyloLikelihood(0))
      throw Exception("Missing phyloLikelihoods.");

    tl_new = (*mPhyl)[0];


    ApplicationTools::displayMessage("");

    // Listing parameters
    string paramNameFile = ApplicationTools::getAFilePath("output.parameter_names.file", bppml.getParams(), false, false, "", true, "none", 1);

    if (paramNameFile != "none")
    {
      ApplicationTools::displayResult("List parameters to", paramNameFile);
      ofstream pnfile(paramNameFile.c_str(), ios::out);

      ParameterList pl = tl_new->getParameters();

      for (size_t i = 0; i < pl.size(); ++i)
      {
        pnfile << pl[i].getName() << endl;
      }
      pnfile.close();
      cout << "BppML's done." << endl;
      exit(0);
    }

    // Output trees
    string treeWIdPath = ApplicationTools::getAFilePath("output.tree_ids.file", bppml.getParams(), false, false, "", true, "none", 1);
    if (treeWIdPath != "none")
    {
      bppml.getParams()["output_ids.tree.file"] = treeWIdPath;

      PhylogeneticsApplicationTools::writePhyloTrees(*SPC, bppml.getParams(), "output_ids.", "", true, true, false, true);

      ApplicationTools::displayResult("Writing tagged tree to", treeWIdPath + "_...");
      cout << "BppML's done." << endl;
      exit(0);
    }


    // Check initial likelihood:

    bppml.fixLikelihood(alphabet, gCode, tl_new);

    tl_new = PhylogeneticsApplicationTools::optimizeParameters(tl_new, bppml.getParams());

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
      
      for (const auto& it2 : mSeqEvol)
      {
        PhylogeneticsApplicationTools::printParameters(*it2.second, out, it2.first);
        out.endLine();
      }

      PhylogeneticsApplicationTools::writePhyloTrees(*SPC, bppml.getParams(), "output.", "", true, true, false, false);
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
    calculateExpectedRatePerSite(*tl_new);

    bppml.done();
  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}

void calculateExpectedRatePerSite(const PhyloLikelihoodInterface& tl) {
  std::vector<double> res;
  auto& lik = dynamic_cast<const SingleProcessPhyloLikelihood&>(tl);
  std::shared_ptr<const SubstitutionProcessInterface> pSP = lik.getSubstitutionProcess();
  auto pDD = pSP->getRateDistribution();
  auto sites = lik.getData();

  for (size_t i = 0; i < sites->getNumberOfSites(); ++i) {
    auto post = lik.getPosteriorProbabilitiesForSitePerClass(i);
    double expectedR = 0;
    for (size_t j = 0; j < post.size(); ++j) {
      auto cat = pDD->getCategory(j);
      expectedR += post[j] * cat;
    }
    cout << "Expected rate for site " << i << ": " << expectedR << endl;
    res.push_back(expectedR);
  }
  normalizeVector(res);
  int i = 0;
  for (const double& entry : res) {
    cout << "Normalized Expected rate for site " << i << ": " << entry << endl;
    i++;
  }
}

void normalizeVector(vector<double>& data) {
    if (data.empty()) {
        throw invalid_argument("Input vector is empty.");
    }

    // Compute mean
    double sum = accumulate(data.begin(), data.end(), 0.0);
    double mean = sum / data.size();

    // Compute standard deviation
    double sq_sum = 0.0;
    for (double val : data) {
        sq_sum += (val - mean) * (val - mean);
    }
    double stdev = sqrt(sq_sum / data.size());

    if (stdev == 0.0) {
        throw runtime_error("Standard deviation is zero. Cannot normalize.");
    }

    // Normalize each element
    for (double& val : data) {
        val = (val - mean) / stdev;
    }
}