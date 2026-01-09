// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

// From the STL:
#include <iostream>
#include <iomanip>

using namespace std;

// From bpp-core:
#include <Bpp/Version.h>
#include <Bpp/Numeric/Prob/DiscreteDistribution.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Text/KeyvalTools.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/Io/IoDistanceMatrixFactory.h>
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/Io/IoDistanceMatrixFactory.h>


// From bpp-phyl:
#include <Bpp/Phyl/Tree/Tree.h>
#include <Bpp/Phyl/PatternTools.h>
#include <Bpp/Phyl/App/BppPhylogeneticsApplication.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Io/Newick.h>
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

int main(int args, char** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*              Bio++ Distance Methods, version " << BPP_VERSION << "             *" << endl;
  cout << "* Author: J. Dutheil                        Created     05/05/07 *" << endl;
  cout << "*                                           Last Modif. " << BPP_REL_DATE << " *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  if (args == 1)
  {
    help();
    return 0;
  }

  try
  {
    Context context;

    BppPhylogeneticsApplication bppdist(args, argv, "BppDist");
    bppdist.startTimer();

    std::shared_ptr<const Alphabet> alphabet = bppdist.getAlphabet();

    /// GeneticCode

    auto gCode(bppdist.getGeneticCode(alphabet));

    ////// Get the map of the sequences

    auto mSitesuniq = bppdist.getConstAlignmentsMap(alphabet, true);

    std::map<size_t, std::shared_ptr<const AlignmentDataInterface >> mSites = PhylogeneticsApplicationTools::uniqueToSharedMap<const TemplateAlignmentDataInterface<string>>(mSitesuniq);

    if (mSites.find(1) == mSites.end())
      throw Exception("Missing number 1 for data.");

    auto align = mSites[1];

    // Build a random tree from data 1, used to construct processes
    map<string, string> unparsedparams;

    map<string, string> args2;
    args2["input.tree1"] = "random(data=1)";

    auto mpTree = PhylogeneticsApplicationTools::getPhyloTrees(args2, mSites, unparsedparams);

    // process & phyloli
    auto SPC = std::shared_ptr<SubstitutionProcessCollection>(bppdist.getCollection(alphabet, gCode, mSites, mpTree, unparsedparams).release());

    /// Only simple processes
    auto mSeqEvoltmp = bppdist.getProcesses(SPC, unparsedparams);

    map<size_t, shared_ptr<SequenceEvolution>> mSeqEvol = PhylogeneticsApplicationTools::uniqueToSharedMap<SequenceEvolution>(mSeqEvoltmp);

    args2["phylo1"] = "Single(process=1, data=1)";
    auto mPhyl = PhylogeneticsApplicationTools::getPhyloLikelihoodContainer(context, SPC, mSeqEvol, mSites, args2);

    //////////////////////////////////////////////////////////
    /// Now get distEstimation
    // Only Model 1 is considered

    if (!SPC->hasSubstitutionProcessNumber(1))
      throw Exception("Missing process 1.");
    auto process = SPC->getSubstitutionProcess(1);

    DistanceEstimation distEstimation(process, align, 1, false); // !process is shared between distEstimation & collection

    /// Dist Method
    string method = ApplicationTools::getStringParameter("method", bppdist.getParams(), "nj");
    ApplicationTools::displayResult("Tree reconstruction method", method);

    unique_ptr<TreeTemplate<Node>> tree;
    unique_ptr<AgglomerativeDistanceMethodInterface> distMethod;
    if (method == "wpgma")
    {
      distMethod = make_unique<PGMA>(true);
    }
    else if (method == "upgma")
    {
      distMethod = make_unique<PGMA>(false);
    }
    else if (method == "nj")
    {
      auto nj = make_unique<NeighborJoining>();
      nj->outputPositiveLengths(true);
      distMethod = std::move(nj);
    }
    else if (method == "bionj")
    {
      auto bionj = make_unique<BioNJ>();
      bionj->outputPositiveLengths(true);
      distMethod = std::move(bionj);
    }
    else
      throw Exception("Unknown tree reconstruction method.");

    string type = ApplicationTools::getStringParameter("optimization.method", bppdist.getParams(), "init");

    ApplicationTools::displayResult("Model parameters estimation method", type);
    if (type == "init")
      type = OptimizationTools::DISTANCEMETHOD_INIT;
    else if (type == "pairwise")
      type = OptimizationTools::DISTANCEMETHOD_PAIRWISE;
    else if (type == "iterations")
      type = OptimizationTools::DISTANCEMETHOD_ITERATIONS;
    else
      throw Exception("Unknown parameter estimation procedure '" + type + "'.");

    //// Optimization

    auto lcp = std::make_shared<LikelihoodCalculationSingleProcess>(context, align, process);
    auto lik = std::make_shared<SingleProcessPhyloLikelihood>(context, lcp);

    auto params = bppdist.getParams();

    params["optimization"] = ApplicationTools::getStringParameter("optimization", params, "D-Brent(derivatives=Newton)");

    OptimizationTools::OptimizationOptions optopt(lik, params);

    tree = OptimizationTools::buildDistanceTree(distEstimation, *distMethod, type, optopt);

    //// Output
    ApplicationTools::warning = ApplicationTools::message;

    string matrixPath = ApplicationTools::getAFilePath("output.matrix.file", bppdist.getParams(), false, false, "", false);
    if (matrixPath != "none")
    {
      ApplicationTools::displayResult("Output matrix file", matrixPath);
      string matrixFormat = ApplicationTools::getStringParameter("output.matrix.format", bppdist.getParams(), IODistanceMatrixFactory::PHYLIP_FORMAT);
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

      auto odm = IODistanceMatrixFactory().createWriter(IODistanceMatrixFactory::PHYLIP_FORMAT, extended);
      odm->writeDistanceMatrix(*distEstimation.getMatrix(), matrixPath, true);
    }

    PhylogeneticsApplicationTools::writeTree(*tree, bppdist.getParams());

    // Output some parameters:
    string parametersFile = ApplicationTools::getAFilePath("output.estimates", bppdist.getParams(), false, false);
    bool withAlias = ApplicationTools::getBooleanParameter("output.estimates.withalias", bppdist.getParams(), true, "", false, 1);

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
    }

    // Bootstrap:
    unsigned int nbBS = ApplicationTools::getParameter<unsigned int>("bootstrap.number", bppdist.getParams(), 0);
    if (nbBS > 0)
    {
      auto sites = dynamic_pointer_cast<const SiteContainerInterface>(align);

      if (sites == 0)
        throw Exception("bppDist: bootstrap yet only for sites. Ask developpers");

      ApplicationTools::displayResult("Number of bootstrap samples", TextTools::toString(nbBS));
      bool approx = ApplicationTools::getBooleanParameter("bootstrap.approximate", bppdist.getParams(), true);
      ApplicationTools::displayResult("Use approximate bootstrap", TextTools::toString(approx ? "yes" : "no"));
      if (approx)
        type = OptimizationTools::DISTANCEMETHOD_INIT;

      uint bootstrapVerbose = ApplicationTools::getBooleanParameter("bootstrap.verbose", bppdist.getParams(), false, "", true, false);
      optopt.verbose = bootstrapVerbose;

      string bsTreesPath = ApplicationTools::getAFilePath("bootstrap.output.file", bppdist.getParams(), false, false);
      shared_ptr<ofstream> out = NULL;
      if (bsTreesPath != "none")
      {
        ApplicationTools::displayResult("Bootstrap trees stored in file", bsTreesPath);
        out = make_shared<ofstream>(bsTreesPath.c_str(), ios::out);
      }
      Newick newick;

      vector<std::unique_ptr<TreeTemplate<Node>>> bsTrees(nbBS);
      ApplicationTools::displayTask("Bootstrapping", true);
      // auto vMod = process->getModelNumbers();

      for (unsigned int i = 0; i < nbBS; i++)
      {
        ApplicationTools::displayGauge(i, nbBS - 1, '=');
        auto sample = SiteContainerTools::bootstrapSites(*sites);
        // if (approx)
        //   for (auto modi:vMod)
        //   {
        //     auto& model = dynamic_pointer_cast<AbstractSubstitutionModel>(process->getModel(modi));
        //     if (model)
        //       model->setFreqFromData(*sample);
        //   }
        distEstimation.setData(shared_ptr<SiteContainerInterface>(sample.release()));
        auto procMb3 = dynamic_pointer_cast<const SubstitutionProcessCollectionMember>(distEstimation.getProcess());

        bsTrees[i].reset(OptimizationTools::buildDistanceTree(distEstimation, *distMethod, type, optopt).release());
        if (out && i == 0)
          newick.writeTree(*bsTrees[i], bsTreesPath, true);
        if (out && i >  0)
          newick.writeTree(*bsTrees[i], bsTreesPath, false);
      }
      if (out)
        out->close();
      ApplicationTools::displayTaskDone();
      ApplicationTools::displayTask("Compute bootstrap values");
      //      TreeTools::computeBootstrapValues(*tree, bsTrees);
      ApplicationTools::displayTaskDone();

      // Write resulting tree:
      PhylogeneticsApplicationTools::writeTree(*tree, bppdist.getParams());
    }

    bppdist.done();
  }

  catch (exception& e)
  {
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}
