// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

// From the STL:
#include <iostream>
#include <iomanip>

using namespace std;

// From bpp-core:
#include <Bpp/Version.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Text/TextTools.h>
#include <Bpp/Numeric/DataTable.h>
#include <Bpp/Numeric/VectorTools.h>

// From bpp-seq:
#include <Bpp/Seq/Alphabet/Alphabet.h>

// From bpp-phyl:
#include <Bpp/Phyl/App/BppPhylogeneticsApplication.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Tree/Tree.h>
#include <Bpp/Phyl/Tree/TreeTemplate.h>
#include <Bpp/Phyl/Model/AbstractBiblioMixedTransitionModel.h>
#include <Bpp/Phyl/Model/MixedTransitionModel.h>
#include <Bpp/Phyl/Model/MixtureOfATransitionModel.h>
#include <Bpp/Phyl/Model/MixtureOfTransitionModels.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>

using namespace bpp;

/******************************************************************************/

int main(int args, char** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*     Bio++ Computation of site likelihoods inside mixed models  *" << endl;
  cout << "*                        Version " << BPP_VERSION << ".                          *" << endl;
  cout << "* Author: L. Guéguen                       Last Modif.: " << BPP_REL_DATE << " *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  try
  {
    BppPhylogeneticsApplication bppmixedlikelihoods(args, argv, "BppMixedLikelihoods");

    if (args == 1)
    {
      bppmixedlikelihoods.help("bppmixedlikelihoods");
      return 0;
    }

    bppmixedlikelihoods.startTimer();

    Context context;

    ///// Alphabet

    shared_ptr<const Alphabet> alphabet(bppmixedlikelihoods.getAlphabet());

    /// GeneticCode

    shared_ptr<const GeneticCode> gCode(bppmixedlikelihoods.getGeneticCode(alphabet));

    // get the data

    auto mSitesuniq = bppmixedlikelihoods.getConstAlignmentsMap(alphabet, true);

    const std::map<size_t, std::shared_ptr<const AlignmentDataInterface >> mSites = PhylogeneticsApplicationTools::uniqueToSharedMap<const TemplateAlignmentDataInterface<string>>(mSitesuniq);

    if (mSites.size() != 1)
      throw Exception("Only one alignment possible.");

    auto sites = mSites.begin()->second;

    /////// Get the map of initial trees

    map<string, string> unparsedParams;

    auto mpTree = bppmixedlikelihoods.getPhyloTreesMap(mSites, unparsedParams);

    /////////////////
    // Computing stuff


    shared_ptr<SubstitutionProcessCollection> SPC(bppmixedlikelihoods.getCollection(alphabet, gCode, mSites, mpTree, unparsedParams));

    auto mSeqEvoltmp = bppmixedlikelihoods.getProcesses(SPC, unparsedParams);

    auto mSeqEvol = PhylogeneticsApplicationTools::uniqueToSharedMap<SequenceEvolution>(mSeqEvoltmp);

    auto mPhyl(bppmixedlikelihoods.getPhyloLikelihoods(context, mSeqEvol, SPC, mSites));

    if (!mPhyl->hasPhyloLikelihood(0))
      throw Exception("Missing phyloLikelihoods.");

    auto tl = dynamic_pointer_cast<AlignedPhyloLikelihoodInterface>((*mPhyl)[0]);

    if (tl == 0)
      throw Exception("Only possible on aligned phyloLikelihood.");

    // Check initial likelihood:

    bppmixedlikelihoods.fixLikelihood(alphabet, gCode, tl);

    // /////////////////////////////////////////////
    // Getting likelihoods per submodel

    string outputFile;
    outputFile = ApplicationTools::getAFilePath("output.likelihoods.file", bppmixedlikelihoods.getParams(), true, false);
    ApplicationTools::displayResult("Output file for likelihoods", outputFile);
    ofstream out(outputFile.c_str(), ios::out);

    size_t nSites = sites->getNumberOfSites();

    // Look for submodel

    size_t modNum(0);

    vector<size_t> vNumMix; // numbers of mixed models in SPC
    auto nMod = SPC->getModelNumbers();
    for (auto n:nMod)
    {
      if (dynamic_pointer_cast<const MixedTransitionModelInterface>(SPC->getModel(n)) != NULL)
        vNumMix.push_back(n);
    }

    // Set model number
    if (vNumMix.size() == 0)
      throw Exception("No mixture models found.");

    // parname: name of the mixed parameter inside the model
    // realparname: name of the parameter as seen by the user (& the process)

    string realparname = ApplicationTools::getStringParameter("likelihoods.parameter_name", bppmixedlikelihoods.getParams(), "", "", true, false);

    string parname = realparname;

    // In case there is a model number at the end of the parameter name
    size_t posund = realparname.find_last_of("_");
    if (posund != string::npos)
    {
      try
      {
        modNum = (size_t)TextTools::toInt(realparname.substr(posund + 1));
        parname = realparname.substr(0, posund);
      }
      catch (exception& e)
      { }
    }


    if (vNumMix.size() == 1)
      modNum = vNumMix[0];
    else
    {
      size_t modNum2 = ApplicationTools::getParameter<size_t>("likelihoods.model_number", bppmixedlikelihoods.getParams(), 0, "", true, true);
      if (modNum && modNum2 != modNum)
        throw Exception("Mismatch numbers between parameter name " + TextTools::toString(modNum) + " and model " + TextTools::toString(modNum2));
      modNum = modNum2;
    }

    if (!modNum)
    // look for model number
    {
      for (auto n:vNumMix)
      {
        auto mod = SPC->getModel(n);
        if (mod->hasParameter(mod->getParameterNameWithoutNamespace(parname)))
        {
          // Check it is a MixtureOfATransitionModel
          bool modok = (dynamic_pointer_cast<const MixtureOfATransitionModel>(mod) != NULL);
          if (!modok)
          {
            auto ptmp = dynamic_pointer_cast<const AbstractBiblioMixedTransitionModel>(mod);
            if (ptmp)
              modok = (dynamic_cast<const MixtureOfATransitionModel*>(&ptmp->mixedModel()) != NULL);
          }

          if (!modok)
            continue;

          if (modNum != 0)
            throw Exception("Ambiguous model numbers for parameter " + parname + ":" + TextTools::toString(modNum) + " & " + TextTools::toString(n));
          else
            modNum = n;
        }
      }
      if (modNum == 0)
        throw Exception("Unknown parameter " + realparname);

      realparname = parname + "_" + TextTools::toString(modNum);
    }


    // at that point model number is known

    // Get the model used to compute the likelihood
    auto model = mPhyl->getCollectionNodes()->getModel(modNum);
    if (!model)
    {
      ApplicationTools::displayError("Unknown number of model " + TextTools::toString(modNum) + ".");
      exit(-1);
    }

    auto mixmodel = dynamic_cast<const MixedTransitionModelInterface*>(model->targetValue().get());

    if (!mixmodel)
    {
      ApplicationTools::displayError("Model " + TextTools::toString(modNum) + " is not a Mixed Model.");
      exit(-1);
    }

    //////////////////////////////////////////////////
    // Case of a MixtureOfTransitionModels

    auto pMSM = dynamic_cast<const MixtureOfTransitionModels*>(mixmodel);

    if (pMSM)
    {
      vector<string> colNames;
      colNames.push_back("Sites");
      colNames.push_back("Ll");

      size_t nummod = pMSM->getNumberOfModels();
      for (unsigned int i = 0; i < nummod; i++)
      {
        colNames.push_back("Ll_" + pMSM->getNModel(i)->getName());
      }
      for (unsigned int i = 0; i < nummod; i++)
      {
        colNames.push_back("Pr_" + pMSM->getNModel(i)->getName());
      }

      auto rates = make_shared<DataTable>(nSites, colNames.size());
      rates->setColumnNames(colNames);

      // output sites
      for (unsigned int i = 0; i < nSites; i++)
      {
        const auto& currentSite = sites->site(i);
        int currentSitePosition = currentSite.getCoordinate();
        (*rates)(i, "Sites") = string(TextTools::toString(currentSitePosition));
      }

      // Likelihoods
      auto vl = tl->getLikelihoodPerSite();
      for (unsigned int j = 0; j < nSites; j++)
      {
        (*rates)(j, "Ll") = TextTools::toString(log(vl[j]));
      }


      Vdouble vprob = pMSM->getProbabilities();
      vector<vector<DataLik>> vvd;

      for (unsigned int i = 0; i < nummod; i++)
      {
        string modname = pMSM->getNModel(i)->getName() + "_" + TextTools::toString(i + 1);

        auto func = [nummod, i](shared_ptr<BranchModelInterface> bmodel){
              auto pAbmtm2 = dynamic_pointer_cast<AbstractBiblioMixedTransitionModel>(bmodel);
              if (pAbmtm2)
              {
                for (unsigned int j = 0; j < nummod; ++j)
                {
                  pAbmtm2->setNProbability(j, (j == i) ? 1 : 0);
                }
              }
              else
              {
                auto pMSM2 = dynamic_pointer_cast<MixtureOfTransitionModels>(bmodel);
                if (!pMSM2)
                  throw Exception("Not mixed model " + bmodel->getName());

                for (unsigned int j = 0; j < nummod; j++)
                {
                  pMSM2->setNProbability(j, (j == i) ? 1 : 0);
                }
              }
            };

        model->modify(func, false);
        auto vd = tl->getLikelihoodPerSite();

        for (unsigned int j = 0; j < nSites; j++)
        {
          (*rates)(j, "Ll_" + modname) = TextTools::toString(log(vd[j]));
        }

        vvd.push_back(vd);

        ApplicationTools::displayMessage("\n");
        ApplicationTools::displayMessage("Model " + modname + ":");
        ApplicationTools::displayResult("Log likelihood", TextTools::toString(tl->getValue(), 15));
        ApplicationTools::displayResult("Probability", TextTools::toString(vprob[i], 15));
      }

      for (size_t j = 0; j < nSites; j++)
      {
        Vdouble vd;
        for (size_t i = 0; i < nummod; i++)
        {
          vd.push_back(std::log(vprob[i]) + log(vvd[i][j]));
        }

        VectorTools::logNorm(vd);
        for (size_t i = 0; i < nummod; i++)
        {
          (*rates)(j, "Pr_" + pMSM->getNModel(i)->getName()) = TextTools::toString(std::exp(vd[i]));
        }
      }

      DataTable::write(*rates, out, "\t");
    }

    //////////////////////////////////////////////////
    // Case of a MixtureOfASubstitutionModel

    else
    {
      // look for parameter name in MixtureOfATransitionModel
      const AbstractBiblioMixedTransitionModel* pAbmtm(0);

      auto pMatm = (pAbmtm = dynamic_cast<const AbstractBiblioMixedTransitionModel*>(mixmodel)) ?
          dynamic_cast<const MixtureOfATransitionModel*>(&pAbmtm->mixedModel())
        : dynamic_cast<const MixtureOfATransitionModel*>(mixmodel);

      // get rid of biblio link
      if (pAbmtm)
      {
        if (parname != "")
        {
          parname = pAbmtm->getPmodelParName(pAbmtm->getParameterNameWithoutNamespace(parname));
          // Remove distribution suffix
          auto pos = parname.rfind("_");
          if (pos != string::npos)
            parname = parname.substr(0, pos);
        }
      }

      // if no parname at that point, look for mixed parameter in model
      size_t nummod = pMatm->getNumberOfModels();
      if (realparname == "")
      {
        ParameterList pl = pMatm->model(0).getParameters();
        for (size_t i2 = 0; i2 < pl.size(); i2++)
        {
          string pl2n = pl[i2].getName();

          if (dynamic_cast<const ConstantDistribution*>(&pMatm->distribution(pl2n)) == NULL)
          {
            parname = pl2n;
            if (parname.size() > 0 && !pMatm->hasDistribution(parname))
              parname = parname.substr(0, parname.rfind("_"));

            if (!pMatm->hasDistribution(parname))
              parname = "";

            if (parname.size() > 0)
            {
              break;
            }
          }
        }
        if (parname != "")
          realparname = parname + "_" + TextTools::toString(modNum);
      }

      // there, all ways to get parname are exhausted
      if (parname == "")
      {
        ApplicationTools::displayError("Argument likelihoods.parameter_name is required.");
        exit(-1);
      }
      else
        ApplicationTools::displayResult("likelihoods.parameter_name", realparname);

      vector< Vuint > vvnmod;
      size_t i2 = 0;

      while (i2 < nummod)
      {
        string par2 = parname + "_" + TextTools::toString(i2 + 1);
        Vuint vnmod = pMatm->getSubmodelNumbers(par2);
        if (vnmod.size() != 0)
          vvnmod.push_back(vnmod);
        i2++;
      }

      size_t nbcl = vvnmod.size();
      if (nbcl <= 1)
        throw Exception("Parameter " + realparname + " is not mixed.");

      Vdouble vprob = pMatm->getProbabilities();

      // vectors of sets of probabilities for each value of parname
      vector<vector<double>> vvprob;

      // vector of total probabilities for each value of parname
      vector<double> vsprob;

      for (size_t i = 0; i < nbcl; i++)
      {
        vector<double> vprob2;
        for (size_t j = 0; j < vvnmod[i].size(); j++)
        {
          vprob2.push_back(vprob[static_cast<size_t>(vvnmod[i][j])]);
        }

        vvprob.push_back(vprob2);
        vsprob.push_back(VectorTools::sum(vvprob[i]));
      }

      vector<string> colNames;
      colNames.push_back("Sites");
      colNames.push_back("Ll");

      Vdouble dval;
      for (size_t i = 0; i < nbcl; i++)
      {
        auto pSM = pMatm->getNModel(static_cast<size_t>(vvnmod[i][0]));
        double valPar = pSM->getParameterValue(pSM->getParameterNameWithoutNamespace(parname));
        dval.push_back(valPar);
        colNames.push_back("Ll_" + realparname + "=" + TextTools::toString(valPar));
      }
      for (size_t i = 0; i < nbcl; i++)
      {
        colNames.push_back("Pr_" + realparname + "=" + TextTools::toString(dval[i]));
      }

      colNames.push_back("mean");

      shared_ptr<DataTable> rates = make_shared<DataTable>(nSites, colNames.size());
      rates->setColumnNames(colNames);

      // output sites
      for (size_t i = 0; i < nSites; i++)
      {
        const auto& currentSite = sites->site(i);
        int currentSitePosition = currentSite.getCoordinate();
        (*rates)(i, "Sites") = TextTools::toString(currentSitePosition);
      }

      // Likelihoods
      auto vl = tl->getLikelihoodPerSite();
      for (unsigned int j = 0; j < nSites; j++)
      {
        (*rates)(j, "Ll") = TextTools::toString(log(vl[j]));
      }


      vector<vector<DataLik>> vvd;
      vector<double> vRates = pMatm->getVRates();

      for (size_t i = 0; i < nbcl; ++i)
      {
        string par2 = realparname + "_" + TextTools::toString(i + 1);

        auto func = [nummod, i, &vvprob, &vvnmod, &vsprob](shared_ptr<BranchModelInterface> cmodel){
              auto pAbmtm2 = dynamic_pointer_cast<AbstractBiblioMixedTransitionModel>(cmodel);
              if (pAbmtm2)
              {
                for (unsigned int j = 0; j < nummod; ++j)
                {
                  pAbmtm2->setNProbability(j, 0);
                }

                for (size_t j = 0; j < vvprob[i].size(); ++j)
                {
                  pAbmtm2->setNProbability(static_cast<size_t>(vvnmod[i][j]), vsprob[i] > NumConstants::TINY() ? vvprob[i][j] / vsprob[i] : 1. / (double)vvprob[i].size());
                }
              }
              else
              {
                auto pMatm2 = dynamic_pointer_cast<MixtureOfATransitionModel>(cmodel);
                if (!pMatm2)
                  throw Exception("Not mixed model " + cmodel->getName());
                for (unsigned int j = 0; j < nummod; ++j)
                {
                  pMatm2->setNProbability(j, 0);
                }

                for (size_t j = 0; j < vvprob[i].size(); ++j)
                {
                  pMatm2->setNProbability(static_cast<size_t>(vvnmod[i][j]),  vsprob[i] > NumConstants::TINY() ? vvprob[i][j] / vsprob[i] : 1. / (double)vvprob[i].size());
                }
              }
            };

        model->modify(func, false);

        auto vd = tl->getLikelihoodPerSite();

        for (unsigned int j = 0; j < nSites; j++)
        {
          (*rates)(j, i + 2) = TextTools::toString(log(vd[j]));
        }

        vvd.push_back(vd);

        ApplicationTools::displayMessage("\n");
        ApplicationTools::displayMessage("Parameter " + par2 + "=" + TextTools::toString(dval[i]) + " with rate=" + TextTools::toString(vRates[i]));

        ApplicationTools::displayResult("Log likelihood", TextTools::toString(tl->getValue(), 15));
        ApplicationTools::displayResult("Probability", TextTools::toString(vsprob[i], 15));
      }


      for (size_t j = 0; j < nSites; j++)
      {
        Vdouble vd;
        for (size_t i = 0; i < nbcl; i++)
        {
          vd.push_back(std::log(vsprob[i]) + log(vvd[i][j]));
        }

        VectorTools::logNorm(vd);
        for (size_t i = 0; i < nbcl; i++)
        {
          (*rates)(j, nbcl + i + 2) = TextTools::toString(std::exp(vd[i]));
        }
        (*rates)(j, 2 * nbcl + 2) = TextTools::toString(VectorTools::sumExp(vd, dval));
      }

      DataTable::write(*rates, out, "\t");
    }

    ApplicationTools::displayMessage("\n");
    bppmixedlikelihoods.done();
  }

  catch (exception& e)
  {
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}
