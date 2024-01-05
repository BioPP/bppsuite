//
// File: bppBranchLik.cpp
// Created by: Laurent Guéguen
// Created on: jeudi 28 septembre 2023, à 09h 11
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
#include <Bpp/Text/KeyvalTools.h>
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
#include <Bpp/Phyl/Likelihood/DataFlow/LikelihoodCalculationOnABranch.h>
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/OnABranchPhyloLikelihood.h>
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/PartitionProcessPhyloLikelihood.h>
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/SingleProcessPhyloLikelihood.h>
#include <Bpp/Phyl/Likelihood/PhyloLikelihoods/OneProcessSequencePhyloLikelihood.h>
#include <Bpp/Phyl/Model/AbstractBiblioMixedTransitionModel.h>
#include <Bpp/Phyl/Model/MixedTransitionModel.h>
#include <Bpp/Phyl/Model/MixtureOfATransitionModel.h>
#include <Bpp/Phyl/Model/MixtureOfTransitionModels.h>
#include <Bpp/Phyl/Model/RateDistribution/ConstantRateDistribution.h>

using namespace bpp;

/******************************************************************************/

struct rangeLik
{
  rangeLik() :
    lik_(0),
    range_()
  {};
  
  std::shared_ptr<LikelihoodCalculationSingleProcess> lik_;   // LikelihoodCalculationSingleProcess used
  std::vector<size_t> range_;    // vector of the positions for this process
};

struct rangeProc
{
  rangeProc() :
    proc_(0),
    range_()
  {};
  
  std::shared_ptr<const SubstitutionProcessInterface> proc_;   // LikelihoodCalculationSingleProcess used
  std::vector<size_t> range_;    // vector of the positions for this process
};

int main(int args, char** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*     Bio++ Computation of likelihoods on branches   *" << endl;
  cout << "*                        Version " << BPP_VERSION << ".                          *" << endl;
  cout << "* Author: L. Guéguen                       Last Modif.: " << BPP_REL_DATE << " *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  try
  {
    BppPhylogeneticsApplication bppbranchlik(args, argv, "BppBranchLik");

    if (args == 1)
    {
      bppbranchlik.help("bppbranchlik");
      return 0;
    }

    bppbranchlik.startTimer();

    Context context;
    
    ///// Alphabet

    shared_ptr<const Alphabet> alphabet(bppbranchlik.getAlphabet());

    /// GeneticCode
    
    shared_ptr<const GeneticCode> gCode(bppbranchlik.getGeneticCode(alphabet));

    // get the data

    auto mSitesuniq = bppbranchlik.getConstAlignmentsMap(alphabet, true);

    const std::map<size_t, std::shared_ptr<const AlignmentDataInterface > > mSites = PhylogeneticsApplicationTools::uniqueToSharedMap<const TemplateAlignmentDataInterface<string>>(mSitesuniq);

    if (mSites.size()!=1)
      throw Exception("Only one alignment possible.");

    auto sites = mSites.begin()->second;
    const size_t nSites = sites->getNumberOfSites();
        
    /////// Get the map of initial trees

    map<string, string> unparsedParams;

    auto mpTree = bppbranchlik.getPhyloTreesMap(mSites, unparsedParams);

    /////////////////
    // Computing stuff


    shared_ptr<SubstitutionProcessCollection> SPC(bppbranchlik.getCollection(alphabet, gCode, mSites, mpTree, unparsedParams));
    
    auto mSeqEvoltmp = bppbranchlik.getProcesses(SPC, unparsedParams);
    
    auto mSeqEvol = PhylogeneticsApplicationTools::uniqueToSharedMap<SequenceEvolution>(mSeqEvoltmp);

    auto mPhyl(bppbranchlik.getPhyloLikelihoods(context, mSeqEvol, SPC, mSites));

    if (!mPhyl->hasPhyloLikelihood(0))
      throw Exception("Missing phyloLikelihoods.");

    auto tl=dynamic_pointer_cast<AlignedPhyloLikelihoodInterface>((*mPhyl)[0]);

    if (tl==0)
      throw Exception("Only possible on aligned phyloLikelihood.");
        
    //Check initial likelihood:
      
    bppbranchlik.fixLikelihood(alphabet, gCode, tl);

    
    // /////////////////////////////////////////////
    // Get the alternate process

    ApplicationTools::displayMessage("\n");
    auto altprocnum = ApplicationTools::getParameter<size_t>("alt_process", bppbranchlik.getParams(), 1, "", true, false);
    ApplicationTools::displayResult("Alternative process", altprocnum);

    std::map<size_t, rangeProc> maltproc;
    
    if (SPC->hasSubstitutionProcessNumber(altprocnum))
    {
      rangeProc rl;
      rl.proc_=SPC->getSubstitutionProcess(altprocnum);
      rl.range_=std::vector<size_t>(nSites);
      std::iota(rl.range_.begin(), rl.range_.end(), 0);
      maltproc[1]=rl;
    }
    else
    {
      if (mSeqEvol.find(altprocnum) == mSeqEvol.end())
        throw BadIntegerException("Unknown alt_process number", static_cast<int>(altprocnum));

      auto seqev = mSeqEvol[altprocnum];

      auto partev = dynamic_pointer_cast<PartitionSequenceEvolution>(seqev);
      if (!partev)
        throw BadIntegerException("Unfit alt_process number: ask developers.", static_cast<int>(altprocnum));

      const auto& mapproc = partev->mapOfProcessSites();

      for (auto& proc : mapproc)
      {
        rangeProc rl;
        rl.proc_ = partev->getSubstitutionProcess(proc.first);
        rl.range_ = proc.second;
        maltproc[proc.first]=rl;
      }
    }
    
    // Up to now only for simple Seq Evolution
    

    // get the branch numbers of tree for alt process: all process must have the same numbers
    std::vector<uint> brnum;
    for (auto prn:maltproc)
    {
      const auto& brn = prn.second.proc_->parametrizablePhyloTree().getAllEdgesIndexes();
      if (brnum.size() == 0)
        brnum = brn;
      else
        if (brn != brnum)
          throw BadIntegerException("Non compatible trees in alternative process ", static_cast<int>(altprocnum));
    }



    // Set up segmentation of data following likelihoods
    std::map<size_t, rangeLik> mlik; // Vector of which likelihood calculation of the main phylo is assessed to each site
    
    auto calc = dynamic_pointer_cast<LikelihoodCalculationSingleProcess>(tl->getAlignedLikelihoodCalculation());

    if (calc)
    {
      rangeLik rl;
      rl.lik_=calc;
      rl.range_=std::vector<size_t>(nSites);
      std::iota(rl.range_.begin(), rl.range_.end(), 0);
      mlik[1]=rl;
    }
    else
    {
      auto parphyl = dynamic_pointer_cast<PartitionProcessPhyloLikelihood>(tl);
      if (!parphyl)
        throw Exception("Up to now, only available non Single Process is Partition. Ask developpers.");

      const auto& mapproc = parphyl->getPartitionSequenceEvolution()->mapOfProcessSites();

      for (auto proc : mapproc)
      {
        rangeLik rl;
        auto pn = parphyl->getPhyloLikelihood(proc.first);
        auto opsp = dynamic_pointer_cast<OneProcessSequencePhyloLikelihood>(pn);
        auto sp = dynamic_pointer_cast<SingleProcessPhyloLikelihood>(pn);
        if (opsp)
          rl.lik_ = opsp->getLikelihoodCalculationSingleProcess();
        else if (sp)
          rl.lik_ = sp->getLikelihoodCalculationSingleProcess();
        else
          throw Exception("LikelihoodSingleProcess not available for phylolikelihood " + TextTools::toString(proc.first));
        rl.range_ = proc.second;
        mlik[proc.first]=rl;
      }
    }
    
    // Output Format

    string outputDesc = ApplicationTools::getStringParameter("output.lik", bppbranchlik.getParams(), "PerBranch(file=outlik.txt)");

    string outputType;
    map<string, string> outputArgs;
    KeyvalTools::parseProcedure(outputDesc, outputType, outputArgs);

    bool perSite=(outputType.find("Site")!=string::npos);
    
    string outputFile = ApplicationTools::getStringParameter("file", outputArgs, "", "", true, 1);
    ofstream out(outputFile.c_str(), ios::out);
    out << std::setprecision(12);

    ApplicationTools::displayResult("Type of output", perSite?"perSitePerBranch":"PerBranch");
    ApplicationTools::displayResult("Output file", outputFile);

    // DataTable for output
    vector<string> colNames;
    shared_ptr<DataTable> rates;

    if (perSite){
      colNames.push_back("Sites");

      for (auto edgeid:brnum)
        colNames.push_back(TextTools::toString(edgeid));

      rates = make_shared<DataTable>(nSites, colNames.size());
      rates->setColumnNames(colNames);

      // output sites
      for (unsigned int i = 0; i < nSites; i++)
      {
        const auto& currentSite = sites->site(i);
        int currentSitePosition = currentSite.getCoordinate();
        (*rates)(i, "Sites") = string(TextTools::toString(currentSitePosition));
      }
    }
    else
    {
      out << "Branch\tLlik\n";
    }

    // Now branch wise computation
    
    for (auto edgeid:brnum)
    {
      double value = 0;
      for (const auto& plik:mlik)
      {
        auto lik = plik.second.lik_;
        const auto& rangelik = plik.second.range_;

        OnABranchPhyloLikelihood obp(context, lik, edgeid);

        for (const auto& pproc:maltproc)
        {
          const auto& rangealt=pproc.second.range_;

          if ((rangealt.back()<rangelik[0]) || (rangealt[0]>rangelik.back()))
            continue;

          auto altmodel = std::unique_ptr<BranchModelInterface>(pproc.second.proc_->model(edgeid,0).clone());
          auto confmodel = ConfiguredParametrizable::createConfigured<BranchModelInterface, ConfiguredModel>(context, std::move(altmodel));
          
          obp.setModel(confmodel);
        
          const auto& vl = obp.getLikelihoodPerSite();
          
          size_t poslik=0;
          while (rangelik[poslik]<rangealt[0])
            poslik++;
          
          size_t posalt=0;
          while (rangealt[posalt]<rangelik[0])
            posalt++;
          
          while (poslik<rangelik.size() && posalt<rangealt.size())
          {
            if (rangealt[posalt]<rangelik[poslik])
            {
              posalt++;
              continue;
            }
            if (rangelik[poslik]<rangealt[posalt])
            {
              poslik++;
              continue;
            }
            
            if (perSite)
              (*rates)(rangealt[posalt], TextTools::toString(edgeid)) = TextTools::toString(log(vl[poslik]));
            else
              value += -log(vl[poslik]);

            posalt++;
            poslik++;
          }
        }
      }
      
      if (!perSite)
      {
        out << edgeid << "\t" << value << "\n";
      }
    }
    
    if (perSite)
    {
      DataTable::write(*rates, out, "\t");
    }
    
    ApplicationTools::displayMessage("\n");
    bppbranchlik.done();
  }
  
  catch (exception& e)
  {
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}

