//
// File: bppTools.h
// Created by: Laurent Guéguen
// Created on: mardi 28 novembre 2017, à 09h 05
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

#include <Bpp/App/ApplicationTools.h>

#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Alphabet/CodonAlphabet.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/Container/AlignedValuesContainer.h>

#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>
#include <Bpp/Phyl/Tree/PhyloTree.h>
#include <Bpp/Phyl/NewLikelihood/SubstitutionProcessCollection.h>
#include <Bpp/Phyl/NewLikelihood/SequenceEvolution.h>
#include <Bpp/Phyl/NewLikelihood/PhyloLikelihoods/PhyloLikelihoodContainer.h>

#include <Bpp/Phyl/NewLikelihood/DataFlow/DataFlow.h>

#ifndef TOOLS_H
#define TOOLS_H

/******************************************************************************/

namespace bpp
{
  
  class bppTools
  {
  public:
  
    /*
     * @brief display the basic help
     *
     */

    static void help(const std::string& program);

    /***************************************
     * @brief Methods to build objects
     *
     * @{
     *
     */
    
    
    /*
     * @brief get the Alphabet
     *
     */
  
    static Alphabet* getAlphabet(const std::map<std::string, std::string>& params);
  
    /*
     * @brief get the GeneticCode
     *
     */
  
    static GeneticCode* getGeneticCode(const std::map<std::string, std::string>& params,
                                       const Alphabet* alphabet);
  
    /*
     *
     * @brief Get the std::map of alignments
     *
     */
    
    static std::map<size_t, AlignedValuesContainer*> getAlignmentsMap(const std::map<std::string, std::string>& params,
                                                                      const Alphabet* alphabet,
                                                                      bool changeGapsToUnknownCharacters = true,
                                                                      bool optionalData = false);
  
    
    /*
     * @brief Get the std::map of initial phylo trees
     *
     */
  
    static std::map<size_t, std::shared_ptr<PhyloTree>> getPhyloTreesMap(const std::map<std::string, std::string>& params,
                                                                         const std::map<size_t, AlignedValuesContainer*>& mSites,
                                                                         std::map<std::string, std::string>& unparsedparams);

    /*
     * @brief get the collection of objects necessary to build
     * substitution processes.
     *
     */
    
    static SubstitutionProcessCollection* getCollection(const std::map<std::string, std::string>& params,
                                                        const Alphabet* alphabet,
                                                        const GeneticCode* gCode,
                                                        const std::map<size_t, AlignedValuesContainer*>& mSites,
                                                        const std::map<size_t, std::shared_ptr<PhyloTree>>& mpTree,
                                                        std::map<std::string, std::string>& unparsedparams);

    static SubstitutionProcessCollection* getCollection(const std::map<std::string, std::string>& params,
                                                        const Alphabet* alphabet,
                                                        const GeneticCode* gCode,
                                                        const std::map<size_t, AlignedValuesContainer*>& mSites,
                                                        std::map<std::string, std::string>& unparsedparams);

    /*
     * @brief get the substitution processes.
     *
     */
    
    static std::map<size_t, SequenceEvolution*> getProcesses(const std::map<std::string, std::string>& params,
                                                             SubstitutionProcessCollection& collection,
                                                             std::map<std::string, std::string>& unparsedparams);
    
    /*
     * @brief get the phylolikelihoods.
     *
     */

    static PhyloLikelihoodContainer* getPhyloLikelihoods(const std::map<std::string, std::string>& params,
                                                         Context& context,
                                                         std::map<size_t, SequenceEvolution*> mSeqEvol, 
                                                         SubstitutionProcessCollection& collection,
                                                         const std::map<size_t, AlignedValuesContainer*>& mSites);


    /*
     * @brief return only final phylo likelihood
     *
     */
    
    static PhyloLikelihood* getResultPhyloLikelihood(const std::map<std::string, std::string>& params,
                                                     Context& context,
                                                     const Alphabet* alphabet,
                                                     const GeneticCode* gCode,
                                                     std::map<std::string, std::string>& unparsedparams);

    /*
     * @}
     *
     */

    /*
     * @brief Method to have a clean likelihood (ie not saturated, nor infinite).
     *
     */
    
    static void fixLikelihood(const std::map<std::string, std::string>& params,
                              const Alphabet* alphabet,
                              const GeneticCode* gCode,
                              PhyloLikelihood* phylolik);
    
    /*
     * @brief Display parameter values.
     * @param tl the phylolikelihood
     * @param displaylL if log-likelihood is displayed (default: true)
     */

    static void displayParameters(const PhyloLikelihood& tl, bool displaylL = true);
    
  };
} // end of namespace


#endif

