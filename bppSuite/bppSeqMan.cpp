//
// File: bppSeqMan.cpp
// Created by: Julien Dutheil
// Created on: Oct Tue 02 9:00 2007
//

/*
Copyright or © or Copr. Bio++ Development Team

This software is a computer program whose purpose is to simulate sequence
data according to a phylogenetic tree and an evolutionary model.

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
#include <fstream>
#include <iomanip>

using namespace std;

#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Io/FileTools.h>
#include <Bpp/Text/KeyvalTools.h>

// From SeqLib:
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/CodonSiteTools.h>
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/Io.all>
#include <Bpp/Seq/Container.all>
#include <Bpp/Seq/SequenceTools.h>
#include <Bpp/Seq/GeneticCode.all>

//From PhylLib:
#include <Bpp/Phyl/Tree/Tree.h>
#include <Bpp/Phyl/App/PhylogeneticsApplicationTools.h>

using namespace bpp;

void help()
{
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
  (*ApplicationTools::message << "bppseqman parameter1_name=parameter1_value").endLine();
  (*ApplicationTools::message << "      parameter2_name=parameter2_value ... param=option_file").endLine();
  (*ApplicationTools::message).endLine();
  (*ApplicationTools::message << "  Refer to the Bio++ Program Suite Manual for a list of available options.").endLine();
  (*ApplicationTools::message << "__________________________________________________________________________").endLine();
}

int main(int args, char** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*           Bio++ Sequence Manipulator, version 0.6              *" << endl;
  cout << "* Author: J. Dutheil                        Last Modif. 21/12/11 *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;
  
  if (args == 1)
  {
    help();
    return 0;
  }
  
  try {

  BppApplication bppseqman(args, argv, "BppSeqMan");
  bppseqman.startTimer();
  
  // Get alphabet
  Alphabet* alphabet = SequenceApplicationTools::getAlphabet(bppseqman.getParams(), "", false, true, true);
  auto_ptr<GeneticCode> gCode;
  CodonAlphabet* codonAlphabet = dynamic_cast<CodonAlphabet*>(alphabet);
  if (codonAlphabet) {
    string codeDesc = ApplicationTools::getStringParameter("genetic_code", bppseqman.getParams(), "Standard", "", true, true);
    gCode.reset(SequenceApplicationTools::getGeneticCode(codonAlphabet->getNucleicAlphabet(), codeDesc));
  }

  // Get sequences:
  SequenceContainer* tmp = SequenceApplicationTools::getSequenceContainer(alphabet, bppseqman.getParams(), "", true, true);
  OrderedSequenceContainer* sequences = new VectorSequenceContainer(*tmp);
  delete tmp;
  ApplicationTools::displayResult("Number of sequences", sequences->getNumberOfSequences());
  
  // Perform manipulations
  
  vector<string> actions = ApplicationTools::getVectorParameter<string>("sequence.manip", bppseqman.getParams(), ',', "", "", false, false);
  
  bool aligned = false;

  for (unsigned int a = 0; a < actions.size(); a++)
  {
    string cmdName;
    map<string, string> cmdArgs;
    KeyvalTools::parseProcedure(actions[a], cmdName, cmdArgs);
    ApplicationTools::displayResult("Performing action", cmdName);

    // +-----------------+
    // | Complementation |
    // +-----------------+
    if (cmdName == "Complement")
    {
      OrderedSequenceContainer* sc = 0;
      if (aligned) sc = new VectorSiteContainer(sequences->getAlphabet());
      else         sc = reinterpret_cast<OrderedSequenceContainer*>(new VectorSequenceContainer(sequences->getAlphabet()));
      for (unsigned int i = 0; i < sequences->getNumberOfSequences(); i++)
      {
        Sequence* seq = SequenceTools::getComplement(sequences->getSequence(i));
        sc->addSequence(*seq, false);
        delete seq;
      }
      delete sequences;
      sequences = sc;
    }
    // +------------------------+
    // | (Reverse)Transcription |
    // +------------------------+
    else if (cmdName == "Transcript")
    {
      if (sequences->getAlphabet()->getAlphabetType() == AlphabetTools::DNA_ALPHABET.getAlphabetType())
      {
        OrderedSequenceContainer* sc = 0;
        if (aligned) sc = new VectorSiteContainer(&AlphabetTools::RNA_ALPHABET);
        else         sc = reinterpret_cast<OrderedSequenceContainer*>(new VectorSequenceContainer(&AlphabetTools::RNA_ALPHABET));
        for (unsigned int i = 0; i < sequences->getNumberOfSequences(); i++)
        {
          Sequence* seq = SequenceTools::transcript(sequences->getSequence(i));
          sc->addSequence(*seq, false);
          delete seq;
        }
        delete sequences;
        sequences = sc;
      }
      else if (sequences->getAlphabet()->getAlphabetType() == AlphabetTools::RNA_ALPHABET.getAlphabetType())
      {
        OrderedSequenceContainer* sc = 0;
        if (aligned) sc = new VectorSiteContainer(&AlphabetTools::DNA_ALPHABET);
        else         sc = reinterpret_cast<OrderedSequenceContainer*>(new VectorSequenceContainer(&AlphabetTools::DNA_ALPHABET));
        for (unsigned int i = 0; i < sequences->getNumberOfSequences(); i++)
        {
          Sequence* seq = SequenceTools::reverseTranscript(sequences->getSequence(i));
          sc->addSequence(*seq, false);
          delete seq;
        }
        delete sequences;
        sequences = sc;
      }
      else throw Exception("Transcription error: input alphabet must be of type 'nucleic'.");
    }
    // +-------------------------------+
    // | Switching nucleotide alphabet |
    // +-------------------------------+
    else if (cmdName == "Switch")
    {
      const Alphabet* alpha = 0;
      if (sequences->getAlphabet()->getAlphabetType() == AlphabetTools::DNA_ALPHABET.getAlphabetType())
      {
        alpha = &AlphabetTools::RNA_ALPHABET;
      }
      else if (sequences->getAlphabet()->getAlphabetType() == AlphabetTools::RNA_ALPHABET.getAlphabetType())
      {
        alpha = &AlphabetTools::DNA_ALPHABET;
      }
      else throw Exception("Cannot switch alphabet type, alphabet is not of type 'nucleic'.");
      OrderedSequenceContainer* sc = 0;
      if (aligned) sc = new VectorSiteContainer(alpha);
      else         sc = reinterpret_cast<OrderedSequenceContainer*>(new VectorSequenceContainer(alpha));
      for (unsigned int i = 0; i < sequences->getNumberOfSequences(); i++)
      {
        const Sequence* old = &sequences->getSequence(i);
        Sequence* seq = new BasicSequence(old->getName(), old->getContent(), old->getComments(), alpha);
        sc->addSequence(*seq, false);
        delete seq;
      }
      delete sequences;
      sequences = sc;
    }
    // +-------------+
    // | Translation |
    // +-------------+
    else if (cmdName == "Translate")
    {
      if (!AlphabetTools::isCodonAlphabet(sequences->getAlphabet()))
        throw Exception("Error in translation: alphabet is not of type 'codon'.");
      GeneticCode* gc = NULL;
      string gcstr = ApplicationTools::getStringParameter("code", cmdArgs, "Standard");
      gc = SequenceApplicationTools::getGeneticCode(dynamic_cast<const CodonAlphabet*>(sequences->getAlphabet())->getNucleicAlphabet(), gcstr);

      OrderedSequenceContainer* sc = 0;
      if (aligned) sc = new VectorSiteContainer(&AlphabetTools::PROTEIN_ALPHABET);
      else         sc = reinterpret_cast<OrderedSequenceContainer*>(new VectorSequenceContainer(&AlphabetTools::PROTEIN_ALPHABET));
      for (unsigned int i = 0; i < sequences->getNumberOfSequences(); i++)
      {
        Sequence* seq = gc->translate(sequences->getSequence(i));
        sc->addSequence(*seq, false);
        delete seq;
      }
      delete sequences;
      sequences = sc;      
    }
    // +-------------+
    // | Remove gaps |
    // +-------------+
    else if (cmdName == "RemoveGaps")
    {
      VectorSequenceContainer* sc = new VectorSequenceContainer(sequences->getAlphabet());
      for (unsigned int i = 0; i < sequences->getNumberOfSequences(); i++)
      {
        auto_ptr<Sequence> seq(sequences->getSequence(i).clone());
        SequenceTools::removeGaps(*seq);
        sc->addSequence(*seq);
      }
      delete sequences;
      sequences = sc;
      aligned = false;
    }
    // +---------------------------+
    // | Change gaps to unresolved |
    // +---------------------------+
    else if (cmdName == "GapToUnknown")
    {
      OrderedSequenceContainer* sc = 0;
      if (aligned) sc = new VectorSiteContainer(sequences->getAlphabet());
      else         sc = reinterpret_cast<OrderedSequenceContainer*>(new VectorSequenceContainer(sequences->getAlphabet()));
      for (unsigned int i = 0; i < sequences->getNumberOfSequences(); i++)
      {
        Sequence* seq = new BasicSequence(sequences->getSequence(i));
        SymbolListTools::changeGapsToUnknownCharacters(*seq);
        sc->addSequence(*seq, false);
        delete seq;
      }
      delete sequences;
      sequences = sc;
    }
    // +---------------------------+
    // | Change unresolved to gaps |
    // +---------------------------+
    else if (cmdName == "UnknownToGap")
    {
      OrderedSequenceContainer* sc = 0;
      if (aligned) sc = new VectorSiteContainer(sequences->getAlphabet());
      else         sc = new VectorSequenceContainer(sequences->getAlphabet());
      for (unsigned int i = 0; i < sequences->getNumberOfSequences(); i++)
      {
        Sequence* seq = new BasicSequence(sequences->getSequence(i));
        SymbolListTools::changeUnresolvedCharactersToGaps(*seq);
        sc->addSequence(*seq, false);
        delete seq;
      }
      delete sequences;
      sequences = sc;
    }
    
    // +--------------+
    // | Remove stops |
    // +--------------+
    else if (cmdName == "RemoveStops")
    {
      SiteContainer* sites = dynamic_cast<SiteContainer*>(sequences);
      if (!sites)
      {
        VectorSequenceContainer* sc = new VectorSequenceContainer(sequences->getAlphabet());
        for (size_t i = 0; i < sequences->getNumberOfSequences(); ++i)
        {
          auto_ptr<Sequence> seq(sequences->getSequence(i).clone());
          SequenceTools::removeStops(*seq, *gCode);
          sc->addSequence(*seq);
        }
        delete sequences;
        sequences = sc;
      } else {
        VectorSiteContainer* sc = new VectorSiteContainer(sequences->getAlphabet());
        for (size_t i = 0; i < sequences->getNumberOfSequences(); ++i)
        {
          auto_ptr<Sequence> seq(sequences->getSequence(i).clone());
          SequenceTools::replaceStopsWithGaps(*seq, *gCode);
          sc->addSequence(*seq);
        }
        delete sequences;
        sequences = sc;
      }
    }

    // +--------------+
    // | Remove stops |
    // +--------------+
    else if (cmdName == "RemoveColumnsWithStop")
    {
      SiteContainer* sites = dynamic_cast<SiteContainer*>(sequences);
      if (!sites)
      {
        throw Exception("'RemoveColumnsWithStop' can only be used on alignment. You may consider using the 'CoerceToAlignment' command.");
      }

      for (size_t i = sites->getNumberOfSites(); i > 0; i--)
      {
        if (CodonSiteTools::hasStop(sites->getSite(i-1), *gCode))
          sites->deleteSite(i - 1);
      }
    }

    // +---------+
    // | Get CDS |
    // +---------+
    else if (cmdName == "GetCDS")
    {
      OrderedSequenceContainer* sc = 0;
      if (aligned) sc = new VectorSiteContainer(sequences->getAlphabet());
      else         sc = reinterpret_cast<OrderedSequenceContainer*>(new VectorSequenceContainer(sequences->getAlphabet()));
      for (unsigned int i = 0; i < sequences->getNumberOfSequences(); i++)
      {
        BasicSequence seq = sequences->getSequence(i);
        size_t len = seq.size();
        SequenceTools::getCDS(seq, *gCode, false, true, true, false);
        if (aligned) {
          for (size_t c = seq.size(); c < len; ++c)
            seq.addElement(seq.getAlphabet()->getGapCharacterCode());
        }
        sc->addSequence(seq, false);
      }
      delete sequences;
      sequences = sc;
    }

    // +--------------------------+
    // | Resolve dotted alignment |
    // +--------------------------+
    else if (actions[a] == "CoerceToAlignment")
    {
      SiteContainer* sites = dynamic_cast<SiteContainer*>(sequences);
      if(! sites)
      {
        sites = new VectorSiteContainer(*sequences);
        delete sequences;
        sequences = sites;
      }
      aligned = true;
    }
    else if (actions[a] == "ResolvedDotted")
    {
      SiteContainer* sites = dynamic_cast<SiteContainer *>(sequences);
      if (!sites)
      {
        throw Exception("'ResolvedDotted' can only be used on alignment. You may consider using the 'CoerceToAlignment' command.");
      }

      const Alphabet* alpha = 0;
      string alphastr = ApplicationTools::getStringParameter("alphabet", cmdArgs, "DNA");
      if (alphastr == "DNA") alpha = &AlphabetTools::DNA_ALPHABET;
      else if (alphastr == "RNA") alpha = &AlphabetTools::RNA_ALPHABET;
      else if (alphastr == "Protein") alpha = &AlphabetTools::PROTEIN_ALPHABET;
      else throw Exception("Resolved alphabet must be one of [DNA|RNA|Protein] for solving dotted alignment.");
      OrderedSequenceContainer* resolvedCont = SiteContainerTools::resolveDottedAlignment(*sites, alpha);
      delete sequences;
      sequences = resolvedCont;
    }
    // +---------------------+
    // | Keep complete sites |
    // +---------------------+
    else if (cmdName == "KeepComplete")
    {
      SiteContainer* sites = dynamic_cast<SiteContainer *>(sequences);
      if (!sites)
      {
        throw Exception("'KeepComplete' can only be used on alignment. You may consider using the 'CoerceToAlignment' command.");
      }

      string maxGapOption = ApplicationTools::getStringParameter("maxGapAllowed", cmdArgs, "100%");
      if (maxGapOption[maxGapOption.size()-1] == '%')
      {
        double gapFreq = TextTools::toDouble(maxGapOption.substr(0, maxGapOption.size()-1)) / 100.;
        for (size_t i = sites->getNumberOfSites(); i > 0; i--)
        {
          map<int, double> freqs;
          SiteTools::getFrequencies(sites->getSite(i - 1), freqs);
          if (freqs[-1] > gapFreq) sites->deleteSite(i - 1);
        }
      }
      else
      {
        size_t gapNum = TextTools::to<size_t>(maxGapOption);
        for (size_t i = sites->getNumberOfSites(); i > 0; i--)
        {
          map<int, size_t> counts;
          SiteTools::getCounts(sites->getSite(i - 1), counts);
          counts[-1]; //Needed in case this entry does not exist in the map. This will set it to 0.
          if (counts[-1] > gapNum) sites->deleteSite(i-1);
        }
      }
    }
    // +-----------------+
    // | Invert sequence |
    // +-----------------+
    else if (cmdName == "Invert")
    {
      OrderedSequenceContainer* sc = 0;
      if (aligned) sc = new VectorSiteContainer(sequences->getAlphabet());
      else         sc = reinterpret_cast<OrderedSequenceContainer*>(new VectorSequenceContainer(sequences->getAlphabet()));
      for (unsigned int i = 0; i < sequences->getNumberOfSequences(); i++)
      {
        const Sequence* old = &sequences->getSequence(i);
        Sequence* seq = SequenceTools::getInvert(*old);
        sc->addSequence(*seq, false);
        delete seq;
      }
      delete sequences;
      sequences = sc;
    }
    // +------------------+
    // | GetCodonPosition |
    // +------------------+
    else if (cmdName == "GetCodonPosition")
    {
      unsigned int pos = ApplicationTools::getParameter<unsigned int>("position", cmdArgs, 3);
      OrderedSequenceContainer* sc = dynamic_cast<OrderedSequenceContainer*>(SequenceContainerTools::getCodonPosition(*sequences, pos - 1));
      delete sequences;
      if (aligned) {
        sequences = new VectorSiteContainer(*sc);
        delete sc;
      } else {
        sequences = sc;
      }
    }
    // +-----------------+
    // | FilterFromTree |
    // +-----------------+
    else if (cmdName == "FilterFromTree")
    {
      auto_ptr<Tree> tree(PhylogeneticsApplicationTools::getTree(cmdArgs, ""));
      vector<string> names = tree->getLeavesNames();
      OrderedSequenceContainer* reorderedSequences = 0;
      if (aligned) {
        reorderedSequences = new VectorSiteContainer(sequences->getAlphabet());
      } else {
        reorderedSequences = new VectorSequenceContainer(sequences->getAlphabet());
      }
      for (size_t i = 0; i < names.size(); ++i) {
        reorderedSequences->addSequence(sequences->getSequence(names[i]), false);
      }
      delete sequences;
      sequences = reorderedSequences;
    }

    else throw Exception("Unknown action: " + cmdName);
  }
  
  // Write sequences
  ApplicationTools::displayBooleanResult("Final sequences are aligned", aligned);
  if (aligned)
  {
    SequenceApplicationTools::writeAlignmentFile(*dynamic_cast<SiteContainer*>(sequences), bppseqman.getParams(), "", true);
  }
  else
  {
    SequenceApplicationTools::writeSequenceFile(*sequences, bppseqman.getParams(), "", true);
  }

  delete alphabet;
  delete sequences;

  bppseqman.done();

  } catch(exception & e) {
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}

