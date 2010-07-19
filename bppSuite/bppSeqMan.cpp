//
// File: bppSeqMan.cpp
// Created by: Julien Dutheil
// Created on: Oct Tue 02 9:00 2007
//

/*
Copyright or © or Copr. CNRS

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

// From SeqLib:
#include <Seq/Alphabet.h>
#include <Seq/VectorSiteContainer.h>
#include <Seq/SequenceApplicationTools.h>
#include <Seq/ioseq>
#include <Seq/containers>
#include <Seq/SequenceTools.h>
#include <Seq/GeneticCode.h>
#include <Seq/StandardGeneticCode.h>
#include <Seq/VertebrateMitochondrialGeneticCode.h>
#include <Seq/InvertebrateMitochondrialGeneticCode.h>
#include <Seq/EchinodermMitochondrialGeneticCode.h>
#include <Seq/SiteTools.h>

// From Utils:
#include <Utils/BppApplication.h>
#include <Utils/ApplicationTools.h>
#include <Utils/FileTools.h>
#include <Utils/KeyvalTools.h>

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
  cout << "*           Bio++ Sequence Manipulator, version 0.3              *" << endl;
  cout << "* Author: J. Dutheil                        Last Modif. 15/04/10 *" << endl;
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
        Sequence* seq = new Sequence(old->getName(), old->getContent(), old->getComments(), alpha);
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
      if (!AlphabetTools::isNucleicAlphabet(sequences->getAlphabet()))
        throw Exception("Error in translation: alphabet is not of type 'nucleic'.");
      GeneticCode* gc = NULL;
      string gcstr = ApplicationTools::getStringParameter("code", cmdArgs, "Standard");
      if (gcstr == "Standard")
        gc = new StandardGeneticCode(dynamic_cast<const NucleicAlphabet *>(sequences->getAlphabet()));
      else if (gcstr == "VerMito")
        gc = new VertebrateMitochondrialGeneticCode(dynamic_cast<const NucleicAlphabet *>(sequences->getAlphabet()));
      else if (gcstr == "InvMito")
        gc = new InvertebrateMitochondrialGeneticCode(dynamic_cast<const NucleicAlphabet *>(sequences->getAlphabet()));
      else if (gcstr == "EchMito")
        gc = new EchinodermMitochondrialGeneticCode(dynamic_cast<const NucleicAlphabet *>(sequences->getAlphabet()));
      else throw Exception("Unknown genetic code: " + gcstr);

      OrderedSequenceContainer* sc = 0;
      if (aligned) sc = new VectorSiteContainer(sequences->getAlphabet());
      else         sc = reinterpret_cast<OrderedSequenceContainer*>(new VectorSequenceContainer(sequences->getAlphabet()));
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
        Sequence* seq = SequenceTools::removeGaps(sequences->getSequence(i));
        sc->addSequence(*seq);
        delete seq;
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
        Sequence* seq = new Sequence(sequences->getSequence(i));
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
        Sequence* seq = new Sequence(sequences->getSequence(i));
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
    if (cmdName == "RemoveStops")
    {
      OrderedSequenceContainer* sc = 0;
      if (aligned) sc = new VectorSiteContainer(sequences->getAlphabet());
      else         sc = reinterpret_cast<OrderedSequenceContainer*>(new VectorSequenceContainer(sequences->getAlphabet()));
      for (unsigned int i = 0; i < sequences->getNumberOfSequences(); i++)
      {
        Sequence seq = sequences->getSequence(i);
        SequenceTools::getCDS(seq, false, true, true, false);
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
        for (unsigned int i = sites->getNumberOfSites(); i > 0; i--)
        {
          map<int, double> freqs;
          SiteTools::getFrequencies(sites->getSite(i - 1), freqs);
          if (freqs[-1] >= gapFreq) sites->deleteSite(i - 1);
        }
      }
      else
      {
        unsigned int gapNum=TextTools::to<unsigned int>(maxGapOption);
        for (unsigned int i = sites->getNumberOfSites(); i > 0; i--)
        {
          map<int, unsigned int> counts;
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

