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
#include <Utils/AttributesTools.h>
#include <Utils/FileTools.h>
#include <Utils/ApplicationTools.h>

using namespace bpp;

void help()
{
  *ApplicationTools::message << "__________________________________________________________________________" << endl;
  *ApplicationTools::message << "       bppseqman arg1=value1 arg2=value2 ..." << endl;
  *ApplicationTools::message << "and/or bppseqman params=argfile" << endl;
  *ApplicationTools::message << "______________________________|___________________________________________" << endl;
}

int main(int args, char ** argv)
{
  cout << "******************************************************************" << endl;
  cout << "*           Bio++ Sequence Manipulator, version 0.1              *" << endl;
  cout << "* Author: J. Dutheil                        Last Modif. 02/03/09 *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;
  
  if(args == 1)
  {
    help();
    exit(0);
  }
  
  try {

  cout << "Parsing options:" << endl;
  
  map<string, string> params = AttributesTools::parseOptions(args, argv);

  // Get alphabet
  Alphabet * alphabet = SequenceApplicationTools::getAlphabet(params, "", false);

  // Get sequences:
  string sequenceFilePath = ApplicationTools::getAFilePath("input.sequence.file", params, true, true);
  string sequenceFormat = ApplicationTools::getStringParameter("input.sequence.format", params, "Fasta");
  ApplicationTools::displayResult("Input sequence format", sequenceFormat);
  ISequence * seqReader = NULL;
  if(sequenceFormat == "Mase")
  {
    seqReader = new Mase();
  }
  else if(sequenceFormat == "Phylip")
  {
    bool sequential = true, extended = true;
    if(params.find("input.sequence.format_phylip.order") != params.end())
    {
           if(params["input.sequence.format_phylip.order"] == "sequential" ) sequential = true;
      else if(params["input.sequence.format_phylip.order"] == "interleaved") sequential = false;
      else ApplicationTools::displayWarning("Argument '" +
             params["input.sequence.format_phylip.order"] +
             "' for parameter 'input.sequence.format_phylip.order' is unknown. " +
             "Default used instead: sequential.");
    }
    else ApplicationTools::displayWarning("Argument 'input.sequence.format_phylip.order' not found. Default used instead: sequential.");
    if(params.find("input.sequence.format_phylip.ext") != params.end())
    {
           if(params["input.sequence.format_phylip.ext"] == "extended") extended = true;
      else if(params["input.sequence.format_phylip.ext"] == "classic" ) extended = false;
      else ApplicationTools::displayWarning("Argument '" +
             params["input.sequence.format_phylip.ext"] +
             "' for parameter 'input.sequence.format_phylip.ext' is unknown. " +
             "Default used instead: extended.");
    }
    else ApplicationTools::displayWarning("Argument 'input.sequence.format_phylip.ext' not found. Default used instead: extended.");
    seqReader = new Phylip(extended, sequential);
  }
  else if(sequenceFormat == "Fasta") seqReader = new Fasta();
  else if(sequenceFormat == "Clustal") seqReader = new Clustal();
  else if(sequenceFormat == "GenBank") seqReader = new GenBank();
  else
  {
    ApplicationTools::displayError("Unknown sequence format.");
    exit(-1);
  }
  ApplicationTools::displayResult("Input sequence file", sequenceFilePath);
  SequenceContainer * tmp = seqReader->read(sequenceFilePath, alphabet);
  OrderedSequenceContainer * sequences = new VectorSequenceContainer(*tmp);
  delete tmp;
  delete seqReader;
  
  // Perform manipulations
  
  vector<string> actions = ApplicationTools::getVectorParameter<string>("sequence.manip", params, ',', "", "", false, false);
  
  for(unsigned int i = 0; i < actions.size(); i++)
  {
    ApplicationTools::displayResult("Performing action", actions[i]);

    // +-----------------+
    // | Complementation |
    // +-----------------+
    if(actions[i] == "complement")
    {
      VectorSequenceContainer * vsc = new VectorSequenceContainer(sequences->getAlphabet());
      for(unsigned int i = 0; i < sequences->getNumberOfSequences(); i++)
      {
        Sequence * tmp = SequenceTools::complement(* sequences->getSequence(i));
        vsc->addSequence(*tmp);
        delete tmp;
      }
      delete sequences;
      sequences = vsc;
    }
    // +------------------------+
    // | (Reverse)Transcription |
    // +------------------------+
    else if(actions[i] == "transcript")
    {
      if(sequences->getAlphabet()->getAlphabetType() == AlphabetTools::DNA_ALPHABET.getAlphabetType())
      {
        VectorSequenceContainer * vsc = new VectorSequenceContainer(&AlphabetTools::RNA_ALPHABET);
        for(unsigned int i = 0; i < sequences->getNumberOfSequences(); i++)
        {
          Sequence * tmp = SequenceTools::transcript(* sequences->getSequence(i));
          vsc->addSequence(*tmp);
          delete tmp;
        }
        delete sequences;
        sequences = vsc;
      }
      else if(sequences->getAlphabet()->getAlphabetType() == AlphabetTools::RNA_ALPHABET.getAlphabetType())
      {
        VectorSequenceContainer * vsc = new VectorSequenceContainer(&AlphabetTools::DNA_ALPHABET);
        for(unsigned int i = 0; i < sequences->getNumberOfSequences(); i++)
        {
          Sequence * tmp = SequenceTools::reverseTranscript(* sequences->getSequence(i));
          vsc->addSequence(*tmp);
          delete tmp;
        }
        delete sequences;
        sequences = vsc;
      }
      else throw Exception("Transcription error: input alphabet must be of type 'nucleic'.");
    }
    // +-------------------------------+
    // | Switching nucleotide alphabet |
    // +-------------------------------+
    else if(actions[i] == "switch")
    {
      const Alphabet * alpha = NULL;
      if(sequences->getAlphabet()->getAlphabetType() == AlphabetTools::DNA_ALPHABET.getAlphabetType())
      {
        alpha = &AlphabetTools::RNA_ALPHABET;
      }
      else if(sequences->getAlphabet()->getAlphabetType() == AlphabetTools::RNA_ALPHABET.getAlphabetType())
      {
        alpha = &AlphabetTools::DNA_ALPHABET;
      }
      else throw Exception("Cannot switch alphabet type, alphabet is not of type 'nucleic'.");
      VectorSequenceContainer * vsc = new VectorSequenceContainer(alpha);
      for(unsigned int i = 0; i < sequences->getNumberOfSequences(); i++)
      {
        const Sequence * old = sequences->getSequence(i);
        Sequence * tmp = new Sequence(old->getName(), old->getContent(), old->getComments(), alpha);
        vsc->addSequence(*tmp);
        delete tmp;
      }
      delete sequences;
      sequences = vsc;
    }
    // +-------------+
    // | Translation |
    // +-------------+
    else if(actions[i] == "translate")
    {
      if(!AlphabetTools::isNucleicAlphabet(sequences->getAlphabet()))
        throw Exception("Error in translation: alphabet is not of type 'nucleic'.");
      GeneticCode * gc = NULL;
      string gcstr = ApplicationTools::getStringParameter("sequence.manip_translate.code", params, "Standard");
      if(gcstr == "Standard")
        gc = new StandardGeneticCode(dynamic_cast<const NucleicAlphabet *>(sequences->getAlphabet()));
      else if(gcstr == "VerMito")
        gc = new VertebrateMitochondrialGeneticCode(dynamic_cast<const NucleicAlphabet *>(sequences->getAlphabet()));
      else if(gcstr == "InvMito")
        gc = new InvertebrateMitochondrialGeneticCode(dynamic_cast<const NucleicAlphabet *>(sequences->getAlphabet()));
      else if(gcstr == "EchMito")
        gc = new EchinodermMitochondrialGeneticCode(dynamic_cast<const NucleicAlphabet *>(sequences->getAlphabet()));
      else throw Exception("Unknown genetic code: " + gcstr);

      VectorSequenceContainer * vsc = new VectorSequenceContainer(&AlphabetTools::PROTEIN_ALPHABET);
      for(unsigned int i = 0; i < sequences->getNumberOfSequences(); i++)
      {
        Sequence * tmp = gc->translate(*sequences->getSequence(i));
        vsc->addSequence(*tmp);
        delete tmp;
      }
      delete sequences;
      sequences = vsc;      
    }
    // +-------------+
    // | Remove gaps |
    // +-------------+
    else if(actions[i] == "remove_gaps")
    {
      VectorSequenceContainer * vsc = new VectorSequenceContainer(sequences->getAlphabet());
      for(unsigned int i = 0; i < sequences->getNumberOfSequences(); i++)
      {
        Sequence * tmp = SequenceTools::removeGaps(* sequences->getSequence(i));
        vsc->addSequence(*tmp);
        delete tmp;
      }
      delete sequences;
      sequences = vsc;
    }
    // +---------------------------+
    // | Change gaps to unresolved |
    // +---------------------------+
    else if(actions[i] == "gap2unknown")
    {
      VectorSequenceContainer * vsc = new VectorSequenceContainer(sequences->getAlphabet());
      for(unsigned int i = 0; i < sequences->getNumberOfSequences(); i++)
      {
        Sequence * tmp = new Sequence(* sequences->getSequence(i));
        SymbolListTools::changeGapsToUnknownCharacters(*tmp);
        vsc->addSequence(*tmp);
        delete tmp;
      }
      delete sequences;
      sequences = vsc;
    }
    // +---------------------------+
    // | Change unresolved to gaps |
    // +---------------------------+
    else if(actions[i] == "unknown2gap")
    {
      VectorSequenceContainer * vsc = new VectorSequenceContainer(sequences->getAlphabet());
      for(unsigned int i = 0; i < sequences->getNumberOfSequences(); i++)
      {
        Sequence * tmp = new Sequence(* sequences->getSequence(i));
        SymbolListTools::changeUnresolvedCharactersToGaps(*tmp);
        vsc->addSequence(*tmp);
        delete tmp;
      }
      delete sequences;
      sequences = vsc;
    }
    // +--------------------------+
    // | Resolve dotted alignment |
    // +--------------------------+
    else if(actions[i] == "resolved_dotted")
    {
      SiteContainer * sites = NULL;
      try { sites = dynamic_cast<SiteContainer *>(sequences); }
      catch(exception & e)
      {
        sites = new VectorSiteContainer(*sequences);
        delete sequences;
        sequences = sites;
      }

      const Alphabet * alpha = NULL;
      string alphastr = ApplicationTools::getStringParameter("sequence.manip_resolve_dotted.alphabet", params, "DNA");
      if(alphastr == "DNA") alpha = &AlphabetTools::DNA_ALPHABET;
      else if(alphastr == "RNA") alpha = &AlphabetTools::RNA_ALPHABET;
      else if(alphastr == "Protein") alpha = &AlphabetTools::PROTEIN_ALPHABET;
      else throw Exception("Resolved alphabet must be one of [DNA|RNA|Protein] for solving dotted alignment.");
      OrderedSequenceContainer * tmp = SiteContainerTools::resolveDottedAlignment(*sites, alpha);
      delete sequences;
      sequences = tmp;
    }
    // +---------------------+
    // | Keep complete sites |
    // +---------------------+
    else if(actions[i] == "keep_complete")
    {
      SiteContainer * sites = NULL;
      try
      {
        sites = dynamic_cast<SiteContainer *>(sequences);
        if(!sites) throw exception();
      }
      catch(exception & e)
      {
        sites = new VectorSiteContainer(*sequences);
        delete sequences;
        sequences = sites;
      }
      string maxGapOption = ApplicationTools::getStringParameter("sequence.manip_keep_complete.max_gap_allowed", params, "100%");
      if(maxGapOption[maxGapOption.size()-1] == '%')
      {
        double gapFreq = TextTools::toDouble(maxGapOption.substr(0, maxGapOption.size()-1)) / 100.;
        for(unsigned int i = sites->getNumberOfSites(); i > 0; i--)
        {
          map<int, double> freqs;
          SiteTools::getFrequencies(*sites->getSite(i-1), freqs);
          if(freqs[-1] >= gapFreq) sites->deleteSite(i-1);
        }
      }
      else
      {
        unsigned int gapNum=TextTools::to<unsigned int>(maxGapOption);
        for(unsigned int i = sites->getNumberOfSites(); i > 0; i--)
        {
          map<int, unsigned int> counts;
          SiteTools::getCounts(*sites->getSite(i-1), counts);
          counts[-1]; //Needed in case this entry does not exist in the map. This will set it to 0.
          if(counts[-1] > gapNum) sites->deleteSite(i-1);
        }
        cout << sites->getNumberOfSites() << endl;
        cout << dynamic_cast<SiteContainer *>(sequences)->getNumberOfSites() << endl;
      }
    }
    // +-----------------+
    // | Invert sequence |
    // +-----------------+
    else if(actions[i] == "invert")
    {
      VectorSequenceContainer * vsc = new VectorSequenceContainer(sequences->getAlphabet());
      for(unsigned int i = 0; i < sequences->getNumberOfSequences(); i++)
      {
        const Sequence * old = sequences->getSequence(i);
        Sequence * tmp = SequenceTools::invert(*old);
        vsc->addSequence(*tmp);
        delete tmp;
      }
      delete sequences;
      sequences = vsc;
    }
    else throw Exception("Unknown action: " + actions[i]);
  }
  
  // Write sequences
  sequenceFilePath = ApplicationTools::getAFilePath("output.sequence.file",params, true, false);
  sequenceFormat = ApplicationTools::getStringParameter("output.sequence.format", params, "Fasta");
  ApplicationTools::displayResult("Output sequence format", sequenceFormat);
  OSequence * seqWriter = NULL;
  if(sequenceFormat == "Mase")
  {
    seqWriter = new Mase();
  }
  else if(sequenceFormat == "Phylip")
  {
    bool sequential = true, extended = true;
    if(params.find("output.sequence.format_phylip.order") != params.end())
    {
           if(params["output.sequence.format_phylip.order"] == "sequential" ) sequential = true;
      else if(params["output.sequence.format_phylip.order"] == "interleaved") sequential = false;
      else ApplicationTools::displayWarning("Argument '" +
             params["output.sequence.format_phylip.order"] +
             "' for parameter 'output.sequence.format_phylip.order' is unknown. " +
             "Default used instead: sequential.");
    }
    else ApplicationTools::displayWarning("Argument 'sequence.format_phylip.order' not found. Default used instead: sequential.");
    if(params.find("output.sequence.format_phylip.ext") != params.end())
    {
           if(params["output.sequence.format_phylip.ext"] == "extended") extended = true;
      else if(params["output.sequence.format_phylip.ext"] == "classic" ) extended = false;
      else ApplicationTools::displayWarning("Argument '" +
             params["output.sequence.format_phylip.ext"] +
             "' for parameter 'output.sequence.format_phylip.ext' is unknown. " +
             "Default used instead: extended.");
    }
    else ApplicationTools::displayWarning("Argument 'output.sequence.format_phylip.ext' not found. Default used instead: extended.");
    seqWriter = new Phylip(extended, sequential);
  }
  else if(sequenceFormat == "Fasta") seqWriter = new Fasta();
  else
  {
    ApplicationTools::displayError("Unknown sequence format.");
    exit(-1);
  }
  ApplicationTools::displayResult("Output sequence file", sequenceFilePath);
  seqWriter->write(sequenceFilePath, *sequences, true);
  delete seqWriter;

  delete alphabet;
  delete sequences;

  } catch(exception & e) {
    cout << e.what() << endl;
    exit(-1);
  }

  return (0);
}

