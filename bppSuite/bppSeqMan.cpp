// SPDX-FileCopyrightText: The Bio++ Development Group
//
// SPDX-License-Identifier: CECILL-2.1

// From the STL:
#include <iostream>
#include <fstream>
#include <iomanip>

using namespace std;

// From bpp-core:
#include <Bpp/Version.h>
#include <Bpp/App/BppApplication.h>
#include <Bpp/App/ApplicationTools.h>
#include <Bpp/Text/KeyvalTools.h>

// From bpp-seq:
#include <Bpp/Seq/SiteTools.h>
#include <Bpp/Seq/CodonSiteTools.h>
#include <Bpp/Seq/Alphabet/Alphabet.h>
#include <Bpp/Seq/Alphabet/AlphabetTools.h>
#include <Bpp/Seq/Container/VectorSiteContainer.h>
#include <Bpp/Seq/Container/SiteContainerTools.h>
#include <Bpp/Seq/Container/SequenceContainerTools.h>
#include <Bpp/Seq/App/SequenceApplicationTools.h>
#include <Bpp/Seq/SequenceTools.h>

// From bpp-phyl:
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
  cout << "*           Bio++ Sequence Manipulator, version " << BPP_VERSION << ".           *" << endl;
  cout << "* Author: J. Dutheil                        Last Modif. " << BPP_REL_DATE << " *" << endl;
  cout << "******************************************************************" << endl;
  cout << endl;

  if (args == 1)
  {
    help();
    return 0;
  }

  try
  {
    BppApplication bppseqman(args, argv, "BppSeqMan");
    bppseqman.startTimer();

    // Get alphabet
    shared_ptr<Alphabet> alphabet = SequenceApplicationTools::getAlphabet(bppseqman.getParams(), "", false, true, true);

    shared_ptr<CodonAlphabet> codonAlphabet = dynamic_pointer_cast<CodonAlphabet>(alphabet);

    // Get sequences:
    bool aligned = ApplicationTools::getBooleanParameter("input.alignment", bppseqman.getParams(), false, "", true, 1);
    shared_ptr<SequenceContainerInterface> sequences = 0;

    if (aligned)
    {
      shared_ptr<VectorSiteContainer> allSites = SequenceApplicationTools::getSiteContainer(alphabet, bppseqman.getParams());
      sequences = SequenceApplicationTools::getSitesToAnalyse(*allSites, bppseqman.getParams(), "", true, false);
    }
    else
    {
      sequences = SequenceApplicationTools::getSequenceContainer(alphabet, bppseqman.getParams(), "", true, true);
    }

    ApplicationTools::displayResult("Number of sequences", sequences->getNumberOfSequences());

    // Perform manipulations

    vector<string> actions = ApplicationTools::getVectorParameter<string>("sequence.manip", bppseqman.getParams(), ',', "", "", false, 1);


    for (size_t a = 0; a < actions.size(); a++)
    {
      auto containerWithKeys = dynamic_pointer_cast<VectorSequenceContainer>(sequences);

      string cmdName;
      map<string, string> cmdArgs;
      KeyvalTools::parseProcedure(actions[a], cmdName, cmdArgs);
      ApplicationTools::displayResult("Performing action", cmdName);

      // +-----------------+
      // | Complementation |
      // +-----------------+
      if (cmdName == "Complement")
      {
        shared_ptr<SequenceContainerInterface> sc = 0;
        if (aligned)
          sc = make_shared<VectorSiteContainer>(sequences->getAlphabet());
        else
          sc = make_shared<VectorSequenceContainer>(sequences->getAlphabet());
        for (size_t i = 0; i < sequences->getNumberOfSequences(); i++)
        {
          auto seq = SequenceTools::getComplement(sequences->sequence(i));
          auto name = containerWithKeys ? containerWithKeys->sequenceKey(i) : "seq_" + TextTools::toString(i);
          sc->addSequence(name, seq);
        }
        sequences = sc;
      }
      // +------------------------+
      // | (Reverse)Transcription |
      // +------------------------+
      else if (cmdName == "Transcript")
      {
        if (sequences->getAlphabet()->getAlphabetType() == AlphabetTools::DNA_ALPHABET->getAlphabetType())
        {
          shared_ptr<SequenceContainerInterface> sc = 0;
          if (aligned)
            sc = make_shared<VectorSiteContainer>(AlphabetTools::RNA_ALPHABET);
          else
            sc = make_shared<VectorSequenceContainer>(AlphabetTools::RNA_ALPHABET);
          for (unsigned int i = 0; i < sequences->getNumberOfSequences(); i++)
          {
            auto seq = SequenceTools::transcript(sequences->sequence(i));
            auto name = containerWithKeys ? containerWithKeys->sequenceKey(i) : "seq_" + TextTools::toString(i);
            sc->addSequence(name, seq);
          }
          sequences = sc;
        }
        else if (sequences->getAlphabet()->getAlphabetType() == AlphabetTools::RNA_ALPHABET->getAlphabetType())
        {
          shared_ptr<SequenceContainerInterface> sc = 0;
          if (aligned)
            sc = make_shared<VectorSiteContainer>(AlphabetTools::DNA_ALPHABET);
          else
            sc = make_shared<VectorSequenceContainer>(AlphabetTools::DNA_ALPHABET);
          for (unsigned int i = 0; i < sequences->getNumberOfSequences(); i++)
          {
            auto seq = SequenceTools::reverseTranscript(sequences->sequence(i));
            auto name = containerWithKeys ? containerWithKeys->sequenceKey(i) : "seq_" + TextTools::toString(i);
            sc->addSequence(name, seq);
          }
          sequences = sc;
        }
        else
          throw Exception("Transcription error: input alphabet must be of type 'nucleic'.");
      }
      // +-------------------------------+
      // | Switching nucleotide alphabet |
      // +-------------------------------+
      else if (cmdName == "Switch")
      {
        std::shared_ptr<const Alphabet> alpha = 0;
        if (sequences->getAlphabet()->getAlphabetType() == AlphabetTools::DNA_ALPHABET->getAlphabetType())
        {
          alpha = AlphabetTools::RNA_ALPHABET;
        }
        else if (sequences->getAlphabet()->getAlphabetType() == AlphabetTools::RNA_ALPHABET->getAlphabetType())
        {
          alpha = AlphabetTools::DNA_ALPHABET;
        }
        else
          throw Exception("Cannot switch alphabet type, alphabet is not of type 'nucleic'.");

        shared_ptr<SequenceContainerInterface> sc = 0;
        if (aligned)
          sc = make_shared<VectorSiteContainer>(alpha);
        else
          sc = make_shared<VectorSequenceContainer>(alpha);

        for (size_t i = 0; i < sequences->getNumberOfSequences(); i++)
        {
          auto old = sequences->sequence(i);
          vector<int> content(old.size());
          for (size_t j = 0; j < old.size(); ++j)
          {
            content[j] = old[j];
          }
          auto seq = make_unique<Sequence>(old.getName(), content, old.getComments(), alpha);
          auto name = containerWithKeys ? containerWithKeys->sequenceKey(i) : "seq_" + TextTools::toString(i);
          sc->addSequence(name, seq);
        }
        sequences = sc;
      }
      // +-------------+
      // | Translation |
      // +-------------+
      else if (cmdName == "Translate")
      {
        if (!AlphabetTools::isCodonAlphabet(sequences->getAlphabet().get()))
          throw Exception("Error in translation: alphabet is not of type 'codon'.");

        string codeDesc = ApplicationTools::getStringParameter("genetic_code", bppseqman.getParams(), "Standard", "", true, 1);

        auto gCode = SequenceApplicationTools::getGeneticCode(codonAlphabet->getNucleicAlphabet(), codeDesc);

        ApplicationTools::displayResult("Genetic Code", codeDesc);

        shared_ptr<SequenceContainerInterface> sc = 0;
        if (aligned)
          sc = make_shared<VectorSiteContainer>(AlphabetTools::PROTEIN_ALPHABET);
        else
          sc = make_shared<VectorSequenceContainer>(AlphabetTools::PROTEIN_ALPHABET);

        for (size_t i = 0; i < sequences->getNumberOfSequences(); ++i)
        {
          auto seq = gCode->translate(sequences->sequence(i));
          auto name = containerWithKeys ? containerWithKeys->sequenceKey(i) : "seq_" + TextTools::toString(i);
          sc->addSequence(name, seq);
        }
        sequences = sc;
      }
      // +-------------+
      // | Remove gaps |
      // +-------------+
      else if (cmdName == "RemoveGaps")
      {
        auto sc = make_shared<VectorSequenceContainer>(sequences->getAlphabet());

        for (size_t i = 0; i < sequences->getNumberOfSequences(); i++)
        {
          unique_ptr<Sequence> seq(sequences->sequence(i).clone());
          SequenceTools::removeGaps(*seq);
          auto name = containerWithKeys ? containerWithKeys->sequenceKey(i) : "seq_" + TextTools::toString(i);
          sc->addSequence(name, seq);
        }
        sequences = sc;
        aligned = false;
      }
      // +---------------------------+
      // | Change gaps to unresolved |
      // +---------------------------+
      else if (cmdName == "GapToUnknown")
      {
        shared_ptr<SequenceContainerInterface> sc = 0;
        if (aligned)
          sc = make_shared<VectorSiteContainer>(sequences->getAlphabet());
        else
          sc = make_shared<VectorSequenceContainer>(sequences->getAlphabet());

        for (size_t i = 0; i < sequences->getNumberOfSequences(); i++)
        {
          auto seq = make_unique<Sequence>(sequences->sequence(i));
          SymbolListTools::changeGapsToUnknownCharacters(*seq);
          auto name = containerWithKeys ? containerWithKeys->sequenceKey(i) : "seq_" + TextTools::toString(i);
          sc->addSequence(name, seq);
        }
        sequences = sc;
      }
      // +---------------------------+
      // | Change unresolved to gaps |
      // +---------------------------+
      else if (cmdName == "UnknownToGap")
      {
        shared_ptr<SequenceContainerInterface> sc = 0;
        if (aligned)
          sc = make_shared<VectorSiteContainer>(sequences->getAlphabet());
        else
          sc = make_shared<VectorSequenceContainer>(sequences->getAlphabet());

        for (size_t i = 0; i < sequences->getNumberOfSequences(); i++)
        {
          auto seq = make_unique<Sequence>(sequences->sequence(i));
          SymbolListTools::changeUnresolvedCharactersToGaps(*seq);
          auto name = containerWithKeys ? containerWithKeys->sequenceKey(i) : "seq_" + TextTools::toString(i);
          sc->addSequence(name, seq);
        }
        sequences = sc;
      }

      // +--------------+
      // | Remove stops |
      // +--------------+
      else if (cmdName == "RemoveStops")
      {
        if (!codonAlphabet)
          throw Exception("RemoveStops: requires a codon alphabet.");

        string codeDesc = ApplicationTools::getStringParameter("genetic_code", bppseqman.getParams(), "Standard", "", true, 1);
        ApplicationTools::displayResult("Genetic Code", codeDesc);
        auto gCode = SequenceApplicationTools::getGeneticCode(codonAlphabet->getNucleicAlphabet(), codeDesc);

        auto sites = dynamic_pointer_cast<SiteContainerInterface>(sequences);
        if (!sites)
        {
          auto sc = make_shared<VectorSequenceContainer>(sequences->getAlphabet());
          for (size_t i = 0; i < sequences->getNumberOfSequences(); ++i)
          {
            unique_ptr<Sequence> seq(sequences->sequence(i).clone());
            SequenceTools::removeStops(*seq, *gCode);
            auto name = containerWithKeys ? containerWithKeys->sequenceKey(i) : "seq_" + TextTools::toString(i);
            sc->addSequence(name, seq);
          }
          sequences = sc;
        }
        else
        {
          auto sc = make_unique<VectorSiteContainer>(sequences->getAlphabet());
          for (size_t i = 0; i < sequences->getNumberOfSequences(); ++i)
          {
            unique_ptr<Sequence> seq(sequences->sequence(i).clone());
            SequenceTools::replaceStopsWithGaps(*seq, *gCode);
            auto name = containerWithKeys ? containerWithKeys->sequenceKey(i) : "seq_" + TextTools::toString(i);
            sc->addSequence(name, seq);
          }
          sequences.reset(sc.release());
        }
      }

      // +--------------+
      // | Remove stops |
      // +--------------+
      else if (cmdName == "RemoveColumnsWithStops")
      {
        auto sites = dynamic_pointer_cast<SiteContainerInterface>(sequences);
        if (!sites)
        {
          throw Exception("'RemoveColumnsWithStops' can only be used on alignment. You may consider using the 'CoerceToAlignment' command.");
        }
        if (!codonAlphabet)
          throw Exception("RemoveColumnsWithStops: requires a codon alphabet.");
        string codeDesc = ApplicationTools::getStringParameter("genetic_code", bppseqman.getParams(), "Standard", "", true, 1);

        ApplicationTools::displayResult("Genetic Code", codeDesc);
        auto gCode = SequenceApplicationTools::getGeneticCode(codonAlphabet->getNucleicAlphabet(), codeDesc);

        for (size_t i = sites->getNumberOfSites(); i > 0; i--)
        {
          if (CodonSiteTools::hasStop(sites->site(i - 1), *gCode))
            sites->deleteSite(i - 1);
        }
      }

      // +---------+
      // | Get CDS |
      // +---------+
      else if (cmdName == "GetCDS")
      {
        if (!codonAlphabet)
          throw Exception("GetCDS: requires a codon alphabet.");
        string codeDesc = ApplicationTools::getStringParameter("genetic_code", bppseqman.getParams(), "Standard", "", true, 1);
        ApplicationTools::displayResult("Genetic Code", codeDesc);
        auto gCode = SequenceApplicationTools::getGeneticCode(codonAlphabet->getNucleicAlphabet(), codeDesc);

        shared_ptr<SequenceContainerInterface> sc = 0;
        if (aligned)
          sc = make_shared<VectorSiteContainer>(sequences->getAlphabet());
        else
          sc = make_shared<VectorSequenceContainer>(sequences->getAlphabet());

        for (size_t i = 0; i < sequences->getNumberOfSequences(); ++i)
        {
          auto seq = unique_ptr<Sequence>(sequences->sequence(i).clone());
          size_t len = seq->size();
          SequenceTools::getCDS(*seq, *gCode, false, true, true, false);
          if (aligned)
          {
            for (size_t c = seq->size(); c < len; ++c)
            {
              seq->addElement(seq->getAlphabet()->getGapCharacterCode());
            }
          }
          auto name = containerWithKeys ? containerWithKeys->sequenceKey(i) : "seq_" + TextTools::toString(i);
          sc->addSequence(name, seq);
        }
        sequences = sc;
      }

      // +--------------------------+
      // | Resolve dotted alignment |
      // +--------------------------+
      else if (actions[a] == "CoerceToAlignment")
      {
        shared_ptr<SiteContainerInterface> sites = dynamic_pointer_cast<SiteContainerInterface>(sequences);
        if (!sites)
        {
          sites = make_shared<VectorSiteContainer>(*sequences);
          sequences = sites;
        }
        aligned = true;
      }
      else if (actions[a] == "ResolvedDotted")
      {
        shared_ptr<SiteContainerInterface> sites = dynamic_pointer_cast<SiteContainerInterface>(sequences);
        if (!sites)
        {
          throw Exception("'ResolvedDotted' can only be used on alignment. You may consider using the 'CoerceToAlignment' command.");
        }

        shared_ptr<const Alphabet> alpha = 0;
        string alphastr = ApplicationTools::getStringParameter("alphabet", cmdArgs, "DNA", "", false, 1);
        if (alphastr == "DNA")
          alpha = AlphabetTools::DNA_ALPHABET;
        else if (alphastr == "RNA")
          alpha = AlphabetTools::RNA_ALPHABET;
        else if (alphastr == "Protein")
          alpha = AlphabetTools::PROTEIN_ALPHABET;
        else
          throw Exception("Resolved alphabet must be one of [DNA|RNA|Protein] for solving dotted alignment.");

        shared_ptr<SiteContainerInterface> resolvedCont = SiteContainerTools::resolveDottedAlignment(*sites, alpha);
        sequences = resolvedCont;
      }
      // +---------------------+
      // | Keep complete sites |
      // +---------------------+
      else if (cmdName == "KeepComplete")
      {
        auto sites = dynamic_pointer_cast<SiteContainerInterface>(sequences);
        if (!sites)
        {
          throw Exception("'KeepComplete' can only be used on alignment. You may consider using the 'CoerceToAlignment' command.");
        }

        string maxGapOption = ApplicationTools::getStringParameter("maxGapAllowed", cmdArgs, "100%", "", false, 1);
        if (maxGapOption[maxGapOption.size() - 1] == '%')
        {
          double gapFreq = TextTools::toDouble(maxGapOption.substr(0, maxGapOption.size() - 1)) / 100.;
          for (size_t i = sites->getNumberOfSites(); i > 0; i--)
          {
            map<int, double> freqs;
            SiteTools::getFrequencies(sites->site(i - 1), freqs);
            if (freqs[-1] > gapFreq)
              sites->deleteSite(i - 1);
          }
        }
        else
        {
          size_t gapNum = TextTools::to<size_t>(maxGapOption);
          for (size_t i = sites->getNumberOfSites(); i > 0; i--)
          {
            map<int, size_t> counts;
            SiteTools::getCounts(sites->site(i - 1), counts);
            counts[-1]; // Needed in case this entry does not exist in the map. This will set it to 0.
            if (counts[-1] > gapNum)
              sites->deleteSite(i - 1);
          }
        }
      }
      // +-----------------+
      // | Invert sequence |
      // +-----------------+
      else if (cmdName == "Invert")
      {
        shared_ptr<SequenceContainerInterface> sc = 0;
        if (aligned)
          sc = make_shared<VectorSiteContainer>(sequences->getAlphabet());
        else
          sc = make_shared<VectorSequenceContainer>(sequences->getAlphabet());

        for (size_t i = 0; i < sequences->getNumberOfSequences(); i++)
        {
          auto seq = make_unique<Sequence>(*SequenceTools::getInvert(sequences->sequence(i)).release());
          auto name = containerWithKeys ? containerWithKeys->sequenceKey(i) : "seq_" + TextTools::toString(i);
          sc->addSequence(name, seq);
        }
        sequences = sc;
      }
      // +------------------+
      // | GetCodonPosition |
      // +------------------+
      else if (cmdName == "GetCodonPosition")
      {
        unsigned int pos = ApplicationTools::getParameter<unsigned int>("position", cmdArgs, 3, "", false, 1);

        shared_ptr<SequenceContainerInterface> sc = SequenceContainerTools::getCodonPosition(*sequences, pos - 1);
//      auto sc = shared_ptr<SequenceContainerInterface>(gp);

        if (aligned)
        {
          sequences = make_shared<VectorSiteContainer>(*sc);
        }
        else
        {
          sequences = sc;
        }
      }
      // +-----------------+
      // | FilterFromTree |
      // +-----------------+
      else if (cmdName == "FilterFromTree")
      {
        unique_ptr<Tree> tree(PhylogeneticsApplicationTools::getTree(cmdArgs, ""));
        vector<string> names = tree->getLeavesNames();

        shared_ptr<SequenceContainerInterface> reorderedSequences = 0;

        if (aligned)
        {
          reorderedSequences = make_shared<VectorSiteContainer>(sequences->getAlphabet());
        }
        else
        {
          reorderedSequences = make_shared<VectorSequenceContainer>(sequences->getAlphabet());
        }

        for (size_t i = 0; i < names.size(); ++i)
        {
          auto seq2 = unique_ptr<Sequence>(sequences->sequence(names[i]).clone());
          reorderedSequences->addSequence(names[i], seq2);
        }
        sequences = reorderedSequences;
      }
      // +----------------------+
      // | RemoveEmptySequences |
      // +----------------------+
      else if (cmdName == "RemoveEmptySequences")
      {
        shared_ptr<SequenceContainerInterface> sc = 0;

        if (aligned)
          sc = make_shared<VectorSiteContainer>(sequences->getAlphabet());
        else
          sc = make_shared<VectorSequenceContainer>(sequences->getAlphabet());
        for (size_t i = 0; i < sequences->getNumberOfSequences(); ++i)
        {
          if (SequenceTools::getNumberOfSites(sequences->sequence(i)) != 0)
          {
            auto name = containerWithKeys ? containerWithKeys->sequenceKey(i) : "seq_" + TextTools::toString(i);
            auto seq2 = unique_ptr<Sequence>(sequences->sequence(i).clone());
            sc->addSequence(name, seq2);
          }
        }
        sequences = sc;
      }

      else
        throw Exception("Unknown action: " + cmdName);
    }

    // Write sequences
    ApplicationTools::displayBooleanResult("Final sequences are aligned", aligned);
    if (aligned)
    {
      auto sites = dynamic_pointer_cast<SiteContainerInterface>(sequences);
      SequenceApplicationTools::writeAlignmentFile(*sites, bppseqman.getParams(), "", true, 1);
    }
    else
    {
      SequenceApplicationTools::writeSequenceFile(*sequences, bppseqman.getParams(), "", true, 1);
    }

    bppseqman.done();
  }
  catch (exception& e)
  {
    cout << e.what() << endl;
    return 1;
  }

  return 0;
}
