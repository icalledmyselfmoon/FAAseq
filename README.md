# `faa_seq.py`
*Filter, alter and analyse sequences.*

## Overview
faa_seq is a Python package performing operations on nucleic acid and amino acid sequences, as well as fastq data. 
### Capacities:
- For nucleic acid sequences:  reverse, transcribe, get complementary sequence, get reversed complementary sequence.
- For peptide sequences: calculate mass of an amino acid chain, get proportions of amino acids grouped by polarity and charge.  
- Filter fastq reads based on quality scores, GC content and sequence length.

##  Modules  
### run_dna_rna_tools
`run_dna_rna_tools('sequence1','sequence2',..., 'operation')`

The last argument is a string indicating one of four possible operations: 'transcribe', 'complement', 'reverse', 'reverse_complement'. 
  You can pass one or multiple sequences in a form of strings to undergo one of these 4 operations. 

####  Inner functions and syntax 
     
   * seq_transcribe('sequence -> transcribed strand (RNA)
     
   * seq_complement('sequence') -> complementary sequence
     
   * seq_reverse('sequence') -> reversed sequence
     
   * seq_reverse_complement('sequence') -> reversed complementary sequence

#### Features of run_dna_rna_tools:
- takes only  DNA and RNA sequences with canonical nucleotide pairs: 'A-T', 'G-C' for DNA, 'A-U', 'G-C' for RNA.
- preserves the letter case of elements in the string with sequence.

#### Examples for run_dna_rna_tools:

* transcribe DNA into RNA:

      run_dna_rna_tools('ATG', 'transcribe') # 'AUG'

* create a complementary strand(s)

        run_dna_rna_tools('AtG', 'complement') # 'TaC'

* create a reverse complementary strand(s) :

         run_dna_rna_tools('ATg', 'reverse_complement') # 'cAT'

* reverse the sequance(s):

      run_dna_rna_tools('ATG', 'aT', 'reverse') # ['GTA', 'Ta']
  
### fastq_check
  `fastq_check(seqs, length_bounds, gc_bounds, quality_threshold) -> dictionary with filtered sequences`
  
  Arguments:
  - seqs: dictionary with keys -strings (names of sequences): values-tuples of strings (sequence, quality scores)
  
  - length_bounds:
    - int or tuple[int], int =  maximum allowed length, tuple = allowed range of lengths 
    - 0.2**32 by default
  
  - gc_bounds:
       - int or tuple[int]: int = maximum allowed gc percentage, tuple =  allowed range of percentage 
       - (0,100) by default
  
  - quality_threshold:
       - int (phred33 quality scores) 
       - 0 by default

    Output:

     A modified dictionary with keys and values only for sequences satisfying given criteria: length, gc proportion, quality threshold.
### aa_tools
`aa_tools('sequence1','sequence2',..., 'operation')`

Sequences of amino acids should be strings of one letter symbols: 'G', 'A', 'V', 'L', 'I', 'P', 'F', 'M', 'W', 'S', 'T', 'C', 'N', 'Q', 'Y', 'E','D', 'K', 'H', 'R'.
####  Inner functions and syntax 
`calculate_protein_mass(seq)`
This function calculates the mass (Da) of a protein based on its amino acid sequence. As input, it takes a string of amino acids and returns the molecular weight in Da.
Usage example:
```python
aa_tools('MARY', 'calculate_protein_mass') #593 (in Da)
```
`group_amino_acids(seq)`
This function returns a dictionary with proportions of hydrophobic, polar, negatively, and positively charged amino acids in the protein. It takes a string of amino acids as an input.
Usage example:
```python
aa_tools('EEKFG', 'group_amino_acids') #{'hydrophobic': 0.4, 'polar': 0.0, '- charged': 0.4, '+ charged': 0.2}
```

### Installation

To get the tool clone the git repository:

```bash
git clone  https://github.com/icalledmyselfmoon/FAAseq.git && cd FAA_seq
```

  ### Usage

1) To run the script call it from the directory where it's located:

```bash
python faa_seq.py
```
2) Import the module in a code interpreter like Jupyter Notebook, PyCharm, GoogleCollab.

```python
import faa_seq
faa_seq.aa_tools('KAALTR', 'calculate_protein_mass')
```
### Feedback
If you have questions, or found any bug in the program, please write to us at icalledmyselfmoon@gmail.com


     
   
