






AA_MASSES: dict = {'A': 89, 'R': 174, 'N': 132, 'D': 133, 'C': 121, 'Q': 146, 'E': 147, 'Z': 147,
                   'G': 75, 'H': 155, 'I': 131, 'L': 131, 'K': 146, 'M': 149, 'F': 165, 'P': 115, 'S': 105,
                   'T': 119, 'W': 204, 'Y': 181, 'V': 117}
AA_GROUPS: dict = {'hydrophobic': ['G', 'A', 'V', 'L', 'I', 'P', 'F', 'M', 'W'],
                   'polar': ['S', 'T', 'C', 'N', 'Q', 'Y'],
                   '- charged': ['E', 'D'], '+ charged': ['K', 'H', 'R']}


def seq_transcribe(seq: str) -> str:
    """

    Returns a transcribed strand (RNA).
    Arguments: DNA sequence as str or a list of str. Input of RNA sequence is returned unchanged.

    """
    rna_seq = seq.replace('t', 'u').replace('T', 'U')
    return rna_seq


def seq_reverse(seq: str) -> str:
    """

    Returns a reversed sequence.
    Arguments: DNA/RNA sequence as str or a list of str

    """
    return seq[::-1]


def seq_complement(seq: str) -> str:
    """

    Returns a complementary sequence.
    Arguments: DNA/RNA sequence as str or a list of str

    """
    if 'T' in set(seq.upper()):
        comp_seq = [COMP_BASES_DNA[nuc] for nuc in seq_list]
    else:
        comp_seq = [COMP_BASES_RNA[nuc] for nuc in seq_list]
    return ''.join(comp_seq)


def seq_reverse_complement(seq: str) -> str:
    """

    Returns a complementary sequence.
    Arguments: DNA/RNA sequence as str or a list of str

    """
    return seq_complement(seq[::-1])


def run_dna_rna_tools(*args: str) -> str:
    """

    Transforms seq in one of four possible ways: transcribed, reversed, complementary, reversed complementary.
    Arguments: DNA/RNA sequence as str or a list of str.
    Returns a sequence as str.

    """
    *seqs, operation = args
    results = []
    for seq in seqs:
        if set(seq.upper()) <= {'A', 'T', 'C', 'G'} or set(seq.upper()) <= {'A', 'U', 'C', 'G'}:
            operations = {'transcribe': seq_transcribe, 'reverse': seq_reverse, 'complement': seq_complement,
                          'reverse_complement': seq_reverse_complement}
            new_seq = operations[operation](seq)
            results.append(new_seq)
        else:
            pass
        if len(results) != 1:
            return results
        else:
            return results[0]
def calculate_protein_mass(seq: str) -> float:
    """

    Calculates the mass (Da) of a protein based on its amino acids sequence.
    Takes a string of amino acids, returns the molecular weight in Da.
    Amino acids in the string should be indicated as one-letter symbols.

    """
    aa_seq = seq.upper()
    mass = 0
    for amino_acid in aa_seq:
        mass += AA_MASSES[amino_acid]
    return mass


def group_amino_acids(seq: str) -> dict:
    """

     Returns a dictionary with proportions of hydrophobic, polar, negatively and positively charged amino acids in the protein.
     Takes a string of amino acids, returns a dictionary.
     Amino acids in the string should be indicated as one-letter symbols.

    """
    aa_seq = seq.upper()
    profile = {'hydrophobic': 0.0, 'polar': 0.0, '- charged': 0.0, '+ charged': 0.0}

    for amino_acid in aa_seq:
        for group_name, group_list in aa_biochemistry.items():
            if amino_acid in group_list:
                profile[group_name] += 1
    for group, count in profile.items():
        profile[group] = round((count / len(seq)), 2)
    return profile


def aa_tools(*args):
    """

    Main function for amino acid sequences processing.
    Parameters: *args (str) - amino acid sequences and operation.
    Returns: -None if non-amino acid chars found
             -Dict with proportions for group_amino_acids operation
             -Int with protein mass in Da for calculate_protein_mass operation

    """
    *seqs, operation = args

    results = []
    for seq in seqs:
        if set(seq.upper()) <= {'G', 'A', 'V', 'L', 'I', 'P', 'F', 'M', 'W', 'S', 'T', 'C', 'N', 'Q', 'Y', 'E',
                                'D', 'K', 'H', 'R'}:
            operations = {'calculate_protein_mass': calculate_protein_mass, 'group_amino_acids': group_amino_acids}
            results = operations[operation](seq)
            return results
        else:
            return None
