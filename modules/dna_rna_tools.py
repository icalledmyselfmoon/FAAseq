COMP_BASES_DNA: dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't',
                        'g': 'c', 'c': 'g', 't': 'a'}
COMP_BASES_RNA: dict = {'A': 'U', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 'u',
                        'g': 'c', 'c': 'g', 'u': 'a'}


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
