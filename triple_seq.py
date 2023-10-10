def fastq_check(seqs: dict, length_bounds: Union[int, tuple[int]] = (0.2 ** 32),
                gc_bounds: Union[int, tuple[int]] = (0, 100), quality_threshold: int = 0) -> dict:
    """

    Returns a dictionary with sequences satisfying given criteria: length, gc proportion, quality threshold.

    Arguments:
    -dictionary with keys -strings (names of sequences): values-tuples of strings (sequence, quality scores)
    -length_bounds: int or tuple[int], int = the maximum allowed length, tuple = allowed range of lengths
                     0.2**32 by default
    -gc_bounds: int or tuple[int]: int indicates = maximum allowed gc percentage, tuple =  allowed range of percentage
                (0,100) by default
    -quality_threshold: int (phred33 quality scores)
                        0 by default


    """
    unsuitable_seqs = (sort_by_length(seqs) + sort_by_gc(seqs) + sort_by_quality(seqs))
    unsuitable_seqs_set = set(unsuitable_seqs)
    selected_seqs = {name: value for name, value in seqs.items() if name not in unsuitable_seqs_set}
    return selected_seqs


def aa_tools(*args):
    """

    Main function for amino acid sequences processing.
    Parameters: *args (str) - amino acid sequences and operation.
    Returns: -None if non-amino acid chars found
             -Dict with proportions for group_amino_acids operation
             -Int with protein mass in Da for calculate_protein_mass operation

    """
    *seqs, operation = args

    answer = []
    for seq in seqs:
        if set(seq.upper()) <= {'G', 'A', 'V', 'L', 'I', 'P', 'F', 'M', 'W', 'S', 'T', 'C', 'N', 'Q', 'Y', 'E',
                                'D', 'K', 'H', 'R'}:
            operations = {'calculate_protein_mass': calculate_protein_mass, 'group_amino_acids': group_amino_acids}
            answer = operations[operation](seq)
            return answer
        else:
            return None


def run_dna_rna_tools(*args):
    *seqs, operation = args
    answer = []
    for seq in seqs:
        if set(seq.upper()) <= {'A', 'T', 'C', 'G'} or set(seq.upper()) <= {'A', 'U', 'C', 'G'}:
            operations = {'transcribe': seq_transcribe, 'reverse': seq_reverse, 'complement': seq_complement,
                          'reverse_complement': seq_reverse_complement}
            new_seq = operations[operation](seq)
            answer.append(new_seq)
        else:
            pass
        if len(answer) != 1:
            return answer
        else:
            return answer[0]
