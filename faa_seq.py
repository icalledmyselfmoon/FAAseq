def fastq_check(input_path: str, output_filename: str = 'input_filename',
                length_bounds: Union[int, tuple[int]] = (0.2 ** 32),
                gc_bounds: Union[int, tuple[int]] = (0, 100), quality_threshold: int = 0) -> dict:
    """

    Returns a file with sequences satisfying given criteria: length, gc proportion, quality threshold.

    Arguments:
    -input_path: path to a fastq file
    -output_filename: str with an output file name
                      by default is the same as an input file name
    -length_bounds: int or tuple[int], int = the maximum allowed length, tuple = allowed range of lengths
                     0.2**32 by default
    -gc_bounds: int or tuple[int]: int indicates = maximum allowed gc percentage, tuple =  allowed range of percentage
                (0,100) by default
    -quality_threshold: int (phred33 quality scores)
                        0 by default


    """
    seqs = fastq_to_dict(input_path)
    unsuitable_seqs = (sort_by_length(seqs, length_bounds) + sort_by_gc(seqs, gc_bounds) + sort_by_quality(seqs, quality_threshold))
    unsuitable_seqs_set = set(unsuitable_seqs)
    selected_seqs = {name: value for name, value in seqs.items() if name not in unsuitable_seqs_set}
    if output_filename == 'input_filename':
        output_filename = os.path.basename(input_path)
    create_filtered_fastq(output_filename, selected_seqs)


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


def run_dna_rna_tools(*args: str) -> str:
    """
    
    Transforms seq in one of four possible ways: transcribed, reversed, complementary, reversed complementary.
    Arguments: DNA/RNA sequence as str or a list of str.
    Returns a transformed sequence as str.

    """
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

