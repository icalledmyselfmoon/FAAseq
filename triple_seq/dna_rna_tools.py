COMP_BASES_DNA: dict = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 't',
                        'g': 'c', 'c': 'g', 't': 'a'}
COMP_BASES_RNA: dict = {'A': 'U', 'C': 'G', 'G': 'C', 'T': 'A', 'a': 'u',
                        'g': 'c', 'c': 'g', 'u': 'a'}


def seq_transcribe(seq):
    rna_seq = seq.replace('t', 'u').replace('T', 'U')
    return rna_seq


def seq_reverse(seq):
    return seq[::-1]


def seq_complement(seq):
    if 'T' in set(seq.upper()):
        comp_seq = [COMP_BASES_DNA[nuc] for nuc in seq_list]
    else:
        comp_seq = [COMP_BASES_RNA[nuc] for nuc in seq_list]
    return ''.join(comp_seq)


def seq_reverse_complement(seq):
    return seq_complement(seq[::-1])


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
