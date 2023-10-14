import os


def convert_multiline_fasta_to_oneline(input_fasta, output_fasta='input_fasta'):
    input_seq = open(input_fasta, 'r')
    if output_fasta == 'input_fasta':
        output_fasta = os.path.basename(input_fasta)
    else:
        output_fasta = output_fasta + '.fasta'
    with open(output_fasta, 'w') as output_seq:
        for line in input_seq:
            output_seq.write(input_seq.readline().strip())
