import os
from typing import List
import string


def convert_multiline_fasta_to_oneline(input_fasta, output_fasta='input_fasta'):
    """
    Create a fasta file with sequences placed in one line from a multiline fasta file.
    Input: fasta file with sequences placed on separate lines
    Output: fasta file with aligned sequences.

    """
    input_seq = open(input_fasta, 'r')
    if output_fasta == 'input_fasta':
        output_fasta = os.path.basename(input_fasta) + '_1line'
    else:
        output_fasta = output_fasta + '.fasta'
    with open(output_fasta, 'w') as output_seq:
        for line in input_seq:
            output_seq.write(input_seq.readline().strip())


def remove_punct(input_string: str) -> str:
    """

    Returns a string without punctuation signs

    """
    translator = str.maketrans("", "", string.punctuation)
    result = input_string.translate(translator)
    return result


def select_genes_from_gbk_to_fasta(input_gbk: str, genes: List[str], n_before: int = 1, n_after: int = 1,
                                   output_fasta='input_fasta'):
    """
    Create a file with sequences for genes of interest
    Include sequences flanking the gene of interest - parse number of genes before and after CDS of interest

    Arguments:
    input_gbk: gbk file
    genes: list of string with gene names
    n_before: number of genes before CDS of interest, 1 by default
    n_after: number of genes before CDS of interest, 1 by default

    Output: fasta file with selected sequences. The output file name is the same as of an input by default


    """
    if output_fasta == 'input_fasta':
        output_fasta = os.path.basename(input_gbk).replace('gbk', 'fasta')
    else:
        output_fasta = output_fasta + '.fasta'
    with open('example.gbk', 'r') as gbk_file:
        gene_lines = []
        seq_lines = []
        counter = 0
        for line in gbk_file:
            if '/gene' in line:
                counter += 1
                gene_lines.append((counter, remove_punct((line.strip()).replace("/gene", ""))))
            elif "/translation=" in line:
                seq_lines.append((counter, remove_punct((line.strip()).replace("/translation=", ""))))
            elif remove_punct(line.strip()).isalpha() or remove_punct(line.strip()).isalpha():
                seq_lines.append((counter, remove_punct(line.strip())))
        genes_seqs_list = []
        for index1, (n1, gene) in enumerate(gene_lines):
            genes_seqs_list.append([gene, ''])
            for (n2, seq) in seq_lines:
                if n1 == n2:
                    genes_seqs_list[index1][1] += seq
        interest_gene_numbers = [n for n, [g, s] in enumerate(genes_seqs_list) if g in genes]
        selected_gene_indexes = list()
        for n in interest_gene_numbers:
            count1 = 0
            while (n - n_before) >= 0 and count1 < n_before:
                if n_before > 0:
                    n -= 1
                    count1 += 1
                selected_gene_indexes.append(n)
        for n in interest_gene_numbers:
            count2 = 0
            while (n + n_after) <= len(genes_seqs_list) and count2 < n_after:
                if n_after > 0:
                    n += 1
                    count2 += 1
                selected_gene_indexes.append(n)
        selected_gene_numbers = selected_gene_indexes + interest_gene_numbers
        selected_gene_numbers = sorted(list(set(selected_gene_numbers)))
        selected_seqs = [s for i, [g, s] in enumerate(genes_seqs_list) if i in selected_gene_numbers]
        with open(output_fasta, 'w') as fasta_file:
            for s in selected_seqs:
                fasta_file.write(s + '\n')
