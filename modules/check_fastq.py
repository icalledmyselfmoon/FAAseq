def sort_by_length(seqs: dict, length_bounds: Union[int, tuple[int]] = (0.2 ** 32)):
    """

    Returns a list of sequence names with length out of the given range().
    Allows to sort out unsuitable sequences based on their length.

    Arguments:
    -dictionary with keys -strings (names of sequences): values-tuples of strings (sequence, quality scores)
    -length_bounds: int or tuple[int], int = the maximum allowed length, tuple = allowed range of lengths
                    0.2**32 by default

    """
    unsuitable_by_length = []
    for name, (sequence, quality) in seqs.items():
        if type(length_bounds) is tuple:
            if length_bounds[0] > len(sequence) > length_bounds[1]:
                unsuitable_by_length.append(name)
        elif type(length_bounds) is int:
            if len(sequence) > length_bounds:
                unsuitable_by_length.append(name)
    return unsuitable_by_length


def sort_by_quality(seq: dict, quality_threshold: int = 0):
    """

    Returns a list of sequence names with the mean of quality scores lower than the threshold value (phred33 score).
    Allows to sort out unsuitable sequences based on the mean of quality scores.

    Arguments:
    -dictionary with keys -strings (names of sequences): values-tuples of strings (sequence, quality scores)
    -quality_threshold: int (phred33 quality scores)
                        0 by default

    """
    unsuitable_by_quality = []
    for name, (sequence, quality) in seq.items():
        q_score = [(ord(sign) - 33) for sign in quality]
        if (sum(q_score) / len(q_score)) < quality_threshold:
            unsuitable_by_quality.append(name)
    return unsuitable_by_quality


def sort_by_gc(seqs: dict, gc_bounds: Union[int, tuple[int]] = (0, 100)):
    """

    Returns a list of sequence names with gc percentage out of the given range (determined by gc_bounds argument).
    Allows to sort out unsuitable sequences based on the percentage of (G + C) in the nucleotide sequence.

    Arguments:
    -seqs: dictionary with keys -strings (names of sequences): values-tuples of strings (sequence, quality scores)
    -qc_bounds: int or tuple[int]: int = the maximum allowed gc percentage, tuple = allowed range of gc percentage
                (0,100) by default

    """
    unsuitable_by_gc = []
    for name, (sequence, quality) in seqs.items():
        gc_score = 100 * ((sequence.count('G') + sequence.count('C')) / len(sequence))
        if type(gc_bounds) is tuple:
            if gc_bounds[0] > gc_score > gc_bounds[1]:
                unsuitable_by_gc.append(name)
        elif type(gc_bounds) is int:
            if gc_score > gc_bounds:
                unsuitable_by_gc.append(name)
    return unsuitable_by_gc


def fastq_to_dict(input_path: str) -> dict:
    """

    Read fastq file and turn it's content into a dictionary.
    Argument: input_path (a path to a fastq file).
    Return type: a dictionary with keys -strings (names of sequences): values-tuples of strings (sequence, quality scores).

    """
    fastq_dict = {}
    fastq_dict = {}
    with open(input_path) as fastq_file:
        for line in fastq_file:
            if line.startswith('@'):
                seq_name = line.strip()
                seq = fastq_file.readline().strip()
                next(fastq_file)
                quality = (fastq_file.readline().strip())
                fastq_dict.update({seq_name: (seq, quality)})
    return fastq_dict


def create_filtered_fastq(output_filename: str, selected_seqs: dict):
    """
    Creates a folder with a fastq file from a dictionary input with filtered reads.
    Arguments:
    -output_filename: str with an output file name
                      by default is the same as an input file name
    -seqs: dictionary with keys -strings (names of sequences): values-tuples of strings (sequence, quality scores)

    """
    output_filename = output_filename
    output_file = os.path.join('fastq_filtrator_results', output_filename)
    input_file = open(input_path, 'r')
    with open(output_file, 'w') as checked_fastq:
        os.makedirs(os.path.dirname(output_file), exist_ok=True)
        lines = input_file.readlines()
        for n, line in enumerate(lines):
            if line.strip() in selected_seqs.keys():
                checked_fastq.write(lines[n])
                checked_fastq.write(lines[n + 1])
                checked_fastq.write(lines[n + 2])
                checked_fastq.write(lines[n + 3])


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
