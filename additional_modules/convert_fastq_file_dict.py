import itertools as it
import os


def convert_fastq_to_dict(input_path: str) -> dict:

    """
    Convert fastq file to dictionary

    :param input_path: path to the fastq file that needs to be filtered
    :return: dictionary with fastq data
    """

    with open(input_path) as fastq_file:
        lines = fastq_file.readlines()

        names = []
        reads = []
        comment = []
        quality = []

        for line in lines[::4]:
            names.append(line.strip())
        for line in lines[1::4]:
            reads.append(line.strip())
        for line in lines[2::4]:
            comment.append(line.strip())
        for line in lines[3::4]:
            quality.append(line.strip())

    seqs = dict(zip(names, it.zip_longest(reads, comment, quality)))
    return seqs


def create_filtr_fastq(output_filename: str, filter_seqs: dict) -> None:

    """
    Convert the filtered dictionary to a fastq file

    :param output_filename: file name with filtered data
    :param filter_seqs: filtered dictionary fastq
    """

    os.makedirs('fastq_filtrator_resuls', exist_ok=True)
    output_dir = 'fastq_filtrator_resuls'

    with open(os.path.join(output_dir, output_filename), mode='w') as file:
        for name, read_com_quality in filter_seqs.items():
            file.write(name + '\n')
            file.write(read_com_quality[0] + '\n')
            file.write(read_com_quality[1] + '\n')
            file.write(read_com_quality[2] + '\n')

    return
