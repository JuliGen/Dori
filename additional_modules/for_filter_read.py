def determine_gc_content(read: str) -> float:
    """
    determine the GC content in percent

    :param read: nucleotide sequence (string)
    :return:  gc content in percent (float)
    """

    count_gc = read.count('C') + read.count('G')
    return (count_gc * 100) / len(read)


def determine_lenght(read: str) -> int:
    """
    determine read lenght

    :param read: nucleotide sequence (string)
    :return: read lenght (int)
    """

    return len(read)


def determine_quality(quality: str) -> float:
    """
    determine the average quality of the read (in Phred33)

    :param quality: quality of the read (in ASCII)
    :return: average read quality (float)
    """

    quality_phred_33 = []
    for quality_single_nucl in quality:
        quality_phred_33.append(ord(quality_single_nucl) - 33)
    return sum(quality_phred_33) / len(quality_phred_33)


def create_report_about_filt(seqs: dict[str, tuple[str, str]], filter_seq: dict[str, tuple[str, str]]) -> str:
    """
    create a filtering report (string)

    :param seqs: dictionary with source fastq (dict)
    :param filter_seq: dictionary with filtered fastq (dict)
    :return: filtering report (str)
    """
    inappropriate_seq = seqs.keys() - filter_seq.keys()

    return (f' number of source sequences = {len(seqs)},'
            f'\n number of sequences after filtering = {len(filter_seq)},'
            f'\n number of sequences that did not pass filtering = {len(seqs) - len(filter_seq)}'
            f'\n sequence names that failed validation and were filtered out = {inappropriate_seq}')
