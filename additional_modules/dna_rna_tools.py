DNA_ALPHABET = set('ATCGatcg')
RNA_ALPHABET = set('AUCGaucg')
COMPLEMENT_NUCL_DICT_DNA = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}
COMPLEMENT_NUCL_DICT_RNA = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G', 'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}


def is_dna_or_rna(seq: str) -> str:
    """
    Determine whether a sequence (string) belongs to DNA or RNA

    :param seq: string with sequence
    :return: string
    """

    if (set(seq).issubset(DNA_ALPHABET)) or (set(seq).issubset(RNA_ALPHABET)):
        if ('u' in seq) or ('U' in seq):
            result = 'RNA'
        else:
            result = 'DNA'
    else:
        result = f'{seq} - is not a DNA or RNA sequence'
    return result


def transcribe(seq: str) -> str:
    """
    Transcribes DNA sequences into RNA,
    if you give RNA as input, the sequence will not be transcribed.

    :param seq: string with sequence
    :return: string
    """

    if is_dna_or_rna(seq) == 'DNA':
        trans_seq = seq.replace('T', 'U').replace('t', 'u')
    else:
        raise ValueError("Provide DNA, RNA sequences cannot be transcribed")
    return trans_seq


def reverse(seq: str) -> str:
    """
    Create a reverse sequence

    :param seq: string with sequence
    :return: string with reverse sequence
    """

    reverse_seq = seq[::-1]
    return reverse_seq


def complement_dna(seq: str) -> str:
    """
    Create a complementary DNA sequence

    :param seq: string with DNA sequence
    :return: string with complementary sequence
    """

    complement_seq = ''
    for nucl in seq:
        complement_seq += COMPLEMENT_NUCL_DICT_DNA[nucl]
    return complement_seq


def complement_rna(seq: str) -> str:
    """
    Create a complementary RNA sequence

    :param seq: string with RNA sequence
    :return: string with complementary sequence
    """
    complement_seq = ''
    for nucl in seq:
        complement_seq += COMPLEMENT_NUCL_DICT_RNA[nucl]
    return complement_seq


def complement(seq: str) -> str:
    """
    Create a complementary sequence

    :param seq: string with sequence DNA or RNA
    :return: string with complementary sequence
    """

    if is_dna_or_rna(seq) == 'DNA':
        res = complement_dna(seq)
    if is_dna_or_rna(seq) == 'RNA':
        res = complement_rna(seq)

    return res


def reverse_complement(seq: str) -> str:
    """
    Create a reverse sequence and then creating its complementary sequence

    :param seq: string with sequence DNA or RNA
    :return: string with complementary sequence for reverse sequence
    """

    reverse_seq = seq[::-1]
    return complement(reverse_seq)
