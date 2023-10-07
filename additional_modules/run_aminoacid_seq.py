AMINO_ACID_WEIGHTS = {
    'A': 89, 'R': 174, 'N': 132, 'D': 133, 'C': 121,
    'E': 147, 'Q': 146, 'G': 75, 'H': 155, 'I': 131,
    'L': 131, 'K': 146, 'M': 149, 'F': 165, 'P': 115,
    'S': 105, 'T': 119, 'W': 204, 'Y': 181, 'V': 117}

MOST_FREQUENT_CODON_FOR_AA_E_COLI = {
    'A': 'GCT', 'R': 'CGT', 'N': 'AAC', 'D': 'GAC', 'C': 'TGC',
    'E': 'GAA', 'Q': 'CAG', 'G': 'GGC', 'H': 'CAC', 'I': 'ATC',
    'L': 'CTG', 'K': 'AAA', 'M': 'ATG', 'F': 'TTC', 'P': 'CCG',
    'S': 'TCT', 'T': 'ACC', 'W': 'TGG', 'Y': 'TAC', 'V': 'GTT',
    'a': 'gct', 'r': 'cgt', 'n': 'aac', 'd': 'gac', 'c': 'tgc',
    'e': 'gaa', 'q': 'cag', 'g': 'ggc', 'h': 'cac', 'i': 'atc',
    'l': 'ctg', 'k': 'aaa', 'm': 'atg', 'f': 'ttc', 'p': 'ccg',
    's': 'tct', 't': 'acc', 'w': 'tgg', 'y': 'tac', 'v': 'gtt'}

DICT_CHARGE_ACID = {
    'negative_charge': ['E', 'D', 'e', 'd'],
    'positive_charge': ['K', 'R', 'H', 'k', 'r', 'h'],
    'neutral_charge': ['V', 'W', 'P', 'w', 'v', 'p', 'i', 'F', 'f', 'm', 'A',
                       'a', 'L', 'M', 'l', 'I', 'S', 's', 'T', 't', 'N', 'n',
                       'Q', 'q', 'C', 'c', 'Y', 'y', 'G', 'g']}

DICT_CLASS_ACID = {
    'hydrophilic': ['t', 'q', 'r', 's', 'y', 'd', 'e', 'g',
                    'c', 'n', 'h', 'k', 'T', 'Q', 'R', 'S',
                    'Y', 'D', 'E', 'G', 'C', 'N', 'H', 'K'],
    'hydrophobic': ['V', 'W', 'P', 'w', 'v', 'p', 'i', 'F',
                    'f', 'm', 'A', 'a', 'L', 'M', 'l', 'I']}

AMINOACID_DICT = {
    'GLY': 'G', 'ALA': 'A', 'VAL': 'V', 'LEU': 'L', 'ILE': 'I', 'MET': 'M',
    'PRO': 'P', 'PHE': 'F', 'TRP': 'W', 'SER': 'S', 'THR': 'T', 'ASN': 'N', 'GLN': 'Q',
    'TYR': 'Y', 'CYS': 'C', 'LYS': 'K', 'ARG': 'R', 'HIS': 'H', 'ASP': 'D', 'GLU': 'E',
    'gly': 'g', 'ala': 'a', 'val': 'v', 'leu': 'l', 'ile': 'i', 'met': 'm',
    'pro': 'p', 'phe': 'f', 'trp': 'w', 'ser': 's', 'thr': 't', 'asn': 'n', 'gln': 'q',
    'tyr': 'y', 'cys': 'c', 'lys': 'k', 'arg': 'r', 'his': 'h', 'asp': 'd', 'glu': 'e',
}

SINGLE_LETTER_ALPHABET = {'G', 'A', 'V', 'L', 'I', 'M', 'P', 'F', 'W', 'S', 'T', 'N', 'Q', 'Y', 'C', 'K', 'R', 'H',
                          'D', 'E',
                          'g', 'a', 'v', 'l', 'i', 'm', 'p', 'f', 'w', 's', 't', 'n', 'q', 'y', 'c', 'k', 'r', 'h',
                          'd', 'e'}
THREE_LETTER_ALPHABET = {'GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PRO', 'PHE', 'TRP', 'SER', 'THR', 'ASN', 'GLN',
                         'TYR', 'CYS', 'LYS', 'ARG', 'HIS', 'ASP', 'GLU',
                         'gly', 'ala', 'val', 'leu', 'ile', 'met', 'pro', 'phe', 'trp', 'ser', 'thr', 'asn', 'gln',
                         'tyr', 'cys', 'lys', 'arg', 'his', 'asp', 'glu',
                         }


def count_molecular_weight(sequence: str) -> int:
    """
    takes an amino acid sequence as input and returns the molecular weight of the protein
    :param sequence: str
    :return: int
    """
    sequence_upper = sequence.upper()
    molecular_weight = sum(AMINO_ACID_WEIGHTS.get(aa, 0) for aa in sequence_upper)
    # Count the molecular weight of protein with using dictionary
    return molecular_weight


def convert_amino_acid_seq_to_dna(sequence: str) -> str:
    """
    takes an amino acid sequence as input and returns the optimal DNA sequence for E.coli
    :param sequence: str
    :return: str
    """

    nucl_str = ''

    for amin_acid in sequence:
        nucl_str += MOST_FREQUENT_CODON_FOR_AA_E_COLI[amin_acid]

    return nucl_str


def determine_charge(amino_seq: str, percent: bool = False) -> dict:
    """
    Takes a string (amino acid sequence),returns the number of positively,
    negatively and neutrally charged amino acids.

    Args:
    - amino_seq - amino acid sequence,
    - percent - optional argument (default False):
    percent = False - output in number of amino acids,
    percent = True - output as a percentage
    """

    charge_amin_acid = []

    for amin_acid in amino_seq:
        for key, values in DICT_CHARGE_ACID.items():
            if amin_acid in values:
                charge_amin_acid.append(key)
    amount_positive = charge_amin_acid.count('positive_charge')
    amount_neutral = charge_amin_acid.count('neutral_charge')
    amount_negative = charge_amin_acid.count('negative_charge')

    if percent:
        result_dict = {"Percentage of positively charged amino acids":
                           (round((amount_positive * 100) / len(amino_seq))),
                       "Percentage of neutrally charged amino acids":
                           (round((amount_neutral * 100) / len(amino_seq))),
                       "Percentage of negatively charged amino acids":
                           (round((amount_negative * 100) / len(amino_seq)))}
    else:
        result_dict = {"Number of positively charged amino acids": amount_positive,
                       "Number of neutrally charged amino acids": amount_neutral,
                       "Number of negatively charged amino acids": amount_negative}

    return result_dict


def determine_polarity(amino_seq: str, percent: bool = False) -> dict:
    """
    Takes a string (amino acid sequence),returns
    a dictionary with the number of hydrophobic and hydrophilic amino acids
    Args:
    - amino_seq - amino acid sequence,
    - percent - optional argument (default False):
    percent = False - output in number of amino acids,
    percent = True - output as a percentage
    """

    class_amin_acid = []

    for amin_acid in amino_seq:
        for key, values in DICT_CLASS_ACID.items():
            if amin_acid in values:
                class_amin_acid.append(key)
    amount_hydrophilic = class_amin_acid.count('hydrophilic')
    amount_hydrophobic = class_amin_acid.count('hydrophobic')

    if percent:
        result_dict = {'Percentage of hydrophilic amino acids':
                           (round((amount_hydrophilic * 100) / len(amino_seq), 2)),
                       'Percentage of hydrophobic amino acids':
                           (round((amount_hydrophobic * 100) / len(amino_seq), 2))}
    else:
        result_dict = {'Number of hydrophilic amino acids': amount_hydrophilic,
                       'Number of hydrophobic amino acids': amount_hydrophobic}

    return result_dict


def translate(sequence: str, record_type: int = 1) -> str:
    """
    Converts one three-letter amino acid sequence notation to one-letter notation
    or one one-letter amino acid sequence notation to three-letter notation

    :param sequence: string with amino acid sequence in three-letter or one-letter notation
    :param record_type: sequence type (one-letter or three-letter). record_type = 3 (three-letter notation),
                  record_type = 1 (one-letter notation, default)
    :return: string with amino acid sequence in one-letter or three-letter notation
    """

    translate_seq = ''

    if record_type == 1:
        for amino_acid in sequence:
            for aa_one_letter, aa_three_letter in AMINOACID_DICT.items():
                if amino_acid == aa_three_letter:
                    translate_seq += aa_one_letter
    else:
        tuple_sequence = [sequence[i:i + 3] for i in range(0, len(sequence), 3)]
        for amino_acid in tuple_sequence:
            translate_seq += AMINOACID_DICT[amino_acid]

    return translate_seq


def count(sequence: str, record_type: int = 1) -> int:
    """
    Counts number of amino acid

    :param sequence: string with amino acid sequence in three-letter or one-letter notation
    :param record_type: sequence type (one-letter or three-letter). record_type = 3 (three-letter notation),
                        record_type = 1 (one-letter notation, default)
    :return: int
    """

    if record_type == 1:
        return len(sequence)
    elif record_type == 3:
        tuple_sequence = [sequence[i:i + 3] for i in range(0, len(sequence), 3)]
        return len(tuple_sequence)


def summary(sequence: str, record_type: int = 1, percent: bool = False) -> str:
    """
     Returns results of all functions
    :param sequence: string with amino acid sequence in three-letter or one-letter notation
    :param record_type: record_type: sequence type (one-letter or three-letter). record_type = 3 (three-letter),
                        record_type = 1 (one-letter notation, default)
    :param percent: optional argument (default False) for 'determine_charge' and 'determine_polarity':
                    percent = False - output in number of amino acids,
                    percent = True - output as a percentage
    :return: dict
    """
    count_seq = count(sequence, record_type)
    translate_seq = translate(sequence, record_type)
    polarity_seq_k = []
    polarity_seq_v = []
    charge_seq_k = []
    charge_seq_v = []

    if record_type == 1:
        for k, v in determine_polarity(sequence, percent).items():
            polarity_seq_k.append(k)
            polarity_seq_v.append(v)
    else:
        for k, v in determine_polarity(translate(sequence, record_type), percent).items():
            polarity_seq_k.append(k)
            polarity_seq_v.append(v)

    if record_type == 1:
        for k, v in determine_charge(sequence, percent).items():
            charge_seq_k.append(k)
            charge_seq_v.append(v)
    else:
        for k, v in determine_charge(translate(sequence, record_type), percent).items():
            charge_seq_k.append(k)
            charge_seq_v.append(v)

    if record_type == 1:
        molecular_weight_seq = count_molecular_weight(sequence)
    else:
        molecular_weight_seq = count_molecular_weight(translate(sequence, record_type))

    if record_type == 1:
        convert_seq_in_dna = convert_amino_acid_seq_to_dna(sequence)
    else:
        convert_seq_in_dna = convert_amino_acid_seq_to_dna(translate(sequence, record_type))

    return (f'Number of amino acids: {count_seq} \n'
            f'Translated amino acid sequence into another format: {translate_seq} \n'
            f'{polarity_seq_k[0]}: {polarity_seq_v[0]} \n'
            f'{polarity_seq_k[1]}: {polarity_seq_v[1]}\n'
            f'{charge_seq_k[0]}:{charge_seq_v[0]} \n'
            f'{charge_seq_k[1]}:{charge_seq_v[1]} \n'
            f'{charge_seq_k[2]}:{charge_seq_v[2]} \n'
            f'Molecular weight: {molecular_weight_seq} \n'
            f'Converted to DNA sequence: {convert_seq_in_dna}')
