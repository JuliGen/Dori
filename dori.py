from typing import Union
import additional_modules.for_filter_read as ffr
import additional_modules.dna_rna_tools as drt
import additional_modules.run_aminoacid_seq as amino


def run_dna_rna_tools(*args: str, action: str) -> Union[str, list]:
    """
    Performs the following list of operations with DNA or RNA:
    -transcribe - transcribes DNA sequences into RNA,
                  if you give RNA as input, the sequence will not be transcribed.
    -reverse - creating a reverse sequence
    -complement - creating a complementary sequence
    -reverse_complement - creating a reverse sequence and then creating its complementary sequence

    The data is first automatically checked to determine whether it belongs to DNA or RNA.
    If any of the sequences are not DNA or RNA, the data will be partially processed.
    In this case, the program output will indicate which of the sequences were not processed.

    :param args: strings with DNA or RNA sequences
    :param action: string with the name of the operation that needs to be performed (the name is stated above)
    :return: string or list with processed sequences
    """

    if action not in ["transcribe", "reverse", "complement", "reverse_complement"]:
        raise ValueError('Please provide the action to be used or check the spelling of the action')

    seqs = []
    wrong_seqs = {}
    res = []

    for seq in args:
        if (drt.is_dna_or_rna(seq) == 'DNA') or (drt.is_dna_or_rna(seq) == 'RNA'):
            seqs.append(seq)
        else:
            wrong_seqs[args.index(seq) + 1] = seq

    for seq in seqs:
        if action == "transcribe":
            res.append(drt.transcribe(seq))
        if action == "reverse":
            res.append(drt.reverse(seq))
        if action == "complement":
            res.append(drt.complement(seq))
        if action == "reverse_complement":
            res.append(drt.reverse_complement(seq))

    if not res:
        raise ValueError('Provide DNA or RNA sequences')

    if len(res) == 1:
        res = res[0]

    if wrong_seqs != {}:
        print(f'Sequences that are not DNA or RNA (number in the input list: seq): {wrong_seqs}')

    print("Result:")
    return res


def filter_read(seqs:  dict[str, tuple[str, str]],
                gc_bounds: Union[tuple[Union[float, int]], Union[float, int]] = (0, 100),
                length_bounds: Union[tuple[int], int] = (0, 2 ** 32),
                quality_threshold: int = 0,
                report: bool = False) -> Union[dict, dict and str]:
    """
    Filter the sequence dictionary (fastq) using the following parameters:
    - GC content
    - sequence length
    - sequence quality

    :param seqs: sequences dictionary (dict)
    :param gc_bounds: gc content range to filtering (in percent);
           sequences that are not included in the range are discarded;
           if you pass one number, this number will be considered the upper limit;
           default = (0, 100) - all reads will be saved, regardless of gc composition
    :param length_bounds: sequence length range for filtering;
           sequences that are not included in the range are discarded;
           if you pass one number, this number will be considered the upper limit;
           default = (0, 2 ** 32)
    :param quality_threshold: average sequence quality threshold for filtering (Phred33 scale);
           sequences with quality below the threshold are discarded
    :param report: report = True (show filtering report), report = True (default)
    :return: A filtered dictionary containing sequences (fastq)
             that have passed all tests (according to specified parameters)
    """

    filter_seqs = {}

    if type(gc_bounds) != tuple:
        gc_bounds = (0, gc_bounds)
    if type(length_bounds) != tuple:
        length_bounds = (0, length_bounds)

    for name, read_quality in seqs.items():
        read = read_quality[0]
        quality = read_quality[1]
        if gc_bounds[0] <= ffr.determine_gc_content(read) <= gc_bounds[1]:
            if length_bounds[0] <= ffr.determine_lenght(read) <= length_bounds[1]:
                if ffr.determine_quality(quality) > quality_threshold:
                    filter_seqs[name] = read_quality

    if report:
        print(ffr.create_report_about_filt(seqs, filter_seqs))

    return filter_seqs


def run_aminoacid_seq(sequence: str, function: str = 'summary', record_type: int = 1,
                      percent: bool = False) -> Union[str, dict, int]:
    """
    Performs the following list of operations:
    count - counts number of amino acid
    translate - converts one record type to another
    determine_charge - counts number or percent of amino acid with different charges
    determine_polarity - counts number or percent of amino acid with different polarity
    convert_amino_acid_seq_to_dna - takes an amino acid sequence as input and returns the optimal DNA sequence for E.coli
    count_molecular_weight - takes an amino acid sequence as input and returns the molecular weight of the protein
    summary - returns results of all functions (default)

    Arguments:
    - sequence:str - sequence for function
    - function:str - name of the function you need to perform. You can use: 'count', 'translate', 'summary'(default)
    - record_type:int - record type of your sequence. You can use 1(default) for single letter type or 3 for three letter type
    - percent: bool - for determine_charge and determine_polarity shows result in percent (True) or in number (False). Default - False

    Return:
    - count - int
    - translate - str
    - determine_charge - dict
    - determine_polarity - dict
    - convert_amino_acid_seq_to_dna - str
    - count_molecular_weight - int
    - summary - str

    """

    if record_type == 1:
        if not set(sequence).issubset(amino.SINGLE_LETTER_ALPHABET):
            raise ValueError('Please provide amino acid sequence'
                             ' or check that the entered parameter record_type is correct')
    elif record_type == 3:
        tuple_sequence = [sequence[i:i + 3] for i in range(0, len(sequence), 3)]
        if not set(tuple_sequence).issubset(amino.THREE_LETTER_ALPHABET):
            raise ValueError('Please provide amino acid sequence '
                             ' or check that the entered parameter record_type is correct')

    if function == 'count':
        return amino.count(sequence, record_type)
    elif function == 'translate':
        return amino.translate(sequence, record_type)
    elif function == 'summary':
        return amino.summary(sequence, record_type, percent)
    elif function == 'determine_polarity':
        if record_type == 3:
            return amino.determine_polarity(amino.translate(sequence, record_type), percent)
        else:
            return amino.determine_polarity(sequence, percent)
    elif function == 'determine_charge':
        if record_type == 3:
            return amino.determine_charge(amino.translate(sequence, record_type), percent)
        else:
            return amino.determine_charge(sequence, percent)
    elif function == 'count_molecular_weight':
        if record_type == 3:
            return amino.count_molecular_weight(amino.translate(sequence, record_type))
        else:
            return amino.count_molecular_weight(sequence)
    elif function == 'convert_amino_acid_seq_to_dna':
        if record_type == 3:
            return amino.convert_amino_acid_seq_to_dna(amino.translate(sequence, record_type))
        else:
            return amino.convert_amino_acid_seq_to_dna(sequence)
