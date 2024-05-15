import datetime
import os
import sys
import time
from abc import ABC, abstractmethod
from dataclasses import dataclass
from io import StringIO
from typing import Union

import numpy as np
import requests
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqUtils import gc_fraction
from bs4 import BeautifulSoup
from dotenv import load_dotenv

# Color settings for text output
BLUE = "\033[36m"
BOLD = "\033[1m"
GREEN = "\033[32m"
START_SETTINGS = "\033[0m"


def filter_read(input_path: str, output_filename: str = None,
                gc_bounds: Union[tuple[Union[float, int]], Union[float, int]] = (0, 100),
                length_bounds: Union[tuple[int], int] = (0, 2 ** 32),
                quality_threshold: int = 0,
                report: bool = False):
    """
    Filter the fastq file using the following parameters:
    - GC content
    - sequence length
    - sequence quality

    :param input_path: path to the fastq file that needs to be filtered
    :param output_filename: file name with filtered data
    :param gc_bounds: gc_bounds: gc content range to filtering (in percent);
           sequences that are not included in the range are discarded;
           if you pass one number, this number will be considered the upper limit;
           default = (0, 100) - all reads will be saved, regardless of gc composition
    :param length_bounds:sequence length range for filtering;
           sequences that are not included in the range are discarded;
           if you pass one number, this number will be considered the upper limit;
           default = (0, 2 ** 32)
    :param quality_threshold:average sequence quality threshold for filtering (Phred33 scale);
           sequences with quality below the threshold are discarded
    :param report:report = True (show filtering report), report = True (default)
    :return:A filtered fastq file containing sequences
             that have passed all tests (according to specified parameters)
    """

    if not isinstance(gc_bounds, tuple):
        gc_bounds = (0, gc_bounds)
    if not isinstance(length_bounds, tuple):
        length_bounds = (0, length_bounds)
    if output_filename is None:
        output_filename = os.path.splitext(os.path.basename(input_path))[0]
    output_filename = f'{output_filename}.fastq'

    os.makedirs('fastq_filtrator_resuls', exist_ok=True)
    output_dir = 'fastq_filtrator_resuls'

    count_filtered_record = 0
    count_record = 0

    with open(os.path.join(output_dir, output_filename), mode='w') as out_file:
        for record in SeqIO.parse("reads.3.fastq", "fastq"):
            count_record += 1
            if gc_bounds[0] <= gc_fraction(record.seq) * 100 <= gc_bounds[1]:
                if length_bounds[0] <= len(record.seq) <= length_bounds[1]:
                    if np.mean(record.letter_annotations['phred_quality']) > quality_threshold:
                        count_filtered_record += 1
                        SeqIO.write(record, out_file, "fastq")

    if report:
        print(f' Number of source sequences = {count_record},'
              f'\n Number of sequences after filtering = {count_filtered_record},')

    return f'{output_filename} file created!'


class BiologicalSequence(ABC, str):

    def __init__(self, seq: str):
        super().__init__()
        self.seq = seq

    @abstractmethod
    def __len__(self):
        """Return the length of the sequence"""
        return super().__len__()

    @abstractmethod
    def __getitem__(self, key: Union[int, slice]) -> object:
        """Return a subsequence in the form of a single letter or a string of several letters"""

        if isinstance(key, int):
            if 0 >= key or key > len(self.seq):
                raise IndexError(f"Index {key} is out of range of the sequence")
            return super().__getitem__(key - 1)
        elif isinstance(key, slice):
            if key.start <= 0 or key.stop > len(self.seq):
                raise IndexError(f"Index {key} is out of range of the sequence")
            new_slice = slice(key.start - 1, key.stop, key.step)
            return super().__getitem__(new_slice)

    @abstractmethod
    def _check_alphabet(self, seq: str):
        """Checks the sequence to match the alphabet"""
        pass

    def __str__(self):
        return f'{self.__class__.__name__}("{self.seq}")'


class NucleicAcidSequence(BiologicalSequence):
    alphabet = set('ATCGUatcgu')
    compliment_dict = None

    def __init__(self, seq: str):
        if self._check_alphabet(seq):
            super().__init__(seq)
        else:
            raise ValueError(f"The sequence does not match the alphabet for the {self.__class__.__name__}")

    def __getitem__(self, key: Union[int, slice]) -> object:
        """Returns a subsequence in the form of a single letter or a string of several letters"""
        return super().__getitem__(key)

    def __len__(self) -> int:
        """Return the length of the sequence"""
        return super().__len__()

    def _check_alphabet(self, seq: str) -> bool:
        """Checks the sequence to match the alphabet"""
        return set(seq).issubset(type(self).alphabet)

    def complement(self):
        """Returns  a complementary sequence"""
        if type(self) is NucleicAcidSequence:
            raise NotImplementedError("The compliment method cannot be called for an object of the NucleicAcidSequence")
        complement_seq = ''
        for nucl in self.seq:
            complement_seq += type(self).compliment_dict[nucl]
        return self.__class__(complement_seq)

    def gc_content(self) -> float:
        """Returns the gc content in percent"""
        nucl_list = list(self.seq.upper())
        gc_count = nucl_list.count("G") + nucl_list.count("C")
        return (gc_count / len(nucl_list)) * 100

    def reverse(self):
        """Create a reverse sequence"""
        return self.__class__(self.seq[::-1])


class DNASequence(NucleicAcidSequence):
    alphabet = set('ATCGatcg')
    compliment_dict = {'A': 'T', 'T': 'A', 'G': 'C', 'C': 'G', 'a': 't', 't': 'a', 'g': 'c', 'c': 'g'}

    def transcribe(self):
        """Returns the transcribed DNA"""
        return RNASequence(self.seq.replace('T', 'U').replace('t', 'u'))


class RNASequence(NucleicAcidSequence):
    alphabet = set('ACUGacug')
    compliment_dict = {'A': 'U', 'U': 'A', 'G': 'C', 'C': 'G', 'a': 'u', 'u': 'a', 'g': 'c', 'c': 'g'}


class AminoAcidSequence(BiologicalSequence):
    one_letter_alphabet = set('GAVLIMPFWSTNQYCKRHDE')
    three_letter_alphabet = {'GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PRO', 'PHE', 'TRP', 'SER', 'THR', 'ASN', 'GLN',
                             'TYR', 'CYS', 'LYS', 'ARG', 'HIS', 'ASP', 'GLU'
                             }
    amino_acid_weights = {
        'A': 89, 'R': 174, 'N': 132, 'D': 133, 'C': 121,
        'E': 147, 'Q': 146, 'G': 75, 'H': 155, 'I': 131,
        'L': 131, 'K': 146, 'M': 149, 'F': 165, 'P': 115,
        'S': 105, 'T': 119, 'W': 204, 'Y': 181, 'V': 117,
        'GLY': 75, 'ALA': 89, 'VAL': 117, 'LEU': 131, 'ILE': 131, 'MET': 149,
        'PRO': 115, 'PHE': 165, 'TRP': 204, 'SER': 105, 'THR': 119, 'ASN': 132, 'GLN': 146,
        'TYR': 181, 'CYS': 121, 'LYS': 146, 'ARG': 174, 'HIS': 155, 'ASP': 133, 'GLU': 147,
    }

    def __new__(cls, seq, record_type='one_letter_sequence'):
        obj = super().__new__(cls, seq)
        return obj

    def __init__(self, seq: str, record_type: str = "one_letter_sequence"):
        """
        Create a AminoAcidSequence object
        :param seq: Amino acid sequence
        :param record_type:str - record type of your sequence.
        You can use "one_letter_sequence"(default) for single letter type
        or "three_letter_sequence" for three letter type
         """
        if self._check_alphabet(seq, record_type):
            super().__init__(seq)
            self.record_type = record_type
        else:
            raise ValueError(f"The sequence does not match the alphabet for the {self.__class__.__name__}")

    def _check_alphabet(self, seq: str, record_type: str) -> bool:
        """Checks the sequence to match the alphabet"""
        if record_type == "one_letter_sequence":
            return set(seq.upper()).issubset(type(self).one_letter_alphabet)
        elif record_type == "three_letter_sequence":
            slice_sequence = [seq.upper()[i:i + 3] for i in range(0, len(seq), 3)]
            return set(slice_sequence).issubset(type(self).three_letter_alphabet)
        else:
            raise ValueError(
                f"Please provide amino acid sequence or check that the entered parameter record_type is correct")

    def __len__(self) -> int:
        """Return the length of the sequence"""
        if self.record_type == "one_letter_sequence":
            return super().__len__()
        if self.record_type == "three_letter_sequence":
            return int(len(self.seq) / 3)

    def __getitem__(self, key: Union[int, slice]) -> object:
        """Returns a subsequence in the form of a single letter or a string of several letters"""
        if self.record_type == "one_letter_sequence":
            return super().__getitem__(key)

        if self.record_type == "three_letter_sequence":
            if isinstance(key, int):
                if 0 >= key or key > len(self.seq):
                    raise IndexError(f"Index {key} is out of range of the sequence")
                return super().__getitem__(slice((key * 3) - 2, key * 3))

            elif isinstance(key, slice):
                try:
                    if key.start <= 0 or key.stop > len(self.seq) / 3:
                        raise IndexError(f"Index {key} is out of range of the sequence")
                    new_slice = slice((key.start * 3) - 2, key.stop * 3)
                    return super().__getitem__(new_slice)

                except TypeError:
                    if key.stop is None and key.start is None:
                        return self.seq
                    if key.start is None:
                        return super().__getitem__(slice(1, key.stop * 3))
                    if key.stop is None:
                        return super().__getitem__(slice((key.start * 3) - 2, int(len(self.seq) / 3) * 3))

    def count_mol_weight(self) -> int:
        """
        Calculate the molecular weight of the protein
        returns: the molecular weight of the protein (int)
        """

        molecular_weight = 0

        if self.record_type == "one_letter_sequence":
            for amino_acid in self.seq.upper():
                molecular_weight += type(self).amino_acid_weights[amino_acid]
        elif self.record_type == "three_letter_sequence":
            for i in range(0, len(self.seq), 3):
                molecular_weight += type(self).amino_acid_weights[self.seq[i:i + 3].upper()]
        return molecular_weight


def telegram_logger(chat_id):
    """
    The decorator function allows you to send a report via the telegram bot
    that the decorated function has completed its work
    :param chat_id: chat_id
    """

    def decorator(func):
        def wrapper(*args, **kwargs):
            old_stdout = sys.stdout
            old_stderr = sys.stdout

            file = StringIO()
            sys.stdout = file
            sys.stderr = file

            start_time = datetime.datetime.now()

            try:
                result = func(*args, **kwargs)
                end_time = datetime.datetime.now()
                execution_time = end_time - start_time
                message = f"\U0001F913 Function `{func.__name__}` successfully finished in `{execution_time}`"

            except RuntimeError as e:
                message = f"\U0001F644 Function `{func.__name__}` failed with an expection: `{type(e).__name__}: {str(e)}`"

            sys.stdout = old_stdout
            sys.stderr = old_stderr

            load_dotenv()
            TOKEN = os.getenv("TOKEN")

            if not file.getvalue():
                url = f"https://api.telegram.org/bot{TOKEN}/sendMessage?"
                data = {
                    "chat_id": chat_id,
                    "text": message,
                    "parse_mode": "Markdown",
                }
                requests.post(url, data=data)

            else:
                url = f"https://api.telegram.org/bot{TOKEN}/sendDocument?"
                data = {
                    "chat_id": chat_id,
                    "caption": message,
                    "parse_mode": "Markdown",
                }

                file.seek(0)
                requests.post(
                    url, data=data, files={"document": (f"{func.__name__}.log", file)}
                )
            return result

        return wrapper

    return decorator


@dataclass
class Intron:
    """
    Data class for intron
    :param number: sequence number
    :param start: starting coordinates
    :param stop: final coordinates
    :param gene_number: belonging to the gene (number)
    """
    number: int
    start: int
    stop: int
    gene_number: int

    def __repr__(self):
        header = f"{BLUE}{BOLD}Intron {self.number} in gene {self.gene_number}:{START_SETTINGS}\nStart: {self.start}\nEnd: {self.stop}\n"
        return f"{header}"


@dataclass
class Exon:
    """
    Data class for exon
    :param number: sequence number
    :param start: starting coordinates
    :param stop: final coordinates
    :param gene_number: belonging to the gene (number)
    """
    number: int
    start: int
    stop: int
    gene_number: int
    type_exon: str

    def __repr__(self):
        header = f"{BLUE}{BOLD}Exon {self.number} in gene {self.gene_number}:{START_SETTINGS}\nStart: {self.start}\nEnd: {self.stop}\nType: {self.type_exon}\n"
        return f"{header}"


@dataclass
class GenscanOutput:
    """
    Data class for genscan output
    :param status: status code
    :param cds_list: list of predicted protein sequences
    :param intron_list: list of found introns
    :param exon_list: list of found exons
    """
    status: str
    cds_list: list[SeqRecord]
    intron_list: list[Intron]
    exon_list: list[Exon]

    def __repr__(self):
        header = f"{BLUE}{BOLD}GenscanOutput:{START_SETTINGS}\nstatus: {self.status}\nnumber cds: {len(self.cds_list)}\nnumber intron: {len(self.intron_list)}\nnumber exon: {len(self.exon_list)}"
        return header


def run_genscan(
        sequence: str = None,
        sequence_file: str = None,
        organism: str = "Vertebrate",
        exon_cutoff: float = 1.00,
        sequence_name: str = "",
) -> GenscanOutput:
    """
    Sends a GenScan request and outputs the results in the format GenscanOutput
    :param sequence: DNA sequence in string format
    :param sequence_file: the path to the file with the DNA sequence
    :param organism: select 'Vertebrate'(default), 'Arabidopsis ' or 'Maize'
    :param exon_cutoff: select 0.01 ,0.02 ,0.05, 0.10, 0.25, 0.50, 1.00 (default)
    :param sequence_name: sequence name
    :return: GenscanOutput
    """

    seq = ""
    cds_list = []
    exon_list = []
    intron_list = []
    url = "http://hollywood.mit.edu//cgi-bin/genscanw_py.cgi"

    form_data = {
        "-o": organism,
        "-e": exon_cutoff,
        "-n": sequence_name,
        "-s": sequence,
        "-p": "Predicted peptides only",
    }

    if sequence_file is None:
        response = requests.post(url, data=form_data)
    else:
        with open(sequence_file, "rb") as file:
            response = requests.post(url, data=form_data, files={"-u": file})

    list_output = BeautifulSoup(response.content, "html.parser").text.split("\n")

    index_start_pred = list_output.index("Predicted genes/exons:") + 10
    index_end_pred = [
        index
        for index in range(len(list_output))
        if list_output[index].startswith("Suboptimal exons with probability")
    ]
    list_exon_line = list_output[index_start_pred: index_end_pred[0] - 3]

    index_start_seq = list_output.index("Predicted peptide sequence(s):")
    list_seq = list_output[index_start_seq + 1: -5]

    list_exon_line, list_seq = list(filter(None, list_exon_line)), list(
        filter(None, list_seq)
    )

    for i in range(len(list_exon_line)):
        if list_exon_line[0].startswith("NO EXONS"):
            return f"No exons have been detected"
        else:
            report_list = list_exon_line[i].split()
            gene_number = report_list[0].split(".")[0]
            exon_number = report_list[0].split(".")[1]
            exon_start = min(int(report_list[3]), int(report_list[4]))
            exon_end = max(int(report_list[3]), int(report_list[4]))
            exon_type = report_list[1]
            exon = Exon(exon_number, exon_start, exon_end, gene_number, exon_type)
            exon_list.append(exon)

            if i < len(list_exon_line) - 1:
                intron_start = exon_end + 1
                intron_end = (
                        min(
                            int(list_exon_line[i + 1].split()[3]),
                            int(list_exon_line[i + 1].split()[4]),
                        )
                        - 1
                )
                intron = Intron(exon_number, intron_start, intron_end, gene_number)
                intron_list.append(intron)

    for line in list_seq:
        if line.startswith(">"):
            if seq != "":
                seq_rec = SeqRecord(Seq(seq), id=id_seq)
                cds_list.append(seq_rec)
                seq += line
            id_seq = line
            seq = ""
        else:
            seq += line

    seq_rec = SeqRecord(Seq(seq), id=id_seq)
    cds_list.append(seq_rec)

    result = GenscanOutput(response.status_code, cds_list, intron_list, exon_list)
    return result


class MeasureTime:
    """Context manager for measuring time"""

    def __enter__(self):
        self.start_time = time.time()

    def __exit__(self, exc_type, exc_value, traceback):
        self.end_time = time.time()
        self.execution_time = self.end_time - self.start_time
        print(f"{BLUE}{BOLD}Время выполнения программы:{START_SETTINGS}{self.execution_time} сек\n"
              f"{GREEN}{BOLD}Старт:{START_SETTINGS}{time.ctime(self.start_time)}\n"
              f"{GREEN}{BOLD}Конец:{START_SETTINGS}{time.ctime(self.end_time)}\n")
