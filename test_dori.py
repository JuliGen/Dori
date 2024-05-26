import os

import pytest

from bio_files_processor import OpenFasta, parse_blast_output
from dori import AminoAcidSequence, DNASequence, NucleicAcidSequence


@pytest.fixture
def input_data_nucl_seq():
    nucl_seq_input = NucleicAcidSequence('AATGCC')
    return nucl_seq_input


def test_complement_nucl_seq(input_data_nucl_seq):
    with pytest.raises(NotImplementedError):
        input_data_nucl_seq.complement()


@pytest.fixture
def input_data_dna_seq():
    dna_seq_input = DNASequence('AATGCC')
    return dna_seq_input


def test_complement_dna_seq(input_data_dna_seq):
    target = DNASequence('TTACGG')
    result = input_data_dna_seq.complement()
    assert target == result


def test_check_alphabet_amino():
    with pytest.raises(ValueError):
        AminoAcidSequence('AAJLERER', record_type='one_letter_sequence')


def test_one_letter_alphabet_am_seq():
    target = set('GAVLIMPFWSTNQYCKRHDE')
    result = set(AminoAcidSequence.one_letter_alphabet)
    assert target == result


def test_three_letter_alphabet_am_seq():
    target = {'GLY', 'ALA', 'VAL', 'LEU', 'ILE', 'MET', 'PRO', 'PHE', 'TRP', 'SER', 'THR', 'ASN', 'GLN',
              'TYR', 'CYS', 'LYS', 'ARG', 'HIS', 'ASP', 'GLU'
              }
    result = set(AminoAcidSequence.three_letter_alphabet)
    assert target == result


@pytest.fixture
def input_fasta_record():
    return ('WP_096906625.1', 'dioxygenase AlkB', 'MTDPLFGVDRAPT')


@pytest.fixture
def input_fasta_file():
    file_path = './data/sequence.fasta'
    return file_path


def test_open_fasta_content(input_fasta_file, input_fasta_record):
    with OpenFasta(input_fasta_file) as file:
        fasta_record = file.read_record()
        id_, description = fasta_record.id, fasta_record.description
        sequence = fasta_record.seq
    result = (id_, description, sequence)
    assert input_fasta_record == result


@pytest.fixture
def blast_output_file():
    file_path = './data/blast_output.txt'
    return file_path


@pytest.fixture
def tmp_file():
    file_path = './tmp_file.txt'
    yield file_path
    if os.path.exists(file_path):
        os.remove(file_path)


def test_parse_blast_output_exists(blast_output_file, tmp_file):
    parse_blast_output(blast_output_file, tmp_file)
    assert os.path.exists(tmp_file)


def test_parse_blast_output_content(blast_output_file, tmp_file):
    parse_blast_output(blast_output_file, tmp_file)
    with open(tmp_file) as file:
        line_1 = file.readline().strip()
        line_2 = file.readline().strip()
        assert line_1 == "PhtB [Rhodococcus sp. TFB]"
        assert line_2 == "alpha-ketoglutarate-dependent dioxygenase AlkB [Dietzia]"
