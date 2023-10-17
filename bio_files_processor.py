import os
import additional_modules.parse_gbk as parse_gbk


def convert_multiline_fasta_to_oneline(input_fasta: str, output_fasta: str = None) -> None:
    """
    Convert multiline fasta format to single line
    :param input_fasta: filename fasta or path to file
    :param output_fasta: name of the processed fasta file
    :return: None
    """
    names = []
    seqs = []
    seq = ''

    with open(input_fasta) as multifasta:
        for line in multifasta:
            line = line.strip()
            if line.startswith('>'):
                if seq != '':
                    seqs.append(seq)
                    names.append(name)
                name = line
                seq = ''
            else:
                seq += line
        seqs.append(seq)
        names.append(name)

    if output_fasta is None:
        output_fasta = f'{os.path.splitext(os.path.basename(input_fasta))[0]}_oneline'

    output_fasta = f'{output_fasta}.fasta'

    with open(output_fasta, mode='w') as file:
        for name, seq in zip(names, seqs):
            file.write(name + '\n')
            file.write(seq + '\n')
    return


def select_genes_from_gbk_to_fasta(input_gbk: str, genes: list = None, n_before: int = 1,
                                   n_after: int = 1, output_fasta: str = None) -> str:

    """
    Write to a file the translation of genes found before and after the gene of interest.

    :param input_gbk: path to gbk file
    :param genes: the desired genes whose neighbors need to be found (list of strings)
    :param n_before: number of genes up to the desired gene that must be included in the file
    :param n_after: number of genes after the desired gene that must be included in the file
    :param output_fasta: output file name
    :return: fasta file
    """
    if output_fasta is None:
        output_fasta = f'{os.path.splitext(os.path.basename(input_gbk))[0]}_select_genes_trans'
    output_fasta = f'{output_fasta}.fasta'

    if genes is None:
        raise ValueError('Please provide the names of the genes whose neighbors need to be found')

    genes_before = []
    genes_after = []

    with open(input_gbk) as gbk:
        lines = gbk.readlines()
        gene_or_locus_tag = parse_gbk.extract_genes(lines)
        translation = parse_gbk.extract_translation(lines)
        dict_gene_or_locus_tag_trans = dict(zip(gene_or_locus_tag, translation))

        for index_gene, gene in enumerate(gene_or_locus_tag):
            for select_gene in genes:
                if select_gene in gene:
                    if (index_gene - n_before) < 0:
                        genes_before.extend(gene_or_locus_tag[:index_gene])
                    if (index_gene + n_after) > len(gene_or_locus_tag):
                        genes_after.extend(gene_or_locus_tag[index_gene + 1: len(gene_or_locus_tag)])
                    else:
                        genes_before.extend(gene_or_locus_tag[index_gene - n_before:index_gene])
                        genes_after.extend(gene_or_locus_tag[index_gene + 1: index_gene + 1 + n_after])

    with open(output_fasta, mode='w') as file:
        for gene_before in genes_before:
            file.write('>' + gene_before + '\n')
            file.write(dict_gene_or_locus_tag_trans[gene_before] + '\n')
        for gene_after in genes_after:
            file.write('>' + gene_after + '\n')
            file.write(dict_gene_or_locus_tag_trans[gene_after] + '\n')

    return f'Created {output_fasta}'


def change_fasta_start_pos(input_fasta: str, shift: int = None, output_fasta: str = None) -> str:

    """
    Write a sequence with a shifted start to a file (by the number of nucleotides specified in the shift parameter)
    :param input_fasta: path to fasta file
    :param shift: the number of nucleotides by which the starting position in the sequence must be shifted
    :param output_fasta: name of the processed fasta file
    :return: str (name of the created file)
    """

    if shift is None:
        raise ValueError('Please provide a parameter "shift"')

    if output_fasta is None:
        output_fasta = f'{os.path.splitext(os.path.basename(input_fasta))[0]}_shifted'

    output_fasta = f'{output_fasta}.fasta'

    with open(input_fasta) as fasta:
        name_seq = fasta.readline().strip()
        seq = fasta.readline().strip()
        seq_shift = seq[shift: len(seq)] + seq[:shift]

    with open(output_fasta, mode='w') as output_file:
        output_file.write(name_seq + '\n')
        output_file.write(seq_shift + '\n')

    return f'Created {output_fasta}'


def parse_blast_output(input_file: str, output_file: str = None) -> None:

    """
    Select the name of the best match from the database for each amino acid sequence
    :param input_file: path to the txt file containing data from blast
    :param output_file: file name with selected protein names
    :return: None
    """

    if output_file is None:
        output_file = f'{os.path.splitext(os.path.basename(input_file))[0]}_select_protein'
    output_file = f'{output_file}.txt'

    list_proteins = []
    with open(input_file) as blast_result:
        for line in blast_result:
            if line.startswith('Description'):
                line = blast_result.readline()
                first_protein = line.replace('...', ']').split(']')[0] + ']'
                list_proteins.append(first_protein)

    list_proteins.sort()

    with open(output_file, mode='w') as output_file:
        for protein in list_proteins:
            output_file.write(protein + '\n')

    return
