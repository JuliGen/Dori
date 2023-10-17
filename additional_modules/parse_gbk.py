def extract_translation(lines: list) -> list:

    """
    export translated sequences from list of lines from gbk file

    :param lines: list of lines from gbk file
    :return: list with sequences
    """
    description_cds = ''
    translation = []
    for number_line, line in enumerate(lines):
        if line.startswith('                     '):
            description_cds += line.strip().replace(' ', '')

    description_cds_lines = description_cds.split('/')
    for line in description_cds_lines:
        if line.startswith('translation'):
            translation.append(line.lstrip('/translation=').replace('"', ''))
    return translation


def extract_genes(lines: list) -> list:

    """
    export gene or locus_tag from list of lines from gbk file

    :param lines: list of lines from gbk file
    :return: list with gene or locus_tag
    """

    gene_or_locus_tag = []

    for number_line, line in enumerate(lines):
        if line.startswith('     CDS             '):
            gene_or_locus_tag.append(
                lines[number_line + 1].rstrip('"\n').
                replace('                     /gene="', '').
                replace('                     /locus_tag="', '')
                )

    return gene_or_locus_tag

