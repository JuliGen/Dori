# 🧬DORi: <img src ='https://clipart-library.com/new_gallery/339-3399184_finding-nemo-facebook-banner.png' width =400 align="right">
## Determine the characteristics of amino acid sequences, Obtaine processed DNA or RNA, Remove low-quality fastq sequences
The idea of this program is to make life easier for experimenters and bioinformatics working with nucleic acids, amino acid sequences (both long and short), as well as fastq sequences. This tool is designed to work with nucleotide, amino acid sequences and fastq sequences.
The tool includes **3 modules:**
- `run_dna_rna_tools` — this module is designed to work with sequences of nucleic acids (DNA and RNA)
- `run_aminoacid_seq` — this module is designed to work with amino acid sequences
- `filter_read` - this module is designed to filter a dictionary with fastq sequences by specified parameters
## Instruction: 
### `run_dna_rna_tools`🧑🏻‍🔬 

### Usage:
The `run_dna_rna_tools` function accepts an arbitrary number of arguments with DNA or RNA sequences (*str*) as input. You also need to specify the argument `action` , which indicates the function to be performed with the sequences (the list below).

### List of actions:
- `transcribe` — conversion of a DNA sequence to a transcribed sequence
- `reverse` — conversion of the original direct sequence to the reverse (inverted)
- `complement` — conversion of the original sequence to a complementary
- `reverse_complement` — conversion of the original sequence to the reverse complementary

#### Example:
```python
run_dna_rna_tools('ATGTTT', action = 'transcribe')
run_dna_rna_tools('ATG', action = 'reverse')
run_dna_rna_tools('AtG', action = 'complement')
run_dna_rna_tools('ATg', action = 'reverse_complement')
run_dna_rna_tools('ATG', 'aT', action = 'reverse')
```

### Module Features:
- The module in the processed sequence saves the character case that was in the original sequences.
- When calling the function only to sequences that are not related to DNA or RNA, the program politely ignores your request and asks you to maintain DNA or RNA sequences.
- At the same time, if sequences related to nucleic acids and sequences of arbitrary letters are simultaneously submitted to the input, then DNA and RNA sequences will be processed.At the same time, if sequences related to nucleic acids and sequences of arbitrary letters are simultaneously submitted to the input, then DNA and RNA sequences will be processed. At the output you will get processed sequences, as well as the number (number in the input list) and sequences that were not processed.

```python
run_dna_rna_tools('TG', 'DF' action = 'reverse')

#Output
#Sequences that are not DNA or RNA (number in the input list: seq): {2: 'DF'}
#Result:
#GT
```
- When calling the transcribe procedure to RNA sequences, the program will display an error message.

```python
run_dna_rna_tools('AUU', action = 'transcribe')

#Output
#ValueError: Provide DNA, RNA sequences cannot be transcribed
```
### `run_aminoacid_seq`👻

### Usage:
The function `run_aminoacid_seq` takes 1 amino acid sequence as input (both in three-letter and one-letter form, str), next you need to specify named arguments: `function`, `record_type`.
- The first argument must specify the amino acid sequence.  The sequence must be written entirely in uppercase or lowercase letters without spaces or commas. The output preserves the case.
- The argument `function` must be passed a string with the name of the action (see possible actions below) that needs to be performed. Default summary.
- The argument `record_type` indicates the form in which you present your sequence. If you are using a three-letter sequence, specify `record_type= 3`. If you are using a one -letter sequence, specify `record_type= 1` (by default). 

#### Example:

```python
run_aminoacid_seq('ALAGLNGLU', function = 'count', record_type = 3)
run_aminoacid_seq('glyvalala', function = 'count', record_type = 3)
run_aminoacid_seq('ASL', function = 'count')
run_aminoacid_seq('alaglyala', function = 'determine_charge', record_type = 3)
run_aminoacid_seq('LLYdD', function = 'determine_charge', record_type = 1, percent=True)
run_aminoacid_seq('alaglyala', function = 'determine_charge', record_type = 3, percent=True)
run_aminoacid_seq('ALAGLYALA', function = 'translate', record_type = 3)
```          
   
Also, if necessary specify a named argument `percent=True` (default False) for actions: determine_charge, determine_polarity (Look in the description of functions).

#### Example:
```python
run_aminoacid_seq('LLYdD', function = 'determine_charge', record_type = 1, percent=True)

#Output
#{'Percentage of positively charged amino acids': 0, 'Percentage of neutrally charged amino acids': 60, 'Percentage of negatively charged amino acids': 40}
```   
### List of functions:
- `translate` - translation of a one-letter amino acid sequence into a three-letter one (for better visual perception), and the reverse operation. 
- `count` - obtaining the length of the amino acid sequence. 
- `count_molecular_weight` - calculating the molecular weight of a protein.
- `determine_charge`- counting the number of positive, negative and neutral amino acids in a protein. To get the output in percent, specify percent=True. 
- `determine_polarity` - counting hybrophobic and hydrophilic amino acids in a protein. To get the output in percent, specify percent=True. 
- `convert_amino_acid_seq_to_dna` - convert an amino acid sequence to the most likely DNA sequence. 
- `summary` - a summary of all information about the sequence (the result of executing all functions). 

## Troubleshooting:
The program works with 20 proteinogenic amino acids {'G', 'A', 'V', 'L', 'I', 'M', 'P', 'F', 'W', 'S', 'T', 'N', 'Q', 'Y', 'C', 'K', 'R', 'H', 'D', 'E'}. If you use a symbol that does not represent an amino acid, or wrong record type, the program will generate an **error** where you can see the first wrong symbol. 
- correct function launch:    
```python
run_aminoacid_seq('ALALEUILE', function = 'count', record_type = 3)
```  
- incorrect function launch:      
```python
run_aminoacid_seq('ALXZEQ', function = 'count', record_type = 3)
run_aminoacid_seq('GluGlu', function = 'count', record_type = 3)
run_aminoacid_seq('ALALEUILE', function = 'count', record_type = 1)
```
### `filter_read`🦖
The `filter_read` function is designed to filter a fastq file by specified parameters (read length, gc content, quality).

### Usage:
The function `filter_read` accepts 5 arguments as input: input_path (str), output_filename (str), gc_bounds (int, float or tuple), length_bounds (int, float или tuple), quality_threshold (int).

- As the first argument, specify the path to the fastq file.
- In argument `output_filename`, specify the name of the file with filtered sequences. If the argument is omitted, the file will be named the same as the source file. 
- `gc_bounds` - GC content interval (in percent) for filtering (default - (0, 100), i.e. all reads are saved). If a single number is passed to the argument, it is assumed that this is the upper bound.
Examples: `gc_bounds = (20, 80)` - save only reads with GC content from 20 to 80%, `gc_bounds = 50.5` - save reads with GC content less than 50.5%.
- `length_bounds` - the length interval for filtering (default - (0, 2**32)).
Examples: `length_bounds = (0, 250)` - sequences from 0 to 250 in length will be saved, `length_bounds = 250` - sequences from 0 to 250 in length will be saved.
- `quality_threshold` - the threshold value of the average quality of the read for filtering (default - 0 (phred33 scale)). Reads with average quality for all nucleotides below the threshold are discarded.
Examples: `quality_threshold = 30` - reads with average quality for all nucleotides below 30 will be discarded.
- The resulting file will be saved to the **"fastq_filtrator_resuls"** folder.

#### Example:

```python
filter_read(input_path='example_fastq.fastq', length_bounds= 21)
filter_read(input_path='example_fastq.fastq', gc_bounds=(20,30), quality_threshold=20)

# Output
# example_fastq.fastq file created!


filter_read(input_path='example_fastq.fastq', output_filename='filt_fastq', quality_threshold=30)
filter_read(input_path='example_fastq.fastq', output_filename='filt_fastq', quality_threshold=30, length_bounds=200)

# Output
# example_fastq.fastq file created!

```   
- If you want to get a report on the number of filtered sequences, add the argument `report=True` (default - `report=False`).

```python
filter_read(input_path='example_fastq.fastq', output_filename='filt_fastq', quality_threshold=30, report=True)
``` 

# 🧬bio_files_processor: 
This tool is designed to work with bioinformatic files, in particular with fasta, gbk, txt (blast output).
The tool includes ** 4 modules:**
- `convert_multiline_fasta_to_oneline` — this module is designed to translate a multi-line sequence entry in a fasta file into a single-line view.
- `select_genes_from_gbk_to_fasta` — this module is designed to isolate a certain number of genes before and after the gene of interest and save their protein sequence (translation) to a fasta file.
-  `change_fasta_start_pos` - this module is designed to shift the starting position in the sequence and write a sequence file with a shifted start.
-  `parse_blast_output` - this module is designed to work with files received after blast processing. The module allows you to get the best match with the database for each sequence that has been analyzed in blast. The output is a file with a list of protein names that are closest to the sequence that was analyzed in blast.

## Instruction: 
### `convert_multiline_fasta_to_oneline`🧑🏻‍🔬

### Usage:
The function `convert_multiline_fasta_to_oneline` takes 2 arguments as input: `input_fasta` (str)  and `output_fasta` (str). 
- As an argument to `input_fasta`, you must pass the path to the fasta file, which contains sequences written in a multiline version.
- The name of the processed file needs to be passed as the `output_fasta` argument. If the argument is not passed to the input, then the processed file will be named the same as the original one with the addition of "_oneline". The extension ".fasta" will be added to the file name.

*At the output, you will receive a fasta file with sequences in a single-line version.*

#### Example:

```python
convert_multiline_fasta_to_oneline(input_fasta='example_multiline_fasta.fasta')

#Output
#file example_multiline_fasta_oneline.fasta

convert_multiline_fasta_to_oneline(input_fasta='example_multiline_fasta.fasta',output_fasta='example_oneline' )

#Output
#file example_oneline.fasta
```

### `select_genes_from_gbk_to_fasta` 🤹🏻‍♂️

### Usage:
The function `select_genes_from_gbk_to_fasta` accepts the following arguments as input:
- `input_gbk`(str) - path to the input GBK file.
- `genes` (list) - genes of interest, next to which neighbors are being searched.  It is necessary to give a list, even if there is only one gene.
- `n_before` (int) - number of genes up to the gene of interest that need to be included in the final file. By default 1.
- `n_after` (int) - number of genes after the gene of interest that need to be included in the final file.By default 1.
- `output_fasta` (str) - name of the output file. If the argument `output_fasta` is not passed to the input, then the processed files will be named the same as the original ones, with the addition of "_select_genes_trans". The extension ".fasta" will also be added to the file name".

*At the output you will receive a file  with protein sequences of genes that were located before and after the genes you are interested in. The resulting file can be used for blast analysis to get data about the environment of the genes you are interested in.*

#### Example:

```python
select_genes_from_gbk_to_fasta(input_gbk='example_gbk.gbk', genes=['dtpD','pxpC'])
select_genes_from_gbk_to_fasta(input_gbk='example_gbk.gbk', genes=['dtpD','pxpC'], n_before=2, n_after=2)
select_genes_from_gbk_to_fasta(input_gbk='example_gbk.gbk', genes=['dtpD'], n_before=2)

#Output
#Created example_gbk_select_genes_trans.fasta

select_genes_from_gbk_to_fasta(input_gbk='example_gbk.gbk', genes=['dtpD'], n_before=2, n_after=2, output_fasta='select_protein')
select_genes_from_gbk_to_fasta(input_gbk='example_gbk.gbk', genes=['dtpD'], n_before=2,  output_fasta='select_protein')
select_genes_from_gbk_to_fasta(input_gbk='example_gbk.gbk', genes=['dtpD','pxpC'], output_fasta='select_protein')

#Output
#Created select_protein.fasta
```

### `change_fasta_start_pos` 🥷🏻

### Usage:
The function `change_fasta_start_pos` accepts the following arguments as input:
- `input_fasta`(str)- path to the input fasta file. A fasta file with a single sequence in a single-line version is accepted as input. If a file with several sequences comes to the input, only the first one will be processed. The multiline fasta format will be processed incorrectly, so if your data is like this, first use the convert_multiline_fasta_to_oneline function.
- `shift`(int) - how much you need to shift the initial position in the file (a positive or negative integer).
- `output_fasta`(str) - name of the output file. If the argument `output_fasta` is not passed to the input, then the processed files will be named the same as the original ones, with the addition of "_select_genes_trans". The extension ".fasta" will also be added to the file name". 

*At the output, you will receive a file with a sequence that has a displaced starting point.*

#### Example:

- correct function launch:

```python
change_fasta_start_pos(input_fasta='example_oneline.fasta', shift=2)
change_fasta_start_pos(input_fasta='example_multiline_fasta.fasta', shift=-3)

#Output
#Created example_oneline_shifted.fasta

change_fasta_start_pos(input_fasta='example_multiline_fasta.fasta', shift=-3, output_fasta='shifted_fasta')
change_fasta_start_pos(input_fasta='example_multiline_fasta.fasta', shift=2, output_fasta='shifted_fasta' )

#Output
#Created shifted_fasta.fasta

```

- incorrect function launch:
  
```python
change_fasta_start_pos(input_fasta='example_multiline_fasta.fasta', output_fasta='shifted_fasta' )

#Output
#ValueError: Please provide a parameter "shift"

```
### `parse_blast_output`🦇

### Usage:
The function `parse_blast_output` takes 2 arguments as input: `input_file` (str)  and `output_file` (str).
- `input_file` - path to the input data file.
- `output_file` - name of the output file. If the argument `output_file` is not passed to the input, then the processed files will be named the same as the original ones, with the addition of "_select_protein".

*At the output, you will receive a txt file with proteins that have the best matches with the database.*

```python
parse_blast_output(input_file='example_blast_results.txt')

#Output
#file example_blast_results_select_protein.txt

parse_blast_output(input_file='example_blast_results.txt', output_file='select_protein')

#Output
#select_protein.txt
```
