# 🧬DORi: <img src ='https://clipart-library.com/new_gallery/339-3399184_finding-nemo-facebook-banner.png' width =400 align="right">
## Determine the characteristics of amino acid sequences, Obtaine processed DNA or RNA, Remove low-quality fastq sequences
**Dori** contains tools that were written while doing home tasks in Python course in **Bioinformatics Institute** (OOP, API, testing, multiprocessing).

## 🐠dori.py:
The idea of this module is to make life easier for experimenters and bioinformatics working with nucleic acids, amino acid sequences (both long and short), as well as fastq sequences. 

**The module includes:**

`RNASequence`/`DNASequence`/`AminoAcidSequence`/`BiologicalSequence`/`NucleicAcidSequence` - designed to manipulate with three main bioinformatics data types - DNA, RNA and proteins.
  Methods:
   - `complement` -  returns a complementary sequence
   - `reverse` - returns a reversed sequence
   - `gc_content` - calculate GC-content of DNAsequence or RNAsequence
   - `transcribe` - returns a transcribed sequence (for DNAsequence only)
   - `count_mol_weight` - calculate count_mol_weight of AminoAcidSequence

`filter_read` - this function is designed to filter a fastq file by specified parameters(read length, gc content, quality)

`telegram_logger` decorator - a decorator that allows you to send information about the function of interest via telegram bot. The information includes the result of run, time of running, stdout and stderr output in logfile. This function is implemented based on [Telegram API](https://core.telegram.org)  without specific libraries.

`run_genscan` - makes queries on the [Genscan website](http://hollywood.mit.edu/GENSCAN.html) to predict possible CDS, exons, and introns in the sequence of interest (input: file or string). Output: object GenscanOutput. Also, `Intron` and `Exon` classes are implemented for more convenient use of the received data.

`MeasureTime` - context manager for measuring time

## 🧬bio_files_processor: 

This module is designed to work with bioinformatic files, in particular with fasta, gbk, blast output.
- `convert_multiline_fasta_to_oneline` — this module is designed to translate a multi-line sequence entry in a fasta file into a single-line view.
- `select_genes_from_gbk_to_fasta` — this module is designed to isolate a certain number of genes before and after the gene of interest and save their protein sequence (translation) to a fasta file.
-  `change_fasta_start_pos` - this module is designed to shift the starting position in the sequence and write a sequence file with a shifted start.
-  `parse_blast_output` - this module is designed to work with files received after blast processing. The module allows you to get the best match with the database for each sequence that has been analyzed in blast. The output is a file with a list of protein names that are closest to the sequence that was analyzed in blast.
-  `OpenFasta` - context manager for iterating on the fasta file. To open fasta file and return individual FASTA records including id, description and sequence. The implementation of the OpenFasta context manager is similar to the built-in open context manager.
- `FastaRecord` - dataclass for storing Fasta

## 🌲custom_random_forest.py: 
- `RandomForestClassifierCustom` - this class represents a custom implementation of the classifier using the random forest algorithm. Additionally, this class supports the use of multiple threads (the n_jobs parameter) to speed up the fitting and prediction process.

## 🛠️test_dori.py: 
The `script test_dori.py` it includes functions for testing the tools presented above.



## Example:
<img src ='https://i.pinimg.com/originals/10/3d/d7/103dd740fe24055d2db8a69056c70fb9.png' width =200 align="right">
Examples of running some programs are presented in Showcases.ipynb.


## 

**I would like to express my deep gratitude to the team of the Institute of Bioinformatics. 
Because thanks to them, this repository was born and began to grow up**❤️





