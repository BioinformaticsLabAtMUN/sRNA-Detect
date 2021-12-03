# sRNA-Detect

sRNA-Detect detects small RNAs from bacterial RNA-seq data.

## Requirements

To run sRNA-Detect, the [HTSeq  python library](https://htseq.readthedocs.io/en/latest/) (version 0.13.5) and Python 3.9 need to be installed. 

## Running sRNA-Detect

After you have installed Python 3 and HTSeq, you can simply run sRNA-Detect in command-line as shown in the example below:

```
python3 sRNADetect.py  -samFile "file1.sam"  -samFile "file2.sam" -samFile "file3.sam" -out "output_sRNAs.gtf"
```

## Script Usage

Executing sRNADetect without parameters prints the following help message:

```
usage: sRNADetect.py [options] -samFile "sam_file1" -samFile "sam_file2"

This script takes alignment files in SAM/BAM format and detect small transcripts
with uniform coverage and no gaps

optional arguments:
  -h, --help          show this help message and exit
  -samFile SAMFILE    REQUIRED: sam/bam file to read.
  -minLen MINLEN      minimum length to consider for small transcripts
  -maxLen MAXLEN      maximum length to consider for small transcripts
  -minH MINH          minimum number of reads across all samples required to
                      detect small transcripts
  -maxPDV MAXPDV      maximum allowed drop in coverage wrt current coverage
  -maxPN MAXPN        maximum allowed negative change (percentage) in mean
                      coverage
  -maxPP MAXPP        maximum allowed positive change (percentage) in mean
                      coverage
  -idPrefix IDPREFIX  Prefix to be used as identifier for sRNAs in gtf file
  -out OUT            Filename of the output gtf file
  -v, --verbose
  --version           show program's version number and exit

Written by Lourdes Pena-Castillo (lourdes at mun dot ca), Memorial University of
Newfoundland. Copyright 2015. "Release under the terms of the GNU General
Public License".
```

## sRNA-Detect Results

sRNA-Detect will create a GTF file with the small transcripts detected as expressed. Below are a few sample lines of a GTF file created by sRNA-Detect.

```
gi|292897734|ref|NC_013971.1|	sRNADetect	sRNA	67	98	5669.9737469	+	.	gene_id "1"; unique_id "sRNA00001"
gi|292897734|ref|NC_013971.1|	sRNADetect	sRNA	122	187	91.0714550309	+	.	gene_id "2"; unique_id "sRNA00002"
gi|292897734|ref|NC_013971.1|	sRNADetect	sRNA	376	425	2632.80376766	+	.	gene_id "3"; unique_id "sRNA00003"
gi|292897734|ref|NC_013971.1|	sRNADetect	sRNA	506	551	26.6309028735	+	.	gene_id "4"; unique_id "sRNA00004"
gi|292897734|ref|NC_013971.1|	sRNADetect	sRNA	598	648	62.1329525246	+	.	gene_id "5"; unique_id "sRNA00005"
```

First column is the sequence name, fourth and fifth columns are the start and end position of the small transcript, respectively; sixth column is the average read depth coverage of the small transcript, seventh column is the strand, and nineth column is the identifier generated by sRNA-Detect.

## Filtering Detected Transcripts

If you are looking for novel sRNAs, we recommend to filter out the transcripts detected by sRNA-Detect by removing those that overlap annotated transcripts. To do that, one can use [BedTools' intersectBed](http://bedtools.readthedocs.io/en/latest/content/tools/intersect.html) as follows:

```
bedtools intersect -wo -a output_sRNAs.gtf -b annotatedTranscripts.bed -sorted -v -s -f 0.8
```
This command won't remove putative asRNAs or putative sRNAs that partially overlap a larger annotated transcript.

## Citing

If you use sRNA-Detect, please cite the following article:

[L. Peña-Castillo, M. Grüll, M.E. Mulligan and A.S. Lang, Detection of bacterial small transcripts from RNA-Seq data: a comparative assessment. Biocomputing 2016: pp. 456-467.](https://doi.org/10.1142/9789814749411_0042)

