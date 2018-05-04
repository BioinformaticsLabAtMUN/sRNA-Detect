# sRNA-Detect

sRNA-Detect detects small RNAs from bacterial RNA-seq data.

## Requirements

To run sRNA-Detect, the [HTSeq  python library](http://htseq.readthedocs.io/en/release_0.9.1/) (version 0.5.4p5) and Python 2.7 need to be installed. 

## Running sRNA-Detect

After you have installed Python 2.7 and HTSeq, you can simply run sRNA-Detect in command-line as shown in the example below:

```
python2.7 sRNADetect.py  -samFile "file1.sam"  -samFile "file2.sam" -samFile "file3.sam" -out "output_sRNAs.gtf"
```

## Script Usage

Executing sRNADetect without parameters prints the following help message:

```
usage: sRNADetect.py [options] -samFile "sam_file1" -samFile "sam_file2"

This script takes alignment files in SAM format and detect small transcripts
with uniform coverage and no gaps

optional arguments:
  -h, --help          show this help message and exit
  -samFile SAMFILE    REQUIRED: sam file to read.
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
## Citing

If you use sRNA-Detect, please cite the following article:

[L. Peña-Castillo, M. Grüll, M.E. Mulligan and A.S. Lang, Detection of bacterial small transcripts from RNA-Seq data: a comparative assessment. Biocomputing 2016: pp. 456-467.](https://doi.org/10.1142/9789814749411_0042)

