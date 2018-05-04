This is a quick start guide to use sRNA-Detect.

---------------------------------------

To run sRNA-Detect,  the HTSeq  python library available at http://www-huber.embl.de/users/anders/HTSeq/doc/index.html and Python 2.7 need to be installed. After you have installed Python 2.7 and HTSeq, you can simply run sRNA-Detect in command-line as shown in the example below:

$ python2.7 sRNADetect.py  -samFile "file1.sam"  -samFile "file2.sam" -samFile "file3.sam" -out "output_sRNAs.gtf"

---------------------------------------

Executing sRNADetect without parameters prints the following help message:

$ python2.7 sRNADetect.py 
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

---------------------------------------

If you use sRNA-Detect, please cite sRNA-Detect as follows:

"L. Peña-Castillo, M. Grüll, M.E. Mulligan and A.S. Lang, Detection of bacterial small transcripts from RNA-Seq data: a comparative assessment. In the proceedings of the Pacific Symposium on Biocomputing (PSB) 2016."

---------------------------------------

sRNA-Detect is developed by Lourdes Peña-Castillo at Memorial University of Newfoundland.  Please do not hesitate to contact me (lourdes at  mun dot ca) if you have any questions or comments.

