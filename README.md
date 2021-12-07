# sRNA-Detect

A Pipeline for small RNA Detection (sRNA-Detect) written in the [Nextflow](http://nextflow.io).

This pipeline goes through the following processes:

* sRNA-Detect process which detects small RNAs from bacterial RNA-seq data.
* Sorting process which sort the sRNAs detected from sRNA-Detect.
* Filtering process which filters out the overlapping transcripts among detected sRNAs.
* Categorizing process which categorizes the overlapping transcripts filtered out from previous process.
##  1 Requisites
### 1.1 Nextflow
Make sure Java 8 or later is already installed on your machine.
```
java -version
```
Install [[nextflow]](http://nextflow.io) with the following command:

```
curl -s https://get.nextflow.io | bash
```
### 1.2 Other software
There are two options here. 
* Install Docker.
* Install native software and dependencies.
### 1.2.1 Docker
It is the most recommended option. 

Install [Docker](https://www.docker.com/) 18.03 (or higher) and then pull the docker image:
```
docker pull penacastillolab/sRNA-Detect
```
To get better understanding of the Docker image, please see the included [Dockerfile](Dockerfile).

### 1.2.2 Natively
If you are not using Docker, you must install the following softwares and packages:
* [Python 3](https://www.python.org/) version 3.8 (or higher)
* [HTSeq](https://htseq.readthedocs.io/en/master/) version 0.13.5 (or higher)
* [Bedtools](http://bedtools.readthedocs.io/en/latest/index.html) version 2.27 (or higher)

## 2 Pipeline usage
After getting the file [detect_filter_sRNA.nf](detect_filter_sRNA.nf), you can check how to use the pipeline by this command:
```
nextflow detect_filter_sRNA.nf --help
```

```
N E X T F L O W  ~  version 21.04.3
Launching `detect_filter_sRNA.nf` [condescending_bernard] - revision: 7b1b4280e0
sRNADetect: sRNA detect and filter pipeline
----------------------------------------------------
Options:
--alignmentDir  directory   the directory of SAM/BAM files [required]
--outputDir  directory   the directory for saving output files
--annotatedGenomeFile  path   the path to the given Genome file
--output  fileName    the name of the sRNA-Detect output gtf file
--idPrefix    prefix  Prefix to be used as identifier for sRNAs in gtf file
--minLength    integer  minimum length to consider for small transcripts
--maxLength    integer  maximum length to consider for small transcripts
--minHeight    integer  minimum number of reads across all samples required to detect small transcripts
--maxPctgDropValue    percentages  maximum allowed drop in coverage wrt current coverage
--maxPctgChangeNeg    maximum allowed negative change (percentage) in mean coverage
--maxPctgChangePos    maximum allowed positive change (percentage) in mean coverage
```
The annotated Genome file is a GFF file with nine columns: Sequence Name, Source, Feature, Start Position, End Position, Score, Strand, Phase, Attributes. For example,
```
Chromosome	RefSeq	region	1	3738958	.	+	.	ID=id0;Dbxref=taxon:272942;Is_circular=true;Name=ANONYMOUS;gbkey=Src;genome=chromosome;mol_type=genomic DNA;strain=SB 1003
Chromosome	RefSeq	gene	351	1724	.	+	.	ID=gene0;Name=RCAP_RS00005;gbkey=Gene;gene_biotype=protein_coding;locus_tag=RCAP_RS00005;old_locus_tag=RCAP_rcc00001
Chromosome	Protein Homology	CDS	351	1724	.	+	0	ID=cds0;Parent=gene0;Dbxref=Genbank:WP_013065748.1;Name=WP_013065748.1;gbkey=CDS;product=chromosomal replication initiator protein DnaA;protein_id=WP_013065748.1;transl_table=11
```

## 3 Running the Pipeline
* By default, sRNA-Detect will run into the [Docker](https://www.docker.com/) container (see [Nextflow config file](nextflow.config)). To use the default settings, you just need to keep your [Docker](https://www.docker.com/) open.
* If you have installed the required software natively, you can comment out the configurations in the [Nextflow config file](nextflow.config). This way, the pipeline will run natively.

You can use the [test data](test_data) to run sRNA-Detect. Go to the root directory of this project, and then run the following command on your terminal:
```
nextflow detect_filter_sRNA.nf --alignmentDir ./test_data --annotatedGenomeFile ./test_data/GCF_000021865.1_ASM2186v1_chromosome.gff --outputDir ./test_data --output sRNADetect_Rcap_output.gtf --idPrefix RCAP_rcs
```
On your terminal, you will see information like this:
```
N E X T F L O W  ~  version 21.04.3
Launching `detect_filter_sRNA.nf` [friendly_carlsson] - revision: 0e432bb127
executor >  local (4)
[ae/ce96e2] process > sRNA_Detect              [100%] 1 of 1 ✔
[bb/5290e0] process > sort_sRNA                [100%] 1 of 1 ✔
[ba/4183cc] process > filter_sorted_sRNA       [100%] 1 of 1 ✔
[49/d5d599] process > separate_overlapped_sRNA [100%] 1 of 1 ✔
Completed at: 07-Dec-2021 13:27:35
Duration    : 1m 53s
CPU hours   : (a few seconds)
Succeeded   : 4
```

## 4 Pipeline Results
The results are saved into the [test data](test_data) directory. 
The following are descriptions of these results:
* ``sRNADetect_Rcap_output.gtf`` -The output from sRNA-Detect process. 
It is a GTF file with the small transcripts detected. 
First column is the sequence name, fourth and fifth columns are the start and end position of the small transcript, respectively; sixth column is the average read depth coverage of the small transcript, seventh column is the strand, and nineth column is the identifier generated by sRNA-Detect.
Below are a few sample lines of the GTF file.
```
pRCB133	sRNADetect	sRNA	72188	72270	10.983726768820178	-	.	gene_id "6"; unique_id "RCAP_rcs00006"
pRCB133	sRNADetect	sRNA	99505	99557	169.09995195961872	-	.	gene_id "7"; unique_id "RCAP_rcs00007"
pRCB133	sRNADetect	sRNA	104684	104798	10.434254695881878	-	.	gene_id "8"; unique_id "RCAP_rcs00008"
Chromosome	sRNADetect	sRNA	91	242	25.91696470201933	+	.	gene_id "9"; unique_id "RCAP_rcs00009"
Chromosome	sRNADetect	sRNA	246	269	12.543877454363388	+	.	gene_id "10"; unique_id "RCAP_rcs00010"
```
* ``non-overlap.gff`` -The output from Filtering process. It is in the same format as ``sRNADetect_Rcap_output.gtf``.
But it contains only novel sRNAs, which are filtered by removing those that overlap with annotated transcripts. 
* ``overlap.gff`` -The output from Filtering process. It is in the same format as ``sRNADetect_Rcap_output.gtf``.
But it contains only transcripts that overlap with annotated transcripts. In the file, each line represents an overlap. 
The left part is a detected transcript, and the middle part is the annotated transcript it overlaps with. The right part is the number of overlapping
bases and the overlap percentage.
```
Chromosome sRNADetect sRNA 35149 35234 12.801680672268908 - . gene_id "755"; unique_id "RCAP_sre00755" Chromosome RefSeq gene 34143 35957 . - . ID=gene26;Name=RCAP_RS00135;gbkey=Gene;gene_biotype=protein_coding;locus_tag=RCAP_RS00135;old_locus_tag=RCAP_rcc00027 86 1
Chromosome sRNADetect sRNA 35149 35234 12.801680672268908 - . gene_id "755"; unique_id "RCAP_sre00755" Chromosome Protein Homology CDS 34143 35957 . - 0 ID=cds26;Parent=gene26;Dbxref=Genbank:WP_013065774.1;Name=WP_013065774.1;gbkey=CDS;product=membrane protein;protein_id=WP_013065774.1;transl_table=11 86 1
Chromosome sRNADetect sRNA 36497 36584 10.02395504610078 - . gene_id "756"; unique_id "RCAP_sre00756" Chromosome RefSeq gene 36047 36577 . - . ID=gene27;Name=RCAP_RS00140;gbkey=Gene;gene_biotype=protein_coding;locus_tag=RCAP_RS00140;old_locus_tag=RCAP_rcc00028 81 0.920455
```
* ``gene_biotype=*.gff`` The output from Categorizing process. These files are results of
categorizing the file overlap.gff. Every line from overlap.gff is picked and then put
into different GFF files according to its genetic biotypes. A single GFF file represents
a single category of genetic biotypes and contains transcripts that belong to this
category. Take the file gene biotype=misc_RNA.gff as an example, every line
from this file has ”misc_RNA” as the value of keyword gene biotype. 
```
Chromosome sRNADetect sRNA 3017608 3017674 318.10821386912625 - . gene_id "1130"; unique_id "RCAP_rcs01130" Chromosome RefSeq gene 3017368 3017698 . - . ID=gene2795;Name=RCAP_RS18220;gbkey=Gene;gene_biotype=misc_RNA;locus_tag=RCAP_RS18220 67 1
```
* ``annotatedTranscripts.gff`` The output from Sorting & Extracting process. It is the sorted version of the given annotated Genome file.

## Citing

If you use sRNA-Detect, please cite the following article:

[L. Peña-Castillo, M. Grüll, M.E. Mulligan and A.S. Lang, Detection of bacterial small transcripts from RNA-Seq data: a comparative assessment. Biocomputing 2016: pp. 456-467.](https://doi.org/10.1142/9789814749411_0042)


