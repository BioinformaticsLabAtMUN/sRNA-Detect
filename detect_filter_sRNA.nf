#!/usr/bin/env nextflow

params.help = false
params.outputDir = '.'
params.alignmentDir = '.'
params.annotatedGenomeFile = ''
params.output = 'sRNAs_1.gtf'
params.idPrefix = 'sRNA'
params.minLength = 20
params.maxLength = 210
params.minHeight = 10
params.maxPctgDropValue = 0.5
params.maxPctgChangeNeg = -0.01
params.maxPctgChangePos = 0.03

// some of the following code was taken from: https://github.com/BioinformaticsLabAtMUN/sRNACharP/blob/master/sRNACharP.nf
// show how to use this file.
if(params.help) {
    log.info 'sRNADetect: sRNA detect and filter pipeline'
    log.info '----------------------------------------------------'
    log.info 'Options:'
    log.info '--alignmentDir  directory   the directory of SAM/BAM files [required]'
    log.info '--outputDir  directory   the directory for saving output files'
    log.info '--annotatedGenomeFile  path   the path to the given Genome file'
    log.info '--output  fileName    the name of the sRNA-Detect output gtf file'
    log.info '--idPrefix    prefix  Prefix to be used as identifier for sRNAs in gtf file'
    log.info '--minLength    integer  minimum length to consider for small transcripts'
    log.info '--maxLength    integer  maximum length to consider for small transcripts'
    log.info '--minHeight    integer  minimum number of reads across all samples required to detect small transcripts'
    log.info '--maxPctgDropValue    percentages  maximum allowed drop in coverage wrt current coverage'
    log.info '--maxPctgChangeNeg    maximum allowed negative change (percentage) in mean coverage'
    log.info '--maxPctgChangePos    maximum allowed positive change (percentage) in mean coverage'
    exit 0
}

alignmentFileList = Channel.fromPath( "${params.alignmentDir}/*.{sam,bam}").collect()
indexFileList = Channel.fromPath( "${params.alignmentDir}/*.bai").collect()

process sRNA_Detect {
    publishDir "${params.outputDir}/sRNADetected"

    input:
    file alignmentFileList
    file indexFileList
    output:
    file params.output into sRNADetected

    script:
    """
    #!/usr/bin/env python3
    import HTSeq
    from multiprocessing import Pool

    # some of the following code cite the reference: https://github.com/BioinformaticsLabAtMUN/sRNA-Detect/blob/master/sRNADetect.py
    def detect_sRNA(cvg):
        # Variables
        start = 0
        end = 0
        meanHeight = 0
        length = 0
        strand = ''
        inSegment = False
        chrom = ''

        # Keep list of sRNAs
        sRNAs = []  # a list of lists, each internal list contains a GenomicInterval and a value
        # Through all intervals in coverage array
        # coverage array is continuos; that is cover the whole genome - an interval starts where the previous ends
        for iv, value in cvg.steps():
            # Changing chromosome or strand, if in a segment check to add it
            if (chrom != iv.chrom or strand != iv.strand) and inSegment:
                if meanHeight > value and ${params.minLength} < length < ${params.maxLength}:
                    # add as a sRNA
                    toAdd = HTSeq.GenomicInterval(chrom, start, end, strand)
                    sRNAs.append([toAdd, meanHeight])
                inSegment = False  # set flag
            # keep chromosome
            chrom = iv.chrom
            if not inSegment:  # Not considering a segment to be a sRNA
                if value < ${params.minHeight} or iv.length > ${params.maxLength}:
                    # height below minimum or too long, skip to next genomicInterval
                    continue
                else:
                    # right height and length, start a segment
                    inSegment = True
                    start = iv.start
                    meanHeight = value
                    end = iv.end
                    length = iv.length
                    strand = iv.strand
            else:  # Considering a segment to be a sRNA
                # check whether there is a change in number of reads above
                newMeanHeight = meanHeight + (value - meanHeight) / float(length + iv.length)
                pctgChange = (newMeanHeight - meanHeight) / float(meanHeight)
                dropValue = (meanHeight - value) / float(meanHeight)
                if pctgChange > ${params.maxPctgChangePos} or pctgChange < ${params.maxPctgChangeNeg} or dropValue > ${params.maxPctgDropValue}:
                    # if drop in value, add to sRNAs if length within range
                    if meanHeight > value and ${params.minLength} < length < ${params.maxLength}:
                        # add as a sRNA
                        toAdd = HTSeq.GenomicInterval(iv.chrom, start, end, strand)
                        sRNAs.append([toAdd, meanHeight])
                        inSegment = False
                    elif meanHeight < value and iv.length < ${params.maxLength}:
                        # if raise in value, restart if current length within allowed range
                        start = iv.start
                        meanHeight = value
                        end = iv.end
                        length = iv.length
                        strand = iv.strand
                    else:  # otherwise discard
                        inSegment = False
                elif iv.strand == strand:
                    # between allowed change, minHeigth and same strand, extend
                    # Update length, meanHeight and end
                    end = iv.end
                    length += iv.length
                    meanHeight = newMeanHeight
                else:  # different strand
                    exit()

        ## Where we in the middle of a segment when the coverage vector ended?
        if inSegment and meanHeight > value and ${params.minLength} < length < ${params.maxLength}:
            # add as a sRNA
            toAdd = HTSeq.GenomicInterval(chrom, start, end, strand)
            sRNAs.append([toAdd, meanHeight])
            inSegment = False

        # print all identified sRNAs in a gtf file

        counter = 1
        myfile = open("${params.output}", "w")
        for s in sRNAs:
            fields = s[0].chrom + "\\tsRNADetect\\tsRNA\\t" + str(s[0].start + 1) + "\\t" + str(s[0].end + 1) +\\
                     "\\t" + str(s[1]) + "\\t" + s[0].strand + "\\t.\\t" + "gene_id \\"" + str(
                counter) + "\\"; unique_id \\"" + "${params.idPrefix}" + '%05d' % counter + "\\"\\n"
            myfile.write(fields)
            counter += 1
        myfile.close()


    def read_Alignment(filename):
        alignment_file = HTSeq.BAM_Reader(filename)
        ivList = []
        for alngt in alignment_file:
            if alngt.aligned:
                ivList.append(alngt.iv)
        return ivList


    if __name__ == '__main__':
        # create a coverage vector
        cvg = HTSeq.GenomicArray('auto', stranded=True, typecode='i')
        alignmentFiles = "$alignmentFileList".split()
        with Pool(processes=len(alignmentFiles)) as pool:
            multiple_results = [pool.apply_async(read_Alignment, (name,)) for name in alignmentFiles]
            for res in multiple_results:
                for iv in res.get():
                    cvg[iv] += 1
        detect_sRNA(cvg)
    """
}
rawGenomeAnnotation = file("${params.annotatedGenomeFile}")
process sort_sRNA {

    input:
    file "sRNA.gtf" from sRNADetected
    file "rawGenomeAnnotation.gff" from rawGenomeAnnotation
    output:
    file "sorted_sRNA.gtf" into sort_sRNA
    file "annotatedTranscripts.gff" into annotatedTranscript
    script:
    """
    bedtools sort -i sRNA.gtf > sorted_sRNA.gtf
    bedtools sort -i rawGenomeAnnotation.gff | awk -F"\\t" '\$3=="gene" {print}{next}' > annotatedTranscripts.gff
    """
}
process filter_sorted_sRNA {
    publishDir "${params.outputDir}/filtered_sRNA"

    input:
    file "sorted_sRNA.gtf" from sort_sRNA
    file "annotatedTranscripts.gff" from annotatedTranscript
    output:
    file "overlap.gff" into overlapped_sRNA
    file "non_overlap.gtf"

    script:
    """
    bedtools intersect -wo -a sorted_sRNA.gtf -b annotatedTranscripts.gff -sorted -v -s -f 0.8 > non_overlap.gtf
    bedtools intersect -wo -a sorted_sRNA.gtf -b annotatedTranscripts.gff -sorted -s -f 0.8 | awk '{\$(NF+1) = \$NF / ((\$5 - \$4) + 1)} 1' > overlap.gff
    """
}

process separate_overlapped_sRNA {
    publishDir "${params.outputDir}/separated_sRNA"

    input:
    file "overlap.gff" from overlapped_sRNA
    file "annotatedTranscripts.gff" from annotatedTranscript
    output:
    file "*.gff"
    file "annotatedTranscripts.gff"

    script:
    """
    (cut -f9 annotatedTranscripts.gff | cut -d";" -f4 | grep biotype | sort -u) | while IFS="" read -r line; do if grep -q "\$line" overlap.gff; then grep "\$line" overlap.gff > "\$line.gff"; fi done
    """
}
