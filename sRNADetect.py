import sys, argparse

import HTSeq


#    Written by Lourdes Pena-Castillo (lourdes at mun dot ca),
#    Memorial University of Newfoundland. Copyright 2015.
#
#    This program is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#        
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#        
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>

#   This program identifies small transcripts from prokaryotic RNA-Seq data
#   Reads a list of sam files (multiple sequencing runs for same organism)
#   Output a GTF file of identified small RNAs


#Read command line arguments
def readCommandLine():
    #create parser
    parser = argparse.ArgumentParser(usage = '%(prog)s [options] -samFile \"sam_file1\" -samFile \"sam_file2\"',
                                     description=
                                     'This script takes alignment files in SAM format and detect small transcripts with uniform coverage and no gaps',
                                     epilog = 'Written by Lourdes Pena-Castillo (lourdes at mun dot ca), Memorial University of Newfoundland. Copyright 2015.\n"Release under the terms of the GNU General Public License".')
    
    #add options
    #Parameters default values
    #minLength = 20
    #maxLength = 210
    #minHeight = 10
    #maxPctgDropValue = .5
    #maxPctgChangeNeg = -.01
    #maxPctgChangePos = 0.03
    
    parser.add_argument('-samFile', action='append', required=True, help = 'REQUIRED: sam file to read.')
    parser.add_argument('-minLen', action='store', default=20, type = int, help = 'minimum length to consider for small transcripts')
    parser.add_argument('-maxLen', action='store', default=210, type = int, help = 'maximum length to consider for small transcripts')
    parser.add_argument('-minH', action='store', default=10, type = int, help = 'minimum number of reads across all samples required to detect small transcripts')
    parser.add_argument('-maxPDV', action='store', default=0.5, type = float, help = 'maximum allowed drop in coverage wrt current coverage')
    parser.add_argument('-maxPN', action='store', default=-0.01, type = float, help = 'maximum allowed negative change (percentage) in mean coverage')
    parser.add_argument('-maxPP', action='store', default=0.03, type = float, help = 'maximum allowed positive change (percentage) in mean coverage')
    parser.add_argument('-idPrefix', action='store', default="sRNA",  help = 'Prefix to be used as identifier for sRNAs in gtf file')
    parser.add_argument('-out', action='store', default="sRNAs_1.gtf",  help = 'Filename of the output gtf file')
    parser.add_argument('-v', '--verbose', action='store_true', default=0)
    parser.add_argument('--version', action='version', version='%(prog)s 2.0')
    
    
    if len( sys.argv ) == 1:
        parser.print_help()
        sys.exit(1)
    
    args = parser.parse_args()
    
    if args.verbose:
        print args
    
    return args



#### MAIN PROGRAM
### Read parameters from command line

args = readCommandLine()
#Parameters
samFiles = args.samFile
minLength = args.minLen
maxLength = args.maxLen
minHeight = args.minH
maxPctgDropValue = args.maxPDV
maxPctgChangeNeg = args.maxPN
maxPctgChangePos = args.maxPP

#Create a coverage vector
cvg = HTSeq.GenomicArray( "auto", stranded=True, typecode="i" )

#For each sam file, add reads to corresponding interval 
print "Reading:", samFiles, "\n"
for file in samFiles:
    alignment_file = HTSeq.SAM_Reader(file)
    #iterate through all the reads and increment the value at the interval to which each read was aligned to
    for alngt in alignment_file:
        if alngt.aligned:
            cvg[ alngt.iv ] += 1


# Cycle through the intervals in the genomic array to identify small RNAs.
# Small RNAs are between minLength and maxLength nt
# Minimum height (ie., number of reads) to be considered is minHeight
# Allow a drop of maxPctgDropValue % from mean height to the current segment value
# Allow an increase of mean height of up to maxPctgChangePos %
# Allow an decrease of mean height of up to maxPctgChangeNeg %
# Calculate mean heigh using Knuth method
# M1 = X1, Mk = Mk-1 + (Xk - Mk-1)/k ---average
# k is number of bps included + the number of bps in the current genomic interval


#Variables
start = 0
end = 0
meanHeight = 0
length = 0
strand = ''
inSegment = False
chrom = ''


print "Identifying sRNAs \n"
#Keep list of sRNAs
sRNAs = [] #a list of lists, each internal list contains a GenomicInterval and a value
#Through all intervals in coverage array
#coverage array is continuos; that is cover the whole genome - an interval starts where the previous ends
for iv, value in cvg.steps():
    if args.verbose:
        print "Working with sequence:", iv.chrom
    #Changing chromosome or strand, if in a segment check to add it
    if (chrom != iv.chrom or strand != iv.strand) and inSegment:
        if meanHeight > value and minLength < length < maxLength: 
            #add as a sRNA
            toAdd = HTSeq.GenomicInterval( chrom, start, end, strand )
            sRNAs.append([toAdd, meanHeight])
            if args.verbose:
                print("Adding sRNA:", toAdd, meanHeight, pctgChange, dropValue, length)
        inSegment = False #set flag
    # keep chromosome
    chrom = iv.chrom
    if not inSegment:#Not considering a segment to be a sRNA
        if value < minHeight or iv.length > maxLength:
            #height below minimum or too long, skip to next genomicInterval
            if args.verbose:
                print("below minHeight or too long, next")
            continue
        else: 
            #right height and length, start a segment
            inSegment = True
            start = iv.start
            meanHeight = value
            end = iv.end
            length = iv.length
            strand = iv.strand
            if args.verbose:
                print("starting a segment\n")
    else: #Considering a segment to be a sRNA
            #check whether there is a change in number of reads above
            newMeanHeight = meanHeight + (value - meanHeight)/float(length + iv.length)
            pctgChange = (newMeanHeight - meanHeight)/float(meanHeight)
            dropValue = (meanHeight - value)/float(meanHeight)
            if pctgChange > maxPctgChangePos or pctgChange < maxPctgChangeNeg or dropValue > maxPctgDropValue:
                #if drop in value, add to sRNAs if length within range
                if meanHeight > value and minLength < length < maxLength: 
                    #add as a sRNA
                    toAdd = HTSeq.GenomicInterval( iv.chrom, start, end, strand )
                    sRNAs.append([toAdd, meanHeight])
                    inSegment = False
                    if args.verbose:
                        print("Adding sRNA:", toAdd, meanHeight, pctgChange, dropValue, length)
                elif meanHeight < value and  iv.length < maxLength: 
                    #if raise in value, restart if current length within allowed range
                    start = iv.start
                    meanHeight = value
                    end = iv.end
                    length = iv.length
                    strand = iv.strand
                    if args.verbose:
                        print("restarting ", meanHeight, pctgChange, dropValue, length)
                else:#otherwise discard
                    inSegment = False
                    if args.verbose:
                        print("discarding", meanHeight, pctgChange, dropValue, length)
            elif iv.strand == strand:
                #between allowed change, minHeigth and same strand, extend
                #Update length, meanHeight and end
                end = iv.end
                length += iv.length
                meanHeight = newMeanHeight
                if args.verbose:
                    print("extending ", meanHeight, pctgChange, dropValue)
            else: #different strand
                if args.verbose:
                    print("Error: Different strand!\n")
                exit()

## Where we in the middle of a segment when the coverage vector ended?
if inSegment and meanHeight > value and minLength < length < maxLength: 
    #add as a sRNA
    toAdd = HTSeq.GenomicInterval( chrom, start, end, strand )
    sRNAs.append([toAdd, meanHeight])
    inSegment = False
    if args.verbose:
        print("Adding sRNA:", toAdd, meanHeight, pctgChange, dropValue, length)

print "Writing gtf file.\n"
#print all identified sRNAs in a gtf file
counter = 1
myfile = open(args.out, "w")
for s in sRNAs:
    fields = s[0].chrom + "\tsRNADetect\tsRNA\t" + str(s[0].start+1) + "\t" + str(s[0].end +1) + \
    "\t" + str(s[1]) + "\t" + s[0].strand + "\t.\t" + "gene_id \"" + str(counter) + "\"; unique_id \"" + args.idPrefix + '%05d' % counter + "\"\n"
    myfile.write(fields)
    counter+=1
myfile.close()





