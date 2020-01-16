import argparse
import os
import logging
from collections import defaultdict
from skbio import TabularMSA, DNA
from itertools import zip_longest


import dispatcher

def processFasta(fasta, read_length):
    '''
    @param fasta: name of file in fasta format of reference sequences to align to
    @return combined_ref: a string with read_length x's between each reference seq
    @return reference_locs: a dictionary with key = ref name and value = position the reference starts, pos ref ends
    '''
    f = open(fasta, 'r')
    #create read_length length stirng of x's
    x = ""
    for i in range(read_length):
        x += 'x'
    #readlines() returns a list with ith line as ith element
    lines = f.readlines()
    combined_ref = ""
    reference_locs = dict()
    x_locs = []
    pos = 0
    for i in range(0, len(lines)-1, 2):
        ref_name = lines[i][1:-1]
        ref_seq = lines[i+1][:-1]
        ref_len = len(ref_seq)
        reference_locs[ref_name] = [pos, pos + ref_len]
        x_locs.append(pos + ref_len)
        x_locs.append(pos + ref_len + read_length)
        combined_ref += ref_seq
        combined_ref += x
        #increment positon by amount of new combined_ref created which is length of reference added plus length of x's added
        pos += read_length + ref_len
    return combined_ref, reference_locs, x_locs

def hashReferences(ref_file, read_length, hash_length=5):
    '''
    @param ref_file: a path to reference file
    @param read_length: read length
    @return hash_table: a dictionary with keys=key_size length seqs of ref and values=positions in combined ref where this seq exist
    @return reference_locs: a dictionary storing where each ref is in combined ref. key=ref name value=[start, end]
    '''
    #create combined reference
    combined_ref, reference_locs, x_locs = processFasta(ref_file, read_length)
    #create hash table, repeated elements have their locations put in a list
    hash_table = defaultdict(list)
    for i in range(len(combined_ref)-(hash_length - 1)):
        hash_table[combined_ref[i:i + hash_length]].append(i)
    return hash_table, reference_locs, x_locs, combined_ref

def getReads(read1_file, read2_file = ''):
    '''
    @param read_list: a list of tuples named Read that contain the sample name followed by an array containing one element, a path to the merged read file, or 2 elements, paths to read1 and read2 files
    @return:
    '''
    #put reads into tabulr multiple sequence alignment structures
    msa1 = TabularMSA.read(read1_file, constructor=DNA, variant='illumina1.8')
    if read2_file != '':
        msa2 = TabularMSA.read(read2_file, constructor=DNA, variant='illumina1.8')
        return msa1, msa2
    return msa1

def in_range(range_list, value, start):
    '''
    @param range_list: a SORTED FROM SMALLER TO LARGER list of form [start, end, start(bigger than last start), end(bigger than last end)...] with start and end representing ranges
    @param value: value to check if is in a range of range_list
    @param start: bool=True when value is start of read and False when end of read
    @return: list with range of where value is in a range of range_list
    '''
    present_range = []
    for i in range(0, len(range_list), 2):

        if range_list[i] <= value:
            if value <= range_list[i + 1]:
                if start:
                    present_range.append(value)
                    present_range.append(range_list[i + 1])
                    #can return because will only be in one range
                    return present_range
                else:
                    present_range.append(range_list[i])
                    present_range.append(value)
                    return present_range
        else:
            #can immedietely return because range_list is increasing
            return present_range
    return present_range

def alignReads(hash_table, combined_ref, x_locs, hash_length, read_length, msa1, msa2 = None, ):
    #paired alignment
    if msa2 != None:
        #for each read from a sample, do in order of read1_1, read2_1, read1_2, read2_2...
        for read in iter(msa1):
            #search every hash_length fragment of read in hash of references
            str_read = str(read)
            for i in range(read_length - hash_length):
                sub_read = str_read[i:i + hash_length]
                #check if the read is in the table
                if sub_read in hash_table:
                    #loop across the locations it was found in the reference
                    for loc in hash_table[sub_read]:
                        #where the complete read starts and ends
                        read_start = loc - i
                        read_end = read_start + read_length
                        #check if read overlaps x's
                        overlaps_start = in_range(x_locs, read_start, True)
                        overlaps_end = in_range(x_locs, read_end, False)
                        #trim locations where read hangs off ref (i.e. overlaps x's)
                        if overlaps_start != []:
                            overlap = overlaps_start[1] - overlaps_start[0]
                            read = str_read[overlap: ]
                            new_start = read_start + overlap
                        else:
                            new_start = read_start
                        if overlaps_end != []:
                            overlap = overlaps_end[1] - overlaps_end[0]
                            read = str_read[: -overlap]
                            new_end = read_end - overlap
                        else:
                            new_end = read_end
                        #align the non-trimmed region
                        diff = 0
                        for base in range(new_start, new_end):
                            if combined_ref[base] != str_read[base-new_start]:
                                diff += 1
    #not paired alignment
    else:
        pass
    return

def main():
    parser = argparse.ArgumentParser()
    required_group = parser.add_argument_group("required arguments")
    required_group.add_argument("-s", "--sample", dest="sample", required=True, help="name of read file(s) sample")
    required_group.add_argument("-r1", "--read_file1", dest="r1", required=True, help="read file of read1s to be aligned")
    required_group.add_argument("-f", "--reference-fasta", dest="ref", required=True, help="fasta describing reference sequences")
    required_group.add_argument("-o", "--out-dir", dest="odir", metavar="DIR", help="directory to write output files to.")
    #required_group.add_argument("--function", dest="function", required=True, help="which function to run")
    optional_group = parser.add_argument_group("optional arguments")
    optional_group.add_argument("-r2", "--read_file2", dest="r2", required=False, default='', help="read file of read2s to be aligned [default: '']")
    optional_group.add_argument("-j", "--job_manager", dest="job_manager", default="SLURM", help="cluster job submitter to use (PBS, SLURM, SGE, none). [default: SLURM]")
    optional_group.add_argument("--hash-length", dest="hash_length", type=int, required=False, default=5, help="reference subsequence length")
    args = parser.parse_args()
    sample = args.sample
    read1_file = args.r1
    read2_file = args.r2
    ref_file = args.ref
    job_manager = args.job_manager
    out_dir = args.odir
    hash_length = args.hash_length
    #get reads from fastq files into skbio multiple sequence alignment objects
    if read2_file != '':
        msa1, msa2 = getReads(read1_file, read2_file)
    else:
        msa1 = getReads(read1_file)
    #get the read length for use in creating reference hash
    read_length = -1
    for seq in iter(msa1):
        read_length = len(seq)
    #create hash table of combined reference sequence and way to decode positon meaning
    hash_table, reference_locs, x_locs, combined_ref = hashReferences(ref_file, read_length, hash_length)
    logging.debug("hashed")
    alignReads(hash_table, combined_ref, x_locs, hash_length, read_length, msa1, msa2)
    logging.debug("finished")
    output_log = os.open(os.path.join(out_dir, "%s_log.out" % sample), os.O_WRONLY | os.O_CREAT)



if __name__ == '__main__':
    main()
