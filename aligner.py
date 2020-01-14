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
    reference_locs = dict()
    f = open(fasta, 'r')
    #create read_length length stirng of x's
    x = ""
    for i in range(read_length):
        x += 'x'
    #readlines() returns a list with ith line as ith element
    lines = f.readlines()
    combined_ref = ""
    pos = 0
    for i in range(0, len(lines)-1, 2):
        ref_name = lines[i][1:-1]
        ref_seq = lines[i+1][:-1]
        reference_locs[ref_name] = [pos, pos + len(ref_seq)]
        combined_ref += ref_seq
        combined_ref += x
        pos += read_length + len(ref_seq)
    return combined_ref, reference_locs

def hashReferences(ref_file, read_length, key_size=5):
    '''
    @param ref_file: a path to reference file
    @param read_length: read length
    @return hash_table: a dictionary with keys=key_size length seqs of ref and values=positions in combined ref where this seq exist
    @return reference_locs: a dictionary storing where each ref is in combined ref. key=ref name value=[start, end]
    '''
    #create combined reference
    combined_ref, reference_locs = processFasta(ref_file, read_length)
    #create hash table, repeated elements have their locations put in a list
    hash_table = defaultdict(list)
    for i in range(len(combined_ref)-(key_size - 1)):
        hash_table[combined_ref[i:i + key_size]].append(i)
    return hash_table, reference_locs

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

def alignReads(hash_table, reference_locs, hash_size, read_length, msa1, msa2 = None, ):
    #paired alignment
    if msa2 != None:
        #for each read from a sample, do in order of read1_1, read2_1, read1_2, read2_2...
        for read in zip_longest(msa1, msa2):
            #search every hash_size fragment of read in hash of references
            for i in range(read_length - hash_size):
                sub_read = read[i:i + hash_size]
                #check if the read is in the table
                if sub_read in hash_table:
                    #loop across the locations it was found in the reference
                    for loc in hash_table[sub_read]:
                        #where the complete read starts and ends
                        read_start = loc - i
                        read_end = read_start + read_length
                        #check if read overlaps x's
                        
                        for base in range(read_start, read_start + read_length):

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
    optional_group.add_argument("--hash-length", dest="hash_size", required=False, default=5, help="reference subsequence length")
    args = parser.parse_args()
    sample = args.sample
    read1_file = args.r1
    read2_file = args.r2
    ref_file = args.ref
    job_manager = args.job_manager
    out_dir = args.odir
    hash_size = args.hash_size
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
    hash_table, reference_locs = hashReferences(ref_file, read_length, hash_size)
    logging.debug("hashed")
    alignReads(hash_table, reference_locs, hash_size, read_length, msa1, msa2)
    logging.debug("finished")
    output_log = os.open(os.path.join(out_dir, "%s_log.out" % sample), os.O_WRONLY | os.O_CREAT)



if __name__ == '__main__':
    main()
