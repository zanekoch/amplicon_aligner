from collections import defaultdict
from skbio import TabularMSA, DNA
import argparse

import dispatcher

def processFasta(fasta):
    '''
    @param fasta: name of file in fasta format of reference sequences to align to
    @return: a dictionary with keys=ref names and values=ref sequences
    '''
    references = dict()
    f = open(fasta, 'r')
    lines = f.readlines()
    for i in range(0, len(lines)-1, 2):
        references[lines[i][1:-1]] = lines[i+1][:-1]
    return references

def hashReferences(ref_file, key_size=5):
    '''
    @param references: a dictionary with keys=ref names and values=ref sequences
    @return: a dictionary with keys=reference name and values=default dict hashing corresponding reference by key_size
    '''
    references = processFasta(ref_file)
    hashes = dict()
    for reference in references:
        hash_table = defaultdict(list)
        reference_seq = references[reference]
        #iterate across the reference sequence hashing 5bp chunks into hashtable, recording their 0-based location
        for i in range(len(reference_seq)-(key_size - 1)):
            hash_table[reference_seq[i:i + key_size]].append(i)
        hashes[reference] = hash_table
    return hashes

def extractReads(read_file):
    '''
    @param read_tuple: sample of reads to be read in
    @return: a list of skbio tabularMSAs of reads
    '''
    MSAs = []
    #iterate across read paths

    MSAs.append(msa)
    return MSAs

def align(sample, read1_file, read2_file, ref_file, out_dir):
    '''
    @param read_list: a list of tuples named Read that contain the sample name followed by an array containing one element, a path to the merged read file, or 2 elements, paths to read1 and read2 files
    @param references: a dictionary with keys=ref names and values=ref sequences
    @return:
    '''
    hashes = hashReferences(ref_file)
    msa1 = TabularMSA.read(read1_file, constructor=DNA, variant='illumina1.8')
    if read2_file != '':
        msa2 = TabularMSA.read(read2_file, constructor=DNA, variant='illumina1.8')
    file1 = out_dir + '/' + sample + "_1"
    with open(file1, 'w+') as f1:
        print(msa1, file=f1)
    file2 = out_dir + '/' + sample + "_2"
    with open(file2, 'w+') as f2:
        print(msa2, file=f2)


def main():
    parser = argparse.ArgumentParser()
    required_group = parser.add_argument_group("required arguments")
    required_group.add_argument("-s", "--sample", dest="sample", required=True, help="name of read file(s) sample")
    required_group.add_argument("-r1", "--read_file1", dest="r1", required=True, help="read file of read1s to be aligned")
    required_group.add_argument("-f", "--reference-fasta", dest="ref", required=True, help="fasta describing reference sequences")
    required_group.add_argument("-o", "--out-dir", dest="odir", metavar="DIR", help="directory to write output files to.")
    optional_group = parser.add_argument_group("optional arguments")
    optional_group.add_argument("-r2", "--read_file2", dest="r2", required=False, default='', help="read file of read2s to be aligned [default: '']")
    optional_group.add_argument("-j", "--job_manager", dest="job_manager", default="SLURM", help="cluster job submitter to use (PBS, SLURM, SGE, none). [default: SLURM]")
    args = parser.parse_args()
    sample = args.sample
    read1_file = args.r1
    read2_file = args.r2
    ref_file = args.ref
    job_manager = args.job_manager
    out_dir = args.odir

    align(sample, read1_file, read2_file, ref_file, out_dir)



if __name__ == '__main__':
    main()
