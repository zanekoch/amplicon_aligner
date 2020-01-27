package main

import(
  "fmt"
	"os"
	"bufio"
	"flag"
	"compress/gzip"
	"strings"

  "github.com/biogo/biogo/alphabet"
  "github.com/biogo/biogo/io/seqio/fastq"
	"github.com/biogo/biogo/seq/linear"

	//"github.com/biogo/biogo/index/kmerindex"
)

func readFastq(read_file1 string, read_file2 string) (reads1 *fastq.Reader, reads2 *fastq.Reader) {
	/*
	@param read_file1: path to first read file
	@param read_file2: path to second read file
	@return reads1, reads2: *linear.NewQseq readers of read files
	*/
	//unzip and read in read_file 1
	if fz, err := os.Open(read_file1); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v.", err)
		os.Exit(1)
	} else {
		if f, err := gzip.NewReader(fz); err != nil{
			fmt.Fprintf(os.Stderr, "Error: %v.", err)
			os.Exit(1)
		} else {
			//returns a fastq reader. Sequences based on linear.NewQSeq template
			reads1 = fastq.NewReader(f, linear.NewQSeq("", nil, alphabet.DNA, alphabet.Sanger))
		}
	}
	//unzip and read file 2
	if fz, err := os.Open(read_file2); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v.", err)
		os.Exit(1)
	} else {
		if f, err := gzip.NewReader(fz); err != nil{
			fmt.Fprintf(os.Stderr, "Error: %v.", err)
			os.Exit(1)
		} else {
			reads2 = fastq.NewReader(f, linear.NewQSeq("", nil, alphabet.DNA, alphabet.Sanger))
		}
	}
	return reads1, reads2
}

func processFasta(ref_fasta string, read_len int) (combined_ref string, reference_locs map[string][]int, x_locs []int) {
	/*
	@param fasta: name of file in fasta format of reference sequences to align to
	@param read_len: length of reads
	@return combined_ref: a string of concatanated reference sequences with read_length x's between each reference sequence
	@return reference_locs: a map with key = ref name and value = [position the reference starts in combined_ref, pos ref ends in combined_ref]
	@return x_locs: slice of form [startxRegion, endxRegion, startXRegion2, ...]
	*/
	reference_locs = make(map[string][]int)
	//read in fasta file
	var scanner *bufio.Scanner
	var all_lines []string
	if f, err := os.Open(ref_fasta); err != nil {
		fmt.Fprintf(os.Stderr, "Error: %v.", err)
		os.Exit(1)
	} else {
		scanner = bufio.NewScanner(f)
		defer f.Close()
	}
	var x string
	for i := 0; i < read_len; i ++ {
		x += "x"
	}
	//create combined reference string with read_len x's between each ref
	for scanner.Scan() {
		line := scanner.Text()
		all_lines = append(all_lines, line)
		if line[0] != '>'{
			combined_ref += line
			combined_ref += x
		}
	}
	//make reference_locs and x_locs
	line_num := len(all_lines)
	pos := 0
	for i := 0; i < line_num; i += 2 {
		ref_name := all_lines[i]
		ref_seq := all_lines[i+1]
		ref_len := len(ref_seq)
		temp := []int {pos, pos+ref_len}
		reference_locs[ref_name] = temp
		x_locs = append(x_locs, pos + ref_len)
		x_locs = append(x_locs, pos + ref_len + read_len)
		//increment positon by amount of new combined_ref created
		//which is length of reference added plus length of x's added
		pos += read_len + ref_len
	}
	return
}

func hashReferences(ref_fasta string, read_len int, hash_length int) (hash_table map[string][]int, reference_locs map[string][]int, x_locs []int, combined_ref string){
	/*
	@param ref_file: a path to reference file
	@param read_length: read length
	@param hash_length: length of keys for hash table
	@return hash_table: a dictionary with keys = key_size length seqs of ref and values = positions in combined ref where this seq exist
	@return reference_locs: a dictionary storing where each ref is in combined ref. key=ref name value=[start, end]
	@return x_locs: slice of form [startxRegion, endxRegion, startXRegion2, ...]
	@return combined_ref: a string of concatanated reference sequences with read_length x's between each reference sequence
	*/
	combined_ref, reference_locs, x_locs = processFasta(ref_fasta, read_len)
	//create hash table, repeated elements have their locations put in a list
  hash_table = make(map[string][]int)
  for i := 0; i < (len(combined_ref)-(hash_length - 1)); i++{
		//ignore any portions with x's because these will never be searched for
		if strings.ContainsAny(combined_ref[i:i + hash_length], "x") {
			continue
		}
		if locs, ok := hash_table[combined_ref[i:i + hash_length]]; ok {
			locs = append(locs, i)
			hash_table[combined_ref[i:i + hash_length]] = locs
		} else {
			locs := []int {i}
			hash_table[combined_ref[i:i + hash_length]] = locs
		}
	}
  return
}

func alignReads(reads1 *fastq.Reader, reads2 *fastq.Reader, hash_table map[string][]int, combined_ref string, x_locs []int, hash_length int, read_len int) {
	/*
	@param reads1, reads2: *linear.NewQseq readers of read files
	@param read_length: read length
	@param hash_length: length of keys for hash table
	@param hash_table: a dictionary with keys = key_size length seqs of ref and values = positions in combined ref where this seq exist
	@param x_locs: slice of form [startxRegion, endxRegion, startXRegion2, ...]
	@param combined_ref: a string of concatanated reference sequences with read_length x's between each reference sequence
	*/
	avg := 0
	n := 0
	for j:=0; j < 100; j++ {
		fmt.Println(j)
		//read in a single *linear.Qseq
		var read1 string
		r1, err := reads1.Read(); if err != nil {
			fmt.Fprintf(os.Stderr, "Error1: %v.", err)
			os.Exit(1)
		} else {
			 read1 = r1.(*linear.QSeq).String()
		}
		var read2 string
		r2, err := reads2.Read(); if err != nil {
			fmt.Fprintf(os.Stderr, "Error2: %v.", err)
			os.Exit(1)
		} else {
		 	read2 = r2.(*linear.QSeq).String()
		}
		var locs1, locs2 []int
		for i := 0; i < (read_len - hash_length); i++ {
			sub_read1 := read1[i:i + hash_length]
			sub_read2 := read2[i:i + hash_length]
			if loc, ok := hash_table[sub_read1]; ok {
				locs1 = loc
			}
			if loc, ok := hash_table[sub_read2]; ok {
				locs2 = loc
			}
			avg = avg + len(locs1) + len(locs2)
			n += 2
		}
	}
	fmt.Println(avg/n)
}

func main() {
	//all flags return a pointer of their type e.g. readFile1 = *string
	readFile1 := flag.String("r1", "", "Filename for read1. Required")
	readFile2 := flag.String("r2", "", "Filename for read2. Default none")
	ref_fasta := flag.String("f", "","fasta describing reference sequences. Required")
	out_dir := flag.String("o", "./", "output directory. Default current directory")
	hash_length := flag.Int("hash-length", 5, "length for hash indexes. Default 5")
 	read_len := flag.Int("read-length", 301, "length for hash indexes. Default 301")
	help := flag.Bool("help", false, "-r1 ReadFile1 -r2 ReadFile2 [optional] -f referenceFasta -o outputDir --hash-length len for hash index [default 5]")
	flag.Parse()
	if *help {
		flag.Usage()
		os.Exit(0)
	}
	fmt.Println(*ref_fasta)
	fmt.Println(*out_dir)
	fmt.Println(*hash_length)

	reads1, reads2 := readFastq(*readFile1, *readFile2)
	fmt.Println(*read_len)
	fmt.Println(*reads1)
	fmt.Println(*reads2)
	hash_table, reference_locs, x_locs, combined_ref := hashReferences(*ref_fasta, *read_len, *hash_length)
	_ = hash_table
	_ = reference_locs
	_ = x_locs
	_ = combined_ref
	alignReads(reads1, reads2, hash_table, combined_ref, x_locs, *hash_length, *read_len)
}
