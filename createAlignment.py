import re
import os
import argparse
import logging

'''
dependencies
-python3

'''

DEBUG = 1
TESTRUN = 0
PROFILE = 0

class CLIError(Exception):
    '''Generic exception to raise and log different fatal errors.'''
    def __init__(self, msg):
        super(CLIError).__init__(type(self))
        self.msg = "E: %s" % msg
    def __str__(self):
        return self.msg
    def __unicode__(self):
        return self.msg

def expandPath(path):
    user_match = re.match('^(~)(.*)$', path)
    if user_match:
        path = os.path.expanduser(path)
    return os.path.abspath(path)

def findReads(path):
    '''
    return: a list of tuples named Read that contain the sample name followed by an array containing one element, a path to the merged read file, or 2 elements, paths to read1 and read2 files
    '''
    import subprocess
    from collections import namedtuple
    read_list = []
    Read = namedtuple('Read', ['sample', 'reads'])
    for file in os.listdir(path):
        is_read = re.search('(.*)(\.f(?:ast)?q(\.gz)?)$', file, re.IGNORECASE)
        if is_read:
            sample_name = is_read.group(1)
            full_file = os.path.join(path, file)
            #check if read file is empty, if so add to read list with empty contents
            if os.path.getsize(full_file) == 0 or (is_read.group(3) and subprocess.getoutput("gzip -l %s | awk 'NR > 1{print $2}'" % full_file) == '0'):
                logging.warning("Read file %s has no data, skipping..." % file)
                read_list.append(Read(sample_name, None))
                continue
            #check if the read1 read2 files have been merged
            is_merged = re.search('^(.*?)(?:[_\.](?:assembled|merged))+$', sample_name, re.IGNORECASE)
            if is_merged:
                sample_name = is_merged.group(1)
                read = Read(sample_name, [os.path.join(path, file)])
                read_list.append(read)
                logging.info(read)
            else:
                is_paired = re.search('^(?:((.*?)(?:_L\d\d\d)?(?:(?:[_\.](?:R(?:ead)?)?)))([12])([_\.])?)(?!.*[_\.](?:R(?:ead)?)?[12][_\.])(.*)$', sample_name, re.IGNORECASE)
                if is_paired:
                    if is_paired.group(3) == '1':  # If paired, only process read 1, so we don't double count the pair, see TODO below
                        sample_name = is_paired.group(2)
                        read1 = file
                        read2 = "%s2%s%s%s" % (is_paired.group(1), is_paired.group(4) or '', is_paired.group(5), is_read.group(2))
                        #print("\t%s\t%s\t%s" % (sample_name, read1, read2))
                        if os.path.exists(os.path.join(path, read2)):
                            read = Read(sample_name, [os.path.join(path, read1), os.path.join(path, read2)])
                            read_list.append(read)
                            logging.info(read)
                        else:
                            # TODO: If only R2 exists, it won't be included
                            logging.warning("Cannot find %s, the matching read to %s. Including as unpaired..." % (read2, read1))
                            read = Read(sample_name, [os.path.join(path, read1)])
                            read_list.append(read)
                            logging.info(read)
                else: #Read is unpaired
                    sample_name = is_read.group(1)
                    read = Read(sample_name, [os.path.join(path, file)])
                    read_list.append(read)
                    logging.info(read)
    return read_list

def main():
    try:
        parser = argparse.ArgumentParser()
        required_group = parser.add_argument_group("required arguments")
        required_group.add_argument("-n", "--name", required=True, help="name for this run. [REQUIRED]")
        required_group.add_argument("-j", "--json", metavar="FILE", required=True, type=argparse.FileType('r'), help="JSON file of assay descriptions. [REQUIRED]")
        required_group.add_argument("-r", "--read-dir", dest="rdir", metavar="rdir", required=True,  help="Directory of reads to be aligned. [REQUIRED]")
        optional_group = parser.add_argument_group("optional arguments")
        optional_group.add_argument("-o", "--out-dir", dest="odir", metavar="DIR", help="directory to write output files to. [default: `pwd`]")
        #parse arguments
        args = parser.parse_args()
        run_name = args.name
        json_filename = expandPath(args.json.name)
        read_dir = args.rdir
        out_dir = args.odir
        #make out_dir and read_dir current directory if not speicified
        if not out_dir:
            out_dir = expandPath(os.getcwd())
        if not read_dir:
            read_dir = expandPath(os.getcwd())
        #if the output path already exists check if should overwrite
        if os.path.exists(out_dir):
            response = input("\nOutput folder %s already exists!\nFiles in it may be overwritten!\nShould we continue anyway [N]? " % out_dir)
            if not re.match('^[Yy]', response):
                print("Operation cancelled!")
                quit()
        else:
            os.makedirs(out_dir)
        #set up logging in output directory
        logfile = os.path.join(out_dir, "TEMP.log")
        logging.basicConfig(level=logging.DEBUG,
                            format='%(asctime)s %(levelname)-8s %(message)s',
                            datefmt='%m/%d/%Y %H:%M:%S',
                            filename=logfile,
                            filemode='w')
        logging.info("Aligning reads in %s to refs in JSON file: %s for run: %s" % (read_dir, json_filename, run_name))
        #get reads into a list of format
        read_list = findReads(read_dir)
        print(read_list)

        return 0
    except KeyboardInterrupt:
        return 0
    except Exception as e:
        if DEBUG or TESTRUN:
            raise(e)
        indent = len(program_name) * " "
        sys.stderr.write(program_name + ": " + repr(e) + "\n")
        sys.stderr.write(indent + "  for help use --help")
        return 2


if __name__ == '__main__':
    main()
