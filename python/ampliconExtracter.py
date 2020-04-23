#!/usr/bin/python3

#This program is for extracting amplicon sequences from given left/right primers and templates.
#The input file must have left primer, right primer and templates informations.
#Usage: python ampliconExtracter.py [--help] [--vesion] [input_file] [output_file] [Column_of_left_primer] [Column_of_right_primer] [Column_of_template]
#Example: python ampliconExtracter.py A_input B_output 3 5 7.
#3: column containing left primers; 5: column containing right primers; 7: column containing templates.
#Version: 1.01
#Author: Zhangdijun  Date: 2020.4.22

import sys  ##import system package, call positional arguments function
from Bio.Seq import Seq  ##import biopython package, call reverse complement function
from Bio.Alphabet import IUPAC

if len(sys.argv) != 6:  ##positional arguments not equal to 6, print usage information
    print("\033[0;31;40mParameter error!!!!\033[0m")
    print("Usage:python ampliconExtracter.py [--help] [--vesion] [input_file] [output_file] [Column_of_left_primer] [Column_of_right_primer] [Column_of_template]")
    print("Author: Zhangdijun  Date: 2020.4.22")
    print("Example: python ampliconExtracter.py A_input B_output 3 5 7")
    print("3: column containing left primers; 5: column containing right primers; 7: column containing templates")
    sys.exit()

elif sys.argv[1].startswith('--'):
    option = sys.argv[1][2:]
    if option == 'version':  ##positional arguments is --vesion, print version information and exit
        print("Version 1.01")
    elif option =='help':  ####positional arguments is --help, print help information and exit
        print("Usage:python ampliconExtracter.py [--help] [--vesion] [input_file] [output_file] [Column_of_left_primer] [Column_of_right_primer] [Column_of_template]")
        print("Author: Zhangdijun  Date: 2020.4.22")
        print("Example: python ampliconExtracter.py A_input B_output 3 5 7")
        print("3: column containing left primers; 5: column containing right primers; 7: column containing templates")
    sys.exit()

else:
    input_file = sys.argv[1]  ##positional arguments 1: input_file
    output_file = sys.argv[2]  ##positional arguments 2: output_file
    column_of_left_primer = sys.argv[3]  ##positional arguments 3:column containing left primers
    column_of_right_primer = sys.argv[4]  ##positional arguments 4:column containing right primers
    column_of_template = sys.argv[5]  ##positional arguments 5:column containing templates

    extract = open(output_file, 'w')
    extract.write("left_primers" + "\t" + "right_primers" + "\t" + "amplicon_seqs" + "\t" + "amplicon_length" + "\n")  ##write table headers of ouuput file

    total_line = 0
    with open(input_file, 'r') as list:
        for line in list:
            left_Primers = line.strip().split('\t')[int(column_of_left_primer)-1]  ##extract left primer sequence
            right_Primers = Seq(line.strip().split('\t')[int(column_of_right_primer)-1],IUPAC.unambiguous_dna)  ##extract right primers sequence and convert to DNA sequence type, which is necessary step for reverse complement
            right_Primers_RC = str(right_Primers.reverse_complement())  ##reverse complement sequence of right primers
            template = line.strip().split('\t')[int(column_of_template)-1]  ##extract template sequense
            f_location = template.find(left_Primers)  ##find the index of left primer in template
            r_location = template.find(right_Primers_RC)  ##find the index of right primer in template
            amplicon = template[f_location:r_location+len(right_Primers)]  ##extract amplicon sequence between the index of left/right primers
            extract.write(left_Primers + "\t" + str(right_Primers) + "\t" + amplicon + "\t" + str(len(amplicon)) + "\n")  ##write extract informations into output file
            total_line += 1
            if total_line%100 == 0:
                if total_line != 0:
                    print(str(total_line) + " " + "amplicons extracted")  ##print extraction progress
    print("=====Extract finished=====")
