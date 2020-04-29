#!/usr/bin/python3

#This program is for convert the fasta file to primer3 input file.
#Usage: python seqs2primer3input.py [--help] [--vesion] [input_file] [output_file] [template_new_ids]
#Example: python ampliconExtracter.py A_input B_output A.
#A_input file must in fasta format. B_output file is the primer3 input file, A is the new id of templates, in output file it will be A_1..A_2.....A_n
#Version: 1.0
#Author: Zhangdijun  Date: 2020.4.29

from Bio import SeqIO
import sys

if len(sys.argv) != 4:
    print("\033[0;31;40mParameter error!!!!\033[0m")
    print("Usage: python ampliconExtracter.py A_input B_output A")
    print("Author: Zhangdijun  Date: 2020.4.29")
    print("Example: python ampliconExtracter.py A_input B_output A")
    print("A_input file must in fasta format. B_output file is the primer3 input file, A_ is the new id of templates, in output file it will be A_1..A_2.....A_n")
    sys.exit()

else:
    input_file = sys.argv[1]
    output_file = sys.argv[2]
    new_id = sys.argv[3]
    extract = open(output_file, 'w')

    i = 1
    for seq_record in SeqIO.parse(input_file,"fasta"):
        seq_record.id = "SEQUENCE_ID=" + new_id + "_" + str(i)  ##modify the name of the seqs
        extract.write(seq_record.id + "\n")  ##write the primer3 information to output_file
        extract.write("SEQUENCE_TEMPLATE=" + str(seq_record.seq) + "\n")
        extract.write("PRIMER_TASK=pick_pcr_primers" + "\n" +
                      "PRIMER_NUM_RETURN=20000" + "\n" +
                      "PRIMER_PRODUCT_SIZE_RANGE=100-400" + "\n" +
                      "PRIMER_MIN_SIZE=18" + "\n" +
                      "PRIMER_OPT_SIZE=22" + "\n" +
                      "PRIMER_MAX_SIZE=26" + "\n" +
                      "PRIMER_WT_SIZE_LT=0.0" + "\n" +
                      "PRIMER_WT_SIZE_GT=0.0" + "\n" +
                      "PRIMER_MIN_GC=30.0" + "\n" +
                      "PRIMER_MAX_GC=70.0" + "\n" +
                      "PRIMER_OPT_GC=50.0" + "\n" +
                      "PRIMER_MAX_END_GC=3" + "\n" +
                      "PRIMER_WT_GC_PERCENT_LT=0.0" + "\n" +
                      "PRIMER_WT_GC_PERCENT_GT=0.0" + "\n" +
                      "PRIMER_MIN_TM=59.0" + "\n" +
                      "PRIMER_OPT_TM=60.0" + "\n" +
                      "PRIMER_MAX_TM=61.0" + "\n" +
                      "PRIMER_PAIR_MAX_DIFF_TM=2" + "\n" +
                      "PRIMER_WT_TM_LT=0.0" + "\n" +
                      "PRIMER_WT_TM_GT=0.0" + "\n" +
                      "PRIMER_SALT_MONOVALENT=50.0" + "\n" +
                      "PRIMER_SALT_DIVALENT=3.0" + "\n" +
                      "PRIMER_DNA_CONC=200.0" + "\n" +
                      "PRIMER_THERMODYNAMIC_OLIGO_ALIGNMENT=1" + "\n" +
                      "PRIMER_MAX_SELF_ANY_TH=37.00" + "\n" +
                      "PRIMER_MAX_SELF_END_TH=37.00" + "\n" +
                      "PRIMER_MAX_HAIRPIN_TH=37.0" + "\n" +
                      "PRIMER_PAIR_MAX_COMPL_ANY_TH=37.00" + "\n" +
                      "PRIMER_PAIR_MAX_COMPL_END_TH=37.0" + "\n" +
                      "PRIMER_MAX_POLY_X=3" + "\n" +
                      "PRIMER_MIN_LEFT_THREE_PRIME_DISTANCE=0" + "\n" +
                      "PRIMER_MIN_RIGHT_THREE_PRIME_DISTANCE=0" + "\n" +
                      "=" + "\n")
        i += 1
