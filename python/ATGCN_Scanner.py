#此脚本用于对收集的序列做质量控制，其主要功能是检查序列的碱基组成，包括ATGC以及其他兼并碱基，特别是N碱基在序列中的数量；
#根据GISAID对于测序结果high coverage(N碱基含量<1%)以及low coverage(N碱基含量>5%)的定义，以评价序列质量；
#一般情况下，优先选择N碱基含量小于5%的序列入主数据库，选择low coverage序列入副数据库，避免引物设计区域存在非ATGC碱基影响统计结果。

import sys
from Bio import SeqIO

if len(sys.argv) != 3:  ##参数数量不等于3的情况下，提示使用信息
    print("脚本用法:python ATGCN_Scanner.py sequences.fasta output_file")
    print("此脚本主要用于序列中N碱基数量和非ATGCN碱基统计，用于排查序列数据质量。")
    print("作者：张迪骏  日期：2020.6.4")
    exit()

if len(sys.argv) == 3:
    sequences_input = sys.argv[1]
    base_count_result_file = sys.argv[2]

    base_count = open(base_count_result_file, "w")
    base_count.write("Sequence ID" + "\t" + "Sequence length" + "\t" +
                     "N_Percentage" + "\t" + "N_judgement" + "\t" +
                     "Number_of_base_not_ATGCN" + "\n")

    for sequence_record in SeqIO.parse(sequences_input, "fasta"):

        A_count = str(sequence_record.seq).count("A")
        T_count = str(sequence_record.seq).count("T")
        G_count = str(sequence_record.seq).count("G")
        C_count = str(sequence_record.seq).count("C")
        N_count = str(sequence_record.seq).count("N")
        N_Percentage = N_count / len(str(sequence_record.seq))

        number_of_base_not_ATGCN = len(str(sequence_record.seq)) - A_count - T_count - G_count - C_count - N_count

        if N_Percentage < 0.01: #高覆盖度序列
            base_count.write(str(sequence_record.id) + "\t" + str(len(str(sequence_record.seq))) + "\t" +
                             str(round(N_Percentage, 3)) + "\t" + "H-C" + "\t" + str(number_of_base_not_ATGCN) + "\n")
        elif N_Percentage < 0.05: #中覆盖度序列
            base_count.write(str(sequence_record.id) + "\t" + str(len(str(sequence_record.seq))) + "\t" +
                             str(round(N_Percentage, 3)) + "\t" + "M-C" + "\t" + str(number_of_base_not_ATGCN) + "\n")
        else: #低覆盖度序列
            base_count.write(str(sequence_record.id) + "\t" + str(len(str(sequence_record.seq))) + "\t" +
                             str(round(N_Percentage, 3)) + "\t" + "L-C" + "\t" + str(number_of_base_not_ATGCN) + "\n")

    base_count.close()
