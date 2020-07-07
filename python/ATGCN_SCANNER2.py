#coding=utf-8

#  1.此程序用于对收集的序列进行质量控制，其主要功能是检查序列的碱基组成,包括GC百分比/兼并碱基含量，特别是N碱基在序列中的占比以及N的分布类型
#  （分散型或是连续型）；
#  2.此程序参考GISAID数据库对测序覆盖度的评价体系，以1%和5%为分界线进行评价，N碱基含量小于1%则为高覆盖度序列（High coverage,H-C),N碱基
#   含量大于等于1%且小于5%则为中覆盖度序列（Moderate coverage, M-C)，N碱基含量大于等于5%则为地覆盖度序列（Low coverage,L-C);
#  3.在数据量足够大的前提下，可只选择H-C的序列入库，可根据实际情况酌情入组M-C的序列，出现数据量较少，且存在大量L-C的序列时，则需要通过判断
#   N碱基是否连续以进行数据选择，若序列中存在连续N碱基（连续N碱基数量约等于总N数量），则该序列的其余部分可用于数据库；若序列中的N碱基呈现随
#   机分布，则该序列的其余部分的质量存在一定缺陷，应谨慎使用该序列。

import sys
from Bio import SeqIO

if len(sys.argv) != 3:  ##参数数量不等于3的情况下，提示使用信息
    print("程序用法:python ATGCN_SCANNER2.py sequences.fasta result_file")
    print("程序主要用于序列序列排查，主要功能包括：")
    print("1.检查序列的碱基组成,包括GC百分比/非ATGCN碱基含量;")
    print("2.N碱基在序列中的占比以及N的分布类型(连续型或分散型）。")
    print("作者：张迪骏  日期：2020.7.6")
    exit()

sequences_input = sys.argv[1]  #输入文件
result_file = sys.argv[2]  #输出文件

result_output = open(result_file, "w")
result_output.write("Seq_ID" + "\t" + "Seq_length" + "\t" + "GC_percentage" + "\t" + "Not_ATGCN_number" + "\t" +
                    "N_number" + "\t" + "N_judgement" + "\t" + "Max_Ns_length" + "\t" + "max_Ns_length/N_number" + "\n")

for sequence_record in SeqIO.parse(sequences_input, "fasta"):
    A_count = sequence_record.seq.count("A")  #A碱基计数
    T_count = sequence_record.seq.count("T")  #T碱基计数
    G_count = sequence_record.seq.count("G")  #G碱基计数
    C_count = sequence_record.seq.count("C")  #C碱基计数
    N_count = sequence_record.seq.count("N")  #N碱基计数
    No_ATGCN_count = len(str(sequence_record.seq)) - A_count - T_count - G_count - C_count - N_count  #非ATGCN碱基计数
    GC_percentage = round((G_count + C_count) / len(str(sequence_record)), 2)  #GC百分比
    N_percentage = N_count / len(str(sequence_record.seq))  #N百分比

    N_location_list = []
    if N_count < 1:  #序列中不包含N碱基的情况
        result_output.write(str(sequence_record.id) + "\t" + str(len(str(sequence_record.seq))) + "\t" +
                            str(GC_percentage) + "\t" + str(No_ATGCN_count) + "\t" + str(N_count) + "\t" +
                            "H-C" + "\t" + "0" + "\t" + "-" + "\n")
    elif N_count == 1:  #序列中包含N碱基的情况
        result_output.write(str(sequence_record.id) + "\t" + str(len(str(sequence_record.seq))) + "\t" +
                            str(GC_percentage) + "\t" + str(No_ATGCN_count) + "\t" + str(N_count) + "\t" +
                            "H-C" + "\t" + "1" + "\t" + "-" + "\n")
    else:  #序列中包含N碱基数量大于1的情况
        index_first_N = str(sequence_record.seq).upper().find("N")  #查找第一个N位置
        N_location_list.append(index_first_N)

        #由于find功能不能查找序列中第二个或以后的N碱基的位置，故采取每确定一个N碱基的位置之后截取其后的所有碱基；
        #再一次find时就使用截取后的序列进行查找，依次类推，最终找到最后一个N碱基的位置。
        while N_count > 1:
            seq_cut_new = str(sequence_record.seq)[index_first_N + 1:len(sequence_record.seq) + 1]
            index_new_N = str(seq_cut_new).upper().find("N")
            index_first_N = index_first_N + 1 + index_new_N
            N_location_list.append(index_first_N)
            N_count -= 1

        #对每两个相邻的N碱基位置进行减法运算，计算每两个N碱基之间的位置差；
        #对这些位置差进行最大连续数统计，如果N碱基为连续，则前后位置差为1；
        #统计1连续出现最大次数，即为N碱基连续出现最大长度。
        N_location_list_minus = []
        for i in range(1, int(len(N_location_list))):
            N_location_list_minus.append(N_location_list[i] - N_location_list[i-1])
        max_time = 0
        cur_time = 1
        pre_element = None
        for j in N_location_list_minus:
            if j == pre_element:
                cur_time += 1
                max_time = max(cur_time, max_time)
            else:
                pre_element = j
                cut_time = j

        if N_percentage < 0.01:
            result_output.write(str(sequence_record.id) + "\t" + str(len(str(sequence_record.seq))) + "\t" +
                                str(GC_percentage) + "\t" + str(No_ATGCN_count) + "\t" +
                                str(str(sequence_record.seq.upper()).count("N")) + "\t" + "H-C" +
                                "\t" + str(max_time) + "\t" + str(round(max_time / int(len(N_location_list)), 2)) + "\n")
        elif N_percentage < 0.05:
            result_output.write(str(sequence_record.id) + "\t" + str(len(str(sequence_record.seq))) + "\t" +
                                str(GC_percentage) + "\t" + str(No_ATGCN_count) + "\t" +
                                str(str(sequence_record.seq.upper()).count("N")) + "\t" + "M-C" +
                                "\t" + str(max_time) + "\t" + str(round(max_time / int(len(N_location_list)), 2)) + "\n")
        else:
            result_output.write(str(sequence_record.id) + "\t" + str(len(str(sequence_record.seq))) + "\t" +
                                str(GC_percentage) + "\t" + str(No_ATGCN_count) + "\t" +
                                str(str(sequence_record.seq.upper()).count("N")) + "\t" + "L-C" +
                                "\t" + str(max_time) + "\t" + str(round(max_time / int(len(N_location_list)), 2)) + "\n")

result_output.close()
