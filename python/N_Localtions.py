#此程序用于检查序列中存在的N是否为连续分布或为分散型分布。
#通过对N碱基在序列中的分布情况判断序列质量。
#作者：张迪骏 日期：2020.06.30

import sys
from Bio import SeqIO

if len(sys.argv) != 3:  ##参数数量不等于3的情况下，提示使用信息
    print("脚本用法:python N_Locations.py sequences.fasta output_file")
    print("此脚本主要用于序列中N碱基的分布情况，判断序列质量。")
    print("作者：张迪骏  日期：2020.6.30")
    exit()

sequences_input = sys.argv[1]  #输入文件
N_analysis_result_file = sys.argv[2]  #输出文件

#写入文件表头
N_count_result = open(N_analysis_result_file, "w")
N_count_result.write("Seq_ID" + "\t" + "N_count" + "\t" + "max_Ns_length" + "\t" + "max_Ns_length/N_count" + "\n")

#依次遍历所有序列
for sequence_record in SeqIO.parse(sequences_input, "fasta"):
    N_location_list = []  #每一个序列开始读取清空N碱基位置列表
    N_count = str(sequence_record.seq).count("N")  #对序列中N碱基数量进行计数

    if N_count < 1:  #序列中不包含N碱基的情况
        N_count_result.write(str(sequence_record.id) + "\t" + "0" + "\t" + "0" + "\t" + "-" + "\n")
    elif N_count == 1:  #序列中包含N碱基的情况
        N_count_result.write(str(sequence_record.id) + "\t" + "1" + "\t" + "1" + "\t" + "-" + "\n")
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
        #结果写入文件
        N_count_result.write(str(sequence_record.id) + "\t" + str(len(N_location_list)) + "\t" + str(max_time) +
                             "\t" + str(round(max_time / int(len(N_location_list)), 2)) + "\n")

N_count_result.close()
