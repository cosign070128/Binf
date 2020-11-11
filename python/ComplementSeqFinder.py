#coding=utf-8

#-----------------------------------------------------------------------------------------------------------------------
#    1.此程序用于判断从NCBI上获取的病毒基因组序列是否为相同方向，是否存在于其他基因组反向互补序列的情况。
#
#    2.此程序主要工作流程如下：
#    2.1.读取包含病毒基因组所有序列的FASTA文件，并将序列ID和序列分别存入列表中；
#    2.2.从列表中依次取出2条序列，并写入文件，以“序列1ID_vs_序列2ID”方式命名文件；
#    2.3.调用muscle读取文件并已“序列1ID_vs_序列2ID_result”方式命名结果文件（采用多进程方式同时运行多个muscle）；
#    2.4.读取比对结果文件，统计相同碱基个数，以两个序列的最大长度为分母进行相似度计算；
#    2.5.将比对结果写入文件，并在屏幕上显示相似度低于85%*的比对信息。
#        *一般情况下，同一种病毒不同基因组序列间的相似度应处于90%-100%之间，若出现低85%的情况，需核对序列信息。
#
#    3.此程序使用方法：
#      python ComplementSeqFinder.py input_file.fasta output_file.result number_of_process
#
#    4.作者：cosign  日期：2020.11.11
#-----------------------------------------------------------------------------------------------------------------------


from Bio import SeqIO
from multiprocessing import Pool
import sys
import subprocess
import os


def muscle_alignment(seq1_ID, seq1_seq, seq2_ID, seq2_seq, result_output_file):
    """
    two sequences global alignment by muscle through shell.
    :param seq1_ID: ID of seq1
    :param seq1_seq: sequence of seq1
    :param seq2_ID: ID of seq2
    :param seq2_seq: sequence of seq2
    :param result_output_file: seq1_ID vs seq2_ID : seq_identity
    :return: nothing
    """
    # create a new file with seq1_ID and seq2_ID as the file name and which contain seq1_seq and seq2_seq.
    file_touch = str(seq1_ID) + "_vs_" + str(seq2_ID)
    file = open(file_touch, "w")
    file.write(">" + str(seq1_ID) + "\n")
    file.write(str(seq1_seq) + "\n")
    file.write(">" + str(seq2_ID) + "\n")
    file.write(str(seq2_seq) + "\n")
    file.close()
    # alignment by muscle through shell.
    file_input = str(cwd) + "/" + str(file_touch)
    file_output = str(cwd) + "/" + str(file_touch) + "_result"
    muscle_input = "muscle -in " + str(file_input) + " -out " + str(file_output) + " -clw -quiet"
    child = subprocess.Popen(str(muscle_input), shell=True)
    child.wait()
    #star count and sequence identity calclation.
    star_count = 0
    with open(file_output, "r") as muscle_result:
        for line in muscle_result:
            star_count += line.count("*")
    seq_identity = star_count / max(len(str(seq1_seq)), len(str(seq2_seq)))
    #print information.
    with open(result_output_file, "a+") as result_output:
        result_output.write(str(seq1_ID) + "\t" + str(seq2_ID) + "\t" + str(round(seq_identity, 4)) + "\n")
        result_output.close()
    if seq_identity <= 0.85:
        print(str(seq1_ID) + "_vs_" + str(seq2_ID) + " : " +
              str(round(seq_identity, 4)) + "\033[31;1m !!! Identity Abnormal !!!\033[0m")
    #remove muscle input file and output file.
    os.remove(file_input)
    os.remove(file_output)
    return

cwd = os.getcwd()

seqs_file = sys.argv[1]
result_output_file = sys.argv[2]
number_of_process = sys.argv[3]

seqs_ID = []
seqs = []
for seq_record in SeqIO.parse(seqs_file, "fasta"):
    seqs_ID.append(seq_record.id)
    seqs.append(seq_record.seq)

pool = Pool(int(number_of_process))
for first_seq_index in range(len(seqs_ID)):
    for second_seq_index in range(first_seq_index + 1, len(seqs_ID)):
        seq1_ID = seqs_ID[first_seq_index]
        seq1_seq = seqs[first_seq_index]
        seq2_ID = seqs_ID[second_seq_index]
        seq2_seq = seqs[second_seq_index]
        pool.apply_async(func=muscle_alignment, args=(seq1_ID, seq1_seq, seq2_ID, seq2_seq, result_output_file))
pool.close()
pool.join()
print("---------Alignments finished---------")
