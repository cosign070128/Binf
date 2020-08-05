#coding=utf-8

#1.此脚本主要功能是查找反向引物在病毒序列中的位置；
#2.此脚本在整个病毒基因组引物设计流程中处于PragmentsEva分析后；
#3.若无此脚本，可手动采用Oligo7，确定反向引物3'端后，进行正向引物设计；
#4.使用此脚本后，可将反向引物3'端输入到primer3中，进行自动化的正向引物设计，通量较高。

import sys
from Bio import SeqIO

if len(sys.argv) != 4:  ##参数数量不等于4的情况下，提示使用信息
    print("Usage:python R_primer_location_finder.py primer_list.fasta refseq_template_file.fasta result_file")
    print("Author: Zhangdijun  Date: 2020.8.5")
    exit()

if len(sys.argv) == 4:
    primer_list = sys.argv[1]
    template = sys.argv[2]
    location_list = sys.argv[3]

    result_output = open(location_list, 'w')
    result_output.write("primer_ID" + "\t" + "primer_seq" + "\t" + "primer_seq_rc" + "\t" + "R_primer_location" + "\n")

    #如下循环获取引物列表中引物的ID/序列/反向序列
    primer_ID = []
    primer_seq = []
    primer_seq_rc = []
    for primer_record in SeqIO.parse(primer_list, "fasta"):
        primer_ID.append(str(primer_record.id))
        primer_seq.append(str(primer_record.seq))
        primer_seq_rc.append((str(primer_record.seq.reverse_complement())))

    #如下循环获取refseq序列信息
    for template_record in SeqIO.parse(template, "fasta"):
        template_seq = str(template_record.seq)

    #如下循环获取每一个反向引物在模版上的位置
    for number in range(len(primer_ID)):
        R_primer_location = template_seq.find(primer_seq_rc[number])
        result_output.write((str(primer_ID[number])) + "\t" + str(primer_seq[number]) + "\t" + str(primer_seq_rc[number]) +
                            "\t" + str(R_primer_location) + "\n")
    result_output.close()
