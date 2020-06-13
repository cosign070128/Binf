# -*- coding: utf-8 -*-

#此程序用于对HPV E6E7基因提取。
#手动提取HPV E6E7基因为以下步骤：
#1.下载HPV序列；2.获取HPV E6E7基因起始位置和终止位置；3.利用efetch工具对E6E7序列进行下载
#手动操作存在的问题：
#1.麻烦；2.部分序列只有包含部分E6E7序列，位置信息不正常，会导致提取失败；3.ID匹配时比较容易出问题。

#此程序实现过程：
#1.读取HPV某一个型别的全CDS文件；2.抓取包含E6E7信息的行并对ID/起始位置/中止位置进行提取；
#3.取E6E7ID交集并提取交集中的E6起始位置和E7终止位置；4.生成efetch命令，以供后续下载。

#作者：张迪骏  日期：2020-06-13

import sys

if len(sys.argv) != 4:  ##参数数量不等于3的情况下，提示使用信息
    print("程序用法:python HPV_E6E7_GENE_EXTRACTER.py hpv_CDS_seqs_file efetch_output_file hpv_type")
    print("程序主要用于从HPV CDS文件中进行E6E7基因的提取。")
    print("作者：张迪骏  日期：2020.6.13")
    exit()

hpv_CDS_seqs = sys.argv[1]  #HPV CDS文件
E6E7_location_list = sys.argv[2]  #最终生成的efetch输入文件
hpv_type = sys.argv[3]  #HPV型别

E6E7_location_list_write = open(E6E7_location_list, 'w')

E6_seqs_ID = []  #建立E6 ID存放用空列表
E6_seqs_start = []  #建立E6基因起始位置存放用空列表
E7_seqs_ID = []  #建立E7 ID存放用空列表
E7_seqs_end = []  #建立E7基因终止位置存放用空列表

with open(hpv_CDS_seqs, "r") as hpv_CDS_seqs_file:
    for line in hpv_CDS_seqs_file:
        '''
        以下代码为E6E7基因在原始CDS文件中的提取过程
        由于在原始文件中对基因命名方式的不统一，故通过多种关键词进行抓取，包括"=E6" " E6" "E6]" "E6 "
        在提取位置时，如果基因完整则为数字，而不完整的基因则会添加">" 或 "<" 或 "complement"等非数字文字
        故在提取时使用try进行异常判断，若出现对位置信息进行整数型转化失败的情况，则直接跳过这组数据。
        '''
        if "=E6" in line or " E6" in line or "E6]" in line or "E6 " in line:
            try:
                if isinstance(int(line.split('location=')[1].split('.')[0]), int) == True:  #整数型判断
                    E6_seqs_ID.append(line.split(' ')[0].split('|')[1].split('_')[0])  #ID提取
                    E6_seqs_start.append(line.split('location=')[1].split('.')[0])  #起始位置提取
            except:  #出现异常跳过，继续循环
                pass

        if "=E7" in line or " E7" in line or "E7]" in line or "E7 " in line:
            try:
                if isinstance(int(line.split('location=')[1].split('.')[2].split(']')[0]), int) == True:  #整数型判断
                    E7_seqs_ID.append(line.split(' ')[0].split('|')[1].split('_')[0])  #ID提取
                    E7_seqs_end.append(line.split('location=')[1].split('.')[2].split(']')[0])  #终止位置提取
            except:  #出现异常跳过，继续循环
                pass

E6E7_seqs_ID = list(set(E6_seqs_ID) & set(E7_seqs_ID))  #E6E7 ID交集提取
E6E7_seqs_start = []  #建立E6E7基因起始位置存放用空列表，该位置即为E6起始位置
E6E7_seqs_end = []  #建立E6E7基因终止位置存放用空列表，该位置即为E7终止位置

for j in range(len(E6E7_seqs_ID)):
    E6E7_seqs_start.append(E6_seqs_start[E6_seqs_ID.index(E6E7_seqs_ID[j])])  #提取E6起始位置
    E6E7_seqs_end.append(E7_seqs_end[E7_seqs_ID.index(E6E7_seqs_ID[j])])  #提取E7终止位置

E6E7_location_list_write.write("#!/bin/bash" + "\n")  #写入文件首行Shebang

for k in range(len(E6E7_seqs_ID)):
    # 对E6E7长度进行判断，在CDS文件中发现存在一个条目中gene=E6 protein=E7的错误注释信息，对结果影响较大
    if int(E6E7_seqs_end[k]) - int(E6E7_seqs_start[k]) >= 600:
        E6E7_location_list_write.write("efetch -db nucleotide -format fasta -id " + E6E7_seqs_ID[k] +
                                       " -seq_start " + E6E7_seqs_start[k] +
                                       " -seq_stop " + E6E7_seqs_end[k] +
                                       " >> " + hpv_type + "\n")

E6E7_location_list_write.close()  #文件流关闭
