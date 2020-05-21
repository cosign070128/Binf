#此脚本主要用于简单的本地primer_blast计算，主要利用引物比对outfmt6格式的比对结果文件，提取引物在基因组序列上的结合位置，并计算产物长度。
#利用这些产物长度可以用于对扩增子稳定评价，也可以弥补单纯对引物进行比对检查的缺陷。
#实现过程：
#1：读取引物列表文件，获得引物ID，并对引物ID进行排序，将引物按照正向在前反向在后的顺序排列，便于遍历引物ID；
#2：遍历引物ID，在比对结果文件中提取包含引物ID的行，构建6个列表，按照相同索引位置储存正向引物命中的基因组编号/命中起始位置/命中终止位置和反向引物命中的基因组编号/命中起始位置/命中终止位置；
#3：将正反向引物命中基因组编号的列表取交集，并形成一个新的列表，遍历这个新列表，通过新列表的值返回去找到基因组编号所在的索引位置，并通过这个索引位置找到正反向引物的四个位置，并将其构建成一个列表；
#4：对这个位置列表取最大值和最小值，相减并加1即为产物长度，将每一对引物的产物长度构建成一个新的列表，计算该列表的最大值/最小值/平均值/长度；
#5：根据产物长度信息评价引物。
#此脚本存在的缺陷：
#1：此脚本进行计算的前提条件是假定所有命中的引物都能够扩增，不进行3'端错配判断，故需要在比对的参数上作出限定。
#2：由于每取一组引物就需要遍历一次比对结果文件，故不太适合用于引物分析前期，分析速度相对较慢。
#3：目前只适用于保守性引物分析，暂不能对非特异性扩增进行分析。

import sys
from Bio import SeqIO
from numpy import *

if len(sys.argv) != 4:  ##参数数量不等于4的情况下，提示使用信息
    print("Usage:python LocalPrimerBlast.py primer_list.fasta blastn_outfmt6_result_file output_file")
    print("Author: Zhangdijun  Date: 2020.5.21")

else:
    primers_list = sys.argv[1]
    blastn_outfmt_6 = sys.argv[2]
    calc_result_output = sys.argv[3]

    calc_output = open(calc_result_output, 'w')
    calc_output.write("F_primer" + "\t" + "R_primer" + "\t" +
                      "max_amplicon_length" + "\t" + "min_amplicon_length" + "\t" +
                      "mean_amplicon_length" + "\t" + "num_amplicon_length" + "\t" + "judgement" + "\n")

    #将引物ID添加至列表中，并通过排序将正反向引物按照正向引物在先反向引物在后的顺序排列，便于遍历引物ID
    primer_ID = []
    for primer_record in SeqIO.parse(primers_list, "fasta"):
        primer_ID.append(str(primer_record.id))
    primer_ID.sort()

    for i in range(int(len(primer_ID)/2)):
        with open(blastn_outfmt_6, "r") as blastn_result:
            #设立6个列表，分别为正向引物的命中基因组/基因组命中起始位置/基因组命中终止位置和反向引物的命中基因组/基因组命中起始位置/基因组命中终止位置
            genome_hit_F = []
            genome_F_location_1 = []
            genome_F_location_2 = []
            genome_hit_R = []
            genome_R_location_1 = []
            genome_R_location_2 = []

            for line in blastn_result:
                if primer_ID[i*2] in line: #取正向引物所在行信息
                    genome_hit_F.append(line.strip().split('\t')[1])
                    genome_F_location_1.append(line.strip().split('\t')[8])
                    genome_F_location_2.append(line.strip().split('\t')[9])
                if primer_ID[i*2+1] in line: #取反向引物所在行信息
                    genome_hit_R.append(line.strip().split('\t')[1])
                    genome_R_location_1.append(line.strip().split('\t')[8])
                    genome_R_location_2.append(line.strip().split('\t')[9])
            #取正反向引物命中基因组交集
            genome_hit_FR = list(set(genome_hit_F) & set(genome_hit_R))

            amplicon_length = []
            for j in range(len(genome_hit_FR)):
                #遍历交集中的基因组，并提取正向引物/反向引物在基因组的上位置
                location_list = []
                location_list.append(int(genome_F_location_1[genome_hit_F.index(genome_hit_FR[j])]))
                location_list.append(int(genome_F_location_2[genome_hit_F.index(genome_hit_FR[j])]))
                location_list.append(int(genome_R_location_1[genome_hit_R.index(genome_hit_FR[j])]))
                location_list.append(int(genome_R_location_2[genome_hit_R.index(genome_hit_FR[j])]))
                amplicon_length.append(int(max(location_list)) - int(min(location_list)) + 1) #产物长度提取
                #以下命令可用于异常数据检查
                #calc_output.write(primer_ID[i*2] + "\t" + genome_hit_FR[j] + "\t" + str(amplicon_length[j]) + "\t" + str(location_list) + "\n")

            #将结果写入输出文件
            if max(amplicon_length) == min(amplicon_length):
                calc_output.write(primer_ID[i*2] + "\t" + primer_ID[i*2+1] + "\t" +
                                  str(max(amplicon_length)) + "\t" + str(min(amplicon_length)) + "\t" + str(round(mean(amplicon_length), 1)) + "\t" +
                                  str(len(amplicon_length)) + "\t" + "0" + "\n")  #产物长度无异常，使用0标注
            else:
                calc_output.write(primer_ID[i*2] + "\t" + primer_ID[i*2+1] + "\t" +
                                  str(max(amplicon_length)) + "\t" + str(min(amplicon_length)) + "\t" + str(round(mean(amplicon_length), 1)) + "\t" +
                                  str(len(amplicon_length)) + "\t" + "1" + "\n")  #产物长度异常，使用1标注
