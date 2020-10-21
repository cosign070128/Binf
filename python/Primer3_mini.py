#coding=utf-8

#-----------------------------------------------------------------------------------------------------------------------
#
# 1.此程序用于病毒靶点的引物设计。
#
# 2.此程序设计原理为：
#    2.1.  读取输入的fasta文件，选择一条不包含N碱基且序列长度最短的序列作为设计参考序列，对应函数reference_sequence_select；
#    2.2.  对设计参考序列进行序列切割（正/反向），以30bp长度为标准切割出每一个3端（类引物），对应函数pragmenter；
#    2.3.  对类引物进行Tm值判断，使类引物截短至设计所需Tm值范围内（默认为59-61），对应函数pragmenters_Tm_calculation;
#    2.4.  对引物的3端5个碱基进行碱基组成进行判断（默认不多于3个GC且大于1个GC），对应函数primer_last_5_base_calculation；
#    2.5.  对引物是否存在超过n个以上连续碱基的情况进行判断(默认为n=3），对应函数primer_poly_X_calculation；
#    2.6.1 对满足上述条件的引物进行自身二聚体判断,只判断全匹配可双向延伸类型，对应函数primer_homodimer_check；
#    2.6.2 对满足上述条件的引物进行自身二聚体判断，采用primer3.py—calcHomodimer计算，对应函数primer_homodimer_check_by_primer3;
#    2.7.  对满足上述条件的引物进行发夹结构判断，对应函数primer_hairpin_check_by_primer3；
#    2.8.  对通过检查的引物进行位置索引；
#    2.9.  配对正反向引物，确保引物长度满足需求；
#
#    ****----****----****----****----****----****----以下功能待验证----****----****----****----****----****----****----***
#    2.10. 检查配对正反向引物之间的相互作用；
#    ****----****----****----****----****----****----****----****----****----****----****----****----****----****----***
#
#    2.11. 输出所有满足要求的配对引物以及产物（产物序列中大写表示引物位置，小写表示扩增序列）。
#
# 3.程序使用方法：
#    python Primer3_mini.py 
#                  [sequence.fasta]   #包含所有序列的FASTA文件
#                  [min_Tm-max_Tm]   #引物Tm值范围
#                  [max_poly_X]   #最大连续碱基个数
#                  [Tm]   #热力学检查Tm值
#                  [min_amplicon_length-max_amplicon_length]   #产物长度
#                  [result_file]   #满足设计条件的引物和产物输出文件
#                  [primer_eliminated_list_file]   #不满足设计条件的引物输出文件
#
# 4.作者：张迪骏     时间：2020.09.30
#-----------------------------------------------------------------------------------------------------------------------

import sys
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.Alphabet import IUPAC
import primer3

def reference_sequence_select(sequence_file):
    """此函数用于设计参考序列的选择。选择条件为：不含N碱基且序列长度最短的序列"""
    sequence_ID = []
    sequence_length = []
    sequence = []

    for sequence_record in SeqIO.parse(sequence_file, "fasta"):
        if sequence_record.seq.count("N") == 0:
            sequence_ID.append(str(sequence_record.id))
            sequence_length.append(str(len(sequence_record.seq)))
            sequence.append(str(sequence_record.seq))

    sequence_min_length = min(sequence_length)  # 获取最短序列长度
    
    # 在索引过程中，可能会碰到有两个均为最短序列长度的序列，由于索引机制，选择第一个被索引到的序列作为设计参考序列
    reference_sequence_selected = {'sequence_ID': sequence_ID[sequence_length.index(sequence_min_length)],
                                   'sequence': sequence[sequence_length.index(sequence_min_length)]}

    # 由于需要返回序列名称和序列碱基，故构建字典返回值
    return reference_sequence_selected

def pragmenter(reference_sequence_selected):
    """此函数用于序列切割，分别对正向和反向互补序列进行切割，切割长度30bp"""

    def pragmenter_inner(sequence_ID, sequence):
        sequence_pragmenters_ID = []
        sequence_pragmenters = []
        for position in range(0, int(len(sequence))-30):
            sequence_pragmenters_ID.append(sequence_ID + "_" + str(position))
            sequence_pragmenters.append(sequence[position:position+30])
        pragmenters_inner_result = [sequence_pragmenters_ID, sequence_pragmenters]
        return pragmenters_inner_result

    sequence_ID = str(reference_sequence_selected['sequence_ID'])  # 后一个sequence_ID为参考序列选择字典返回值中的键
    sequence = str(reference_sequence_selected['sequence'])  #后一个sequence为参考序列选择字典返回值中的键
    sequence_reverse_complement_ID = str(sequence_ID + "_RC_")
    sequence_reverse_complement = str(Seq(sequence, IUPAC.ambiguous_dna).reverse_complement())

    sequence_pragmenters_result = pragmenter_inner(sequence_ID, sequence)
    sequence_reverse_complement_pragmenters_result = pragmenter_inner(sequence_reverse_complement_ID,
                                                                      sequence_reverse_complement)

    pragmenters_result = [sequence_pragmenters_result, sequence_reverse_complement_pragmenters_result]
    #pragmenters_result索引方式:
    #pragment_result[0][0][i] = sequence_pragments_ID_i;
    #pragment_result[0][1][i] = sequence_pragmenters_i;
    #pragment_result[1][0][i] = sequence_reverse_complement_pragmenters_ID_i;
    #pragment_result[1][1][i] = sequence_reverse_complement_pragmenters_i;

    return pragmenters_result

def pragmenters_Tm_calculation(pragmenters, min_Tm = 59, max_Tm = 61):
    """此函数用于对切割后的类引物进行Tm值计算，通过碱基的删减，输出所有满足和不满足预设Tm值的类引物"""
    sequence_pragmenters_id = pragmenters[0]
    sequence_pragmenters = pragmenters[1]
    sequence_RC_pragments_id = pragmenters[2]
    sequence_RC_pragments = pragmenters[3]

    # 以下为正向类引物Tm值计算
    sequence_pragmenters_id_satisfied = []  # 列表用于收录Tm值满足最小Tm值和最大Tm值类引物的ID
    sequence_pragmenters_Tm_satisfied = []  # 列表用于收录Tm值满足最小Tm值和最大Tm值类引物的序列
    sequence_pragmenters_id_not_satisfied = []  #列 表用于收录Tm值不满足最小Tm值和最大Tm值类引物的ID
    sequence_pragmenters_Tm_not_satisfied = []  # 列表用于收录Tm值不满足最小Tm值和最大Tm值类引物的序列

    for id, pragmenter in zip(sequence_pragmenters_id, sequence_pragmenters):
        for split in range(9):
            pragmenter_Tm = primer3.calcTm(str(pragmenter[split:]), mv_conc=50, dv_conc=3, dna_conc=200)
            if min_Tm <= pragmenter_Tm <= max_Tm:
                sequence_pragmenters_id_satisfied.append(str(id) + "_" + str(split) + "@" +
                                                         str(round(pragmenter_Tm, 2)))
                sequence_pragmenters_Tm_satisfied.append(pragmenter[split:])
            else:
                sequence_pragmenters_id_not_satisfied.append(str(id) + "_" + str(split) + "@" +
                                                             str(round(pragmenter_Tm, 2)))
                sequence_pragmenters_Tm_not_satisfied.append(pragmenter[split:])

    # 以下为反向类引物Tm值计算
    sequence_RC_pragmenters_id_satisfied = []  # 列表用于收录Tm值满足最小Tm值和最大Tm值类引物的ID
    sequence_RC_pragmenters_Tm_satisfied = []  # 列表用于收录Tm值满足最小Tm值和最大Tm值类引物的序列
    sequence_RC_pragmenters_id_not_satisfied = []  # 列表用于收录Tm值不满足最小Tm值和最大Tm值类引物的ID
    sequence_RC_pragmenters_Tm_not_satisfied = []  # 列表用于收录Tm值不满足最小Tm值和最大Tm值类引物的序列

    for RC_id, RC_pragmenter in zip(sequence_RC_pragments_id, sequence_RC_pragments):
        for RC_split in range(9):
            RC_pragmenter_Tm = primer3.calcTm(str(RC_pragmenter[RC_split:]), mv_conc=50, dv_conc=3, dna_conc=200)
            if min_Tm <= RC_pragmenter_Tm <= max_Tm:
                sequence_RC_pragmenters_id_satisfied.append(str(RC_id) + "_" + str(RC_split) + "@" +
                                                            str(round(RC_pragmenter_Tm, 2)))
                sequence_RC_pragmenters_Tm_satisfied.append(RC_pragmenter[RC_split:])
            else:
                sequence_RC_pragmenters_id_not_satisfied.append(str(RC_id) + "_" + str(RC_split) + "@" +
                                                                str(round(RC_pragmenter_Tm, 2)))
                sequence_RC_pragmenters_Tm_not_satisfied.append(RC_pragmenter[RC_split:])

    pragmenters_Tm_calculation_result = [sequence_pragmenters_id_satisfied,
                                         sequence_pragmenters_Tm_satisfied,
                                         sequence_RC_pragmenters_id_satisfied,
                                         sequence_RC_pragmenters_Tm_satisfied,
                                         sequence_pragmenters_id_not_satisfied,
                                         sequence_pragmenters_Tm_not_satisfied,
                                         sequence_RC_pragmenters_id_not_satisfied,
                                         sequence_RC_pragmenters_Tm_not_satisfied]

    return pragmenters_Tm_calculation_result

def primer_last_5_base_calculation(pragmenters_Tm_calculation_result):
    """此函数用于检查引物3端最后5个碱基的组成情况，一般情况下，最后5个碱基允许出现的GC碱基个数应小于等于3个且大于等于1个，若出现过多GC碱
    基，会导致3端引物结合能力过强，错配概率大幅度增加，若没有GC碱基，则会导致引物3端结合能力过低，影响引导效率，一般而言，引物设计应追求碱
    基的平均分布，尽可能不要出现GC/AT碱基连续出现的情况。"""
    primer_id = pragmenters_Tm_calculation_result[0]
    primer = pragmenters_Tm_calculation_result[1]
    primer_RC_id = pragmenters_Tm_calculation_result[2]
    primer_RC = pragmenters_Tm_calculation_result[3]

    # 以下为正向引物3端5个碱基GC含量判断
    primer_id_last_5_base_GC_count_satisfied = []  # 列表用于收录3端5个碱基GC含量满足要求的引物的ID
    primer_last_5_base_GC_count_satisfied = []  # 列表用于收录3端5个碱基GC含量满足要求的引物的序列
    primer_id_last_5_base_GC_count_not_satisfied = []  # 列表用于收录3端5个碱基GC含量不满足要求的引物的ID
    primer_last_5_base_GC_count_not_satisfied = []  # 列表用于收录3端5个碱基GC含量不满足要求的引物的序列

    for id, primer_seq in zip(primer_id, primer):
        primer_last_5_base = str(primer_seq[-5:])
        GC_count = primer_last_5_base.count("G") + primer_last_5_base.count("C")
        if 1 <= GC_count <= 3:
            primer_id_last_5_base_GC_count_satisfied.append(id)
            primer_last_5_base_GC_count_satisfied.append(primer_seq)
        else:
            primer_id_last_5_base_GC_count_not_satisfied.append(id)
            primer_last_5_base_GC_count_not_satisfied.append(primer_seq)

    # 以下为反向引物3端5个碱基GC含量判断
    primer_RC_id_last_5_base_GC_count_satisfied = []  # 列表用于收录3端5个碱基GC含量满足要求的引物的ID
    primer_RC_last_5_base_GC_count_satisfied = []  # 列表用于收录3端5个碱基GC含量满足要求的引物的序列
    primer_RC_id_last_5_base_GC_count_not_satisfied = []  # 列表用于收录3端5个碱基GC含量不满足要求的引物的ID
    primer_RC_last_5_base_GC_count_not_satisfied = []  # 列表用于收录3端5个碱基GC含量不满足要求的引物的序列

    for RC_id, primer_RC_seq in zip(primer_RC_id, primer_RC):
        primer_RC_last_5_base = str(primer_RC_seq[-5:])
        GC_RC_count = primer_RC_last_5_base.count("G") + primer_RC_last_5_base.count("C")
        if 1 <= GC_RC_count <= 3:
            primer_RC_id_last_5_base_GC_count_satisfied.append(RC_id)
            primer_RC_last_5_base_GC_count_satisfied.append(primer_RC_seq)
        else:
            primer_RC_id_last_5_base_GC_count_not_satisfied.append(RC_id)
            primer_RC_last_5_base_GC_count_not_satisfied.append(primer_RC_seq)

    primer_last_5_base_calculation_result = [primer_id_last_5_base_GC_count_satisfied,
                                             primer_last_5_base_GC_count_satisfied,
                                             primer_RC_id_last_5_base_GC_count_satisfied,
                                             primer_RC_last_5_base_GC_count_satisfied,
                                             primer_id_last_5_base_GC_count_not_satisfied,
                                             primer_last_5_base_GC_count_not_satisfied,
                                             primer_RC_id_last_5_base_GC_count_not_satisfied,
                                             primer_RC_last_5_base_GC_count_not_satisfied]

    return primer_last_5_base_calculation_result

def primer_poly_X_calculation(primer_last_5_base_calculation_result, max_poly_X_number = 3):
    """此函数用于检查引物连续碱基数量，在PCR过程中应尽可能避免连续相同碱基（大于4个）的出现，若出现连续碱基，可能会出现Stutter峰，会影响
    结果判读，故在引物设计过程中，应避免引物中存在的连续相同碱基，函数默认不超过3个连续碱基"""
    primer_id = primer_last_5_base_calculation_result[0]
    primer_seq = primer_last_5_base_calculation_result[1]
    primer_RC_id = primer_last_5_base_calculation_result[2]
    primer_RC_seq = primer_last_5_base_calculation_result[3]

    # 以下为正向引物连续碱基数量判断
    primer_id_poly_X_satisfied = []  # 列表用于收录连续碱基数量满足要求的引物的ID
    primer_poly_X_satisfied = []  # 列表用于收录连续碱基数量满足要求的引物的序列
    primer_id_poly_X_not_satisfied = []  # 列表用于收录连续碱基数量不满足要求的引物的ID
    primer_poly_X_not_satisfied = []  # 列表用于收录连续碱基数量不满足要求的引物的序列

    for id, seq in zip(primer_id, primer_seq):
        poly_A_count = seq.count("A" * (max_poly_X_number + 1))
        poly_T_count = seq.count("T" * (max_poly_X_number + 1))
        poly_C_count = seq.count("C" * (max_poly_X_number + 1))
        poly_G_count = seq.count("G" * (max_poly_X_number + 1))
        poly_X_count = poly_A_count + poly_T_count + poly_G_count + poly_C_count

        if poly_X_count == 0:
            primer_id_poly_X_satisfied.append(id)
            primer_poly_X_satisfied.append(seq)
        else:
            primer_id_poly_X_not_satisfied.append(id)
            primer_poly_X_not_satisfied.append(seq)

    # 以下为反向引物连续碱基数量判断
    primer_RC_id_poly_X_satisfied = []  # 列表用于收录连续碱基数量满足要求的引物的ID
    primer_RC_poly_X_satisfied = []  # 列表用于收录连续碱基数量满足要求的引物的序列
    primer_RC_id_poly_X_not_satisfied = []  # 列表用于收录连续碱基数量不满足要求的引物的ID
    primer_RC_poly_X_not_satisfied = []  # 列表用于收录连续碱基数量不满足要求的引物的序列

    for RC_id, RC_seq in zip(primer_RC_id, primer_RC_seq):
        RC_poly_A_count = RC_seq.count("A" * (max_poly_X_number + 1))
        RC_poly_T_count = RC_seq.count("T" * (max_poly_X_number + 1))
        RC_poly_C_count = RC_seq.count("C" * (max_poly_X_number + 1))
        RC_poly_G_count = RC_seq.count("G" * (max_poly_X_number + 1))
        RC_poly_X_count = RC_poly_A_count + RC_poly_T_count + RC_poly_G_count + RC_poly_C_count

        if RC_poly_X_count == 0:
            primer_RC_id_poly_X_satisfied.append(RC_id)
            primer_RC_poly_X_satisfied.append(RC_seq)
        else:
            primer_RC_id_poly_X_not_satisfied.append(RC_id)
            primer_RC_poly_X_not_satisfied.append(RC_seq)

    primer_poly_X_calculation_result = [primer_id_poly_X_satisfied,
                                        primer_poly_X_satisfied,
                                        primer_RC_id_poly_X_satisfied,
                                        primer_RC_poly_X_satisfied,
                                        primer_id_poly_X_not_satisfied,
                                        primer_poly_X_not_satisfied,
                                        primer_RC_id_poly_X_not_satisfied,
                                        primer_RC_poly_X_not_satisfied]

    return primer_poly_X_calculation_result

def primer_homodimer_check(primer_poly_X_calclation_result):
    """此函数用于自身二聚体检查（严格模式），主要筛查全匹配可双向延伸类型自身引物二聚体，即自身3端序列为回文序列
    例如：5-nnnnnnnnnATGCGCAT-3
                    ||||||||
                  3-TACGCGTAnnnnnnnnn-5
    检查通过序列本身是否等于序列本身反向互补进行判断"""

    primer_id = primer_poly_X_calclation_result[0]
    primer_seq = primer_poly_X_calclation_result[1]
    primer_RC_id = primer_poly_X_calclation_result[2]
    primer_RC_seq = primer_poly_X_calclation_result[3]

    # 以下为正向引物自身二聚体检查
    primer_id_homodimer_satisfied = []  # 列表用于收录通过自身二聚体严格模式检查的引物的ID
    primer_homodimer_satisfied = []  # 列表用于收录通过自身二聚体严格模式检查的引物序列
    primer_id_homodimer_not_satisfied = []  # 列表用于收录未通过自身二聚体严格模式检查的引物的ID
    primer_homodimer_not_satisfied = []  # 列表用于未收录通过自身二聚体严格模式检查的引物序列

    for id, seq in zip(primer_id, primer_seq):
        for length in range(4, len(seq) + 1, 2):  # 长度加1为确保引物全部都能够被检查到
            primer_sub_seq = str(seq[-length:])  # 截取3端序列
            primer_sub_seq_RC = str(Seq(primer_sub_seq, IUPAC.ambiguous_dna).reverse_complement())
            primer_sub_seq_AT_count = primer_sub_seq.count("A") + primer_sub_seq.count("T")
            if primer_sub_seq.strip() is primer_sub_seq_RC.strip():
                if length == 4 and primer_sub_seq_AT_count == 4:  # 3端最后4个碱基均为AT时结合能力叫差，较难形成引物二聚体，排除
                    continue
                else:
                    primer_id_homodimer_not_satisfied.append(str(id) + "_" + "HoD" + "@" + str(length))
                    primer_homodimer_not_satisfied.append(seq)
                    break
        primer_id_homodimer_satisfied.append(id)
        primer_homodimer_satisfied.append(seq)

    # 以下为反向引物自身二聚体检查
    primer_RC_id_homodimer_satisfied = []  # 列表用于收录通过自身二聚体严格模式检查的引物的ID
    primer_RC_homodimer_satisfied = []  # 列表用于收录通过自身二聚体严格模式检查的引物序列
    primer_RC_id_homodimer_not_satisfied = []  # 列表用于收录未通过自身二聚体严格模式检查的引物的ID
    primer_RC_homodimer_not_satisfied = []  # 列表用于未收录通过自身二聚体严格模式检查的引物序列

    for RC_id, RC_seq in zip(primer_RC_id, primer_RC_seq):
        for RC_length in range(4, len(RC_seq) + 1, 2):  # 长度加1为确保引物全部都能够被检查到
            primer_RC_sub_seq = str(RC_seq[-RC_length:])  # 截取3端序列
            primer_RC_sub_seq_RC = str(Seq(primer_RC_sub_seq, IUPAC.ambiguous_dna).reverse_complement())
            primer_RC_sub_seq_AT_count = primer_RC_sub_seq.count("A") + primer_RC_sub_seq.count("T")
            if primer_RC_sub_seq.strip() is primer_RC_sub_seq_RC.strip():
                if RC_length == 4 and primer_RC_sub_seq_AT_count == 4:  # 3端最后4个碱基均为AT时结合能力叫差，较难形成引物二聚体，排除
                    continue
                else:
                    primer_RC_id_homodimer_not_satisfied.append(str(RC_id) + "_" + "HoD" + "@" + str(RC_length))
                    primer_RC_homodimer_not_satisfied.append(RC_seq)
                    break
        primer_RC_id_homodimer_satisfied.append(RC_id)
        primer_RC_homodimer_satisfied.append(RC_seq)

    primer_homodimer_check_result = [primer_id_homodimer_satisfied,
                                     primer_homodimer_satisfied,
                                     primer_RC_id_homodimer_satisfied,
                                     primer_RC_homodimer_satisfied,
                                     primer_id_homodimer_not_satisfied,
                                     primer_homodimer_not_satisfied,
                                     primer_RC_id_homodimer_not_satisfied,
                                     primer_RC_homodimer_not_satisfied]

    return primer_homodimer_check_result

def primer_homodimer_check_by_primer3(primer_poly_X_calclation_result, Tm = 47):
    """此函数用于自身二聚体检查（primer3.py模式），通过热力学计算方式检查自身二聚体，在参数上较严格模式宽松，会杀掉一些非可延伸的二聚体，
    目前对于这些非可延伸的二聚体，手动操作过程中采取保留的方式。这项检查可以与严格模式相互补充。
    例如：    5-CCAGAGCTTAAGCTCTTTAGAAAT-3
                 ||||||||||||||
        3-TAAAGATTCTCGAATTCGAGACC-5
    默认二聚体Tm值需低于47摄氏度。"""

    primer_id = primer_poly_X_calclation_result[0]
    primer_seq = primer_poly_X_calclation_result[1]
    primer_RC_id = primer_poly_X_calclation_result[2]
    primer_RC_seq = primer_poly_X_calclation_result[3]

    # 以下为正向引物自身二聚体检查
    primer_id_homodimer_satisfied = []  # 列表用于收录通过自身二聚体严格模式检查的引物的ID
    primer_homodimer_satisfied = []  # 列表用于收录通过自身二聚体严格模式检查的引物序列
    primer_id_homodimer_not_satisfied = []  # 列表用于收录未通过自身二聚体严格模式检查的引物的ID
    primer_homodimer_not_satisfied = []  # 列表用于未收录通过自身二聚体严格模式检查的引物序列

    for id, seq in zip(primer_id, primer_seq):
        homodimer = primer3.calcHomodimer((str(seq)), mv_conc=50.0, dv_conc=3, dna_conc=200, temp_c=25)
        homodimer_Tm = str(homodimer).split('tm=')[1].split(',')[0]
        if float(homodimer_Tm) <= float(Tm):
            primer_id_homodimer_satisfied.append(id)
            primer_homodimer_satisfied.append(seq)
        else:
            primer_id_homodimer_not_satisfied.append(str(id) + "_" + "HoD" + "@" + str(homodimer_Tm))
            primer_homodimer_not_satisfied.append(seq)

    # 以下为反向引物自身二聚体检查
    primer_RC_id_homodimer_satisfied = []  # 列表用于收录通过自身二聚体严格模式检查的引物的ID
    primer_RC_homodimer_satisfied = []  # 列表用于收录通过自身二聚体严格模式检查的引物序列
    primer_RC_id_homodimer_not_satisfied = []  # 列表用于收录未通过自身二聚体严格模式检查的引物的ID
    primer_RC_homodimer_not_satisfied = []  # 列表用于未收录通过自身二聚体严格模式检查的引物序列

    for RC_id, RC_seq in zip(primer_RC_id, primer_RC_seq):
        RC_homodimer = primer3.calcHomodimer((str(RC_seq)), mv_conc=50.0, dv_conc=3, dna_conc=200, temp_c=25)
        RC_homodimer_Tm = str(RC_homodimer).split('tm=')[1].split(',')[0]
        if float(RC_homodimer_Tm) <= float(Tm):
            primer_RC_id_homodimer_satisfied.append(RC_id)
            primer_RC_homodimer_satisfied.append(RC_seq)
        else:
            primer_RC_id_homodimer_not_satisfied.append(str(RC_id) + "_" + "HoD" + "@" + str(RC_homodimer_Tm))
            primer_RC_homodimer_not_satisfied.append(RC_seq)

    primer_homodimer_check_by_primer3_result = [primer_id_homodimer_satisfied,
                                                primer_homodimer_satisfied,
                                                primer_RC_id_homodimer_satisfied,
                                                primer_RC_homodimer_satisfied,
                                                primer_id_homodimer_not_satisfied,
                                                primer_homodimer_not_satisfied,
                                                primer_RC_id_homodimer_not_satisfied,
                                                primer_RC_homodimer_not_satisfied]

    return primer_homodimer_check_by_primer3_result

def primer_hairpin_check_by_primer3(primer_homodimer_check_result, Tm = 47):
    """此函数用于检查引物发夹结构，采用primer3.py进行热力学计算，默认发夹结构Tm值47摄氏度。此方法会有一定程度的误杀。"""
    primer_id = primer_homodimer_check_result[0]
    primer_seq = primer_homodimer_check_result[1]
    primer_RC_id = primer_homodimer_check_result[2]
    primer_RC_seq = primer_homodimer_check_result[3]

    # 以下为正向引物发夹结构检查
    primer_id_hairpin_satisfied = []  # 列表用于收录通过发夹结构检查的引物的ID
    primer_hairpin_satisfied = []  # 列表用于收录通过发夹结构检查的引物序列
    primer_id_hairpin_not_satisfied = []  # 列表用于收录未通过发夹结构检查的引物的ID
    primer_hairpin_not_satisfied = []  # 列表用于收录未通过发夹结构检查的引物序列

    for id, seq in zip(primer_id, primer_seq):
        hairpin = primer3.calcHairpin(str(seq), mv_conc=50.0, dv_conc=3, dna_conc=200, temp_c=25)
        if float(hairpin.tm) <= float(Tm):
            primer_id_hairpin_satisfied.append(id)
            primer_hairpin_satisfied.append(seq)
        else:
            primer_id_hairpin_not_satisfied.append(str(id) + "_" + "hairpin" + "@" + str(round(hairpin.tm, 2)))
            primer_hairpin_not_satisfied.append(seq)

    # 以下为反向引物发夹结构检查
    primer_RC_id_hairpin_satisfied = []  # 列表用于收录通过发夹结构检查的引物的ID
    primer_RC_hairpin_satisfied = []  # 列表用于收录通过发夹结构检查的引物序列
    primer_RC_id_hairpin_not_satisfied = []  # 列表用于收录未通过发夹结构检查的引物的ID
    primer_RC_hairpin_not_satisfied = []  # 列表用于收录未通过发夹结构检查的引物序列

    for RC_id, RC_seq in zip(primer_RC_id, primer_RC_seq):
        RC_hairpin = primer3.calcHairpin(str(RC_seq), mv_conc=50.0, dv_conc=3, dna_conc=200, temp_c=25)
        if float(RC_hairpin.tm) <= float(Tm):
            primer_RC_id_hairpin_satisfied.append(RC_id)
            primer_RC_hairpin_satisfied.append(RC_seq)
        else:
            primer_RC_id_hairpin_not_satisfied.append(str(RC_id) + "_" + "hairpin" + "@" + str(round(RC_hairpin.tm, 2)))
            primer_RC_hairpin_not_satisfied.append(RC_seq)

    primer_hairpin_check_by_primer3_result = [primer_id_hairpin_satisfied,
                                              primer_hairpin_satisfied,
                                              primer_RC_id_hairpin_satisfied,
                                              primer_RC_hairpin_satisfied,
                                              primer_id_hairpin_not_satisfied,
                                              primer_hairpin_not_satisfied,
                                              primer_RC_id_hairpin_not_satisfied,
                                              primer_RC_hairpin_not_satisfied]

    return primer_hairpin_check_by_primer3_result





#-----------------------------------------------------------------------------------------------------------------------

file_test = sys.argv[1]

j = reference_sequence_select(file_test)
k = pragmenter(j)
for i in range(len(k[0][0])):
    print(str(">") + k[0][0][i])
    print(k[0][1][i])
for j in range(len(k[1][0])):
    print(str(">") + k[1][0][i])
    print(k[1][1][i])


#l = pragmenters_Tm_calculation(k, 59, 61)
#e = primer_last_5_base_calculation(l)
#o = primer_poly_X_calculation(e, 3)
#q = primer_homodimer_check(o)
#w = primer_homodimer_check_by_primer3(q)
'''
for kk in range(len(w[0])):
    print(">" + w[0][kk])
    print(str(w[1][kk]))
for ll in range(len(w[2])):
    print(">" + w[2][ll])
    print(str(w[3][ll]))

for jj in range(len(w[4])):
    print(">" + w[4][jj])
    print(str(w[5][jj]))
for qq in range(len(w[6])):
    print(">" + w[6][qq])
    print(str(w[7][qq]))
'''
#e = primer_hairpin_check_by_primer3(w)

#for kk in range(len(e[0])):
#    print(">" + e[0][kk])
#    print(str(e[1][kk]))
#for ll in range(len(e[2])):
#    print(">" + e[2][ll])
#    print(str(e[3][ll]))
'''
for jj in range(len(e[4])):
    print(">" + e[4][jj])
    print(str(e[5][jj]))
for qq in range(len(e[6])):
    print(">" + e[6][qq])
    print(str(e[7][qq]))
'''
