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
#    2.8.  对通过检查的引物进行位置索引,并输出满足产物长度且通过相互作用检查的配对引物，对应函数primer_pairing；
#    2.9.  输出所有满足要求的配对引物以及产物（产物序列中大写表示引物位置，小写表示扩增序列）。
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

    for sequence_record in SeqIO.parse(sequence_file, "fasta"):  # 将不含N碱基的序列添加进列表中，并统计长度
        if sequence_record.seq.count("N") == 0:
            sequence_ID.append(str(sequence_record.id))
            sequence_length.append(str(len(sequence_record.seq)))
            sequence.append(str(sequence_record.seq))

    sequence_min_length = min(sequence_length)  # 获取最短序列长度
    
    # 在索引过程中，可能会碰到有两个均为最短序列长度的序列，由于索引机制，选择第一个被索引到的序列作为设计参考序列
    reference_sequence_selected = {'sequence_ID': sequence_ID[sequence_length.index(sequence_min_length)],
                                   'sequence': sequence[sequence_length.index(sequence_min_length)]}

    return reference_sequence_selected

def pragmenter(reference_sequence_selected):
    """此函数用于序列切割，分别对正向和反向互补序列进行切割，切割长度30bp"""

    def pragmenter_inner(sequence_ID, sequence):
        sequence_pragmenters_ID = []
        sequence_pragmenters = []
        for position in range(0, int(len(sequence))-30):
            sequence_pragmenters_ID.append(sequence_ID + "_" + str(position))
            sequence_pragmenters.append(sequence[position:position+30]) #将序列切割成30bp的类引物
        pragmenters_inner_result = [sequence_pragmenters_ID, sequence_pragmenters]
        return pragmenters_inner_result

    sequence_ID = str(reference_sequence_selected['sequence_ID'])  # 后一个sequence_ID为参考序列选择字典返回值中的键
    sequence = str(reference_sequence_selected['sequence'])  # 后一个sequence为参考序列选择字典返回值中的键
    sequence_reverse_complement_ID = str(sequence_ID + "_RC")  # 反向互补序列命名
    sequence_reverse_complement = str(Seq(sequence, IUPAC.ambiguous_dna).reverse_complement()) # 反向互补序列

    sequence_pragmenters_result = pragmenter_inner(sequence_ID, sequence) # 正向类引物生成
    sequence_reverse_complement_pragmenters_result = pragmenter_inner(sequence_reverse_complement_ID,
                                                                      sequence_reverse_complement) # 反向类引物生成

    pragmenters_result = [sequence_pragmenters_result, sequence_reverse_complement_pragmenters_result]
    """
    pragmenters_result索引方式:
    pragment_result[0][0][i] = sequence_pragments_ID_i;
    pragment_result[0][1][i] = sequence_pragmenters_i;
    pragment_result[1][0][i] = sequence_reverse_complement_pragmenters_ID_i;
    pragment_result[1][1][i] = sequence_reverse_complement_pragmenters_i;
    """
    return pragmenters_result

def pragmenters_Tm_calculation(pragmenters_result, min_Tm = 59, max_Tm = 61):
    """此函数用于对切割后的类引物进行Tm值计算，通过碱基的删减，输出所有满足和不满足预设Tm值的类引物"""
    def pragmenters_Tm_calculation_inner(pragmenters_inner_ID, pragmenters, min_Tm_inner = 59, max_Tm_inner = 61):
        sequence_pragmenters_Tm_satisfied_ID = []
        sequence_pragmenters_Tm_satisfied = []
        sequence_pragmenters_Tm_not_satisfied_ID = []
        sequence_pragmenters_Tm_not_satisfied = []

        for ID, pragmenter in zip(pragmenters_inner_ID, pragmenters):
            for split in range(9):
            # 从30bp长度切到20bp长度，原因在于，如果序列大于30bp才能够达到预设最大Tm,表明序列中含有大量的AT;
            # 如果序列小于20bp就能够达到预设最小Tm,表明序列中含有大量的GC;以上两种情况均不满足引物中ATGC尽可能均匀分布的需求。
                pragmenter_Tm = primer3.calcTm(str(pragmenter[split:]), mv_conc=50, dv_conc=3, dna_conc=200)
                if float(min_Tm_inner) <= float(pragmenter_Tm) <= float(max_Tm_inner):  # 判断Tm值是否在最小最大Tm值之间
                    sequence_pragmenters_Tm_satisfied_ID.append(str(ID) + "_" + str(30 - split) + "bp@" +
                                                                str(round(pragmenter_Tm, 2)))
                    sequence_pragmenters_Tm_satisfied.append(str(pragmenter[split:]))
                else:
                    sequence_pragmenters_Tm_not_satisfied_ID.append(str(ID) + "_" + str(30 - split) + "bp@" +
                                                                    str(round(pragmenter_Tm, 2)))
                    sequence_pragmenters_Tm_not_satisfied.append(str(pragmenter[split:]))

        pragmenters_Tm_calculation_inner_result = [sequence_pragmenters_Tm_satisfied_ID,
                                                   sequence_pragmenters_Tm_satisfied,
                                                   sequence_pragmenters_Tm_not_satisfied_ID,
                                                   sequence_pragmenters_Tm_not_satisfied]

        return pragmenters_Tm_calculation_inner_result

    sequence_pragments_ID = pragmenters_result[0][0]
    sequence_pragments = pragmenters_result[0][1]

    sequence_pragmenters_Tm_calculation = \
        pragmenters_Tm_calculation_inner(sequence_pragments_ID,
                                         sequence_pragments,
                                         min_Tm, max_Tm) #正向类引物判断

    sequence_reverse_complement_pragmenters_ID = pragmenters_result[1][0]
    sequence_reverse_complement_pragmenters = pragmenters_result[1][1]

    sequence_reverse_complement_pragmenters_Tm_calculation = \
        pragmenters_Tm_calculation_inner(sequence_reverse_complement_pragmenters_ID,
                                         sequence_reverse_complement_pragmenters,
                                         min_Tm, max_Tm) #反向类引物判断

    pragmenters_Tm_calculation_result = [sequence_pragmenters_Tm_calculation,
                                         sequence_reverse_complement_pragmenters_Tm_calculation]
    """
    pragmenters_Tm_calculation_result索引方式：
    pragmenters_Tm_calculation[0][0][i] = sequence_pragmenters_Tm_satisfied_ID_i
    pragmenters_Tm_calculation[0][1][i] = sequence_pragmenters_Tm_satisfied_i
    pragmenters_Tm_calculation[0][2][i] = sequence_pragmenters_Tm_not_satisfied_ID_i
    pragmenters_Tm_calculation[0][3][i] = sequence_pragmenters_Tm_not_satisfied_i
    pragmenters_Tm_calculation[1][0][i] = sequence_reverse_complement_pragmenters_Tm_satisfied_ID_i
    pragmenters_Tm_calculation[1][1][i] = sequence_reverse_complement_pragmenters_Tm_satisfied_i
    pragmenters_Tm_calculation[1][2][i] = sequence_reverse_complement_pragmenters_Tm_not_satisfied_ID_i
    pragmenters_Tm_calculation[1][3][i] = sequence_reverse_complement_pragmenters_Tm_not_satisfied_i
    """

    return pragmenters_Tm_calculation_result

def primer_last_5_base_calculation(pragmenters_Tm_calculation_result):
    """此函数用于检查引物3端最后5个碱基的组成情况，一般情况下，最后5个碱基允许出现的GC碱基个数应小于等于3个且大于等于1个，若出现过多GC碱
    基，会导致3端引物结合能力过强，错配概率大幅度增加，若没有GC碱基，则会导致引物3端结合能力过低，影响引导效率，一般而言，引物设计应追求碱
    基的平均分布，尽可能不要出现GC/AT碱基连续出现的情况。至此筛选步骤，类引物已转化为正常引物。"""
    def primer_last_5_base_calculation_inner(primer_ID_inner, primers_inner):
        primer_last_5_base_GC_count_satisfied_ID = []
        primer_last_5_base_GC_count_satisfied = []
        primer_last_5_base_GC_count_not_satisfied_ID = []
        primer_last_5_base_GC_count_not_satisfied = []

        for ID, primer in zip(primer_ID_inner, primers_inner):
            primer_last_5_base = str(primer[-5:]) # 切出引物最后5个碱基
            GC_count = primer_last_5_base.count("G") + primer_last_5_base.count("C")
            if 1 <= GC_count <= 3:
                primer_last_5_base_GC_count_satisfied_ID.append(ID)
                primer_last_5_base_GC_count_satisfied.append(primer)
            else:
                primer_last_5_base_GC_count_not_satisfied_ID.append(str(ID) + "@GC_" + str(GC_count))
                primer_last_5_base_GC_count_not_satisfied.append(primer)

        primer_last_5_base_calculation_inner_result = [primer_last_5_base_GC_count_satisfied_ID,
                                                       primer_last_5_base_GC_count_satisfied,
                                                       primer_last_5_base_GC_count_not_satisfied_ID,
                                                       primer_last_5_base_GC_count_not_satisfied]

        return primer_last_5_base_calculation_inner_result

    primer_ID = pragmenters_Tm_calculation_result[0][0]
    primer = pragmenters_Tm_calculation_result[0][1]

    primer_last_5_base_calculation_inner_result = \
        primer_last_5_base_calculation_inner(primer_ID, primer) # 正向引物判断

    primer_reverse_complement_ID = pragmenters_Tm_calculation_result[1][0]
    primer_reverse_complement = pragmenters_Tm_calculation_result[1][1]

    primer_reverse_complement_last_5_base_calculation_inner_result = \
        primer_last_5_base_calculation_inner(primer_reverse_complement_ID, primer_reverse_complement) # 反向引物判断

    primer_last_5_base_calculation_result = [primer_last_5_base_calculation_inner_result,
                                             primer_reverse_complement_last_5_base_calculation_inner_result]

    """
    primer_last_5_base_calculation_result索引方式：
    primer_last_5_base_calculation_result[0][0][i] = primer_last_5_base_GC_count_satisfied_ID_i
    primer_last_5_base_calculation_result[0][1][i] = primer_last_5_base_GC_count_satisfied_i
    primer_last_5_base_calculation_result[0][2][i] = primer_last_5_base_GC_count_not_satisfied_ID_i
    primer_last_5_base_calculation_result[0][3][i] = primer_last_5_base_GC_count_not_satisfied_i
    primer_last_5_base_calculation_result[1][0][i] = primer_reverse_complement_last_5_base_GC_count_satisfied_ID_i
    primer_last_5_base_calculation_result[1][1][i] = primer_reverse_complement_last_5_base_GC_count_satisfied_i
    primer_last_5_base_calculation_result[1][2][i] = primer_reverse_complement_last_5_base_GC_count_not_satisfied_ID_i
    primer_last_5_base_calculation_result[1][3][i] = primer_reverse_complement_last_5_base_GC_count_not_satisfied_i
    """


    return primer_last_5_base_calculation_result

def primer_poly_X_calculation(primer_last_5_base_calculation_result, max_poly_X_number = 3):
    """此函数用于检查引物连续碱基数量，在PCR过程中应尽可能避免连续相同碱基（大于4个）的出现，若出现连续碱基，可能会出现Stutter峰，会影响
    结果判读，故在引物设计过程中，应避免引物中存在的连续相同碱基，函数默认不超过3个连续碱基"""
    def primer_poly_X_calculation_inner(primer_ID_inner, primer_inner):
        primer_poly_X_satisfied_ID = []
        primer_poly_X_satisfied = []
        primer_poly_X_not_satisfied_ID = []
        primer_poly_X_not_satisfied = []

        for ID, primer in zip(primer_ID_inner, primer_inner):
            # 统计poly_X总的个数
            poly_A_count = primer.count("A" * (max_poly_X_number + 1))
            poly_T_count = primer.count("T" * (max_poly_X_number + 1))
            poly_C_count = primer.count("C" * (max_poly_X_number + 1))
            poly_G_count = primer.count("G" * (max_poly_X_number + 1))
            poly_X_count = poly_A_count + poly_T_count + poly_G_count + poly_C_count 
            if poly_X_count == 0:
                primer_poly_X_satisfied_ID.append(ID)
                primer_poly_X_satisfied.append(primer)
            else:
                primer_poly_X_not_satisfied_ID.append(str(ID) + "@poly_X_" + str(poly_X_count))
                primer_poly_X_not_satisfied.append(primer)

        primer_poly_X_calculation_inner_result = [primer_poly_X_satisfied_ID,
                                                  primer_poly_X_satisfied,
                                                  primer_poly_X_not_satisfied_ID,
                                                  primer_poly_X_not_satisfied]

        return primer_poly_X_calculation_inner_result

    primer_ID = primer_last_5_base_calculation_result[0][0]
    primer = primer_last_5_base_calculation_result[0][1]

    primer_poly_X_calculation_inner_result = \
        primer_poly_X_calculation_inner(primer_ID, primer) # 正向引物判断

    primer_reverse_complement_ID = primer_last_5_base_calculation_result[1][0]
    primer_reverse_complement = primer_last_5_base_calculation_result[1][1]

    primer_reverse_complement_poly_X_calculation_inner_result = \
        primer_poly_X_calculation_inner(primer_reverse_complement_ID, primer_reverse_complement) # 反向引物判断

    primer_poly_X_calculation_result = [primer_poly_X_calculation_inner_result,
                                        primer_reverse_complement_poly_X_calculation_inner_result]
    """
    primer_poly_X_calculation_result索引方式：
    primer_poly_X_calculation_result[0][0][i] = primer_poly_X_satisfied_ID_i
    primer_poly_X_calculation_result[0][1][i] = primer_poly_X_satisfied_i
    primer_poly_X_calculation_result[0][2][i] = primer_poly_X_not_satisfied_ID_i
    primer_poly_X_calculation_result[0][3][i] = primer_poly_X_not_satisfied_i
    primer_poly_X_calculation_result[1][0][i] = primer_reverse_complement_poly_X_satisfied_ID_i
    primer_poly_X_calculation_result[1][1][i] = primer_reverse_complement_poly_X_satisfied_i
    primer_poly_X_calculation_result[1][2][i] = primer_reverse_complement_poly_X_not_satisfied_ID_i
    primer_poly_X_calculation_result[1][3][i] = primer_reverse_complement_poly_X_not_satisfied_i
    """

    return primer_poly_X_calculation_result

def primer_homodimer_check(primer_poly_X_calclation_result):
    """此函数用于自身二聚体检查（严格模式），主要筛查全匹配可双向延伸类型自身引物二聚体，即自身3端序列为回文序列
    例如：5'-nnnnnnnnnATGCGCAT-3'
                      ||||||||
                   3'-TACGCGTAnnnnnnnnn-5'
    检查通过序列本身是否等于序列本身反向互补进行判断"""
    def primer_homodimer_check_inner(primer_ID, primer):
        primer_homodimer_satisfied_ID = []
        primer_homodimer_satisfied = []
        primer_homodimer_not_satisfied_ID = []
        primer_homodimer_not_satisfied = []

        for ID, seq in zip(primer_ID, primer):
            time = 0 #统计homodimer次数
            for split_length in range(4, len(seq) + 1, 2): # length+1为确保引物全长在双数是被全长被取到，从3端4个碱基开始计算
                primer_sub_seq = str(seq[-split_length:])
                primer_sub_seq_RC = str(Seq(primer_sub_seq, IUPAC.ambiguous_dna).reverse_complement())
                primer_sub_seq_AT_count = primer_sub_seq.count("A") + primer_sub_seq.count("T")
                if primer_sub_seq == primer_sub_seq_RC:
                    if int(len(primer_sub_seq)) == int(primer_sub_seq_AT_count) == 4:
                        continue #如果homodimer为4个碱基且均为AT，其结合能力较差，故接受此种情况
                    else:
                        primer_homodimer_not_satisfied_ID.append(str(ID) + "@HoD_" + str(split_length))
                        primer_homodimer_not_satisfied.append(str(seq))
                        time += 1
            if time == 0:
                primer_homodimer_satisfied_ID.append(str(ID))
                primer_homodimer_satisfied.append(str(seq))

        primer_homodimer_check_inner_result = [primer_homodimer_satisfied_ID,
                                              primer_homodimer_satisfied,
                                              primer_homodimer_not_satisfied_ID,
                                              primer_homodimer_not_satisfied]

        return primer_homodimer_check_inner_result

    primer_ID = primer_poly_X_calclation_result[0][0]
    primer = primer_poly_X_calclation_result[0][1]

    primer_homodimer_check_inner_result = \
        primer_homodimer_check_inner(primer_ID, primer) # 正向引物判断

    primer_reverse_complement_ID = primer_poly_X_calclation_result[1][0]
    primer_reverse_complement = primer_poly_X_calclation_result[1][1]

    primer_reverse_complement_homodimer_check_inner_result = \
        primer_homodimer_check_inner(primer_reverse_complement_ID, primer_reverse_complement) # 反向引物判断

    primer_homodimer_check_result = [primer_homodimer_check_inner_result,
                                     primer_reverse_complement_homodimer_check_inner_result]
    """
    primer_homodimer_check_result索引方式：
    primer_homodimer_check_result[0][0][i] = primer_homodimer_satisfied_ID_i
    primer_homodimer_check_result[0][1][i] = primer_homodimer_satisfied_i
    primer_homodimer_check_result[0][2][i] = primer_homodimer_not_satisfied_ID_i
    primer_homodimer_check_result[0][3][i] = primer_homodimer_not_satisfied_i
    primer_homodimer_check_result[1][0][i] = primer_reverse_complement_homodimer_satisfied_ID_i
    primer_homodimer_check_result[1][1][i] = primer_reverse_complement_homodimer_satisfied_i
    primer_homodimer_check_result[1][2][i] = primer_reverse_complement_homodimer_not_satisfied_ID_i
    primer_homodimer_check_result[1][3][i] = primer_reverse_complement_homodimer_not_satisfied_i
    
    """
    return primer_homodimer_check_result

def primer_homodimer_check_by_primer3(primer_poly_X_calclation_result, Tm = 47):
    """此函数用于自身二聚体检查（primer3.py模式），通过热力学计算方式检查自身二聚体，在参数上较严格模式宽松，会杀掉一些非可延伸的二聚体，
    目前对于这些非可延伸的二聚体，手动操作过程中采取保留的方式。这项检查可以与严格模式相互补充。
    例如：    5‘-CCAGAGCTTAAGCTCTTTAGAAAT-3’
                   ||||||||||||||
         3‘-TAAAGATTCTCGAATTCGAGACC-5’
    默认二聚体Tm值需低于47摄氏度。"""
    def primer_homodimer_check_by_primer3_inner(primer_ID, primers, Tm):
        primer_homodimer_satisfied_ID = []
        primer_homodimer_satisfied = []
        primer_homodimer_not_satisfied_ID = []
        primer_homodimer_not_satisfied = []

        for ID, primer in zip(primer_ID, primers):
            homodimer = primer3.calcHomodimer((str(primer)), mv_conc=50.0, dv_conc=3, dna_conc=200, temp_c=25) # homodimer检查
            homodimer_Tm = str(homodimer).split('tm=')[1].split(',')[0] # Tm值提取
            if float(homodimer_Tm) <= float(Tm):
                primer_homodimer_satisfied_ID.append(ID)
                primer_homodimer_satisfied.append(primer)
            else:
                primer_homodimer_not_satisfied_ID.append(str(ID) + "@HoD_" + str(homodimer_Tm))
                primer_homodimer_not_satisfied.append(primer)

        primer_homodimer_check_by_primer3_inner_result = [primer_homodimer_satisfied_ID,
                                                          primer_homodimer_satisfied,
                                                          primer_homodimer_not_satisfied_ID,
                                                          primer_homodimer_not_satisfied]

        return primer_homodimer_check_by_primer3_inner_result

    primer_ID = primer_poly_X_calclation_result[0][0]
    primer = primer_poly_X_calclation_result[0][1]

    primer_homodimer_check_by_primer3_inner_result = \
        primer_homodimer_check_by_primer3_inner(primer_ID, primer, Tm) # 正向引物判断

    primer_reverse_complement_ID = primer_poly_X_calclation_result[1][0]
    primer_reverse_complement = primer_poly_X_calclation_result[1][1]

    primer_reverse_complement_homodimer_check_by_primer3_inner_result = \
        primer_homodimer_check_by_primer3_inner(primer_reverse_complement_ID, primer_reverse_complement, Tm) # 反向引物判断

    primer_homodimer_check_by_primer3_result = [primer_homodimer_check_by_primer3_inner_result,
                                                primer_reverse_complement_homodimer_check_by_primer3_inner_result]

    """
    primer_homodimer_check_by_primer3_result索引方式：
    primer_homodimer_check_by_primer3_result[0][0][i] = primer_homodimer_satisfied_ID
    primer_homodimer_check_by_primer3_result[0][1][i] = primer_homodimer_satisfied
    primer_homodimer_check_by_primer3_result[0][2][i] = primer_homodimer_not_satisfied_ID
    primer_homodimer_check_by_primer3_result[0][3][i] = primer_homodimer_not_satisfied
    primer_homodimer_check_by_primer3_result[1][0][i] = primer_reverse_complement_homodimer_satisfied_ID
    primer_homodimer_check_by_primer3_result[1][1][i] = primer_reverse_complement_homodimer_satisfied
    primer_homodimer_check_by_primer3_result[1][2][i] = primer_reverse_complement_homodimer_not_satisfied_ID
    primer_homodimer_check_by_primer3_result[1][3][i] = primer_reverse_complement_homodimer_not_satisfied
    """
    return primer_homodimer_check_by_primer3_result

def primer_hairpin_check_by_primer3(primer_homodimer_check_result, Tm = 47):
    """此函数用于检查引物发夹结构，采用primer3.py进行热力学计算，默认发夹结构Tm值47摄氏度。此方法会有一定程度的误杀。"""
    def primer_hairpin_check_by_primer3_inner(primers_ID, primers, Tm):
        primer_hairpin_satisfied_ID = []
        primer_hairpin_satisfied = []
        primer_hairpin_not_satisfied_ID = []
        primer_hairpin_not_satisfied = []

        for ID, primer in zip(primers_ID, primers):
            hairpin = primer3.calcHairpin(str(primer), mv_conc=50.0, dv_conc=3, dna_conc=200, temp_c=25 # 发夹结构检查
            if float(hairpin.tm) <= float(Tm):
                primer_hairpin_satisfied_ID.append(str(ID))
                primer_hairpin_satisfied.append(str(primer))
            else:
                primer_hairpin_not_satisfied_ID.append(str(ID) + "_hairpin@" + str(round(hairpin.tm, 2)))
                primer_hairpin_not_satisfied.append(str(primer))

        primer_hairpin_check_by_primer3_inner_result = [primer_hairpin_satisfied_ID,
                                                        primer_hairpin_satisfied,
                                                        primer_hairpin_not_satisfied_ID,
                                                        primer_hairpin_not_satisfied]

        return primer_hairpin_check_by_primer3_inner_result

    primer_ID = primer_homodimer_check_result[0][0]
    primer = primer_homodimer_check_result[0][1]

    primer_hairpin_check_by_primer3_inner_result = \
        primer_hairpin_check_by_primer3_inner(primer_ID, primer, Tm) # 正向引物筛选

    primer_reverse_complement_ID = primer_homodimer_check_result[1][0]
    primer_reverse_complement = primer_homodimer_check_result[1][1]

    primer_reverse_complement_hairpin_check_by_primer3_inner_result = \
        primer_hairpin_check_by_primer3_inner(primer_reverse_complement_ID, primer_reverse_complement, Tm) # 反向引物筛选

    primer_hairpin_check_by_primer3_result = [primer_hairpin_check_by_primer3_inner_result,
                                              primer_reverse_complement_hairpin_check_by_primer3_inner_result]

    """
    primer_hairpin_check_by_primer3_result索引方式：
    primer_hairpin_check_by_primer3_result[0][0][i] = primer_hairpin_satisfied_ID_i
    primer_hairpin_check_by_primer3_result[0][1][i] = primer_hairpin_satisfied_i
    primer_hairpin_check_by_primer3_result[0][2][i] = primer_hairpin_not_satisfied_ID_i
    primer_hairpin_check_by_primer3_result[0][3][i] = primer_hairpin_not_satisfied_i
    primer_hairpin_check_by_primer3_result[1][0][i] = primer_reverse_complement_hairpin_satisfied_ID_i
    primer_hairpin_check_by_primer3_result[1][1][i] = primer_reverse_complement_hairpin_satisfied_i
    primer_hairpin_check_by_primer3_result[1][2][i] = primer_reverse_complement_hairpin_not_satisfied_ID_i
    primer_hairpin_check_by_primer3_result[1][3][i] = primer_reverse_complement_hairpin_not_satisfied_i
    
    """
    return primer_hairpin_check_by_primer3_result

def primer_pairing(primer_list, reference_sequence, min_amplicon_length, max_amplicon_length, Tm = 47):
    """此函数用于筛选过滤后的引物配对，并对配对后的产物长度进行筛选，并检查配对引物相互作用，并输出引物组及扩增产物序列。"""
    def primer_position_index(primers_ID, primers_list, reference_sequence, primer_direction):
        primer_ID = []
        primer = []
        primer_position = []
       
        for ID, primer_seq in zip(primers_ID, primers_list):
            primer_ID.append(str(ID))
            primer.append(str(primer_seq))
            if str(primer_direction) == str("forward"): # 正向引物位置索引
                primer_position.append(int(reference_sequence.find(primer_seq)))
            if str(primer_direction) == str("reverse"): # 反向引物位置索引
                primer_seq_RC = str(Seq(primer_seq, IUPAC.ambiguous_dna).reverse_complement())
                primer_position.append(int(reference_sequence.find(primer_seq_RC) + len(primer_seq_RC)))

        primer_position_index_result = [primer_ID, primer, primer_position]

        return primer_position_index_result

    primer_ID = primer_list[0][0]
    primer = primer_list[0][1]
    ref_seq = str(reference_sequence['sequence'])

    primer_position_index_result = \
        primer_position_index(primer_ID, primer, ref_seq, "forward") # 正向引物位置索引

    primer_reverse_complement_ID = primer_list[1][0]
    primer_reverse_complement = primer_list[1][1]

    primer_reverse_complement_position_index_result = \
        primer_position_index(primer_reverse_complement_ID, primer_reverse_complement, ref_seq, "reverse") # 反向引物位置索引

    def primer_pairing_inner(primer_position_index_result, primer_reverse_complement_position_index_result,
                             min_amplicon_length, max_amplicon_length, reference_sequence, Tm = 47):
        primer_ID = []
        primer = []
        primer_reverse_complement_ID = []
        primer_reverse_complement = []
        amplicons_length = []
        amplicon = []

        for ID, primer_seq, position in zip(
                primer_position_index_result[0],
                primer_position_index_result[1],
                primer_position_index_result[2]
        ): # 正向引物信息提取
            for reverse_complement_ID, reverse_complement_primer_seq, reverse_complement_position in zip(
                primer_reverse_complement_position_index_result[0],
                primer_reverse_complement_position_index_result[1],
                primer_reverse_complement_position_index_result[2]
            ): # 反向引物信息提取
                amplicon_length = int(reverse_complement_position) - int(position) # 产物长度计算
                if min_amplicon_length <= amplicon_length <= max_amplicon_length:
                    # 在产物长度满足限制条件时，检查正反向引物相互作用
                    heterodimer = primer3.calcHeterodimer(
                        primer_seq, reverse_complement_primer_seq, mv_conc=50.0, dv_conc=3, dna_conc=200, temp_c=25)
                    heterodimer_Tm = str(heterodimer).split('tm=')[1].split(',')[0]
                    if float(heterodimer_Tm) < float(Tm):
                        primer_ID.append(str(ID))
                        primer.append(str(primer_seq))
                        primer_reverse_complement_ID.append(str(reverse_complement_ID))
                        primer_reverse_complement.append(str(reverse_complement_primer_seq))
                        amplicons_length.append(str(amplicon_length))
                        amplicon.append(reference_sequence[position:reverse_complement_position])

        primer_pairing_inner_result = [primer_ID, primer,
                                       primer_reverse_complement_ID, primer_reverse_complement,
                                       amplicon, amplicons_length]

        return primer_pairing_inner_result


    primer_pairing_result = primer_pairing_inner(primer_position_index_result,
                                                 primer_reverse_complement_position_index_result,
                                                 min_amplicon_length, max_amplicon_length,
                                                 ref_seq, Tm)

    return primer_pairing_result


#-----------------------------------------------------------------------------------------------------------------------

file_test = sys.argv[1]

j = reference_sequence_select(file_test)
k = pragmenter(j)
l = pragmenters_Tm_calculation(k, 59, 61)
e = primer_last_5_base_calculation(l)
o = primer_poly_X_calculation(e, 3)
q = primer_homodimer_check(o)
w = primer_homodimer_check_by_primer3(q)
e = primer_hairpin_check_by_primer3(w)
r = primer_pairing(e, j, 190, 200)

for i in range(len(r[0])):
    print(str(r[0][i]) + "\t" + str(r[1][i]) + "\t" + str(r[2][i]) + "\t" +
          str(r[3][i]) + "\t" + str(r[4][i]) + "\t" + str(r[5][i]))



