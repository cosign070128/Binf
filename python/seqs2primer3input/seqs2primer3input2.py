#!/usr/bin/python3
#coding=utf-8

#*********************************************************************************************************************
#                                                                                                                     
#  1.此程序主要用于primer3输入文件准备，可以将序列文件转化成primer3格式的输入文件。                                     
#                                                                                                                   
#  2.较1.0版本主要改进内容：从参数内嵌式(每次使用软件需要修改源代码)改成参数外挂式(每次使用不涉及源代码修改)。       
#                                                                                                                    
#  3.软件使用方式：                                                                                                  
#      python seqs2primer3input.py SEQUENCE.fasta PARAMETERS.primer3input TEMPLATE_NEW_ID PRIMER3_INPUT_FILE         
#        SEQUENCE.fasta             包含所有需要设计序列的原始文件                                                   
#        PARAMETERS.primer3input    包含此次引物设计参数的文件                                                       
#        TEMPLATE_NEW_ID            为明确设计模版名称，需自定义模版ID，输出文件会按照ID_1，ID_2...形式对模版重新命名
#        PRIMER3_INPUT_FILE         primer3_core输入文件                                                             
#                                                                                                                    
#  4.软件信息：                                                                                                      
#      作者：张迪骏  编写日期：2020.09.24   版本号：2.0                                                              
#                                                                                                                    
#*********************************************************************************************************************

from Bio import SeqIO
import sys

#通过读取PARAMETERS.primer3input文件获取参数，利用字典键和值记录参数，并最终返回字典
def parameter_dict(parameterFile):
    parameter_file_open = open(parameterFile, 'r')
    parameter_dict_extracter = {}  #构建空字典
    for line in parameter_file_open:
        key = str(line.split('=')[0].strip())  #获取键
        value = str(line.split('=')[1].strip())   #获取值
        parameter_dict_extracter[key] = value   #为字典添加键：值
    parameter_file_open.close()
    return parameter_dict_extracter   #返回字典

if len(sys.argv) != 5:
    print("python seqs2primer3input.py SEQUENCE.fasta PARAMETERS.primer3input TEMPLATE_NEW_ID PRIMER3_INPUT_FILE")
    sys.exit()

sequenceFile = sys.argv[1]
parametersPrimer3Input = sys.argv[2]
templateNewID = sys.argv[3]
primer3InputFile = sys.argv[4]
primer3InputFileOpen = open(primer3InputFile, 'w')
parameters = parameter_dict(parametersPrimer3Input)   #用parameters接走parameter_dict返回的字典

i = 1
for sequenceRecord in SeqIO.parse(sequenceFile, "fasta"):
    sequenceRecord.id = "SEQUENCE_ID=" + templateNewID + "_" + str(i)  #修改模版ID
    primer3InputFileOpen.write(sequenceRecord.id + "\n")
    primer3InputFileOpen.write("SEQUENCE_TEMPLATE=" + str(sequenceRecord.seq) + "\n")
    for parameters_key, parameters_value in parameters.items():
        primer3InputFileOpen.write(parameters_key + '=' + parameters_value + "\n")
    primer3InputFileOpen.write("=" + "\n")  #输出字典中的参数内容
    i += 1
primer3InputFileOpen.close()
