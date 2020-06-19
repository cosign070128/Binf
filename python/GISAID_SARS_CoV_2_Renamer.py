#此程序主要用于GISAID SARS_CoV_2序列ID调整。
#GISAID命名规则不遵循NCBI命名规则，在数据库构建过程中会出现序列ID过长的情况，导致建库失败。
#此程序可以将GISAID命名规格调整为NCBI命名规则，且速度明显快于linux-sed功能。
#举例：
#GISAID:    >hCoV-19/USA/WA-UW-1729/2020|EPI_ISL_424220|2020-03-21
#NCBI:      >EPI_ISL_424220 hCoV-19/USA/WA-UW-1729/2020|EPI_ISL_424220|2020-03-21

from Bio import SeqIO
import sys

if len(sys.argv) != 3:  ##参数数量不等于3的情况下，提示使用信息
    print("程序用法:python GISAID_SARS_CoV_2_Renamer.py GISAID_Download_file File_after_renamed")
    print("此程序用于GISAID SARS_CoV_2序列ID调整成符合NCBI建库规则的序列ID。")
    print("转化前：>hCoV-19/USA/WA-UW-1729/2020|EPI_ISL_424220|2020-03-21")
    print("转化后：>EPI_ISL_424220 hCoV-19/USA/WA-UW-1729/2020|EPI_ISL_424220|2020-03-21")
    print("作者：张迪骏  日期：2020.6.19")
    exit()

GISAID_download_file = sys.argv[1]
file_after_renamed = sys.argv[2]
file_after_renamed_write = open(file_after_renamed, "w")

for seq_record in SeqIO.parse(GISAID_download_file, "fasta"):
#使用异常模块的原因：
#   在GISAID命名规则中，会加入国家名称，有一些国家的名称中包含空格，例如
#   >hCoV-19/Bosnia and Herzegovina/03_Tuzla/2020|EPI_ISL_463893|2020-05-25
#   这会导致SeqIO模块将 hCoV-19/Bosnia 作为seq.id;
#   and Herzegovina/03_Tuzla/2020|EPI_ISL_463893|2020-05-25作为seq.description；
#   若不存在空格，则EPI_ISL_******包含在id中，反之，则包含在description中。
    try:
        seq_record_new_id = str(seq_record.id).split("|")[1]
        seq_record.description = seq_record.id
    except:
        seq_record_new_id = str(seq_record.description).split("|")[1]

    seq_record.id = seq_record_new_id
    file_after_renamed_write.write(">" + seq_record.id + " " + seq_record.description + "\n")
    file_after_renamed_write.write(str(seq_record.seq) + "\n")

file_after_renamed_write.close()
