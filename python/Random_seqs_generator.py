#coding=utf-8
#此程序主要用于生成一定数量的随机核酸序列。
#主要参数如下：
#1.需要生成的随机序列数量；
#2.需要生成的随机序列最小长度；
#3.需要生成的随机序列最大长度。

import sys
import random

#定义随机序列生成器函数
def random_seqs_generator(min_seq_length, max_seq_length):
    base_dict = {'1': 'A', '2': 'T', '3': 'G', '4': 'C'}  #定义ATGC字典，键设定为1-4，方便随机索引，值设定为ATGC
    seq = []  #定义一个空列表，用于存在序列碱基
    seq_length = random.randint(min_seq_length, max_seq_length)  #在最小序列长度和最大序列长度之间随机选取一个序列长度
    for i in range(seq_length):
        seq.append((base_dict[str(random.randint(1, 4))]))  #在指定序列长度下，随机生成相同数量的随机碱基并添加到列表中
    random_seq = ''.join(map(str, seq))  #将列表转化为字符串
#    print(random_seq)
    return random_seq

if len(sys.argv) != 5:
    print("Usage: python Random_seqs_generator.py number_of_random_seqs min_seq_length max_seq_length output_file")
    print("This program is used to generate random sequences.")
    print("Author: Dijun Zhang   Date: 2020.08.19")

seqs_number = sys.argv[1]
seq_length_min = sys.argv[2]
seq_length_max = sys.argv[3]
result_file = sys.argv[4]

f = open(result_file, "w")

for number in range(int(seqs_number)):
    f.write(str(">" + "random_seq_" + str(number+1) + "\n"))
    f.write(str(random_seqs_generator(int(seq_length_min), int(seq_length_max))) + "\n")
f.close()
