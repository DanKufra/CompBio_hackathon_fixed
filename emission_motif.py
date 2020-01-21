import numpy as np
import csv

ABC = ['A','C','G','T','^']
M_E = [0.21692973973638763, 0.30088008759486623, 0.2920972724433625, 0.1900929002253837]
I_E = [0.2372296613851884, 0.2456412115178716, 0.25376575262676015, 0.2633633744701799]

with open("e_exon_start.txt",'r') as reader:
    ms = reader.readlines()
for i, line in enumerate(ms):
    line = line[:-1]
    ms[i] = line.split(" ")

with open("e_exon_end.txt",'r') as reader:
    me = reader.readlines()
for i, line in enumerate(me):
    line = line[:-1]
    me[i] = line.split(" ")

with open("e_intron_start.txt",'r') as reader:
    ist = reader.readlines()
for i, line in enumerate(ist):
    line = line[:-1]
    ist[i] = line.split(" ")

with open("e_intron_end.txt",'r') as reader:
    ie = reader.readlines()
for i, line in enumerate(ie):
    line = line[:-1]
    ie[i] = line.split(" ")

def create_emiss(start_m, middle_m, end_m, start_i, middle_i, end_i, abc = ABC):
    states_num = start_m+middle_m+end_m+start_i+middle_i+end_i+1
    my_e = np.zeros((states_num,len(abc)))
    states = [0]*states_num
    states_so_far = 0
    for i in range(states_so_far, states_so_far+start_m):
        my_e[i][:-1] = np.array(ms[i-states_so_far])
        states[i] = 'Ms{num}'.format(num = i-states_so_far)
    states_so_far+=start_m
    for i in range(states_so_far,states_so_far+middle_m):
        my_e[i][:-1] = np.array(M_E)
        states[i] = 'Mm{num}'.format(num=i - states_so_far)
    states_so_far += middle_m
    me_index = end_m - 1
    for i in range(states_so_far,states_so_far+end_m):

        update = me[me_index]
        my_e[i][:-1] = np.array(update)
        me_index-=1
        states[i] = 'Me{num}'.format(num=i - states_so_far)
    states_so_far += end_m


    for i in range(states_so_far, states_so_far+start_i):
        my_e[i][:-1] = np.array(ist[i-states_so_far])
        states[i] = 'Is{num}'.format(num = i-states_so_far)
    states_so_far+=start_i
    for i in range(states_so_far,states_so_far+middle_i):
        my_e[i][:-1] = np.array(I_E)
        states[i] = 'Im{num}'.format(num=i - states_so_far)
    states_so_far += middle_i
    ie_index = end_i - 1
    for i in range(states_so_far,states_so_far+end_i):

        update = ie[ie_index]
        my_e[i][:-1] = np.array(update)
        ie_index-=1
        states[i] = 'Ie{num}'.format(num=i - states_so_far)
    states_so_far += end_i

    states[-1] = "end"
    my_e[-1,-1] = 1.


    with open('emi.tsv', 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter=' ')
        tsv_writer.writerow(['index']+list(abc))
        for i in range(states_num):
            tsv_writer.writerow([states[i]] + list(my_e[i]))








