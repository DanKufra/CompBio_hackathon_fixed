import numpy as np
import csv

AVG_NUM_M = 6
AVG_LEN_M = 160
SUM_LEN_M = 160*6
AVG_LEN_I = 18459.46



def create_trans(start_m, middle_m, end_m, start_i, middle_i, end_i):
    p_m_m = 1 / (AVG_LEN_M - start_m-end_m)
    p_m_m = p_m_m **(1/middle_m)
    p_i_m = AVG_NUM_M/(AVG_LEN_I-((AVG_NUM_M+1)*(end_i+start_i)))
    p_i_m = p_i_m ** (1/middle_i)
    p_i_i = 1-p_i_m
    states_num = start_m+middle_m+end_m+start_i+middle_i+end_i+1
    states = [0] * states_num
    my_t = np.zeros((states_num,states_num))
    states_so_far = 0
    my_t[-1,-1]=1
    for i in range(start_m):
        my_t[i,i+1] = 1
        states[i] = 'Ms{num}'.format(num = i-states_so_far)
    states_so_far +=start_m
    for i in range(states_so_far, states_so_far + middle_m):
        my_t[i,i+1] = 1-p_m_m
        my_t[i,i] = p_m_m
        states[i] = 'Mm{num}'.format(num=i - states_so_far)
    states_so_far+=middle_m
    for i in range(states_so_far, states_so_far+end_m):
        my_t[i, i + 1] = 1
        states[i] = 'Me{num}'.format(num=i - states_so_far)
    states_so_far += end_m


    for i in range(states_so_far, states_so_far+start_i):
        my_t[i,i+1] = 1
        states[i] = 'Is{num}'.format(num=i - states_so_far)
    states_so_far +=start_i
    for i in range(states_so_far, states_so_far + middle_i):
        my_t[i,i+1] = 1-p_i_i
        my_t[i,i] = p_i_i
        states[i] = 'Im{num}'.format(num=i - states_so_far)
    states_so_far+=middle_i
    for i in range(states_so_far, states_so_far+end_i):
        my_t[i, i + 1] = 1
        states[i] = 'Ie{num}'.format(num=i - states_so_far)
    states_so_far += end_i

    my_t[-2,0] =1-1/6
    my_t[-2,-1]= 1/6
    states[-1]='end'


    with open('trans.tsv', 'wt') as out_file:
        tsv_writer = csv.writer(out_file, delimiter=' ')
        tsv_writer.writerow(['index']+list(states))
        for i in range(states_num):
            tsv_writer.writerow([states[i]] + list(my_t[i]))

