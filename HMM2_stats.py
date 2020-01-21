import emission_motif
import our_gen
import transition_motifs
import markov_model_runner_guy
import pandas as pd

STATES = ["a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s"
          ,"t","u","v","w","x","y","z","A","B","C","D","H","I","J","K","L","M","N","O","P",
          "Q","R","S","T","U","V","W","X","Y","Z","1","2","3","4","5","6","7","8","9",
          "{","}","_","-","[","]","(",")"]
# transition_motifs.create_trans(15,1,15,15,1,15)

# for i in range(1, 16):
k=2
l=4
j=15
i=8
transition_motifs.create_trans(k, 1, l, j, 1, i)
emission_motif.create_emiss(k, 1, l, j, 1, i)


i_start = 6
i_end = 16
for i in range(i_start, i_end):
    emission_df = pd.read_csv("emi.tsv", sep=' ', index_col=0)
    transition_df = pd.read_csv("trans.tsv", sep=' ', index_col=0)
    pos_vals = STATES[:2 * i + 1] + ['M']
    neg_vals = STATES[2*i+1:4*i+2] + ['I','E']
    gen = our_gen.getSeq("hg19.genes.NR.chr19.exonCount2_29.bed", "hg19.2bit")
    loss = markov_model_runner_guy.run_and_stat(emission_df=emission_df, transition_df=transition_df, start_state=STATES[2*i+1],
                                                generator=gen, pos_vals=pos_vals, neg_vals=neg_vals, max_amnt=100)
    print(markov_model_runner_guy.total_acc)



