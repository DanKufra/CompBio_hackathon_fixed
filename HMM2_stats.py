import emission_motif
import our_gen
import transition_motifs
import markov_model_runner
import pandas as pd

# transition_motifs.create_trans(15,1,15,15,1,15)

# for i in range(1, 16):
    #transition_motifs.create_trans(i, 1, i, i, 1, i)
    #emission_motif.create_emiss(i, 1, i, i, 1, i)

    emission_df = pd.read_csv("emi.tsv", sep=' ', index_col=0)
    transition_df = pd.read_csv("trans.tsv", sep=' ', index_col=0)

    gen = our_gen.getSeq("hg19.genes.NR.chr19.exonCount2_29.bed", "hg19.2bit")
    loss = markov_model_runner.run_and_stat(emission_df=emission_df, transition_df=transition_df, start_state="d",
                                            generator=gen,pos_vals=['a','b','c'], neg_vals=['d','e','f'])



