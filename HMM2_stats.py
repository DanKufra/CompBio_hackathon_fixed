import emission_motif
import our_gen
import transition_motifs
import markov_model_runner

#transition_motifs.create_trans(15,1,15,15,1,15)

for i in range(1,16):
    print(i)
    transition_motifs.create_trans(i,1,i,i,1,i)
    emission_motif.create_emiss(i,1,i,i,1,i)
    gen = our_gen.getSeq("hg19.genes.NR.chr19.exonCount2_29.bed", "hg19.2bit")
    loss = markov_model_runner("trans.tsv","emi.tsv",'')