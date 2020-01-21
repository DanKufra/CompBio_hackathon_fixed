import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
from Bio.HMM.MarkovModel import MarkovModelBuilder, HiddenMarkovModel
from Bio import Alphabet
from Bio.Seq import MutableSeq

if __name__ == '__main__':
    emission_stats_df = pd.read_csv('./emission_stats_df.tsv', sep='\t')
    new_df = pd.DataFrame(columns=['name', 'state', 'count', 'letter'])
    for l in ['A', 'C', 'G', 'T']:
        tmp_df = emission_stats_df[['name', 'state', l]]
        tmp_df.insert(2, 'letter', l)
        tmp_df = tmp_df.rename(columns={'name': 'name', 'state':'state', l: 'count', 'letter': 'letter'})
        new_df = new_df.append(tmp_df)

    # g = sns.FacetGrid(new_df, col="state", hue='letter')
    # g = g.map(plt.hist, "count", histtype='barstacked', bins=np.arange(0, 4000, 50))
    # g.add_legend()
    # plt.show()

    # Intron emissions
    intron_emissions = emission_stats_df[emission_stats_df['state'] == 'I'][['A', 'C', 'G', 'T']]
    intron_emissions = np.sum(intron_emissions) / np.sum(intron_emissions.values)
    # First extron emissions
    first_exon_emissions = emission_stats_df[emission_stats_df['state'] == 'F'][['A', 'C', 'G', 'T']]
    first_exon_emissions = np.sum(first_exon_emissions) / np.sum(first_exon_emissions.values)
    # Middle extron emissions
    middle_exon_emissions = emission_stats_df[emission_stats_df['state'] == 'M'][['A', 'C', 'G', 'T']]
    middle_exon_emissions = np.sum(middle_exon_emissions) / np.sum(middle_exon_emissions.values)
    # First extron emissions
    last_exon_emissions = emission_stats_df[emission_stats_df['state'] == 'L'][['A', 'C', 'G', 'T']]
    last_exon_emissions = np.sum(last_exon_emissions) / np.sum(last_exon_emissions.values)

    emission_matrix = np.stack([first_exon_emissions, middle_exon_emissions, last_exon_emissions, intron_emissions])
    emission_matrix = pd.DataFrame(emission_matrix, columns=['A', 'C', 'G', 'T'])
    emission_matrix.insert(4, '^', [0.0, 0.0, 0.0, 0.0])
    emission_matrix = emission_matrix.append(pd.Series({'A': 0.0, 'C': 0.0, 'G': 0.0, 'T': 0.0, '^': 1.0}), ignore_index=True)
    emission_matrix.insert(0, 'index', ['F', 'M', 'L', 'I', 'E'])
    emission_matrix.to_csv('./initial_emissions_matrix.tsv', sep=' ', index=False)
