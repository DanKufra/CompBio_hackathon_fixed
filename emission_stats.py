import pandas as pd
import numpy as np
import seaborn
import argparse
from seq_generator import getSeq, GENE_TSV_HEADER, REVERSE_COMPLEMENT_MAP
import twobitreader
from tqdm import tqdm


def get_next_exon_intron(gene_tsv_path, twobitpath, intron=False):
    genes_df = pd.read_csv(gene_tsv_path, sep='\t', names=GENE_TSV_HEADER)
    genome_reader = twobitreader.TwoBitFile(twobitpath)
    for i, gene in genes_df.iterrows():
        exon_starts = np.array(gene['exonStart'].split(',')[:-1], dtype=np.int64) + int(gene['startTranscription'])
        exon_sizes = np.array(gene['exonSize'].split(',')[:-1], dtype=np.int64)
        if intron:
            starts = [start + size for start, size in zip(exon_starts[:-1], exon_sizes[:-1])]
            sizes = [exon_start - intron_start for exon_start, intron_start in zip(exon_starts[1:], starts[:-1])]
            if len(starts) == 1:
                sizes = [exon_starts[-1] - starts[0]]
        else:
            starts = exon_starts
            sizes = exon_sizes
        for idx, (start, size) in enumerate(zip(starts, sizes)):
            if size <= 0:
                continue
            seq = np.array(list(genome_reader[gene['chrom']][start:start + size].upper()))
            if gene['strand'] == '-':
                seq = np.flip(seq)
                seq = np.vectorize(REVERSE_COMPLEMENT_MAP.get)(seq)
            if intron:
                state = 'I'
            elif idx == 0:
                state = 'F'
            elif idx == len(exon_starts) - 1:
                state = 'L'
            else:
                state = 'M'
            yield seq, state, gene['name']


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument('twobit_path', help='2bit file we use to get sequences)')
    parser.add_argument('gene_tsv_path', help='gene tsv path')
    return parser.parse_args()


if __name__ == '__main__':
    args = parse_args()
    # usage example
    exon_gen = get_next_exon_intron(args.gene_tsv_path, args.twobit_path, intron=False)
    intron_gen = get_next_exon_intron(args.gene_tsv_path, args.twobit_path, intron=True)
    emission_stats_df = pd.DataFrame(columns=['name', 'state', 'A', 'C', 'G', 'T'])

    for i, (seq, lbl, gene_name) in tqdm(enumerate(intron_gen)):
        emission_stats_df = emission_stats_df.append(pd.Series({"name": gene_name, "state": lbl,
                                                                "A": np.sum(seq == 'A'),
                                                                "C": np.sum(seq == 'C'),
                                                                "G": np.sum(seq == 'G'),
                                                                "T": np.sum(seq == 'T')}), ignore_index=True)


    for i, (seq, lbl, gene_name) in tqdm(enumerate(exon_gen)):
        emission_stats_df = emission_stats_df.append(pd.Series({"name" : gene_name, "state" : lbl,
                                                                "A": np.sum(seq == 'A'),
                                                                "C": np.sum(seq == 'C'),
                                                                "G": np.sum(seq == 'G'),
                                                                "T": np.sum(seq == 'T')}), ignore_index=True)



    emission_stats_df.to_csv("./emission_stats_df.tsv", sep=" ")
    print('end')