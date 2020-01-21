import twobitreader
import argparse
import pandas as pd
import numpy as np
import collections
import matplotlib.pyplot as plt


GENE_PATH = "hg19.genes.NR.chr19.exonCount2_29.bed"
TWO_BIT_PATH = "hg19.2bit"

ABC = {"AA":'a', "AC":'b', "AG":'c', "AT":'d',
       "CA":'e', "CC":'f', "CG":'g', "CT":'h',
       "GA":'i', "GC":'j',"GG":'k',"GT":'l',
       "TA":'m', "TC":'n', "TG":'o', "TT":'p'}
# GENE_TSV_HEADER = ['name', 'chrom', 'strand', 'txStart', 'txEnd', 'cdsStart', 'cdsEnd', 'exonCount', 'exonStarts', 'exonEnds', 'proteinId', 'alignID']
GENE_TSV_HEADER = ['chrom', 'startTranscription', 'endTranscription', 'name', 'unimportant', 'strand', 'startTranslation', 'endTranslation', 'unimportant2', 'exonCount', 'exonSize', 'exonStart']

REVERSE_COMPLEMENT_MAP = {"A": "T", "T": "A", "C": "G", "G": "C"}


def getSeq(gene_tsv_path, twobitpath):
    genes_df = pd.read_csv(gene_tsv_path, sep='\t', names=GENE_TSV_HEADER)
    genome_reader = twobitreader.TwoBitFile(twobitpath)
    for i, gene in genes_df.iterrows():
        seq = genome_reader[gene['chrom']][int(gene['startTranscription']):int(gene['endTranscription'])]
        seq = seq.upper()
        seq = np.array(list(seq))

        if gene['strand'] == '-':
            seq = np.flip(seq)
            seq = np.vectorize(REVERSE_COMPLEMENT_MAP.get)(seq)
        exon_starts = np.array(gene['exonStart'].split(',')[:-1], dtype=np.int64)
        exon_sizes = np.array(gene['exonSize'].split(',')[:-1], dtype=np.int64)
        intron_exon_lbl = np.full(len(seq), "I")
        for exon_idx, (start, size) in enumerate(zip(exon_starts, exon_sizes)):
            if exon_idx == 0:
                exon_lbl = 'F'
            elif exon_idx == len(exon_starts) - 1:
                exon_lbl = 'L'
            else:
                exon_lbl = 'M'
            intron_exon_lbl[start: start + size] = exon_lbl
        first_intron_idx = np.where(intron_exon_lbl == 'I')[0][0]
        first_L_idx = np.where(intron_exon_lbl == 'L')[0][0]
        seq = seq[first_intron_idx:first_L_idx]
        label = intron_exon_lbl[first_intron_idx:first_L_idx]
        seq_ret = np.array([seq[0]]+[ABC[seq[i]+seq[i+1]] for i in range(len(seq)-1)])


        yield seq_ret[1:], label[1:]


# def parse_args():
#     parser = argparse.ArgumentParser()
#     parser.add_argument('twobit_path', help='2bit file we use to get sequences)')
#     parser.add_argument('gene_tsv_path', help='gene tsv path')
#     return parser.parse_args()


b=next(getSeq(GENE_PATH, TWO_BIT_PATH))
a=2