import pandas as pd
import numpy as np
import argparse
import seaborn as sns
import matplotlib.pyplot as plt
from Bio.HMM.MarkovModel import MarkovModelBuilder
from Bio import Alphabet
from Bio.Seq import MutableSeq
from seq_generator import getSeq
from tqdm import tqdm
from sklearn.metrics import confusion_matrix
from confusion_matrix_pretty_print import plot_confusion_matrix_from_data


def build_mm(emission_df, transition_df, start_state):
    state_alphabet = Alphabet.Alphabet()
    state_alphabet.letters = transition_df.columns
    emission_alphabet = Alphabet.Alphabet()
    emission_alphabet.letters = emission_df.columns
    m = MarkovModelBuilder(state_alphabet, emission_alphabet)
    m.set_initial_probabilities({start_state: 1.0})
    for state in state_alphabet.letters:
        for state2 in state_alphabet.letters:
            prob = transition_df.loc[state][state2]
            m.allow_transition(state, state2, prob)
    for state in state_alphabet.letters:
        for letter in emission_alphabet.letters:
                prob = emission_df.loc[state][letter]
                m.set_emission_score(state, letter, prob)

    mm = m.get_markov_model()
    return mm


def run_viterbi(markov_model, seq, emission_array, state_array, blocks=False):
    state_alphabet = Alphabet.Alphabet()
    state_alphabet.letters = state_array
    emission_alphabet = Alphabet.Alphabet()
    emission_alphabet.letters = emission_array

    mutable_seq = MutableSeq("", emission_alphabet)
    # generate the sequence
    last_best = None
    for i, letter in enumerate(seq):
        # add on a new roll to the sequence
        if blocks:
            try:
                if len(letter) == 0:
                    mutable_seq.append(last_best)
                else:
                    mutable_seq.append(letter[0])
                    last_best = letter[0]
            except TypeError:
                mutable_seq.append(last_best)
        else:
            if letter == 'N':
                continue
            mutable_seq.append(letter[0].upper())
    mutable_seq = mutable_seq.toseq()

    predicted_states, prob = markov_model.viterbi(mutable_seq, state_alphabet)
    return predicted_states, prob


def calc_statistics(transition_df, pred, lbl, pos_vals=None, neg_vals=None, i=None, verbose=False):
    # get overall accuracy
    overall_accuracy = np.sum(pred == lbl) / len(pred)

    # TPR
    if pos_vals  is not None:
        pred_is_pos = np.zeros(len(pred))
        for val in pos_vals:
            pred_is_pos[pred == val] = 1.0
        lbl_is_pos = np.zeros(len(lbl))
        for val in pos_vals:
            lbl_is_pos[lbl == val] = 1.0
        TPR = np.sum((pred_is_pos == lbl_is_pos) & (lbl_is_pos == 1.0)) / np.sum(lbl_is_pos)

    # TNR
    if neg_vals  is not None:
        pred_is_neg = np.zeros(len(pred))
        for val in neg_vals:
            pred_is_neg[pred == val] = 1.0
        lbl_is_neg = np.zeros(len(lbl))
        for val in neg_vals:
            lbl_is_neg[lbl == val] = 1.0
        TNR = np.sum((pred_is_neg == lbl_is_neg) & (lbl_is_neg == 1.0)) / np.sum(lbl_is_neg)

    if verbose:
        print("ACC: %f TPR: %f TNR: %f iter: %d" % (overall_accuracy, TPR, TNR, i))
    # get confusion matrix of preds vs labels (Multiclass)
    confusion = confusion_matrix(lbl, pred, labels=transition_df.columns)
    return overall_accuracy, TPR, TNR, confusion

def map_block_to_emission(block, blocks):
    new_emission = ''
    block = np.array(block)
    for l in ['A', 'C', 'G', 'T']:
        new_emission += str(np.sum(block == l))
    return new_emission, np.where(blocks == new_emission)[0]


def create_all_blocks():
    # dummy_strings = ['0005', '0014', '0023', '0122', '0113', '1112']
    dummy_strings = ['0003', '0012', '0111']
    from itertools import permutations
    blocks = []
    for string in dummy_strings:
        blocks.extend([''.join(p) for p in permutations(string)])
    blocks = np.unique(blocks)
    unicode_blocks = np.array([chr(i) for i in np.arange(len(blocks))])
    return blocks, unicode_blocks


def run_on_blocks():
    blocks, unicode_blocks = create_all_blocks()
    # TODO create real emission matrix
    emission_dict = {}
    for b in unicode_blocks:
        emission_dict[b] = [1.0 / len(blocks)]
    emission_df = pd.DataFrame(emission_dict, columns=unicode_blocks)
    emission_dict = {}
    for b in unicode_blocks:
        emission_dict[b] = 1.0 / len(blocks)
    emission_df = emission_df.append(pd.Series(emission_dict), ignore_index=True)
    emission_df.insert(0, 'index', ['e', 'i'])
    emission_df.set_index('index', inplace=True)

    transition_df = pd.DataFrame([{'e': 0.01, 'i': 0.99},
                                  {'e': 0.99, 'i': 0.01}], columns=['e', 'i'])
    transition_df.insert(0, 'index', ['e', 'i'])
    transition_df.set_index('index', inplace=True)
    mm = build_mm(emission_df, transition_df, 'e')
    gene_gen = getSeq('./hg19.genes.NR.chr19.exonCount2_29.bed', './hg19.2bit')
    for i, (seq, lbl) in enumerate(gene_gen):
        seq_mapped = []
        for block in range(0, len(seq), 3):
            cur_block, cur_unicode_block_idx = map_block_to_emission(seq[block: block + 3].tolist(), np.array(blocks))
            seq_mapped.append(unicode_blocks[cur_unicode_block_idx])
        predicted_states, prob  = run_viterbi(mm, seq_mapped, emission_df.columns, transition_df.columns, blocks=True)
        print(predicted_states)


def run_and_stat(emission_df, transition_df, start_state, generator, pos_vals=None, neg_vals=None, max_amnt=None):
    mm = build_mm(emission_df, transition_df, start_state)

    acc_arr = []
    tpr_arr = []
    tnr_arr = []
    total_confusion = np.zeros((len(transition_df.columns), len(transition_df.columns)))
    for i, (seq, lbl) in enumerate(generator):
        seq = seq.tolist()
        seq.append('^')
        predicted_states, prob = run_viterbi(mm, seq, emission_df.columns, transition_df.columns)
        # print(predicted_states)
        lbl_str = ''
        for l, c in zip(lbl, seq):
            lbl_str += l
        lbl_str += 'E'
        # print(lbl_str)
        intron_only = ''
        for j in range(len(lbl_str)):
            intron_only += 'I'

        intron_only = np.array(list(intron_only))
        lbl_arr = np.array(list(lbl_str))
        pred_arr = np.array(list(predicted_states))
        acc, tpr, tnr, confusion = calc_statistics(transition_df, pred_arr, lbl_arr,
                                                   pos_vals=pos_vals, neg_vals=neg_vals,
                                                   i=i, verbose=True)
        acc_arr.append(acc)
        tpr_arr.append(tpr)
        tnr_arr.append(tnr)
        total_confusion += confusion
        if max_amnt is not None and i > max_amnt:
            break

    acc_arr = np.array(acc_arr)
    tpr_arr = np.array(tpr_arr)
    tnr_arr = np.array(tnr_arr)
    print(acc_arr.mean(), tpr_arr.mean(), tnr_arr.mean())
    # print(total_confusion[0])
    # print(np.sum(total_confusion, axis=1)[0])
    # print(np.sum(total_confusion, axis=0)[0])
    # # print(np.round(total_confusion / np.sum(total_confusion, axis=1), decimals=4))
    # normalized_confusion = total_confusion / np.sum(total_confusion, axis=1)
    # print(normalized_confusion)
    # # y is label
    # # x is pred
    # confusion_matrix = pd.DataFrame(normalized_confusion, columns=transition_df.columns)
    # ax = sns.heatmap(confusion_matrix, annot=True, cbar=False)
    # plt.show()
    plot_confusion_matrix_from_data(lbl_arr, pred_arr, transition_df.columns, True,
                                    'Oranges', '.3f', 12, 0.5, False, [9, 9], 2, 'y', confusion=total_confusion)


if __name__ == '__main__':
    # run_on_blocks()
    # args = parse_args()
    emission_df = pd.read_csv('./initial_emissions_matrix.tsv', sep=' ', index_col=0)
    transition_df = pd.read_csv('./initial_transitions_matrix.tsv', sep=' ', index_col=0)


    gene_gen = getSeq('./hg19.genes.NR.chr19.exonCount2_29.bed', './hg19.2bit')
    run_and_stat(emission_df, transition_df, 'F', gene_gen,  ['F', 'M', 'L'], ['I'], max_amnt=150)