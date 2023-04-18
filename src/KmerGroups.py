import pandas as pd
import numpy as np
import textdistance as td
import sklearn.cluster
from skbio import RNA
from skbio.alignment import global_pairwise_align_nucleotide
from itertools import combinations
import warnings
import seqlogo
import logomaker
import matplotlib.pyplot as plt
import seaborn as sns
warnings.filterwarnings("ignore", message="You're using skbio's python implementation of Needleman-Wunsch alignment.")

########## Get kmer clusters ##########
def get_unique_string_tokens(test_str):
    tokens = [test_str[i: j] for i in range(len(test_str)) for j in range(i + 1, len(test_str) + 1)]
    # Filter tokens that correspond to test_str
    tokens = [t for t in tokens if t != test_str]
    ugca_tokens = []
    unique_tokens = sorted(list(set(tokens)))
    token_counts = [tokens.count(i) for i in unique_tokens]
    for tup in list(zip(unique_tokens, token_counts)):
        for i in range(tup[1]):
            ugca_tokens.append(f'{tup[0]}{i+1}')
    return ugca_tokens

def get_unique_string_tokensNoOccurrence(test_str):
    tokens = [test_str[i: j] for i in range(len(test_str)) for j in range(i + 1, len(test_str) + 1)]
    # Filter tokens that correspond to test_str
    tokens = [t for t in tokens if t != test_str]
    unique_tokens = sorted(list(set(tokens)))
    return unique_tokens


def SortTokensByLength(tokens):
    lenDict = {}
    for length in range(np.min([len(t) for t in tokens]), np.max([len(t) for t in tokens]) + 1):
        lenDict[length-1] = [t for t in tokens if len(t) == length]
    return lenDict

def ComputePairwiseKmerSimilarity(k1, k2):
    dict1 = SortTokensByLength(get_unique_string_tokens(k1))
    dict2 = SortTokensByLength(get_unique_string_tokens(k2))
    jaccardSimilarities = []
    lengths = [l for l in dict1.keys() if l in dict2.keys()]
    for tokenLength in lengths:
        jaccardSimilarities.append(td.jaccard.similarity(dict1[tokenLength], dict2[tokenLength]))
    # Mean Jaaccard is between 0 and 1
    MeanJaccard = np.mean(jaccardSimilarities)
    return MeanJaccard

def ComputePairwiseKmerSimilarityNoOccurrence(k1, k2, minLen=1):
    dict1 = SortTokensByLength(get_unique_string_tokensNoOccurrence(k1))
    dict2 = SortTokensByLength(get_unique_string_tokensNoOccurrence(k2))
    jaccardSimilarities = []
    lengths = [l for l in dict1.keys() if l in dict2.keys()]
    for tokenLength in lengths:
        if tokenLength >= minLen:
            jaccardSimilarities.append(td.jaccard.similarity(dict1[tokenLength], dict2[tokenLength]))
    # Mean Jaaccard is between 0 and 1
    MeanJaccard = np.mean(jaccardSimilarities)
    return MeanJaccard


def GetKmerSimilarityMatrix(words):
    words = np.asarray(words)
    SimilarityMatrix = np.array(
        [[ComputePairwiseKmerSimilarity(w1, w2) for w1 in words] for w2 in words]
        )
    return SimilarityMatrix

def GetKmerSimilarityMatrixNoOccurrence(words, minLen=1):
    words = np.asarray(words)
    SimilarityMatrix = np.array(
        [[ComputePairwiseKmerSimilarityNoOccurrence(w1, w2, minLen) for w1 in words] for w2 in words]
        )
    return SimilarityMatrix


def get_unique_string_tokens_w_occurrence(test_str, minSubstrLength=1):
    tokens = [test_str[i: j] for i in range(len(test_str)) for j in range(i + 1, len(test_str) + 1)]
    tokens = [t for t in tokens if   len(t)>=minSubstrLength]
    ugca_tokens = []
    unique_tokens = sorted(list(set(tokens)))
    token_counts = [tokens.count(i) for i in unique_tokens]
    for tup in list(zip(unique_tokens, token_counts)):
        for i in range(tup[1]):
            ugca_tokens.append(f'{tup[0]}{i+1}')
    return ugca_tokens

def get_distances(words, minSubstrLength=1):
    words = np.asarray(words) #So that indexing with a list will work
    similarity_ugca = np.array([[td.jaccard.similarity(get_unique_string_tokens_w_occurrence(w1, minSubstrLength), get_unique_string_tokens_w_occurrence(w2, minSubstrLength)) for w1 in words] for w2 in words])
    return similarity_ugca

def get_similarity_weighed(words, minSubstrLength=1):
    words = np.asarray(words) #So that indexing with a list will work
    similarity_ugca = np.array([[td.jaccard.similarity(get_unique_string_tokens_weighed(w1, minSubstrLength), get_unique_string_tokens_weighed(w2, minSubstrLength)) for w1 in words] for w2 in words])
    return similarity_ugca

def get_clustering(similarity_matrix, d):
    affprop = sklearn.cluster.AffinityPropagation(affinity="precomputed", damping=d, max_iter=1000, convergence_iter=200, random_state=0)
    affprop.fit(similarity_matrix)
    return affprop.labels_, len(np.unique(affprop.labels_))

def get_clustered_df(motifs, d=0.5):
    """
    Returns a pandas dataframe with k-mers and the cluster to which kmer was assigned.
    Distance matrix is generated from jaccard index of tokenized k-mers.
    Clustering of k-mers is performed with affinity propagation algorithm.
    -------
    Inputs_
    motifs - list of kmers to cluster
    d : damping parameter for affinity propafation, defaoult 0.5 (accepts vals between 0.5 an 1)
    """
    df_cl = pd.DataFrame(columns=['motif', 'cluster'])
    df_cl['motif'] = motifs
    similarity_ugca = get_distances(motifs, minSubstrLength=1)
    labels = get_clustering(similarity_ugca, d)[0]
    df_cl['cluster'] = labels
    df_cl.sort_values(by='cluster', ascending=True, inplace=True)
    return df_cl

def GetClusteredDfFixedNClusters(motifs, n=3):
    """
    Returns a pandas dataframe with k-mers and the cluster to which kmer was assigned.
    Distance matrix is generated from jaccard index of tokenized k-mers.
    Clustering of k-mers is performed with kMeans algorithm.
    """
    distanceMatrix = 1 - get_distances(motifs, minSubstrLength=1)
    kmeans = sklearn.cluster.KMeans(n_clusters=n)
    X = kmeans.fit(distanceMatrix)
    df = pd.DataFrame([motifs, X.labels_], index=['kmers', 'cluster']).T.sort_values(by='cluster')
    return df

########## Get alignment of k-mers ##########

def get_prefix_suffix(j, motif):
    prefix = j.split(motif)[0]
    suffix = j.split(motif)[-1]
    return prefix, suffix

def get_alignments(motif_pairs_list, match_s, missm_s):
    df_cons = pd.DataFrame(columns=['consensus', 'score', 'aligned_mots', 'motif_pair'])
    for i, pair in enumerate(motif_pairs_list):
        alignment, score, start_end_positions = global_pairwise_align_nucleotide(
            RNA(pair[0]), RNA(pair[1]),
            match_score=match_s, mismatch_score=missm_s
            )
        df_cons.loc[i, :] = str(alignment.consensus()), score, tuple(str(alignment).split('\n')[-2:]), pair
        df_cons = df_cons.sort_values(by='score', ascending=False)
    return df_cons

def get_alignment(mots, match_s=2, missm_s=-1):
    motif_pairs = list(combinations(mots, r=2))
    motif_pairs = [tuple(sorted(m)) for m in motif_pairs]
    #Get pairwise alignments
    df_cons = get_alignments(motif_pairs, match_s, missm_s)
    # Get best alignment out of all as a tuple of aligned motifs
    aligned = list(df_cons.sort_values(by='score', ascending=False)['motif_pair'].values.tolist()[0])
    aligned_sequences = list(df_cons.sort_values(by='score', ascending=False)['aligned_mots'].values.tolist()[0])
    # remove aligned motifs from the list of all motifs
    remaining_mots = [m for m in mots if m not in aligned]
    while len(remaining_mots)>0:
        #Find the best alignment related to one of the aligned motifs (reference motifs)
        aligned_pairs = list(combinations(aligned, r=2))
        aligned_pairs = [tuple(sorted(m)) for m in aligned_pairs]
        lookup_pairs = [m for m in motif_pairs if m not in aligned_pairs]
        lookup_pairs = [m for m in lookup_pairs if (m[0] in aligned) or (m[1] in aligned)]
        df_select = df_cons.loc[df_cons['motif_pair'].isin(lookup_pairs)].sort_values(by='score', ascending=False)
        second_alignment = df_select['aligned_mots'].values.tolist()[0]
        second_pair = df_select['motif_pair'].values.tolist()[0]
        # Set reference motif
        y = [m for m in aligned if m in second_pair][0]
        ref1 = [m for m in aligned_sequences if m.replace('-','')==y][0]
        ref2 = [m for m in second_alignment if m.replace('-','')==y][0]
        # Define other motif to align
        seq2 = [m for m in second_alignment if m.replace('-','')!=y][0]
        # All other aligned motifs must be changed accordingly
        al_seqs = [m for m in aligned_sequences if m.replace('-', '')!=y]
        #The ref2 and seq2 will get this prefix and suffix
        p2, s2 = get_prefix_suffix(ref1, y)
        #The ref1 and seq1 will get this prefix and suffix
        p1, s1 = get_prefix_suffix(ref2, y)
        #New sequences
        ref1 = p1 + ref1 + s1
        al_seqs = [p1 + seq1 + s1 for seq1 in al_seqs]
        seq2 = p2 + seq2 + s2
        aligned_sequences = al_seqs + [ref1, seq2]
        aligned = [s.replace('-', '') for s in aligned_sequences]
        remaining_mots = [m for m in mots if m not in aligned]
    return aligned_sequences

########## Get consensus from motifs ##########

def get_consensus(aligned, df_iupac_codes='../../data/Other/IUPAC_nt_code.csv', l=5):
    #Generate consensus
    df_consensus = pd.DataFrame(0, columns=[i for i in range(len(aligned[0]))], index=['U', 'G', 'C', 'A'])
    for s in aligned:
        for pos, letter in enumerate(s):
            if letter in df_consensus.index.tolist():
                df_consensus.loc[letter, pos] += 1
    #Make a full scored consensus
    df_sequence = pd.DataFrame(index=df_consensus.columns, columns=['nt', 'score'])
    for pos in df_consensus.columns.tolist():
        s = df_consensus[pos].where(df_consensus[pos] == df_consensus[pos].max()).dropna()
        #Only one max nucleotide
        if len(s)==1:
            df_sequence.loc[pos, 'nt'] = s.index.tolist()[0]
            df_sequence.loc[pos, 'score'] = s.values.tolist()[0]
        #If there is a tie - write both as sequence
        else:
            bases_tied = ', '.join(sorted(s.index.tolist()))
            code = df_iupac_codes.loc[df_iupac_codes['bases_sorted']==bases_tied, 'IUPAC_code'].values.tolist()[0]
            df_sequence.loc[pos, 'score'] = s.sum()
            df_sequence.loc[pos, 'nt'] = code
    #use a sliding window of kmer length to determine consensus sequence
    df_sequence['rolling_sum'] = df_sequence['score'].rolling(l).sum()
    end = df_sequence['rolling_sum'].idxmax()
    start = end - (l-1)
    consensus = ''.join(df_sequence.loc[start:end, 'nt'].values.tolist())
    return consensus, df_sequence

def get_consensus_from_motifs(mots, iupac_path='../../data/Other/IUPAC_nt_code.csv', l=4, match_s=2, missm_s=-1):
    aligned_sequences = get_alignment(mots, match_s, missm_s)
    #Import iupac code to make consensuses
    df_iupac_codes = pd.read_csv(iupac_path, sep='\t')
    df_iupac_codes['bases_sorted'] = [', '.join(sorted(s.split(' or '))) for s in df_iupac_codes['Base'].values.tolist()]
    consensus = get_consensus(aligned_sequences, df_iupac_codes, l)[0]
    return consensus, aligned_sequences

########## Plot a k-mer logo from alignment ##########
def GetPfm(aligned):
    """
    Input: aligned sequences from get_alignment.
    """
    # Generate position frequency matrix from alignment
    idx = [i for i in range(len(aligned[0]))]
    df_pfm = pd.DataFrame(0, index=idx, columns=['A', 'C', 'G', 'U', 'N', '-'])
    # print('Alignment:')
    for s in aligned:
        # print(s)
        for pos, letter in enumerate(s):
            df_pfm.loc[pos, letter] += 1
    # Remove positions in PFM, where no nucleotides
    df_pfm = df_pfm.loc[(df_pfm[['A', 'C', 'G', 'U']] != 0).any(axis=1)]
    df_pfm.reset_index(drop=True, inplace=True)
    return df_pfm


def PlotKmerLogoSeqlogo(kmer_list, outFileName, match_s=2, missm_s=-1):
    """
    Plots a weblogo from k-mer alignment.
    inputs: kmer list (with uracils instead of thymine)
    output: sequence logo
    """
    aligned = get_alignment(kmer_list, match_s, missm_s)
    df_pfm = GetPfm(aligned)
    # convert pfm to ppm and pwm
    cpm = seqlogo.CompletePm(pfm=df_pfm, alphabet_type='reduced RNA', background=0)
    ppm = seqlogo.Ppm(cpm.ppm, alphabet_type='reduced RNA', background=0)
    #plot seq logos
    try:
        seqlogo.seqlogo(ppm, filename=outFileName, format='pdf', size='medium', ic_scale=False)
    except ValueError: #Needed to add this to prevent value error: pdf not supported for plotting in console. It saves the pdfs anyway
        pass
    return df_pfm

def PlotKmerLogoLogomaker(kmer_list, match_s=2, missm_s=-1, colorScheme='default_RNA'):
    if colorScheme == 'default_RNA':
        colorScheme = {
        'A': '#48AD7E',
        'U': '#FC5252',
        'G': '#F5D741',
        'C': '#55BDF5',
        '-': 'white',
        'N': 'white',
        }
    if len(kmer_list) > 1:
        aligned = get_alignment(kmer_list, match_s, missm_s)
    elif len(kmer_list) == 1:
        aligned = kmer_list
    df_pfm = GetPfm(aligned)
    # convert pfm to ppm and pwm
    cpm = seqlogo.CompletePm(pfm=df_pfm, alphabet_type='reduced RNA', background=0)
    # convert pfm to position probability matrix
    ppm = cpm.ppm
    nts = ['A', 'U', 'G', 'C']
    ppm = ppm[nts]
    fig, ax = plt.subplots()
    logomaker.Logo(ppm, color_scheme=colorScheme, ax=ax)
    return fig, ax


def get_coverage(sequence, kmer_group, kmer_length):
    "Get sequence positions that are covered by k-mers in a given k-mer group."
    pos_dict = {pos:0 for pos in range(0, len(sequence))}
    for pos, nt in enumerate(sequence):
        if pos <= len(sequence)-kmer_length:
            kmer = sequence[pos: pos+kmer_length]
            if kmer in kmer_group:
                for p in range(pos, pos+kmer_length):
                    pos_dict[p] = 1
    return pos_dict

def get_percNtCovered(seq, kmer_group, kmer_length):
    "Returns the percentage of nt covered by a given k-mer group in a sequence."
    coverage = get_coverage(seq, kmer_group, kmer_length)
    covered_nt = sum(list(coverage.values()))
    total_nt = len(seq)
    perc_covered = covered_nt * 100 / total_nt
    return perc_covered
