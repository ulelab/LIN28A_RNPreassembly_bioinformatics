import pandas as pd
from itertools import product
import numpy as np
import pybedtools as pbt
import GtfBedLibrary as gbl

def getName(f, method='eclip'):
    """
    Get a name that represents a dataset (protein name and cell line) from file path.
    """
    methods = ['eclip', 'iclip', 'parclip']
    if not method:
        name = f.split('/')[-1].split('.')[0]
    elif (method == 'eclip') or (method == 'parclip'):
        name = f.split('/')[-1].split('.')[0].split('_')[0].replace('-merged', '')
    elif method == 'iclip':
        name = f.split('/')[-1].split('.')[0].split('-')[0].upper().split('_')[0]
    else:
        print(method, 'not in', methods)
        print('wrong method specified, exiting...')
        return
    return name

def GetPEKAScoresFromTsv(file_list):
    """
    Combines PEKA-scores columns from a list of tsv files produced by PEKA into one table.
    Columns are filenames and rows are k-mers.
    """
    df = pd.DataFrame()
    for f in file_list:
        name = f.split('/')[-1]
        dft = pd.read_csv(f, sep='\t', index_col=0)
        df[name] = dft['PEKA-score']
    return df

def GetSignificantKmersFileList(file_list, p=0.05):
    """
    Returns a dict with filename as key and a list pf significant k-mers (p-value < p).
    """
    out = {}
    for f in file_list:
        # Get filename
        name = f.split('/')[-1]
        dft = pd.read_csv(f, sep='\t', index_col=0)
        significantKmers = dft.loc[dft['p-value'] < p].index.tolist()
        out[name] = significantKmers
    return out

def GetKmerRanksFromPEKAScore(PekaScoreSeriesOrDf):
    """
    Ranks scores in form of pandas series or dataframe (scores in columns).
    Highest score gets the lowest rank.
    NaN values are given the highest ranks (bottom).
    For a group of indices that have the same value, an average rank is assigned.
    """
    ranked = PekaScoreSeriesOrDf.rank(axis='index', ascending=False, na_option='bottom', method='average')
    return ranked

def GetNtComposition(kmer_list):
    """
    Returns a dictionary of nucleotide composition for a list of k-mers.
    -----
    input: k-mer list
    -----
    output: dict of nt content per k-mer
    """
    kmer_list = [k.upper() for k in kmer_list]
    U_counts = sum([m.count('U') for m in kmer_list])/len(kmer_list)
    A_counts = sum([m.count('A') for m in kmer_list])/len(kmer_list)
    C_counts = sum([m.count('C') for m in kmer_list])/len(kmer_list)
    G_counts = sum([m.count('G') for m in kmer_list])/len(kmer_list)
    nt_comp = {
        'U' : U_counts,
        'A' : A_counts,
        'C' : C_counts,
        'G' : G_counts,
    }
    return nt_comp

def ParseNtxnNoxnFromOutfile(f):
    """
    Reads the number of thresholded crosslinks and reference crosslinks in each region analyzed with PEKA from
    standard output file.
    -----
    Input: file which contains standard output from PEKA.
    -----
    Returns: A dict with transcript regions as keys and ntxn and noxn as values
    """
    outDict = {'ntxn': {}, 'noxn': {}}
    with open(f, 'r') as file:
        for line in file:
            line = line.strip('\n')
            ntxn = None
            noxn = None
            if ('ntxn' in line) and ('on' in line):
                ntxn = int(line.split(' ')[1])
                region = line.split(' ')[-1]
                outDict['ntxn'][region] = ntxn
            if ('noxn' in line) and ('on' in line):
                noxn = int(line.split(' ')[1])
                region = line.split(' ')[-1]
                outDict['noxn'][region] = noxn
    return outDict


def ConvertScoresKmerLen(DfScores, currentLen, newLen):
    """
    Converts PEKA-scores between k-mers of different lengths.
    -----
    Input: Table with PEKA-scores, kmers in index and dataset name in column.
    -----
    Output: Table with PEKA-scores for k-mers of desired length.
    """
    # Get possible k-mers of a desired length
    possible_kmers = []
    for i in product('ACGU', repeat=newLen):
        possible_kmers.append("".join(i))
    DfOut = pd.DataFrame(columns=[c for c in DfScores.columns], index=possible_kmers)
    # If the new length is lower than old length, take an average across all longer kmers,
    # that contain a shorter k-mer.
    if newLen < currentLen:
        for kmer in possible_kmers:
            df_t = DfScores.filter(like=f'{kmer}', axis=0)
            longKmers = df_t.index.tolist()
            longKmers2 = []
            for k in longKmers:
                splitLong = [k[end-newLen:end] for end in range(newLen, len(k) + 1)]
                longKmers2.extend([k for _ in range(splitLong.count(kmer))])
            vals = [DfScores.loc[k, :].values.tolist() for k in longKmers2]
            df_t = pd.DataFrame(vals, columns=DfScores.columns)
            DfOut.loc[kmer, :] = df_t.mean(axis=0)
    elif newLen == currentLen:
        print('No difference in k-mer length, exiting.')
        return
    else:
        # Kmers are short and we want to convert scores to longer k-mer.
        # We take an arithmetic average of all shorter k-mers that comprise a longer k-mer.
        for kmer in possible_kmers:
            shortKmers = [kmer[end-currentLen:end] for end in range(currentLen, len(kmer) + 1)]
            vals = [DfScores.loc[k, :].values.tolist() for k in shortKmers]
            df_t = pd.DataFrame(vals, columns=DfScores.columns)
            DfOut.loc[kmer, :] = df_t.mean(axis=0)
    return DfOut

def GetRecall(ClipPekaScores, InVitroRanks, topCLIP=50, topInVitro=20):
    """
    Computes recall for a table of CLIP datasets, that contains PEKA score.
    -----
    Input:
    - Table of peka scores across k-mers. Names of Clip datasets are in columns, k-mers are in index.
    - Table of in-vitro k-mer ranks. Names of proteins are in columns, k-mers are in index.
    - topCLIP: number of top k-mers from CLIP to take as top.
    - topInVitro: number of top k-mers from in-vitro data to take as top.
    -----
    Output:
    Table of recall values for a particular dataset.
    If a dataset does not have corresponding in vitro data available, it gets a NaN value for recall.
    """
    Recall = pd.Series()
    InVitroRanks.columns = [c.upper() for c in InVitroRanks.columns]
    for dataset in ClipPekaScores.columns:
        RBP = dataset.replace('HepG2-', '').replace('K562-', '').upper().split('_')[0]
        if RBP in InVitroRanks.columns:
            topkmersCLIP = ClipPekaScores[dataset].sort_values(ascending=False).head(topCLIP).index.tolist()
            # We have in vitro ranks, so I sort in ascending order
            topkmersinVitro = InVitroRanks[RBP].sort_values(ascending=True).head(topInVitro).index.tolist()
            #  Calculate recall
            recall = len([k for k in topkmersinVitro if k in topkmersCLIP]) / topInVitro
            Recall[dataset] = recall
        else:
            Recall[dataset] = np.nan
    return Recall

def CountXlsInRegion(files, regions_file, outpath, prefix='', method=None):
    """
    Counts the number of crosslink positions and cDNA counts that fall into a particular genomic region.
    Regions are defined with a gtf file that contains a region as a feature in the second column.
    Such region file can be generated from GENCODE / Ensembl annotation by running iCount segment.
    -----
    Input:
    - list of crosslink files
    - region file (can be genome level segmentation or a file that contains other features - i.e. repeats ...).
    - path to a folder where output tables should be saved
    - prefix for output filenames, defaults to empty string
    - method {eclip, iclip, parclip} - for which files we are generating names, default is None.
    -----
    Output: Saves two tables, one with summed cDNA counts (scores) for a particular region and
    one with the number of positions in each region (equal to number of entries/lines for a
    particular feature).
    """
    regions = gbl.Gtf2Bed(regions_file)

    cdnaRows = []
    ntxnRows = []
    for i, f in enumerate(files):
        name = getName(f, method=method)
        txn = pbt.BedTool(f)
        txn = gbl.AnnotateBed(txn, regions)
        df = gbl.ParseBedToolToDf(txn)
        df.rename(columns={6: 'region'}, inplace=True)
        df_summedCna = df.groupby('region')['score'].sum()
        df_countedTxnPositions = df.groupby('region')['score'].count()
        Cdna = df_summedCna.to_dict()
        Cdna['dataset'] = name
        Ntxn = df_countedTxnPositions.to_dict()
        Ntxn['dataset'] = name
        # Append
        cdnaRows.append(Cdna)
        ntxnRows.append(Ntxn)
        print(i+1, '/', len(files))
    CdnaCounts = pd.DataFrame(data=cdnaRows)
    NtxnPositions = pd.DataFrame(data=ntxnRows)
    NtxnPositions.set_index('dataset', inplace=True)
    CdnaCounts.set_index('dataset', inplace=True)
    CdnaCounts.to_csv(f'{outpath}/{prefix}CdnaCountTxnPerRegion.tsv', sep='\t')
    NtxnPositions.to_csv(f'{outpath}/{prefix}NxnPositionsPerRegion.tsv', sep='\t')
    return

def FilterTxnInRepeats(bt_object, fasta):
    """Gets all sequences (crosslinks) from bed file and drops those which fall on a lowercase letter in masked genome."""
    df_sites = pd.read_csv(bt_object.fn,
                names=['chrom', 'start', 'end', 'name', 'score', 'strand'],
                sep='\t',
                header=None,
                dtype={'chrom': str, 'start': int, 'end': int, 'name': str, 'score': int, 'strand': str})
    seq_tab = bt_object.sequence(s=True, fi=fasta, tab=True)
    df_sites['seq'] = [line.split("\t")[1].strip() for line in open(seq_tab.seqfn)]
    df_sites = df_sites[df_sites['seq'].str.isupper() == True]
    df_sites.drop(['seq'], axis=1, inplace=True)
    return pbt.BedTool.from_dataframe(df_sites)

def ComputeSimilarity(dfRanks, ntop):
    df_similarity = pd.DataFrame(columns=['similarity_score'], index=dfRanks.columns)
    # Where rank is lower or equal to ntop, give score of 1, else 0
    df_sim = dfRanks.apply(lambda x: np.where(x <= ntop, 1, 0))
    for c in dfRanks.columns:
        # Find kmers where rank is less than 50 for a dataset, get all other columns except for a given dataset.
        # Across other columns, sum the 1nes (this represents the number of k-mers that are in top 50 for other datasets) and divide by 50 to get proportion
        pairwise_sim_scores = df_sim.copy().loc[df_sim[c]==1, [i for i in dfRanks.columns if i!=c]].sum(axis='index').divide(ntop)
        # Get a mean of proportions across all datasets
        similarity_score = pairwise_sim_scores.mean()
        df_similarity.loc[c, 'similarity_score'] = similarity_score
    return df_similarity
