import pandas as pd
import pybedtools as pbt
from glob import glob
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import os
import numpy as np
import math
import dask.bag as db
from scipy import stats


def readBed(f):
    return pd.read_csv(f, sep='\t', header=None)

def createOutputDir(outpath):
    os.makedirs(outpath, exist_ok=True)
    return

def ExplodeXls(df, use_scores):
    if use_scores:
        df = df.loc[df.index.repeat(df[4])]
        df[4] = 1
    else:
        df[4] = 1
    return df

def getMatrix(df):
    cols = df.columns.tolist()
    # Split into positive and negative entries
    dfpos = df.loc[df[5]=='+'].copy()
    dfneg = df.loc[df[5]=='-'].copy()
    valspos = dfpos.groupby(cols[:-2])[cols[-1]].apply(list)
    # Reverse values on negative entries so they are three to five prime
    valsneg = dfneg.groupby(cols[:-2])[cols[-1]].apply(list).apply(lambda lst: lst[::-1])
    array = np.array(valspos.values.tolist() + valsneg.values.tolist())
    return array.astype(float)

def CoverageOnChunk(chunk, bt, libsize, norm):
    # Convert to bedtool
    bedRegions = pbt.BedTool.from_dataframe(chunk)
    # get raw coverages
    Coverage = bedRegions.coverage(bt, d=True, s=True, nonamecheck=True).to_dataframe(disable_auto_names=True, header=None)
    # normalise
    if norm == 'libsize':
    # normalise by libsize to get CPM in region
        Coverage.iloc[:,-1] = Coverage.iloc[:,-1] * (10**6 / libsize)
    elif norm == 'by_reg':
        # normalise by libsize and region value
        Coverage.iloc[:,-1] = Coverage.iloc[:,-1] * (10**6 / libsize)
        Coverage.iloc[:,-1] = Coverage.iloc[:,-1] / Coverage.iloc[:,4]
    else:
        print('Invalid argument for norm. Exiting.')
        sys.exit()
    coverage_array = getMatrix(Coverage)
    return coverage_array

def getCoverage(xls, dfRegions, chunk_size, use_scores, norm, pos_limits):
    libsizes = {}
    chunks_d = {}
    for f in xls:
        # Get filename to construct outfile name
        outname = f.split('/')[-1]
        chunks_d[outname] = []
        df = ExplodeXls(readBed(f), use_scores)
        # get library-size (n features in file)
        ls = df[4].sum()
        libsizes[outname] = ls
        # Convert to bedtool
        bt = pbt.BedTool.from_dataframe(df).sort()
        # Process in parallel
        chunks = np.array_split(dfRegions, math.ceil(len(dfRegions) / chunk_size))
        bag = db.from_sequence(chunks)
        coverage = bag.map(lambda x: CoverageOnChunk(x, bt=bt, libsize=ls, norm=norm))
        chunks_d[outname] = np.vstack(coverage.compute())
    # Collate - compute mean coverage across regions and 95% ci
    dfList = [pd.DataFrame()]
    for k, v in chunks_d.items():
        mean = np.mean(v, axis=0) # Compute a mean across region
        std = np.std(v, axis=0, ddof=1) # use ddof=1 for unbiased estimate of standard deviation
        # compute the 95% confidence interval for each column
        n = v.shape[0]
        ci = stats.t.interval(0.95, n-1, loc=mean, scale=std/np.sqrt(n))
        df = pd.DataFrame()
        df[f'mean'] = list(mean)
        df[f'ciMin'] = ci[0]
        df[f'ciMax'] = ci[1]
        df['Sample'] = k
        if pos_limits != None:
            df['Position'] = [i for i in range(pos_limits[0], pos_limits[1]+1)]
        else:
            df['Position'] = df.index.tolist()
        dfList.append(df)
    df_out = pd.concat(dfList)
    return df_out


def smoothDf(df, sw):
    cols = [c for c in df.columns.tolist() if c in ['ciMax', 'ciMin', 'mean']]
    df[cols] = df[cols].rolling(sw, center=True,  win_type='triang', axis=0).mean()\
                    .fillna(axis=0, method='ffill')\
                    .fillna(axis=0, method='bfill')
    return df


def main(regions, xls, outpath='./results', use_scores=True, sw=5, norm='libsize', chunk_size=2000, pos_limits=None):
    """
    norm: values with which to normalise xls. Options are 'libsize' and 'by_reg'. Default is 'libsize'.
    If option 'by_reg' is selected, the user should provide normalization value for each region in the score column.
    This value will be used in addition to normalization by library size.
    """
    # Create output directory
    createOutputDir(outpath)
    # Read regions into dataframe
    df_regions = readBed(regions)
    # Drop duplicated regions by chrom, start, end, strand
    df_regions = df_regions.drop_duplicates(subset=[0,1,2,5])
    # Check regions length
    length = (df_regions[2] - df_regions[1])
    if length.nunique() == 1:
        print(f"All regions have the same length: {length.iloc[0]}")
        region_sites = pbt.BedTool.from_dataframe(df_regions).sort()
    else:
        print("Regions don't have the same length. Exiting.")
        sys.exit()
    n_regions = len(region_sites)
    # Get library-normalised xl coverage within regions
    dfPlot = getCoverage(xls, df_regions, chunk_size, use_scores, norm, pos_limits)
    # Impute missing ci limits
    dfPlot[['ciMax', 'ciMin']] = dfPlot[['ciMax', 'ciMin']].apply(lambda col: col.fillna(dfPlot['mean'], axis='index'))
    # Save unsmoothed coverage values
    dfPlot.to_csv(f'{outpath}/coverage_unsmoothed_norm-{norm}.tsv', sep='\t', index=False)
    # Smooth and fill nan values with the closest valid observation
    dfPlot = dfPlot.groupby('Sample').apply(smoothDf, sw=sw)
    # Save smoothed data
    dfPlot.to_csv(f'{outpath}/coverage_smoothed_norm-{norm}.tsv', sep='\t', index=False)
    # Plot metaprofiles
    sampleNames = dfPlot['Sample'].unique()
    colorMap = {sample: sns.color_palette('husl', len(sampleNames))[i] for i, sample in enumerate(sampleNames)}

    fig, ax = plt.subplots()
    for s, df in dfPlot.groupby('Sample'):
        sns.lineplot(data=df, x='Position', y='mean', hue='Sample', palette=[colorMap[s]], ax=ax)
        ax.fill_between(df.Position, df.ciMin, df.ciMax, alpha=0.2, color=colorMap[s])
        if norm == 'libsize':
            ax.set_ylabel('mean CPM across regions')
        elif norm == 'by_reg':
            ax.set_ylabel('mean CPM per RegVal across regions')
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5))
    fig.savefig(f'{outpath}/metaprofile_norm-{norm}_Coverage.pdf', bbox_inches='tight')
    pbt.helpers.cleanup()
    return
