# GTF and Bed Library
import pandas as pd
import pybedtools as pbt
import csv

def ReadGtf(segmentation):
    df_segment = pd.read_csv(segmentation,
                             sep='\t',
                             names=['chrom', 'source', 'feature', 'start', 'end', 'name', 'strand', 'name2', 'annotations'],
                             header=None,
                             comment='#',
                             dtype={
                                    "chrom": str,
                                    "source": str,
                                    "feature": str,
                                    "start": int,
                                    "end": int,
                                    "name": str,
                                    "strand": str,
                                    "name2": str,
                                    "annotations": str,
                                    },
                                    quoting=csv.QUOTE_NONE)
    return df_segment

def ReadBed(f):
    """Parse BedFilewith 6 or more columns to pandas.DataFrame."""
    cols = {0: "chrom", 1: "start", 2: "end", 3: "name", 4: "score", 5: "strand"}
    df = pd.read_csv(
        f,
        sep="\t",
        header=None,
        dtype={0: str, 1: int, 2: int, 3: str, 5: str},
    )
    df.rename(columns=cols, inplace=True)
    return df

def Gtf2Bed(gtf, segment=True, fromFile=True):
    # Get pandas dataframe
    if fromFile:
        df_gtf = ReadGtf(gtf)
    else:
        df_gtf = gtf
    # Conert to bed format
    if segment:
        bed_gtf = df_gtf.assign(start=df_gtf['start']-1, score=0)[['chrom', 'start', 'end', 'feature', 'score','strand', 'annotations']]
    else:
        bed_gtf = df_gtf.assign(start=df_gtf['start']-1, score=0)[['chrom', 'start', 'end', 'name', 'score','strand']]
    bed_gtf = pbt.BedTool.from_dataframe(bed_gtf, quoting=csv.QUOTE_NONE).sort()
    return bed_gtf


def Fai2Bed(fai):
    df_chromosomes = pd.read_csv(fai, sep='\t', header=None, names=['chr', 'end', 'offset', 'linebases', 'linewidth'])
    df_chromosomes = df_chromosomes[['chr', 'end']].assign(start=0, name='.', score=0)
    df_chromosomes_p = df_chromosomes.copy()
    df_chromosomes_p['strand'] = '+'
    df_chromosomes_p = df_chromosomes_p[['chr', 'start', 'end', 'name', 'score', 'strand']]
    df_chromosomes_m = df_chromosomes.copy()
    df_chromosomes_m['strand'] = '-'
    df_chromosomes_m = df_chromosomes_m[['chr', 'start', 'end', 'name', 'score', 'strand']]
    df_chromosomes = pd.concat([df_chromosomes_p, df_chromosomes_m], ignore_index=True)
    bed_chr = pbt.BedTool.from_dataframe(df_chromosomes).sort()
    return(bed_chr)

def MapXlCdnaCounts(IntervalBed, XlBed):
    IntervalBed = IntervalBed.sort()
    XlBed = XlBed.sort()
    BedCounts = IntervalBed.map(XlBed, c=5, o='sum', null=0, s=True, nonamecheck=True).sort()
    return BedCounts

def ParseBedToolToDf(f):
    """Parse BedTool with 6 or more columns to pandas.DataFrame."""
    cols = {0: "chrom", 1: "start", 2: "end", 3: "name", 4: "score", 5: "strand"}
    df = pd.read_csv(
        f.fn,
        sep="\t",
        header=None,
        dtype={0: str, 1: int, 2: int, 3: str, 5: str},
    )
    df.rename(columns=cols, inplace=True)
    return df

def AnnotateBed(BedTool, BedAnnotation, MapCols=[4]):
    """Annotates BedTool with another BedTool, which has feature in the 4th column
    and more detailed annotations in the """
    BedTool = BedTool.sort()
    BedAnnotation = BedAnnotation.sort()
    return BedTool.map(BedAnnotation, s=True, c=MapCols, o='collapse', nonamecheck=True).sort()

# def BedCombine(f1, f2, bedtool=False, mergedist=0):
#     """
#   Combine two bed files or two BedTool objects into one BedTool object.
#   This merges bookended features if mergedist is 0.
#   """
#     if not bedtool:
#         r1 = pbt.BedTool(f1).sort()
#         r2 = pbt.BedTool(f2).sort()
#     else:
#         r1 = f1.sort()
#         r2 = f2.sort()
#     concat = r1.cat(r2, s=True, force_truncate=False, postmerge=False).sort()
#     merged = concat.merge(s=True, c=[4, 5, 6], o=['distinct', 'sum', 'distinct'], d=mergedist).sort()
#     return merged

def CombineReplicates(f1, f2, bedtool=False):
    """
    Combines crosslinks from replicates.
    """
    if not bedtool:
        r1 = pbt.BedTool(f1).sort()
        r2 = pbt.BedTool(f2).sort()
    else:
        r1 = f1.sort()
        r2 = f2.sort()
    # Concatenate bedtools
    concat = r1.cat(r2, s=True, force_truncate=False, postmerge=False).sort()
    # Groupby chrom, start, end, strand
    # Sum scores, and keep first name
    combined = concat.groupby(g=[1, 2, 3, 6], c=[5, 4], o=['sum', 'first'])
    # Convert to df to rearrange columns in correct order (chrom, start, end, name, score, strand)
    df = pd.read_csv(
        combined.fn,
        sep="\t",
        header=None
    )
    df = df[[0, 1, 2, 5, 4, 3]]
    # Return bedtool
    btOut = pbt.BedTool.from_dataframe(df).sort()
    return btOut

def ExtractSequencesToList(df_sites, fasta, slop=False, fai=None, slop_bp=3):
    """Get genome sequences in peaks."""
    sites_extended = pbt.BedTool.from_dataframe(df_sites).sort()
    if slop:
        sites_extended = sites_extended.slop(b=slop_bp, g=fai)
    seq_tab = sites_extended.sequence(s=True, fi=fasta, tab=True)# If s=True, the sequence is reverse complemented.
    return [line.split("\t")[1].strip() for line in open(seq_tab.seqfn)]

def Bedgraph2Bed(f, outputFolder, skiprows=[0], header=None):
    fname = f.split('/')[-1].split('.')[0]
    if skiprows:
        df = pd.read_csv(f, sep='\t', skiprows=skiprows, header=header)
    else:
        df = pd.read_csv(f, sep='\t', header=header)
    df.columns = ['chrom', 'start', 'end', 'BgScore']
    df = df.assign(
        name = '.',
        strand = ['-' if v < 0 else '+' for v in df.BgScore.values.tolist()],
        score = df['BgScore'].abs()
    )
    df = df[['chrom', 'start', 'end', 'name', 'score', 'strand']]
    df.to_csv(f'{outputFolder}/{fname}.bed.gz', compression='gzip', index=False, header=False, sep='\t')
    return df

def Bed2BedGraph(f, outputFolder):
    fname = f.split('/')[-1].split('.')[0]
    df = ReadBed(f)
    print(df.head())
    df['BgScore'] = df.apply(lambda row: -(row.score) if row.strand =='-' else row.score, axis='columns')
    df = df[['chrom', 'start', 'end', 'BgScore']].sort_values(by=['chrom', 'start'])
    df.to_csv(f'{outputFolder}/{fname}.bedgraph.gz', compression='gzip', index=False, header=False, sep='\t', quoting=csv.QUOTE_NONE)
    return df




