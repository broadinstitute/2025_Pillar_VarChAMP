import pandas as pd
import numpy as np
from random import choice
from scipy.stats import norm


wt_gfp_threshold = 100


def get_wt_variability_d(pDEST_DUAL_avg_df, pDEST_DUAL_df):
    '''
    get mean WT, per plate to determine WT/WT variability to infer Allele/WT score
    '''
    mean_pDEST_empty = np.mean(pDEST_DUAL_df[pDEST_DUAL_df.well.isin(['A07','A08'])].GFP_mCherry_ratio)
    std_pDEST_empty = np.std(pDEST_DUAL_df[pDEST_DUAL_df.well.isin(['A07','A08'])].GFP_mCherry_ratio)

    empty_threshold = mean_pDEST_empty + 3*std_pDEST_empty

    wt_df = pDEST_DUAL_avg_df[(pDEST_DUAL_avg_df.mut_id == 0) & # get wildtype cells
                            (pDEST_DUAL_avg_df.avg_gfp >= wt_gfp_threshold) & # GFP higher than some absolute signal strength
                            (pDEST_DUAL_avg_df.avg_GFP_mCherry_ratio >= empty_threshold) # GFP:mCherry signal is higher than empty wells
                            ].copy()

    # Get WT GFP/mCherry ratios per ORF ID per plate
    wt_d = {}
    for _, row in wt_df.iterrows():
        key = (row.orf_id, row.pla)
        value = row.avg_GFP_mCherry_ratio

        if key not in wt_d:
            wt_d[key] = [value]
        else:
            wt_d[key].append(value)

    # Compute wt:wt variability by randomly pairing two wells with the same ORF on the same plate
    wt_ratio_l = []
    l = []
    for (orf_id, pla), avg_ratio in wt_d.items():
        l.append([orf_id, pla, np.mean(avg_ratio), len(avg_ratio)])
        while len(avg_ratio) >= 2:
            x1 = choice(avg_ratio)
            avg_ratio.remove(x1)
            x2 = choice(avg_ratio)
            avg_ratio.remove(x2)
            
            wt_ratio_l.append(x1/x2)

    wt_avg_df = pd.DataFrame(l, columns=['orf_id','pla','avg_ratio','n_wt'])
    wt_avg_df.n_wt.value_counts()

    wt_avg_d = {(orf_id, pla):avg_ratio for i, (orf_id, pla, avg_ratio, n_wt) in wt_avg_df.query("n_wt >= 2").iterrows()}

    return wt_avg_d, wt_ratio_l


def get_wt_variability_d_median(pDEST_DUAL_median_df,pDEST_DUAL_df):
    '''
    get median WT, per plate to determine WT/WT variability to infer Allele/WT score
    '''
    median_pDEST_empty = np.median(pDEST_DUAL_df[pDEST_DUAL_df.well.isin(['A07','A08'])].GFP_mCherry_ratio)
    std_pDEST_empty = np.std(pDEST_DUAL_df[pDEST_DUAL_df.well.isin(['A07','A08'])].GFP_mCherry_ratio)
    empty_threshold = median_pDEST_empty + 3*std_pDEST_empty

    wt_df = pDEST_DUAL_median_df[(pDEST_DUAL_median_df.mut_id == 0) &
                                 (pDEST_DUAL_median_df.median_gfp >= 100) &
                                 (pDEST_DUAL_median_df.median_GFP_mCherry_ratio >= empty_threshold)
                                 ].copy()

    # Get WT GFP/mCherry ratios per ORF ID per plate
    wt_d = {}
    for _, row in wt_df.iterrows():
        key = (row.orf_id, row.pla)
        value = row.median_GFP_mCherry_ratio

        if key not in wt_d:
            wt_d[key] = [value]
        else:
            wt_d[key].append(value)

    # Compute wt:wt variability by randomly pairing two wells with the same ORF on the same plate
    wt_ratio_l_median = []
    l = []
    for (orf_id, pla), median_ratio in wt_d.items():
        l.append([orf_id, pla, np.mean(median_ratio), len(median_ratio)])
        while len(median_ratio) >= 2:
            x1 = choice(median_ratio)
            median_ratio.remove(x1)
            x2 = choice(median_ratio)
            median_ratio.remove(x2)
            
            wt_ratio_l_median.append(x1/x2)

    wt_median_df = pd.DataFrame(l, columns=['orf_id','pla','median_ratio','n_wt'])
    wt_median_df.n_wt.value_counts()

    wt_median_d = {(orf_id, pla):median_ratio for i, (orf_id, pla, median_ratio, n_wt) in wt_median_df.query("n_wt >= 2").iterrows()}
    return wt_median_d, wt_ratio_l_median


def wt_log2fc_variability(pDEST_DUAL_avg_df, pDEST_DUAL_df):
    '''
    Maxime addition
    Instead of using Georges' approach to compute the assay's variability, which uses a step with random pairings, 
    Luke suggested to compute a STD from the Log2FC of all individual WT measurements, relative to the mean of the WT of each gene.
    This is a more robust approach, as it does not rely on random pairings.
    This function computes the STD of the Log2FC of all individual WT measurements, relative to the mean of the WT of each gene.
    Returns a tuple with the mean and the STD of the Log2FC of all individual WT measurements, relative to the mean of the WT of each gene.
    '''
    
    # Extract empty wells fluorescence measurements
    mean_pDEST_empty = np.mean(pDEST_DUAL_df[pDEST_DUAL_df.well.isin(['A07','A08'])].GFP_mCherry_ratio)
    std_pDEST_empty = np.std(pDEST_DUAL_df[pDEST_DUAL_df.well.isin(['A07','A08'])].GFP_mCherry_ratio)

    empty_threshold = mean_pDEST_empty + 3*std_pDEST_empty

    # Extract and filter WT data
    wt_df = pDEST_DUAL_avg_df[(pDEST_DUAL_avg_df.mut_id == 0.0) 
                          & (pDEST_DUAL_avg_df.avg_GFP_mCherry_ratio >= empty_threshold)
                          & (pDEST_DUAL_avg_df.avg_gfp >= wt_gfp_threshold)
                        ].copy()
    
    # Compute the WT mean for each orf and plate
    wt_mean = wt_df.groupby(['orf_id', 'pla'])['avg_GFP_mCherry_ratio'].mean().rename('WT_gene_avg_GFP_mCherry_ratio').reset_index()
    
    wt_df = wt_df.merge(wt_mean, on=['orf_id', 'pla'])
    
    # Compute the log2FC for each WT well/ log2 or not?
    wt_log2FC = np.log2(wt_df['avg_GFP_mCherry_ratio'] / wt_df['WT_gene_avg_GFP_mCherry_ratio'])
    
    # Compute the standard deviation and mean of the log2FC values
    wt_std = np.std(wt_log2FC)
    wt_mean = np.mean(wt_log2FC)
    
    return wt_std, wt_mean


def get_pDEST_DUAL_avg_allele_df(pDEST_DUAL_avg_df, wt_avg_d, wt_ratio_l, wt_std, wt_mean):
    pDEST_DUAL_avg_allele_df = pDEST_DUAL_avg_df[(pDEST_DUAL_avg_df.mut_id > 0)].copy()

    pDEST_DUAL_avg_allele_df['wt_GFP_mCherry_ratio'] = [wt_avg_d[orf_id, pla] if (orf_id, pla) in wt_avg_d else np.nan for orf_id, pla in zip(pDEST_DUAL_avg_allele_df.orf_id, pDEST_DUAL_avg_allele_df.pla)]
    pDEST_DUAL_avg_allele_df['allele_wt_ratio'] = pDEST_DUAL_avg_allele_df.avg_GFP_mCherry_ratio / pDEST_DUAL_avg_allele_df.wt_GFP_mCherry_ratio
    
    pDEST_DUAL_avg_allele_df['zscore'] = pDEST_DUAL_avg_allele_df.allele_wt_ratio.apply(lambda x: (x - np.mean(wt_ratio_l)) / np.std(wt_ratio_l))
    pDEST_DUAL_avg_allele_df['zcat'] = pDEST_DUAL_avg_allele_df.zscore.apply(lambda x: -2 if x <= -2 else -1 if x <= -1 else 2 if x >= 2 else 1 if x >= 1 else 0)

    # Maxime add
    pDEST_DUAL_avg_allele_df['allele_wt_log2fc'] = np.log2(pDEST_DUAL_avg_allele_df['allele_wt_ratio'])
    pDEST_DUAL_avg_allele_df['allele_wt_log2fc_zscore'] = pDEST_DUAL_avg_allele_df.allele_wt_log2fc.apply(
        lambda x: (x - wt_mean) / wt_std
    )
    pDEST_DUAL_avg_allele_df['allele_wt_log2fc_pvalue'] = pDEST_DUAL_avg_allele_df['allele_wt_log2fc_zscore'].apply(
        lambda z: 2 * norm.sf(np.abs(z))
    )
    
    # Jess add
    wt_log = np.log2(wt_ratio_l)
    wt_log_mean = np.mean(wt_log)
    wt_log_sd = np.std(wt_log)
    pDEST_DUAL_avg_allele_df['zscore_log2'] = pDEST_DUAL_avg_allele_df.allele_wt_ratio.apply(lambda x: (np.log2(x) - wt_log_mean) / wt_log_sd)

    return pDEST_DUAL_avg_allele_df


def get_pDEST_DUAL_median_allele_df(pDEST_DUAL_median_df, wt_median_d, wt_ratio_l_median):
    pDEST_DUAL_median_allele_df = pDEST_DUAL_median_df[(pDEST_DUAL_median_df.mut_id > 0)].copy()

    pDEST_DUAL_median_allele_df['wt_GFP_mCherry_ratio_median'] = [wt_median_d[orf_id, pla] if (orf_id, pla) in wt_median_d else np.nan for orf_id, pla in zip(pDEST_DUAL_median_allele_df.orf_id, pDEST_DUAL_median_allele_df.pla)]
    pDEST_DUAL_median_allele_df['allele_wt_ratio_median'] = pDEST_DUAL_median_allele_df.median_GFP_mCherry_ratio / pDEST_DUAL_median_allele_df.wt_GFP_mCherry_ratio_median

    pDEST_DUAL_median_allele_df['zscore_median'] = pDEST_DUAL_median_allele_df.allele_wt_ratio_median.apply(lambda x: (x - np.mean(wt_ratio_l_median)) / np.std(wt_ratio_l_median))
    pDEST_DUAL_median_allele_df['zcat_median'] = pDEST_DUAL_median_allele_df.zscore_median.apply(lambda x: -2 if x <= -2 else -1 if x <= -1 else 2 if x >= 2 else 1 if x >= 1 else 0)

    # Jess add
    wt_log = np.log2(wt_ratio_l_median)
    wt_log_mean = np.mean(wt_log)
    wt_log_sd = np.std(wt_log)
    pDEST_DUAL_median_allele_df['zscore_log2'] = pDEST_DUAL_median_allele_df.allele_wt_ratio_median.apply(lambda x: (np.log2(x) - wt_log_mean) / wt_log_sd)

    return pDEST_DUAL_median_allele_df