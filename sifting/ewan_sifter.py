import sys
import os
import json
import pandas as pd
import numpy as np
import time
import glob
import shutil
from sqlalchemy import create_engine


sys.path.append("/hercules/results/jjawor/MMGPS_sifter")
from MMGPS_sifter.sifting.sifters.ewanfilter import cluster_cands, spatial_rfi,\
    filtering, known_filter


def a_to_pdot(P_s, acc_ms2):
    LIGHT_SPEED = 2.99792458e8                 # Speed of Light in SI
    return P_s * acc_ms2 /LIGHT_SPEED


def period_modified(p0,pdot,no_of_samples,tsamp,fft_size):
    if (fft_size==0.0):
        return p0 - pdot*float(1<<(int(no_of_samples).bit_length()-1))*tsamp/2
    else:
        return p0 - pdot*float(fft_size)*tsamp/2


def excise_birdies(df, birdie_list_file, obs_meta_data, nharmonics, p_tol):
    # Get a temporary ordered array to work with
    df = df.reset_index()
    
    # Get candidate periods and dms
    cand_mid_periods = df['period'].to_numpy()
    cand_dms = df['dm'].to_numpy()
    cand_snrs = df['snr'].to_numpy()
    cand_accs = df['acc'].to_numpy()

    # Modify periods to starting epoch reference
    print("Calculating modified period based on starting epoch reference")
    tsamp = obs_meta_data['tsamp']
    fft_size = obs_meta_data['fft_size']
    nsamples = obs_meta_data['nsamples']

    mod_periods=[]
    pdots = []
    for i in range(len(cand_mid_periods)):
        pdot = a_to_pdot(cand_mid_periods[i], cand_accs[i])
        mod_periods.append(period_modified(cand_mid_periods[i], pdot, nsamples, tsamp, fft_size))
        pdots.append(pdot)

    cand_periods = np.asarray(mod_periods, dtype=float)
    cand_freqs = 1 / cand_periods
    cand_mid_freqs = 1 / cand_mid_periods

    known_rfi_indices = known_filter.get_known_rfi(
        cand_mid_freqs, birdie_list_file, nharmonics, p_tol)
    print("Number of RFI instances: %d"%len(known_rfi_indices))
    df.iloc[known_rfi_indices,:].to_csv("known_rfi_cands.csv")
    df_excised = df[~df.index.isin(known_rfi_indices)]
    return df_excised


def filter_clusters_ewan(df, cluster_df=None):
    """Remove all cands from clusters that are spatially associated with RFI and have 
    low nassoc. Also remove all cand that aren't the strongest in their cluster"""

    # Select only the top cluster candidates
    df = df[(df['strongest_in_cluster'] == True)]

    # Remove candidates that are spatially associated with RFI and/or have low
    # nassoc
    condition = ((df['spatial_rfi'] == False) & (df['low_nassoc'] == False))
    df['passed_sifting'] =  0
    df.loc[condition, 'passed_sifting'] = 1

    # Other logic can be implemented to filter out candidates (DM thresh, etc)


    return df

def score_ewan_pointing(cluster_df):
    """Score clusters in cands_df. Used to score clusters in the whole pointing"""

    sorted_indices_cand = cluster_df.sort_values("max_snr", ascending=False).index
    cluster_df.loc[sorted_indices_cand, 'pointing_rank'] = range(len(sorted_indices_cand))

    return cluster_df


def score_ewan_beam(df, beam_id, cluster_df=None):
    """Score clusters in a target beam (Ewan sifter version)
    """

    # Select only the target beam
    target_beam_df = df[(df['beam_id'] == beam_id)]

    # Assign the rank of the cluster to the candidates
    target_beam_df = target_beam_df.sort_values('snr', ascending=False)\
                    .reset_index(drop=True)
    target_beam_df["beam_rank"] = target_beam_df.index

    return target_beam_df


def run_ewansifter(df_cands_ini, obs_meta_data, threshold, birdies=None, config=baseconfig,
         debug=False, ptol=5e-4, harmonics=16, output=os.getcwd(),
         psr_filter=False):
    """Run the candidate filtering pipeline.
    Args:
        df_cands_ini (DataFrame): DataFrame containing initial candidates
        obs_meta_data (dict): The metadata for the observation that the
            candidates were found in
        threshold (float): S/N threshold to apply when reading candidates from
            input files
        birdies (str): Known birdie list name
        config (str): The path to the config file containing the filtering
            arguments
        debug (bool): If True, write out intermediate files
        ptol (float): The period tolerance for harmonic clustering
        harmonics (int): Harmonic number to search upto
        output (str): The root output directory
        psr_filter (bool): If True, filter out known PSRs. NOTE: This is not
            currently implemented!
    """


    # Load config
    with open(config) as json_data_file:
        config = json.load(json_data_file)

    if debug:
        df_cands_ini.to_csv("all_candidates.csv")
     
    # Remove DM candidates below 2 pc cm^-3 - 0 DM + too low a DM to probably be real
    print("Removing candidates below DM of 2 pc cm^-3") 
    df_cands_ini = df_cands_ini[df_cands_ini['dm'] > 2.0]
    
    # Apply S/N threshold
    print("Removing candidates below S/N {}".format(threshold))
    df_cands_ini = df_cands_ini[df_cands_ini['snr'] >= threshold]

    # Reindex the data frame to allow interoperability of numpy and
    # pandas indexing 
    df_cands_ini.reset_index(inplace=True)

    if debug:
        df_cands_ini.to_csv("sn_and_dm_filtered_candidates.csv")

    # We don't use a birdies list anymore
    # if birdies is not None:
    #     df_cands_ini = excise_birdies(df_cands_ini, birdies, obs_meta_data,
    #                                   harmonics, p_tol)

    # Create clusters from remaining candidates
    df_cands_clustered = cluster_cands.cluster_cand_df(
        df_cands_ini, obs_meta_data, config)

    # Find spatial RFI and write out details about clusters
    df_clusters = spatial_rfi.label_spatial_rfi(df_cands_clustered, config)

    # Label bad clusters
    df_cands_filtered, df_clusters_filtered = filtering.filter_clusters(
        df_cands_clustered, df_clusters, config, output)

    if psr_filter:
        # Filter away known PSRs. WARNING: DOSEN'T WORK R/N
        # Write out known pulsar file
        df_cands_ini.iloc[known_psr_indices,:].to_csv(f"{output}_known_psr_cands.csv")

        # Write out known pulsar harmonics file
        df_cands_ini.iloc[known_ph_indices,:].to_csv(f"{output}_known_ph_cands.csv")

    # Write out reduced candidate list
    # df_cands_filtered.to_csv(f"{output}_cands.csv")
    # Write out cluster list
    # df_clusters_filtered.to_csv(f"{output}_clusters.csv")


    return df_cands_filtered, df_clusters_filtered


def main(args):

    # Load the config file
    outfile_cluster = os.path.join(args.outpath, f"cluster_score.csv")
    outfile = os.path.join(args.outpath, f"cand_score.csv")

    base_cand_df = pd.read_csv(args.cand_file)
    obs_meta_df = pd.read_csv(args.metafile)

    try: 
        # Launch sifter
        cands_df, cluster_df = run_ewansifter(base_cand_df, obs_meta_df,
                                              threshold=args.threshold,
                                              config=args.config)

        # Apply thresholding to remove bad candidates
        cands_df = filter_clusters_ewan(cands_df)

        # Save the scored dfs
        cands_df.to_csv(outfile)


        # # Create the pointing-wise score
        # scored_cluster_df = score_ewan_pointing(cluster_df)

        # # Create the beam-wise score
        # label1_beam_IDs = cands_df[cands_df['label'] == 1]['beam_id'].unique()
        # scores = []
        # for beam_id in label1_beam_IDs:
        #     scores.append(score_ewan_beam(cands_df, beam_id))
        # scored_cands_df = pd.concat(scores)

        # scored_cluster_df.to_csv(outfile_cluster)
        # # Save the scored dfs
        # scored_cands_df.to_csv(outfile)

    except Exception as e:
        logging.error(f"Error in Ewan sifter. Pointing {pid}: {e}")


if __name__ == '__main__':

    parser = argparse.ArgumentParser(
    description='Sift candidates specified in a file. The candidates must also'
                ' exist in the sifterDB.')
    parser.add_argument('-o', '--outpath', help='Output path to save results',
                        default=os.getcwd(), type=str)
    parser.add_argument('-c', '--config', help='Path to config file',
                        required=True, type=str)
    parser.add_argument('-m', '--metafile', help='Path to metafile', required=True, type=str)
    parser.add_argument('-t', '--threshold', help='S/N threshold to apply',
                        default=7.0, type=float)
    parser.add_argument('cand_file', help='Path to candidate file', type=str)
    args = parser.parse_args()

    main(args)