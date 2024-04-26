import xml.etree.ElementTree as ET
import sys, os, subprocess
import argparse, errno
import pandas as pd
import numpy as np

parser = argparse.ArgumentParser(description='Fold Peasoup candidates')
parser.add_argument('-i', '--input_file', help='Name of the input xml file', type=str)
parser.add_argument('-ps', '--period_start', help='Spin Period start in seconds', type=float)
parser.add_argument('-pe', '--period_end', help='Spin Period end in seconds', type=float)
parser.add_argument('-t', '--target_name', help='Target name', type=str)


args = parser.parse_args()

xml_file = args.input_file
tree = ET.parse(xml_file)
root = tree.getroot()
header_params = root[1]
search_params = root[2]
candidates = root[6]
tstart = float(header_params.find("tstart").text)
fft_size = float(search_params.find("size").text)
tsamp = float(header_params.find("tsamp").text)
tobs = tsamp * fft_size
tobs_days = tobs/(3600. * 24.)
new_epoch = tstart + tobs_days/2
period_start = args.period_start
period_end = args.period_end
target_name = args.target_name
ignored_entries = ['candidate', 'opt_period', 'folded_snr', 'byte_offset', 'is_adjacent', 'is_physical', 'ddm_count_ratio', 'ddm_snr_ratio']
rows = []
for candidate in candidates:
    cand_dict = {}
    for cand_entry in candidate.iter():
        if not cand_entry.tag in ignored_entries:
            cand_dict[cand_entry.tag] = cand_entry.text

    cand_dict['cand_id_in_file'] = candidate.attrib.get("id")
    rows.append(cand_dict)

df = pd.DataFrame(rows)
# Select the top 50 candidates to fold
df.sort_values(by='snr', ascending=False, inplace=True)
df = df.head(50)
df['Epoch_Mid'] = new_epoch
df = df.astype({"snr": float, "dm": float, "period": float, "nh":int, "acc": float, "nassoc": int, "cand_id_in_file": int})
cands_selected = df.loc[(df['period'] >= period_start) & (df['period'] <= period_end)]

for index, row in cands_selected.iterrows():
    with open('predictor_candidate_%d' %row['cand_id_in_file'], 'w') as f:
        f.write("SOURCE: " + target_name + '\n' + \
                "EPOCH: " + str(new_epoch) + '\n' + \
                "PERIOD: " + str(row['period']) + 's' + '\n' + \
                "DM: " + str(row['dm']) + '\n' + \
                "ACC: " +  str(row['acc']) + " (m/s/s)" + '\n')
