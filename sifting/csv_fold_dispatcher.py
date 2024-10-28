import numpy as np
import glob
import xml.etree.ElementTree as ET
import os
import optparse
import re
import subprocess
import itertools
import logging
import sys
import pandas as pd
import pika_process



log = logging.getLogger('manual_presto_fold')
FORMAT = "[%(levelname)s - %(asctime)s - %(filename)s:%(lineno)s] %(message)s"
logging.basicConfig(format=FORMAT)

log.setLevel('INFO')



parser = optparse.OptionParser()
pika_process.add_pika_producer_opts(parser)
parser.add_option('-C',type=str,help='Path to CSV file post candidate filtering to fold',dest="csv_path")
parser.add_option('-B',type=int,help='Batches',dest="batch_no",default=15)
parser.add_option('-M',type=str,help='Mask Path',dest="mask_path",default="/beegfs/PROCESSING/TRAPUM/RFIFIND_masks/Ter5_16apr20_frankenmask/frankenbeam_mask_rfifind.mask")

opts,args = parser.parse_args()




def a_to_pdot(P_s, acc_ms2):
    LIGHT_SPEED = 2.99792458e8                 # Speed of Light in SI
    return P_s * acc_ms2 /LIGHT_SPEED


def period_modified(p0,pdot,no_of_samples,tsamp,fft_size):
    if (fft_size==0.0):
        return p0 - pdot*float(1<<(no_of_samples.bit_length()-1))*tsamp/2
    else:
        return p0 - pdot*float(fft_size)*tsamp/2

def sortKeyFunc(s):
    return int(re.search('f1_(.*)_',s).group(1).split('_')[0])


def get_params_from_csv_and_fold(opts):
    csv_file = opts.csv_path

    # Read file
    all_data = pd.read_csv(opts.csv_path)

    cand_periods = all_data['period'].to_numpy()
    cand_accs = all_data['acc'].to_numpy()
    cand_dms = all_data['dm'].to_numpy()
    cand_ids = all_data['cand_id_in_file'].to_numpy()
    filenames = all_data['file'].to_numpy()

    sample_xml_file = filenames[0]
    tree = ET.parse(sample_xml_file)
    root = tree.getroot()
  
    tsamp = float(root.find("header_parameters/tsamp").text)
    fft_size = float(root.find('search_parameters/size').text)
    no_of_samples = int(root.find("header_parameters/nsamples").text)
  
    mod_periods=[]
    pdots= []
    for i in range(len(cand_periods)):
        Pdot = a_to_pdot(cand_periods[i],cand_accs[i])
        mod_periods.append(period_modified(cand_periods[i],Pdot,no_of_samples,tsamp,fft_size))
        pdots.append(Pdot) 

    cand_mod_periods = np.asarray(mod_periods,dtype=float)


    no_of_cands = len(cand_mod_periods)
    batch_no = opts.batch_no
    
    # Fold all csv file candidates in batches
    extra = no_of_cands%batch_no
    batches = int(no_of_cands/batch_no) +1
    for x in range(batches):
       start = x*batch_no
       if(x==batches-1):
           end = x*batch_no+extra
       else:
           end = (x+1)*batch_no
       for i in range(start,end):
           folding_packet={}
           folding_packet['period'] = cand_mod_periods[i]
           folding_packet['acc'] = cand_accs[i]
           folding_packet['pdot'] = pdots[i]
           folding_packet['dm'] = cand_dms[i]
           output_name= "%s_cand_no_%03d_dm_%.2f_acc_%.2f"%(os.path.basename(opts.csv_path)[:-4],cand_ids[i],folding_packet['dm'],folding_packet['acc'])
           output_path = filenames[i].split('overview.xml')[0][:-1]

           # Get input filterbank
           tree = ET.parse(filenames[i])
           root = tree.getroot()
           input_name = root.find('search_parameters/infilename').text
           

           script = "prepfold -ncpus 15 -nsub 64 -mask %s -noxwin -topo -p %s -pd %s -dm %s %s -o %s"%(opts.mask_path,str(folding_packet['period']),str(folding_packet['pdot']),str(folding_packet['dm']),input_name,output_name)
           # Wrap a message with output path and command
           info={}
           info['output_path'] = output_path
           info['script'] = script 
 
           # Publish to Rabbit
           #print(info)
           log.debug("Ready to publish to RabbitMQ")
           pika_process.publish_info(opts,info)
 


if __name__ == "__main__":
    
    get_params_from_csv_and_fold(opts)
    

