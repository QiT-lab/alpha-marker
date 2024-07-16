import csv
import numpy as np
import pandas as pd
import sys
import math
import collections
import multiprocessing
from multiprocessing import Process, Queue, Pool, Manager
import time, random, cv2

#################################################



def split_file(file_path, process_number):
    with open(file_path, 'r') as file:
        lines = file.readlines()

    total_lines = len(lines)
    lines_per_part = total_lines // process_number

    parts = []
    for i in range(process_number):
        start_idx = i * lines_per_part
        end_idx = (i + 1) * lines_per_part - 1 if i != process_number - 1 else total_lines
        #part = lines[start_idx:end_idx]
        part = [start_idx,end_idx]
        parts.append(part)
    return parts


def get_region_dict(file,part):
    """
    retrieve the block regions from text file
    """
    regions_dict = collections.defaultdict(list)  

    with open(file, "r") as input_file:
        regions_file = csv.reader(input_file, delimiter="\t")
        for i, line in enumerate(regions_file):
           if i >= part[0] and i <= part[1]:
             chromo, start, end, CpGstart, CpGend = line[0], line[1], line[2], int(line[3]), int(line[4])
             regions_dict[(CpGstart,CpGend)] = [0]
    return regions_dict

 
def in_intervals(num, intervals):
    for interval in intervals:
        if num >= interval[0] and num <= interval[1]:
            return interval
    return None


def get_mean_alpha(regions_dict):
  df = pd.DataFrame(columns=['cpg_sites', 'alpha_mean'])
  i=0
  for key, value in regions_dict.items():
    if value != [0] and len(value) > 4:
      v = value[1:]
      v_array = np.array(v)
      alpha_mean = np.mean(v_array)
      df.loc[i] = [str(key),alpha_mean]
      i=i+1
  return df


def get_methylation_counts(cpg_file, regions_file,part,result):
    """
    calculate the mean methylation values for all reads in the selected region
    """
    regions_dict = get_region_dict(regions_file,part) 
    intervals = list(regions_dict.keys())
    start_index = intervals[0][0]
    end_index = intervals[len(intervals)-1][1]
  
    with open(cpg_file, "r") as input_file:
        cpg = csv.reader(input_file, delimiter="\t")
        for line in cpg:
            chromo, start, number, alpha = line[0], int(line[1]), int(line[3]), float(line[4])
            if start >= start_index and start <= end_index:
              key = in_intervals(start, intervals)
              if key != None:
                regions_dict[key] = regions_dict[key] + [alpha]*number
    df = get_mean_alpha(regions_dict)
    result.append(df)
    #return df




if __name__ == "__main__":

    regions_file = sys.argv[1]   # the block file
    alpha_cpg_file = sys.argv[2] # the input file
    threads = int(sys.argv[3])
    output_file_name = sys.argv[4]
    
    #regions_file='/mnt/data3/qiting/cfDNA/methylation/cell_type_marker/uxm_marker/plasma_marker/marker/block.bed'
    #alpha_cpg_file = '/mnt/data3/qiting/cfDNA/methylation/cell_type_marker/own_method/GSM5652315_Blood-Granulocytes_alpha.bed'

    #block = pd.read_csv(regions_file,sep='\t',header=None)
    parts = split_file(regions_file,threads)
    #chr = block[0].unique()
    #df = get_methylation_counts(alpha_cpg_file, regions_file,parts[0])  # get methylation read counts of Cpgs within region
    
    time_start = time.time()
    processes = []
    manager = multiprocessing.Manager()
    result = manager.list()
    for i in range(threads):
       # queue = multiprocessing.Queue()
       # queues.append(queue)
       part = parts[i]
       #df = get_methylation_counts(alpha_cpg_file, regions_file,part,result)
       #process = multiprocessing.Process(target=lambda queue: queue.put(get_methylation_counts(alpha_cpg_file, regions_file,part)), args=(queue,))
       process = multiprocessing.Process(target=get_methylation_counts, args=(alpha_cpg_file, regions_file,part,result))
       processes.append(process)
       process.start()
       
    for process in processes:
      process.join()
    
    time_end = time.time()-time_start
    print(time_end)
    df = pd.concat(result, ignore_index=True)
    df.to_csv(output_file_name, index=False,quoting=None)



