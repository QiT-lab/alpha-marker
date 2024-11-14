import sys
import os
import csv
import ast
import pandas as pd
import numpy as np
from scipy.stats import ranksums
#from scipy.stats import wilcoxon
import multiprocessing
from multiprocessing import Process, Queue, Pool, Manager



def merge_alpha_files(alpha_path):
  
  files = os.listdir(alpha_path,)
  all_df = pd.read_csv(os.path.join(alpha_path, files[0]))
  name = files[0].split('.')[0]
  all_df = all_df.rename(columns={'alpha_mean':name})
  
  for f in files[1:]:
    if os.path.isfile(os.path.join(alpha_path, f)):
      df = pd.read_csv(os.path.join(alpha_path, f))
      all_df = pd.merge(all_df, df, on = 'cpg_sites')
      name = f.split('.')[0]
      all_df = all_df.rename(columns={'alpha_mean':name})
      
  return all_df



def split_df(all_df, process_number):
  
    total_lines = len(all_df)
    lines_per_part = total_lines // process_number

    parts = []
    for i in range(process_number):
        start_idx = i * lines_per_part
        end_idx = (i + 1) * lines_per_part if i != process_number - 1 else total_lines
        #part = start_idx + ':' + end_idx
        part = [start_idx,end_idx]
        parts.append(part)
    return parts




def test_mean_alpha(target,df,result):
  target_columns = [col for col in df.columns if target in col]
  tg = df[target_columns]
  bg = df.drop(target_columns, axis=1)
  
  c_df = df.copy()
  c_df = c_df.assign(delta_mean=['nan'] * len(c_df))
  c_df = c_df.assign(p_value=['nan'] * len(c_df))
  for i in range(len(c_df)):
    A = tg.iloc[i,:]
    B = bg.iloc[i,1:]
    #A = list(A)
    #B = list(B)
    stat, p_value = ranksums(A, B)
    #stat, p_value = wilcoxon(A, B)
    mean = np.mean(A) - np.mean(B)
    c_df.iloc[i,-2] = mean
    c_df.iloc[i,-1] = p_value
  df = c_df.copy()
  result.append(df)
    


if __name__ == "__main__":
  
  
  alpha_path = sys.argv[1]
  marker_path = sys.argv[2]
  threads = int(sys.argv[3])
  all_type = ast.literal_eval(sys.argv[4])
  
  all_df = merge_alpha_files(alpha_path)
  #all_df.replace([np.inf, -np.inf], np.nan, inplace=True)
  #all_df = all_df.dropna()
  
  #all_type = ['Vein-Endothel','Hepatocytes','Erythrocyte','Monocytes','Granulocytes']
  #all_type = ['cfDNA','Breast']
  parts = split_df(all_df,threads)
  

  for s in all_type:
    print(s)
    processes = []
    manager = multiprocessing.Manager()
    result = manager.list()
    for j in range(threads):
       part = parts[j]
       df = all_df.iloc[part[0]:part[1],:]
       #df = test_mean_alpha(target,df)
       process = multiprocessing.Process(target=test_mean_alpha, args=(s, df, result))
       processes.append(process)
       process.start()
       
    for process in processes:
      process.join()
    
    diff_df = pd.concat(result, ignore_index=True)
    s = s + '_marker.txt'
    diff_df.to_csv(os.path.join(marker_path, s),sep='\t',index=False)
    
  

      
      
  
