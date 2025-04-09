import argparse
import os
import subprocess
from datetime import datetime

#python3 /work/users/d/t/dturissi/lutzomyia/msmc2/run_msmc2.py JB_LF11-lutzomyia_female_78_021921_TAGGCATG-TCGACTAG_S2


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("sample")
    parser.add_argument("time_pattern")

    return parser.parse_args()


def main():
  args = parse_args()
  sample = args.sample  
  time_pattern = args.time_pattern  
  
  time_pattern_formatted = time_pattern.replace('*', '-')
  
  chroms = ['chr_1', 'chr_2', 'chr_3', 'chr_4']  
  input_dir = '/work/users/d/t/dturissi/lutzomyia/msmc2/inputs'
  results_dir = '/work/users/d/t/dturissi/lutzomyia/msmc2/results'


  #process chromosomes separately
  input_file_str = ''
  for chrom in chroms:  
    input_file = os.path.join(input_dir, sample + '_' + chrom + '_multihetsep.txt')
    input_file_str += os.path.join(input_dir, sample + '_' + chrom + '_multihetsep.txt') + ' '
      
  #run msmc2
  results_file = os.path.join(results_dir, sample + '_genome_' + time_pattern_formatted + '.msmc2')
  os.system(f"msmc2 -t 1 -p {time_pattern} -o {results_file} {input_file_str}")
    

if __name__ == '__main__':
  main()
