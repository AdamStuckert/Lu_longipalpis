# Import packages.
import sys
import sqlite3
import re
import os
import glob
import tempfile


def main():
  db_file = '/work/users/d/t/dturissi/lutzomyia/poly_windows/lutzomyia_poly_win.db'
  pi_table = 'poly_win_pops_sobral_1s_subpop_10000_vcftools'  
  
  #establish db connection and create lutzo_sweepfinder_sobral_1s_subpop table
  conn = sqlite3.connect(db_file)  

  conn.execute(f"drop table if exists {pi_table}")
  conn.execute(f"""create table {pi_table}
                  (pwpv_id int primary key,
                   pop varchar(20),
                   chrom varchar(20),
                   start int,
                   end int,
                   num_sites int,
                   pi float)""")


  conn.execute(f"create index idx_pwpv_pop on poly_win_pops_sobral_1s_subpop_10000_vcftools(pop)")
  conn.execute(f"create index idx_pwpv_start_end on poly_win_pops_sobral_1s_subpop_10000_vcftools(start, end)")
  conn.execute(f"create index idx_pwpv_end on poly_win_pops_sobral_1s_subpop_10000_vcftools(end)")
  
  conn.close()

  results_dir = 'pi_vcftools_results'
  results_file_str = os.path.join(results_dir, '*.pi')
  
  pwpv_id = 0
  with tempfile.NamedTemporaryFile(mode='w') as t:
    for results_file in glob.glob(results_file_str): 
      m = re.search(r'.+/pi_vcftools_10000_(.+)_chr_\d+.windowed.pi', results_file)
      pop = m.group(1)
      
      with open(results_file, 'r') as r: 
        next(r)
        for line in r:
          pwpv_id += 1
          line = line.strip()

          (chrom, start, end, num_sites, pi) = line.split('\t')

          t.write(f"{pwpv_id}\t{pop}\t{chrom}\t{start}\t{end}\t{num_sites}\t{pi}\n")   
                         
    t.flush()      
    os.system(f"""sqlite3 {db_file} ".mode tabs" ".import {t.name} {pi_table}" """)




if __name__ == '__main__':
  main()
