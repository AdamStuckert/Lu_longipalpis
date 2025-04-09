import sqlite3
import os
import tempfile
import glob
import re


#python3 load_msmc2_results.py 



def main():
  db_file = 'msmc2_results.db'
  results_dir = '/work/users/d/t/dturissi/lutzomyia/msmc2/results'
  results_file_str = os.path.join(results_dir, '*.msmc2.final.txt')
  
  conn = sqlite3.connect(db_file)
    
  cur = conn.cursor()                            

  cur.execute("drop table if exists msmc2_results")
  cur.execute("""create table msmc2_results
                 (mr_id int auto_increment primary key,
                  pattern varchar(50),
                  sample varchar(50),
                  chrom varchar(10),
                  time_index int,
                  left_time float,
                  right_time float,
                  lambda float)""")
   
  cur.execute("create index idx_mr_sample on msmc2_results(sample)")
  cur.execute("create index idx_mr_chrom on msmc2_results(chrom)")
  cur.execute("create index idx_mr_time on msmc2_results(left_time, right_time)")
  conn.close()

  
  mr_id = 0
  with tempfile.NamedTemporaryFile(mode='w') as t:
    for results_file in glob.glob(results_file_str): 
      m = re.search(r'.+/(.+)_(chr_.+|genome)_(.+).msmc2.final.txt', results_file)
      sample, chrom, pattern = m.group(1), m.group(2), m.group(3)
      
      if chrom[:4] == 'chr_':
        chrom = chrom[4]
      
      with open(results_file, 'r') as r: 
        next(r)
        for line in r:
          mr_id += 1
          line = line.strip()
          (time_index, left_time, right_time, lamb) = line.split('\t')

          if right_time == 'inf':
            right_time = 0
            
          t.write(f"{mr_id}\t{pattern}\t{sample}\t{chrom}\t{time_index}\t{left_time}\t{right_time}\t{lamb}\n")   
                         
    t.flush()      
    os.system(f"""sqlite3 {db_file} ".mode tabs" ".import {t.name} msmc2_results" """)




if __name__ == '__main__':
  main()
