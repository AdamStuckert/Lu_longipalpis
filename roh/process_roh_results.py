# Import packages.
import sys
import sqlite3
import re
import pandas as pd
import glob


def main():
  #establish db connection and create lutzo_sweepfinder table
  db_file = 'lutzo_roh.db'
  
  conn = sqlite3.connect(db_file)  
  conn.execute(f"attach '/work/users/d/t/dturissi/lutzomyia/lutzomyia_metadata.db' as m")
  

  conn.execute(f"drop table if exists lutzo_roh")
  conn.execute(f"""create table lutzo_roh
                  (lr_id int primary key,
                   min_kb int,
                   pop varchar(20),
                   sample varchar(50),
                   num_roh_regions int,
                   sum_roh_kb decimal(10, 3))""")

  conn.execute(f"create index idx_lr_sample on lutzo_roh(sample)")
  conn.execute(f"create index idx_lr_pop on lutzo_roh(pop)")
  
  
  conn.execute(f"drop table if exists lutzo_roh_regions")
  conn.execute(f"""create table lutzo_roh_regions
                  (lrr_id int primary key,
                   min_kb int,
                   pop varchar(20),
                   sample varchar(50),
                   chrom varchar(20),
                   start int,
                   end int,
                   num_snps int,
                   prop_homo decimal(4,3),
                   prop_het decimal(4,3))""")

  conn.execute(f"create index idx_lrr_sample on lutzo_roh_regions(sample)")
  conn.execute(f"create index idx_lrr_pop on lutzo_roh_regions(pop)")
  conn.execute(f"create index idx_lrr_start on lutzo_roh_regions(start)")
  conn.execute(f"create index idx_lrr_end on lutzo_roh_regions(end)")
  
  
  #define primary key value since annoyingly sqlite doesn't have a functional auto increment for primary keys
  lr_id = 0
  roh_dict = {'lr_id': [],
              'min_kb': [],
              'pop': [],
              'sample': [],
              'num_roh_regions': [],
              'sum_roh_kb': []}

  lrr_id = 0
  roh_region_dict = {'lrr_id': [],
                     'min_kb': [],
                     'pop': [],
                     'sample': [],
                     'chrom': [],
                     'start': [],
                     'end': [],
                     'num_snps': [],
                     'prop_homo': [],
                     'prop_het': []}
  
  
  site_query = conn.execute(f"""select distinct site
                                from lutzomyia_samples""")
                              
                              
  for site, in site_query:      
    pop = site.split(',')[0].replace(' ', '_')
    if pop == 'Ricaurte':
      pop = 'Colombia'
  
    #loop over all results files
    for results_file in glob.glob(f"results/{pop}/*kb.hom.indiv"):
      #get values from name of result file
      min_kb = re.match(r'results/.+_(\d+)kb.hom.indiv', results_file).groups()[0]    
      print(results_file)
    
    
      #process results file                      
      with open(results_file, 'r') as v: 
        next(v)
        for line in v:
          lr_id += 1        
          values = line.strip().split()
          
          roh_dict['lr_id'].append(lr_id)
          roh_dict['min_kb'].append(min_kb)
          roh_dict['pop'].append(pop)
          roh_dict['sample'].append(values[1])
          roh_dict['num_roh_regions'].append(values[3])
          roh_dict['sum_roh_kb'].append(values[4])
              
    for region_file in glob.glob(f"results/{pop}/*kb.hom"):
      #get values from name of result file
      min_kb = re.match(r'results/.+_(\d+)kb.hom', region_file).groups()[0]    
      print(region_file)
    
      #process results file                      
      with open(region_file, 'r') as v: 
        next(v)
        for line in v:
          lrr_id += 1        
          values = line.strip().split()
          
          chrom = values[3][4:]
          
          roh_region_dict['lrr_id'].append(lrr_id)
          roh_region_dict['min_kb'].append(min_kb)
          roh_region_dict['pop'].append(pop)
          roh_region_dict['sample'].append(values[1])
          roh_region_dict['chrom'].append(chrom)
          roh_region_dict['start'].append(values[6])
          roh_region_dict['end'].append(values[7])
          roh_region_dict['num_snps'].append(values[9])
          roh_region_dict['prop_homo'].append(values[11])
          roh_region_dict['prop_het'].append(values[12])


 
              
#directly insert pandas dataframe into db table          
  roh_df = pd.DataFrame(roh_dict)  
  roh_df.to_sql('lutzo_roh', if_exists = 'append', index=False, con=conn)    
  
  roh_region_df = pd.DataFrame(roh_region_dict)  
  roh_region_df.to_sql('lutzo_roh_regions', if_exists = 'append', index=False, con=conn)    
  


if __name__ == '__main__':
  main()
