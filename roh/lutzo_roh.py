import sqlite3
import os
import sys


#python3 lutzo_roh 50



def main():  
  min_kb = sys.argv[1]
  
  base_dir = '/work/users/d/t/dturissi/lutzomyia'  
  roh_dir = os.path.join(base_dir, 'roh')
  results_dir = os.path.join(roh_dir, 'results')
    
  db_file = os.path.join(base_dir, 'lutzomyia_metadata.db')
  plink_stub = os.path.join(base_dir, 'lutzo_all_male')
 
  conn = sqlite3.connect(db_file)    
  site_query = conn.execute(f"""select distinct site
                                from lutzomyia_samples""")
                              
                              
  total_cmds = 0
  for site, in site_query:      
    pop = site.split(',')[0].replace(' ', '_')
    if pop == 'Ricaurte':
      pop = 'Colombia'
    
    print(pop)
    
    pop_results_dir = os.path.join(results_dir, pop)
    pop_id_file = os.path.join(pop_results_dir, pop + '_ids.txt')
    
    os.system(f"mkdir -p {pop_results_dir}")
    
    os.system(f"""sqlite3 -separator $'\t' {db_file} "select 0, sample from lutzomyia_samples where site = '{site}'" > {pop_id_file} """)
    os.system(f"""plink --allow-extra-chr --keep {pop_id_file} -bfile {plink_stub} --homozyg-kb {min_kb} --out {pop_results_dir}/{pop}_{min_kb}kb""")
            
       
if __name__ == '__main__':
  main()
