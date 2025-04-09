import sqlite3
import os
import tempfile


#python3 load_lutzomyia_samples.py 



def main():
  db_file = 'lutzomyia_metadata.db'
  metadata_file = 'Lu_longipalpis_metadata.tsv'
  
  conn = sqlite3.connect(db_file)
    

  conn.execute("drop table if exists lutzomyia_samples")
  conn.execute("""create table lutzomyia_samples
                 (sample varchar(50) primary key,
                  site varchar(30),
                  country varchar(30),
                  sex varchar(5),
                  pop varchar(30))""")
   
  conn.execute("create index idx_ls_site on lutzomyia_samples(site)")
  conn.execute("create index idx_ls_pop on lutzomyia_samples(pop)")

  
  conn.execute("drop table if exists lutzomyia_pops")
  conn.execute("""create table lutzomyia_pops
                 (pop varchar(30) primary key,
                  pop_manuscript varchar(30),
                  pop_col varchar(10),
                  country_col varchar(10))""")
   
  
  conn.execute("""insert into lutzomyia_pops
                  values
                  ('Colombia', 'Colombia', '#0000AA', '#0000AA'),
                  ('Jacobina', 'Jacobina', '#0080FF', 'maroon'),
                  ('Jacobina_L_S1', 'Jacobina S1', 'grey', 'maroon'),
                  ('Jacobina_M_S2', 'Jacobina S2', 'tan4', 'maroon'),
                  ('Laphina', 'Laphina', '#00DCDC', 'maroon'),
                  ('Marajo', 'Marajo', '#00AB00', 'maroon'),
                  ('Sobral_1S', 'Sobral 1S', '#FF8000', 'maroon'),
                  ('Sobral_1S_A', 'Sobral 1S A', 'darkgoldenrod1', 'maroon'),
                  ('Sobral_1S_B', 'Sobral 1S B', 'pink', 'maroon'),
                  ('Sobral_2S', 'Sobral 2S', '#AA0000', 'maroon')""")
  

  conn.commit()
  conn.close()

  
  with tempfile.NamedTemporaryFile(mode='w') as t:
    with open(metadata_file, 'r') as r: 
      next(r)
      for line in r:
        line = line.strip()
        (sample, site, sex, path) = line.split('\t')
        place, country = site.split(', ')
        
        pop = site.split(', ')[0].replace(' ', '_')
        if pop == 'Ricaurte':
          pop = 'Colombia'
    
        t.write(f"{sample}\t{site}\t{country}\t{sex}\t{pop}\n")   
                         
    t.flush()      
    os.system(f"""sqlite3 {db_file} ".mode tabs" ".import {t.name} lutzomyia_samples" """)




if __name__ == '__main__':
  main()
