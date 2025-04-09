import argparse
import sqlite3
import os
import sys
import tempfile


#python3 /proj/matutelb/projects/gwas/scripts/make_anno_db_lutzomyia.py /proj/matutelb/data_share/longipalpis/Lu_long_male.gff3 /proj/matutelb/projects/gwas/genotype_datasets/lutzomyia_male/lutzomyia_male_anno.db

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("gff_file")
    parser.add_argument("db_file")
 
    return parser.parse_args()


def main():
  args = parse_args()
  gff_file = args.gff_file  
  db_file = args.db_file  

  conn = sqlite3.connect(db_file)
    
  cur = conn.cursor()                            


  cur.execute("drop table if exists gene_info")
  cur.execute("""create table gene_info
                 (gene_id varchar(30) primary key,
                  chrom varchar(20),
                  start int,
                  end int,
                  strand varchar(1))""")
   
  cur.execute("create index idx_gi_se on gene_info(start,end)")
 
  cur.execute("drop table if exists transcript_info")
  cur.execute("""create table transcript_info
                 (transcript_id varchar(30) primary key,
                  gene_id varchar(30))""")

  cur.execute("create index idx_ti_ngi on transcript_info(gene_id)")


  cur.execute("drop table if exists cds_info")
  cur.execute("""create table cds_info
                 (cds_id varchar(20) primary key,
                  transcript_id varchar(30),
                  chrom varchar(20),
                  start int,
                  end int,
                  strand varchar(1),
                  phase int)""")


  cur.execute("create index idx_ci_ti on cds_info(transcript_id)")
  cur.execute("create index idx_ci_se on cds_info(start,end)")


  cur.execute("drop table if exists intron_info")
  cur.execute("""create table intron_info
                 (intron_id int primary key,
                  transcript_id varchar(30),
                  chrom varchar(20),
                  start int,
                  end int,
                  strand varchar(1))""")

  cur.execute("create index idx_ii_ti on intron_info(transcript_id)")
  cur.execute("create index idx_ii_se on intron_info(start,end)")

  cur.execute("drop table if exists utr_info")
  cur.execute("""create table utr_info
                 (utr_id int primary key,
                  transcript_id varchar(30),
                  chrom varchar(20),
                  start int,
                  end int,
                  strand varchar(1),
                  utr_type varchar(6))""")

  cur.execute("create index idx_ui_ti on utr_info(transcript_id)")
  cur.execute("create index idx_ui_se on utr_info(start,end)")

  conn.close()

  
  counter = 0
  intron_id = 0
  
  with tempfile.NamedTemporaryFile(mode='w') as gene_file:
    with tempfile.NamedTemporaryFile(mode='w') as tran_file:
      with tempfile.NamedTemporaryFile(mode='w') as cds_file:
        with tempfile.NamedTemporaryFile(mode='w') as intron_file:
          with open(gff_file, 'r') as gff: 
            for line in gff:
              counter += 1
              if counter % 100000 == 0:
                print(counter)
              
              if line[0] != "#":
                line = line.strip()
                gff_vals = line.split('\t')
                
                if len(gff_vals) == 9:                 
                  (chrom, source, feature_type, start, end, score, strand, phase, attributes) = gff_vals
                  
                  if chrom[:3] == 'chr':
                    chrom = chrom[3:]
                  
                  if feature_type == 'gene': 
                    attribute_dict = process_attributes(attributes)
                    gene_file.write(f"{attribute_dict['ID']}\t{chrom}\t{start}\t{end}\t{strand}\n")                  
                  elif feature_type == 'mRNA': 
                    attribute_dict = process_attributes(attributes)
                    tran_file.write(f"{attribute_dict['ID']}\t{attribute_dict['Parent']}\n")                  
                  elif feature_type == 'CDS':
                    attribute_dict = process_attributes(attributes)
                    for parent in attribute_dict['Parent'].split(','):
                      cds_file.write(f"{attribute_dict['ID']}\t{attribute_dict['Parent']}\t{chrom}\t{start}\t{end}\t{strand}\t{phase}\n")                  
                  elif feature_type == 'intron':
                    attribute_dict = process_attributes(attributes)
                    for parent in attribute_dict['Parent'].split(','):
                      intron_file.write(f"{attribute_dict['ID']}\t{attribute_dict['Parent']}\t{chrom}\t{start}\t{end}\t{strand}\n")                  
          
          gene_file.flush() 
          tran_file.flush() 
          cds_file.flush() 
          intron_file.flush() 
          
          os.system(f"""sqlite3 {db_file} ".mode tabs" ".import {gene_file.name} gene_info" """)
          os.system(f"""sqlite3 {db_file} ".mode tabs" ".import {tran_file.name} transcript_info" """)
          os.system(f"""sqlite3 {db_file} ".mode tabs" ".import {cds_file.name} cds_info" """)
          os.system(f"""sqlite3 {db_file} ".mode tabs" ".import {intron_file.name} intron_info" """)

  
  #deduce utrs
  conn = sqlite3.connect(db_file)
  t_cur = conn.cursor()                            
  cds_cur = conn.cursor()                            
  transcript_query = t_cur.execute(f"""select transcript_id, start, end, strand
                                       from transcript_info t, gene_info g
                                       where t.gene_id = g.gene_id""")
  
  utr_id = 0                            
  with tempfile.NamedTemporaryFile(mode='w') as utr_file:
    for (transcript_id, gene_start, gene_end, strand) in transcript_query:
      cds_query = cds_cur.execute(f"""select chrom, start, end
                                      from cds_info
                                      where transcript_id = '{transcript_id}'
                                      order by start""")
      cds_count = 0
      last_cds_end = 0
      for (chrom, cds_start, cds_end) in cds_query:
        cds_count += 1
        if cds_count == 1:
          if cds_start != gene_start:
            utr_start = gene_start
            utr_end = cds_start - 1

            if strand == '+':
              utr_type = '5prime'
            else:
              utr_type = '3prime'
                            
            utr_id += 1
            utr_file.write(f"{utr_id}\t{transcript_id}\t{chrom}\t{utr_start}\t{utr_end}\t{strand}\t{utr_type}\n")  
          
        last_cds_end = cds_end
      
      if last_cds_end != gene_end:
          utr_start = last_cds_end + 1
          utr_end = gene_end
          
          if strand == '+':
            utr_type = '5prime'
          else:
            utr_type = '3prime'
            
          utr_id += 1
          utr_file.write(f"{utr_id}\t{transcript_id}\t{chrom}\t{utr_start}\t{utr_end}\t{strand}\t{utr_type}\n")  
    
    conn.close()
    utr_file.flush() 
    os.system(f"""sqlite3 {db_file} ".mode tabs" ".import {utr_file.name} utr_info" """)



def process_attributes(attributes):
  attributes_dict = {}
  for attribute in attributes.split(';'):
    if attribute != '':
      (key, val) = attribute.split('=')
      attributes_dict[key] = val

  return attributes_dict  
  



if __name__ == '__main__':
  main()
