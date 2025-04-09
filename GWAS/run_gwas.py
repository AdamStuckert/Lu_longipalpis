import argparse
import json                                                   
import os
import sys
import re
from datetime import datetime

def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("input_json_file", help="json file containing gwas parameters. Example file: /proj/matutelb/projects/gwas/gwas_template.json")

    return parser.parse_args()


def main():
  #process json
  args = parse_args()
  input_json_file = args.input_json_file  
  jf = open(input_json_file)
  gwas_info = json.load(jf)
  
  
  #define dirs and files
  base_dir = '/proj/matutelb/projects/gwas'
  scripts_dir = os.path.join(base_dir, 'scripts')
  gwas_dir = os.path.join(base_dir, 'gwas_results', gwas_info['gwas_name'])
  input_dir = os.path.join(base_dir, 'gwas_results', gwas_info['gwas_name'], 'input')
  gwas_sh_file = os.path.join(input_dir, gwas_info['gwas_name'] + '_gwas_plink.sh')
  
  snp_report_suffix = gwas_info['genotype_dataset']    
  snp_report_file = os.path.join(scripts_dir, 'gwas_snp_reports_' + snp_report_suffix + '.R')
  
  
  #prep gwas
  print('prepping gwas', datetime.now(), "\n")
  prep_gwas(gwas_info, input_json_file, base_dir)

  #run gwas  
  print('running gwas', datetime.now(), "\n")
  os.system(f"{gwas_sh_file}")

  #gwas plots  
  print("\n", 'making gwas plots', datetime.now(), "\n")
  os.system(f"Rscript {base_dir}/scripts/gwas_plots.R {gwas_info['gwas_name']}")

  print('done', datetime.now(), "\n") 
  print('snp reports can be generated after determining an appropriate p cutoff:', "\n", f"Rscript {snp_report_file} {gwas_info['gwas_name']} 1e-5")
  


def prep_gwas(gwas_info, input_json_file, base_dir):
  #define dirs and files
  gwas_dir = os.path.join(base_dir, 'gwas_results', gwas_info['gwas_name'])
  input_dir = os.path.join(base_dir, 'gwas_results', gwas_info['gwas_name'], 'input')
  results_dir = os.path.join(base_dir, 'gwas_results', gwas_info['gwas_name'], 'results')
  
  plink_stub = os.path.join(base_dir, 'genotype_datasets', gwas_info['genotype_dataset'], gwas_info['genotype_dataset'])

  gwas_json_file = os.path.join(input_dir, 'gwas_' + gwas_info['gwas_name'] + '.json')
  gwas_covar_file = os.path.join(input_dir, gwas_info['gwas_name'] + '_covar.txt')
  gwas_phenotype_file = os.path.join(input_dir, gwas_info['gwas_name'] + '_phenotype.txt')
  unformatted_phenotype_file = os.path.join(input_dir, gwas_info['gwas_name'] + '_phenotype_unformatted.txt')

  run_gwas_sh_file = os.path.join(input_dir, gwas_info['gwas_name'] + '_gwas_plink.sh')
  gwas_out_stub = os.path.join(results_dir, gwas_info['gwas_name'] + '_gwas')
  gwas_results_file = os.path.join(results_dir, gwas_info['gwas_name'] + '_gwas.' + gwas_info['phenotype_name'] + '.glm.linear')
  gwas_results_file_raw = gwas_results_file + '.raw'
  gwas_results_logistic_file = os.path.join(results_dir, gwas_info['gwas_name'] + '_gwas.' + gwas_info['phenotype_name'] + '.glm.logistic')
  gwas_results_logistic_file_raw = gwas_results_logistic_file + '.raw'
  
  
  #check if gwas_dir already exists
  if os.path.isdir(gwas_dir):
    sys.exit(f"ERROR: {gwas_dir} already exists. Either remove that directory or specifiy a new gwas name.")
  
  
  #ensure phenotype_name is alphanumeric (underscroes ok)
  if not re.match(r'^\w+$', gwas_info['phenotype_name']):
    sys.exit(f"ERROR: phenotype_name '{gwas_info['phenotype_name']}' contains spaces or punctuation. Please use only alphanumerics or underscores.") 
  
  
  #make gwas dirs
  os.system(f"mkdir {gwas_dir}")
  os.system(f"mkdir {input_dir}")
  os.system(f"mkdir {results_dir}")
  
  
  #archive copy of json
  os.system(f"cp {input_json_file} {gwas_json_file}")
  
  
  #define covar param and archive file
  covar_param = ""
  if gwas_info['covariate_file'] != '':
    os.system(f"cp {gwas_info['covariate_file']} {gwas_covar_file}")
    covar_param = '--covar ' + gwas_covar_file
  
  
  #archive copy of phenotype_file
  os.system(f"cp {gwas_info['phenotype_file']} {unformatted_phenotype_file}")
  
  
  #process phenotype file
  with open(gwas_phenotype_file, 'w') as po:  
    with open(unformatted_phenotype_file, 'r') as pi: 
      header = pi.readline().strip()
      (id_name, old_phenotype_name) = header.split('\t')      
      po.write(f"#FID\tIID\t{gwas_info['phenotype_name']}\n")

      for line in pi.readlines():
        line = line.strip()
        (id, phenotype) = line.split('\t')
        po.write(f"0\t{id}\t{phenotype}\n")
  

  #make gwas plink sh file
  with open(run_gwas_sh_file, 'w') as gp:  
    gp.write(f"#!/bin/bash\n")
    gp.write(f"plink2 -bfile {plink_stub} --glm {gwas_info['plink_parameters']} {covar_param} --pheno {gwas_phenotype_file} --out {gwas_out_stub}\n")

    #filter gwas results if covariates were used so R can load a smaller dataframe
    if gwas_info['covariate_file'] != '':
      gp.write(f"""if test -f {gwas_results_logistic_file}; then\n""")
      gp.write(f"""  mv {gwas_results_logistic_file} {gwas_results_logistic_file_raw}\n""")
      gp.write(f"""  cat <(head -1 {gwas_results_logistic_file_raw}) <(awk '$6 == "ADD" {{print $0}}' {gwas_results_logistic_file_raw}) > {gwas_results_logistic_file}\n""")
      gp.write(f"""else\n""")
      gp.write(f"""  mv {gwas_results_file} {gwas_results_file_raw}\n""")
      gp.write(f"""  cat <(head -1 {gwas_results_file_raw}) <(awk '$6 == "ADD" {{print $0}}' {gwas_results_file_raw}) > {gwas_results_file}\n""")
      gp.write(f"""fi\n""")
      
  #make gwas sh file executable
  os.system(f"chmod +x {run_gwas_sh_file}")  
  

if __name__ == '__main__':
  main()
