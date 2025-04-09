import argparse
import os
import subprocess
from datetime import datetime


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("sample")

    return parser.parse_args()


def main():
  args = parse_args()
  sample = args.sample  
  
  chroms = ['chr_1', 'chr_2', 'chr_3', 'chr_4']
  fasta_file = '/work/users/d/t/dturissi/lutzomyia/msmc2/genomes/Lutzomyia_longipalpis_male_1.0.fa'
  mask_file = '/work/users/d/t/dturissi/lutzomyia/msmc2/genmap_map/genmap_map.bedgraph'
  
  bam_dir = '/work/users/d/t/dturissi/lutzomyia/msmc2/bams'
  bed_dir = '/work/users/d/t/dturissi/lutzomyia/msmc2/beds'
  vcf_dir = '/work/users/d/t/dturissi/lutzomyia/msmc2/vcfs'
  input_dir = '/work/users/d/t/dturissi/lutzomyia/msmc2/inputs'

  bam_file = os.path.join(bam_dir, sample + '.dedupd.bam')
  
    
#index bam
  os.system(f"samtools index {bam_file}")



  #process chromosomes separately
  for chrom in chroms:
    print(chrom, datetime.now())   
    bed_file = os.path.join(bed_dir, sample + '_' + chrom + '.bed.gz')
    vcf_file = os.path.join(vcf_dir, sample + '_' + chrom + '.vcf.gz')
    input_file = os.path.join(input_dir, sample + '_' + chrom + '_multihetsep.txt')
    
    #get mean depth
    depth = subprocess.run(f"samtools depth -r {chrom} {bam_file} | awk '{{sum += $3}} END {{print sum / NR}}'", shell=True, capture_output=True, text=True).stdout.strip()    
    
    #make bed and vcf
    os.system(f"bcftools mpileup -q 20 -Q 20 -C 50 -r {chrom} -f {fasta_file} {bam_file} | bcftools call -c -V indels | /proj/matutelb/software/msmc-tools-master/bamCaller.py {depth} {bed_file} | gzip -c > {vcf_file}")
    
    #prepare input file
    os.system(f"/proj/matutelb/software/msmc-tools-master/generate_multihetsep.py --mask={bed_file} --mask={mask_file} {vcf_file} > {input_file}")
    

if __name__ == '__main__':
  main()
