#calculate runs of homozygosity

python3 lutzo_roh.py


#process roh results and save in database

python3 process_roh_results.py

#generate figure for manuscript

Rscript lutzo_roh_manuscript.R
