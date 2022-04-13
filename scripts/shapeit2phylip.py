#!/usr/bin/python

'''
convert Shapeit phased genotype ".haps" file to .geno format file.
A.A. Comeault
16 July 2018

shapeit2phylip.py

eg.
python shapeit2phylip.py /proj/matutelb/projects/drosophila/zaprionus/shapeit/ zindianus-clade_50x_scaffold9_HaplotypeData.sample zindianus-clade_50x_scaffold9_HaplotypeData.haps 250
'''

########################################################
# ---- python packages and user-defined arguments ---- #
########################################################
import sys, glob, fileinput, re, math, os
from pprint import pprint
from subprocess import call
from collections import defaultdict
import numpy as np
from itertools import imap

# ---- parse arguments ---- #
if len(sys.argv[1:]) == 4:
    vals = sys.argv[1:]
    print "arguments successfully loaded!"
    print "path:   ", vals[0]
    print "samples:   ", vals[1]
    print "haplo:   ", vals[2]
    print "n_snps:   ", vals[3]
else:
    print "\nUsage: <shapeit2phylip.py> [/path/to/data/] [sample file; eg. shapeit.sample] [haplotype file; eg. shapeit.haps] [number of snps per window; eg. 250]" 


my_shapit_dir  = str(vals[0])
my_sample_file = str(vals[1]) 
my_phase_dat   = str(vals[2])
n_snps         = int(vals[3])

#my_shapit_dir  = '/proj/matutelb/projects/drosophila/zaprionus/shapeit/'
#my_sample_file = 'zindianus-clade_50x_scaffold9_HaplotypeData.sample'
#my_phase_dat   = 'zindianus-clade_50x_scaffold9_HaplotypeData.haps'
#n_snps = 250 



# ---- get sample information ---- #
my_ind_dat = open( my_shapit_dir + my_sample_file )
my_ind_ids = []

for i, line in enumerate(my_ind_dat):
    if i > 0 and not re.match('^0 0', line):
        line = line.rstrip('\n')        
        line = re.split(r' ', line)
        my_ind_ids.append(line[0])
        


# ---- prepare haplotype dictionary ---- #

ind_haplos = [val for val in my_ind_ids for __ in (0,1)]

for index, object in enumerate(ind_haplos):
	if index % 2 == 0:
		ind_haplos[index] = ind_haplos[index]+'-1'
	elif index % 2 == 1:
		ind_haplos[index] = ind_haplos[index]+'-2'

my_seq_dic = {}
for haplotype in ind_haplos:
	my_seq_dic[haplotype] = ''


# ---- open and write header to window summary '.data.tsv' file ---- #
header = ['scaffold', 'start', 'end', 'mid', 'sites']
data_sum = open(my_shapit_dir + re.split('\.',my_sample_file)[0] +'.' + str(n_snps) +'winds.data.tsv', "a")
data_sum.write('\t'.join(header) +'\n')


       
#############################################################################################
# ---- loop through .haplo file and convert numeric-coded haplotypes into nucleotides  ---- #
#############################################################################################
snp_positions = []
phase_dat = open( my_shapit_dir + my_phase_dat )

for i, line in enumerate(phase_dat):
	line = line.rstrip('\n')
	line = re.split(r' ', line)
	
	my_map_info = [ line[0], line[2], line[3], line[4] ]
	my_haplo_genos = [my_map_info[2:4][geno] for geno in map(int, line[5:])]
	
	snp_positions.append(int(my_map_info[1]))
	
	for id_index, ind in enumerate(ind_haplos):
		my_seq_dic[ind]=  my_seq_dic[ind]+my_haplo_genos[id_index]
	
	if len(snp_positions) > 0 and len(snp_positions) % n_snps == 0:
		wind_summary = [np.min(snp_positions), np.max(snp_positions), np.min(snp_positions)+((np.max(snp_positions)-np.min(snp_positions)) / 2), len(snp_positions)]
		data_sum.write( my_map_info[0] + '\t' + '\t'.join(map(str, wind_summary)) +'\n')
		
		snp_positions = []

# write data for the last window:
if len(snp_positions) > 0:
        wind_summary = [np.min(snp_positions), np.max(snp_positions), np.min(snp_positions)+((np.max(snp_positions)-np.min(snp_positions)) / 2), len(snp_positions)]
        data_sum.write( my_map_info[0] + '\t' + '\t'.join(map(str, wind_summary)) +'\n')
        data_sum.close()
else:
        data_sum.close()



####################################################################################
# ---- write .phylp alignments for widows containing 'n_snps' number of snps  ---- #
####################################################################################

seq_len  = len( my_seq_dic[ind_haplos[0]] )
wind_num = range(seq_len/n_snps + 1)              

for window in wind_num:
	m_seq_len = len(my_seq_dic[ind][window*n_snps:(window+1)*n_snps])
	
	out = open(my_shapit_dir + re.split('\.',my_sample_file)[0] + '_w' + str(window) + '.phylp', "a")
	out.write(' '+ str( len(ind_haplos) ) + ' ' + str( m_seq_len ) +'\n')
	
	
	for id_index, ind in enumerate(ind_haplos):
		out.write( "{:<25}".format(ind) + '\t' + my_seq_dic[ind][window*n_snps:(window+1)*n_snps] + '\n' )
	
	
	# close output file:
	out.close()
