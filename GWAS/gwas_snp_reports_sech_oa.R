library(rjson)
library("RSQLite")  

args <- commandArgs(trailingOnly = TRUE)
gwas_name <- args[1]
p_cutoff <- as.numeric(args[2])


#process json and name files
base_dir <- '/proj/matutelb/projects/gwas'

gwas_json <- paste(base_dir, '/gwas_results/', gwas_name, '/input/gwas_', gwas_name, '.json', sep='')
gwas_info <- fromJSON(file=gwas_json)

anno_db_file <-  paste(base_dir, '/genotype_datasets/', gwas_info$genotype_dataset, '/', gwas_info$genotype_dataset, '_anno.db', sep='')
results_dir <- paste(base_dir, '/gwas_results/', gwas_info$gwas_name, '/results', sep='')
gwas_snp_db <- paste(gwas_info$gwas_name, '_snp.db', sep='')
gwas_results_file <- paste(gwas_info$gwas_name, '_gwas.', gwas_info$phenotype_name, '.glm.linear', sep='')
gwas_logistic_results_file <- paste(gwas_info$gwas_name, '_gwas.', gwas_info$phenotype_name, '.glm.logistic', sep='')
chrom_name_file <- paste(base_dir, '/genotype_datasets/', gwas_info$genotype_dataset, '/', gwas_info$genotype_dataset, '_chrom_names.txt', sep='')


gene_file <- paste(gwas_info$gwas_name, '_snp_genes_', p_cutoff, '.txt', sep='')
go_file <- paste(gwas_info$gwas_name, '_snp_go_', p_cutoff, '.txt', sep='')
intergenic_file <- paste(gwas_info$gwas_name, '_snp_intergenic_', p_cutoff, '.txt', sep='')
seq_type_file <- paste(gwas_info$gwas_name, '_snp_seq_type_', p_cutoff, '.txt', sep='')

setwd(results_dir)

#check for logistic regression results_file
if (file.exists(gwas_logistic_results_file))
  {gwas_results_file = gwas_logistic_results_file}


#load chromosome names
chrom_names <- scan(chrom_name_file, 'character', sep=',')

#establish database connections which also creates a sqlite database to store gwas results
gwas_snp_con <- dbConnect(dbDriver("SQLite"), gwas_snp_db)
dbSendQuery(gwas_snp_con, paste("attach database '", anno_db_file, "' as a", sep=''))

#connect orthodb
orthodb_file <-  '/proj/matutelb/projects/gwas/orthodb/odb11v0_subsetted_sim_mel.db'
mel_species_id <- '7227_0'
dbSendQuery(gwas_snp_con, paste("attach database '", orthodb_file, "' as o", sep=''))
 

#create table to store gwas results for significant snps
dbSendQuery(gwas_snp_con, "drop table if exists gwas_snp")
dbSendQuery(gwas_snp_con, "create table gwas_snp
                           (snp_id varchar(30) primary key,
                            chrom varchar(10),
                            pos int,
                            p float)")


#process gwas results
gwas <- read.table(gwas_results_file, comment.char='', header=T)
names(gwas)[1] <- 'CHROM'


#insert significant snps into db
for (i in which(gwas$P < p_cutoff))
  {    
  if (gwas$ID[i] == '.')
    {gwas$ID[i] <- paste(chrom_names[gwas$CHROM[i]], gwas$POS[i], sep='_')}
  
  dbSendQuery(gwas_snp_con, paste("insert into gwas_snp
                                   values
                                   ('", gwas$ID[i], "', '", chrom_names[gwas$CHROM[i]], "', ", gwas$POS[i], ", ", gwas$P[i], ")", sep='')) 
  }


#make dataframe of significant snps located in genes
gs_gene <- dbGetQuery(gwas_snp_con, paste("select snp_id, gs.chrom, gs.pos, P, ncbi_gene_id, group_concat(distinct synonyms order by synonyms) synonyms, group_concat(distinct desc order by desc) desc
                                           from gwas_snp gs, gene_info gi, gene_xrefs gx, og2genes og, og2genes og2, genes g
                                           where pos between start and end
                                           and gs.chrom = gi.chrom
                                           and gi.ncbi_gene_id = external_id
                                           and gx.orthodb_gene_id = og.orthodb_gene_id
                                           and og.orthogroup_id = og2.orthogroup_id
                                           and og.orthodb_gene_id != og2.orthodb_gene_id
                                           and og2.orthodb_gene_id = g.orthodb_gene_id
                                           and orthodb_species_id = '", mel_species_id, "'
                                           group by snp_id, gs.chrom, gs.pos, P, ncbi_gene_id
                                           order by gs.chrom, gs.pos, ncbi_gene_id", sep='')) 
                                 

#make dataframe of significant snps and go terms
gs_go <- dbGetQuery(gwas_snp_con, paste("select snp_id, gs.chrom, gs.pos, P, ncbi_gene_id, group_concat(distinct synonyms order by synonyms) synonyms, group_concat(distinct desc order by desc) desc, go_id, name
                                         from gwas_snp gs, gene_info gi, gene_xrefs gx, og2genes og, og2genes og2, genes g, gene_xrefs gx2, go
                                         where pos between start and end
                                         and gs.chrom = gi.chrom
                                         and gi.ncbi_gene_id = gx.external_id
                                         and gx.orthodb_gene_id = og.orthodb_gene_id
                                         and og.orthogroup_id = og2.orthogroup_id
                                         and og.orthodb_gene_id != og2.orthodb_gene_id
                                         and og2.orthodb_gene_id = g.orthodb_gene_id
                                         and orthodb_species_id = '", mel_species_id, "'
                                         and g.orthodb_gene_id = gx2.orthodb_gene_id
                                         and gx2.external_source = 'GOterm' 
                                         and gx2.external_id = go.go_id
                                         and go.namespace = 'biological_process'
                                         group by snp_id, gs.chrom, gs.pos, P, ncbi_gene_id, go_id, name
                                         order by gs.chrom, gs.pos, ncbi_gene_id", sep='')) 


#make dataframe of significant snps located in intergenic regions and identify nearest upstream and downstream genes
gs_intergenic <- dbGetQuery(gwas_snp_con, "select snp_id, gs.chrom, gs.pos, P
                                           from gwas_snp gs
                                           where not exists (select 'x'
                                                             from gene_info gi
                                                             where pos between start and end
                                                             and gs.chrom = gi.chrom)               
                                           order by gs.chrom, pos") 

nearest_gemes <- data.frame()                                           
for (i in 1:nrow(gs_intergenic)) 
  {
  upstream_gene <- dbGetQuery(gwas_snp_con, paste("select gi.ncbi_gene_id upstream_ncbi_gene_id, ", gs_intergenic$pos[i], " - gi.end upstream_dist, group_concat(distinct synonyms order by synonyms) upstream_synonyms
                                                   from gene_info gi, gene_xrefs gx, og2genes og, og2genes og2, genes g
                                                   where gi.chrom = '", gs_intergenic$chrom[i], "'
                                                   and gi.end < ", gs_intergenic$pos[i], "
                                                   and gi.ncbi_gene_id = external_id
                                                   and gx.orthodb_gene_id = og.orthodb_gene_id
                                                   and og.orthogroup_id = og2.orthogroup_id
                                                   and og.orthodb_gene_id != og2.orthodb_gene_id
                                                   and og2.orthodb_gene_id = g.orthodb_gene_id
                                                   and orthodb_species_id = '", mel_species_id, "'
                                                   group by gi.ncbi_gene_id, ", gs_intergenic$pos[i], " - gi.end
                                                   order by ", gs_intergenic$pos[i], " - gi.end
                                                   limit 1", sep='')) 

  downstream_gene <- dbGetQuery(gwas_snp_con, paste("select gi.ncbi_gene_id downstream_ncbi_gene_id, gi.start - ", gs_intergenic$pos[i], " downstream_dist, group_concat(distinct synonyms order by synonyms) downstream_synonyms
                                                     from gene_info gi, gene_xrefs gx, og2genes og, og2genes og2, genes g
                                                     where gi.chrom = '", gs_intergenic$chrom[i], "'
                                                     and gi.start > ", gs_intergenic$pos[i], "
                                                     and gi.ncbi_gene_id = external_id
                                                     and gx.orthodb_gene_id = og.orthodb_gene_id
                                                     and og.orthogroup_id = og2.orthogroup_id
                                                     and og.orthodb_gene_id != og2.orthodb_gene_id
                                                     and og2.orthodb_gene_id = g.orthodb_gene_id
                                                     and orthodb_species_id = '", mel_species_id, "'
                                                     group by gi.ncbi_gene_id, gi.start - ", gs_intergenic$pos[i], "
                                                     order by gi.start - ", gs_intergenic$pos[i], "
                                                     limit 1", sep='')) 
                                                   
  nearest_gemes <- rbind(nearest_gemes, cbind(upstream_gene, downstream_gene))                                                    
  }                                  

gs_intergenic <- cbind(gs_intergenic, nearest_gemes)


#make dataframe of seq types (CDS, intron, 3prime_utr, 5prime utr) for significant genic snps
gs_seq_type <- dbGetQuery(gwas_snp_con, "select snp_id, 'cds' seq_type, gs.chrom, gs.pos, P, ncbi_gene_id, group_concat(distinct ti.transcript_id order by ti.transcript_id) transcripts
                                         from gwas_snp gs, cds_info ci, transcript_info ti
                                         where pos between start and end
                                         and gs.chrom = ci.chrom
                                         and ci.transcript_id = ti.transcript_id
                                         group by snp_id, gs.chrom, pos, P, ncbi_gene_id
                                         union all
                                         select snp_id, 'intron', gs.chrom, gs.pos, P, ncbi_gene_id, group_concat(distinct ti.transcript_id order by ti.transcript_id) transcripts
                                         from gwas_snp gs, intron_info ii, transcript_info ti
                                         where pos between start and end
                                         and gs.chrom = ii.chrom
                                         and ii.transcript_id = ti.transcript_id
                                         group by snp_id, gs.chrom, gs.pos, P, ncbi_gene_id
                                         union all
                                         select snp_id, concat(utr_type, ' utr'), gs.chrom, gs.pos, P, ncbi_gene_id, group_concat(distinct ti.transcript_id order by ti.transcript_id) transcripts
                                         from gwas_snp gs, utr_info ui, transcript_info ti
                                         where pos between start and end
                                         and gs.chrom = ui.chrom
                                         and ui.transcript_id = ti.transcript_id
                                         group by snp_id, gs.chrom, gs.pos, P, ncbi_gene_id, utr_type
                                         order by gs.chrom, pos, seq_type") 

#write dataframes to text files in gwas results dir
write.table(gs_gene, gene_file, quote=F, sep='\t', row.names=F)
write.table(gs_go, go_file, quote=F, sep='\t', row.names=F)
write.table(gs_intergenic, intergenic_file, quote=F, sep='\t', row.names=F)
write.table(gs_seq_type, seq_type_file, quote=F, sep='\t', row.names=F)
