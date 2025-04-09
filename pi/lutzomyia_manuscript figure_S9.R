library("RSQLite")  
library(rstatix)


win_size <- 10000

#define file paths
base_dir <- '/work/users/d/t/dturissi/lutzomyia/poly_windows'
db_file <- paste(base_dir, '/lutzomyia_poly_win.db', sep='')

setwd(base_dir)

#db connection
conn <- dbConnect(dbDriver("SQLite"), db_file)


dbSendQuery(conn, "attach database '/work/users/d/t/dturissi/lutzomyia/sweep/sweepfinder/sobral_1s_subpop_sweepfinder/lutzo_sweepfinder_sobral_1s_subpop.db' as s")
dbSendQuery(conn, "attach database '/work/users/d/t/dturissi/lutzomyia/lutzomyia_metadata.db' as m")


chroms <- dbGetQuery(conn, "select chrom, chrom_len 
                            from chrom_lens")

chrom_label_pos <- c(chroms$chrom_len[chroms$chrom == chroms$chrom[1]] / 2)
for (i in 2:(nrow(chroms)))
  {
  chrom_label_pos <- c(chrom_label_pos, sum(chroms$chrom_len[1:(i - 1)]) + chroms$chrom_len[chroms$chrom == chroms$chrom[i]] / 2)
  }



#load window data from db
poly_table <- paste('poly_win_pops_sobral_1s_subpop_', win_size, '_vcftools', sep='')
pops <- dbGetQuery(conn, paste("select p.pop, pop_manuscript, round(avg(pi), 4) avg_pi
                                from ", poly_table, " s, lutzomyia_pops p
                                where s.pop = p.pop
                                and s.pop in ('Colombia', 'Laphina', 'Sobral_1S', 'Sobral_2S')
                                group by p.pop, pop_manuscript", sep=''))

num_pops <- nrow(pops)

pi <- dbGetQuery(conn, paste("select *, (start + end) / 2.0 pos, 
                                     case when chrom in ('chr_1', 'chr_3') then 'black' else 'blue' end chrom_col
                              from ", poly_table, "
                              where pop in ('Colombia', 'Laphina', 'Sobral_1S', 'Sobral_2S')", sep=''))

pi_adj_pos <- pi$pos
for (i in 1:(nrow(chroms) - 1))
  {
  pi_adj_pos[pi$chrom > chroms$chrom[i]] <- pi_adj_pos[pi$chrom > chroms$chrom[i]] + chroms$chrom_len[i]
  }






pdf_file <- paste('lutzo_manuscript_figure_s9_pi.pdf', sep='')
pdf(pdf_file, height=12, width=12)
layout(matrix(1:4, ncol=2, byrow=T))
for (i in 1:nrow(pops))
  {
  pop <- pops$pop[i]
  pop_manuscript <- pops$pop_manuscript[i]
  plot(1, type='n', xlim=c(0, max(pi_adj_pos)), cex.main=2, cex.axis=1.5, cex.lab=1.75, ylim=c(0, max(pi$pi)), xaxt='n', xlab='Chromosome', ylab=bquote(pi), main=bquote(.(pop_manuscript) ~ pi))
  mtext(paste(LETTERS[i], ')', sep=''), side=3, line=2, at= -.06 * max(pi_adj_pos), cex=2.2)

  for (chr_i in 1:4)
    {
    chr <- paste('chr_', chr_i, sep='')
    points(pi_adj_pos[pi$pop == pop & pi$chrom == chr], pi$pi[pi$pop == pop & pi$chrom == chr], type='l', col=pi$chrom_col[pi$pop == pop & pi$chrom == chr])
    mtext(chr_i, 1, at=chrom_label_pos[chr_i], line=1, cex=1.5)
    }
  }
layout(1)
dev.off()


