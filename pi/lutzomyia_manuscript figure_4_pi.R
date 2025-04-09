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
                                and s.pop in ('Colombia', 'Jacobina', 'Laphina', 'Marajo', 'Sobral_1S', 'Sobral_2S')
                                group by p.pop, pop_manuscript", sep=''))

num_pops <- nrow(pops)

pi <- dbGetQuery(conn, paste("select *, (start + end) / 2.0 pos, 
                                     case when chrom in ('chr_1', 'chr_3') then 'black' else 'blue' end chrom_col
                              from ", poly_table, "
                              where pop in ('All', 'Colombia', 'Jacobina', 'Laphina', 'Marajo', 'Sobral_1S', 'Sobral_2S')", sep=''))

pi_adj_pos <- pi$pos
for (i in 1:(nrow(chroms) - 1))
  {
  pi_adj_pos[pi$chrom > chroms$chrom[i]] <- pi_adj_pos[pi$chrom > chroms$chrom[i]] + chroms$chrom_len[i]
  }






pdf_file <- paste('lutzo_manuscript_figure_4_pi.pdf', sep='')
pdf(pdf_file, height=24, width=12)
layout(matrix(1:8, ncol=2, byrow=T))
par(mar=c(5.1, 4.6, 4.1, 2.1))
###########
#A pi all
###########
hist(pi$pi[pi$pop == 'All'], breaks=seq(0, max(pi$pi) + 0.001, 0.001), cex.main=2, cex.axis=1.5, cex.lab=1.75, col='black', xlab=bquote(pi), ylab='10kb windows', main=bquote("All" ~ pi))
mtext('A)', side=3, line=1.5, at= -.06 * max(pi$pi) + 0.001, cex=2.2)



###########
#B pi Colombia
###########
hist(pi$pi[pi$pop == 'Colombia'], breaks=seq(0, max(pi$pi) + 0.001, 0.001), cex.main=2, cex.axis=1.5, cex.lab=1.75, col='black', xlab=bquote(pi), ylab='10kb windows', main=bquote("Colombia" ~ pi))
mtext('B)', side=3, line=1.5, at= -.06 * max(pi$pi) + 0.001, cex=2.2)


###########
#C pi Jacobina
###########
hist(pi$pi[pi$pop == 'Jacobina'], breaks=seq(0, max(pi$pi) + 0.001, 0.001), cex.main=2, cex.axis=1.5, cex.lab=1.75, col='black', xlab=bquote(pi), ylab='10kb windows', main=bquote("Jacobina" ~ pi))
mtext('C)', side=3, line=2, at= -.06 * max(pi$pi) + 0.001, cex=2.2)


###########
#D pi Marajo
###########
hist(pi$pi[pi$pop == 'Marajo'], breaks=seq(0, max(pi$pi) + 0.001, 0.001), cex.main=2, cex.axis=1.5, cex.lab=1.75, col='black', xlab=bquote(pi), ylab='10kb windows', main=bquote("Marajo" ~ pi))
mtext('D)', side=3, line=2, at= -.06 * max(pi$pi) + 0.001, cex=2.2)


###########
#E pi Sobral 2S
###########
hist(pi$pi[pi$pop == 'Sobral_2S'], breaks=seq(0, max(pi$pi) + 0.001, 0.001), cex.main=2, cex.axis=1.5, cex.lab=1.75, col='black', xlab=bquote(pi), ylab='10kb windows', main=bquote("Sobral 2S" ~ pi))
mtext('E)', side=3, line=2, at= -.06 * max(pi$pi) + 0.001, cex=2.2)


###########
#F pi along the genome in Jacobina. 
###########
pop <- 'Jacobina'
plot(1, type='n', xlim=c(0, max(pi_adj_pos)), cex.main=2, cex.axis=1.5, cex.lab=1.75, ylim=c(0, max(pi$pi)), xaxt='n', xlab='Chromosome', ylab=bquote(pi), main=bquote(.(pop) ~ pi))
mtext('F)', side=3, line=2, at= -.06 * max(pi_adj_pos), cex=2.2)

for (chr_i in 1:4)
  {
  chr <- paste('chr_', chr_i, sep='')
  points(pi_adj_pos[pi$pop == pop & pi$chrom == chr], pi$pi[pi$pop == pop & pi$chrom == chr], type='l', col=pi$chrom_col[pi$pop == pop & pi$chrom == chr])
  mtext(chr_i, 1, at=chrom_label_pos[chr_i], line=1, cex=1.5)
  }




###########
#G pi along the genome in Marajo
###########
pop <- 'Marajo'
plot(1, type='n', xlim=c(0, max(pi_adj_pos)), cex.main=2, cex.axis=1.5, cex.lab=1.75, ylim=c(0, max(pi$pi)), xaxt='n', xlab='Chromosome', ylab=bquote(pi), main=bquote(.(pop) ~ pi))
mtext('G)', side=3, line=2, at= -.06 * max(pi_adj_pos), cex=2.2)

for (chr_i in 1:4)
  {
  chr <- paste('chr_', chr_i, sep='')
  points(pi_adj_pos[pi$pop == pop & pi$chrom == chr], pi$pi[pi$pop == pop & pi$chrom == chr], type='l', col=pi$chrom_col[pi$pop == pop & pi$chrom == chr])
  mtext(chr_i, 1, at=chrom_label_pos[chr_i], line=1, cex=1.5)
  }



###########
#H pi pairwise comparisons
###########
par(mar=c(2.1, 7.1, 7.1, 2.1))
plot(1, type='n', bty='n', xaxt='n', yaxt='n', xlim=c(0.5, num_pops - 0.5), ylim=c(1.5, num_pops + 0.5), xlab='', ylab='', main='')
axis(2, at=2:num_pops, labels=pops$pop_manuscript[2:num_pops], las=2, cex.axis=1.5)
axis(2, at=-0.2 + 2:num_pops, tick=F, line=F, labels=as.expression(sapply(pops$avg_pi[2:num_pops], function(var) bquote(pi == .(var)))), las=2, cex.axis=1.5)
axis(3, at=1:(num_pops - 1), labels=pops$pop_manuscript[1:(num_pops - 1)], las=2, cex.axis=1.5)
axis(3, at=0.2 + 1:(num_pops - 1), tick=F, line=F, labels=as.expression(sapply(pops$avg_pi[1:(num_pops - 1)], function(var) bquote(pi == .(var)))), las=2, cex.axis=1.5)


mtext('H)', side=3, line=2, at= 0.5 - .06 * (num_pops + 0.5), cex=2.2)
legend(3.5, 3, c('No effect', 'Small effect', 'Moderate effect', 'Large effect'), fill=adjustcolor(c('gray', 'red', 'green', 'blue'), .2), border=NA, cex=1.5)

#mann-whitney u test
for (i in 1:num_pops)
  {
  for (j in (i + 1):num_pops)
    {
    if (j <= num_pops & i < j)
      {
      cat(i, j, pops$pop[i], pops$pop[j], "\n")
      poly_wins <- dbGetQuery(conn, paste("select pi, 'pop_1' pop, chrom, start
                                    from ", poly_table, " p
                                    where pop = '", pops$pop[i], "'
                                    and exists (select 'x'
                                                from ", poly_table, " p2
                                                where p2.pop = '", pops$pop[j], "'
                                                and p.chrom = p2.chrom
                                                and p.start = p2.start)
                                    union all
                                    select pi, 'pop_2' pop, chrom, start
                                    from ", poly_table, " p
                                    where pop = '", pops$pop[j], "'
                                    and exists (select 'x'
                                                from ", poly_table, " p2
                                                where p2.pop = '", pops$pop[i], "'
                                                and p.chrom = p2.chrom
                                                and p.start = p2.start)
                                    order by pop, chrom, start", sep=''))

      poly_wins_xy <- dbGetQuery(conn, paste("select p.pi pi_1, p2.pi pi_2
                                           from ", poly_table, " p, ", poly_table, " p2
                                           where p.chrom = p2.chrom
                                           and p.start = p2.start
                                           and p.pop = '", pops$pop[i], "'
                                           and p2.pop = '", pops$pop[j], "'", sep=''))

      effect_size <- wilcox_effsize(poly_wins, pi ~ pop, paired=T)$effsize
      p <- wilcox.test(Pair(pi_1, pi_2) ~ 1, data=poly_wins_xy)$p.value
      
      if (effect_size < 0.1)
        {
        effect_col <- adjustcolor('grey', 0.2)
        } else if (effect_size < 0.3) {
        effect_col <- adjustcolor('red', 0.2)         
        } else if (effect_size < 0.5) {
        effect_col <- adjustcolor('green', 0.2)         
        } else {
        effect_col <- adjustcolor('blue', 0.2)
        }
            
      rect(i - 0.5, j - 0.5, i + 0.5, j + 0.5, col=effect_col, border='black')
      text(i, j + .1, paste('Effect size =', round(effect_size, 2)), cex=.95)
      text(i, j - .1, paste('p =', signif(p, 2)), cex=.95)
      }
    }
  }
par(mar=c(5.1, 4.1, 4.1, 2.1))
layout(1)
dev.off()


