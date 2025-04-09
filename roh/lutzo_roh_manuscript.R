library("RSQLite")  


#define file paths
base_dir <- '/work/users/d/t/dturissi/lutzomyia/roh'

db_file <-  paste('lutzo_roh.db', sep='')

setwd(base_dir)

#create db connection
conn <- dbConnect(dbDriver("SQLite"), db_file)


dbSendQuery(conn, "attach '/work/users/d/t/dturissi/lutzomyia/lutzomyia_metadata.db' as m")


pop_cols <- dbGetQuery(conn, "select distinct pop_manuscript, pop_col
                              from lutzo_roh r, lutzomyia_pops p
                              where min_kb = 50
                              and r.pop = p.pop
                              order by pop_manuscript")



lutzo_roh <- dbGetQuery(conn, "select r.*, pop_manuscript, pop_col
                               from lutzo_roh r, lutzomyia_pops p
                               where min_kb = 50
                               and r.pop = p.pop")

lutzo_roh_regions <- dbGetQuery(conn, "select r.*, (end - start + 1) roh_len,(end + start) / 2 mid, pop_manuscript, pop_col
                                       from lutzo_roh_regions r, lutzomyia_pops p
                                       where min_kb = 50
                                       and r.pop = p.pop")
                                             



pdf('lutzo_roh_manuscript.pdf', height=8, width=12)
par(mfrow=c(2,3))
for (pop in pop_cols$pop_manuscript)
  {
  pop_filter <- lutzo_roh$pop_manuscript == pop
  hist(lutzo_roh$sum_roh_kb[pop_filter], breaks=seq(0, max(lutzo_roh$sum_roh_kb) + 500, 500), col=lutzo_roh$pop_col[pop_filter], border=NA, xlab='Total ROH kb', ylab='Individuals', main=pop)
  }
par(mfrow=c(1,1))


plot(1, type='n', xlim=c(0, max(lutzo_roh$num_roh_regions)), ylim=c(0, max(lutzo_roh$sum_roh_kb)), xlab='ROH regions', ylab='Total ROH kb', main=c('Lutzomyia ROH', 'min kb = 50'))
legend('topleft', pop_cols$pop_manuscript, fill = pop_cols$pop_col, border=NA)
for (pop in pop_cols$pop_manuscript)
  {
  pop_filter <- lutzo_roh$pop_manuscript == pop
  points(lutzo_roh$num_roh_regions[pop_filter], lutzo_roh$sum_roh_kb[pop_filter], pch=20, col=lutzo_roh$pop_col[pop_filter])
  } 


par(mfrow=c(2,3))
for (pop in pop_cols$pop_manuscript)
  {
  pop_filter <- lutzo_roh_regions$pop_manuscript == pop
  hist(lutzo_roh_regions$roh_len[pop_filter], breaks=seq(0, max(lutzo_roh_regions$roh_len) + 10000, 10000), col=lutzo_roh_regions$pop_col[pop_filter], border=NA, xlab='ROH length (kb)', ylab='Num ROH regions', main=pop)
  }
par(mfrow=c(1,1))   

par(mfrow=c(2,3))
for (chrom in sort(unique(lutzo_roh_regions$chrom)))
  {
  for (pop in pop_cols$pop_manuscript)
    {
    pop_chrom_filter <- lutzo_roh_regions$pop_manuscript == pop & lutzo_roh_regions$chrom == chrom
    hist(lutzo_roh_regions$mid[pop_chrom_filter], breaks=seq(0, max(lutzo_roh_regions$start[lutzo_roh_regions$chrom == chrom]) + 1e6, 1e6), col=lutzo_roh_regions$pop_col[pop_chrom_filter], border=NA, xlab=paste('Chromosome', chrom), ylab='Num ROH regions', main=pop)
    }
  }
par(mfrow=c(1,1))  
dev.off()  
