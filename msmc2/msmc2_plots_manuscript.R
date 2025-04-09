library(RSQLite)

setwd('/work/users/d/t/dturissi/lutzomyia/msmc2')
#generation time of 7 weeks
gen <- (7 * 7) / 365
mu <- 1.25e-8
time_pattern <- '1-6+8-1+1-4'

db_file <- '/work/users/d/t/dturissi/lutzomyia/msmc2/msmc2_results.db'


conn <- dbConnect(dbDriver("SQLite"), db_file)
dbSendQuery(conn, "attach '/work/users/d/t/dturissi/lutzomyia/lutzomyia_metadata.db' as m")




pop_cols <- dbGetQuery(conn, paste("select distinct pop_manuscript, pop_col
                                    from msmc2_results m, lutzomyia_samples s, lutzomyia_pops p
                                    where m.sample = s.sample
                                    and s.pop = p.pop
                                    and pattern = '", time_pattern, "'
                                    and chrom = 'genome'
                                    order by pop_manuscript", sep=''))

country_cols <- dbGetQuery(conn, paste("select distinct country, country_col
                                        from msmc2_results m, lutzomyia_samples s, lutzomyia_pops p
                                        where m.sample = s.sample
                                        and s.pop = p.pop
                                        and pattern = '", time_pattern, "'
                                        and chrom = 'genome'
                                        order by country", sep=''))

msmc2 <- dbGetQuery(conn, paste("select m.sample, chrom, time_index, left_time, right_time, lambda, country, pop_manuscript, sex, pop_col, country_col
                                  from msmc2_results m, lutzomyia_samples s, lutzomyia_pops p
                                  where m.sample = s.sample
                                  and s.pop = p.pop
                                  and pattern = '", time_pattern, "'
                                  and chrom = 'genome'
                                  and left_time > 0
                                  order by m.sample, chrom, time_index", sep=''))
                                  

                      
samples <- unique(msmc2$sample)
br_samples <- unique(msmc2$sample[msmc2$country == 'Brazil'])

pdf('lutzomyia_msmc2_manuscript.pdf', height=8, width=10.5)
#plot countries
plot(msmc2$left_time / mu * gen, (1 / msmc2$lambda) / (2 * mu), type='n', log='x', xlab="Years ago", ylab=expression("Effective population size (N"['e']*")"), main='')
legend('topright', country_cols[,1], fill=country_cols[,2], border=NA)

for (sample in samples)
  {
  df_filter <- msmc2$sample == sample
  lines(msmc2$left_time[df_filter] / mu * gen, (1 / msmc2$lambda[df_filter]) / (2 * mu), type="s", col=msmc2$country_col[df_filter])  
  }


#plot all pops
plot(msmc2$left_time / mu * gen, (1 / msmc2$lambda) / (2 * mu), type='n', log='x', xlab="Years ago", ylab=expression("Effective population size (N"['e']*")"), main='')
legend('topright', pop_cols$pop_manuscript, fill=pop_cols$pop_col, border=NA)
for (sample in samples)
  {
  df_filter <- msmc2$sample == sample
  lines(msmc2$left_time[df_filter] / mu * gen, (1 / msmc2$lambda[df_filter]) / (2 * mu), type="s", col=msmc2$pop_col[df_filter])  
  }


  #plot only Brazil pops
plot(msmc2$left_time / mu * gen, (1 / msmc2$lambda) / (2 * mu), type='n', log='x', xlab="Years ago", ylab=expression("Effective population size (N"['e']*")"), main='')
legend('topright', pop_cols$pop_manuscript[pop_cols$pop_manuscript != 'Colombia'], fill=pop_cols$pop_col[pop_cols$pop_manuscript != 'Colombia'], border=NA)
for (sample in br_samples)
  {
  df_filter <- msmc2$sample == sample
  lines(msmc2$left_time[df_filter] / mu * gen, (1 / msmc2$lambda[df_filter]) / (2 * mu), type="s", col=msmc2$pop_col[df_filter])  
  }
dev.off()

