# 2021 01 21 I.Zliobaite

set.seed(1981)

# parameters
n_years <- 1200 #also number of species
max_duration <- 100
max_peak <- 100

all_correlations <- c()
for (ff in seq(-3,2,0.1)){
  fossilization_ratio <- 10^ff
  print(fossilization_ratio)
  
  real_species <- matrix(0, n_years, n_years*2)
  fossilized_species <- matrix(0, n_years, n_years*2)
  
  spcies_statistics <- c()
  
  for (sk in 1:n_years){
    
    dur <- sample(3:max_duration, 1)
    pek <- sample(1:max_peak, 1)
    
    spcies_statistics <- rbind(spcies_statistics,c(dur,pek))
    
    if (dur %% 2 == 0){
      half_life <- dur/2
      step <- pek/half_life  
      species_all <- 1
      species_now <- 1
      for (sk2 in 2:(half_life)){
        species_now <- species_now + step
        species_all <- c(species_all,species_now)
      }  
      species_all <- c(species_all,species_now)
      for (sk2 in (half_life+2):dur){
        species_now <- species_now - step
        species_all <- c(species_all,species_now)
      }
    }else{
      half_life <- ceiling(dur/2)
      step <- pek/half_life  
      species_all <- 1
      species_now <- 1
      for (sk2 in 2:(half_life)){
        species_now <- species_now + step
        species_all <- c(species_all,species_now)
      }  
      for (sk2 in (half_life+1):dur){
        species_now <- species_now - step
        species_all <- c(species_all,species_now)
      }
    }
    
    start <- sk
    end <- sk+dur-1
    real_species[sk,start:end] <- species_all
  }
  
  real_diversity <- c()
  for (sk in 1:(n_years*2)){
    div <- sum(real_species[,sk]>0)
    real_diversity <- rbind(real_diversity,div)
  }

  sample_diversity <- c()
  for (sk in 1:(n_years+50)){
    n_fossils <- floor(real_abundance[sk]*fossilization_ratio)
    pp <- real_species[,sk]
    ss <- sample(1:n_years,n_fossils,replace=TRUE,prob=pp)
    un_sample <- unique(ss)
    for (sk2 in 1:length(un_sample)){
      u_now <- un_sample[sk2]
      fossilized_species[u_now,sk] <- sum(ss==u_now)
    }
    sample_diversity <- rbind(sample_diversity,length(un_sample))
  }

  stats <- c()
  
  for (sk in 1:n_years){
    ind <- which(fossilized_species[sk,]>0)
    mn <- min(ind)
    mx <- max(ind)
    dur <- mx-mn+1
    pk <- max(fossilized_species[sk,])
    ind <- which(fossilized_species[sk,]==pk)
    ind <- ind[1]
    stats <- rbind(stats,c(dur,pk,mn,mx,ind))
  }
  

  spcies_statistics <- cbind(spcies_statistics,stats)
  spcies_statistics <- spcies_statistics[101:1100,]
  
  ind_sp <- which(spcies_statistics[,3]>0)
  spcies_statistics <- spcies_statistics[ind_sp,]
  real_diversity <- real_diversity[ind_sp]
  
  cor_dur_pek <- cor(spcies_statistics[,3],spcies_statistics[,4])
  
  all_correlations <- rbind(all_correlations, cbind(ff,fossilization_ratio,cor_dur_pek,length(ind_sp), mean(real_diversity),mean(sample_diversity)))
  
}

pdf('fig_correlations.pdf',height = 6, width = 6)
plot(all_correlations[,2],all_correlations[,3], pch = 16, log="x", ylim = c(0,0.5), xlab = "Fossilization ratio", ylab = "Correlation between sp. duration and abundance", xaxt = "n")
axis(1, at=c(0.001,0.01,0.1,1,10,100),labels = c("0.001",'0.01',"0.1",'1',"10",'100'))
dev.off()

pdf('fig_loop_diversity.pdf',height = 6, width = 6)
plot(all_correlations[,2],all_correlations[,5], type = 'l', log="x", ylim = c(0,55), xlab = "Fossilization ratio", ylab = "Species count", xaxt = "n",lwd = 2)
lines(all_correlations[,2],all_correlations[,6], col = 'brown', lwd = 2)
axis(1, at=c(0.001,0.01,0.1,1,10,100),labels = c("0.001",'0.01',"0.1",'1',"10",'100'))
legend(0.3, 7, legend=c("Simulated species", "Species in the 'fossil' record"),col=c("black", "brown"), lty=1, lwd = 2, cex = 0.8)
dev.off()