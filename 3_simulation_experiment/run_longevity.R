# 2021 01 21 I.Zliobaite

set.seed(1981)

# parameters
n_years <- 1200 #also number of species
max_duration <- 100
max_peak <- 100
fossilization_ratio <- 0.05

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

real_abundance <- apply(real_species,2,sum)

pdf('fig_real_abundance.pdf',height = 6, width = 10)
plot(real_abundance, type = 'l')
dev.off()

real_diversity <- c()
for (sk in 1:(n_years*2)){
  div <- sum(real_species[,sk]>0)
  real_diversity <- rbind(real_diversity,div)
}

pdf('fig_real_diversity.pdf',height = 6, width = 10)
plot(real_diversity, type = 'l')
dev.off()

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

sample_diversity2 <- c()
for (sk in 1:(n_years+50)){
  pp <- real_species[,sk]
  ss <- sample(1:n_years,50,replace=TRUE,prob=pp)
  un_sample <- unique(ss)
  sample_diversity2 <- rbind(sample_diversity2,length(un_sample))
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

pdf('fig_hats.pdf',height = 4, width = 10)
plot(NA,NA,xlim = c(101,300),ylim = c(0,100),xlab = 'Time step', ylab = 'Relative abundance', xaxt='n', yaxt='n', bty="n")
for (sk in seq(1,n_years,15)){
  ind <- which(real_species[sk,]>0)
  lines(ind,real_species[sk,ind],lwd = 2)
}
dev.off()

spcies_statistics <- cbind(spcies_statistics,stats)
spcies_statistics <- spcies_statistics[101:1100,]

pdf('fig_fossils.pdf',height = 6, width = 6)
plot(spcies_statistics[,3],spcies_statistics[,4])
dev.off()

pdf('fig_fossils_vs_real.pdf',height = 7, width = 7)
plot(NA,NA,xlim=c(0,100),ylim=c(0,100), xlab = 'True duration' , ylab = 'Duration in the "fossil" record')
abline(0,1)
points(spcies_statistics[,1],spcies_statistics[,3],pch = 16)
dev.off()

pdf('fig_diversity_both.pdf',height = 6, width = 10)
plot(real_diversity, type = 'l')
lines(sample_diversity, col='blue')
lines(sample_diversity2, col='red')
dev.off()

spcies_statistics <- spcies_statistics[spcies_statistics[,3]>0,]
print(cor(spcies_statistics[,3],spcies_statistics[,4]))

