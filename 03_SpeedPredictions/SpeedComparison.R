# my speed estimates
# author: Alexander Stein

library("tidyverse")

#setwd("~/ETH/Master_Thesis/Code/speedcomparison")
source("speedestimates.R")

# Define parameters for the "expected dispersal speed for migration model" function
r1 <- 1
r2 <- 1.1
s <- r2-r1
K <- 256
m <- 0.01
migration_type <- 0
symmetric <- FALSE
two_dim <- FALSE
d <- NA
approx <- FALSE

disp_rate(r1, r2, K, m, migration_type, symmetric, two_dim, d, approx)


### Plot over migration rate (m in [0,1])
ms <- 1:50
our_speed <- 1:50
fisher_speed <- 1:50
mf_speed <- 1:50
noise_speed <- 1:50

# Speed values found in the demon simulation for s=0.01, K=256
mig <- c(0.001,0.005,0.01,0.05,0.1,0.2,0.5)
speed <- c(3.21,8.30,12.7,33.5,54.0,90.1,189.)/sqrt(256)
std <- c(0.66,0.54,0.69,1.39,2.52,2.32,6.62)/sqrt(256)
demon_m <- data.frame(mig, speed, std)

for (k in 1:1000) {
  ms[k] <- k/1000
  our_speed[k] <- disp_rate(r1, r2, K, ms[k], migration_type, symmetric, two_dim, d, approx)
  fisher_speed[k] <- 2*sqrt(ms[k]*K*s) # using c = 2*sqrt(Ds) with D=ml^2/2=mK/2, but m=2m
  mf_speed[k] <- 2*sqrt(ms[k]*K*s*(s+1)) # mean field from Houchmandzadeh and Vallade, but m=2m
  #noise_speed[k] <- 2*K*ms[k]*(r2-r1) # Halaltschek and Koroloev
}

speed_over_migration <- data.frame(ms, our_speed, fisher_speed, mf_speed, noise_speed)

ggplot()+geom_line(data=speed_over_migration, aes(x=ms,y=our_speed), color="blue", linetype = "solid", size=1)+
  geom_line(data=speed_over_migration, aes(x=ms,y=fisher_speed), color="orange", linetype = "dashed", size=1)+
  geom_line(data=speed_over_migration, aes(x=ms,y=mf_speed), color="red", linetype = "dotted", size=1)+
  #geom_line(data=speed_over_migration, aes(x=ms,y=noise_speed), color="yellow")+
  geom_point(data=demon_m, aes(x=mig, y=speed), color="black", size=2.5) + #scale_x_log10() + scale_y_log10() +
  geom_errorbar(data=demon_m, aes(x=mig,ymin=speed-std, ymax=speed+std), width=0.0, size=1) +
  xlab("Migration Rate m") + #+ theme(axis.title.x = element_text(face="bold",size=15)) +
  ylab("Speed of Propagation c") + #+ theme(axis.title.y = element_text(face="bold",size=15))
  theme_bw(base_size=22,base_family = "bold", base_line_size = 1.0, base_rect_size = 1.0)
ggsave(filename="Speed_over_migration.png")



### Plot over selection stength s=r2-r1 (r2 in [1,2])
r2s <- 1:100
ss <- 1:100
ms <- 1:100
our_speed <- 1:100
fisher_speed <- 1:100
mf_speed <- 1:100
noise_speed <- 1:100

# Speed values found in the demon simulation for m=0.01, K=256
sel <- c(0.001, 0.005, 0.01, 0.05, 0.1, 0.2, 0.5, 1.0)
speed <- c(2.30,3.24,4.01,8.29,12.7,18.0,29.0,39.6)/sqrt(256)
std <- c(0.464, 0.660, 0.829, 0.639, 0.691, 0.466, 0.601, 1.11)/sqrt(256)
demon_s <- data.frame(sel,speed,std)


for (k in 1:1000) {
  r2s[k] <- 1+(k-1)/1000
  ss[k] <- r2s[k]-r1
  our_speed[k] <- disp_rate(r1, r2s[k], K, m, migration_type, symmetric, two_dim, d, approx)
  fisher_speed[k] <- 2*sqrt(m*K*ss[k]) # using c = 2*sqrt(Ds) with D=ml^2/2=mK/2
  mf_speed[k] <- 2*sqrt(m*K*ss[k]*(ss[k]+1)) # mean field from Houchmandzadeh and Vallade but m=2m
  #noise_speed[k] <- 2*K*m*(r2-r1) # Halaltschek and Koroloev but m=2m
}

speed_over_selection <- data.frame(s, our_speed, fisher_speed, mf_speed, noise_speed)

ggplot()+geom_line(data=speed_over_selection, aes(x=ss,y=our_speed), color="blue", linetype="solid", size=1.0)+
  geom_line(data=speed_over_selection, aes(x=ss,y=fisher_speed), color="orange",linetype = "dashed", size=1.0)+
  geom_line(data=speed_over_selection, aes(x=ss,y=mf_speed), color="red",linetype = "dotted", size=1.0)+
  #geom_line(data=speed_over_selection, aes(x=ss,y=noise_speed), color="yellow")+
  geom_point(data=demon_s, aes(x=sel, y=speed), color="black", size=2.5) +
  geom_errorbar(data=demon_s, aes(x=sel,ymin=speed-std, ymax=speed+std), width=0.0, size=1.0) +
  xlab("Selection strength s") + ylab("Speed of Propagation c") + scale_x_log10() + scale_y_log10() +
  theme_bw(base_size=22,base_family = "bold", base_line_size = 1.0, base_rect_size = 1.0)
ggsave(filename="Speed_over_selection_loglog.png")


### Plot over deme size K = 10,20,...,500
Ks <- 1:500
our_speed <- 1:500
fisher_speed <- 1:500
sfisher_speed <- 1:500
mf_speed <- 1:500
noise_speed <- 1:500


# Data points from simulation
K_2 <- c(2,4,8,16,32,64,128,256,512,1024)
speed <- c(0.01291,0.03387,0.09319,0.24730,0.78790,2.0400,5.395,12.73,27.5,59.41)/sqrt(K_2)
std <- c(0.001727,0.010712,0.022277,0.0333,0.0954,0.2136,0.2896,0.6912,1.3581,1.7272)/sqrt(K_2)
  
  
demon_K <- data.frame(K_2, speed, std)


for (k in 1:500) {
  Ks[k] <- 2*k
  our_speed[k] <- disp_rate(r1, r2, Ks[k], m, migration_type, symmetric, two_dim, d, approx)
  fisher_speed[k] <- 2*sqrt(m*Ks[k]*s) # using c = 2*sqrt(Ds) with D=ml^2/2=mK/2
  sfisher_speed[k] <- ( 2-(pi/(log(Ks[k])))^2 ) *sqrt(m*Ks[k]*s)
  mf_speed[k] <- 2*sqrt(m*Ks[k]*s*(s+1)) # mean field from Houchmandzadeh and Vallade
  noise_speed[k] <- 2*Ks[k]*m*(r2-r1) # Halaltschek and Koroloev with c=2
}


speed_over_size <- data.frame(Ks, our_speed, fisher_speed, sfisher_speed ,mf_speed, noise_speed)

ggplot()+geom_line(data=speed_over_size, aes(x=Ks,y=our_speed), color="blue", linetype = "solid", size=1)+
  geom_line(data=speed_over_size, aes(x=Ks,y=fisher_speed), color="orange", linetype = "dashed", size=1)+
  geom_line(data=speed_over_size, aes(x=Ks,y=sfisher_speed), color="purple", linetype = "dotdash", size=1)+
  geom_line(data=speed_over_size, aes(x=Ks,y=mf_speed), color="red", linetype = "dotted", size=1)+
  #geom_line(data=speed_over_size, aes(x=Ks,y=noise_speed), color="yellow")+
  geom_point(data=demon_K, aes(x=K_2, y=speed), color="black", size=2.5) +
  geom_errorbar(data=demon_K, aes(x=K_2,ymin=speed-std, ymax=speed+std), width=0.0, size=1) +
  xlab("Deme size K")+ylab("Speed of Propagation c") + #scale_x_log10() + scale_y_log10() +
  theme_bw(base_size=22,base_family = "bold", base_line_size = 1.0, base_rect_size = 1.0) +
  ylim(0,2)
ggsave(filename="Speed_over_size.png")




