#Functions updated 

library(popbio)
library(MASS)
library(boot)
library(ggplot2)
library(patchwork)

CV <- function(x){
  (sd(x,na.rm = T)/mean(x,na.rm = T))*100
}

#########################################################################
## VITAL RATE MODELS AND MPM SETUP 

source("/Users/maria/Dropbox/MSCA2019/theoretical_model/demographic_rate_functions.R")
# NULL MODEL

source("/Users/maria/Dropbox/MSCA2019/theoretical_model/parameters.R")

################  
### ONCE ALL VITAL RATES ARE PARAMETERIZED, CREATE A SIMULATION

# Joint model

##Set initial conditions

years=1050 # years of simulations

# simulate winter and summer environment from multivariate normal (in base: no covariance and relatively small variance)

cov=matrix(c(0.5,0,0,0.1),2,2)

set.seed(123)
env=mvrnorm(years, mu=c(-1,1), cov)

envW=env[,1]
envS=env[,2]

## data frame to hold results in
pdf(file="/Users/maria/Dropbox/MSCA2019/theoretical_model/plots/PreyPredCompManip.pdf", width = 14.5, height =13,pointsize=12)

for(k in 1:length(scen)){

  
  # Parameters  prey
  A.mean.vec <- c(scen.Sr_w[k],scen.D[k],scen.Sd[k],scen.B[k], scen.Snr_w[k],
                  scen.Sj_w[k],scen.Tj[k],scen.Tnr[k],scen.Tr[k],scen.Sr_s[k],scen.Snr_s[k])

  AB.diff.vec <- c(AB.scen.Sr_w[k]*abs(A.mean.vec[1]),0,AB.scen.Sd[k]*abs(A.mean.vec[3]),
                   0,AB.scen.Snr_w[k]*abs(A.mean.vec[5]),
                   AB.scen.Sj_w[k]*abs(A.mean.vec[6]),0,
                   0,0,
                   AB.scen.Sr_s[k]*abs(A.mean.vec[10]),AB.scen.Snr_s[k]*abs(A.mean.vec[11]))

  slope.env.vec <- c(env.scen.Sr_w[k]*abs(A.mean.vec[1]),env.scen.D[k]*abs(A.mean.vec[2]),env.scen.Sd[k]*abs(A.mean.vec[3]),
                     env.scen.B[k]*abs(A.mean.vec[4]),env.scen.Snr_w[k]*abs(A.mean.vec[5]),
                     env.scen.Sj_w[k]*abs(A.mean.vec[6]),env.scen.Tj[k]*abs(A.mean.vec[7]),
                     env.scen.Tnr[k]*abs(A.mean.vec[8]),env.scen.Tr[k]*abs(A.mean.vec[9]),
                     env.scen.Sr_s[k]*abs(A.mean.vec[10]),env.scen.Snr_s[k]*abs(A.mean.vec[11]))
  slope.dens.vec <- c(-0.02*abs(A.mean.vec[1]),0.02*abs(A.mean.vec[2]),-0.2*abs(A.mean.vec[3]),-0.03*abs(A.mean.vec[4]),
                      -0.02*abs(A.mean.vec[5]),0.1*abs(A.mean.vec[6]),-0.02*abs(A.mean.vec[7]),-0.1*abs(A.mean.vec[8]),-0.02*abs(A.mean.vec[9]),-0.01*abs(A.mean.vec[10]),-0.01*abs(A.mean.vec[11]))
  slope.bi.vec <- c(bi.scen.Sr_w[k]*abs(slope.dens.vec)[1],bi.scen.D[k]*abs(slope.dens.vec)[2],bi.scen.Sd[k]*abs(slope.dens.vec)[3],bi.scen.B[k]*abs(slope.dens.vec)[4],
                    bi.scen.Snr_w[k]*abs(slope.dens.vec)[5],bi.scen.Sj_w[k]*abs(slope.dens.vec)[6],bi.scen.Tj[k]*abs(slope.dens.vec)[7],bi.scen.Tnr[k]*abs(slope.dens.vec)[8],
                    bi.scen.Tr[k]*abs(slope.dens.vec)[9],bi.scen.Sr_s[k]*abs(slope.dens.vec)[10],bi.scen.Snr_s[k]*abs(slope.dens.vec)[11])
  
  slope.biS2.vec <- c(biS2.scen.Sr_w[k]*abs(slope.dens.vec)[1],biS2.scen.D[k]*abs(slope.dens.vec)[2],biS2.scen.Sd[k]*abs(slope.dens.vec)[3],biS2.scen.B[k]*abs(slope.dens.vec)[4],
                    biS2.scen.Snr_w[k]*abs(slope.dens.vec)[5],biS2.scen.Sj_w[k]*abs(slope.dens.vec)[6],biS2.scen.Tj[k]*abs(slope.dens.vec)[7],biS2.scen.Tnr[k]*abs(slope.dens.vec)[8],
                    biS2.scen.Tr[k]*abs(slope.dens.vec)[9],biS2.scen.Sr_s[k]*abs(slope.dens.vec)[10],biS2.scen.Snr_s[k]*abs(slope.dens.vec)[11])
  
  AB.env.diff.vec<- c(AB.env.scen.Sr_w[k]*abs(slope.env.vec[1]),AB.env.scen.D[k]*abs(slope.env.vec[2]),AB.env.scen.Sd[k]*abs(slope.env.vec[3]),
                      AB.env.scen.B[k]*abs(slope.env.vec[4]),AB.env.scen.Snr_w[k]*abs(slope.env.vec[5]),
                      AB.env.scen.Sj_w[k]*abs(slope.env.vec[6]),AB.env.scen.Tj[k]*abs(slope.env.vec[7]),
                      AB.env.scen.Tnr[k]*abs(slope.env.vec[8]),AB.env.scen.Tr[k]*abs(slope.env.vec[9]),
                      AB.env.scen.Sr_s[k]*abs(slope.env.vec[10]),AB.env.scen.Snr_s[k]*abs(slope.env.vec[11]))

  env.dens.diff.vec<- c(0.01,0.03,0.01,0.01,0.01,0.01,0,0,0,0.01,0.01) # density has + effect at bad env and - effect at good env
  dens.dens.bi.diff.vec<- c(dens.bi.scen.Sr_w[k]*abs(slope.dens.vec)[1],dens.bi.scen.D[k]*abs(slope.dens.vec)[2],dens.bi.scen.Sd[k]*abs(slope.dens.vec)[3],dens.bi.scen.B[k]*abs(slope.dens.vec)[4],
                            dens.bi.scen.Snr_w[k]*abs(slope.dens.vec)[5],dens.bi.scen.Sj_w[k]*abs(slope.dens.vec)[6],dens.bi.scen.Tj[k]*abs(slope.dens.vec)[7],dens.bi.scen.Tnr[k]*abs(slope.dens.vec)[8],
                            dens.bi.scen.Tr[k]*abs(slope.dens.vec)[9],dens.bi.scen.Sr_s[k]*abs(slope.dens.vec)[10],dens.bi.scen.Snr_s[k]*abs(slope.dens.vec)[11])
  AB.dens.bi.diff.vec<- rep(0,length(vr))

  # Parameters predator

  A.mean.vecS2 <- c(S2.scen.Sr_w[k],S2.scen.D[k],S2.scen.Sd[k],S2.scen.B[k], S2.scen.Snr_w[k],
                    S2.scen.Sj_w[k],S2.scen.Tj[k],S2.scen.Tnr[k],S2.scen.Tr[k],S2.scen.Sr_s[k],S2.scen.Snr_s[k])
  AB.diff.vecS2 <- c(S2.AB.scen.Sr_w[k]*abs(A.mean.vecS2[1]),0,S2.AB.scen.Sd[k]*abs(A.mean.vecS2[3]),
                   0,S2.AB.scen.Snr_w[k]*abs(A.mean.vecS2[5]),
                   S2.AB.scen.Sj_w[k]*abs(A.mean.vecS2[6]),0,
                   0,0,
                   S2.AB.scen.Sr_s[k]*abs(A.mean.vecS2[10]),S2.AB.scen.Snr_s[k]*abs(A.mean.vecS2[11]))


  slope.env.vecS2 <- c(S2.env.scen.Sr_w[k]*abs(A.mean.vecS2[1]),S2.env.scen.D[k]*abs(A.mean.vecS2[2]),S2.env.scen.Sd[k]*abs(A.mean.vecS2[3]),
                       S2.env.scen.B[k]*abs(A.mean.vecS2[4]),S2.env.scen.Snr_w[k]*abs(A.mean.vecS2[5]),
                       S2.env.scen.Sj_w[k]*abs(A.mean.vecS2[6]),S2.env.scen.Tj[k]*abs(A.mean.vecS2[7]),
                       S2.env.scen.Tnr[k]*abs(A.mean.vecS2[8]),S2.env.scen.Tr[k]*abs(A.mean.vecS2[9]),
                       S2.env.scen.Sr_s[k]*abs(A.mean.vecS2[10]),S2.env.scen.Snr_s[k]*abs(A.mean.vecS2[11]))
  slope.dens.vecS2 <- c(-0.05*abs(A.mean.vecS2[1]),5*abs(A.mean.vecS2[2]),-0.05*abs(A.mean.vecS2[3]),-0.05*abs(A.mean.vecS2[4]),
                        -0.05*abs(A.mean.vecS2[5]),-2*abs(A.mean.vecS2[6]),-0.01*abs(A.mean.vecS2[7]),-2*abs(A.mean.vecS2[8]),-0.01*abs(A.mean.vecS2[9]),-0.01*abs(A.mean.vecS2[10]),-0.01*abs(A.mean.vecS2[11]))

  slope.bi.vecS2 <-  c(S2.bi.scen.Sr_w[k]*abs(slope.dens.vecS2)[1],S2.bi.scen.D[k]*abs(slope.dens.vecS2)[2],S2.bi.scen.Sd[k]*abs(slope.dens.vecS2)[3],
                       S2.bi.scen.B[k]*abs(slope.dens.vecS2)[4],S2.bi.scen.Snr_w[k]*abs(slope.dens.vecS2)[5],S2.bi.scen.Sj_w[k]*abs(slope.dens.vecS2)[6],
                       S2.bi.scen.Tj[k]*abs(slope.dens.vecS2)[7],S2.bi.scen.Tnr[k]*abs(slope.dens.vecS2)[8],
                       S2.bi.scen.Tr[k]*abs(slope.dens.vecS2)[9],S2.bi.scen.Sr_s[k]*abs(slope.dens.vecS2)[10],S2.bi.scen.Snr_s[k]*abs(slope.dens.vecS2)[11])
  
  slope.biS2.vecS2 <-  c(S2.biS2.scen.Sr_w[k]*abs(slope.dens.vecS2)[1],S2.biS2.scen.D[k]*abs(slope.dens.vecS2)[2],S2.biS2.scen.Sd[k]*abs(slope.dens.vecS2)[3],
                       S2.biS2.scen.B[k]*abs(slope.dens.vecS2)[4],S2.biS2.scen.Snr_w[k]*abs(slope.dens.vecS2)[5],S2.biS2.scen.Sj_w[k]*abs(slope.dens.vecS2)[6],
                       S2.biS2.scen.Tj[k]*abs(slope.dens.vecS2)[7],S2.biS2.scen.Tnr[k]*abs(slope.dens.vecS2)[8],
                       S2.biS2.scen.Tr[k]*abs(slope.dens.vecS2)[9],S2.biS2.scen.Sr_s[k]*abs(slope.dens.vecS2)[10],S2.biS2.scen.Snr_s[k]*abs(slope.dens.vecS2)[11])
  
  slope.bi2.vecS2 <- c(-0.001,0,-0.001,S2.bi2.scen.B[k],-0.001,-0.001,-0.001,-0.001,-0.001,-0.001,-0.001)
  AB.env.diff.vecS2<- c(S2.AB.env.scen.Sr_w[k]*abs(slope.env.vecS2[1]),S2.AB.env.scen.D[k]*abs(slope.env.vecS2[2]),S2.AB.env.scen.Sd[k]*abs(slope.env.vecS2[3]),
                        S2.AB.env.scen.B[k]*abs(slope.env.vecS2[4]),S2.AB.env.scen.Snr_w[k]*abs(slope.env.vecS2[5]),
                        S2.AB.env.scen.Sj_w[k]*abs(slope.env.vecS2[6]),S2.AB.env.scen.Tj[k]*abs(slope.env.vecS2[7]),
                        S2.AB.env.scen.Tnr[k]*abs(slope.env.vecS2[8]),S2.AB.env.scen.Tr[k]*abs(slope.env.vecS2[9]),
                        S2.AB.env.scen.Sr_s[k]*abs(slope.env.vecS2[10]),S2.AB.env.scen.Snr_s[k]*abs(slope.env.vecS2[11]))

  dens.dens.bi.diff.vecS2<- c(S2.dens.bi.scen.Sr_w[k]*abs(slope.dens.vecS2)[1],S2.dens.bi.scen.D[k]*abs(slope.dens.vecS2)[2],S2.dens.bi.scen.Sd[k]*abs(slope.dens.vecS2)[3],
                              S2.dens.bi.scen.B[k]*abs(slope.dens.vecS2)[4],S2.dens.bi.scen.Snr_w[k]*abs(slope.dens.vecS2)[5],S2.dens.bi.scen.Sj_w[k]*abs(slope.dens.vecS2)[6],
                              S2.dens.bi.scen.Tj[k]*abs(slope.dens.vecS2)[7],S2.dens.bi.scen.Tnr[k]*abs(slope.dens.vecS2)[8],
                              S2.dens.bi.scen.Tr[k]*abs(slope.dens.vecS2)[9],S2.dens.bi.scen.Sr_s[k]*abs(slope.dens.vecS2)[10],S2.dens.bi.scen.Snr_s[k]*abs(slope.dens.vecS2)[11])
  env.dens.diff.vecS2<- c(0.01,0.01,0,0.01,0.01,0.01,0,0,0,0.01,0.01) # density has + effect at bad env and - effect at good env

  # Parameters predator 2
  
  A.mean.vecS3 <- c(S3.scen.Sr_w[k],S3.scen.D[k],S3.scen.Sd[k],S3.scen.B[k], S3.scen.Snr_w[k],
                    S3.scen.Sj_w[k],S3.scen.Tj[k],S3.scen.Tnr[k],S3.scen.Tr[k],S3.scen.Sr_s[k],S3.scen.Snr_s[k])
  AB.diff.vecS3 <- c(S3.AB.scen.Sr_w[k]*abs(A.mean.vecS3[1]),0,S3.AB.scen.Sd[k]*abs(A.mean.vecS3[3]),
                     0,S3.AB.scen.Snr_w[k]*abs(A.mean.vecS3[5]),
                     S3.AB.scen.Sj_w[k]*abs(A.mean.vecS3[6]),0,
                     0,0,
                     S3.AB.scen.Sr_s[k]*abs(A.mean.vecS3[10]),S3.AB.scen.Snr_s[k]*abs(A.mean.vecS3[11]))
  
  
  slope.env.vecS3 <- c(S3.env.scen.Sr_w[k]*abs(A.mean.vecS3[1]),S3.env.scen.D[k]*abs(A.mean.vecS3[2]),S3.env.scen.Sd[k]*abs(A.mean.vecS3[3]),
                       S3.env.scen.B[k]*abs(A.mean.vecS3[4]),S3.env.scen.Snr_w[k]*abs(A.mean.vecS3[5]),
                       S3.env.scen.Sj_w[k]*abs(A.mean.vecS3[6]),S3.env.scen.Tj[k]*abs(A.mean.vecS3[7]),
                       S3.env.scen.Tnr[k]*abs(A.mean.vecS3[8]),S3.env.scen.Tr[k]*abs(A.mean.vecS3[9]),
                       S3.env.scen.Sr_s[k]*abs(A.mean.vecS3[10]),S3.env.scen.Snr_s[k]*abs(A.mean.vecS3[11]))
  slope.dens.vecS3 <- c(-0.05*abs(A.mean.vecS3[1]),5*abs(A.mean.vecS3[2]),-0.05*abs(A.mean.vecS3[3]),-0.05*abs(A.mean.vecS3[4]),
                        -0.05*abs(A.mean.vecS3[5]),-2*abs(A.mean.vecS3[6]),-0.01*abs(A.mean.vecS3[7]),-2*abs(A.mean.vecS3[8]),-0.01*abs(A.mean.vecS3[9]),-0.01*abs(A.mean.vecS3[10]),-0.01*abs(A.mean.vecS3[11]))
  
  slope.bi.vecS3 <-  c(S3.bi.scen.Sr_w[k]*abs(slope.dens.vecS3)[1],S3.bi.scen.D[k]*abs(slope.dens.vecS3)[2],S3.bi.scen.Sd[k]*abs(slope.dens.vecS3)[3],
                       S3.bi.scen.B[k]*abs(slope.dens.vecS3)[4],S3.bi.scen.Snr_w[k]*abs(slope.dens.vecS3)[5],S3.bi.scen.Sj_w[k]*abs(slope.dens.vecS3)[6],
                       S3.bi.scen.Tj[k]*abs(slope.dens.vecS3)[7],S3.bi.scen.Tnr[k]*abs(slope.dens.vecS3)[8],
                       S3.bi.scen.Tr[k]*abs(slope.dens.vecS3)[9],S3.bi.scen.Sr_s[k]*abs(slope.dens.vecS3)[10],S3.bi.scen.Snr_s[k]*abs(slope.dens.vecS3)[11])
  
  slope.biS2.vecS3 <-  c(S3.biS2.scen.Sr_w[k]*abs(slope.dens.vecS2)[1],S3.biS2.scen.D[k]*abs(slope.dens.vecS2)[2],S3.biS2.scen.Sd[k]*abs(slope.dens.vecS2)[3],
                         S3.biS2.scen.B[k]*abs(slope.dens.vecS2)[4],S3.biS2.scen.Snr_w[k]*abs(slope.dens.vecS2)[5],S3.biS2.scen.Sj_w[k]*abs(slope.dens.vecS2)[6],
                         S3.biS2.scen.Tj[k]*abs(slope.dens.vecS2)[7],S3.biS2.scen.Tnr[k]*abs(slope.dens.vecS2)[8],
                         S3.biS2.scen.Tr[k]*abs(slope.dens.vecS2)[9],S3.biS2.scen.Sr_s[k]*abs(slope.dens.vecS2)[10],S3.biS2.scen.Snr_s[k]*abs(slope.dens.vecS2)[11])
  
  slope.bi2.vecS3 <- c(-0.001,0,-0.001,S3.bi2.scen.B[k],-0.001,-0.001,-0.001,-0.001,-0.001,-0.001,-0.001)
  
  AB.env.diff.vecS3<- c(S3.AB.env.scen.Sr_w[k]*abs(slope.env.vecS3[1]),S3.AB.env.scen.D[k]*abs(slope.env.vecS3[2]),S3.AB.env.scen.Sd[k]*abs(slope.env.vecS3[3]),
                        S3.AB.env.scen.B[k]*abs(slope.env.vecS3[4]),S3.AB.env.scen.Snr_w[k]*abs(slope.env.vecS3[5]),
                        S3.AB.env.scen.Sj_w[k]*abs(slope.env.vecS3[6]),S3.AB.env.scen.Tj[k]*abs(slope.env.vecS3[7]),
                        S3.AB.env.scen.Tnr[k]*abs(slope.env.vecS3[8]),S3.AB.env.scen.Tr[k]*abs(slope.env.vecS3[9]),
                        S3.AB.env.scen.Sr_s[k]*abs(slope.env.vecS3[10]),S3.AB.env.scen.Snr_s[k]*abs(slope.env.vecS3[11]))
  
  dens.dens.bi.diff.vecS3<- c(S3.dens.bi.scen.Sr_w[k],S3.dens.bi.scen.D[k],S3.dens.bi.scen.Sd[k],
                              S3.dens.bi.scen.B[k],S3.dens.bi.scen.Snr_w[k],
                              S3.dens.bi.scen.Sj_w[k],S3.dens.bi.scen.Tj[k],
                              S3.dens.bi.scen.Tnr[k],S3.dens.bi.scen.Tr[k],
                              S3.dens.bi.scen.Sr_s[k],S3.dens.bi.scen.Snr_s[k])
  
    
  env.dens.diff.vecS3<- c(0.01,0.01,0,0.01,0.01,0.01,0,0,0,0.01,0.01) # density has + effect at bad env and - effect at good env
  
  

  ## Initial prey

  densA=c(0,30,20)
  densB=c(0,30,20)

  densMet=c(densA,densB)

  dens.biA=sum(densA)
  dens.biB=sum(densB)

  # Initial predator 1

  S2_densA=c(0,15,7)
  S2_densB=c(0,15,7)

  S2_densMet=c(S2_densA,S2_densB)


  S2_dens.biA=sum(S2_densA)
  S2_dens.biB=sum(S2_densB)
  
  # Initial predator 2
  
  S3_densA=c(0,20,20)
  S3_densB=c(0,25,15)
  
  S3_densMet=c(S3_densA,S3_densB)
  
  
  S3_dens.biA=sum(S3_densA)
  S3_dens.biB=sum(S3_densB)

  output=data.frame(dens=c(densA,densB,S2_densA,S2_densB,S3_densA,S3_densB),
                    species=rep(c("Prey","Predator1", "Predator2"),each=6),
                    site=rep(rep(c("A","B"),each=3),3),
                    season="summer",
                    stage=rep(c("J","N","R"),6),
                    year=1)

  matsS1=array(0,c(3,3,50))
  matsS2=array(0,c(3,3,50))
  matsS3=array(0,c(3,3,50))

  for(y in 1:years){

    #summer

    #Prey
    summerA=matrix(c(0,0,0,
                     0,Snr_s("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA),0,
                     Sr_s("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*B("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA),0,Sr_s("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)),3,3,byrow = F)

    summerB=matrix(c(0,0,0,
                     0,Snr_s("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB),0,
                     Sr_s("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*B("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB),0,Sr_s("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)),3,3,byrow = F)



    Bs[1:3,1:3]=summerA
    Bs[4:6,4:6]=summerB


    Ms[3,3]=1-D("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)
    Ms[4,4]=1-D("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)
    Ms[6,3]=D("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Ds("A",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)
    Ms[5,4]=D("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Ds("B",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)



    densW=t(Ps)%*%Ms%*%Ps%*%Bs%*%densMet

    output=rbind(output,data.frame(dens=c(round(densW[1:3],2),round(densW[4:6],2)),
                                   species="Prey",
                                   site=rep(c("A","B"),each=3),
                                   season=rep("winter",6),stage=rep(c("J","N","R"),2),year=rep(y,6)))

    densA=densW[1:3]
    densB=densW[4:6]
    densMet=c(densA,densB)

    #predator1 
    S2_summerA=matrix(c(0,0,0,
                        0,S2_Snr_s("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA),0,
                        S2_Sr_s("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_B("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA),0,S2_Sr_s("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)),3,3,byrow = F)

    S2_summerB=matrix(c(0,0,0,
                        0,S2_Snr_s("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB),0,
                        S2_Sr_s("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)*S2_B("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB),0,S2_Sr_s("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)),3,3,byrow = F)



    S2_Bs[1:3,1:3]=S2_summerA
    S2_Bs[4:6,4:6]=S2_summerB


    S2_Ms[3,3]=1-S2_D("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)
    S2_Ms[4,4]=1-S2_D("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)
    S2_Ms[6,3]=S2_D("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Ds("A",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)
    S2_Ms[5,4]=S2_D("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)*S2_Ds("B",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)

    S2_densW=t(S2_Ps)%*%S2_Ms%*%S2_Ps%*%S2_Bs%*%S2_densMet

    output=rbind(output,data.frame(dens=c(round(S2_densW[1:3],2),round(S2_densW[4:6],2)),
                                   species="Predator1",
                                   site=rep(c("A","B"),each=3),
                                   season=rep("winter",6),stage=rep(c("J","N","R"),2),year=rep(y,6)))

    S2_densA=S2_densW[1:3]
    S2_densB=S2_densW[4:6]
    S2_densMet=c(S2_densA,S2_densB)

    #predator 2 
    S3_summerA=matrix(c(0,0,0,
                        0,S3_Snr_s("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA),0,
                        S3_Sr_s("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_B("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA),0,S3_Sr_s("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)),3,3,byrow = F)
    
    S3_summerB=matrix(c(0,0,0,
                        0,S3_Snr_s("B",envS[y],sum(S3_densB[1:3]),dens.biB,S2_dens.biB),0,
                        S3_Sr_s("B",envS[y],sum(S3_densB[1:3]),dens.biB,S2_dens.biB)*S3_B("B",envS[y],sum(S3_densB[1:3]),dens.biB,S2_dens.biB),0,S3_Sr_s("B",envS[y],sum(S3_densB[1:3]),dens.biB,S2_dens.biB)),3,3,byrow = F)
    
    
    
    S3_Bs[1:3,1:3]=S3_summerA
    S3_Bs[4:6,4:6]=S3_summerB
    
    
    S3_Ms[3,3]=1-S3_D("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)
    S3_Ms[4,4]=1-S3_D("B",envS[y],sum(S3_densB[1:3]),dens.biB,S2_dens.biB)
    S3_Ms[6,3]=S3_D("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Ds("A",envS[y],sum(S3_densB[1:3]),dens.biB,S2_dens.biB)
    S3_Ms[5,4]=S3_D("B",envS[y],sum(S3_densB[1:3]),dens.biB,S2_dens.biB)*S3_Ds("B",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)
    
    S3_densW=t(S3_Ps)%*%S3_Ms%*%S3_Ps%*%S3_Bs%*%S3_densMet
    
    output=rbind(output,data.frame(dens=c(round(S3_densW[1:3],2),round(S3_densW[4:6],2)),
                                   species="Predator2",
                                   site=rep(c("A","B"),each=3),
                                   season=rep("winter",6),stage=rep(c("J","N","R"),2),year=rep(y,6)))
    
    S3_densA=S3_densW[1:3]
    S3_densB=S3_densW[4:6]
    S3_densMet=c(S3_densA,S3_densB)
    
    ### Update bi densities

    dens.biA=sum(densA)
    dens.biB=sum(densB)

    S2_dens.biA=sum(S2_densA)
    S2_dens.biB=sum(S2_densB)

    S3_dens.biA=sum(S3_densA)
    S3_dens.biB=sum(S3_densB)
    
    # WINTER

    # Prey

    winterA=matrix(c(0,(1-Tj("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA))*Sj_w("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA),Tj("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Sj_w("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA),
                     0,(1-Tnr("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA))*Snr_w("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA),Tnr("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Snr_w("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA),
                     0, (1-Tr("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA))*Sr_w("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA),Tr("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Sr_w("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)),3,3,byrow = F)

    winterB=matrix(c(0,(1-Tj("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB))*Sj_w("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB),Tj("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Sj_w("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB),
                     0,(1-Tnr("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB))*Snr_w("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB),Tnr("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Snr_w("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB),
                     0, (1-Tr("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB))*Sr_w("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB),Tr("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Sr_w("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)),3,3,byrow = F)


    Bw[1:3,1:3]=winterA
    Bw[4:6,4:6]=winterB

    densS=t(Pw)%*%Mw%*%Pw%*%Bw%*%densMet

    output=rbind(output,data.frame(dens=c(round(densS[1:3],2),round(densS[4:6],2)),
                                   species="Prey",
                                   site=rep(c("A","B"),each=3),
                                   season=rep("summer",6),stage=rep(c("J","N","R"),2),year=rep(y+1,6)))

    # Predator 1
    S2_winterA=matrix(c(0,(1-S2_Tj("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA))*S2_Sj_w("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA),S2_Tj("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Sj_w("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA),
                        0,(1-S2_Tnr("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA))*S2_Snr_w("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA),S2_Tnr("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Snr_w("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA),
                        0, (1-S2_Tr("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA))*S2_Sr_w("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA),S2_Tr("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Sr_w("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)),3,3,byrow = F)

    S2_winterB=matrix(c(0,(1-S2_Tj("B",envW[y],sum(S2_densA[1:3]),dens.biB,S3_dens.biB))*S2_Sj_w("B",envW[y],sum(S2_densA[1:3]),dens.biB,S3_dens.biB),S2_Tj("B",envW[y],sum(S2_densA[1:3]),dens.biB,S3_dens.biB)*S2_Sj_w("B",envW[y],sum(S2_densA[1:3]),dens.biB,S3_dens.biB),
                        0,(1-S2_Tnr("B",envW[y],sum(S2_densA[1:3]),dens.biB,S3_dens.biB))*S2_Snr_w("B",envW[y],sum(S2_densA[1:3]),dens.biB,S3_dens.biB),S2_Tnr("B",envW[y],sum(S2_densA[1:3]),dens.biB,S3_dens.biB)*S2_Snr_w("B",envW[y],sum(S2_densA[1:3]),dens.biB,S3_dens.biB),
                        0, (1-S2_Tr("B",envW[y],sum(S2_densA[1:3]),dens.biB,S3_dens.biB))*S2_Sr_w("B",envW[y],sum(S2_densA[1:3]),dens.biB,S3_dens.biB),S2_Tr("B",envW[y],sum(S2_densA[1:3]),dens.biB,S3_dens.biB)*S2_Sr_w("B",envW[y],sum(S2_densA[1:3]),dens.biB,S3_dens.biB)),3,3,byrow = F)


    S2_Bw[1:3,1:3]=S2_winterA
    S2_Bw[4:6,4:6]=S2_winterB

    S2_densS=t(S2_Pw)%*%S2_Mw%*%S2_Pw%*%S2_Bw%*%S2_densMet

    output=rbind(output,data.frame(dens=c(round(S2_densS[1:3],2),round(S2_densS[4:6],2)),
                                   species="Predator1",
                                   site=rep(c("A","B"),each=3),
                                   season=rep("summer",6),stage=rep(c("J","N","R"),2),year=rep(y+1,6)))

    # Predator 2
    S3_winterA=matrix(c(0,(1-S3_Tj("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA))*S3_Sj_w("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA),S3_Tj("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Sj_w("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA),
                        0,(1-S3_Tnr("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA))*S3_Snr_w("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA),S3_Tnr("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Snr_w("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA),
                        0, (1-S3_Tr("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA))*S3_Sr_w("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA),S3_Tr("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Sr_w("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)),3,3,byrow = F)
    
    S3_winterB=matrix(c(0,(1-S3_Tj("B",envW[y],sum(S3_densA[1:3]),dens.biB,S3_dens.biB))*S3_Sj_w("B",envW[y],sum(S3_densA[1:3]),dens.biB,S3_dens.biB),S3_Tj("B",envW[y],sum(S3_densA[1:3]),dens.biB,S3_dens.biB)*S3_Sj_w("B",envW[y],sum(S3_densA[1:3]),dens.biB,S3_dens.biB),
                        0,(1-S3_Tnr("B",envW[y],sum(S3_densA[1:3]),dens.biB,S3_dens.biB))*S3_Snr_w("B",envW[y],sum(S3_densA[1:3]),dens.biB,S3_dens.biB),S3_Tnr("B",envW[y],sum(S3_densA[1:3]),dens.biB,S3_dens.biB)*S3_Snr_w("B",envW[y],sum(S3_densA[1:3]),dens.biB,S3_dens.biB),
                        0, (1-S3_Tr("B",envW[y],sum(S3_densA[1:3]),dens.biB,S3_dens.biB))*S3_Sr_w("B",envW[y],sum(S3_densA[1:3]),dens.biB,S3_dens.biB),S3_Tr("B",envW[y],sum(S3_densA[1:3]),dens.biB,S3_dens.biB)*S3_Sr_w("B",envW[y],sum(S3_densA[1:3]),dens.biB,S3_dens.biB)),3,3,byrow = F)
    
    
    S3_Bw[1:3,1:3]=S3_winterA
    S3_Bw[4:6,4:6]=S3_winterB
    
    S3_densS=t(S3_Pw)%*%S3_Mw%*%S3_Pw%*%S3_Bw%*%S3_densMet
    
    output=rbind(output,data.frame(dens=c(round(S3_densS[1:3],2),round(S3_densS[4:6],2)),
                                   species="Predator2",
                                   site=rep(c("A","B"),each=3),
                                   season=rep("summer",6),stage=rep(c("J","N","R"),2),year=rep(y+1,6)))
    
    if(y>1000){
      ## save MPMs to calculate life histories:

      ann.matS1A=matrix(c(0,(1-Tj("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA))*Sj_w("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Snr_s("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA),Tj("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Sj_w("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Sr_s("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA),
                          0,(1-Tnr("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA))*Snr_w("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Snr_s("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA),Tnr("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Snr_w("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Sr_s("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA),
                          Sr_w("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Sr_s("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*B("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA), (1-Tr("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA))*Sr_w("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Snr_s("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA),Tr("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Sr_w("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Sr_s("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)),3,3,byrow = F)

      ann.matS1B=matrix(c(0,(1-Tj("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB))*Sj_w("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Snr_s("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB),Tj("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Sj_w("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Sr_s("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB),
                          0,(1-Tnr("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB))*Snr_w("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Snr_s("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB),Tnr("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Snr_w("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Sr_s("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB),
                          Sr_w("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Sr_s("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*B("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB), (1-Tr("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB))*Sr_w("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Snr_s("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB),Tr("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Sr_w("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Sr_s("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)),3,3,byrow = F)

      ann.matS2A=matrix(c(0,(1-S2_Tj("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA))*S2_Sj_w("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Snr_s("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA),S2_Tj("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Sj_w("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Sr_s("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA),
                          0,(1-S2_Tnr("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA))*S2_Snr_w("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Snr_s("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA),S2_Tnr("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Snr_w("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Sr_s("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA),
                          S2_Sr_w("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Sr_s("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_B("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA), (1-S2_Tr("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA))*S2_Sr_w("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Snr_s("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA),S2_Tr("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Sr_w("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Sr_s("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)),3,3,byrow = F)

      ann.matS2B=matrix(c(0,(1-S2_Tj("B",envW[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB))*S2_Sj_w("B",envW[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)*S2_Snr_s("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB),S2_Tj("B",envW[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)*S2_Sj_w("B",envW[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)*S2_Sr_s("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB),
                          0,(1-S2_Tnr("B",envW[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB))*S2_Snr_w("B",envW[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)*S2_Snr_s("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB),S2_Tnr("B",envW[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)*S2_Snr_w("B",envW[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)*S2_Sr_s("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB),
                          S2_Sr_w("B",envW[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)*S2_Sr_s("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)*S2_B("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB), (1-S2_Tr("B",envW[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB))*S2_Sr_w("B",envW[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)*S2_Snr_s("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB),S2_Tr("B",envW[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)*S2_Sr_w("B",envW[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)*S2_Sr_s("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)),3,3,byrow = F)

      ann.matS3A=matrix(c(0,(1-S3_Tj("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA))*S3_Sj_w("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Snr_s("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA),S3_Tj("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Sj_w("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Sr_s("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA),
                          0,(1-S3_Tnr("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA))*S3_Snr_w("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Snr_s("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA),S3_Tnr("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Snr_w("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Sr_s("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA),
                          S3_Sr_w("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Sr_s("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_B("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA), (1-S3_Tr("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA))*S3_Sr_w("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Snr_s("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA),S3_Tr("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Sr_w("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Sr_s("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)),3,3,byrow = F)
      
      ann.matS3B=matrix(c(0,(1-S3_Tj("B",envW[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB))*S3_Sj_w("B",envW[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB)*S3_Snr_s("B",envS[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB),S3_Tj("B",envW[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB)*S3_Sj_w("B",envW[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB)*S3_Sr_s("B",envS[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB),
                          0,(1-S3_Tnr("B",envW[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB))*S3_Snr_w("B",envW[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB)*S3_Snr_s("B",envS[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB),S3_Tnr("B",envW[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB)*S3_Snr_w("B",envW[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB)*S3_Sr_s("B",envS[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB),
                          S3_Sr_w("B",envW[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB)*S3_Sr_s("B",envS[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB)*S3_B("B",envS[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB), (1-S3_Tr("B",envW[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB))*S3_Sr_w("B",envW[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB)*S3_Snr_s("B",envS[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB),S3_Tr("B",envW[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB)*S3_Sr_w("B",envW[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB)*S3_Sr_s("B",envS[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB)),3,3,byrow = F)
      
      matsS1[,,y-1000]=(ann.matS1A+ann.matS1B)/2
      matsS2[,,y-1000]=(ann.matS2A+ann.matS2B)/2
      matsS3[,,y-1000]=(ann.matS3A+ann.matS3B)/2
    }

    # Update densities
    densA=densS[1:3]
    densB=densS[4:6]
    densMet=c(densA,densB)

    S2_densA=S2_densS[1:3]
    S2_densB=S2_densS[4:6]
    S2_densMet=c(S2_densA,S2_densB)
    
    S3_densA=S3_densS[1:3]
    S3_densB=S3_densS[4:6]
    S3_densMet=c(S3_densA,S3_densB)

    ### Update bi densities

    dens.biA=sum(densA)
    dens.biB=sum(densB)

    S2_dens.biA=sum(S2_densA)
    S2_dens.biB=sum(S2_densB)
    
    S3_dens.biA=sum(S3_densA)
    S3_dens.biB=sum(S3_densB)


  }

  output$year[output$season%in%"winter"]=output$year[output$season%in%"winter"]+0.5
  output$int=paste(output$species,output$site)

  a=ggplot(output,aes(year,dens,col=species))+
    geom_line()+
    facet_grid(stage~site,scales = "free")+
    scale_color_manual(name="",values=c("#D55E00","orange","#009E73"))+
    xlab("Simulation year")+ylab("Total density")+theme_bw(base_size=20)+
    theme(panel.grid = element_blank())+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black"))+
    geom_vline(xintercept = 1000, col="red")+
    ggtitle(paste("Scenario",k))

  base=output[output$year>1000,]
  base$year=base$year-1000
  base$sim="base"

  b=ggplot(base,aes(year,dens,col=species))+
    geom_line()+
    facet_grid(stage~site,scales = "free")+
    scale_color_manual(name="",values=c("#D55E00","orange","#009E73"))+
    xlab("Simulation year")+ylab("Total density")+theme_bw(base_size=20)+
    theme(panel.grid = element_blank())+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black"))+
    ggtitle(scen[k])


  GT_S1=data.frame(gt=apply(matsS1,3,generation.time)[apply(matsS1,3,generation.time)<10],
                   species="Prey")
  GT_S2=data.frame(gt=apply(matsS2,3,generation.time)[apply(matsS2,3,generation.time)<10],
                   species="Predator1")
  
  GT_S3=data.frame(gt=apply(matsS3,3,generation.time)[apply(matsS3,3,generation.time)<10],
                   species="Predator2")

  GT=rbind(GT_S1,GT_S2,GT_S3)


  c=ggplot(GT,aes(gt,col=species,fill=species))+
    geom_histogram(aes(y=..density..), position="identity", alpha=0.3)+
    guides(fill=F,color=F)+
    scale_color_manual(name="",values=c("#D55E00","orange","#009E73"))+
    scale_fill_manual(name="",values=c("#D55E00","#E69F00","#009E73"))+
    xlab("Generation time")+theme_bw(base_size=20)+
    theme(panel.grid = element_blank())+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black"),
          legend.position = c(0.81,.83))

  NR_S1=data.frame(nr=apply(matsS1,3,net.reproductive.rate),
                   species="Prey")
  NR_S2=data.frame(nr=apply(matsS2,3,net.reproductive.rate),
                   species="Predator1")
  
  NR_S3=data.frame(nr=apply(matsS3,3,net.reproductive.rate),
                   species="Predator2")

  NR=rbind(NR_S1,NR_S2,NR_S3)


  d=ggplot(NR,aes(nr,col=species,fill=species))+
    geom_histogram(aes(y=..density..), position="identity", alpha=0.3)+

    scale_color_manual(name="",values=c("#D55E00","#E69F00","#009E73"))+
    scale_fill_manual(name="",values=c("#D55E00","#E69F00","#009E73"))+
    # guides(fill=F,color=F)+
    xlab("Net reproductive rate")+theme_bw(base_size=20)+
    theme(panel.grid = element_blank())+
    theme(panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.background = element_blank(),
          panel.border = element_rect(colour = "black"),
          legend.position = c(0.81,.83))


  print(( b / ( c + d ) ) )

}

dev.off()

#############################################################

## Base simulations

years=1050 # years of simulations
simul=100

# simulate winter and summer environment from multivariate normal
cov=matrix(c(0.5,0,0,0.1),2,2)

envir.sim=NULL
base.sim=NULL
base.sim.stage=NULL
remove.sim=NULL

init.dens=NULL # save starting densities for clim change simulations 

cor.ab=array(NA,c(length(scen),simul,3,3))

for(k in 1:length(scen)){
  
  print(paste("Running scenario ", k))
  
  # Parameters  prey
  A.mean.vec <- c(scen.Sr_w[k],scen.D[k],scen.Sd[k],scen.B[k], scen.Snr_w[k],
                  scen.Sj_w[k],scen.Tj[k],scen.Tnr[k],scen.Tr[k],scen.Sr_s[k],scen.Snr_s[k])
  
  AB.diff.vec <- c(AB.scen.Sr_w[k]*abs(A.mean.vec[1]),0,AB.scen.Sd[k]*abs(A.mean.vec[3]),
                   0,AB.scen.Snr_w[k]*abs(A.mean.vec[5]),
                   AB.scen.Sj_w[k]*abs(A.mean.vec[6]),0,
                   0,0,
                   AB.scen.Sr_s[k]*abs(A.mean.vec[10]),AB.scen.Snr_s[k]*abs(A.mean.vec[11]))
  
  slope.env.vec <- c(env.scen.Sr_w[k]*abs(A.mean.vec[1]),env.scen.D[k]*abs(A.mean.vec[2]),env.scen.Sd[k]*abs(A.mean.vec[3]),
                     env.scen.B[k]*abs(A.mean.vec[4]),env.scen.Snr_w[k]*abs(A.mean.vec[5]),
                     env.scen.Sj_w[k]*abs(A.mean.vec[6]),env.scen.Tj[k]*abs(A.mean.vec[7]),
                     env.scen.Tnr[k]*abs(A.mean.vec[8]),env.scen.Tr[k]*abs(A.mean.vec[9]),
                     env.scen.Sr_s[k]*abs(A.mean.vec[10]),env.scen.Snr_s[k]*abs(A.mean.vec[11]))
  slope.dens.vec <- c(-0.02*abs(A.mean.vec[1]),0.02*abs(A.mean.vec[2]),-0.2*abs(A.mean.vec[3]),-0.03*abs(A.mean.vec[4]),
                      -0.02*abs(A.mean.vec[5]),0.1*abs(A.mean.vec[6]),-0.02*abs(A.mean.vec[7]),-0.1*abs(A.mean.vec[8]),-0.02*abs(A.mean.vec[9]),-0.01*abs(A.mean.vec[10]),-0.01*abs(A.mean.vec[11]))
  slope.bi.vec <- c(bi.scen.Sr_w[k]*abs(slope.dens.vec)[1],bi.scen.D[k]*abs(slope.dens.vec)[2],bi.scen.Sd[k]*abs(slope.dens.vec)[3],bi.scen.B[k]*abs(slope.dens.vec)[4],
                    bi.scen.Snr_w[k]*abs(slope.dens.vec)[5],bi.scen.Sj_w[k]*abs(slope.dens.vec)[6],bi.scen.Tj[k]*abs(slope.dens.vec)[7],bi.scen.Tnr[k]*abs(slope.dens.vec)[8],
                    bi.scen.Tr[k]*abs(slope.dens.vec)[9],bi.scen.Sr_s[k]*abs(slope.dens.vec)[10],bi.scen.Snr_s[k]*abs(slope.dens.vec)[11])
  
  slope.biS2.vec <- c(biS2.scen.Sr_w[k]*abs(slope.dens.vec)[1],biS2.scen.D[k]*abs(slope.dens.vec)[2],biS2.scen.Sd[k]*abs(slope.dens.vec)[3],biS2.scen.B[k]*abs(slope.dens.vec)[4],
                      biS2.scen.Snr_w[k]*abs(slope.dens.vec)[5],biS2.scen.Sj_w[k]*abs(slope.dens.vec)[6],biS2.scen.Tj[k]*abs(slope.dens.vec)[7],biS2.scen.Tnr[k]*abs(slope.dens.vec)[8],
                      biS2.scen.Tr[k]*abs(slope.dens.vec)[9],biS2.scen.Sr_s[k]*abs(slope.dens.vec)[10],biS2.scen.Snr_s[k]*abs(slope.dens.vec)[11])
  
  AB.env.diff.vec<- c(AB.env.scen.Sr_w[k]*abs(slope.env.vec[1]),AB.env.scen.D[k]*abs(slope.env.vec[2]),AB.env.scen.Sd[k]*abs(slope.env.vec[3]),
                      AB.env.scen.B[k]*abs(slope.env.vec[4]),AB.env.scen.Snr_w[k]*abs(slope.env.vec[5]),
                      AB.env.scen.Sj_w[k]*abs(slope.env.vec[6]),AB.env.scen.Tj[k]*abs(slope.env.vec[7]),
                      AB.env.scen.Tnr[k]*abs(slope.env.vec[8]),AB.env.scen.Tr[k]*abs(slope.env.vec[9]),
                      AB.env.scen.Sr_s[k]*abs(slope.env.vec[10]),AB.env.scen.Snr_s[k]*abs(slope.env.vec[11]))
  
  env.dens.diff.vec<- c(0.01,0.03,0.01,0.01,0.01,0.01,0,0,0,0.01,0.01) # density has + effect at bad env and - effect at good env
  dens.dens.bi.diff.vec<- c(dens.bi.scen.Sr_w[k]*abs(slope.dens.vec)[1],dens.bi.scen.D[k]*abs(slope.dens.vec)[2],dens.bi.scen.Sd[k]*abs(slope.dens.vec)[3],dens.bi.scen.B[k]*abs(slope.dens.vec)[4],
                            dens.bi.scen.Snr_w[k]*abs(slope.dens.vec)[5],dens.bi.scen.Sj_w[k]*abs(slope.dens.vec)[6],dens.bi.scen.Tj[k]*abs(slope.dens.vec)[7],dens.bi.scen.Tnr[k]*abs(slope.dens.vec)[8],
                            dens.bi.scen.Tr[k]*abs(slope.dens.vec)[9],dens.bi.scen.Sr_s[k]*abs(slope.dens.vec)[10],dens.bi.scen.Snr_s[k]*abs(slope.dens.vec)[11])
  AB.dens.bi.diff.vec<- rep(0,length(vr))
  
  # Parameters predator
  
  A.mean.vecS2 <- c(S2.scen.Sr_w[k],S2.scen.D[k],S2.scen.Sd[k],S2.scen.B[k], S2.scen.Snr_w[k],
                    S2.scen.Sj_w[k],S2.scen.Tj[k],S2.scen.Tnr[k],S2.scen.Tr[k],S2.scen.Sr_s[k],S2.scen.Snr_s[k])
  AB.diff.vecS2 <- c(S2.AB.scen.Sr_w[k]*abs(A.mean.vecS2[1]),0,S2.AB.scen.Sd[k]*abs(A.mean.vecS2[3]),
                     0,S2.AB.scen.Snr_w[k]*abs(A.mean.vecS2[5]),
                     S2.AB.scen.Sj_w[k]*abs(A.mean.vecS2[6]),0,
                     0,0,
                     S2.AB.scen.Sr_s[k]*abs(A.mean.vecS2[10]),S2.AB.scen.Snr_s[k]*abs(A.mean.vecS2[11]))
  
  
  slope.env.vecS2 <- c(S2.env.scen.Sr_w[k]*abs(A.mean.vecS2[1]),S2.env.scen.D[k]*abs(A.mean.vecS2[2]),S2.env.scen.Sd[k]*abs(A.mean.vecS2[3]),
                       S2.env.scen.B[k]*abs(A.mean.vecS2[4]),S2.env.scen.Snr_w[k]*abs(A.mean.vecS2[5]),
                       S2.env.scen.Sj_w[k]*abs(A.mean.vecS2[6]),S2.env.scen.Tj[k]*abs(A.mean.vecS2[7]),
                       S2.env.scen.Tnr[k]*abs(A.mean.vecS2[8]),S2.env.scen.Tr[k]*abs(A.mean.vecS2[9]),
                       S2.env.scen.Sr_s[k]*abs(A.mean.vecS2[10]),S2.env.scen.Snr_s[k]*abs(A.mean.vecS2[11]))
  slope.dens.vecS2 <- c(-0.05*abs(A.mean.vecS2[1]),5*abs(A.mean.vecS2[2]),-0.05*abs(A.mean.vecS2[3]),-0.05*abs(A.mean.vecS2[4]),
                        -0.05*abs(A.mean.vecS2[5]),-2*abs(A.mean.vecS2[6]),-0.01*abs(A.mean.vecS2[7]),-2*abs(A.mean.vecS2[8]),-0.01*abs(A.mean.vecS2[9]),-0.01*abs(A.mean.vecS2[10]),-0.01*abs(A.mean.vecS2[11]))
  
  slope.bi.vecS2 <-  c(S2.bi.scen.Sr_w[k]*abs(slope.dens.vecS2)[1],S2.bi.scen.D[k]*abs(slope.dens.vecS2)[2],S2.bi.scen.Sd[k]*abs(slope.dens.vecS2)[3],
                       S2.bi.scen.B[k]*abs(slope.dens.vecS2)[4],S2.bi.scen.Snr_w[k]*abs(slope.dens.vecS2)[5],S2.bi.scen.Sj_w[k]*abs(slope.dens.vecS2)[6],
                       S2.bi.scen.Tj[k]*abs(slope.dens.vecS2)[7],S2.bi.scen.Tnr[k]*abs(slope.dens.vecS2)[8],
                       S2.bi.scen.Tr[k]*abs(slope.dens.vecS2)[9],S2.bi.scen.Sr_s[k]*abs(slope.dens.vecS2)[10],S2.bi.scen.Snr_s[k]*abs(slope.dens.vecS2)[11])
  
  slope.biS2.vecS2 <-  c(S2.biS2.scen.Sr_w[k]*abs(slope.dens.vecS2)[1],S2.biS2.scen.D[k]*abs(slope.dens.vecS2)[2],S2.biS2.scen.Sd[k]*abs(slope.dens.vecS2)[3],
                         S2.biS2.scen.B[k]*abs(slope.dens.vecS2)[4],S2.biS2.scen.Snr_w[k]*abs(slope.dens.vecS2)[5],S2.biS2.scen.Sj_w[k]*abs(slope.dens.vecS2)[6],
                         S2.biS2.scen.Tj[k]*abs(slope.dens.vecS2)[7],S2.biS2.scen.Tnr[k]*abs(slope.dens.vecS2)[8],
                         S2.biS2.scen.Tr[k]*abs(slope.dens.vecS2)[9],S2.biS2.scen.Sr_s[k]*abs(slope.dens.vecS2)[10],S2.biS2.scen.Snr_s[k]*abs(slope.dens.vecS2)[11])
  
  slope.bi2.vecS2 <- c(-0.001,0,-0.001,S2.bi2.scen.B[k],-0.001,-0.001,-0.001,-0.001,-0.001,-0.001,-0.001)
  AB.env.diff.vecS2<- c(S2.AB.env.scen.Sr_w[k]*abs(slope.env.vecS2[1]),S2.AB.env.scen.D[k]*abs(slope.env.vecS2[2]),S2.AB.env.scen.Sd[k]*abs(slope.env.vecS2[3]),
                        S2.AB.env.scen.B[k]*abs(slope.env.vecS2[4]),S2.AB.env.scen.Snr_w[k]*abs(slope.env.vecS2[5]),
                        S2.AB.env.scen.Sj_w[k]*abs(slope.env.vecS2[6]),S2.AB.env.scen.Tj[k]*abs(slope.env.vecS2[7]),
                        S2.AB.env.scen.Tnr[k]*abs(slope.env.vecS2[8]),S2.AB.env.scen.Tr[k]*abs(slope.env.vecS2[9]),
                        S2.AB.env.scen.Sr_s[k]*abs(slope.env.vecS2[10]),S2.AB.env.scen.Snr_s[k]*abs(slope.env.vecS2[11]))
  
  dens.dens.bi.diff.vecS2<- c(S2.dens.bi.scen.Sr_w[k]*abs(slope.dens.vecS2)[1],S2.dens.bi.scen.D[k]*abs(slope.dens.vecS2)[2],S2.dens.bi.scen.Sd[k]*abs(slope.dens.vecS2)[3],
                              S2.dens.bi.scen.B[k]*abs(slope.dens.vecS2)[4],S2.dens.bi.scen.Snr_w[k]*abs(slope.dens.vecS2)[5],S2.dens.bi.scen.Sj_w[k]*abs(slope.dens.vecS2)[6],
                              S2.dens.bi.scen.Tj[k]*abs(slope.dens.vecS2)[7],S2.dens.bi.scen.Tnr[k]*abs(slope.dens.vecS2)[8],
                              S2.dens.bi.scen.Tr[k]*abs(slope.dens.vecS2)[9],S2.dens.bi.scen.Sr_s[k]*abs(slope.dens.vecS2)[10],S2.dens.bi.scen.Snr_s[k]*abs(slope.dens.vecS2)[11])
  env.dens.diff.vecS2<- c(0.01,0.01,0,0.01,0.01,0.01,0,0,0,0.01,0.01) # density has + effect at bad env and - effect at good env
  
  # Parameters predator 2
  
  A.mean.vecS3 <- c(S3.scen.Sr_w[k],S3.scen.D[k],S3.scen.Sd[k],S3.scen.B[k], S3.scen.Snr_w[k],
                    S3.scen.Sj_w[k],S3.scen.Tj[k],S3.scen.Tnr[k],S3.scen.Tr[k],S3.scen.Sr_s[k],S3.scen.Snr_s[k])
  AB.diff.vecS3 <- c(S3.AB.scen.Sr_w[k]*abs(A.mean.vecS3[1]),0,S3.AB.scen.Sd[k]*abs(A.mean.vecS3[3]),
                     0,S3.AB.scen.Snr_w[k]*abs(A.mean.vecS3[5]),
                     S3.AB.scen.Sj_w[k]*abs(A.mean.vecS3[6]),0,
                     0,0,
                     S3.AB.scen.Sr_s[k]*abs(A.mean.vecS3[10]),S3.AB.scen.Snr_s[k]*abs(A.mean.vecS3[11]))
  
  
  slope.env.vecS3 <- c(S3.env.scen.Sr_w[k]*abs(A.mean.vecS3[1]),S3.env.scen.D[k]*abs(A.mean.vecS3[2]),S3.env.scen.Sd[k]*abs(A.mean.vecS3[3]),
                       S3.env.scen.B[k]*abs(A.mean.vecS3[4]),S3.env.scen.Snr_w[k]*abs(A.mean.vecS3[5]),
                       S3.env.scen.Sj_w[k]*abs(A.mean.vecS3[6]),S3.env.scen.Tj[k]*abs(A.mean.vecS3[7]),
                       S3.env.scen.Tnr[k]*abs(A.mean.vecS3[8]),S3.env.scen.Tr[k]*abs(A.mean.vecS3[9]),
                       S3.env.scen.Sr_s[k]*abs(A.mean.vecS3[10]),S3.env.scen.Snr_s[k]*abs(A.mean.vecS3[11]))
  slope.dens.vecS3 <- c(-0.05*abs(A.mean.vecS3[1]),5*abs(A.mean.vecS3[2]),-0.05*abs(A.mean.vecS3[3]),-0.05*abs(A.mean.vecS3[4]),
                        -0.05*abs(A.mean.vecS3[5]),-2*abs(A.mean.vecS3[6]),-0.01*abs(A.mean.vecS3[7]),-2*abs(A.mean.vecS3[8]),-0.01*abs(A.mean.vecS3[9]),-0.01*abs(A.mean.vecS3[10]),-0.01*abs(A.mean.vecS3[11]))
  
  slope.bi.vecS3 <-  c(S3.bi.scen.Sr_w[k]*abs(slope.dens.vecS3)[1],S3.bi.scen.D[k]*abs(slope.dens.vecS3)[2],S3.bi.scen.Sd[k]*abs(slope.dens.vecS3)[3],
                       S3.bi.scen.B[k]*abs(slope.dens.vecS3)[4],S3.bi.scen.Snr_w[k]*abs(slope.dens.vecS3)[5],S3.bi.scen.Sj_w[k]*abs(slope.dens.vecS3)[6],
                       S3.bi.scen.Tj[k]*abs(slope.dens.vecS3)[7],S3.bi.scen.Tnr[k]*abs(slope.dens.vecS3)[8],
                       S3.bi.scen.Tr[k]*abs(slope.dens.vecS3)[9],S3.bi.scen.Sr_s[k]*abs(slope.dens.vecS3)[10],S3.bi.scen.Snr_s[k]*abs(slope.dens.vecS3)[11])
  
  slope.biS2.vecS3 <-  c(S3.biS2.scen.Sr_w[k]*abs(slope.dens.vecS2)[1],S3.biS2.scen.D[k]*abs(slope.dens.vecS2)[2],S3.biS2.scen.Sd[k]*abs(slope.dens.vecS2)[3],
                         S3.biS2.scen.B[k]*abs(slope.dens.vecS2)[4],S3.biS2.scen.Snr_w[k]*abs(slope.dens.vecS2)[5],S3.biS2.scen.Sj_w[k]*abs(slope.dens.vecS2)[6],
                         S3.biS2.scen.Tj[k]*abs(slope.dens.vecS2)[7],S3.biS2.scen.Tnr[k]*abs(slope.dens.vecS2)[8],
                         S3.biS2.scen.Tr[k]*abs(slope.dens.vecS2)[9],S3.biS2.scen.Sr_s[k]*abs(slope.dens.vecS2)[10],S3.biS2.scen.Snr_s[k]*abs(slope.dens.vecS2)[11])
  
  slope.bi2.vecS3 <- c(-0.001,0,-0.001,S3.bi2.scen.B[k],-0.001,-0.001,-0.001,-0.001,-0.001,-0.001,-0.001)
  
  AB.env.diff.vecS3<- c(S3.AB.env.scen.Sr_w[k]*abs(slope.env.vecS3[1]),S3.AB.env.scen.D[k]*abs(slope.env.vecS3[2]),S3.AB.env.scen.Sd[k]*abs(slope.env.vecS3[3]),
                        S3.AB.env.scen.B[k]*abs(slope.env.vecS3[4]),S3.AB.env.scen.Snr_w[k]*abs(slope.env.vecS3[5]),
                        S3.AB.env.scen.Sj_w[k]*abs(slope.env.vecS3[6]),S3.AB.env.scen.Tj[k]*abs(slope.env.vecS3[7]),
                        S3.AB.env.scen.Tnr[k]*abs(slope.env.vecS3[8]),S3.AB.env.scen.Tr[k]*abs(slope.env.vecS3[9]),
                        S3.AB.env.scen.Sr_s[k]*abs(slope.env.vecS3[10]),S3.AB.env.scen.Snr_s[k]*abs(slope.env.vecS3[11]))
  
  dens.dens.bi.diff.vecS3<- c(S3.dens.bi.scen.Sr_w[k],S3.dens.bi.scen.D[k],S3.dens.bi.scen.Sd[k],
                              S3.dens.bi.scen.B[k],S3.dens.bi.scen.Snr_w[k],
                              S3.dens.bi.scen.Sj_w[k],S3.dens.bi.scen.Tj[k],
                              S3.dens.bi.scen.Tnr[k],S3.dens.bi.scen.Tr[k],
                              S3.dens.bi.scen.Sr_s[k],S3.dens.bi.scen.Snr_s[k])
  
  
  env.dens.diff.vecS3<- c(0.01,0.01,0,0.01,0.01,0.01,0,0,0,0.01,0.01) # density has + effect at bad env and - effect at good env
  
  
  for(x in 1:simul){
    
    env=mvrnorm(years, mu=c(-1,1), cov)
    
    envW=env[,1]
    envS=env[,2]
    
    # Initial prey 
    densA=c(0,30,20)
    densB=c(0,30,20)
    
    densMet=c(densA,densB)
    
    dens.biA=sum(densA)
    dens.biB=sum(densB)
    
    # Initial predator 1
    
    S2_densA=c(0,15,7)
    S2_densB=c(0,15,7)
    
    S2_densMet=c(S2_densA,S2_densB)
    
    
    S2_dens.biA=sum(S2_densA)
    S2_dens.biB=sum(S2_densB)
    
    # Initial predator 2
    
    S3_densA=c(0,20,20)
    S3_densB=c(0,25,15)
    
    S3_densMet=c(S3_densA,S3_densB)
    
    
    S3_dens.biA=sum(S3_densA)
    S3_dens.biB=sum(S3_densB)
    
   temp=data.frame(dens=c(densA,densB,S2_densA,S2_densB,S3_densA,S3_densB),
                      species=rep(c("Prey","Predator1", "Predator2"),each=6),
                      site=rep(rep(c("A","B"),each=3),3),
                      season="summer",
                      stage=rep(c("J","N","R"),6),
                      year=1)
    
    matsS1=array(0,c(3,3,50))
    matsS2=array(0,c(3,3,50))
    matsS3=array(0,c(3,3,50))
   
    
    for(y in 1:years){
      
      #summer
      
      #Prey
      summerA=matrix(c(0,0,0,
                       0,Snr_s("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA),0,
                       Sr_s("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*B("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA),0,Sr_s("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)),3,3,byrow = F)
      
      summerB=matrix(c(0,0,0,
                       0,Snr_s("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB),0,
                       Sr_s("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*B("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB),0,Sr_s("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)),3,3,byrow = F)
      
      
      
      Bs[1:3,1:3]=summerA
      Bs[4:6,4:6]=summerB
      
      
      Ms[3,3]=1-D("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)
      Ms[4,4]=1-D("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)
      Ms[6,3]=D("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Ds("A",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)
      Ms[5,4]=D("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Ds("B",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)
      
      
      
      densW=t(Ps)%*%Ms%*%Ps%*%Bs%*%densMet
      
      if(y>1000){
        
        temp=rbind(temp,data.frame(dens=c(round(densW[1:3],2),round(densW[4:6],2)),
                                   species="Prey",
                                   site=rep(c("A","B"),each=3),
                                   season=rep("winter",6),stage=rep(c("J","N","R"),2),year=rep(y,6)))
      }
      
      
      densA=densW[1:3]
      densB=densW[4:6]
      densMet=c(densA,densB)
      
      #predator1 
      S2_summerA=matrix(c(0,0,0,
                          0,S2_Snr_s("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA),0,
                          S2_Sr_s("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_B("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA),0,S2_Sr_s("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)),3,3,byrow = F)
      
      S2_summerB=matrix(c(0,0,0,
                          0,S2_Snr_s("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB),0,
                          S2_Sr_s("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)*S2_B("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB),0,S2_Sr_s("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)),3,3,byrow = F)
      
      
      
      S2_Bs[1:3,1:3]=S2_summerA
      S2_Bs[4:6,4:6]=S2_summerB
      
      
      S2_Ms[3,3]=1-S2_D("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)
      S2_Ms[4,4]=1-S2_D("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)
      S2_Ms[6,3]=S2_D("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Ds("A",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)
      S2_Ms[5,4]=S2_D("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)*S2_Ds("B",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)
      
      S2_densW=t(S2_Ps)%*%S2_Ms%*%S2_Ps%*%S2_Bs%*%S2_densMet
      
      if(y>1000){
        temp=rbind(temp,data.frame(dens=c(round(S2_densW[1:3],2),round(S2_densW[4:6],2)),
                                   species="Predator1",
                                   site=rep(c("A","B"),each=3),
                                   season=rep("winter",6),stage=rep(c("J","N","R"),2),year=rep(y,6)))

      }
      
      S2_densA=S2_densW[1:3]
      S2_densB=S2_densW[4:6]
      S2_densMet=c(S2_densA,S2_densB)
      
      #predator 2 
      S3_summerA=matrix(c(0,0,0,
                          0,S3_Snr_s("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA),0,
                          S3_Sr_s("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_B("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA),0,S3_Sr_s("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)),3,3,byrow = F)
      
      S3_summerB=matrix(c(0,0,0,
                          0,S3_Snr_s("B",envS[y],sum(S3_densB[1:3]),dens.biB,S2_dens.biB),0,
                          S3_Sr_s("B",envS[y],sum(S3_densB[1:3]),dens.biB,S2_dens.biB)*S3_B("B",envS[y],sum(S3_densB[1:3]),dens.biB,S2_dens.biB),0,S3_Sr_s("B",envS[y],sum(S3_densB[1:3]),dens.biB,S2_dens.biB)),3,3,byrow = F)
      
      
      
      S3_Bs[1:3,1:3]=S3_summerA
      S3_Bs[4:6,4:6]=S3_summerB
      
      
      S3_Ms[3,3]=1-S3_D("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)
      S3_Ms[4,4]=1-S3_D("B",envS[y],sum(S3_densB[1:3]),dens.biB,S2_dens.biB)
      S3_Ms[6,3]=S3_D("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Ds("A",envS[y],sum(S3_densB[1:3]),dens.biB,S2_dens.biB)
      S3_Ms[5,4]=S3_D("B",envS[y],sum(S3_densB[1:3]),dens.biB,S2_dens.biB)*S3_Ds("B",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)
      
      S3_densW=t(S3_Ps)%*%S3_Ms%*%S3_Ps%*%S3_Bs%*%S3_densMet
      
      if(y>1000){
        temp=rbind(temp,data.frame(dens=c(round(S3_densW[1:3],2),round(S3_densW[4:6],2)),
                                   species="Predator2",
                                   site=rep(c("A","B"),each=3),
                                   season=rep("winter",6),stage=rep(c("J","N","R"),2),year=rep(y,6)))
        
      }
      
      S3_densA=S3_densW[1:3]
      S3_densB=S3_densW[4:6]
      S3_densMet=c(S3_densA,S3_densB)
      
      ### Update bi densities
      
      dens.biA=sum(densA)
      dens.biB=sum(densB)
      
      S2_dens.biA=sum(S2_densA)
      S2_dens.biB=sum(S2_densB)
      
      S3_dens.biA=sum(S3_densA)
      S3_dens.biB=sum(S3_densB)
      
      # WINTER
      
      # Prey
      
      winterA=matrix(c(0,(1-Tj("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA))*Sj_w("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA),Tj("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Sj_w("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA),
                       0,(1-Tnr("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA))*Snr_w("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA),Tnr("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Snr_w("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA),
                       0, (1-Tr("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA))*Sr_w("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA),Tr("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Sr_w("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)),3,3,byrow = F)
      
      winterB=matrix(c(0,(1-Tj("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB))*Sj_w("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB),Tj("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Sj_w("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB),
                       0,(1-Tnr("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB))*Snr_w("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB),Tnr("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Snr_w("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB),
                       0, (1-Tr("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB))*Sr_w("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB),Tr("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Sr_w("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)),3,3,byrow = F)
      
      
      Bw[1:3,1:3]=winterA
      Bw[4:6,4:6]=winterB
      
      densS=t(Pw)%*%Mw%*%Pw%*%Bw%*%densMet
      
      if(y>1000){
        
        temp=rbind(temp,data.frame(dens=c(round(densS[1:3],2),round(densS[4:6],2)),
                                   species="Prey",
                                   site=rep(c("A","B"),each=3),
                                   season=rep("summer",6),stage=rep(c("J","N","R"),2),year=rep(y+1,6)))
        
      }
      
      
      
      # Predator 1
      S2_winterA=matrix(c(0,(1-S2_Tj("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA))*S2_Sj_w("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA),S2_Tj("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Sj_w("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA),
                          0,(1-S2_Tnr("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA))*S2_Snr_w("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA),S2_Tnr("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Snr_w("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA),
                          0, (1-S2_Tr("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA))*S2_Sr_w("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA),S2_Tr("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Sr_w("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)),3,3,byrow = F)
      
      S2_winterB=matrix(c(0,(1-S2_Tj("B",envW[y],sum(S2_densA[1:3]),dens.biB,S3_dens.biB))*S2_Sj_w("B",envW[y],sum(S2_densA[1:3]),dens.biB,S3_dens.biB),S2_Tj("B",envW[y],sum(S2_densA[1:3]),dens.biB,S3_dens.biB)*S2_Sj_w("B",envW[y],sum(S2_densA[1:3]),dens.biB,S3_dens.biB),
                          0,(1-S2_Tnr("B",envW[y],sum(S2_densA[1:3]),dens.biB,S3_dens.biB))*S2_Snr_w("B",envW[y],sum(S2_densA[1:3]),dens.biB,S3_dens.biB),S2_Tnr("B",envW[y],sum(S2_densA[1:3]),dens.biB,S3_dens.biB)*S2_Snr_w("B",envW[y],sum(S2_densA[1:3]),dens.biB,S3_dens.biB),
                          0, (1-S2_Tr("B",envW[y],sum(S2_densA[1:3]),dens.biB,S3_dens.biB))*S2_Sr_w("B",envW[y],sum(S2_densA[1:3]),dens.biB,S3_dens.biB),S2_Tr("B",envW[y],sum(S2_densA[1:3]),dens.biB,S3_dens.biB)*S2_Sr_w("B",envW[y],sum(S2_densA[1:3]),dens.biB,S3_dens.biB)),3,3,byrow = F)
      
      
      S2_Bw[1:3,1:3]=S2_winterA
      S2_Bw[4:6,4:6]=S2_winterB
      
      S2_densS=t(S2_Pw)%*%S2_Mw%*%S2_Pw%*%S2_Bw%*%S2_densMet
      
      if(y>1000){
        
        temp=rbind(temp,data.frame(dens=c(round(S2_densS[1:3],2),round(S2_densS[4:6],2)),
                                   species="Predator1",
                                   site=rep(c("A","B"),each=3),
                                   season=rep("summer",6),stage=rep(c("J","N","R"),2),year=rep(y+1,6)))
        
      }
      
      # Predator 2
      S3_winterA=matrix(c(0,(1-S3_Tj("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA))*S3_Sj_w("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA),S3_Tj("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Sj_w("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA),
                          0,(1-S3_Tnr("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA))*S3_Snr_w("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA),S3_Tnr("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Snr_w("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA),
                          0, (1-S3_Tr("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA))*S3_Sr_w("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA),S3_Tr("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Sr_w("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)),3,3,byrow = F)
      
      S3_winterB=matrix(c(0,(1-S3_Tj("B",envW[y],sum(S3_densA[1:3]),dens.biB,S3_dens.biB))*S3_Sj_w("B",envW[y],sum(S3_densA[1:3]),dens.biB,S3_dens.biB),S3_Tj("B",envW[y],sum(S3_densA[1:3]),dens.biB,S3_dens.biB)*S3_Sj_w("B",envW[y],sum(S3_densA[1:3]),dens.biB,S3_dens.biB),
                          0,(1-S3_Tnr("B",envW[y],sum(S3_densA[1:3]),dens.biB,S3_dens.biB))*S3_Snr_w("B",envW[y],sum(S3_densA[1:3]),dens.biB,S3_dens.biB),S3_Tnr("B",envW[y],sum(S3_densA[1:3]),dens.biB,S3_dens.biB)*S3_Snr_w("B",envW[y],sum(S3_densA[1:3]),dens.biB,S3_dens.biB),
                          0, (1-S3_Tr("B",envW[y],sum(S3_densA[1:3]),dens.biB,S3_dens.biB))*S3_Sr_w("B",envW[y],sum(S3_densA[1:3]),dens.biB,S3_dens.biB),S3_Tr("B",envW[y],sum(S3_densA[1:3]),dens.biB,S3_dens.biB)*S3_Sr_w("B",envW[y],sum(S3_densA[1:3]),dens.biB,S3_dens.biB)),3,3,byrow = F)
      
      
      S3_Bw[1:3,1:3]=S3_winterA
      S3_Bw[4:6,4:6]=S3_winterB
      
      S3_densS=t(S3_Pw)%*%S3_Mw%*%S3_Pw%*%S3_Bw%*%S3_densMet
      
      if(y>1000){
        
        temp=rbind(temp,data.frame(dens=c(round(S3_densS[1:3],2),round(S3_densS[4:6],2)),
                                   species="Predator2",
                                   site=rep(c("A","B"),each=3),
                                   season=rep("summer",6),stage=rep(c("J","N","R"),2),year=rep(y+1,6)))
        
        
        if(y==1001){
          
          init.dens=rbind(init.dens,data.frame(densA=round(densS[1:3],2),densB=round(densS[4:6],2),
                                               S2_densA=round(S2_densS[1:3],2),S2_densB=round(S2_densS[4:6],2),
                                               S3_densA=round(S3_densS[1:3],2),S3_densB=round(S3_densS[4:6],2),
                                               sim=x,
                                               scen=k))
        }
        ## save MPMs to calculate life histories:
        
        ann.matS1A=matrix(c(0,(1-Tj("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA))*Sj_w("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Snr_s("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA),Tj("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Sj_w("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Sr_s("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA),
                            0,(1-Tnr("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA))*Snr_w("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Snr_s("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA),Tnr("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Snr_w("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Sr_s("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA),
                            Sr_w("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Sr_s("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*B("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA), (1-Tr("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA))*Sr_w("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Snr_s("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA),Tr("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Sr_w("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Sr_s("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)),3,3,byrow = F)
        
        ann.matS1B=matrix(c(0,(1-Tj("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB))*Sj_w("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Snr_s("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB),Tj("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Sj_w("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Sr_s("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB),
                            0,(1-Tnr("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB))*Snr_w("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Snr_s("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB),Tnr("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Snr_w("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Sr_s("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB),
                            Sr_w("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Sr_s("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*B("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB), (1-Tr("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB))*Sr_w("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Snr_s("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB),Tr("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Sr_w("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Sr_s("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)),3,3,byrow = F)
        
        ann.matS2A=matrix(c(0,(1-S2_Tj("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA))*S2_Sj_w("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Snr_s("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA),S2_Tj("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Sj_w("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Sr_s("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA),
                            0,(1-S2_Tnr("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA))*S2_Snr_w("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Snr_s("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA),S2_Tnr("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Snr_w("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Sr_s("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA),
                            S2_Sr_w("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Sr_s("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_B("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA), (1-S2_Tr("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA))*S2_Sr_w("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Snr_s("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA),S2_Tr("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Sr_w("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Sr_s("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)),3,3,byrow = F)
        
        ann.matS2B=matrix(c(0,(1-S2_Tj("B",envW[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB))*S2_Sj_w("B",envW[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)*S2_Snr_s("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB),S2_Tj("B",envW[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)*S2_Sj_w("B",envW[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)*S2_Sr_s("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB),
                            0,(1-S2_Tnr("B",envW[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB))*S2_Snr_w("B",envW[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)*S2_Snr_s("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB),S2_Tnr("B",envW[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)*S2_Snr_w("B",envW[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)*S2_Sr_s("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB),
                            S2_Sr_w("B",envW[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)*S2_Sr_s("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)*S2_B("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB), (1-S2_Tr("B",envW[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB))*S2_Sr_w("B",envW[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)*S2_Snr_s("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB),S2_Tr("B",envW[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)*S2_Sr_w("B",envW[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)*S2_Sr_s("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)),3,3,byrow = F)
        
        ann.matS3A=matrix(c(0,(1-S3_Tj("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA))*S3_Sj_w("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Snr_s("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA),S3_Tj("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Sj_w("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Sr_s("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA),
                            0,(1-S3_Tnr("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA))*S3_Snr_w("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Snr_s("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA),S3_Tnr("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Snr_w("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Sr_s("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA),
                            S3_Sr_w("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Sr_s("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_B("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA), (1-S3_Tr("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA))*S3_Sr_w("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Snr_s("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA),S3_Tr("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Sr_w("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Sr_s("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)),3,3,byrow = F)
        
        ann.matS3B=matrix(c(0,(1-S3_Tj("B",envW[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB))*S3_Sj_w("B",envW[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB)*S3_Snr_s("B",envS[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB),S3_Tj("B",envW[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB)*S3_Sj_w("B",envW[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB)*S3_Sr_s("B",envS[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB),
                            0,(1-S3_Tnr("B",envW[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB))*S3_Snr_w("B",envW[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB)*S3_Snr_s("B",envS[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB),S3_Tnr("B",envW[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB)*S3_Snr_w("B",envW[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB)*S3_Sr_s("B",envS[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB),
                            S3_Sr_w("B",envW[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB)*S3_Sr_s("B",envS[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB)*S3_B("B",envS[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB), (1-S3_Tr("B",envW[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB))*S3_Sr_w("B",envW[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB)*S3_Snr_s("B",envS[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB),S3_Tr("B",envW[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB)*S3_Sr_w("B",envW[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB)*S3_Sr_s("B",envS[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB)),3,3,byrow = F)
        
        matsS1[,,y-1000]=(ann.matS1A+ann.matS1B)/2
        matsS2[,,y-1000]=(ann.matS2A+ann.matS2B)/2
        matsS3[,,y-1000]=(ann.matS3A+ann.matS3B)/2
        
        }
      
      # Update densities
      densA=densS[1:3]
      densB=densS[4:6]
      densMet=c(densA,densB)
      
      S2_densA=S2_densS[1:3]
      S2_densB=S2_densS[4:6]
      S2_densMet=c(S2_densA,S2_densB)
      
      S3_densA=S3_densS[1:3]
      S3_densB=S3_densS[4:6]
      S3_densMet=c(S3_densA,S3_densB)
      
      ### Update bi densities
      
      dens.biA=sum(densA)
      dens.biB=sum(densB)
      
      S2_dens.biA=sum(S2_densA)
      S2_dens.biB=sum(S2_densB)
      
      S3_dens.biA=sum(S3_densA)
      S3_dens.biB=sum(S3_densB)

     
    }
    
    # get summary metrics
    
    # environment 
    
    envir.sim=rbind(envir.sim,data.frame(envW=env[1000:years,1],
                                         envS=env[1000:years,2],
                                         sim=x,
                                         scen=k))
    temp=temp[!temp$year%in%1,]
    temp$year=temp$year-1000
    
    
    # total density per year and species
    densS1=aggregate(dens~year,data = temp[temp$species=="Prey"&temp$season%in%"winter",],sum)
    densS2=aggregate(dens~year,data = temp[temp$species=="Predator1"&temp$season%in%"winter",],sum)
    densS3=aggregate(dens~year,data = temp[temp$species=="Predator2"&temp$season%in%"winter",],sum)
    
    if(length(unique(densS1$dens))<2|length(unique(densS2$dens))<2|length(unique(densS3$dens))<2){
      
      print("Issue with Abundance")
      
      remove.sim=rbind(remove.sim,data.frame(sim=x,
                                             scen=k))
      next
      
      }else{
      
      # Extniction probability
      ext1=ext2=ext3=NA
      
      if(any(densS1$dens<1)) ext1 <- 1 else ext1 <- 0
      
      if(any(densS2$dens<1)) ext2 <- 1 else ext2 <- 0
      
      if(any(densS3$dens<1)) ext3 <- 1 else ext3 <- 0
      
      # generation time and net reproductive rate
      
      GTS1=mean(apply(matsS1[,,densS1$year[densS1$year%in%densS1$year[densS1$dens>0]&!densS1$year%in%(101)]],3,generation.time))
      RRS1=mean(apply(matsS1[,,densS1$year[densS1$year%in%densS1$year[densS1$dens>0]&!densS1$year%in%(101)]],3,net.reproductive.rate))
      
      GTS2=mean(apply(matsS2[,,densS2$year[densS2$year%in%densS2$year[densS2$dens>0]&!densS2$year%in%(101)]],3,generation.time))
      RRS2=mean(apply(matsS2[,,densS2$year[densS2$year%in%densS2$year[densS2$dens>0]&!densS2$year%in%(101)]],3,net.reproductive.rate))
      
      GTS3=mean(apply(matsS3[,,densS3$year[densS3$year%in%densS3$year[densS3$dens>0]&!densS3$year%in%(101)]],3,generation.time))
      RRS3=mean(apply(matsS3[,,densS3$year[densS3$year%in%densS3$year[densS3$dens>0]&!densS3$year%in%(101)]],3,net.reproductive.rate))
      
      # Mean total abundance (before extinction)
      
      mu.densS1=mean(densS1$dens[densS1$dens>0])
      mu.densS2=mean(densS2$dens[densS2$dens>0])
      mu.densS3=mean(densS3$dens[densS3$dens>0])
      
      # CV in total abundance (before extinction)
      
      cv.densS1=CV(densS1$dens[densS1$dens>0])
      cv.densS2=CV(densS2$dens[densS2$dens>0])
      cv.densS3=CV(densS3$dens[densS3$dens>0])
      
      # Covariation in total abudance
      if(ext1==0 & ext2==0 &ext3==0) {
        
        cor.ab[k,x,,]= cor(cbind(densS1$dens,densS2$dens,densS3$dens))
        
      
      }
      # mean stage-specific abundance
      
      densS1_stage=aggregate(dens~stage+year,data = temp[temp$species=="Prey"&temp$year%in%densS1$year[densS1$dens>0]&temp$season%in%"winter",],sum)
      mu.densS1_stage=aggregate(dens~stage,data = densS1_stage,mean)
      # CV
      cv.densS1_stage=aggregate(dens~stage,data = densS1_stage,CV)
      
      densS2_stage=aggregate(dens~stage+year,data = temp[temp$species=="Predator1"&temp$year%in%densS2$year[densS2$dens>0&temp$season%in%"winter"],],sum)
      mu.densS2_stage=aggregate(dens~stage,data = densS2_stage,mean)
      # CV
      cv.densS2_stage=aggregate(dens~stage,data = densS2_stage,CV)
      
      densS3_stage=aggregate(dens~stage+year,data = temp[temp$species=="Predator2"&temp$year%in%densS3$year[densS3$dens>0&temp$season%in%"winter"],],sum)
      mu.densS3_stage=aggregate(dens~stage,data = densS3_stage,mean)
      # CV
      cv.densS3_stage=aggregate(dens~stage,data = densS3_stage,CV)
      
      # seasonality in total abundance
      
      densS1v2=aggregate(dens~season+year,data = temp[temp$species=="Prey",],sum)
      
      densS1v2=densS1v2[densS1v2$year%in%densS1$year[densS1$dens>0]&!densS1v2$year%in%(years+1),]
      
      densS2v2=aggregate(dens~season+year,data = temp[temp$species=="Predator1",],sum)
      
      densS2v2=densS2v2[densS2v2$year%in%densS2$year[densS2$dens>0]&!densS2v2$year%in%(years+1),]
      
      densS3v2=aggregate(dens~season+year,data = temp[temp$species=="Predator2",],sum)
      
      densS3v2=densS3v2[densS3v2$year%in%densS3$year[densS3$dens>0]&!densS3v2$year%in%(years+1),]
      
      
      ts.S1=ts(densS1v2$dens,frequency = 2,start = c(1,1))
      
      seas.S1=decompose(ts.S1)$figure[2]
      
      ts.S2=ts(densS2v2$dens,frequency = 2,start = c(1,1))
      
      seas.S2=decompose(ts.S2)$figure[2]
      
      ts.S3=ts(densS3v2$dens,frequency = 2,start = c(1,1))
      
      seas.S3=decompose(ts.S3)$figure[2]
      
      
      # Seasonality in stage-specific abundance
      
      densS1v2_R=temp[temp$species=="Prey"&temp$stage%in%"R"&temp$year%in%densS1$year[densS1$dens>0]&!temp$year%in%(years+1),]
      densS1v2_N=temp[temp$species=="Prey"&temp$stage%in%"N"&temp$year%in%densS1$year[densS1$dens>0]&!temp$year%in%(years+1),]
      
      ts.S1=ts(densS1v2_R$dens,frequency = 2,start = c(1,1))
      
      seas.S1_R=decompose(ts.S1)$figure[2]
      
      ts.S1=ts(densS1v2_N$dens,frequency = 2,start = c(1,1))
      
      seas.S1_N=decompose(ts.S1)$figure[2]
      
      densS2v2_R=temp[temp$species=="Predator1"&temp$stage%in%"R"&temp$year%in%densS2$year[densS2$dens>0]&!temp$year%in%(years+1),]
      densS2v2_N=temp[temp$species=="Predator1"&temp$stage%in%"N"&temp$year%in%densS2$year[densS2$dens>0]&!temp$year%in%(years+1),]
      
      ts.S2=ts(densS2v2_R$dens,frequency = 2,start = c(1,1))
      
      seas.S2_R=decompose(ts.S2)$figure[2]
      
      ts.S2=ts(densS2v2_N$dens,frequency = 2,start = c(1,1))
      
      seas.S2_N=decompose(ts.S2)$figure[2]
      
      densS3v2_R=temp[temp$species=="Predator2"&temp$stage%in%"R"&temp$year%in%densS3$year[densS3$dens>0]&!temp$year%in%(years+1),]
      densS3v2_N=temp[temp$species=="Predator2"&temp$stage%in%"N"&temp$year%in%densS3$year[densS3$dens>0]&!temp$year%in%(years+1),]
      
      ts.S3=ts(densS3v2_R$dens,frequency = 2,start = c(1,1))
      
      seas.S3_R=decompose(ts.S3)$figure[2]
      
      ts.S3=ts(densS3v2_N$dens,frequency = 2,start = c(1,1))
      
      seas.S3_N=decompose(ts.S3)$figure[2]
      
      base.sim=rbind(base.sim,data.frame(species=c("Prey","Predator1","Predator2"),
                                         ext=c(ext1,ext2,ext3),
                                         Dens.mu=c(mu.densS1,mu.densS2,mu.densS3),
                                         Dens.cv=c(cv.densS1,cv.densS2,cv.densS3),
                                         seas=c(seas.S1,seas.S2,seas.S3),
                                         GT=c(GTS1,GTS2,GTS3),
                                         RR=c(RRS1,RRS2,RRS3),
                                         sim=x,
                                         scen=k))
      
      base.sim.stage=rbind(base.sim.stage,data.frame(species=rep(c("Prey","Predator1","Predator2"),each=3),
                                                     Dens.mu=c(mu.densS1_stage,mu.densS2_stage,mu.densS3_stage),
                                                     Dens.cv=c(cv.densS1_stage,cv.densS2_stage,cv.densS3_stage),
                                                     seas=c(0,seas.S1_N,seas.S1_R,0,seas.S2_N,seas.S2_R,0,seas.S3_N,seas.S3_R),
                                                     sim=x,
                                                     scen=k))
      
      
    } 
    
    
  }
}

base.sim[base.sim$ext%in%1,] # check for extinction

hist(base.sim$Dens.mu[base.sim$species%in%"Prey"])
hist(base.sim$Dens.mu[base.sim$species%in%"Predator1"])
hist(base.sim$Dens.mu[base.sim$species%in%"Predator2"])

# save output
write.csv(envir.sim,"envir.sim.csv",row.names = F)
write.csv(base.sim,"base.sim.csv",row.names = F)
write.csv(base.sim.stage,"base.sim.stage.csv",row.names = F)
write.csv(init.dens,"init.dens.csv",row.names = F)
write.csv(remove.sim,"remove.sim.csv",row.names = F)

save(cor.ab,file = "corr.sim.Rdata" )

#############################################################

## Perturb climate - covariance and higher variance

init.dens=read.csv("init.dens.csv")
remove.sim=read.csv("remove.sim.csv")

years=50 # years of simulations

# simulate winter and summer environment from multivariate normal

cov=c(0,0.2)
var=c(0,0.1)
cov.change=c("no","yes")

means=rbind(cbind(-1,1),cbind(-1.5,1),cbind(-1,0.5),cbind(-1.5,0.5),cbind(-0.5,0.5),cbind(-1.5,1.5))

clim.change.scen=c("baseline","-W","-S","-W-S","+W-S","-W+S")

var.change=c("no","yes")

clim.change.sim=NULL
clim.change.sim.stage=NULL
envir.clim.change.sim=NULL
remove.sim.clim.change=NULL

for(cp in 1:2){
  for(cpx in 1:2){
    
    covariance= matrix(c(0.5+var[cp],cov[cpx],cov[cpx],0.1+var[cp]),2,2)
    
    for(cc in 1:length(clim.change.scen)){
      
      for(k in 1:48){
        
        print(paste("Running scenario ", k, "for climate change ", clim.change.scen[cc],"cov change", cov.change[cpx],"and var change", var.change[cp] ))
        
        # Parameters  prey
        A.mean.vec <- c(scen.Sr_w[k],scen.D[k],scen.Sd[k],scen.B[k], scen.Snr_w[k],
                        scen.Sj_w[k],scen.Tj[k],scen.Tnr[k],scen.Tr[k],scen.Sr_s[k],scen.Snr_s[k])
        
        AB.diff.vec <- c(AB.scen.Sr_w[k]*abs(A.mean.vec[1]),0,AB.scen.Sd[k]*abs(A.mean.vec[3]),
                         0,AB.scen.Snr_w[k]*abs(A.mean.vec[5]),
                         AB.scen.Sj_w[k]*abs(A.mean.vec[6]),0,
                         0,0,
                         AB.scen.Sr_s[k]*abs(A.mean.vec[10]),AB.scen.Snr_s[k]*abs(A.mean.vec[11]))
        
        slope.env.vec <- c(env.scen.Sr_w[k]*abs(A.mean.vec[1]),env.scen.D[k]*abs(A.mean.vec[2]),env.scen.Sd[k]*abs(A.mean.vec[3]),
                           env.scen.B[k]*abs(A.mean.vec[4]),env.scen.Snr_w[k]*abs(A.mean.vec[5]),
                           env.scen.Sj_w[k]*abs(A.mean.vec[6]),env.scen.Tj[k]*abs(A.mean.vec[7]),
                           env.scen.Tnr[k]*abs(A.mean.vec[8]),env.scen.Tr[k]*abs(A.mean.vec[9]),
                           env.scen.Sr_s[k]*abs(A.mean.vec[10]),env.scen.Snr_s[k]*abs(A.mean.vec[11]))
        slope.dens.vec <- c(-0.02*abs(A.mean.vec[1]),0.02*abs(A.mean.vec[2]),-0.2*abs(A.mean.vec[3]),-0.03*abs(A.mean.vec[4]),
                            -0.02*abs(A.mean.vec[5]),0.1*abs(A.mean.vec[6]),-0.02*abs(A.mean.vec[7]),-0.1*abs(A.mean.vec[8]),-0.02*abs(A.mean.vec[9]),-0.01*abs(A.mean.vec[10]),-0.01*abs(A.mean.vec[11]))
        slope.bi.vec <- c(bi.scen.Sr_w[k]*abs(slope.dens.vec)[1],bi.scen.D[k]*abs(slope.dens.vec)[2],bi.scen.Sd[k]*abs(slope.dens.vec)[3],bi.scen.B[k]*abs(slope.dens.vec)[4],
                          bi.scen.Snr_w[k]*abs(slope.dens.vec)[5],bi.scen.Sj_w[k]*abs(slope.dens.vec)[6],bi.scen.Tj[k]*abs(slope.dens.vec)[7],bi.scen.Tnr[k]*abs(slope.dens.vec)[8],
                          bi.scen.Tr[k]*abs(slope.dens.vec)[9],bi.scen.Sr_s[k]*abs(slope.dens.vec)[10],bi.scen.Snr_s[k]*abs(slope.dens.vec)[11])
        
        slope.biS2.vec <- c(biS2.scen.Sr_w[k]*abs(slope.dens.vec)[1],biS2.scen.D[k]*abs(slope.dens.vec)[2],biS2.scen.Sd[k]*abs(slope.dens.vec)[3],biS2.scen.B[k]*abs(slope.dens.vec)[4],
                            biS2.scen.Snr_w[k]*abs(slope.dens.vec)[5],biS2.scen.Sj_w[k]*abs(slope.dens.vec)[6],biS2.scen.Tj[k]*abs(slope.dens.vec)[7],biS2.scen.Tnr[k]*abs(slope.dens.vec)[8],
                            biS2.scen.Tr[k]*abs(slope.dens.vec)[9],biS2.scen.Sr_s[k]*abs(slope.dens.vec)[10],biS2.scen.Snr_s[k]*abs(slope.dens.vec)[11])
        
        AB.env.diff.vec<- c(AB.env.scen.Sr_w[k]*abs(slope.env.vec[1]),AB.env.scen.D[k]*abs(slope.env.vec[2]),AB.env.scen.Sd[k]*abs(slope.env.vec[3]),
                            AB.env.scen.B[k]*abs(slope.env.vec[4]),AB.env.scen.Snr_w[k]*abs(slope.env.vec[5]),
                            AB.env.scen.Sj_w[k]*abs(slope.env.vec[6]),AB.env.scen.Tj[k]*abs(slope.env.vec[7]),
                            AB.env.scen.Tnr[k]*abs(slope.env.vec[8]),AB.env.scen.Tr[k]*abs(slope.env.vec[9]),
                            AB.env.scen.Sr_s[k]*abs(slope.env.vec[10]),AB.env.scen.Snr_s[k]*abs(slope.env.vec[11]))
        
        env.dens.diff.vec<- c(0.01,0.03,0.01,0.01,0.01,0.01,0,0,0,0.01,0.01) # density has + effect at bad env and - effect at good env
        dens.dens.bi.diff.vec<- c(dens.bi.scen.Sr_w[k]*abs(slope.dens.vec)[1],dens.bi.scen.D[k]*abs(slope.dens.vec)[2],dens.bi.scen.Sd[k]*abs(slope.dens.vec)[3],dens.bi.scen.B[k]*abs(slope.dens.vec)[4],
                                  dens.bi.scen.Snr_w[k]*abs(slope.dens.vec)[5],dens.bi.scen.Sj_w[k]*abs(slope.dens.vec)[6],dens.bi.scen.Tj[k]*abs(slope.dens.vec)[7],dens.bi.scen.Tnr[k]*abs(slope.dens.vec)[8],
                                  dens.bi.scen.Tr[k]*abs(slope.dens.vec)[9],dens.bi.scen.Sr_s[k]*abs(slope.dens.vec)[10],dens.bi.scen.Snr_s[k]*abs(slope.dens.vec)[11])
        AB.dens.bi.diff.vec<- rep(0,length(vr))
        
        # Parameters predator
        
        A.mean.vecS2 <- c(S2.scen.Sr_w[k],S2.scen.D[k],S2.scen.Sd[k],S2.scen.B[k], S2.scen.Snr_w[k],
                          S2.scen.Sj_w[k],S2.scen.Tj[k],S2.scen.Tnr[k],S2.scen.Tr[k],S2.scen.Sr_s[k],S2.scen.Snr_s[k])
        AB.diff.vecS2 <- c(S2.AB.scen.Sr_w[k]*abs(A.mean.vecS2[1]),0,S2.AB.scen.Sd[k]*abs(A.mean.vecS2[3]),
                           0,S2.AB.scen.Snr_w[k]*abs(A.mean.vecS2[5]),
                           S2.AB.scen.Sj_w[k]*abs(A.mean.vecS2[6]),0,
                           0,0,
                           S2.AB.scen.Sr_s[k]*abs(A.mean.vecS2[10]),S2.AB.scen.Snr_s[k]*abs(A.mean.vecS2[11]))
        
        
        slope.env.vecS2 <- c(S2.env.scen.Sr_w[k]*abs(A.mean.vecS2[1]),S2.env.scen.D[k]*abs(A.mean.vecS2[2]),S2.env.scen.Sd[k]*abs(A.mean.vecS2[3]),
                             S2.env.scen.B[k]*abs(A.mean.vecS2[4]),S2.env.scen.Snr_w[k]*abs(A.mean.vecS2[5]),
                             S2.env.scen.Sj_w[k]*abs(A.mean.vecS2[6]),S2.env.scen.Tj[k]*abs(A.mean.vecS2[7]),
                             S2.env.scen.Tnr[k]*abs(A.mean.vecS2[8]),S2.env.scen.Tr[k]*abs(A.mean.vecS2[9]),
                             S2.env.scen.Sr_s[k]*abs(A.mean.vecS2[10]),S2.env.scen.Snr_s[k]*abs(A.mean.vecS2[11]))
        slope.dens.vecS2 <- c(-0.05*abs(A.mean.vecS2[1]),5*abs(A.mean.vecS2[2]),-0.05*abs(A.mean.vecS2[3]),-0.05*abs(A.mean.vecS2[4]),
                              -0.05*abs(A.mean.vecS2[5]),-2*abs(A.mean.vecS2[6]),-0.01*abs(A.mean.vecS2[7]),-2*abs(A.mean.vecS2[8]),-0.01*abs(A.mean.vecS2[9]),-0.01*abs(A.mean.vecS2[10]),-0.01*abs(A.mean.vecS2[11]))
        
        slope.bi.vecS2 <-  c(S2.bi.scen.Sr_w[k]*abs(slope.dens.vecS2)[1],S2.bi.scen.D[k]*abs(slope.dens.vecS2)[2],S2.bi.scen.Sd[k]*abs(slope.dens.vecS2)[3],
                             S2.bi.scen.B[k]*abs(slope.dens.vecS2)[4],S2.bi.scen.Snr_w[k]*abs(slope.dens.vecS2)[5],S2.bi.scen.Sj_w[k]*abs(slope.dens.vecS2)[6],
                             S2.bi.scen.Tj[k]*abs(slope.dens.vecS2)[7],S2.bi.scen.Tnr[k]*abs(slope.dens.vecS2)[8],
                             S2.bi.scen.Tr[k]*abs(slope.dens.vecS2)[9],S2.bi.scen.Sr_s[k]*abs(slope.dens.vecS2)[10],S2.bi.scen.Snr_s[k]*abs(slope.dens.vecS2)[11])
        
        slope.biS2.vecS2 <-  c(S2.biS2.scen.Sr_w[k]*abs(slope.dens.vecS2)[1],S2.biS2.scen.D[k]*abs(slope.dens.vecS2)[2],S2.biS2.scen.Sd[k]*abs(slope.dens.vecS2)[3],
                               S2.biS2.scen.B[k]*abs(slope.dens.vecS2)[4],S2.biS2.scen.Snr_w[k]*abs(slope.dens.vecS2)[5],S2.biS2.scen.Sj_w[k]*abs(slope.dens.vecS2)[6],
                               S2.biS2.scen.Tj[k]*abs(slope.dens.vecS2)[7],S2.biS2.scen.Tnr[k]*abs(slope.dens.vecS2)[8],
                               S2.biS2.scen.Tr[k]*abs(slope.dens.vecS2)[9],S2.biS2.scen.Sr_s[k]*abs(slope.dens.vecS2)[10],S2.biS2.scen.Snr_s[k]*abs(slope.dens.vecS2)[11])
        
        slope.bi2.vecS2 <- c(-0.001,0,-0.001,S2.bi2.scen.B[k],-0.001,-0.001,-0.001,-0.001,-0.001,-0.001,-0.001)
        AB.env.diff.vecS2<- c(S2.AB.env.scen.Sr_w[k]*abs(slope.env.vecS2[1]),S2.AB.env.scen.D[k]*abs(slope.env.vecS2[2]),S2.AB.env.scen.Sd[k]*abs(slope.env.vecS2[3]),
                              S2.AB.env.scen.B[k]*abs(slope.env.vecS2[4]),S2.AB.env.scen.Snr_w[k]*abs(slope.env.vecS2[5]),
                              S2.AB.env.scen.Sj_w[k]*abs(slope.env.vecS2[6]),S2.AB.env.scen.Tj[k]*abs(slope.env.vecS2[7]),
                              S2.AB.env.scen.Tnr[k]*abs(slope.env.vecS2[8]),S2.AB.env.scen.Tr[k]*abs(slope.env.vecS2[9]),
                              S2.AB.env.scen.Sr_s[k]*abs(slope.env.vecS2[10]),S2.AB.env.scen.Snr_s[k]*abs(slope.env.vecS2[11]))
        
        dens.dens.bi.diff.vecS2<- c(S2.dens.bi.scen.Sr_w[k]*abs(slope.dens.vecS2)[1],S2.dens.bi.scen.D[k]*abs(slope.dens.vecS2)[2],S2.dens.bi.scen.Sd[k]*abs(slope.dens.vecS2)[3],
                                    S2.dens.bi.scen.B[k]*abs(slope.dens.vecS2)[4],S2.dens.bi.scen.Snr_w[k]*abs(slope.dens.vecS2)[5],S2.dens.bi.scen.Sj_w[k]*abs(slope.dens.vecS2)[6],
                                    S2.dens.bi.scen.Tj[k]*abs(slope.dens.vecS2)[7],S2.dens.bi.scen.Tnr[k]*abs(slope.dens.vecS2)[8],
                                    S2.dens.bi.scen.Tr[k]*abs(slope.dens.vecS2)[9],S2.dens.bi.scen.Sr_s[k]*abs(slope.dens.vecS2)[10],S2.dens.bi.scen.Snr_s[k]*abs(slope.dens.vecS2)[11])
        env.dens.diff.vecS2<- c(0.01,0.01,0,0.01,0.01,0.01,0,0,0,0.01,0.01) # density has + effect at bad env and - effect at good env
        
        # Parameters predator 2
        
        A.mean.vecS3 <- c(S3.scen.Sr_w[k],S3.scen.D[k],S3.scen.Sd[k],S3.scen.B[k], S3.scen.Snr_w[k],
                          S3.scen.Sj_w[k],S3.scen.Tj[k],S3.scen.Tnr[k],S3.scen.Tr[k],S3.scen.Sr_s[k],S3.scen.Snr_s[k])
        AB.diff.vecS3 <- c(S3.AB.scen.Sr_w[k]*abs(A.mean.vecS3[1]),0,S3.AB.scen.Sd[k]*abs(A.mean.vecS3[3]),
                           0,S3.AB.scen.Snr_w[k]*abs(A.mean.vecS3[5]),
                           S3.AB.scen.Sj_w[k]*abs(A.mean.vecS3[6]),0,
                           0,0,
                           S3.AB.scen.Sr_s[k]*abs(A.mean.vecS3[10]),S3.AB.scen.Snr_s[k]*abs(A.mean.vecS3[11]))
        
        
        slope.env.vecS3 <- c(S3.env.scen.Sr_w[k]*abs(A.mean.vecS3[1]),S3.env.scen.D[k]*abs(A.mean.vecS3[2]),S3.env.scen.Sd[k]*abs(A.mean.vecS3[3]),
                             S3.env.scen.B[k]*abs(A.mean.vecS3[4]),S3.env.scen.Snr_w[k]*abs(A.mean.vecS3[5]),
                             S3.env.scen.Sj_w[k]*abs(A.mean.vecS3[6]),S3.env.scen.Tj[k]*abs(A.mean.vecS3[7]),
                             S3.env.scen.Tnr[k]*abs(A.mean.vecS3[8]),S3.env.scen.Tr[k]*abs(A.mean.vecS3[9]),
                             S3.env.scen.Sr_s[k]*abs(A.mean.vecS3[10]),S3.env.scen.Snr_s[k]*abs(A.mean.vecS3[11]))
        slope.dens.vecS3 <- c(-0.05*abs(A.mean.vecS3[1]),5*abs(A.mean.vecS3[2]),-0.05*abs(A.mean.vecS3[3]),-0.05*abs(A.mean.vecS3[4]),
                              -0.05*abs(A.mean.vecS3[5]),-2*abs(A.mean.vecS3[6]),-0.01*abs(A.mean.vecS3[7]),-2*abs(A.mean.vecS3[8]),-0.01*abs(A.mean.vecS3[9]),-0.01*abs(A.mean.vecS3[10]),-0.01*abs(A.mean.vecS3[11]))
        
        slope.bi.vecS3 <-  c(S3.bi.scen.Sr_w[k]*abs(slope.dens.vecS3)[1],S3.bi.scen.D[k]*abs(slope.dens.vecS3)[2],S3.bi.scen.Sd[k]*abs(slope.dens.vecS3)[3],
                             S3.bi.scen.B[k]*abs(slope.dens.vecS3)[4],S3.bi.scen.Snr_w[k]*abs(slope.dens.vecS3)[5],S3.bi.scen.Sj_w[k]*abs(slope.dens.vecS3)[6],
                             S3.bi.scen.Tj[k]*abs(slope.dens.vecS3)[7],S3.bi.scen.Tnr[k]*abs(slope.dens.vecS3)[8],
                             S3.bi.scen.Tr[k]*abs(slope.dens.vecS3)[9],S3.bi.scen.Sr_s[k]*abs(slope.dens.vecS3)[10],S3.bi.scen.Snr_s[k]*abs(slope.dens.vecS3)[11])
        
        slope.biS2.vecS3 <-  c(S3.biS2.scen.Sr_w[k]*abs(slope.dens.vecS2)[1],S3.biS2.scen.D[k]*abs(slope.dens.vecS2)[2],S3.biS2.scen.Sd[k]*abs(slope.dens.vecS2)[3],
                               S3.biS2.scen.B[k]*abs(slope.dens.vecS2)[4],S3.biS2.scen.Snr_w[k]*abs(slope.dens.vecS2)[5],S3.biS2.scen.Sj_w[k]*abs(slope.dens.vecS2)[6],
                               S3.biS2.scen.Tj[k]*abs(slope.dens.vecS2)[7],S3.biS2.scen.Tnr[k]*abs(slope.dens.vecS2)[8],
                               S3.biS2.scen.Tr[k]*abs(slope.dens.vecS2)[9],S3.biS2.scen.Sr_s[k]*abs(slope.dens.vecS2)[10],S3.biS2.scen.Snr_s[k]*abs(slope.dens.vecS2)[11])
        
        slope.bi2.vecS3 <- c(-0.001,0,-0.001,S3.bi2.scen.B[k],-0.001,-0.001,-0.001,-0.001,-0.001,-0.001,-0.001)
        
        AB.env.diff.vecS3<- c(S3.AB.env.scen.Sr_w[k]*abs(slope.env.vecS3[1]),S3.AB.env.scen.D[k]*abs(slope.env.vecS3[2]),S3.AB.env.scen.Sd[k]*abs(slope.env.vecS3[3]),
                              S3.AB.env.scen.B[k]*abs(slope.env.vecS3[4]),S3.AB.env.scen.Snr_w[k]*abs(slope.env.vecS3[5]),
                              S3.AB.env.scen.Sj_w[k]*abs(slope.env.vecS3[6]),S3.AB.env.scen.Tj[k]*abs(slope.env.vecS3[7]),
                              S3.AB.env.scen.Tnr[k]*abs(slope.env.vecS3[8]),S3.AB.env.scen.Tr[k]*abs(slope.env.vecS3[9]),
                              S3.AB.env.scen.Sr_s[k]*abs(slope.env.vecS3[10]),S3.AB.env.scen.Snr_s[k]*abs(slope.env.vecS3[11]))
        
        dens.dens.bi.diff.vecS3<- c(S3.dens.bi.scen.Sr_w[k],S3.dens.bi.scen.D[k],S3.dens.bi.scen.Sd[k],
                                    S3.dens.bi.scen.B[k],S3.dens.bi.scen.Snr_w[k],
                                    S3.dens.bi.scen.Sj_w[k],S3.dens.bi.scen.Tj[k],
                                    S3.dens.bi.scen.Tnr[k],S3.dens.bi.scen.Tr[k],
                                    S3.dens.bi.scen.Sr_s[k],S3.dens.bi.scen.Snr_s[k])
        
        
        env.dens.diff.vecS3<- c(0.01,0.01,0,0.01,0.01,0.01,0,0,0,0.01,0.01) # density has + effect at bad env and - effect at good env
        
        remove=remove.sim$sim[remove.sim$scen%in%k]
        simul=unique(init.dens$sim[!init.dens$sim%in%remove])
        
        for(x in simul){
          
          
          env=mvrnorm(years, mu=means[cc,], covariance)
          
          envW=env[,1]
          envS=env[,2]
          
          # Initial prey 
          densA=init.dens$densA[init.dens$sim%in%x&init.dens$scen%in%k]
          densB=init.dens$densB[init.dens$sim%in%x&init.dens$scen%in%k]
          
          densMet=c(densA,densB)
          
          dens.biA=sum(densA)
          dens.biB=sum(densB)
          
          # Initial predator 1
          
          S2_densA=init.dens$S2_densA[init.dens$sim%in%x&init.dens$scen%in%k]
          S2_densB=init.dens$S2_densB[init.dens$sim%in%x&init.dens$scen%in%k]
          
          S2_densMet=c(S2_densA,S2_densB)
          
          
          S2_dens.biA=sum(S2_densA)
          S2_dens.biB=sum(S2_densB)
          
          # Initial predator 2
          
          S3_densA=init.dens$S3_densA[init.dens$sim%in%x&init.dens$scen%in%k]
          S3_densB=init.dens$S3_densB[init.dens$sim%in%x&init.dens$scen%in%k]
          
          S3_densMet=c(S3_densA,S3_densB)
          
          
          S3_dens.biA=sum(S3_densA)
          S3_dens.biB=sum(S3_densB)
          
          
          ## data frame to hold results in
          
          temp=data.frame(dens=c(densA,densB,S2_densA,S2_densB,S3_densA,S3_densB),
                          species=rep(c("Prey","Predator1", "Predator2"),each=6),
                          site=rep(rep(c("A","B"),each=3),3),
                          season="summer",
                          stage=rep(c("J","N","R"),6),
                          year=1)
          
          matsS1=array(0,c(3,3,50))
          matsS2=array(0,c(3,3,50))
          matsS3=array(0,c(3,3,50))
          
          for(y in 1:years){
            
            #summer
            
            #Prey
            summerA=matrix(c(0,0,0,
                             0,Snr_s("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA),0,
                             Sr_s("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*B("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA),0,Sr_s("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)),3,3,byrow = F)
            
            summerB=matrix(c(0,0,0,
                             0,Snr_s("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB),0,
                             Sr_s("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*B("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB),0,Sr_s("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)),3,3,byrow = F)
            
            
            
            Bs[1:3,1:3]=summerA
            Bs[4:6,4:6]=summerB
            
            
            Ms[3,3]=1-D("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)
            Ms[4,4]=1-D("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)
            Ms[6,3]=D("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Ds("A",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)
            Ms[5,4]=D("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Ds("B",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)
            
            
            
            densW=t(Ps)%*%Ms%*%Ps%*%Bs%*%densMet
            
            
            temp=rbind(temp,data.frame(dens=c(round(densW[1:3],2),round(densW[4:6],2)),
                                       species="Prey",
                                       site=rep(c("A","B"),each=3),
                                       season=rep("winter",6),stage=rep(c("J","N","R"),2),year=rep(y,6)))
            
            
            
            densA=densW[1:3]
            densB=densW[4:6]
            densMet=c(densA,densB)
            
            
            #predator1 
            S2_summerA=matrix(c(0,0,0,
                                0,S2_Snr_s("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA),0,
                                S2_Sr_s("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_B("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA),0,S2_Sr_s("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)),3,3,byrow = F)
            
            S2_summerB=matrix(c(0,0,0,
                                0,S2_Snr_s("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB),0,
                                S2_Sr_s("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)*S2_B("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB),0,S2_Sr_s("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)),3,3,byrow = F)
            
            
            
            S2_Bs[1:3,1:3]=S2_summerA
            S2_Bs[4:6,4:6]=S2_summerB
            
            
            S2_Ms[3,3]=1-S2_D("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)
            S2_Ms[4,4]=1-S2_D("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)
            S2_Ms[6,3]=S2_D("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Ds("A",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)
            S2_Ms[5,4]=S2_D("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)*S2_Ds("B",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)
            
            S2_densW=t(S2_Ps)%*%S2_Ms%*%S2_Ps%*%S2_Bs%*%S2_densMet
            
            
            temp=rbind(temp,data.frame(dens=c(round(S2_densW[1:3],2),round(S2_densW[4:6],2)),
                                       species="Predator1",
                                       site=rep(c("A","B"),each=3),
                                       season=rep("winter",6),stage=rep(c("J","N","R"),2),year=rep(y,6)))
            
            
            
            S2_densA=S2_densW[1:3]
            S2_densB=S2_densW[4:6]
            S2_densMet=c(S2_densA,S2_densB)
            
            #predator 2 
            S3_summerA=matrix(c(0,0,0,
                                0,S3_Snr_s("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA),0,
                                S3_Sr_s("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_B("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA),0,S3_Sr_s("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)),3,3,byrow = F)
            
            S3_summerB=matrix(c(0,0,0,
                                0,S3_Snr_s("B",envS[y],sum(S3_densB[1:3]),dens.biB,S2_dens.biB),0,
                                S3_Sr_s("B",envS[y],sum(S3_densB[1:3]),dens.biB,S2_dens.biB)*S3_B("B",envS[y],sum(S3_densB[1:3]),dens.biB,S2_dens.biB),0,S3_Sr_s("B",envS[y],sum(S3_densB[1:3]),dens.biB,S2_dens.biB)),3,3,byrow = F)
            
            
            
            S3_Bs[1:3,1:3]=S3_summerA
            S3_Bs[4:6,4:6]=S3_summerB
            
            
            S3_Ms[3,3]=1-S3_D("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)
            S3_Ms[4,4]=1-S3_D("B",envS[y],sum(S3_densB[1:3]),dens.biB,S2_dens.biB)
            S3_Ms[6,3]=S3_D("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Ds("A",envS[y],sum(S3_densB[1:3]),dens.biB,S2_dens.biB)
            S3_Ms[5,4]=S3_D("B",envS[y],sum(S3_densB[1:3]),dens.biB,S2_dens.biB)*S3_Ds("B",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)
            
            S3_densW=t(S3_Ps)%*%S3_Ms%*%S3_Ps%*%S3_Bs%*%S3_densMet
            
            
            temp=rbind(temp,data.frame(dens=c(round(S3_densW[1:3],2),round(S3_densW[4:6],2)),
                                       species="Predator2",
                                       site=rep(c("A","B"),each=3),
                                       season=rep("winter",6),stage=rep(c("J","N","R"),2),year=rep(y,6)))
            
            
            S3_densA=S3_densW[1:3]
            S3_densB=S3_densW[4:6]
            S3_densMet=c(S3_densA,S3_densB)
            
            ### Update bi densities
            
           
            # 
            dens.biA=sum(densA)
            dens.biB=sum(densB)
            
            S2_dens.biA=sum(S2_densA)
            S2_dens.biB=sum(S2_densB)
            
            S3_dens.biA=sum(S3_densA)
            S3_dens.biB=sum(S3_densB)
            
            # WINTER
            
            # Prey
            
            winterA=matrix(c(0,(1-Tj("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA))*Sj_w("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA),Tj("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Sj_w("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA),
                             0,(1-Tnr("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA))*Snr_w("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA),Tnr("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Snr_w("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA),
                             0, (1-Tr("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA))*Sr_w("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA),Tr("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Sr_w("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)),3,3,byrow = F)
            
            winterB=matrix(c(0,(1-Tj("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB))*Sj_w("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB),Tj("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Sj_w("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB),
                             0,(1-Tnr("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB))*Snr_w("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB),Tnr("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Snr_w("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB),
                             0, (1-Tr("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB))*Sr_w("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB),Tr("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Sr_w("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)),3,3,byrow = F)
            
            
            Bw[1:3,1:3]=winterA
            Bw[4:6,4:6]=winterB
            
            densS=t(Pw)%*%Mw%*%Pw%*%Bw%*%densMet
            
            
            temp=rbind(temp,data.frame(dens=c(round(densS[1:3],2),round(densS[4:6],2)),
                                       species="Prey",
                                       site=rep(c("A","B"),each=3),
                                       season=rep("summer",6),stage=rep(c("J","N","R"),2),year=rep(y+1,6)))
            
            
            # Predator 1
            S2_winterA=matrix(c(0,(1-S2_Tj("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA))*S2_Sj_w("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA),S2_Tj("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Sj_w("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA),
                                0,(1-S2_Tnr("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA))*S2_Snr_w("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA),S2_Tnr("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Snr_w("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA),
                                0, (1-S2_Tr("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA))*S2_Sr_w("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA),S2_Tr("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Sr_w("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)),3,3,byrow = F)
            
            S2_winterB=matrix(c(0,(1-S2_Tj("B",envW[y],sum(S2_densA[1:3]),dens.biB,S3_dens.biB))*S2_Sj_w("B",envW[y],sum(S2_densA[1:3]),dens.biB,S3_dens.biB),S2_Tj("B",envW[y],sum(S2_densA[1:3]),dens.biB,S3_dens.biB)*S2_Sj_w("B",envW[y],sum(S2_densA[1:3]),dens.biB,S3_dens.biB),
                                0,(1-S2_Tnr("B",envW[y],sum(S2_densA[1:3]),dens.biB,S3_dens.biB))*S2_Snr_w("B",envW[y],sum(S2_densA[1:3]),dens.biB,S3_dens.biB),S2_Tnr("B",envW[y],sum(S2_densA[1:3]),dens.biB,S3_dens.biB)*S2_Snr_w("B",envW[y],sum(S2_densA[1:3]),dens.biB,S3_dens.biB),
                                0, (1-S2_Tr("B",envW[y],sum(S2_densA[1:3]),dens.biB,S3_dens.biB))*S2_Sr_w("B",envW[y],sum(S2_densA[1:3]),dens.biB,S3_dens.biB),S2_Tr("B",envW[y],sum(S2_densA[1:3]),dens.biB,S3_dens.biB)*S2_Sr_w("B",envW[y],sum(S2_densA[1:3]),dens.biB,S3_dens.biB)),3,3,byrow = F)
            
            
            S2_Bw[1:3,1:3]=S2_winterA
            S2_Bw[4:6,4:6]=S2_winterB
            
            S2_densS=t(S2_Pw)%*%S2_Mw%*%S2_Pw%*%S2_Bw%*%S2_densMet
            
            
            
            temp=rbind(temp,data.frame(dens=c(round(S2_densS[1:3],2),round(S2_densS[4:6],2)),
                                       species="Predator1",
                                       site=rep(c("A","B"),each=3),
                                       season=rep("summer",6),stage=rep(c("J","N","R"),2),year=rep(y+1,6)))
            
            # Predator 2
            S3_winterA=matrix(c(0,(1-S3_Tj("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA))*S3_Sj_w("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA),S3_Tj("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Sj_w("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA),
                                0,(1-S3_Tnr("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA))*S3_Snr_w("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA),S3_Tnr("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Snr_w("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA),
                                0, (1-S3_Tr("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA))*S3_Sr_w("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA),S3_Tr("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Sr_w("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)),3,3,byrow = F)
            
            S3_winterB=matrix(c(0,(1-S3_Tj("B",envW[y],sum(S3_densA[1:3]),dens.biB,S3_dens.biB))*S3_Sj_w("B",envW[y],sum(S3_densA[1:3]),dens.biB,S3_dens.biB),S3_Tj("B",envW[y],sum(S3_densA[1:3]),dens.biB,S3_dens.biB)*S3_Sj_w("B",envW[y],sum(S3_densA[1:3]),dens.biB,S3_dens.biB),
                                0,(1-S3_Tnr("B",envW[y],sum(S3_densA[1:3]),dens.biB,S3_dens.biB))*S3_Snr_w("B",envW[y],sum(S3_densA[1:3]),dens.biB,S3_dens.biB),S3_Tnr("B",envW[y],sum(S3_densA[1:3]),dens.biB,S3_dens.biB)*S3_Snr_w("B",envW[y],sum(S3_densA[1:3]),dens.biB,S3_dens.biB),
                                0, (1-S3_Tr("B",envW[y],sum(S3_densA[1:3]),dens.biB,S3_dens.biB))*S3_Sr_w("B",envW[y],sum(S3_densA[1:3]),dens.biB,S3_dens.biB),S3_Tr("B",envW[y],sum(S3_densA[1:3]),dens.biB,S3_dens.biB)*S3_Sr_w("B",envW[y],sum(S3_densA[1:3]),dens.biB,S3_dens.biB)),3,3,byrow = F)
            
            
            S3_Bw[1:3,1:3]=S3_winterA
            S3_Bw[4:6,4:6]=S3_winterB
            
            S3_densS=t(S3_Pw)%*%S3_Mw%*%S3_Pw%*%S3_Bw%*%S3_densMet
            
            temp=rbind(temp,data.frame(dens=c(round(S3_densS[1:3],2),round(S3_densS[4:6],2)),
                                       species="Predator2",
                                       site=rep(c("A","B"),each=3),
                                       season=rep("summer",6),stage=rep(c("J","N","R"),2),year=rep(y+1,6)))
            
            
            ## save MPMs to calculate life histories:
            
            ann.matS1A=matrix(c(0,(1-Tj("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA))*Sj_w("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Snr_s("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA),Tj("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Sj_w("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Sr_s("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA),
                                0,(1-Tnr("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA))*Snr_w("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Snr_s("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA),Tnr("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Snr_w("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Sr_s("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA),
                                Sr_w("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Sr_s("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*B("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA), (1-Tr("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA))*Sr_w("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Snr_s("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA),Tr("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Sr_w("A",envW[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Sr_s("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)),3,3,byrow = F)
            
            ann.matS1B=matrix(c(0,(1-Tj("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB))*Sj_w("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Snr_s("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB),Tj("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Sj_w("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Sr_s("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB),
                                0,(1-Tnr("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB))*Snr_w("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Snr_s("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB),Tnr("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Snr_w("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Sr_s("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB),
                                Sr_w("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Sr_s("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*B("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB), (1-Tr("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB))*Sr_w("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Snr_s("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB),Tr("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Sr_w("B",envW[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Sr_s("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)),3,3,byrow = F)
            
            ann.matS2A=matrix(c(0,(1-S2_Tj("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA))*S2_Sj_w("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Snr_s("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA),S2_Tj("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Sj_w("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Sr_s("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA),
                                0,(1-S2_Tnr("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA))*S2_Snr_w("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Snr_s("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA),S2_Tnr("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Snr_w("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Sr_s("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA),
                                S2_Sr_w("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Sr_s("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_B("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA), (1-S2_Tr("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA))*S2_Sr_w("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Snr_s("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA),S2_Tr("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Sr_w("A",envW[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)*S2_Sr_s("A",envS[y],sum(S2_densA[1:3]),dens.biA,S3_dens.biA)),3,3,byrow = F)
            
            ann.matS2B=matrix(c(0,(1-S2_Tj("B",envW[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB))*S2_Sj_w("B",envW[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)*S2_Snr_s("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB),S2_Tj("B",envW[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)*S2_Sj_w("B",envW[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)*S2_Sr_s("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB),
                                0,(1-S2_Tnr("B",envW[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB))*S2_Snr_w("B",envW[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)*S2_Snr_s("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB),S2_Tnr("B",envW[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)*S2_Snr_w("B",envW[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)*S2_Sr_s("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB),
                                S2_Sr_w("B",envW[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)*S2_Sr_s("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)*S2_B("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB), (1-S2_Tr("B",envW[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB))*S2_Sr_w("B",envW[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)*S2_Snr_s("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB),S2_Tr("B",envW[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)*S2_Sr_w("B",envW[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)*S2_Sr_s("B",envS[y],sum(S2_densB[1:3]),dens.biB,S3_dens.biB)),3,3,byrow = F)
            
            ann.matS3A=matrix(c(0,(1-S3_Tj("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA))*S3_Sj_w("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Snr_s("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA),S3_Tj("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Sj_w("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Sr_s("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA),
                                0,(1-S3_Tnr("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA))*S3_Snr_w("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Snr_s("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA),S3_Tnr("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Snr_w("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Sr_s("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA),
                                S3_Sr_w("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Sr_s("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_B("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA), (1-S3_Tr("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA))*S3_Sr_w("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Snr_s("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA),S3_Tr("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Sr_w("A",envW[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)*S3_Sr_s("A",envS[y],sum(S3_densA[1:3]),dens.biA,S2_dens.biA)),3,3,byrow = F)
            
            ann.matS3B=matrix(c(0,(1-S3_Tj("B",envW[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB))*S3_Sj_w("B",envW[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB)*S3_Snr_s("B",envS[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB),S3_Tj("B",envW[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB)*S3_Sj_w("B",envW[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB)*S3_Sr_s("B",envS[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB),
                                0,(1-S3_Tnr("B",envW[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB))*S3_Snr_w("B",envW[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB)*S3_Snr_s("B",envS[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB),S3_Tnr("B",envW[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB)*S3_Snr_w("B",envW[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB)*S3_Sr_s("B",envS[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB),
                                S3_Sr_w("B",envW[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB)*S3_Sr_s("B",envS[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB)*S3_B("B",envS[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB), (1-S3_Tr("B",envW[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB))*S3_Sr_w("B",envW[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB)*S3_Snr_s("B",envS[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB),S3_Tr("B",envW[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB)*S3_Sr_w("B",envW[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB)*S3_Sr_s("B",envS[y],sum(S3_densB[1:3]),dens.biB,S3_dens.biB)),3,3,byrow = F)
            
            
            matsS1[,,y]=(ann.matS1A+ann.matS1B)/2
            matsS2[,,y]=(ann.matS2A+ann.matS2B)/2
            matsS3[,,y]=(ann.matS3A+ann.matS3B)/2
            
            
            
            # Update densities
            densA=densS[1:3]
            densB=densS[4:6]
            densMet=c(densA,densB)
            
            S2_densA=S2_densS[1:3]
            S2_densB=S2_densS[4:6]
            S2_densMet=c(S2_densA,S2_densB)
            
            S3_densA=S3_densS[1:3]
            S3_densB=S3_densS[4:6]
            S3_densMet=c(S3_densA,S3_densB)
            
            ### Update bi densities
            
            dens.biA=sum(densA)
            dens.biB=sum(densB)
            
            
            S2_dens.biA=sum(S2_densA)
            S2_dens.biB=sum(S2_densB)
            
            S3_dens.biA=sum(S3_densA)
            S3_dens.biB=sum(S3_densB)
            
            
            
          }
          
          # get summary metrics
          
          # environment 
          
          envir.clim.change.sim=rbind(envir.clim.change.sim,data.frame(envW=env[,1],
                                                                       envS=env[,2],
                                                                       sim=x,
                                                                       scen=k,
                                                                       clim.change=clim.change.scen[cc],
                                                                       var.change=var.change[cp],
                                                                       cov.change=cov.change[cpx]))
          
          # total density per year and species
          densS1=aggregate(dens~year,data = temp[temp$species=="Prey"&temp$season%in%"winter",],sum)
          densS2=aggregate(dens~year,data = temp[temp$species=="Predator1"&temp$season%in%"winter",],sum)
          densS3=aggregate(dens~year,data = temp[temp$species=="Predator2"&temp$season%in%"winter",],sum)
          
          if(length(unique(densS1))<2|length(unique(densS2))<2|length(unique(densS3))<2){
            
            remove.sim.clim.change=rbind(remove.sim.clim.change,data.frame(sim=x,
                                                                           scen=k,
                                                                           clim.change=clim.change.scen[cc],
                                                                           var.change=var.change[cp],
                                                                           cov.change=cov.change[cpx]))
            print("Issue with Abundance")
            next
          }else{
            
            # Extniction probability
            ext1=ext2=ext3=NA
            
            if(any(densS1$dens<1)) ext1 <- 1 else ext1 <- 0
            
            if(any(densS2$dens<1)) ext2 <- 1 else ext2 <- 0
            
            if(any(densS3$dens<1)) ext3 <- 1 else ext3 <- 0
            
            # generation time and net reproductive rate
            
            GTS1=mean(apply(matsS1[,,densS1$year[densS1$year%in%densS1$year[densS1$dens>0]&!densS1$year%in%(101)]],3,generation.time))
            RRS1=mean(apply(matsS1[,,densS1$year[densS1$year%in%densS1$year[densS1$dens>0]&!densS1$year%in%(101)]],3,net.reproductive.rate))
            
            GTS2=mean(apply(matsS2[,,densS2$year[densS2$year%in%densS2$year[densS2$dens>0]&!densS2$year%in%(101)]],3,generation.time))
            RRS2=mean(apply(matsS2[,,densS2$year[densS2$year%in%densS2$year[densS2$dens>0]&!densS2$year%in%(101)]],3,net.reproductive.rate))
            
            GTS3=mean(apply(matsS3[,,densS3$year[densS3$year%in%densS3$year[densS3$dens>0]&!densS3$year%in%(101)]],3,generation.time))
            RRS3=mean(apply(matsS3[,,densS3$year[densS3$year%in%densS3$year[densS3$dens>0]&!densS3$year%in%(101)]],3,net.reproductive.rate))
            
            # Mean total abundance (before extinction)
            
            mu.densS1=mean(densS1$dens[densS1$dens>0])
            mu.densS2=mean(densS2$dens[densS2$dens>0])
            mu.densS3=mean(densS3$dens[densS3$dens>0])
            
            # CV in total abundance (before extinction)
            
            cv.densS1=CV(densS1$dens[densS1$dens>0])
            cv.densS2=CV(densS2$dens[densS2$dens>0])
            cv.densS3=CV(densS3$dens[densS3$dens>0])
            
            # mean stage-specific abudance
            
            densS1_stage=aggregate(dens~stage+year,data = temp[temp$species=="Prey"&temp$year%in%densS1$year[densS1$dens>0]&temp$season%in%"winter",],sum)
            mu.densS1_stage=aggregate(dens~stage,data = densS1_stage,mean)
            # CV
            cv.densS1_stage=aggregate(dens~stage,data = densS1_stage,CV)
            
            densS2_stage=aggregate(dens~stage+year,data = temp[temp$species=="Predator1"&temp$year%in%densS2$year[densS2$dens>0]&temp$season%in%"winter",],sum)
            mu.densS2_stage=aggregate(dens~stage,data = densS2_stage,mean)
            # CV
            cv.densS2_stage=aggregate(dens~stage,data = densS2_stage,CV)
            
            densS3_stage=aggregate(dens~stage+year,data = temp[temp$species=="Predator2"&temp$year%in%densS3$year[densS3$dens>0]&temp$season%in%"winter",],sum)
            mu.densS3_stage=aggregate(dens~stage,data = densS3_stage,mean)
            # CV
            cv.densS3_stage=aggregate(dens~stage,data = densS3_stage,CV)
            
            # seasonality in total abundance
            
            densS1v2=aggregate(dens~season+year,data = temp[temp$species=="Prey",],sum)
            
            densS1v2=densS1v2[densS1v2$year%in%densS1$year[densS1$dens>0]&!densS1v2$year%in%(years+1),]
            
            densS2v2=aggregate(dens~season+year,data = temp[temp$species=="Predator1",],sum)
            
            densS2v2=densS2v2[densS2v2$year%in%densS2$year[densS2$dens>0]&!densS2v2$year%in%(years+1),]
            
            densS3v2=aggregate(dens~season+year,data = temp[temp$species=="Predator2",],sum)
            
            densS3v2=densS3v2[densS3v2$year%in%densS3$year[densS3$dens>0]&!densS3v2$year%in%(years+1),]
            
            
            ts.S1=ts(densS1v2$dens,frequency = 2,start = c(1,1))
            
            seas.S1=decompose(ts.S1)$figure[2]
            
            ts.S2=ts(densS2v2$dens,frequency = 2,start = c(1,1))
            
            seas.S2=decompose(ts.S2)$figure[2]
            
            ts.S3=ts(densS3v2$dens,frequency = 2,start = c(1,1))
            
            seas.S3=decompose(ts.S3)$figure[2]
            
            
            # Seasonality in stage-specific abundance
            
            densS1v2_R=temp[temp$species=="Prey"&temp$stage%in%"R"&temp$year%in%densS1$year[densS1$dens>0]&!temp$year%in%(years+1),]
            densS1v2_N=temp[temp$species=="Prey"&temp$stage%in%"N"&temp$year%in%densS1$year[densS1$dens>0]&!temp$year%in%(years+1),]
            
            ts.S1=ts(densS1v2_R$dens,frequency = 2,start = c(1,1))
            
            seas.S1_R=decompose(ts.S1)$figure[2]
            
            ts.S1=ts(densS1v2_N$dens,frequency = 2,start = c(1,1))
            
            seas.S1_N=decompose(ts.S1)$figure[2]
            
            densS2v2_R=temp[temp$species=="Predator1"&temp$stage%in%"R"&temp$year%in%densS2$year[densS2$dens>0]&!temp$year%in%(years+1),]
            densS2v2_N=temp[temp$species=="Predator1"&temp$stage%in%"N"&temp$year%in%densS2$year[densS2$dens>0]&!temp$year%in%(years+1),]
            
            ts.S2=ts(densS2v2_R$dens,frequency = 2,start = c(1,1))
            
            seas.S2_R=decompose(ts.S2)$figure[2]
            
            ts.S2=ts(densS2v2_N$dens,frequency = 2,start = c(1,1))
            
            seas.S2_N=decompose(ts.S2)$figure[2]
            
            densS3v2_R=temp[temp$species=="Predator2"&temp$stage%in%"R"&temp$year%in%densS3$year[densS3$dens>0]&!temp$year%in%(years+1),]
            densS3v2_N=temp[temp$species=="Predator2"&temp$stage%in%"N"&temp$year%in%densS3$year[densS3$dens>0]&!temp$year%in%(years+1),]
            
            ts.S3=ts(densS3v2_R$dens,frequency = 2,start = c(1,1))
            
            seas.S3_R=decompose(ts.S3)$figure[2]
            
            ts.S3=ts(densS3v2_N$dens,frequency = 2,start = c(1,1))
            
            seas.S3_N=decompose(ts.S3)$figure[2]
            
            clim.change.sim=rbind(clim.change.sim,data.frame(species=c("Prey","Predator1","Predator2"),
                                                             ext=c(ext1,ext2,ext3),
                                                             Dens.mu=c(mu.densS1,mu.densS2,mu.densS3),
                                                             Dens.cv=c(cv.densS1,cv.densS2,cv.densS3),
                                                             seas=c(seas.S1,seas.S2,seas.S3),
                                                             GT=c(GTS1,GTS2,GTS3),
                                                             RR=c(RRS1,RRS2,RRS3),
                                                             sim=x,
                                                             scen=k,
                                                             clim.change=clim.change.scen[cc],
                                                             var.change=var.change[cp],
                                                             cov.change=cov.change[cpx]))
            
            clim.change.sim.stage=rbind(clim.change.sim.stage,data.frame(species=rep(c("Prey","Predator1","Predator2"),each=3),
                                                                         Dens.mu=c(mu.densS1_stage,mu.densS2_stage,mu.densS3_stage),
                                                                         Dens.cv=c(cv.densS1_stage,cv.densS2_stage,cv.densS3_stage),
                                                                         seas=c(0,seas.S1_N,seas.S1_R,0,seas.S2_N,seas.S2_R,0,seas.S3_N,seas.S3_R),
                                                                         sim=x,
                                                                         scen=k,
                                                                         clim.change=clim.change.scen[cc],
                                                                         var.change=var.change[cp],
                                                                         cov.change=cov.change[cpx]))
            
            
            
          }
          
        }
      }
    }
  }
}


nrow(clim.change.sim[clim.change.sim$ext%in%1,]) # not much extinction

# save output
write.csv(envir.clim.change.sim,"envir.clim.change.sim.csv",row.names = F)

write.csv(clim.change.sim,"clim.change.sim.csv",row.names = F)
write.csv(clim.change.sim.stage,"clim.change.sim.stage.csv",row.names = F)
write.csv(remove.sim.clim.change,"remove.sim.csv",row.names = F)

