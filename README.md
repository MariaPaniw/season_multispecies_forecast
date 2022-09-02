# season_multispecies_forecast

---
title: "Multispecies seasonal demographic models"

author: "Maria Paniw"

date: "01/09/2022"

---

This framework allows users to construct and project seasonal demographic metapopulation models for thee interacting species. 

## 1. Species

**S1** - a simulated prey species

**S2** - a simulated predator (Pred1)

**S3** - a simulated predator (Pred2) and competitor with S2

S1 (prey) is a relatively short-lived prey species and responds stronger to the environment than the other two species. S2 (predator 1) S3 (predator 2) predate on S1 and compete for this resource. S2 has a longer generation time than S3, is relatively more dependent on S1; is a stronger competitor (has a more negative effect on vital rates of S3 than vice versa); and exploits the prey more (stronger negative effect on S1 than in the case of S3). The environment affects all vital rates positively, except for Tj. The effects of predation are negative, except for Tj and Tn - speeding up transitions. All demographic rates of S2 and S3 are positively affected by prey density - except dispersal (less likely when food is abundant).


## 2. Demographic rates
For each of the three species, we parameterize a seasonal stage-structured population model at two sites (A and B). We determined demographic rates (survival, recruitment, and dispersal) for three stages across two discrete growing seasons: unproductive (winter) and productive (summer) season. The three stages are juveniles (J), non-reproductive adults (NR) and reproductive adults (R). Individuals in each stage survive (S), transition between stages (T), and reproductive adults also produce offspring (B). Non-reproductive adults also have a probability of dispersing (D) in order to reproduce, and the ones that disperse must survive (SD) to become established at the new site. 

We modeled each demographic rate in each season as a function of environmental states (env) parameterized as standardized deviations from seasonal means;  intraspecific density (dens); and interspecific density (dens.bi). 

These are demographic-rate functions for the three species:

```{r}

## S1 - PREY:

Sr_w <- function(pop,env,dens,dens.bi,dens.biS2){
  
  if(pop%in%"A"){
    
    mean=exp(A.mean.vec[1]+slope.env.vec[1]*env+slope.dens.vec[1]*dens+slope.bi.vec[1]*dens.bi+slope.biS2.vec[1]*dens.biS2+env.dens.diff.vec[1]*env*dens+dens.dens.bi.diff.vec[1]*dens*dens.bi)
  }else if(pop%in%"B"){
    
    mean=exp(A.mean.vec[1]+AB.diff.vec[1]+slope.env.vec[1]*env+slope.dens.vec[1]*dens+slope.bi.vec[1]*dens.bi+slope.biS2.vec[1]*dens.biS2+AB.env.diff.vec[1]*env+env.dens.diff.vec[1]*env*dens+dens.dens.bi.diff.vec[1]*dens*dens.bi+AB.dens.bi.diff.vec[1]*dens.bi)
  }
  
  return( mean/ (1 + mean))
}

Snr_w <- function(pop,env,dens,dens.bi,dens.biS2){
  
  if(pop%in%"A"){
    
    mean=exp(A.mean.vec[5]+slope.env.vec[5]*env+slope.dens.vec[5]*dens+slope.bi.vec[5]*dens.bi+slope.biS2.vec[5]*dens.biS2+env.dens.diff.vec[5]*env*dens+dens.dens.bi.diff.vec[5]*dens*dens.bi)
  }else if(pop%in%"B"){
    
    mean=exp(A.mean.vec[5]+AB.diff.vec[5]+slope.env.vec[5]*env+slope.dens.vec[5]*dens+slope.bi.vec[5]*dens.bi+slope.biS2.vec[5]*dens.biS2+AB.env.diff.vec[5]*env+env.dens.diff.vec[5]*env*dens+dens.dens.bi.diff.vec[5]*dens*dens.bi+AB.dens.bi.diff.vec[5]*dens.bi)
  }
  
  return( mean/ (1 + mean))
}

Sj_w <- function(pop,env,dens,dens.bi,dens.biS2){
  
  if(pop%in%"A"){
    
    mean=exp(A.mean.vec[6]+slope.env.vec[6]*env+slope.dens.vec[6]*dens+slope.bi.vec[6]*dens.bi+slope.biS2.vec[6]*dens.biS2+env.dens.diff.vec[6]*env*dens+dens.dens.bi.diff.vec[6]*dens*dens.bi)
  }else if(pop%in%"B"){
    
    mean=exp(A.mean.vec[6]+AB.diff.vec[6]+slope.env.vec[6]*env+slope.dens.vec[6]*dens+slope.bi.vec[6]*dens.bi+slope.biS2.vec[6]*dens.biS2+AB.env.diff.vec[6]*env+env.dens.diff.vec[6]*env*dens+dens.dens.bi.diff.vec[6]*dens*dens.bi+AB.dens.bi.diff.vec[6]*dens.bi)
  }
  
  return( mean/ (1 + mean))
}

Tj <- function(pop,env,dens,dens.bi,dens.biS2){
  
  if(pop%in%"A"){
    
    mean=exp(A.mean.vec[7]+slope.env.vec[7]*env+slope.dens.vec[7]*dens+slope.bi.vec[7]*dens.bi+slope.biS2.vec[7]*dens.biS2+env.dens.diff.vec[7]*env*dens+dens.dens.bi.diff.vec[7]*dens*dens.bi)
  }else if(pop%in%"B"){
    
    mean=exp(A.mean.vec[7]+AB.diff.vec[7]+slope.env.vec[7]*env+slope.dens.vec[7]*dens+slope.bi.vec[7]*dens.bi+slope.biS2.vec[7]*dens.biS2+AB.env.diff.vec[7]*env+env.dens.diff.vec[7]*env*dens+dens.dens.bi.diff.vec[7]*dens*dens.bi+AB.dens.bi.diff.vec[7]*dens.bi)
  }
  
  return( mean/ (1 + mean))
}

Tnr <- function(pop,env,dens,dens.bi,dens.biS2){
  
  if(pop%in%"A"){
    
    mean=exp(A.mean.vec[8]+slope.env.vec[8]*env+slope.dens.vec[8]*dens+slope.bi.vec[8]*dens.bi+slope.biS2.vec[8]*dens.biS2+env.dens.diff.vec[8]*env*dens+dens.dens.bi.diff.vec[8]*dens*dens.bi)
  }else if(pop%in%"B"){
    
    mean=exp(A.mean.vec[8]+AB.diff.vec[8]+slope.env.vec[8]*env+slope.dens.vec[8]*dens+slope.bi.vec[8]*dens.bi+slope.biS2.vec[8]*dens.biS2+AB.env.diff.vec[8]*env+env.dens.diff.vec[8]*env*dens+dens.dens.bi.diff.vec[8]*dens*dens.bi+AB.dens.bi.diff.vec[8]*dens.bi)
  }
  
  return( mean/ (1 + mean))
}

Tr <- function(pop,env,dens,dens.bi,dens.biS2){
  
  if(pop%in%"A"){
    
    mean=exp(A.mean.vec[9]+slope.env.vec[9]*env+slope.dens.vec[9]*dens+slope.bi.vec[9]*dens.bi+slope.biS2.vec[9]*dens.biS2+env.dens.diff.vec[9]*env*dens+dens.dens.bi.diff.vec[9]*dens*dens.bi)
  }else if(pop%in%"B"){
    
    mean=exp(A.mean.vec[9]+AB.diff.vec[9]+slope.env.vec[9]*env+slope.dens.vec[9]*dens+slope.bi.vec[9]*dens.bi+slope.biS2.vec[9]*dens.biS2+AB.env.diff.vec[9]*env+env.dens.diff.vec[9]*env*dens+dens.dens.bi.diff.vec[9]*dens*dens.bi+AB.dens.bi.diff.vec[9]*dens.bi)
  }
  
  return( mean/ (1 + mean))
}

Sr_s <- function(pop,env,dens,dens.bi,dens.biS2){
  
  if(pop%in%"A"){
    
    mean=exp(A.mean.vec[10]+slope.env.vec[10]*env+slope.dens.vec[10]*dens+slope.bi.vec[10]*dens.bi+slope.biS2.vec[10]*dens.biS2+env.dens.diff.vec[10]*env*dens+dens.dens.bi.diff.vec[10]*dens*dens.bi)
  }else if(pop%in%"B"){
    
    mean=exp(A.mean.vec[10]+AB.diff.vec[10]+slope.env.vec[10]*env+slope.dens.vec[10]*dens+slope.bi.vec[10]*dens.bi+slope.biS2.vec[10]*dens.biS2+AB.env.diff.vec[10]*env+env.dens.diff.vec[10]*env*dens+dens.dens.bi.diff.vec[10]*dens*dens.bi+AB.dens.bi.diff.vec[10]*dens.bi)
  }
  
  return( mean/ (1 + mean))
}

Snr_s <- function(pop,env,dens,dens.bi,dens.biS2){
  
  if(pop%in%"A"){
    
    mean=exp(A.mean.vec[11]+slope.env.vec[11]*env+slope.dens.vec[11]*dens+slope.bi.vec[11]*dens.bi+slope.biS2.vec[11]*dens.biS2+env.dens.diff.vec[11]*env*dens+dens.dens.bi.diff.vec[11]*dens*dens.bi)
  }else if(pop%in%"B"){
    
    mean=exp(A.mean.vec[11]+AB.diff.vec[11]+slope.env.vec[11]*env+slope.dens.vec[11]*dens+slope.bi.vec[11]*dens.bi+slope.biS2.vec[11]*dens.biS2+AB.env.diff.vec[11]*env+env.dens.diff.vec[11]*env*dens+dens.dens.bi.diff.vec[11]*dens*dens.bi+AB.dens.bi.diff.vec[11]*dens.bi)
  }
  
  return( mean/ (1 + mean))
}

D <- function(pop,env,dens,dens.bi,dens.biS2){
  
  if(pop%in%"A"){
    
    mean=exp(A.mean.vec[2]+slope.env.vec[2]*env+slope.dens.vec[2]*dens+slope.bi.vec[2]*dens.bi+slope.biS2.vec[2]*dens.biS2+env.dens.diff.vec[2]*env*dens+dens.dens.bi.diff.vec[2]*dens*dens.bi)
  }else if(pop%in%"B"){
    
    mean=exp(A.mean.vec[2]+AB.diff.vec[2]+slope.env.vec[2]*env+slope.dens.vec[2]*dens+slope.bi.vec[2]*dens.bi+slope.biS2.vec[2]*dens.biS2+AB.env.diff.vec[2]*env+env.dens.diff.vec[2]*env*dens+dens.dens.bi.diff.vec[2]*dens*dens.bi+AB.dens.bi.diff.vec[2]*dens.bi)
  }
  
  return( mean/ (1 + mean))
}

Ds <- function(pop,env,dens,dens.bi,dens.biS2){
  
  if(pop%in%"A"){
    
    mean=exp(A.mean.vec[3]+slope.env.vec[3]*env+slope.dens.vec[3]*dens+slope.bi.vec[3]*dens.bi+slope.biS2.vec[3]*dens.biS2+env.dens.diff.vec[3]*env*dens+dens.dens.bi.diff.vec[3]*dens*dens.bi)
  }else if(pop%in%"B"){
    
    mean=exp(A.mean.vec[3]+AB.diff.vec[3]+slope.env.vec[3]*env+slope.dens.vec[3]*dens+slope.bi.vec[3]*dens.bi+slope.biS2.vec[3]*dens.biS2+AB.env.diff.vec[3]*env+env.dens.diff.vec[3]*env*dens+dens.dens.bi.diff.vec[3]*dens*dens.bi+AB.dens.bi.diff.vec[3]*dens.bi)
  }
  
  return( mean/ (1 + mean))
}

B <- function(pop,env,dens,dens.bi,dens.biS2){
  
  if(pop%in%"A"){
    
    mean=exp(A.mean.vec[4]+slope.env.vec[4]*env+slope.dens.vec[4]*dens+slope.bi.vec[4]*dens.bi+slope.biS2.vec[4]*dens.biS2+env.dens.diff.vec[4]*env*dens+dens.dens.bi.diff.vec[4]*dens*dens.bi)
  }else if(pop%in%"B"){
    
    mean=exp(A.mean.vec[4]+AB.diff.vec[4]+slope.env.vec[4]*env+slope.dens.vec[4]*dens+slope.bi.vec[4]*dens.bi+slope.biS2.vec[4]*dens.biS2+AB.env.diff.vec[4]*env+env.dens.diff.vec[4]*env*dens+dens.dens.bi.diff.vec[4]*dens*dens.bi+AB.dens.bi.diff.vec[4]*dens.bi)
  }
  
  
  return(mean)
}


# vec permutation approach
Bs=matrix(0,3*2,3*2)
Ms=matrix(0,3*2,3*2)

diag(Ms)=1

Ps=matrix(0,3*2,3*2)

for(i in 1:3){
  for(j in 1:2){
    E=matrix(0,3,2)
    E[i,j]=1
    Ps=Ps+(E%x%t(E))
  }
} 

Bw=matrix(0,3*2,3*2)

Mw=matrix(0,3*2,3*2)
diag(Mw)=1

Pw=matrix(0,3*2,3*2)

for(i in 1:3){
  for(j in 1:2){
    E=matrix(0,3,2)
    E[i,j]=1
    Pw=Pw+(E%x%t(E))
  }
}

## S2 - PREDATOR 1

### vital rate functions:
S2_Sr_w <- function(pop,env,dens,dens.bi,dens.biS2){
  
  dens.bi2=dens.bi^2
  if(pop%in%"A"){
    
    mean=exp(A.mean.vecS2[1]+slope.env.vecS2[1]*env+slope.dens.vecS2[1]*dens+slope.bi.vecS2[1]*dens.bi+slope.biS2.vecS2[1]*dens.biS2+slope.bi2.vecS2[1]*dens.bi2+env.dens.diff.vecS2[1]*env*dens+dens.dens.bi.diff.vecS2[1]*dens*dens.bi)
  }else if(pop%in%"B"){
    
    mean=exp(A.mean.vecS2[1]+AB.diff.vecS2[1]+slope.env.vecS2[1]*env+slope.dens.vecS2[1]*dens+slope.bi.vecS2[1]*dens.bi+slope.biS2.vecS2[1]*dens.biS2+slope.bi2.vecS2[1]*dens.bi2+AB.env.diff.vecS2[1]*env+env.dens.diff.vecS2[1]*env*dens+dens.dens.bi.diff.vecS2[1]*dens*dens.bi)
  }
  
  return( mean/ (1 + mean))
}

S2_Snr_w <- function(pop,env,dens,dens.bi,dens.biS2){
  dens.bi2=dens.bi^2
  if(pop%in%"A"){
    
    mean=exp(A.mean.vecS2[5]+slope.env.vecS2[5]*env+slope.dens.vecS2[5]*dens+slope.bi.vecS2[5]*dens.bi+slope.biS2.vecS2[5]*dens.biS2+slope.bi2.vecS2[5]*dens.bi2+env.dens.diff.vecS2[5]*env*dens+dens.dens.bi.diff.vecS2[5]*dens*dens.bi)
  }else if(pop%in%"B"){
    
    mean=exp(A.mean.vecS2[5]+AB.diff.vecS2[5]+slope.env.vecS2[5]*env+slope.dens.vecS2[5]*dens+slope.bi.vecS2[5]*dens.bi+slope.biS2.vecS2[5]*dens.biS2+slope.bi2.vecS2[5]*dens.bi2+AB.env.diff.vecS2[5]*env+env.dens.diff.vecS2[5]*env*dens +dens.dens.bi.diff.vecS2[5]*dens*dens.bi)
  }
  
  return( mean/ (1 + mean))
}

S2_Sj_w <- function(pop,env,dens,dens.bi,dens.biS2){
  dens.bi2=dens.bi^2
  if(pop%in%"A"){
    
    mean=exp(A.mean.vecS2[6]+slope.env.vecS2[6]*env+slope.dens.vecS2[6]*dens+slope.bi.vecS2[6]*dens.bi+slope.biS2.vecS2[6]*dens.biS2+slope.bi2.vecS2[6]*dens.bi2+env.dens.diff.vecS2[6]*env*dens+dens.dens.bi.diff.vecS2[6]*dens*dens.bi)
  }else if(pop%in%"B"){
    
    mean=exp(A.mean.vecS2[6]+AB.diff.vecS2[6]+slope.env.vecS2[6]*env+slope.dens.vecS2[6]*dens+slope.bi.vecS2[6]*dens.bi+slope.biS2.vecS2[6]*dens.biS2+slope.bi2.vecS2[6]*dens.bi2+AB.env.diff.vecS2[6]*env+env.dens.diff.vecS2[6]*env*dens+dens.dens.bi.diff.vecS2[6]*dens*dens.bi)
  }
  
  return( mean/ (1 + mean))
}

S2_Tj <- function(pop,env,dens,dens.bi,dens.biS2){
  dens.bi2=dens.bi^2
  if(pop%in%"A"){
    
    mean=exp(A.mean.vecS2[7]+slope.env.vecS2[7]*env+slope.dens.vecS2[7]*dens+slope.bi.vecS2[7]*dens.bi+slope.biS2.vecS2[7]*dens.biS2+slope.bi2.vecS2[7]*dens.bi2+env.dens.diff.vecS2[7]*env*dens+dens.dens.bi.diff.vecS2[7]*dens*dens.bi)
  }else if(pop%in%"B"){
    
    mean=exp(A.mean.vecS2[7]+AB.diff.vecS2[7]+slope.env.vecS2[7]*env+slope.dens.vecS2[7]*dens+slope.bi.vecS2[7]*dens.bi+slope.biS2.vecS2[7]*dens.biS2+slope.bi2.vecS2[7]*dens.bi2+AB.env.diff.vecS2[7]*env+env.dens.diff.vecS2[7]*env*dens+dens.dens.bi.diff.vecS2[7]*dens*dens.bi)
  }
  
  return( mean/ (1 + mean))
}

S2_Tnr <- function(pop,env,dens,dens.bi,dens.biS2){
  dens.bi2=dens.bi^2
  if(pop%in%"A"){
    
    mean=exp(A.mean.vecS2[8]+slope.env.vecS2[8]*env+slope.dens.vecS2[8]*dens+slope.bi.vecS2[8]*dens.bi+slope.biS2.vecS2[8]*dens.biS2+slope.bi2.vecS2[8]*dens.bi2+env.dens.diff.vecS2[8]*env*dens+dens.dens.bi.diff.vecS2[8]*dens*dens.bi)
  }else if(pop%in%"B"){
    
    mean=exp(A.mean.vecS2[8]+AB.diff.vecS2[8]+slope.env.vecS2[8]*env+slope.dens.vecS2[8]*dens+slope.bi.vecS2[8]*dens.bi+slope.biS2.vecS2[8]*dens.biS2+slope.bi2.vecS2[8]*dens.bi2+AB.env.diff.vecS2[8]*env+env.dens.diff.vecS2[8]*env*dens+dens.dens.bi.diff.vecS2[8]*dens*dens.bi)
  }
  
  return( mean/ (1 + mean))
}

S2_Tr <- function(pop,env,dens,dens.bi,dens.biS2){
  dens.bi2=dens.bi^2
  if(pop%in%"A"){
    
    mean=exp(A.mean.vecS2[9]+slope.env.vecS2[9]*env+slope.dens.vecS2[9]*dens+slope.bi.vecS2[9]*dens.bi+slope.biS2.vecS2[9]*dens.biS2+slope.bi2.vecS2[9]*dens.bi2+env.dens.diff.vecS2[9]*env*dens+dens.dens.bi.diff.vecS2[9]*dens*dens.bi)
  }else if(pop%in%"B"){
    
    mean=exp(A.mean.vecS2[9]+AB.diff.vecS2[9]+slope.env.vecS2[9]*env+slope.dens.vecS2[9]*dens+slope.bi.vecS2[9]*dens.bi+slope.biS2.vecS2[9]*dens.biS2+slope.bi2.vecS2[9]*dens.bi2+AB.env.diff.vecS2[9]*env+env.dens.diff.vecS2[9]*env*dens+dens.dens.bi.diff.vecS2[9]*dens*dens.bi)
  }
  
  return( mean/ (1 + mean))
}

S2_Sr_s <- function(pop,env,dens,dens.bi,dens.biS2){
  dens.bi2=dens.bi^2
  if(pop%in%"A"){
    
    mean=exp(A.mean.vecS2[10]+slope.env.vecS2[10]*env+slope.dens.vecS2[10]*dens+slope.bi.vecS2[10]*dens.bi+slope.biS2.vecS2[10]*dens.biS2+slope.bi2.vecS2[10]*dens.bi2+env.dens.diff.vecS2[10]*env*dens+dens.dens.bi.diff.vecS2[10]*dens*dens.bi)
  }else if(pop%in%"B"){
    
    mean=exp(A.mean.vecS2[10]+AB.diff.vecS2[10]+slope.env.vecS2[10]*env+slope.dens.vecS2[10]*dens+slope.bi.vecS2[10]*dens.bi+slope.biS2.vecS2[10]*dens.biS2+slope.bi2.vecS2[10]*dens.bi2+AB.env.diff.vecS2[10]*env+env.dens.diff.vecS2[10]*env*dens+dens.dens.bi.diff.vecS2[10]*dens*dens.bi)
  }
  
  return( mean/ (1 + mean))
}

S2_Snr_s <- function(pop,env,dens,dens.bi,dens.biS2){
  dens.bi2=dens.bi^2
  if(pop%in%"A"){
    
    mean=exp(A.mean.vecS2[11]+slope.env.vecS2[11]*env+slope.dens.vecS2[11]*dens+slope.bi.vecS2[11]*dens.bi+slope.biS2.vecS2[11]*dens.biS2+slope.bi2.vecS2[11]*dens.bi2+env.dens.diff.vecS2[11]*env*dens+dens.dens.bi.diff.vecS2[11]*dens*dens.bi)
  }else if(pop%in%"B"){
    
    mean=exp(A.mean.vecS2[11]+AB.diff.vecS2[11]+slope.env.vecS2[11]*env+slope.dens.vecS2[11]*dens+slope.bi.vecS2[11]*dens.bi+slope.biS2.vecS2[11]*dens.biS2+slope.bi2.vecS2[11]*dens.bi2+AB.env.diff.vecS2[11]*env+env.dens.diff.vecS2[11]*env*dens+dens.dens.bi.diff.vecS2[11]*dens*dens.bi)
  }
  
  return( mean/ (1 + mean))
}

S2_D <- function(pop,env,dens,dens.bi,dens.biS2){
  dens.bi2=dens.bi^2
  if(pop%in%"A"){
    
    mean=exp(A.mean.vecS2[2]+slope.env.vecS2[2]*env+slope.dens.vecS2[2]*dens+slope.bi.vecS2[2]*dens.bi+slope.biS2.vecS2[2]*dens.biS2+slope.bi2.vecS2[2]*dens.bi2+env.dens.diff.vecS2[2]*env*dens+dens.dens.bi.diff.vecS2[2]*dens*dens.bi)
  }else if(pop%in%"B"){
    
    mean=exp(A.mean.vecS2[2]+AB.diff.vecS2[2]+slope.env.vecS2[2]*env+slope.dens.vecS2[2]*dens+slope.bi.vecS2[2]*dens.bi+slope.biS2.vecS2[2]*dens.biS2+slope.bi2.vecS2[2]*dens.bi2+AB.env.diff.vecS2[2]*env+env.dens.diff.vecS2[2]*env*dens+dens.dens.bi.diff.vecS2[2]*dens*dens.bi)
  }
  
  return( mean/ (1 + mean))
}

S2_Ds <- function(pop,env,dens,dens.bi,dens.biS2){
  dens.bi2=dens.bi^2
  if(pop%in%"A"){
    
    mean=exp(A.mean.vecS2[3]+slope.env.vecS2[3]*env+slope.dens.vecS2[3]*dens+slope.bi.vecS2[3]*dens.bi+slope.biS2.vecS2[3]*dens.biS2+slope.bi2.vecS2[3]*dens.bi2+env.dens.diff.vecS2[3]*env*dens+dens.dens.bi.diff.vecS2[3]*dens*dens.bi)
  }else if(pop%in%"B"){
    
    mean=exp(A.mean.vecS2[3]+AB.diff.vecS2[3]+slope.env.vecS2[3]*env+slope.dens.vecS2[3]*dens+slope.bi.vecS2[3]*dens.bi+slope.biS2.vecS2[3]*dens.biS2+slope.bi2.vecS2[3]*dens.bi2+AB.env.diff.vecS2[3]*env+env.dens.diff.vecS2[3]*env*dens+dens.dens.bi.diff.vecS2[3]*dens*dens.bi)
  }
  
  return( mean/ (1 + mean))
}

S2_B <- function(pop,env,dens,dens.bi,dens.biS2){
  dens.bi2=dens.bi^2
  if(pop%in%"A"){
    
    mean=exp(A.mean.vecS2[4]+slope.env.vecS2[4]*env+slope.dens.vecS2[4]*dens+slope.bi.vecS2[4]*dens.bi+slope.biS2.vecS2[4]*dens.biS2+slope.bi2.vecS2[4]*dens.bi2+env.dens.diff.vecS2[4]*env*dens+dens.dens.bi.diff.vecS2[4]*dens*dens.bi)
  }else if(pop%in%"B"){
    
    mean=exp(A.mean.vecS2[4]+AB.diff.vecS2[4]+slope.env.vecS2[4]*env+slope.dens.vecS2[4]*dens+slope.bi.vecS2[4]*dens.bi+slope.biS2.vecS2[4]*dens.biS2+slope.bi2.vecS2[4]*dens.bi2+AB.env.diff.vecS2[4]*env+env.dens.diff.vecS2[4]*env*dens+dens.dens.bi.diff.vecS2[4]*dens*dens.bi)
  }
  
  if(mean>10) mean<-10 #set cap on birth
  return(mean)
}


# vec permutation approach
S2_Bs=matrix(0,3*2,3*2)
S2_Ms=matrix(0,3*2,3*2)

diag(S2_Ms)=1

S2_Ps=matrix(0,3*2,3*2)

for(i in 1:3){
  for(j in 1:2){
    E=matrix(0,3,2)
    E[i,j]=1
    S2_Ps=S2_Ps+(E%x%t(E))
  }
} 

S2_Bw=matrix(0,3*2,3*2)

S2_Mw=matrix(0,3*2,3*2)
diag(S2_Mw)=1

S2_Pw=matrix(0,3*2,3*2)

for(i in 1:3){
  for(j in 1:2){
    E=matrix(0,3,2)
    E[i,j]=1
    S2_Pw=S2_Pw+(E%x%t(E))
  }
}


#### S3 - PREDATOR 2

### vital rate functions:
S3_Sr_w <- function(pop,env,dens,dens.bi,dens.biS2){
  
  dens.bi2=dens.bi^2
  if(pop%in%"A"){
    
    mean=exp(A.mean.vecS3[1]+slope.env.vecS3[1]*env+slope.dens.vecS3[1]*dens+slope.bi.vecS3[1]*dens.bi+slope.biS2.vecS3[1]*dens.biS2+slope.bi2.vecS3[1]*dens.bi2+env.dens.diff.vecS3[1]*env*dens+dens.dens.bi.diff.vecS3[1]*dens*dens.bi)
  }else if(pop%in%"B"){
    
    mean=exp(A.mean.vecS3[1]+AB.diff.vecS3[1]+slope.env.vecS3[1]*env+slope.dens.vecS3[1]*dens+slope.bi.vecS3[1]*dens.bi+slope.biS2.vecS3[1]*dens.biS2+slope.bi2.vecS3[1]*dens.bi2+AB.env.diff.vecS3[1]*env+env.dens.diff.vecS3[1]*env*dens+dens.dens.bi.diff.vecS3[1]*dens*dens.bi)
  }
  
  return( mean/ (1 + mean))
}

S3_Snr_w <- function(pop,env,dens,dens.bi,dens.biS2){
  dens.bi2=dens.bi^2
  if(pop%in%"A"){
    
    mean=exp(A.mean.vecS3[5]+slope.env.vecS3[5]*env+slope.dens.vecS3[5]*dens+slope.bi.vecS3[5]*dens.bi+slope.biS2.vecS3[5]*dens.biS2+slope.bi2.vecS3[5]*dens.bi2+env.dens.diff.vecS3[5]*env*dens+dens.dens.bi.diff.vecS3[5]*dens*dens.bi)
  }else if(pop%in%"B"){
    
    mean=exp(A.mean.vecS3[5]+AB.diff.vecS3[5]+slope.env.vecS3[5]*env+slope.dens.vecS3[5]*dens+slope.bi.vecS3[5]*dens.bi+slope.biS2.vecS3[5]*dens.biS2+slope.bi2.vecS3[5]*dens.bi2+AB.env.diff.vecS3[5]*env+env.dens.diff.vecS3[5]*env*dens +dens.dens.bi.diff.vecS3[5]*dens*dens.bi)
  }
  
  return( mean/ (1 + mean))
}

S3_Sj_w <- function(pop,env,dens,dens.bi,dens.biS2){
  dens.bi2=dens.bi^2
  if(pop%in%"A"){
    
    mean=exp(A.mean.vecS3[6]+slope.env.vecS3[6]*env+slope.dens.vecS3[6]*dens+slope.bi.vecS3[6]*dens.bi+slope.biS2.vecS3[6]*dens.biS2+slope.bi2.vecS3[6]*dens.bi2+env.dens.diff.vecS3[6]*env*dens+dens.dens.bi.diff.vecS3[6]*dens*dens.bi)
  }else if(pop%in%"B"){
    
    mean=exp(A.mean.vecS3[6]+AB.diff.vecS3[6]+slope.env.vecS3[6]*env+slope.dens.vecS3[6]*dens+slope.bi.vecS3[6]*dens.bi+slope.biS2.vecS3[6]*dens.biS2+slope.bi2.vecS3[6]*dens.bi2+AB.env.diff.vecS3[6]*env+env.dens.diff.vecS3[6]*env*dens+dens.dens.bi.diff.vecS3[6]*dens*dens.bi)
  }
  
  return( mean/ (1 + mean))
}

S3_Tj <- function(pop,env,dens,dens.bi,dens.biS2){
  dens.bi2=dens.bi^2
  if(pop%in%"A"){
    
    mean=exp(A.mean.vecS3[7]+slope.env.vecS3[7]*env+slope.dens.vecS3[7]*dens+slope.bi.vecS3[7]*dens.bi+slope.biS2.vecS3[7]*dens.biS2+slope.bi2.vecS3[7]*dens.bi2+env.dens.diff.vecS3[7]*env*dens+dens.dens.bi.diff.vecS3[7]*dens*dens.bi)
  }else if(pop%in%"B"){
    
    mean=exp(A.mean.vecS3[7]+AB.diff.vecS3[7]+slope.env.vecS3[7]*env+slope.dens.vecS3[7]*dens+slope.bi.vecS3[7]*dens.bi+slope.biS2.vecS3[7]*dens.biS2+slope.bi2.vecS3[7]*dens.bi2+AB.env.diff.vecS3[7]*env+env.dens.diff.vecS3[7]*env*dens+dens.dens.bi.diff.vecS3[7]*dens*dens.bi)
  }
  
  return( mean/ (1 + mean))
}

S3_Tnr <- function(pop,env,dens,dens.bi,dens.biS2){
  dens.bi2=dens.bi^2
  if(pop%in%"A"){
    
    mean=exp(A.mean.vecS3[8]+slope.env.vecS3[8]*env+slope.dens.vecS3[8]*dens+slope.bi.vecS3[8]*dens.bi+slope.biS2.vecS3[8]*dens.biS2+slope.bi2.vecS3[8]*dens.bi2+env.dens.diff.vecS3[8]*env*dens+dens.dens.bi.diff.vecS3[8]*dens*dens.bi)
  }else if(pop%in%"B"){
    
    mean=exp(A.mean.vecS3[8]+AB.diff.vecS3[8]+slope.env.vecS3[8]*env+slope.dens.vecS3[8]*dens+slope.bi.vecS3[8]*dens.bi+slope.biS2.vecS3[8]*dens.biS2+slope.bi2.vecS3[8]*dens.bi2+AB.env.diff.vecS3[8]*env+env.dens.diff.vecS3[8]*env*dens+dens.dens.bi.diff.vecS3[8]*dens*dens.bi)
  }
  
  return( mean/ (1 + mean))
}

S3_Tr <- function(pop,env,dens,dens.bi,dens.biS2){
  dens.bi2=dens.bi^2
  if(pop%in%"A"){
    
    mean=exp(A.mean.vecS3[9]+slope.env.vecS3[9]*env+slope.dens.vecS3[9]*dens+slope.bi.vecS3[9]*dens.bi+slope.biS2.vecS3[9]*dens.biS2+slope.bi2.vecS3[9]*dens.bi2+env.dens.diff.vecS3[9]*env*dens+dens.dens.bi.diff.vecS3[9]*dens*dens.bi)
  }else if(pop%in%"B"){
    
    mean=exp(A.mean.vecS3[9]+AB.diff.vecS3[9]+slope.env.vecS3[9]*env+slope.dens.vecS3[9]*dens+slope.bi.vecS3[9]*dens.bi+slope.biS2.vecS3[9]*dens.biS2+slope.bi2.vecS3[9]*dens.bi2+AB.env.diff.vecS3[9]*env+env.dens.diff.vecS3[9]*env*dens+dens.dens.bi.diff.vecS3[9]*dens*dens.bi)
  }
  
  return( mean/ (1 + mean))
}

S3_Sr_s <- function(pop,env,dens,dens.bi,dens.biS2){
  dens.bi2=dens.bi^2
  if(pop%in%"A"){
    
    mean=exp(A.mean.vecS3[10]+slope.env.vecS3[10]*env+slope.dens.vecS3[10]*dens+slope.bi.vecS3[10]*dens.bi+slope.biS2.vecS3[10]*dens.biS2+slope.bi2.vecS3[10]*dens.bi2+env.dens.diff.vecS3[10]*env*dens+dens.dens.bi.diff.vecS3[10]*dens*dens.bi)
  }else if(pop%in%"B"){
    
    mean=exp(A.mean.vecS3[10]+AB.diff.vecS3[10]+slope.env.vecS3[10]*env+slope.dens.vecS3[10]*dens+slope.bi.vecS3[10]*dens.bi+slope.biS2.vecS3[10]*dens.biS2+slope.bi2.vecS3[10]*dens.bi2+AB.env.diff.vecS3[10]*env+env.dens.diff.vecS3[10]*env*dens+dens.dens.bi.diff.vecS3[10]*dens*dens.bi)
  }
  
  return( mean/ (1 + mean))
}

S3_Snr_s <- function(pop,env,dens,dens.bi,dens.biS2){
  dens.bi2=dens.bi^2
  if(pop%in%"A"){
    
    mean=exp(A.mean.vecS3[11]+slope.env.vecS3[11]*env+slope.dens.vecS3[11]*dens+slope.bi.vecS3[11]*dens.bi+slope.biS2.vecS3[11]*dens.biS2+slope.bi2.vecS3[11]*dens.bi2+env.dens.diff.vecS3[11]*env*dens+dens.dens.bi.diff.vecS3[11]*dens*dens.bi)
  }else if(pop%in%"B"){
    
    mean=exp(A.mean.vecS3[11]+AB.diff.vecS3[11]+slope.env.vecS3[11]*env+slope.dens.vecS3[11]*dens+slope.bi.vecS3[11]*dens.bi+slope.biS2.vecS3[11]*dens.biS2+slope.bi2.vecS3[11]*dens.bi2+AB.env.diff.vecS3[11]*env+env.dens.diff.vecS3[11]*env*dens+dens.dens.bi.diff.vecS3[11]*dens*dens.bi)
  }
  
  return( mean/ (1 + mean))
}

S3_D <- function(pop,env,dens,dens.bi,dens.biS2){
  dens.bi2=dens.bi^2
  if(pop%in%"A"){
    
    mean=exp(A.mean.vecS3[2]+slope.env.vecS3[2]*env+slope.dens.vecS3[2]*dens+slope.bi.vecS3[2]*dens.bi+slope.biS2.vecS3[2]*dens.biS2+slope.bi2.vecS3[2]*dens.bi2+env.dens.diff.vecS3[2]*env*dens+dens.dens.bi.diff.vecS3[2]*dens*dens.bi)
  }else if(pop%in%"B"){
    
    mean=exp(A.mean.vecS3[2]+AB.diff.vecS3[2]+slope.env.vecS3[2]*env+slope.dens.vecS3[2]*dens+slope.bi.vecS3[2]*dens.bi+slope.biS2.vecS3[2]*dens.biS2+slope.bi2.vecS3[2]*dens.bi2+AB.env.diff.vecS3[2]*env+env.dens.diff.vecS3[2]*env*dens+dens.dens.bi.diff.vecS3[2]*dens*dens.bi)
  }
  
  return( mean/ (1 + mean))
}

S3_Ds <- function(pop,env,dens,dens.bi,dens.biS2){
  dens.bi2=dens.bi^2
  if(pop%in%"A"){
    
    mean=exp(A.mean.vecS3[3]+slope.env.vecS3[3]*env+slope.dens.vecS3[3]*dens+slope.bi.vecS3[3]*dens.bi+slope.biS2.vecS3[3]*dens.biS2+slope.bi2.vecS3[3]*dens.bi2+env.dens.diff.vecS3[3]*env*dens+dens.dens.bi.diff.vecS3[3]*dens*dens.bi)
  }else if(pop%in%"B"){
    
    mean=exp(A.mean.vecS3[3]+AB.diff.vecS3[3]+slope.env.vecS3[3]*env+slope.dens.vecS3[3]*dens+slope.bi.vecS3[3]*dens.bi+slope.biS2.vecS3[3]*dens.biS2+slope.bi2.vecS3[3]*dens.bi2+AB.env.diff.vecS3[3]*env+env.dens.diff.vecS3[3]*env*dens+dens.dens.bi.diff.vecS3[3]*dens*dens.bi)
  }
  
  return( mean/ (1 + mean))
}

S3_B <- function(pop,env,dens,dens.bi,dens.biS2){
  dens.bi2=dens.bi^2
  if(pop%in%"A"){
    
    mean=exp(A.mean.vecS3[4]+slope.env.vecS3[4]*env+slope.dens.vecS3[4]*dens+slope.bi.vecS3[4]*dens.bi+slope.biS2.vecS3[4]*dens.biS2+slope.bi2.vecS3[4]*dens.bi2+env.dens.diff.vecS3[4]*env*dens+dens.dens.bi.diff.vecS3[4]*dens*dens.bi)
  }else if(pop%in%"B"){
    
    mean=exp(A.mean.vecS3[4]+AB.diff.vecS3[4]+slope.env.vecS3[4]*env+slope.dens.vecS3[4]*dens+slope.bi.vecS3[4]*dens.bi+slope.biS2.vecS3[4]*dens.biS2+slope.bi2.vecS3[4]*dens.bi2+AB.env.diff.vecS3[4]*env+env.dens.diff.vecS3[4]*env*dens+dens.dens.bi.diff.vecS3[4]*dens*dens.bi)
  }
  
  if(mean>20) mean<-20 #set cap on birth
  return(mean)
}


# vec permutation approach
S3_Bs=matrix(0,3*2,3*2)
S3_Ms=matrix(0,3*2,3*2)

diag(S3_Ms)=1

S3_Ps=matrix(0,3*2,3*2)

for(i in 1:3){
  for(j in 1:2){
    E=matrix(0,3,2)
    E[i,j]=1
    S3_Ps=S3_Ps+(E%x%t(E))
  }
} 

S3_Bw=matrix(0,3*2,3*2)

S3_Mw=matrix(0,3*2,3*2)
diag(S3_Mw)=1

S3_Pw=matrix(0,3*2,3*2)

for(i in 1:3){
  for(j in 1:2){
    E=matrix(0,3,2)
    E[i,j]=1
    S3_Pw=S3_Pw+(E%x%t(E))
  }
}

```


## 3. Constructing communities

By perturbing the parameters in the demographic models above, 48 different communities can be created (Fig. 1).

![communities](outline_communities_V2.png)


In all communities, effect of interspecific density is expressed relative to intraspecific density; as percent change in intraspecific density. Intraspecific density effects are negative, expect for dispersal. Survival or non-reproductive and reproductive is lower in unproductive season.

The communities differ as follows: 
 
**Life-history differences** (LH):

(a) Relatively higher average adult survival and lower recruitment for S2 (-LH)

(b) Relatively higher recruitment for S1 (+ B prey)

**Seasonality of average demography** (LS vs HS): 

Survival or non-reproductive and reproductive is lower in unproductive season, but difference are higher in “ high-seasonality” (HS) simulations compared to baseline low seasonality (LS)


**Interspecific interaction strength & seasonality**:

(a) Stronger effect of S1 density on S2 (higher specialization on S1) (+ Spec pred1)

(b) More negative effects of predators (S2 and S3) for summer survival in S1 + (a) (+ Spec pred1 + HS pred)

(c) More neagative effect of competition for summer survival in S2 and S3 + (a) + (b) (+ Spec pred1 + HS pred + HS comp)

Below is the code with the relevant parameters to construct the communities:

```{r}

###  Prey species

vr=c("Sr_w","D","Sd","B","Snr_w","Sj_w","Tj","Tnr","Tr","Sr_s","Snr_s")

scen=c("LS: Base","LS: - LH","LS: + B in prey","HS: Base","HS: - LH","HS: + B in prey",
       paste(c("LS: Base","LS: - LH","LS: + B in prey","HS: Base","HS: - LH","HS: + B in prey"),"+ Spec pred1"),
       paste(c("LS: Base","LS: - LH","LS: + B in prey","HS: Base","HS: - LH","HS: + B in prey"),"+ Spec pred1 + HS pred"),
       paste(c("LS: Base","LS: - LH","LS: + B in prey","HS: Base","HS: - LH","HS: + B in prey"),"+ Spec pred1 + HS pred/comp"))

scen=c(scen,paste(scen,"HH"))

# mean vital rates (LS - HS differences)
scen.Sr_w=rep(c(4,4,2,3,3,1),8)
scen.D=rep(-1,48)
scen.Sd=rep(0.3,48)
scen.B=rep(c(1.7,1.7,3,1.7,1.7,3),8)
scen.Snr_w=rep(c(3.5,3.5,1.8,3,3,1.5),8)
scen.Sj_w=rep(c(0.6,0.6,0.8,0.6,0.6,0.8),8)
scen.Tj=rep(c(1,1,1.2,1,1,1.2),8)
scen.Tnr=rep(c(1,1,1.2,1,1,1.2),8)
scen.Tr=rep(1.6,48)
scen.Sr_s=rep(c(5,5,3,6,6,8),8)
scen.Snr_s=rep(c(4,4,3,4.5,4.5,3.3),8)


# environment slope
env.scen.Sr_w=rep(c(0.15,0.15,0.2,0.15,0.15,0.2),8)
env.scen.D=rep(c(0.05,0.05,0.05,0.05,0.05,0.05),8)
env.scen.Sd=rep(c(0.15,0.15,0.2,0.15,0.15,0.2),8)
env.scen.B=rep(c(0.05,0.05,0.1,0.05,0.05,0.1),8)
env.scen.Snr_w=rep(c(0.2,0.2,0.2,0.2,0.2,0.2),8)
env.scen.Sj_w=rep(c(1.2,1.2,1.4,1.2,1.2,1.4),8)
env.scen.Tj=rep(c(-0.05,-0.05,-0.1,-0.05,-0.05,-0.1),8)
env.scen.Tnr=rep(c(0.1,0.1,0.1,0.1,0.1,0.1),8)
env.scen.Tr=rep(c(0.3,0.3,0.3,0.3,0.3,0.3),8)
env.scen.Sr_s=rep(c(0.1,0.1,0.2,0.1,0.1,0.2),8)
env.scen.Snr_s=rep(c(0.1,0.1,0.2,0.1,0.1,0.2),8)

# biotic interactions Predator 1

bi.scen.Sr_w=c(rep(-1.1,12),rep(-0.8,12),rep(-1.1,12),rep(-0.8,12))
bi.scen.D=rep(-2,48)
bi.scen.Sd=rep(-2,48)
bi.scen.B=rep(-1,48)
bi.scen.Snr_w=c(rep(-1,12),rep(-0.7,12),rep(-1,12),rep(-0.7,12))
bi.scen.Sj_w=rep(-1,48)
bi.scen.Tj=rep(1,48)
bi.scen.Tnr=rep(0.5,48)
bi.scen.Tr=rep(-1,48)
bi.scen.Sr_s=c(rep(-1,12),rep(-1.4,12),rep(-1,12),rep(-1.4,12))
bi.scen.Snr_s=c(rep(-1,12),rep(-1.3,12),rep(-1,12),rep(-1.3,12))

# biotic interactions Predator 2 (less effect than predotr 1)

biS2.scen.Sr_w=c(rep(-0.8,12),rep(-0.5,12),rep(-0.8,12),rep(-0.5,12))
biS2.scen.D=rep(-1,48)
biS2.scen.Sd=rep(-1,48)
biS2.scen.B=rep(-0.8,48)
biS2.scen.Snr_w=c(rep(-0.8,12),rep(-0.5,12),rep(-0.8,12),rep(-0.5,12))
biS2.scen.Sj_w=rep(-0.8,48)
biS2.scen.Tj=rep(0.8,48)
biS2.scen.Tnr=rep(0.01,48)
biS2.scen.Tr=rep(-0.8,48)
biS2.scen.Sr_s=c(rep(-0.8,12),rep(-1.1,12),rep(-0.8,12),rep(-1.1,12))
biS2.scen.Snr_s=c(rep(-0.8,12),rep(-1.1,12),rep(-0.8,12),rep(-1.1,12))

# dens dens interactions

dens.bi.scen.Sr_w=rep(0,48)
dens.bi.scen.D=rep(0,48)
dens.bi.scen.Sd=rep(0,48)
dens.bi.scen.B=rep(0,48)
dens.bi.scen.Snr_w=rep(0,48)
dens.bi.scen.Sj_w=rep(0,48)
dens.bi.scen.Tj=rep(0,48)
dens.bi.scen.Tnr=rep(0,48)
dens.bi.scen.Tr=rep(0,48)
dens.bi.scen.Sr_s=rep(0,48)
dens.bi.scen.Snr_s=rep(0,48)

## spatial differences

AB.scen.Sr_w=c(rep(0,24),rep(-0.1,24))
AB.scen.Sd=c(rep(0,24),rep(-0.1,24))
AB.scen.Snr_w=c(rep(0,24),rep(-0.1,24))
AB.scen.Sj_w=c(rep(0,24),rep(-0.1,24))
AB.scen.Sr_s=c(rep(0,24),rep(-0.1,24))
AB.scen.Snr_s=c(rep(0,24),rep(-0.1,24))

# spatial differences in environment slope

AB.env.scen.Sr_w=c(rep(0,24),rep(-0.1,24))
AB.env.scen.D=c(rep(0,24),rep(-0.1,24))
AB.env.scen.Sd=c(rep(0,24),rep(-0.1,24))
AB.env.scen.B=c(rep(0,24),rep(-0.1,24))
AB.env.scen.Snr_w=c(rep(0,24),rep(-0.1,24))
AB.env.scen.Sj_w=c(rep(0,24),rep(-0.1,24))
AB.env.scen.Tj=c(rep(0,24),rep(-0.1,24))
AB.env.scen.Tnr=c(rep(0,24),rep(-0.1,24))
AB.env.scen.Tr=c(rep(0,24),rep(-0.1,24))
AB.env.scen.Sr_s=c(rep(0,24),rep(-0.1,24))
AB.env.scen.Snr_s=c(rep(0,24),rep(-0.1,24))

###  Predator species 1 (bigger and more specialized)

# mean vital rates
S2.scen.Sr_w=rep(c(1.3,1.4,1.3,0.8,0.9,0.8),8)
S2.scen.D=rep(0.01,48)
S2.scen.Sd=rep(c(1.5,1,1.5,1.5,1,1.5),8)
S2.scen.B=rep(c(-1.5,-0.8,-1.5,-1.5,-0.8,-1.5),8)
S2.scen.Snr_w=rep(c(1.8,1.6,1.8,1.6,1.4,1.6),8)
S2.scen.Sj_w=rep(c(0.01,0.01,0.01,0.01,0.01,0.01),8)
S2.scen.Tj=rep(c(-2.8,0.0007,-2.8,-2.8,0.0007,-2.8),8)
S2.scen.Tnr=rep(c(0.5,0.5,0.5,0.5,0.5,0.5),8)
S2.scen.Tr=rep(c(1.5,1.5,1.5,1.5,1.5,1.5),8)
S2.scen.Sr_s=rep(c(2,2.1,2,2.5,2.6,2.5),8)
S2.scen.Snr_s=rep(c(1,2,1,1.2,2.2,1.2),8)

# environment slope
S2.env.scen.Sr_w=rep(c(0.005,0.08,0.005,0.005,0.08,0.005),8)
S2.env.scen.D=rep(c(1.1,1.1,1.1,1.1,1.1,1.1),8)
S2.env.scen.Sd=rep(c(0.01,0.1,0.01,0.01,0.1,0.01),8)
S2.env.scen.B=rep(c(0.1,0.06,0.1,0.1,0.06,0.1),8)
S2.env.scen.Snr_w=rep(c(0.005,0.1,0.005,0.005,0.1,0.005),8)
S2.env.scen.Sj_w=rep(c(10,5,10,10,5,10),8)
S2.env.scen.Tj=rep(c(-0.05,-0.05,-0.05,-0.05,-0.05,-0.05),8)
S2.env.scen.Tnr=rep(c(-1.1,-0.5,-1.1,-1.1,-0.5,-1.1),8)
S2.env.scen.Tr=rep(c(-0.005,0.01,-0.005,-0.005,0.01,-0.005),8)
S2.env.scen.Sr_s=rep(c(0.005,0.05,0.005,0.005,0.05,0.005),8)
S2.env.scen.Snr_s=rep(c(0.01,0.05,0.01,0.01,0.05,0.01),8)

# cap on positive density effects (for reproduction only):

S2.bi2.scen.B=rep(c(rep(c(rep(-0.001,2),-0.0012),2),rep(c(rep(-0.001,2),-0.0017),2),rep(c(rep(-0.001,2),-0.0017),4)),2)

# biotic interactions Prey 

S2.bi.scen.Sr_w=c(rep(1.2,6),rep(1.5,6),rep(1.5,12),rep(1.2,6),rep(1.5,6),rep(1.5,12))
S2.bi.scen.D=c(rep(-1.2,6),rep(-1.5,6),rep(-1.5,12),rep(-1.2,6),rep(-1.5,6),rep(-1.5,12))
S2.bi.scen.Sd=c(rep(1.5,6),rep(1.8,6),rep(1.8,12),rep(1.5,6),rep(1.8,6),rep(1.8,12))
S2.bi.scen.B=c(rep(2,6),rep(2.3,6),rep(2.3,12),rep(2,6),rep(2.3,6),rep(2.3,12))
S2.bi.scen.Snr_w=c(rep(8,6),rep(8.3,6),rep(8.3,12),rep(8,6),rep(8.3,6),rep(8.3,12))
S2.bi.scen.Sj_w=c(rep(2.8,6),rep(3.0,6),rep(3.0,12),rep(2.8,6),rep(3.0,6),rep(3.0,12))
S2.bi.scen.Tj=c(rep(3.2,6),rep(3.5,6),rep(3.5,12),rep(3.2,6),rep(3.5,6),rep(3.5,12))
S2.bi.scen.Tnr=c(rep(1.5,6),rep(1.8,6),rep(1.8,12),rep(1.5,6),rep(1.8,6),rep(1.8,12))
S2.bi.scen.Tr=c(rep(1.4,6),rep(1.7,6),rep(1.7,12),rep(1.4,6),rep(1.7,6),rep(1.7,12))
S2.bi.scen.Sr_s=c(rep(2.3,6),rep(2.6,6),rep(2.6,12),rep(2.3,6),rep(2.6,6),rep(2.6,12))
S2.bi.scen.Snr_s=c(rep(3.2,6),rep(3.5,6),rep(3.5,12),rep(3.2,6),rep(3.5,6),rep(3.5,12))

# biotic interactions competitor (predator 2) 

S2.biS2.scen.Sr_w=c(rep(-0.07,18),rep(-0.04,6),rep(-0.07,18),rep(-0.04,6))
S2.biS2.scen.D=rep(0.1,48)
S2.biS2.scen.Sd=rep(0,48)
S2.biS2.scen.B=rep(-0.5,48)
S2.biS2.scen.Snr_w=c(rep(-0.07,18),rep(-0.04,6),rep(-0.07,18),rep(-0.04,6))
S2.biS2.scen.Sj_w=rep(-0.07,48)
S2.biS2.scen.Tj=rep(0,48)
S2.biS2.scen.Tnr=rep(0,48)
S2.biS2.scen.Tr=rep(0,48)
S2.biS2.scen.Sr_s=c(rep(-0.07,18),rep(-0.1,6),rep(-0.07,18),rep(-0.1,6))
S2.biS2.scen.Snr_s=c(rep(-0.07,18),rep(-0.1,6),rep(-0.07,18),rep(-0.1,6))


# dens dens interactions

S2.dens.bi.scen.Sr_w=c(rep(0,6),rep(0.0001,18),rep(0,6),rep(0.0001,18))
S2.dens.bi.scen.D=c(rep(0,6),rep(0.0001,18),rep(0,6),rep(0.0001,18))
S2.dens.bi.scen.Sd=c(rep(0,6),rep(0.0001,18),rep(0,6),rep(0.0001,18))
S2.dens.bi.scen.B=c(rep(0,6),rep(0.0001,18),rep(0,6),rep(0.0001,18))
S2.dens.bi.scen.Snr_w=c(rep(0,6),rep(0.0001,18),rep(0,6),rep(0.0001,18))
S2.dens.bi.scen.Sj_w=c(rep(0,6),rep(0.0001,18),rep(0,6),rep(0.0001,18))
S2.dens.bi.scen.Tj=c(rep(0,6),rep(0.0001,18),rep(0,6),rep(0.0001,18))
S2.dens.bi.scen.Tnr=c(rep(0,6),rep(0.0001,18),rep(0,6),rep(0.0001,18))
S2.dens.bi.scen.Tr=c(rep(0,6),rep(0.0001,18),rep(0,6),rep(0.0001,18))
S2.dens.bi.scen.Sr_s=c(rep(0,6),rep(0.0001,18),rep(0,6),rep(0.0001,18))
S2.dens.bi.scen.Snr_s=c(rep(0,6),rep(0.0001,18),rep(0,6),rep(0.0001,18))

## spatial differences

S2.AB.scen.Sr_w=c(rep(0,24),rep(-0.1,24))
S2.AB.scen.Sd=c(rep(0,24),rep(-0.1,24))
S2.AB.scen.Snr_w=c(rep(0,24),rep(-0.1,24))
S2.AB.scen.Sj_w=c(rep(0,24),rep(-0.1,24))
S2.AB.scen.Sr_s=c(rep(0,24),rep(-0.1,24))
S2.AB.scen.Snr_s=c(rep(0,24),rep(-0.1,24))

# spatial differences in environment slope

S2.AB.env.scen.Sr_w=c(rep(0,24),rep(-0.1,24))
S2.AB.env.scen.D=c(rep(0,24),rep(-0.1,24))
S2.AB.env.scen.Sd=c(rep(0,24),rep(-0.1,24))
S2.AB.env.scen.B=c(rep(0,24),rep(-0.1,24))
S2.AB.env.scen.Snr_w=c(rep(0,24),rep(-0.1,24))
S2.AB.env.scen.Sj_w=c(rep(0,24),rep(-0.1,24))
S2.AB.env.scen.Tj=c(rep(0,24),rep(-0.1,24))
S2.AB.env.scen.Tnr=c(rep(0,24),rep(-0.1,24))
S2.AB.env.scen.Tr=c(rep(0,24),rep(-0.1,24))
S2.AB.env.scen.Sr_s=c(rep(0,24),rep(-0.1,24))
S2.AB.env.scen.Snr_s=c(rep(0,24),rep(-0.1,24))

####################

###  Predator species 2

# mean vital rates
S3.scen.Sr_w=rep(c(2,2,2,1.5,1.5,1.5),8)
S3.scen.D=rep(0.1,48)
S3.scen.Sd=rep(c(1,0.5,1,1,0.5,1),8)
S3.scen.B=rep(c(-0.1,0.01,-0.1,-0.1,0.01,-0.1),8)
S3.scen.Snr_w=rep(c(2.1,2,2.1,2,1.8,2),8)
S3.scen.Sj_w=rep(c(0.04,0.04,0.04,0.04,0.04,0.04),8)
S3.scen.Tj=rep(c(-1,0.001,-1,-1,0.001,-1),8)
S3.scen.Tnr=rep(c(0.7,0.7,0.7,0.7,0.7,0.7),8)
S3.scen.Tr=rep(c(1,1,1,1,1,1),8)
S3.scen.Sr_s=rep(c(2.5,2.7,2.5,3,3.2,3),8)
S3.scen.Snr_s=rep(c(0.5,0.7,0.5,0.6,1,0.6),8)

# environment slope
S3.env.scen.Sr_w=rep(c(0.009,0.15,0.009,0.009,0.15,0.009),8)
S3.env.scen.D=rep(c(0.1,0.1,0.1,0.1,0.1,0.1),8)
S3.env.scen.Sd=rep(c(0.07,0.7,0.07,0.07,0.7,0.07),8)
S3.env.scen.B=rep(c(0.07,0.05,0.07,0.07,0.05,0.07),8)
S3.env.scen.Snr_w=rep(c(0.01,0.1,0.01,0.01,0.1,0.01),8)
S3.env.scen.Sj_w=rep(c(8,3,8,8,3,8),8)
S3.env.scen.Tj=rep(c(-0.05,-0.05,-0.05,-0.05,-0.05,-0.05),8)
S3.env.scen.Tnr=rep(c(-0.8,-0.5,-0.8,-0.8,-0.5,-0.8),8)
S3.env.scen.Tr=rep(c(-0.005,0.01,-0.005,-0.005,0.01,-0.005),8)
S3.env.scen.Sr_s=rep(c(0.005,0.05,0.005,0.005,0.05,0.005),8)
S3.env.scen.Snr_s=rep(c(0.01,0.05,0.01,0.01,0.05,0.01),8)

# cap on positive density effects (for reproduction only):

S3.bi2.scen.B=rep(rep(c(rep(-0.001,2),-0.0005),8),2)

# biotic interactions Prey (a bit lower than fro predator 1 - less specialized)

S3.bi.scen.Sr_w=rep(0.7,48)
S3.bi.scen.D=rep(-1.5,48)
S3.bi.scen.Sd=rep(1.1,48)
S3.bi.scen.B=rep(2.2,48)
S3.bi.scen.Snr_w=rep(6,48)
S3.bi.scen.Sj_w=rep(2.5,48)
S3.bi.scen.Tj=rep(1.2,48)
S3.bi.scen.Tnr=rep(0.9,48)
S3.bi.scen.Tr=rep(0.9,48)
S3.bi.scen.Sr_s=rep(1.8,48)
S3.bi.scen.Snr_s=rep(2,48)

# biotic interactions competitor (predator 1) - more negatively affected than predator 1 is by predator 2

S3.biS2.scen.Sr_w=c(rep(-0.6,18),rep(-0.3,6),rep(-0.6,18),rep(-0.3,6))
S3.biS2.scen.D=rep(0.3,48)
S3.biS2.scen.Sd=rep(0.2,48)
S3.biS2.scen.B=rep(-0.7,48)
S3.biS2.scen.Snr_w=c(rep(-0.7,18),rep(-0.4,6),rep(-0.7,18),rep(-0.4,6))
S3.biS2.scen.Sj_w=rep(-0.7,48)
S3.biS2.scen.Tj=rep(-0.1,48)
S3.biS2.scen.Tnr=rep(0,48)
S3.biS2.scen.Tr=rep(0,48)
S3.biS2.scen.Sr_s=c(rep(-0.6,18),rep(-0.9,6),rep(-0.6,18),rep(-0.9,6))
S3.biS2.scen.Snr_s=c(rep(-0.6,18),rep(-0.9,6),rep(-0.6,18),rep(-0.9,6))


## spatial differences

S3.AB.scen.Sr_w=c(rep(0,24),rep(0.05,24))
S3.AB.scen.Sd=c(rep(0,24),rep(0.05,24))
S3.AB.scen.Snr_w=c(rep(0,24),rep(0.05,24))
S3.AB.scen.Sj_w=c(rep(0,24),rep(0.05,24))
S3.AB.scen.Sr_s=c(rep(0,24),rep(0.05,24))
S3.AB.scen.Snr_s=c(rep(0,24),rep(0.05,24))

# dens dens interactions

S3.dens.bi.scen.Sr_w=rep(c(0.0001,0,0.001,0.0001,0,0.001),8)
S3.dens.bi.scen.D=rep(c(0.0001,0,0.001,0.0001,0,0.001),8)
S3.dens.bi.scen.Sd=rep(c(0.0001,0,0.001,0.0001,0,0.001),8)
S3.dens.bi.scen.B=rep(c(0.0001,0,0.001,0.0001,0,0.001),8)
S3.dens.bi.scen.Snr_w=rep(c(0.0001,0,0.001,0.0001,0,0.001),8)
S3.dens.bi.scen.Sj_w=rep(c(0.0001,0,0.001,0.0001,0,0.001),8)
S3.dens.bi.scen.Tj=rep(c(0.0001,0,0.001,0.0001,0,0.001),8)
S3.dens.bi.scen.Tnr=rep(c(0.0001,0,0.001,0.0001,0,0.001),8)
S3.dens.bi.scen.Tr=rep(c(0.0001,0,0.001,0.0001,0,0.001),8)
S3.dens.bi.scen.Sr_s=rep(c(0.0001,0,0.001,0.0001,0,0.001),8)
S3.dens.bi.scen.Snr_s=rep(c(0.0001,0,0.001,0.0001,0,0.001),8)


# spatial differences in environment slope

S3.AB.env.scen.Sr_w=c(rep(0,24),rep(0.2,24))
S3.AB.env.scen.D=c(rep(0,24),rep(0.2,24))
S3.AB.env.scen.Sd=c(rep(0,24),rep(0.2,24))
S3.AB.env.scen.B=c(rep(0,24),rep(0.2,24))
S3.AB.env.scen.Snr_w=c(rep(0,24),rep(0.2,24))
S3.AB.env.scen.Sj_w=c(rep(0,24),rep(0.2,24))
S3.AB.env.scen.Tj=c(rep(0,24),rep(0.2,24))
S3.AB.env.scen.Tnr=c(rep(0,24),rep(0.2,24))
S3.AB.env.scen.Tr=c(rep(0,24),rep(0.2,24))
S3.AB.env.scen.Sr_s=c(rep(0,24),rep(0.2,24))
S3.AB.env.scen.Snr_s=c(rep(0,24),rep(0.2,24))

```

## 4. Simulating baseline environmental fluctuations

Once the demographic rate models are created, and a range of parameters for these models are defined, which allow to parameterize different communities, we can simulate baseline metapopulaiton dynamics of the three interacting species. These baseline simulations allow us to set up the climate simulations later on and are meant as a control to ensure everything runs smoothly.

The raw R script in this repository (_seasonal_mltispecies_demo.R_), shows how to run one simulation (fixing the environmental variation with set_seed() to obtain relevant plots of abundances in each community. Here, we skip this step and go right to the core simulations.

We run each simulation for 1050 years, retaining the last 50 years (using 1000 years as burn-in for the simulations).

```{r}
## Base simulations

years=1050 # years of simulations
```

We create environmental variation, by running 100 simulations, and in each simulation, creating an environmental variable for the length of the years simulated. The values of this variable are sampled from a multivariate normal distribution with $\mu$ = -1 for the unproductive season (winter) and $\mu$ = 1 for the productive season (summer) and a covariance matrix between the two seasons as follows: 

```{r}
simul=100

# simulate winter and summer environment from multivariate normal
cov=matrix(c(0.5,0,0,0.1),2,2)

```

So, initially, we simulate a higher variance in the unproductive season, but no covariance.

We create empty objects to fill with simulation results:

```{r}
envir.sim=NULL  # save the environmental values
base.sim=NULL   # save total abundances (mean and CV across years)
base.sim.stage=NULL # save stage-specific abundances (mean and CV across years)
remove.sim=NULL # if extinction occurred in the 1000 burn-in years, remove these simulations (which we need to get initial values for climate-change simulations)

init.dens=NULL # save starting densities for climate change simulations 

cor.ab=array(NA,c(length(scen),simul,3,3)) # correlattion between mean abundances of the three species
```

And then, proceed to simulate metapopulation dynamics for the three interacting species, starting with the summer season. 

```{r, eval = FALSE}
library(popbio)
library(MASS)
library(boot)
library(ggplot2)
library(patchwork)

CV <- function(x){
  (sd(x,na.rm = T)/mean(x,na.rm = T))*100
}

# Start with the communities (scen)
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
  
  # Once the parameters for a community are set, go through each of 100 simulations
  for(x in 1:simul){
    
    # Define environmental variation
    env=mvrnorm(years, mu=c(-1,1), cov)
    
    # Environmental vectors for winter and summer
    envW=env[,1]
    envS=env[,2]
    
    # Set initial densities
    
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
   
    # Initial data frame 
   temp=data.frame(dens=c(densA,densB,S2_densA,S2_densB,S3_densA,S3_densB),
                      species=rep(c("Prey","Predator1", "Predator2"),each=6),
                      site=rep(rep(c("A","B"),each=3),3),
                      season="summer",
                      stage=rep(c("J","N","R"),6),
                      year=1)
    
   # Place holder to calculate life history traits 
   
    matsS1=array(0,c(3,3,50))
    matsS2=array(0,c(3,3,50))
    matsS3=array(0,c(3,3,50))
   
    # Go through years: 
    for(y in 1:years){
      
      #summer
      
      #### LOCAL DEMOGRAPHY
      #Prey
      summerA=matrix(c(0,0,0,
                       0,Snr_s("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA),0,
                       Sr_s("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*B("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA),0,Sr_s("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)),3,3,byrow = F)
      
      summerB=matrix(c(0,0,0,
                       0,Snr_s("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB),0,
                       Sr_s("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*B("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB),0,Sr_s("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)),3,3,byrow = F)
      
      
      
      Bs[1:3,1:3]=summerA
      Bs[4:6,4:6]=summerB
      
      ### DISPERSAL
      
      Ms[3,3]=1-D("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)
      Ms[4,4]=1-D("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)
      Ms[6,3]=D("A",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)*Ds("A",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)
      Ms[5,4]=D("B",envS[y],sum(densB[1:3]),S2_dens.biB,S3_dens.biB)*Ds("B",envS[y],sum(densA[1:3]),S2_dens.biA,S3_dens.biA)
      
      
      # VEC PERMUTATION APPROACH TO LINK LOCAL DEMOGRAPHY TO DISPERSAL AND GET POPULATION VECTOR ACROSS TWO SITES AT THE BEGINNING OF WINDER SEASON
      
      densW=t(Ps)%*%Ms%*%Ps%*%Bs%*%densMet
      
      # SAVE OUTPUT FROM LAST 50 YEARS
      if(y>1000){
        
        temp=rbind(temp,data.frame(dens=c(round(densW[1:3],2),round(densW[4:6],2)),
                                   species="Prey",
                                   site=rep(c("A","B"),each=3),
                                   season=rep("winter",6),stage=rep(c("J","N","R"),2),year=rep(y,6)))
      }
      
      # UPDATE DENSITIES 
      densA=densW[1:3]
      densB=densW[4:6]
      densMet=c(densA,densB)
      
      # SAME APPROACH AS FOR S1 
      
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
      
      # ONCE WE GO THROUGH ALL SPECIES, WE UPDATE THE DENSITIES USED AS PREDICTORS DESCRIBING BIOTIC INTERACTIONS 
      
      ### Update bi densities
      
      dens.biA=sum(densA)
      dens.biB=sum(densB)
      
      S2_dens.biA=sum(S2_densA)
      S2_dens.biB=sum(S2_densB)
      
      S3_dens.biA=sum(S3_densA)
      S3_dens.biB=sum(S3_densB)
      
      ############ REPEAT FOR WINTER 
      
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
    
    # GET SUMMARY METRICS
    
    # environment 
    
    envir.sim=rbind(envir.sim,data.frame(envW=env[1000:years,1],
                                         envS=env[1000:years,2],
                                         sim=x,
                                         scen=k))
    temp=temp[!temp$year%in%1,]
    temp$year=temp$year-1000
    
    
    # TOTAL DENSITY PER YEAR AND SPECIES
    densS1=aggregate(dens~year,data = temp[temp$species=="Prey"&temp$season%in%"winter",],sum)
    densS2=aggregate(dens~year,data = temp[temp$species=="Predator1"&temp$season%in%"winter",],sum)
    densS3=aggregate(dens~year,data = temp[temp$species=="Predator2"&temp$season%in%"winter",],sum)
    
    if(length(unique(densS1$dens))<2|length(unique(densS2$dens))<2|length(unique(densS3$dens))<2){
      
      print("Issue with Abundance")
      
      remove.sim=rbind(remove.sim,data.frame(sim=x,
                                             scen=k))
      next
      
      }else{
      
      # EXTINCTION PROBABILITY 
      ext1=ext2=ext3=NA
      
      if(any(densS1$dens<1)) ext1 <- 1 else ext1 <- 0
      
      if(any(densS2$dens<1)) ext2 <- 1 else ext2 <- 0
      
      if(any(densS3$dens<1)) ext3 <- 1 else ext3 <- 0
      
      # GENERATION TIME AND NET REPRODUCTIVE RATE 
      
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

# Quick test that everything is correct (there should be few extinctions)

base.sim[base.sim$ext%in%1,] # check for extinction

hist(base.sim$Dens.mu[base.sim$species%in%"Prey"])
hist(base.sim$Dens.mu[base.sim$species%in%"Predator1"])
hist(base.sim$Dens.mu[base.sim$species%in%"Predator2"])

# Save output if necessary

write.csv(envir.sim,"envir.sim.csv",row.names = F)
write.csv(base.sim,"base.sim.csv",row.names = F)
write.csv(base.sim.stage,"base.sim.stage.csv",row.names = F)
write.csv(init.dens,"init.dens.csv",row.names = F)
write.csv(remove.sim,"remove.sim.csv",row.names = F)

save(cor.ab,file = "corr.sim.Rdata" )


```


## 5. Simulating climate change

Once the baseline model is run, and we plot the results, ensuring that extinctions are rates, and abundance distributions are as expected, we proceed to run the climate change simulation, which are then used to obtain our main results. 

The climate-cahnge simulations perturb the seasonal mean environments: 

**Baseline** - $\mu$ = -1 and $\mu$ = 1 for unproductive and productive season, respectively. 

**-W** - $\mu$ = -1.5 in unproductive season

**-S** - $\mu$ = 0.5 in productive season

**-W-S** - $\mu$ decreases in both seasons

**+W-S** - $\mu$ = -0.5 and $\mu$ = 0.5 for unproductive and productive season, respectively. 

**-W+S** - $\mu$ = -1.5 and $\mu$ = 1.5 for unproductive and productive season, respectively. 


Variance in environmental states is also perturbed by **increasing variance in both seasons by 0.1**

Covariance in environmental states is also perturbed by **setting the covariance between seasons to 0.2**

The script to run the simulations mimics the one for section 4, except that simulations are only run for 50 years, and the environment is changed based on the above scenarios.

```{r, eval=FALSE}

# Load if necessary 

init.dens=read.csv("init.dens.csv")
remove.sim=read.csv("remove.sim.csv")

years=50 # years of simulations

# simulate winter and summer environment from multivariate normal

cov=c(0,0.2)
var=c(0,0.1)

means=rbind(cbind(-1,1),cbind(-1.5,1),cbind(-1,0.5),cbind(-1.5,0.5),cbind(-0.5,0.5),cbind(-1.5,1.5))

clim.change.scen=c("baseline","-W","-S","-W-S","+W-S","-W+S")

var.change=c("no","yes")

cov.change=c("no","yes")

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

```

You can save the output of these simulations as they are needed for plotting:

```{r,eval=FALSE}
nrow(clim.change.sim[clim.change.sim$ext%in%1,]) # not much extinction

# save output
write.csv(envir.clim.change.sim,"envir.clim.change.sim.csv",row.names = F)

write.csv(clim.change.sim,"clim.change.sim.csv",row.names = F)
write.csv(clim.change.sim.stage,"clim.change.sim.stage.csv",row.names = F)
write.csv(remove.sim.clim.change,"remove.sim.csv",row.names = F)


```




