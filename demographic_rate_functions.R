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
