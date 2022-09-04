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
