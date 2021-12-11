# Data Kemiskinan #
library('readxl')
library(writexl)
dataset <- read_excel("E:/Dokumenku/2019_Kemiskinan_Var.xlsx")


attach(dataset)
x0 <- array(1, dim = length(y))                  ##Intercecpt##
X  <- cbind(x0,x1,x2,x3,x4,x5,x6,x7,x8,x9)
k  <- 9
n  <- length(y)

library('sp')
Y  <- scale(y)
X1 <- scale(x1)
X2 <- scale(x2)
X3 <- scale(x3)
X4 <- scale(x4)
X5 <- scale(x5)
X6 <- scale(x6)
X7 <- scale(x7)
X8 <- scale(x8)
X9 <- scale(x9)

#regresi parametrik
reg=lm(y~x1+x2+x3+x4+x5+x6+x7+x8+x9)
summary(reg)



#longlat
library('rgdal')
Koordinat    <- cbind(b,l)
Z<-as.matrix(Koordinat)

#data
library('GWmodel')
dtgab <- data.frame(y,x1,x2,x3,x4,x5,x6,x7,x8,x9)
m     <- matrix(1,nrow(dataset),1)
Xm    <- cbind(m,x1,x2,x3,x4,x5,x6,x7,x8,x9)
X     <- cbind(x1,x2,x3,x4,x5,x6,x7,x8,x9)
Y     <- cbind(Y)
MatriksKoordinat     <- as.matrix(Cordinat)
d.mat <- gw.dist(dp.locat=MatriksKoordinat,p=2, focus=0,theta=0,longlat=T)
dtnew<-SpatialPointsDataFrame(Z,as.data.frame(dtgab))

Matriksdij <- read_excel("E:/Dokumenku/dij.xlsx")
Matriksdij<-as.matrix(Matriksdij)
datapakailonglat <- SpatialPointsDataFrame(Matriksdij,data=dtgab, match.ID=TRUE)


###############################################################
####----Model ROBUST GEOGRAPHICALLY WEIGHTED REGRESSION----####
###############################################################
library('SparseM')
library('quantreg')

##Mendefinisikan Bandwidth
rgwr.b=function (bw, X, Y, kernel, adaptive = F, dp.locat,
                 p = 2, theta = 0, longlat = F, dMat, verbose = T)
{
  dp.n = length(dp.locat[, 1])
  dim.dMat = dim(dMat)
  ACV = numeric(dp.n)
  m = matrix(1,length(dp.locat[,1]),1)
  for (i in 1:dp.n) {
    dist.vi = gw.dist(dp.locat = dp.locat, focus = i, p = p, theta = theta, longlat = longlat)
    W.i = gw.weight(dist.vi, bw, kernel, adaptive)
    W.i[i] = 0
    W.i = c(W.i)
    rq.f = function(X,Y,W){
      Xwi = X*W
      Ywi = Y*W
      Xwicsr = as.matrix.csr(cbind(m,Xwi))
      m.rq = rq.fit.sfn(Xwicsr,Ywi,tau=0.5)
      beta = cbind(m.rq$coef)
      return(beta)
    }
    hasil = try(rq.f(X, Y, W.i))
    xnew = cbind(m,X)
    yhat.noi = xnew[i, ] %*% hasil
    ACV[i] = Y[i] - yhat.noi
  }
  ACV.score = sum(abs(ACV))
  if (verbose) {
    if (adaptive)
      cat("Adaptive bandwidth:", bw, "ACV score:", ACV.score, "\n")
    else cat("Fixed bandwidth:", bw, "ACV score:", ACV.score, "\n")
  }
  ACV.score
}
rgwr.bo=function (formula, data, approach = "ACV", kernel,
                  adaptive = F, p = 2, theta = 0, longlat = F, dMat)
{
  dp.locat = coordinates(data)
  data = as(data, "data.frame")
  data = as.matrix(as.data.frame(data))
  mf = match.call(expand.dots = FALSE)
  m = match(c("formula", "data"), names(mf), 0L)
  mf = mf[c(1L, m)]
  mf$drop.unused.levels =TRUE
  mf[[1L]] = as.name("model.frame")
  mf = eval(mf, parent.frame())
  mt = attr(mf, "terms")
  y = model.extract(mf, "response")
  x = model.matrix(mt, mf)
  dp.n = nrow(data)
  if (dp.n > 1500) {
    cat("Take a cup of tea and have a break, it will take a few minutes.\n")
    cat("-----A kind suggestion from GWmodel development group\n")
  }
  if (missing(dMat)) {
    DM.given = F
    if (dp.n + dp.n <= 10000) {
      dMat = gw.dist(dp.locat = dp.locat, rp.locat = dp.locat,
                     p = p, theta = theta, longlat = longlat)
      DM.given = T
    }}
  else {
    DM.given = T
    dim.dMat = dim(dMat)
    if (dim.dMat[1] != dp.n || dim.dMat[2] != dp.n)
      stop("Dimensions of dMat are not correct")
  }
  if (adaptive) {
    upper = dp.n
    lower = 20
  }
  else {
    if (DM.given) {
      upper = range(dMat)[2]
      lower = upper/5000
    }
    else {
      dMat = NULL
      if (p == 2) {
        b.box = bbox(dp.locat)
        upper = sqrt((b.box[1, 2] - b.box[1, 1])^2 +
                       (b.box[2, 2] - b.box[2, 1])^2)
        lower = upper/5000
      }
      else {
        upper = 0
        for (i in 1:dp.n) {
          dist.vi = gw.dist(dp.locat = dp.locat, focus = i, p = p, theta = theta, longlat = longlat)
          upper = max(upper, range(dist.vi)[2])
        }
        lower = upper/5000
      } } }
  bw = NA
  bw = gold(rgwr.b, lower, upper, adapt.bw = adaptive, x, y, kernel, adaptive, dp.locat, p, theta, longlat,
            dMat)
  bw
}
##Bandwidth RGWR##
bwo.AI <- rgwr.bo(y~x1+x2+x3+x4+x5+x6+x7+x8+x9, data=dtnew, kernel="boxcar", dMat=Matriksdij, adaptive=T)
bwo.FI <- rgwr.bo(y~x1+x2+x3+x4+x5+x6+x7+x8+x9, data=dtnew, kernel="boxcar", dMat=Matriksdij, adaptive=F)
bwo.AG <- rgwr.bo(y~x1+x2+x3+x4+x5+x6+x7+x8+x9, data=dtnew, kernel="gaussian", dMat=Matriksdij, adaptive=T)
bwo.FG <- rgwr.bo(y~x1+x2+x3+x4+x5+x6+x7+x8+x9, data=dtnew, kernel="gaussian", dMat=Matriksdij, adaptive=F)
bwo.AE <- rgwr.bo(y~x1+x2+x3+x4+x5+x6+x7+x8+x9, data=dtnew, kernel="exponential", dMat=Matriksdij, adaptive=T)
bwo.FE <- rgwr.bo(y~x1+x2+x3+x4+x5+x6+x7+x8+x9, data=dtnew, kernel="exponential", dMat=Matriksdij, adaptive=F)
bwo.AB <- rgwr.bo(y~x1+x2+x3+x4+x5+x6+x7+x8+x9, data=dtnew, kernel="bisquare", dMat=Matriksdij, adaptive=T)
bwo.FB <- rgwr.bo(y~x1+x2+x3+x4+x5+x6+x7+x8+x9, data=dtnew, kernel="bisquare", dMat=Matriksdij, adaptive=F)
bwo.AT <- rgwr.bo(y~x1+x2+x3+x4+x5+x6+x7+x8+x9, data=dtnew, kernel="tricube", dMat=Matriksdij, adaptive=T)
bwo.FT <- rgwr.bo(y~x1+x2+x3+x4+x5+x6+x7+x8+x9, data=dtnew, kernel="tricube", dMat=Matriksdij, adaptive=F)


###RGWR ADAPTIVE BISQUARE###
##Bobot Adaptive Bisquare
W.RGWR.AB <- gw.weight(Matriksdij,bw=bwo.AB,kernel="bisquare",adaptive=TRUE)
W.RGWR.AB<-as.data.frame(W.RGWR.AB)
library(writexl)
write_xlsx(W.RGWR, "E:/Dokumenku/BobotRGWR.xlsx")

##model RGWR Lokal FI
SDFRGWRLokal_FI <- matrix(0,nrow(X),(ncol(X)+1))
for (i in 1:nrow(X)){
  dmat                   <- Matriksdij[,i]
  W_RGWR_FI              <- gw.weight(dmat,bw=bwo.FI,kernel="boxcar",adaptive=FALSE)
  rgwr_invA             <- rq(y~x1+x2+x3+x4+x5+x6+x7+x8+x9,data=dtnew,tau=0.5,weight=W_RGWR_FI)
  SDFRGWRLokal_FI[i,]    <- rgwr_invA[[1]]
}
RGWRLokal_FI<-as.data.frame(SDFRGWRLokal_FI)
write_xlsx(RGWRLokal_FI, "E:/Dokumenku/SDF/RGWRLokal_FI.xlsx")
###################

##model RGWR Lokal AI
SDFRGWRLokal_AI <- matrix(0,nrow(X),(ncol(X)+1))
for (i in 1:nrow(X)){
  dmat                   <- Matriksdij[,i]
  W_RGWR_AI              <- gw.weight(dmat,bw=bwo.AI,kernel="boxcar",adaptive=TRUE)
  rgwr_invB             <- rq(y~x1+x2+x3+x4+x5+x6+x7+x8+x9,data=dtnew,tau=0.5,weight=W_RGWR_AI)
  SDFRGWRLokal_AI[i,]    <- rgwr_invB[[1]]
}
RGWRLokal_AI<-as.data.frame(SDFRGWRLokal_AI)
write_xlsx(RGWRLokal_AI, "E:/Dokumenku/SDF/RGWRLokal_AI.xlsx")
###################

##model RGWR Lokal FG
SDFRGWRLokal_FG <- matrix(0,nrow(X),(ncol(X)+1))
for (i in 1:nrow(X)){
  dmat                   <- Matriksdij[,i]
  W_RGWR_FG              <- gw.weight(dmat,bw=bwo.FG,kernel="gaussian",adaptive=FALSE)
  rgwr_gaussA             <- rq(y~x1+x2+x3+x4+x5+x6+x7+x8+x9,data=dtnew,tau=0.5,weight=W_RGWR_FG)
  SDFRGWRLokal_FG[i,]    <- rgwr_gaussA[[1]]
}
RGWRLokal_FG<-as.data.frame(SDFRGWRLokal_FG)
write_xlsx(RGWRLokal_FG, "E:/Dokumenku/SDF/RGWRLokal_FG.xlsx")
###################

##model RGWR Lokal AG
SDFRGWRLokal_AG <- matrix(0,nrow(X),(ncol(X)+1))
for (i in 1:nrow(X)){
  dmat                   <- Matriksdij[,i]
  W_RGWR_AG              <- gw.weight(dmat,bw=bwo.AG,kernel="gaussian",adaptive=TRUE)
  rgwr_gaussB             <- rq(y~x1+x2+x3+x4+x5+x6+x7+x8+x9,data=dtnew,tau=0.5,weight=W_RGWR_AG)
  SDFRGWRLokal_AG[i,]    <- rgwr_gaussB[[1]]
}
RGWRLokal_AG<-as.data.frame(SDFRGWRLokal_AG)
write_xlsx(RGWRLokal_AG, "E:/Dokumenku/SDF/RGWRLokal_AG.xlsx")
###################

##model RGWR Lokal FE
SDFRGWRLokal_FE <- matrix(0,nrow(X),(ncol(X)+1))
for (i in 1:nrow(X)){
  dmat                   <- Matriksdij[,i]
  W_RGWR_FE              <- gw.weight(dmat,bw=bwo.FE,kernel="exponential",adaptive=FALSE)
  rgwr_expA             <- rq(y~x1+x2+x3+x4+x5+x6+x7+x8+x9,data=dtnew,tau=0.5,weight=W_RGWR_FE)
  SDFRGWRLokal_FE[i,]    <- rgwr_expA[[1]]
}
RGWRLokal_FE<-as.data.frame(SDFRGWRLokal_FE)
write_xlsx(RGWRLokal_FE, "E:/Dokumenku/SDF/RGWRLokal_FE.xlsx")
###################

##model RGWR Lokal AE
SDFRGWRLokal_AE <- matrix(0,nrow(X),(ncol(X)+1))
for (i in 1:nrow(X)){
  dmat                   <- Matriksdij[,i]
  W_RGWR_AE              <- gw.weight(dmat,bw=bwo.AE,kernel="exponential",adaptive=TRUE)
  rgwr_expB             <- rq(y~x1+x2+x3+x4+x5+x6+x7+x8+x9,data=dtnew,tau=0.5,weight=W_RGWR_AE)
  SDFRGWRLokal_AE[i,]    <- rgwr_expB[[1]]
}
RGWRLokal_AE<-as.data.frame(SDFRGWRLokal_AE)
write_xlsx(RGWRLokal_AE, "E:/Dokumenku/SDF/RGWRLokal_AE.xlsx")
###################


##model RGWR Lokal FB
SDFRGWRLokal_FB <- matrix(0,nrow(X),(ncol(X)+1))
for (i in 1:nrow(X)){
  dmat                   <- Matriksdij[,i]
  W_RGWR_FB              <- gw.weight(dmat,bw=bwo.FB,kernel="bisquare",adaptive=FALSE)
  rgwr_bisqA             <- rq(y~x1+x2+x3+x4+x5+x6+x7+x8+x9,data=dtnew,tau=0.5,weight=W_RGWR_FB)
  SDFRGWRLokal_FB[i,]    <- rgwr_bisqA[[1]]
}
RGWRLokal_FB<-as.data.frame(SDFRGWRLokal_FB)
write_xlsx(RGWRLokal_FB, "E:/Dokumenku/SDF/RGWRLokal_FB.xlsx")
###################

##model RGWR Lokal AB
SDFRGWRLokal_AB <- matrix(0,nrow(X),(ncol(X)+1))
for (i in 1:nrow(X)){
  dmat                   <- Matriksdij[,i]
  W_RGWR_AB              <- gw.weight(dmat,bw=bwo.AB,kernel="bisquare",adaptive=TRUE)
  rgwr_bisqB             <- rq(y~x1+x2+x3+x4+x5+x6+x7+x8+x9,data=dtnew,tau=0.5,weight=W_RGWR_AB)
  SDFRGWRLokal_AB[i,]    <- rgwr_bisqB[[1]]
}
RGWRLokal_AB<-as.data.frame(SDFRGWRLokal_AB)
write_xlsx(RGWRLokal_AB, "E:/Dokumenku/SDF/RGWRLokal_AB.xlsx")
###################

##model RGWR Lokal FT
SDFRGWRLokal_FT <- matrix(0,nrow(X),(ncol(X)+1))
for (i in 1:nrow(X)){
  dmat                   <- Matriksdij[,i]
  W_RGWR_FT              <- gw.weight(dmat,bw=bwo.FT,kernel="tricube",adaptive=FALSE)
  rgwr_tricA             <- rq(y~x1+x2+x3+x4+x5+x6+x7+x8+x9,data=dtnew,tau=0.5,weight=W_RGWR_FT)
  SDFRGWRLokal_FT[i,]    <- rgwr_tricA[[1]]
}
RGWRLokal_FT<-as.data.frame(SDFRGWRLokal_FT)
write_xlsx(RGWRLokal_FT, "E:/Dokumenku/SDF/RGWRLokal_FT.xlsx")
###################

##model RGWR Lokal AT
SDFRGWRLokal_AT <- matrix(0,nrow(X),(ncol(X)+1))
for (i in 1:nrow(X)){
  dmat                   <- Matriksdij[,i]
  W_RGWR_AT              <- gw.weight(dmat,bw=bwo.AT,kernel="tricube",adaptive=TRUE)
  rgwr_tricB             <- rq(y~x1+x2+x3+x4+x5+x6+x7+x8+x9,data=dtnew,tau=0.5,weight=W_RGWR_AT)
  SDFRGWRLokal_AT[i,]    <- rgwr_tricB[[1]]
}
RGWRLokal_AT<-as.data.frame(SDFRGWRLokal_AT)
write_xlsx(RGWRLokal_AT, "E:/Dokumenku/SDF/RGWRLokal_AT.xlsx")
###################

#menghitung error RGWR
error.rgwr <- Y-yhat