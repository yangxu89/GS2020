dir <- "/home/xuyang/Data/phe"
setwd(dir)

library(ggplot2)
library(ggthemes)
library(ggpubr)
library(reshape2)
library(doParallel)
library(foreach)
library(parallel)

source(file = "cv.R")
source(file = "cv_fast.R")
source(file = "predict.R")

phe <- read.csv(file = "Input/RIL-Phenotypes.csv", header = T)
snp <- read.csv(file = "Input/RIL-Genotypes.csv", header = T)
exp <- read.csv(file = "Input/RIL-Expressions.csv", header = T)
met <- read.csv(file = "Input/RIL-Metabolites.csv", header = T)
n <- nrow(phe)
gn <- phe$gn
tp <- phe$tp
yd <- phe$yd
kg <- phe$kgw


phe1 <- cbind(yd, tp, gn, kg)
phe0 <- scale(phe1)
phef <- read.csv(file = "Input/IMF2-Phenotypes.csv", header = T)

snp0 <- t(snp[, -1])
exp0 <- scale(t(exp[, -1]))
met0 <- scale(t(met[, -1]))


phe.r11 <- phe0[match(phef[, 2], phe[, 1]), ]
phe.r22 <- phe0[match(phef[, 3], phe[, 1]), ]

snp.r11 <- snp0[match(phef[, 2], rownames(snp0)), ]
snp.r22 <- snp0[match(phef[, 3], rownames(snp0)), ]

exp.r11 <- exp0[match(phef[, 2], rownames(exp0)), ]
exp.r22 <- exp0[match(phef[, 3], rownames(exp0)), ]

met.r11 <- met0[match(phef[, 2], rownames(met0)), ]
met.r22 <- met0[match(phef[, 3], rownames(met0)), ]


########### IMF2 gene phe
snpa <- 0.5 * (snp.r11 + snp.r22)
snpd <- 0.5 * abs(snp.r11 - snp.r22)

expa <- 0.5 * (exp.r11 + exp.r22)
expd <- 0.5 * abs(exp.r11 - exp.r22)

meta <- 0.5 * (met.r11 + met.r22)
metd <- 0.5 * abs(met.r11 - met.r22)

phea <- 0.5 * (phe.r11 + phe.r22)
phed <- 0.5 * abs(phe.r11 - phe.r22)



yield <- phef$yield
kgw <- phef$kgw
grain <- phef
tiller <- phef$tiller
phe12 <- cbind(yield, tiller, grain, kgw)



k1 <- tcrossprod(snpa)/ncol(snpa)
k2 <- tcrossprod(snpd)/ncol(snpd)
k3 <- tcrossprod(expa)/ncol(expa)
dim(k1)
k4 <- tcrossprod(expd)/ncol(expd)
k5 <- tcrossprod(meta)/ncol(meta)
dim(k1)
k6 <- tcrossprod(metd)/ncol(metd)
k7 <- tcrossprod(phea)/ncol(phea)
k8 <- tcrossprod(phed)/ncol(phed)


kin1 <- list()
kin1[[1]] <- k1

kin2 <- list()
kin2[[1]] <- k3

kin3 <- list()
kin3[[1]] <- k5

kin4 <- list()
kin4[[1]] <- k7

kin5 <- list()
kin5[[1]] <- k1
kin5[[2]] <- k2

kin6 <- list()
kin6[[1]] <- k3
kin6[[2]] <- k4

kin7 <- list()
kin7[[1]] <- k5
kin7[[2]] <- k6

kin8 <- list()
kin8[[1]] <- k7
kin8[[2]] <- k8

#### combination
kin9 <- list()
kin9[[1]] <- k1
kin9[[2]] <- k3

kin10 <- list()
kin10[[1]] <- k1
kin10[[2]] <- k5


kin11 <- list()
kin11[[1]] <- k1
kin11[[2]] <- k7

kin12 <- list()
kin12[[1]] <- k3
kin12[[2]] <- k5

kin13 <- list()
kin13[[1]] <- k3
kin13[[2]] <- k7

kin14 <- list()
kin14[[1]] <- k5
kin14[[2]] <- k7

kin15 <- list()
kin15[[1]] <- k1
kin15[[2]] <- k3
kin15[[3]] <- k5

kin16 <- list()
kin16[[1]] <- k1
kin16[[2]] <- k3
kin16[[3]] <- k7

kin17 <- list()
kin17[[1]] <- k1
kin17[[2]] <- k5
kin17[[3]] <- k7

kin18 <- list()
kin18[[1]] <- k3
kin18[[2]] <- k5
kin18[[3]] <- k7

kin19 <- list()
kin19[[1]] <- k1
kin19[[2]] <- k3
kin19[[3]] <- k5
kin19[[4]] <- k7



kin20 <- list()
kin20[[1]] <- k1
kin20[[2]] <- k2
kin20[[3]] <- k3
kin20[[4]] <- k4


kin21 <- list()
kin21[[1]] <- k1
kin21[[2]] <- k2
kin21[[3]] <- k5
kin21[[4]] <- k6

kin22 <- list()
kin22[[1]] <- k1
kin22[[2]] <- k2
kin22[[3]] <- k7
kin22[[4]] <- k8

kin23 <- list()
kin23[[1]] <- k3
kin23[[2]] <- k4
kin23[[3]] <- k5
kin23[[4]] <- k6

kin24 <- list()
kin24[[1]] <- k3
kin24[[2]] <- k4
kin24[[3]] <- k7
kin24[[4]] <- k8


kin25 <- list()
kin25[[1]] <- k5
kin25[[2]] <- k6
kin25[[3]] <- k7
kin25[[4]] <- k8


kin26 <- list()
kin26[[1]] <- k1
kin26[[2]] <- k2
kin26[[3]] <- k3
kin26[[4]] <- k4
kin26[[5]] <- k5
kin26[[6]] <- k6

kin27 <- list()
kin27[[1]] <- k1
kin27[[2]] <- k2
kin27[[3]] <- k3
kin27[[4]] <- k4
kin27[[5]] <- k7
kin27[[6]] <- k8

kin28 <- list()
kin28[[1]] <- k1
kin28[[2]] <- k2
kin28[[3]] <- k5
kin28[[4]] <- k6
kin28[[5]] <- k7
kin28[[6]] <- k8


kin29 <- list()
kin29[[1]] <- k3
kin29[[2]] <- k4
kin29[[3]] <- k5
kin29[[4]] <- k6
kin29[[5]] <- k7
kin29[[6]] <- k8

kin30 <- list()
kin30[[1]] <- k1
kin30[[2]] <- k2
kin30[[3]] <- k3
kin30[[4]] <- k4
kin30[[5]] <- k5
kin30[[6]] <- k6
kin30[[7]] <- k7
kin30[[8]] <- k8

date()
res1 <- list()
resfold <- list()

for (jj in 1:50) {
    cat(jj)
    for (ii in 1:30) {
        cat(ii)
        kk <- get(paste("kin", ii, sep = ""))
        r <- NULL
        for (i in 1:4) {
            y <- as.matrix(phe12[, i])
            # n<-length(y) x<-matrix(1,n,1)
            pred <- cv(y = y, kk = kk, nfold = 10, seed = jj)
            r <- rbind(r, pred)
        }
        res1[[ii]] <- r
    }
    resfold[[jj]] <- do.call(cbind, res1)
}
save(x = resfold, file = " resfold.Rdata")


resfold10 <- resfold
mean0 <- list()
sd0 <- list()
for (i in 1:4) {
    ma <- lapply(resfold, function(x) {
        x[i, ]
    })
    ma1 <- do.call(rbind, ma)
    mean0[[i]] <- apply(ma1, 2, mean)
    sd0[[i]] <- apply(ma1, 2, sd)
}

foldsd <- do.call(rbind, sd0)
foldmean <- Reduce("+", resfold)/50

date()
res2 <- list()
reshat <- list()

for (jj in 1:50) {
    cat(jj)
    for (ii in 1:30) {
        cat(ii)
        kk <- get(paste("kin", ii, sep = ""))
        r <- NULL
        for (i in 1:4) {
            y <- as.matrix(phe12[, i])
            # n<-length(y) x<-matrix(1,n,1)
            pred <- cv_fast(y = y, kk = kk, nfold = 10, seed = jj)
            r <- rbind(r, pred)
        }
        res2[[ii]] <- r
    }
    reshat[[jj]] <- do.call(cbind, res2)
}

sd1 <- list()
for (i in 1:4) {
    ma <- lapply(reshat, function(x) {
        x[i, ]
    })
    ma1 <- do.call(rbind, ma)
    sd1[[i]] <- apply(ma1, 2, sd)
}

hatmean <- Reduce("+", reshat)/50
hatsd <- do.call(rbind, sd1)
######### resnfold
res11 <- list()
for (ii in 1:30) {
    cat(ii)
    kk <- get(paste("kin", ii, sep = ""))
    r <- NULL
    for (i in 1:4) {
        y <- as.matrix(phe12[, i])
        n<-length(y) 
        pred <- cv(y = y, kk = kk, nfold = n)
        r <- rbind(r, pred)
    }
    res11[[ii]] <- r
}
resnfold <- do.call(cbind, res11)

date()
res22 <- list()

for (ii in 1:30) {
    cat(ii)
    kk <- get(paste("kin", ii, sep = ""))
    r <- NULL
    for (i in 1:4) {
        y <- as.matrix(phe12[, i])
        n<-length(y) 
        pred <- cv_fast(y = y, kk = kk, nfold = n)
        r <- rbind(r, pred)
    }
    res22[[ii]] <- r
}
resnhat <- do.call(cbind, res22)

save(resfold, reshat, resnfold, resnhat, file = "predcv.Rdata")




############## add P
king1 <- list()
king1[[1]] <- k1
king1[[2]] <- k2

king2 <- list()
king2[[1]] <- k3
king2[[2]] <- k4

king3 <- list()
king3[[1]] <- k5
king3[[2]] <- k6


king4 <- list()
king4[[1]] <- k1
king4[[2]] <- k2
king4[[3]] <- k3
king4[[4]] <- k4


king5 <- list()
king5[[1]] <- k1
king5[[2]] <- k2
king5[[3]] <- k5
king5[[4]] <- k6

king6 <- list()
king6[[1]] <- k3
king6[[2]] <- k4
king6[[3]] <- k5
king6[[4]] <- k6

king7 <- list()
king7[[1]] <- k1
king7[[2]] <- k2
king7[[3]] <- k3
king7[[4]] <- k4
king7[[5]] <- k5
king7[[6]] <- k6


# phe.r11 phe.r22

phea1 <- 0.5 * (phe.r11 + phe.r22)
phed1 <- 0.5 * abs(phe.r11 - phe.r22)


res22 <- list()
for (ii in 1:7) {
    cat(ii)
    kk <- get(paste("king", ii, sep = ""))
    r1 <- NULL
    r2 <- NULL
    r3 <- NULL
    r4 <- NULL
    r5 <- NULL
    r6 <- NULL
    r7 <- NULL
    r8 <- NULL
    
    for (i in 1:4) {
        cat(i)
        y <- as.matrix(phe12[, i])
        n <- length(y)
        
        x <- cbind(matrix(1, n, 1), phea1, phed1)
        r1[i] <- cv_fast(fix=x, y=y, kk=kk,nfold=n)
        
        
        x <- cbind(matrix(1, n, 1), phea1[, i], phed1[, i])
        r2[i] <- cv_fast(fix=x, y=y, kk=kk, nfold=n)
        
        x <- cbind(matrix(1, n, 1), phe.r11, phe.r22)
        r3[i] <-cv_fast(fix=x, y=y, kk=kk, nfold=n)
        
        x <- cbind(matrix(1, n, 1), phe.r11[, i], phe.r22[, i])
        r4[i] <- cv_fast(fix=x, y=y, kk=kk, nfold=n)
        
        x <- cbind(matrix(1, n, 1), phea1)
        r5[i] <- cv_fast(fix=x, y=y, kk=kk, nfold=n)
        
        x <- cbind(matrix(1, n, 1), phea1[, i, drop = F])
        r6[i] <- cv_fast(fix=x, y=y, kk=kk, nfold=n)
        
        x <- cbind(matrix(1, n, 1), phed1)
        r7[i] <- cv_fast(fix=x, y=y, kk=kk, nfold=n)
        
        x <- cbind(matrix(1, n, 1), phed1[, i, drop = F])
        r8[i] <- cv_fast(fix=x, y=y, kk=kk, nfold=n)
    }
    
    res22[[ii]] <- cbind(r1, r2, r3, r4, r5, r6, r7, r8)
    
}
save(res22, file = "res22.Rdata")

## 
knew <- function(parent, train, add = TRUE) {
    parent <- as.matrix(parent)
    train <- as.matrix(train)
    kk <- list()
    if (add == TRUE) {
        for (i in 1:nrow(parent)) {
            p0 <- (parent[i, ] + parent[-(1:i), ])/2
            kk[[i]] <- p0 %*% t(train)/ncol(train)
        }
    } else {
        for (i in 1:nrow(parent)) {
            p0 <- abs(parent[i, ] - parent[-(1:i), ])/2
            kk[[i]] <- p0 %*% t(train)/ncol(train)
        }
    }
    kknew <- do.call(rbind, kk)
    return(kknew)
}

pp <- list()
for (i in 1:209) {
    name <- row.names(snp0)
    pp[[i]] <- paste(name[i], name[-(1:i)], sep = "/")
}
nameall <- do.call(c, pp)


traita <- list()
for (i in 1:209) {
    traita[[i]] <- (phe0[i, ] + phe0[-(1:i), ])/2
}

pa <- do.call(rbind, traita)
row.names(pa) <- nameall

traitd <- list()
for (i in 1:209) {
    traitd[[i]] <- abs(phe0[i, ] - phe0[-(1:i), ])/2
}
pd <- do.call(rbind, traitd)
row.names(pd) <- nameall


nrow(pa)
kn1 <- knew(parent = snp0, train = snpa, add = TRUE)
kn2 <- knew(parent = snp0, train = snpd, add = FALSE)

kn3 <- knew(parent = exp0, train = expa, add = TRUE)
kn4 <- knew(parent = exp0, train = expd, add = FALSE)

kn5 <- knew(parent = met0, train = meta, add = TRUE)
kn6 <- knew(parent = met0, train = metd, add = FALSE)

n2 <- nrow(pa)
######################### 

ki1 <- list()
ki1[[1]] <- k1
ki1[[2]] <- k2


ke1 <- list()
ke1[[1]] <- kn1
ke1[[2]] <- kn2


ki2 <- list()
ki2[[1]] <- k3
ki2[[2]] <- k4

ke2 <- list()
ke2[[1]] <- kn3
ke2[[2]] <- kn4


ki3 <- list()
ki3[[1]] <- k5
ki3[[2]] <- k6

ke3 <- list()
ke3[[1]] <- kn5
ke3[[2]] <- kn6

ki4 <- list()
ki4[[1]] <- k1
ki4[[2]] <- k2
ki4[[3]] <- k3
ki4[[4]] <- k4

ke4 <- list()
ke4[[1]] <- kn1
ke4[[2]] <- kn2
ke4[[3]] <- kn3
ke4[[4]] <- kn4


ki5 <- list()
ki5[[1]] <- k1
ki5[[2]] <- k2
ki5[[3]] <- k5
ki5[[4]] <- k6

ke5 <- list()
ke5[[1]] <- kn1
ke5[[2]] <- kn2
ke5[[3]] <- kn5
ke5[[4]] <- kn6


ki6 <- list()
ki6[[1]] <- k3
ki6[[2]] <- k4
ki6[[3]] <- k5
ki6[[4]] <- k6

ke6 <- list()
ke6[[1]] <- kn3
ke6[[2]] <- kn4
ke6[[3]] <- kn5
ke6[[4]] <- kn6

ki7 <- list()
ki7[[1]] <- k1
ki7[[2]] <- k2
ki7[[3]] <- k3
ki7[[4]] <- k4
ki7[[5]] <- k5
ki7[[6]] <- k6

ke7 <- list()
ke7[[1]] <- kn1
ke7[[2]] <- kn2
ke7[[3]] <- kn3
ke7[[4]] <- kn4
ke7[[5]] <- kn5
ke7[[6]] <- kn6

# phe.r11 phe.r22 pred
pred <- list()
for (ii in 1:7) {
    cat(ii)
    ktrain <- get(paste("ki", ii, sep = ""))
    ktest <- get(paste("ke", ii, sep = ""))
    r1 <- list()
    for (i in 1:4) {
        cat(i)
        y <- as.matrix(phe12[, i])
        n <- length(y)
        # x<- cbind(matrix(1,n,1),phea1,phed1)
        x <- matrix(1, n, 1)
        xnew <- matrix(1, n2, 1)
        r1[[i]] <- predict(x, y, kk = ktrain, xnew = xnew, knew = ktest)
        max(r1[[i]])
    }
    rr <- do.call(cbind, r1)
    pred[[ii]] <- rr
    
}
save(x = pred, file = "output\\pred.Rdata")

apply(phe12, 2, mean)
predp <- list()
for (ii in 1:7) {
    cat(ii)
    ktrain <- get(paste("ki", ii, sep = ""))
    ktest <- get(paste("ke", ii, sep = ""))
    r1 <- list()
    for (i in 1:4) {
        cat(i)
        y <- as.matrix(phe12[, i])
        n <- length(y)
        x <- cbind(matrix(1, n, 1), phea, phed)
        # x<-matrix(1,n,1)
        xnew <- cbind(matrix(1, n2, 1), pa, pd)
        r1[[i]] <- predict(x, y, kk = ktrain, xnew = xnew, knew = ktest)
        max(r1[[i]])
    }
    rr <- do.call(cbind, r1)
    predp[[ii]] <- rr
    
}
save(x = predp, file = "output\\predp.Rdata")
################################################## output

pre1 <- list()
for (i in 1:4) {
    pp <- lapply(pred, function(x) {
        a <- x[, i]
        z <- sort(a, decreasing = T)[c(1:200)]
        f <- sort(a, decreasing = F)[c(1:200)]
        res <- c(mean(z), sd(z), mean(f), sd(f))
        return(res)
    })
    pre1[[i]] <- do.call(rbind, pp)
}

pre <- do.call(cbind, pre1)

pre2 <- list()
for (i in 1:4) {
    pp <- lapply(predp, function(x) {
        a <- x[, i]
        z <- sort(a, decreasing = T)[c(1:200)]
        f <- sort(a, decreasing = F)[c(1:200)]
        res <- c(mean(z), sd(z), mean(f), sd(f))
        return(res)
    })
    pre2[[i]] <- do.call(rbind, pp)
}
prep <- do.call(cbind, pre2)
predres <- rbind(pre, prep)
row.names(predres) <- c("G", "T", "M", "GT", "GM", "TM", "GTM", "GP", "TP", "MP", "GTP", "GMP", "TMP", 
    "GTMP")
write.csv(x = predres, file = "output\\predres.csv")


