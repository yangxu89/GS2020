dir <- "/home/xuyang/Data/phe"
setwd(dir)
# load(file='mixed.RData')
source(file = "cross.R")
source(file = "all.R")

library(ggplot2)
library(ggthemes)
library(ggpubr)
library(reshape2)


phe <- read.csv(file = "Input\\RIL-Phenotypes.csv", header = T)
snp <- read.csv(file = "Input\\RIL-Genotypes.csv", header = T)
exp <- read.csv(file = "Input\\RIL-Expressions.csv", header = T)
met <- read.csv(file = "Input\\RIL-Metabolites.csv", header = T)


head(met)
n <- nrow(phe)

gn <- (phe$gn97x + phe$gn98x + phe$gn98h + phe$gn99h)/4
tp <- (phe$tp97x + phe$tp98x + phe$tp98h + phe$tp99h)/4
yd <- (phe$yd97x + phe$yd98x + phe$yd98h + phe$yd99h)/4
kg <- (phe$kgw97x + phe$kgw98x + phe$kgw98h + phe$kgw99h)/4


phe1 <- cbind(yd, tp, gn, kg)
phe0 <- scale(phe1)
row.names(phe0) <- row.names(phe1)
phef <- read.csv(file = "Input\\IMF2-Phenotypes.csv", header = T)

# y<-yd y<- as.matrix(y) x<-rep(1,length(y)) x<-as.matrix(x)
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

# snp.r11<-scale(t(snp.r1)) snp.r22<-scale(t(snp.r2)) exp.r11<-scale(t(exp.r1))
# exp.r22<-scale(t(exp.r2)) met.r11<-scale(t(met.r1)) met.r22<-scale(t(met.r2))
# phe.r11<-scale(phe.r1) phe.r22<-scale(phe.r2)


########### IMF2 gene phe snpa<-0.5*(snp.r1+snp.r2) snpd<-0.5*abs(snp.r1-snp.r2)

# meta<-0.5*(met.r1+met.r2) metd<-0.5*abs(met.r1-met.r2) expa<-0.5*(exp.r1+exp.r2)
# expd<-0.5*abs(exp.r1-exp.r2) phea<-0.5*(phe.r1+phe.r2) phed<-0.5*abs(phe.r1-phe.r2) IMF2 gene
# phe
snpa <- 0.5 * (snp.r11 + snp.r22)
snpd <- 0.5 * abs(snp.r11 - snp.r22)

expa <- 0.5 * (exp.r11 + exp.r22)
expd <- 0.5 * abs(exp.r11 - exp.r22)

meta <- 0.5 * (met.r11 + met.r22)
metd <- 0.5 * abs(met.r11 - met.r22)

phea <- 0.5 * (phe.r11 + phe.r22)
phed <- 0.5 * abs(phe.r11 - phe.r22)



yield <- (phef$yield98 + phef$yield99)/2
kgw <- (phef$kgw98 + phef$kgw99)/2
grain <- (phef$grain98 + phef$grain99)/2
tiller <- (phef$tiller98 + phef$tiller99)/2
phe12 <- cbind(yield, tiller, grain, kgw)



k1 <- tcrossprod(snpa)/ncol(snpa)
dim(k1)
k2 <- tcrossprod(snpd)/ncol(snpd)

k3 <- tcrossprod(expa)/ncol(expa)
dim(k1)
k4 <- tcrossprod(expd)/ncol(expd)


k5 <- tcrossprod(meta)/ncol(meta)
dim(k1)
k6 <- tcrossprod(metd)/ncol(metd)

k7 <- tcrossprod(phea)/ncol(phea)
k8 <- tcrossprod(phed)/ncol(phed)



z1 <- as.matrix(snp.r11)
k11 <- tcrossprod(z1)/nrow(z1)
dim(k11)
z2 <- as.matrix(snp.r22)
k22 <- tcrossprod(z2)/nrow(z2)

z3 <- as.matrix(exp.r11)
k33 <- tcrossprod(z3)/nrow(z3)

z4 <- as.matrix(exp.r22)
k44 <- tcrossprod(z4)/nrow(z4)

z5 <- as.matrix(met.r11)
k55 <- tcrossprod(z5)/nrow(z5)
dim(k55)
z6 <- as.matrix(met.r22)
k66 <- tcrossprod(z6)/nrow(z6)

z7 <- as.matrix(phe.r11)
k77 <- tcrossprod(z7)/nrow(z7)

z8 <- as.matrix(phe.r22)
k88 <- tcrossprod(z8)/nrow(z8)
dim(k88)
##################################################### all hybrid
nrow(snp0)

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


############################################################################# 

# name1<-list() top1<-list() for (i in 1:4){ pp<-lapply(pred, function(x){ a<-x[,i] names(a)<-
# nameall z<-sort(a,decreasing =T)[c(1:200)] #f<-sort(a,decreasing =F)[c(1:200)]
# res<-data.frame(names(z),z) return(res) } ) top1[[i]]<-do.call(cbind,pp) #top1[[i]]<-pp }
# #pre<-do.call(cbind,pre1) name2<-list() top2<-list() for (i in 1:4){ pp<-lapply(predp,
# function(x){ a<-x[,i] names(a)<- nameall z<-sort(a,decreasing =T)[c(1:200)]
# res<-data.frame(names(z),z) return(res) } ) top2[[i]]<-do.call(cbind,pp) }
# #########################################################


name1 <- list()
top1 <- list()

for (i in 1:4) {
    pp <- lapply(pred, function(x) {
        a <- x[, i]
        names(a) <- nameall
        z <- sort(a, decreasing = T)[c(1:200)]
        # f<-sort(a,decreasing =F)[c(1:200)]
        res <- data.frame(names(z), z)
        return(res)
    })
    top1[[i]] <- pp
    # top1[[i]]<-pp
}

# pre<-do.call(cbind,pre1)
name2 <- list()
top2 <- list()
for (i in 1:4) {
    pp <- lapply(predp, function(x) {
        a <- x[, i]
        names(a) <- nameall
        z <- sort(a, decreasing = T)[c(1:200)]
        res <- data.frame(names(z), z)
        return(res)
    })
    top2[[i]] <- pp
}


omic <- c("G", "T", "M", "GT", "GM", "TM", "GTM", "GP", "TP", "MP", "GTP", "GMP", "TMP", "GTMP")

trait <- c("YIELD", "TILLER", "GRAIN", "KGW")
for (i in 1:4) {
    aa <- top2[[i]]
    for (j in 1:7) {
        cat(j)
        row.names(aa[[j]]) <- c(1:200)
        colnames(aa[[j]]) <- c("IMF2", trait[i])
        write.xlsx(aa[[j]], paste("output\\", trait[i], ".xlsx", sep = ""), omic[[j + 7]], append = TRUE)
    }
    
}
############################################ all

name1 <- list()
all1 <- list()

for (i in 1:4) {
    pp <- lapply(pred, function(x) {
        a <- x[, i]
        names(a) <- nameall
        z <- sort(a, decreasing = T)
        res <- data.frame(names(z), z)
        return(res)
    })
    all1[[i]] <- pp
    # top1[[i]]<-pp
}

pred

# pre<-do.call(cbind,pre1)
name2 <- list()
all2 <- list()
for (i in 1:4) {
    pp <- lapply(predp, function(x) {
        a <- x[, i]
        names(a) <- nameall
        z <- sort(a, decreasing = T)
        res <- data.frame(names(z), z)
        return(res)
    })
    all2[[i]] <- pp
}

omic <- c("G", "T", "M", "GT", "GM", "TM", "GTM", "GP", "TP", "MP", "GTP", "GMP", "TMP", "GTMP")

trait <- c("YIELD", "TILLER", "GRAIN", "KGW")
for (i in 1:4) {
    aa <- all1[[i]]
    for (j in 1:7) {
        cat(j)
        row.names(aa[[j]]) <- c(1:21945)
        colnames(aa[[j]]) <- c("Cross", trait[i])
        write.xlsx(aa[[j]], paste("output\\", trait[i], "nop.xlsx", sep = ""), paste("sheet", j, 
            sep = ""), append = TRUE)
    }
}

nrow(aa[[j]])
write.xlsx(incometsforecasts2, "/Users/ap/Desktop/incometsforecast.xlsx", )

############################################## name
name1 <- list()
for (i in 1:4) {
    pp <- lapply(pred, function(x) {
        a <- x[, i]
        names(a) <- nameall
        z <- sort(a, decreasing = T)[c(1:200)]
        # f<-sort(a,decreasing =F)[c(1:200)]
        res <- names(z)
        return(res)
    })
    name1[[i]] <- do.call(cbind, pp)
}

# pre<-do.call(cbind,pre1)
name2 <- list()
top2 <- list()
for (i in 1:4) {
    pp <- lapply(predp, function(x) {
        a <- x[, i]
        names(a) <- nameall
        z <- sort(a, decreasing = T)[c(1:200)]
        res <- names(z)
        return(res)
    })
    name2[[i]] <- do.call(cbind, pp)
}




for (i in 1:4) {
    aa <- name1[[i]]
    bb <- name2[[i]]
    for (ii in 1:7) {
        d <- intersect(aa[, ii], bb[, ii])
        cat(length(d), sep = "\n")
    }
}



