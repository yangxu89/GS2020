
cv <- function(fix = NULL, y, kk, nfold = 5, seed = 123) {
    n <- length(y)
    y <- as.matrix(y)
    if (is.null(fix)) {
        fix <- matrix(1, n, 1)
    } else {
        fix <- as.matrix(fix)
    }
    set.seed(seed)
    foldid <- sample(rep(1:nfold, ceiling(n/nfold))[1:n])
    yobs <- NULL
    yhat <- NULL
    fold <- NULL
    id <- NULL
    k11 <- list()
    k21 <- list()
    g <- length(kk)
    mixed <- function(fix = NULL, y, kk) {
        n <- length(y)
        y <- as.matrix(y)
        if (is.null(fix)) {
            fix <- matrix(1, n, 1)
        } else {
            fix <- as.matrix(fix)
        }
        g <- length(kk)
        nloglik_REML <- function(pm) {
            v_phi <- 0
            for (p in 1:g) {
                v_phi <- v_phi + kk[[p]] * pm[p]
            }
            v_sigma <- diag(n) * pm[g + 1]
            v <- v_phi + v_sigma + diag(1e-09, n)
            v_i <- solve(v, tol = -50)
            beta <- solve(t(fix) %*% v_i %*% fix, t(fix) %*% v_i %*% y)
            nloglik <- 0.5 * (unlist(determinant(v))[[1]] + unlist(determinant(t(fix) %*% v_i %*% 
                fix))[[1]] + t(y - fix %*% beta) %*% v_i %*% (y - fix %*% beta))
            return(nloglik)
        }
        parm0 <- rep(1, g + 1)
        parm <- optim(par = parm0, fn = nloglik_REML, method = "L-BFGS-B", hessian = FALSE, lower = 0)
        v_phi <- 0
        for (p in 1:g) {
            v_phi <- v_phi + kk[[p]] * parm$par[p]
        }
        v_sigma <- diag(n) * parm$par[g + 1]
        v <- v_phi + v_sigma
        v_i <- solve(v, tol = -50)
        beta <- solve(t(fix) %*% v_i %*% fix, t(fix) %*% v_i %*% y)
        res <- list(v_i = v_i, var = parm$par[-(g + 1)], ve = parm$par[g + 1], beta = beta)
        return(res)
    }
    
    cl.cores <- detectCores()
    if (cl.cores <= 2) {
        cl.cores <- 1
    } else if (cl.cores > 2) {
        if (cl.cores > 10) {
            cl.cores <- 10
        } else {
            cl.cores <- detectCores() - 1
        }
    }
    cl <- makeCluster(cl.cores)
    registerDoParallel(cl)
    respar <- foreach(k = 1:nfold, .multicombine = TRUE, .combine = "rbind") %dopar% {
        # for (k in 1:nfold) {
        i1 <- which(foldid != k)
        i2 <- which(foldid == k)
        x1 <- fix[i1, , drop = F]
        y1 <- y[i1, , drop = F]
        x2 <- fix[i2, , drop = F]
        y2 <- y[i2, , drop = F]
        for (i in 1:g) {
            k11[[i]] <- kk[[i]][i1, i1]
            k21[[i]] <- kk[[i]][i2, i1]
        }
        parm <- mixed(fix = x1, y = y1, kk = k11)
        G21 <- 0
        for (i in 1:g) {
            G21 <- G21 + k21[[i]] * parm$var[i]
        }
        v_i <- parm$v_i
        beta <- parm$beta
        y3 <- x2 %*% beta + G21 %*% v_i %*% (y1 - x1 %*% beta)
        pred <- data.frame(i2, yobs = y2, yhat = y3)
    }
    stopCluster(cl)
    R2 <- cor(respar$yobs, respar$yhat)^2
    return(R2)
}

