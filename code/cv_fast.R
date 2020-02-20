cv_fast <- function(fix = NULL, y, kk, nfold = 5, seed = 123) {
    n <- length(y)
    y <- as.matrix(y)
    if (is.null(fix)) {
        fix <- matrix(1, n, 1)
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
        nloglik <- 0.5 * (unlist(determinant(v))[[1]] + unlist(determinant(t(fix) %*% v_i %*% fix))[[1]] + 
            t(y - fix %*% beta) %*% v_i %*% (y - fix %*% beta))
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
    v_i <- solve(v)
    beta <- solve(t(fix) %*% v_i %*% fix, t(fix) %*% v_i %*% y)
    H <- v_phi %*% v_i
    r <- y - fix %*% beta
    SSP <- var(r) * (n - 1)
    rhat <- H %*% r
    ehat <- r - rhat
    if (nfold == n) {
        foldid <- c(1:n)
    } else {
        set.seed(seed)
        foldid <- sample(rep(1:nfold, ceiling(n/nfold))[1:n])
    }
    PRESS <- 0
    for (i in 1:nfold) {
        indx <- which(foldid == i)
        nk <- length(indx)
        Hkk <- H[indx, indx]
        e <- solve(diag(nk) - Hkk) %*% ehat[indx]
        PRESS <- PRESS + sum(e^2)
    }
    R2 <- 1 - PRESS/SSP
    return(R2)
}
