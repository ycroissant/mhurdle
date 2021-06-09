gbdiag <- function(...){
  if (nargs() == 1)
    x <- as.list(...)
  else
    x <- list(...)
  n <- length(x)
  if(n==0) return(NULL)
  x <- lapply(x, function(y) if(length(y)) as.matrix(y) else
              stop("Zero-length component in x"))
  d <- array(unlist(lapply(x, dim)), c(2, n))
  rr <- d[1,]
  cc <- d[2,]
  rsum <- sum(rr)
  csum <- sum(cc)
  out <- array(0, c(rsum, csum))
  ind <- array(0, c(4, n))
  rcum <- cumsum(rr)
  ccum <- cumsum(cc)
  ind[1,-1] <- rcum[-n]
  ind[2,] <- rcum
  ind[3,-1] <- ccum[-n]
  ind[4,] <- ccum
  imat <- array(1:(rsum * csum), c(rsum, csum))
  iuse <- apply(ind, 2, function(y, imat) imat[(y[1]+1):y[2],
                                               (y[3]+1):y[4]], imat=imat)
  iuse <- as.vector(unlist(iuse))
  out[iuse] <- unlist(x)
  return(out)
}

sqrtm <- function(x){
    eigenx <- eigen(x)
#    if (any(eigenx$values < 0)) stop("negative eigen values")
#    eigenx$vectors %*% diag(sqrt(abs(eigenx$values))) %*% t(eigenx$vectors)
    eigenx$vectors %*% diag(sqrt(eigenx$values)) %*% t(eigenx$vectors)
}



#' Vuoung test for non-nested models
#' 
#' The Vuong test is suitable to discriminate between two non-nested models.
#' 
#' 
#' @aliases ndvuongtest
#' @param x a first fitted model of class \code{"mhurdle"},
#' @param y a second fitted model of class \code{"mhurdle"},
#' @param size the size of the test,
#' @param pval should the p-value be computed ?
#' @param type the kind of test to be computed,
#' @param ndraws the number of draws for the simulations,
#' @param diffnorm a creuser,
#' @param seed the seed,
#' @param print.level the level of details to be printed.
#' @return an object of class \code{"htest"}
#' @seealso \code{vuong} in package \code{pscl}.
#' @references Vuong Q.H. (1989) Likelihood ratio tests for model selection and
#' non-nested hypothesis, Econometrica, vol.57(2), pp.307-33.
#' @keywords htest
#' @examples
#' 
#' data("Interview", package = "mhurdle")
#' # dependent double hurdle model
#' dhm <- mhurdle(vacations ~ car + size | linc + linc2 | 0, Interview,
#'               dist = "ln", h2 = TRUE, method = "bhhh", corr = TRUE)
#' 
#' # a double hurdle p-tobit model
#' ptm <- mhurdle(vacations ~ 0 | linc + linc2 | car + size, Interview,
#'               dist = "ln", h2 = TRUE, method = "bhhh", corr = TRUE)
#' vuongtest(dhm, ptm)
#' @export
ndvuongtest <- function(x, y, size = 0.05, pval = TRUE,
                        type = c("non-nested", "nested", "overlapping"),
                        ndraws = 1E04, diffnorm = 0.1, seed = 1,
                        print.level = 0){
    type <- match.arg(type)
    eigentol <- 1E-12
    data.name <- c(
        paste(deparse(substitute(x))),
        paste(deparse(substitute(y)))
    )
    data.name <- paste(data.name, collapse = "-")
    set.seed(seed)

    if (!is.null(x$gradientObs)) x$gradient <- x$gradientObs
    if (!is.null(y$gradientObs)) y$gradient <- y$gradientObs
    if (!is.null(x$lnl)) x$logLik <- x$lnl 
    if (!is.null(y$lnl)) y$logLik <- y$lnl 
    
    fmodel <- x                    ;  gmodel <- y    
    hessian.f <- fmodel$hessian    ;  hessian.g <- gmodel$hessian
    # this version is uncorrect where ther
    hessian.f <- - solve(vcov(fmodel))    ;  hessian.g <- - solve(vcov(gmodel))
    grad.f <- fmodel$gradient      ;  grad.g <- gmodel$gradient
    
    N <- nrow(grad.f)
    lnl.f <- fmodel$logLik         ;  lnl.g <- gmodel$logLik
    K.f <- ncol(grad.f)            ;  K.g <- ncol(grad.g)
    K <- K.f + K.g
    A <- bdiag(hessian.f, - hessian.g) / N
    solveA <- bdiag(solve(hessian.f / N), - solve(hessian.g / N))
    grad.tot <- cbind(grad.f, - grad.g)
    grad.tot <- t(t(grad.tot) - apply(grad.tot, 2, mean))
    B <- crossprod(grad.tot) / N + diag(rep(eigentol, K))
    sqrtB <- sqrtm(B)
    V <- sqrtB %*% solveA %*% sqrtB
    V <- eigen(V, symmetric = TRUE)$values
    V <- round(V, 5)
    Vf <- - sum(eigen(solve(hessian.f / N) %*%
                      (crossprod( t(t(grad.f) - apply(grad.f, 2, mean))) / N),
                      only.values = TRUE)$values)
    Vg <- - sum(eigen(solve(hessian.g / N) %*%
                      (crossprod( t(t(grad.g) - apply(grad.g, 2, mean))) / N),
                      only.values = TRUE)$values)
    if (type != "nested"){
        rho <- as.numeric((abs(V) - max(abs(V))) == 0)
        rho <- rho / sqrt(sum(rho))
    }
    else rho <- rep(0, length(V))
    SIGMA <- rbind(c(1 , rho),
                   cbind(rho, diag(rep(1, K))))
    Z <- matrix(rnorm(ndraws * (K + 1)), ndraws)
#    Z <- Z %*% sqrtm(SIGMA + diag(rep(eigentol, K + 1)))
    Z <- Z %*% sqrtm(SIGMA)
    ZL <- Z[,   1]
    ZP <- Z[, - 1]
    trVsq <- sum(V ^ 2)
    Vnmlzd <- V / sqrt(trVsq)
    A1 <- ZL
    A2 <- as.numeric(apply(ZP, 1, function(x) sum(Vnmlzd * x ^ 2) / 2) -
                     sum(Vnmlzd) / 2)
    A3 <- as.numeric(apply(ZP, 1, function(x) sum(Vnmlzd * rho * x)))
    A4 <- as.numeric(apply(ZP, 1, function(x) sum(x ^ 2 * Vnmlzd ^ 2)))
    
    LR <- mean(lnl.f - lnl.g)
    w2 <- mean( (lnl.f - lnl.g) ^ 2) - LR ^ 2
    Tvuong <- sqrt(N) * LR / sqrt(w2)
    
    Tmod <- function(sigma, cst){
        # R draws in the distribution of the Vuong statistic
        num <- sigma * A1 - A2
        denom <- sigma ^ 2 - 2 * sigma * A3 + A4 + cst
        num / sqrt(denom)
    }

    quant <- function(sigma, cst, size)
        # compute the empirical quantile of level 1 - size in the
        # distribution of the Vuong statistic
        as.numeric(quantile(abs(Tmod(sigma, cst)), 1 - size))
        
    sigstar <- function(cst, size)
        # compute the value of sigma which maximize the empirical
        # quantile of level 1 - size of the distribution of the Vuong
        # statistic
        optimize(function(x) quant(x, cst, size), c(0, 5), maximum = TRUE)$maximum

    seekpval <- function(size){
        # for a given size, compute the constant so that the critical
        # value equals the target
        c.value <- 0
        cv.value <- quant(sigstar(0, size), 0, size)
        cv.normal <- qnorm(1 - size / 2)
        cv.target <- cv.normal + diffnorm
       
        if (cv.value < cv.target){
            cv.value <- max(cv.value, cv.normal)
        }
        else {
            froot <- function(cst) quant(sigstar(cst, size), cst, size) - cv.target
            zo <- uniroot(froot, c(0, 10))
            c.value <- zo$root
            cv.value <- quant(sigstar(c.value, size), c.value, size)
        }
        LRmod <- LR + sum(V) / (2 * N)
        w2mod <- w2 + c.value * sum(V ^ 2) / N
        Tnd <- sqrt(N) * (LRmod) / sqrt(w2mod)      
        pvalue <- 1 - ecdf(abs(Tmod(sigstar(c.value, size), c.value)))(abs(Tnd))
        list(pvalue = pvalue, stat = Tnd,
             constant = c.value, cv = cv.value)#,
    }

    if (type == "nested"){
        LRmod <- LR + sum(V) / (2 * N)
        Tnd <- sqrt(N) * (LRmod) / sqrt(w2)
        if (Tnd < 0) pvalue <- ecdf(Tmod(0, 0))(Tnd)
        else pvalue <- 1 - ecdf(Tmod(0, 0))(Tnd)
        TndVuong <- structure(list(statistic = c(z = Tnd),
                                   method = "Non-degenerate Vuong test",
                                   p.value = pvalue,
                                   data.name = data.name,
                                   alternative = "different models"),
                              class = "htest")
        trad <- vuongtest(x, y, type = "nested")
        structure(list(trad = trad, nd = TndVuong),
                  class = "ndvuongtest")
    }
    else{
        if (pval){
            froot <- function(alpha) seekpval(alpha)$pvalue - alpha           
            pvalue <- uniroot(froot, c(1E-7, 1-1E-7))
            results <- seekpval(pvalue$root)
            Tndvuong <- structure(list(statistic = c(z = as.numeric(results$stat)),
                                       method = "Non-degenerate Vuong test",
                                       p.value = results$pvalue,
                                       data.name = data.name,
                                       alternative = "different models",
                                       parameters = c(constant = results$constant,
                                                      'crit-value' = results$cv,
                                                      'sum e.v.' = sum(V),
                                                      'tracef' = Vf,
                                                      'traceg' = Vg)),
                                  class = "htest")
        }
        else{
            results <- seekpval(size)
            Tndvuong <- structure(list(statistic = c(z = as.numeric(results$stat)),
                                       method = "Non-degenerate Vuong test",
                                       p.value = results$pvalue,
                                       data.name = data.name,
                                       size = size,
                                       alternative = "different models",
                                       parameters = c(constant = results$constant,
                                                      'crit. value' = results$cv,
                                                      'sum e.v.' = sum(V))),
                                  class = "htest")
        }
        Tvuong <- structure(list(statistic = c(z = Tvuong),
                                 method = "Vuong test",
                                 p.value = 2 * pnorm(abs(Tvuong), lower.tail = FALSE),
                                 data.name = data.name,
                                 alternative = "different models"),
                            class = "htest")
        
        structure(list(trad = Tvuong, nd = Tndvuong), class = "ndvuongtest")
    }
}

print.ndvuongtest <- function(x, model = c("nd", "trad", "both"), ...){
    model <- match.arg(model)
    if (model != "nd") print(x$trad, ...)
    if (model != "trad") print(x$nd, ...)
}
