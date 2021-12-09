is_negative_definite <- function(x){
    vls <- eigen(-x)$values
    ifelse(all(vls > 0) & min(abs(vls)) > 1E-08, TRUE, FALSE)
}

compare_gradient <- function(f, start){
    fs <- function(x) sum(as.numeric(f(x)))
    fstart <- f(start)
    ngrad <- numDeriv::grad(fs, start)
    agrad <- apply(attr(fstart, "gradient"), 2, sum)
    value <- fs(start)
    parts <- attr(fstart, "parts")
    grad <- data.frame(param = start,
                       analytic = agrad,
                       numerical = ngrad,
                       rel_diff = abs(agrad - ngrad) / (agrad + ngrad) * 2)
    return(list(value = value,
                gradient = grad,
                parts = parts))
}


log2 <- function(x) ifelse(x > 0,log(x), 0)

mills <- function(x) exp(dnorm(x, log = TRUE) - pnorm(x, log.p = TRUE))
dmills <- function(x) - mills(x) * (x + mills(x))

# a function to construct block-diagonal matrix (can't remember from
# which package it is borrowed)
bdiag <- function(...){
  if (nargs() == 1)
    x <- as.list(...)
  else
    x <- list(...)
  n <- length(x)
  if(n == 0) return(NULL)
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

# print a summary of the non-linear opimization
print.est.stat <- function(x, ...){
  et <- x$elaps.time[3]
  i <- x$nb.iter[1]
  halton <- x$halton
  method <- x$method
  if (!is.null(x$type) && x$type != "simple"){
    R <- x$nb.draws
    cat(paste("Simulated maximum likelihood with", R, "draws\n"))
  }
  s <- round(et,0)
  h <- s%/%3600
  s <- s-3600*h
  m <- s%/%60
  s <- s-60*m
  cat(paste(method, "method\n"))
  tstr <- paste(h, "h:", m, "m:", s, "s", sep="")
  cat(paste(i,"iterations,",tstr,"\n"))
  if (!is.null(halton)) cat("Halton's sequences used\n")
  if (!is.null(x$eps)) cat(paste("g'(-H)^-1g =", sprintf("%5.3G", as.numeric(x$eps)),"\n"))
  if (is.numeric(x$code)){
    msg <- switch(x$code,
                  "1" = "gradient close to zero",
                  "2" = "successive fonction values within tolerance limits",
                  "3" = "last step couldn't find higher value",
                  "4" = "iteration limit exceeded"
                  )
    cat(paste(msg, "\n"))
  }
  else cat(paste(x$code, "\n"))
}

# Compute the estimation of one-equation models. Those are only
# relevant for distributions that admit negative values of y^*, namely
# ln2, bc2, ihs and n


# Compute the "naive" model, i.e. a model with no explanatory
# variables.  Full version with correlation ; not used because
# identification problems

lnl.naive <- function(param, dist = c("ln", "tn", "n", "ln2"), moments,
                      h1 = TRUE, h3 = FALSE){
    dist <- match.arg(dist)
    n <- moments[1]
    ym <- moments[2]
    s2 <- moments[3]
    if (h1){
        alpha1 <- param[1]
        alpha2 <- param[2]
        param <- param[- c(1,2)]
    }
    else{
        alpha2 <- param[1]
        param <- param[- 1]
    }
    if (h3){
        alpha3 <- param[1]
        param <- param[- 1]
    }
    sigma <- param[1]
    
    if (h1){
        Phi1 <- pnorm(alpha1)
        phi1 <- dnorm(alpha1)
    }
    else Phi1 <- 1
    if (h3){
        Phi3 <- pnorm(alpha3)
        phi3 <- dnorm(alpha3)
    }
    else Phi3 <- 1
    Phi2 <- pnorm(alpha2 / sigma)
    phi2 <- dnorm(alpha2 / sigma)
    scr <- ifelse(dist == "ln",
                  s2 + (ym + log(Phi3) - alpha2) ^ 2,
                  Phi3 ^ 2 * (s2 + (ym - alpha2 / Phi3) ^ 2)
                  )
    Pbiv <- Phi1 * Phi2
    Phi1bis <- Phi1
    s2term <- 0
    P0 <- switch(dist,
                 "ln" = 1 - Phi1 * Phi3,
                 "ln2" = 1 - Phi1 * Phi3,
                 "tn" = 1 - Pbiv / Phi2 * Phi3,
                 "n" = 1 - Pbiv * Phi3
                 )
    lnPos <-
        -log(sigma) - 0.5 * log(2 * pi) -
            scr / (2 * sigma ^ 2) +
                log(Phi3) +
                    (log(Phi1bis) + s2term) * h1 -
                        ym * (dist == "ln") +
                            log(Phi3) * (dist != "ln") -
                                log(Phi2) * (dist == "tn")
    
    lnL <- n * log(P0) + (1 - n) * lnPos
    #YC 20180110 extract the parts in order to compute the new R2s
    attr(lnL, "parts") <- c(lnLNull = log(P0), lnLOne = log(1 - P0), lnLPos = lnPos)
    lnL
}


lnl.naive <- function(param, dist = c("ln", "tn", "n", "ln2"), moments,
                      h1 = TRUE, h3 = FALSE){
    dist <- match.arg(dist)
    n <- moments[1]
    ym <- moments[2]
    s2 <- moments[3]
    if (h1){
        alpha1 <- param[1]
        alpha2 <- param[2]
        param <- param[- c(1,2)]
    }
    else{
        alpha2 <- param[1]
        param <- param[- 1]
    }
    if (h3){
        alpha3 <- param[1]
        param <- param[- 1]
    }
    sigma <- param[1]
    
    if (h1){
        Phi1 <- pnorm(alpha1)
        phi1 <- dnorm(alpha1)
    }
    else Phi1 <- 1
    if (h3){
        Phi3 <- pnorm(alpha3)
        phi3 <- dnorm(alpha3)
    }
    else Phi3 <- 1
    Phi2 <- pnorm(alpha2 / sigma)
    phi2 <- dnorm(alpha2 / sigma)
    scr <- ifelse(dist == "ln",
                  s2 + (ym + log(Phi3) - alpha2) ^ 2,
                  Phi3 ^ 2 * (s2 + ym ^ 2) + alpha2 ^ 2 - 2 * alpha2 * Phi3 * ym
                  )
    Pbiv <- Phi1 * Phi2
    Phi1bis <- Phi1
    s2term <- 0
    P0 <- switch(dist,
                 "ln" = 1 - Phi1 * Phi3,
                 "ln2" = 1 - Phi1 * Phi3,
                 "tn" = 1 - Pbiv / Phi2 * Phi3,
                 "n" = 1 - Pbiv * Phi3
                 )
    lnPos <-
        - log(sigma) - 0.5 * log(2 * pi) -
            scr / (2 * sigma ^ 2) +
                log(Phi3) +
                    log(Phi1) -
                        ym * (dist == "ln") +
                            log(Phi3) * (dist != "ln") -
                                log(Phi2) * (dist == "tn")
    
    lnL <- n * log(P0) + (1 - n) * lnPos
    #YC 20180110 extract the parts in order to compute the new R2s
    attr(lnL, "parts") <- c(lnLNull = log(P0), lnLOne = log(1 - P0), lnLPos = lnPos)
    lnL
}



## Ostart.mhurdle <- function(X1, X2, X3, y, dist){
##     h1 <- ! is.null(X1)
##     h3 <- ! is.null(X3)
##     # for models (2ld, 2td), estimates of models (2li, 2ti) which can be
##     # estimated simply in two parts using seperate.mhurdle() are used
##     # as starting values
##     if (h1 && ! h3 &&  dist %in% c("ln", "tn")){
##         result <- seperate.mhurdle(X1, X2, y, dist = dist)
##         start <- result$coefficients
##     }
##     # model (3) tobit : use linear model as starting values
##     if (!h1 && !h3 &&  dist == "n"){
##         lin <- lm(y ~ X2 - 1)
##         sdest <- summary(lin)$sigma
##         start <- c(coef(lin), sigma = sdest)
##     }
##     # model (4, 7) h3 whithout h1
##     if (!h1 && h3){
##         probit <- glm( (y != 0) ~ X3 - 1, family = binomial(link = 'probit'))
##         bX3 <- as.numeric(crossprod(coef(probit), t(X3)))
##         Phi3 <- pnorm(bX3)
##         yPP <- y * Phi3
##         lin <- switch(dist,
##                       "ln" = lm(log(yPP)  ~ X2 - 1, subset = y != 0),
##                       "n" = lm(yPP        ~ X2 - 1, subset = y != 0),
##                       "tn" = truncreg::truncreg(yPP ~ X2 - 1, subset = y != 0, scaled = TRUE)
##                       )
##         if (dist %in% c("n", "ln")) start <- c(coef(lin), coef(probit), sigma = summary(lin)$sigma)
##         if (dist == "tn") start <- c(coef(lin)[- length(coef(lin))], coef(probit), sigma = coef(lin)[length(coef(lin))])
##     }
##     # model (5), double hurdle use model (3i) as starting values
##     if (h1 && !h3 && dist == "n"){
##         result <- seperate.mhurdle(X1, X2, y, dist = dist)
##         start <- result$coefficients
##     }
##     # model (6 and 8)
##     if (h1 && h3){
##         probit.h3 <- glm( (y != 0) ~ X3 - 1, family = binomial(link = 'probit'))
##         probit.h1 <- glm( (y != 0) ~ X1 - 1, family = binomial(link = 'probit'))
##         beta3 <- coef(probit.h3)
##         beta1 <- coef(probit.h1)
##         bX3 <- as.numeric(crossprod(beta3, t(X3)))
##         Phi3 <- pnorm(bX3)
##         bX1 <- as.numeric(crossprod(beta1, t(X1)))
##         P0 <- mean(y > 0)
##         yPP <- y * Phi3
##         lin <- switch(dist,
##                       "ln" = lm(log(yPP)  ~ X2 - 1, subset = y!= 0),
##                       "n" = lm(yPP       ~ X2 - 1, subset = y!= 0),
##                       "tn" = truncreg::truncreg(yPP ~ X2 - 1, subset = y!= 0, scaled = TRUE)
##                       )
##         if(dist != "tn")
##             start <- c(beta1, coef(lin), beta3, sigma = summary(lin)$sigma)
##         else
##             start <- c(beta1, coef(lin)[- length(coef(lin))],
##                        beta3, sigma = coef(lin)[length(coef(lin))])
##     }
##     start
## }


start.mhurdle <- function(X1, X2, X3, y, dist){
    if (! is.null(X1)) beta1 <- coef(glm( (y != 0) ~ X1 - 1, family = binomial(link = 'probit'))) else beta1 <- NULL
    if (! is.null(X3)) beta3 <- coef(glm( (y != 0) ~ X3 - 1, family = binomial(link = 'probit'))) else beta3 <- NULL
    if (dist == "n") lin <- lm(y ~ X2 - 1)
    if (dist == "tn") lin <- lm(y ~ X2 - 1, subset = y > 0); truncreg::truncreg(y ~ X2 - 1, subset = y > 0)
    if (dist == "ln") lin <- lm(log(y) ~ X2 - 1, subset = y > 0)
    beta2 <- coef(lin)[1:ncol(X2)]
#    beta2 <- coef(lin[1:ncol(X2)])#YC20171204 wrong parenthesis
#    if (dist %in% c("n", "tn")) sigma <- summary(lin)$sigma
#    else sigma <- coef(lin)[ncol(X2) + 1]
    sigma <- summary(lin)$sigma
    c(beta1, beta2, beta3, sigma = sigma)
}
