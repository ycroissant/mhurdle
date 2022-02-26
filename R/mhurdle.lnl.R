mhurdle.lnl <- function(param, X1, X2, X3, X4, y, gradient = FALSE,
                        fitted = FALSE, dist = NULL, corr = FALSE,
                        robust = TRUE){

    #----------------------------------------------------------------
    # 1. Declare the transformation functions for the parameters and
    # the associated gradients
    #----------------------------------------------------------------

    # if robust, compute the transformation of the parameter and its
    # gradient in order to get the structural parameter, ie, 0 is
    # entered for sigma and is actually log(sigma), so that sigma =
    # exp(0) = 1

    if (robust){
        frho <- function(x) atan(x) * 2 / pi
        grho <- function(x) 2 / pi / (1 +  x ^ 2)
        fmu <- function(x) exp(x)
        gmu <- function(x) exp(x)
        fsd <- function(x) exp(x)
        gsd <- function(x) exp(x)
    }
    else{
        frho <- function(x) x
        grho <- function(x) rep(1, length(x))
        fmu <- function(x) x
        gmu <- function(x) 1
        fsd <- function(x) x
        gsd <- function(x) 1
    }

    #-----------------------------------------------------
    # 2. Extract the coefficients and compute the relevant
    # crossproducts
    #-----------------------------------------------------

    ##.........................................................
    ## 2.a  Dummies for the existing equations hi and number of
    ##  coefficients Ki
    ##.........................................................

    h1 <- ! is.null(X1) ;  K1 <- ifelse(h1, ncol(X1), 0)
    h3 <- ! is.null(X3) ;  K3 <- ifelse(h3, ncol(X3), 0)
    h4 <- ! is.null(X4) ;  K4 <- ifelse(h4, ncol(X4), 0)
    K2 <- ncol(X2)

    ##.............................
    ##   2.b Equation 1 (selection)
    ##.............................

    if (h1){
        beta1 <- param[1:K1]
        bX1 <- as.numeric(crossprod(t(X1), beta1))
        Phi1 <- pnorm(bX1) ; phi1 <- dnorm(bX1)
        Phi1 <- sanitize(Phi1, string = c("Phi1", "mhurdle.lnl"),
                         verbal = FALSE, replace = TRUE)
    }
    else{
        bX1 <- Inf ; beta1 <- NULL
        Phi1 <- 1 ; phi1 <- 0;
    }

    ##.........................
    ##   2.c Equation 2 (level)
    ##.........................

    beta2 <- param[(K1 + 1):(K1 + K2)]
    bX2 <- as.numeric(crossprod(t(X2), beta2))

    ##............................
    ##   2.d Equation 3 (purchase)
    ##............................

    if (h3){
        beta3 <- param[(K1 + K2 + 1):(K1 + K2 + K3)]
        bX3 <- as.numeric(crossprod(t(X3), beta3))
        Phi3 <- pnorm(bX3) ; phi3 <- dnorm(bX3)
        Phi3 <- sanitize(Phi3, string = c("Phi3", "mhurdle.lnl"),
                         verbal = FALSE, replace = TRUE)
    }
    else{
        bX3 <- Inf ; beta3 <- NULL
        Phi3 <- 1 ; phi3 <- 0
    }
    
    ##.........................
    ##   2.e Standard deviation
    ##.........................

    sd <- param[K1 + K2 + K3 + 1]
    gradientsd <- gsd(sd)
    sd <- fsd(sd)
    
    ##.......................................
    ##   2.f Equation 4 (heterosckedasticity)
    ##.......................................

    if (h4){
        beta4 <- param[(K1 + K2 + K3 + 2):(K1 + K2 + K3 + K4 + 1)]
        bX4 <- as.numeric(crossprod(t(X4), beta4))
        pbX4 <- pnorm(bX4)
        #PKG sigma <- sd * pbX4
    }
    else{
        sigma <- sd ; pbX4 <- 1
    }

    ##...............................
    ##   2.g Correlation coefficients
    ##...............................

    # KR is the length of rho, either 1 (h1 or h3) or 3 (h1 and h3)
    KR <- ifelse(corr, h1 + h3 + h1 * h3, 0)
    if (corr){
        posrho <- (K1 + K2 + K3 + K4 + 2):(K1 + K2 + K3 + K4 + 1 + KR)
        rho <- param[posrho]
        # In case of only one correlation coefficient, use the whole
        # vector with only one non-null component
        if (h1 & ! h3) rho <- c(rho, 0, 0)
        if (! h1 & h3) rho <- c(0, 0, rho)
        gradientrho <- grho(rho)
        rho <- frho(rho)
        if (h1 & h3){
            # In case of a trivariate normal distribution, check the
            # joint relation of the three coefficients
            rho3 <- function(rho) rho[1] * rho[2] + c(-1, 1) *
                sqrt(1 + rho[1] ^ 2 * rho[2] ^ 2 - rho[1] ^ 2 - rho[2] ^ 2)
            if (rho[3] < rho3(rho)[1]) rho[3] <- rho3(rho)[1] + 1E-04 
            if (rho[3] > rho3(rho)[2]) rho[3] <- rho3(rho)[2] - 1E-04
        }
    }
    else rho <- rep(0, 3)

    ##........................................................
    ##   2.h Transformation coefficient (Box-Cox and ish only)
    ##........................................................

    if (dist %in% c("bc", "bc2", "ihs")) lambda <- param[K1 + K2 + K3 + 1 + K4 + KR + 1]
    if (dist %in% c("bc", "bc2")) sgn <- sign(lambda) else sgn <- + 1

    ##...............................................
    ##   2.i Position parameter (ln and Box-Cox only)
    ##...............................................

    if (dist %in% c("ln2", "bc2")){
        posmu <- K1 + K2 + K3 + 1 + K4 + KR + ifelse(dist == "ln2", 1, 2)
        mu <- fmu(param[posmu])
        gradientmu <- gmu(param[posmu])
    }

    #----------------------------------------------
    # 3. Compute the elements of the log likelihood
    #----------------------------------------------

    ##................................................
    ##   3.a. Transformation of the dependent variable
    ##................................................

    Ty <- switch(dist,
                 "ln"  = log2(y) + log(Phi3),
                 "ln2" = log2(y * Phi3 + mu),
                 "bc"  = (exp(lambda * log(y * Phi3)) - 1) / lambda,
                 "bc2" = (exp(lambda * log(y * Phi3 + mu)) - 1) / lambda,
                 "ihs" = log(lambda * y * Phi3 + sqrt(1 + (lambda  * y * Phi3) ^ 2)) / lambda,
                 y * Phi3
                 )
    
    ##................................
    ##   3.b Logarithm of the jacobian
    ##................................

    lnJ <- switch(dist,
                  "ln"  = - log2(y),
                  "ln2" = log(Phi3) - log(mu + Phi3 * y),
                  "bc"  = (lambda - 1) * log2(y) + lambda * log(Phi3),
                  "bc2" = (lambda - 1) * log(Phi3 * y + mu) + log(Phi3),
                  "ihs" = - 0.5 * log(1 + (lambda * Phi3 * y) ^ 2) + log(Phi3),
                  log(Phi3)
                  )
    
    ##............................................
    ##   3.c  Residual of the consumption equation
    ##............................................

    resid <- (Ty - bX2)
    # problem with bc and lambda < 0, for y = 0, resid = -inf and
    # log(dnorm(resid)) = -inf
    resid[y == 0] <- 0

    ##..............................................................
    ##  3.d Computation of the (opposite of) the range of z and of z
    ##  for y = 0
    ##..............................................................
    
    # The minimum value of z is - Inf, except for bc/bc2 with lambda >
    # 0 and for tn (left truncation)
    mzn <- + Inf
    if (dist %in% c("bc", "bc2") && lambda > 0) mzn <- (1 / lambda + bX2) / sigma
    if (dist == "tn") mzn <- bX2 / sigma

    # The maximum value of z is + Inf, except for bc/bc2 with lambda <
    # 0 (right truncation)
    mzx <- - Inf
    if (dist %in% c("bc", "bc2") && lambda < 0) mzx <- (1 / lambda + bX2) / sigma

    # (T(0) - b2x2) / sigma
    mz0 <- + Inf
    if (dist == "bc2") mz0 <- ( - (mu ^ lambda - 1) / lambda + bX2) / sigma
    if (dist == "bc" && lambda > 0)  mz0 <- (1 / lambda + bX2) / sigma
    if (dist == "ln2") mz0 <- (bX2 - log(mu)) / sigma
    if (dist %in% c("n", "tn")) mz0 <- bX2 / sigma
    
    ##.................................................
    ## 3.e Trivariate and bivariate cummulative normals
    ##.................................................

    Pr123A <- PHI3(bX1, mz0, bX3, rho)
    Pr123B <- PHI3(bX1, mzx, bX3, rho)
    if (h1) arg1 <- (bX1 + rho[1] * resid / sigma) / sqrt(1 - rho[1] ^ 2) else arg1 <- Inf
    if (h3) arg3 <- (bX3 + rho[3] * resid / sigma) / sqrt(1 - rho[3] ^ 2) else arg3 <- Inf
    Pr13 <- PHI2(arg1, arg3,
                 (rho[2] - rho[1] * rho[3]) / sqrt(1 - rho[1] ^ 2) / sqrt(1 - rho[3] ^ 2)
                 )
    Pr13$f <- sanitize(Pr13$f, string = c("Pr13", "mhurdle.ln"),
                       verbal = FALSE, replace = TRUE)
    # PI is the correction of the truncature
    PI <- punorm(mzn) - punorm(mzx)
    ##......................................
    ## 3.e Computation of the log-likelihhod   
    ##......................................

    Pplus <- (Pr123A$f - Pr123B$f) / PI
    Numerator <- PI - Pr123A$f + Pr123B$f
    lnL_null <- log(Numerator) - log(PI)
    lnL_one <- log(PI - Numerator) - log(PI)
    lnL_null[y != 0] <- 0
    lnL_pos <-
        - log(sigma) +
        dnorm(resid / sigma, log = TRUE) +
        log(Pr13$f) +
        lnJ - log(PI)
    lnL_pos[y == 0] <- 0    
    lnL <- lnL_null * (y == 0) + lnL_pos * (y != 0)
    if (any(is.na(lnL) | any(is.infinite(lnL)))){
        lnL[is.na(lnL)] <- 0
        warnings("infinite or missing values of lnLi")
    }
    attr(lnL, "parts") <- c(lnLNull = sum(lnL_null[y == 0]),
                            lnLOne  = sum(lnL_one[y != 0]),
                            lnLPos  = sum(lnL_pos[y != 0]))

    #-------------------------------
    # 4. Computation of the gradient
    #-------------------------------

    if (gradient){
        gradi <- c()
        
        ##......................................
        ##   4.a Derivatives respective to beta1
        ##......................................

        if (h1){
            lnL_beta1 <- (y == 0) * ( (- Pr123A$d1 + Pr123B$d1) / Numerator) +
                (y != 0) * ( Pr13$d1 / sqrt(1 - rho[1] ^ 2) / Pr13$f)
            gradi <- cbind(gradi, lnL_beta1 * X1)
        }

        ##......................................
        ##   4.b Derivatives respective to beta2
        ##......................................

        PI2 <-  (dnorm(mzn) - dnorm(mzx) ) / sigma
        lnL_beta2 <- (y == 0) * ( (PI2 - Pr123A$d2 / sigma + Pr123B$d2 / sigma) / Numerator -
                                     PI2 / PI) +
            (y != 0) * ( resid / sigma ^ 2 -
                            ( Pr13$d1 * rho[1] / sigma / sqrt(1 - rho[1] ^ 2) +
                                 Pr13$d2 * rho[3] / sigma / sqrt(1 - rho[3] ^ 2) ) / Pr13$f -
                                     PI2 / PI)
        gradi <- cbind(gradi, lnL_beta2 * X2)
        
        ##......................................
        ##   4.c Derivatives respective to beta3
        ##......................................

        if (h3){
            # derivative of T(Phi3 y) with bX3
            odist <- dist
            if (dist == "ln2" && mu == 0) dist <- "ln"
            Ty_bX3 <- switch(dist,
                             "ln" = mills(bX3),
                             "ln2" = y * phi3 / (mu + y * Phi3),
                             "bc" = exp(lambda * log(y * Phi3)) * mills(bX3),
                             "bc2" = exp((lambda - 1) * log(y * Phi3 + mu)) * phi3 * y,
                             "ihs" = y * phi3 / sqrt( 1 + (lambda * y * Phi3) ^ 2),
                             y * phi3
                             )
            # derivative of lnJ with bX3
            lnJ_bX3 <- switch(dist,
                              "ln" = 0,
                              "ln2" = mills(bX3) - y * phi3 / (mu + y * Phi3),
                              "bc" = lambda * mills(bX3),
                              "bc2" = (lambda - 1) * phi3 * y / (Phi3 * y + mu) + mills(bX3),
                              "ihs" = - phi3 * Phi3 * lambda ^ 2 * y ^ 2 /
                                  (1 + (lambda * y * Phi3) ^ 2) + mills(bX3),
                              mills(bX3)
                              )
            
            lnL_beta3 <- (y == 0) * (- Pr123A$d3 + Pr123B$d3) / Numerator + 
                (y != 0) * (- resid / sigma ^ 2 * Ty_bX3 +
                            ( Pr13$d1 * Ty_bX3 * rho[1] / sigma / sqrt(1 - rho[1] ^ 2) +
                              Pr13$d2 * (1 + Ty_bX3 * rho[3] / sigma) /
                              sqrt(1 - rho[3] ^ 2) ) / Pr13$f +
                            lnJ_bX3
                )
            dist <- odist
            gradi <- cbind(gradi, lnL_beta3 * X3)
        }
        
        ##   4.d  Derivative respective to sigma
        PI_sigma <- 0
        Pr123B_sigma <- 0
        Pr123A_sigma <- ifelse(is.infinite(mz0), 0, Pr123A$d2 * (- mz0 / sigma))
        if (dist == "tn" | (dist %in% c("bc", "bc2") && lambda > 0)) PI_sigma <- PI_sigma + dnorm(mzn) * (- mzn / sigma)
        if (dist %in% c("bc", "bc2") && lambda < 0){
            PI_sigma <- PI_sigma - dnorm(mzx) * (- mzx / sigma)
            Pr123B_sigma <- Pr123B$d2 * (- mzx / sigma)
        }
        lnL_sigma <- (y == 0) * ( (PI_sigma - Pr123A_sigma + Pr123B_sigma) / Numerator - PI_sigma / PI) +
            (y != 0) * (- 1 / sigma + resid ^ 2 / sigma ^ 3 +
                        ( - Pr13$d1 * rho[1] * resid / sigma ^ 2 / sqrt(1 - rho[1] ^ 2) -
                          Pr13$d2 * rho[3] * resid / sigma ^ 2 /
                          sqrt(1 - rho[3] ^ 2))  / Pr13$f -
                        PI_sigma / PI)
        gradi <- cbind(gradi, sigma = lnL_sigma * pbX4  * gradientsd)
        
        ##......................................
        ##   4.e Derivatives respective to beta4
        ##......................................

        if (! is.null(X4)){
            gradi <- cbind(gradi,  lnL_sigma * sd * dnorm(bX4) * X4)
        }
        
        ##....................................
        ##   4.f Derivatives respective to rho
        ##....................................

        if (corr){
            if (h1 & ! h3){
                Pr123A$dr <- cbind(Pr123A$dr, 0, 0)
                Pr123B$dr <- cbind(Pr123B$dr, 0, 0)
            }
            if (! h1 & h3){
                Pr123A$dr <- cbind(0, 0, Pr123A$dr)
                Pr123B$dr <- cbind(0, 0, Pr123B$dr)
            }
            if ( (h1 & h3) & !is.matrix(Pr123A$dr) ) Pr123A$dr <- cbind(0, Pr123A$dr, 0)

            if (length(Pr123B$dr) == 1 && Pr123B$dr == 0) Pr123B$dr <- cbind(0, 0, 0)
            lnL_rho <- c()
            if (h1){
                Drho12 <- (rho[1] * rho[2] - rho[3]) /
                    (1 - rho[1] ^ 2) ^ 1.5  / sqrt(1 - rho[3] ^ 2)
                lnL_rho12 <- (y == 0) * (- Pr123A$dr[, 1] + Pr123B$dr[, 1]) / Numerator  +
                    (y != 0) * ( Pr13$d1 * (resid / sigma + rho[1] / (1 - rho[1] ^ 2) *
                                                (bX1 + rho[1] * resid / sigma) ) /
                                    (1 - rho[1] ^ 2) ^ 0.5 +
                                        Pr13$dr * Drho12) / Pr13$f
                lnL_rho <- cbind(lnL_rho, rho12 = lnL_rho12 * gradientrho[1])
            }
            if (h1 & h3){
                Drho13 <- 1 / sqrt(1 - rho[1] ^ 2) / sqrt(1 - rho[3] ^ 2)
                lnL_rho13 <- (y == 0) * (- Pr123A$dr[, 2] + Pr123B$dr[, 2]) / Numerator  +
                    (y != 0) * ( Pr13$dr * Drho13) / Pr13$f
                lnL_rho <- cbind(lnL_rho, rho13 = lnL_rho13 * gradientrho[2])
            }
            if (h3){
                Drho23 <- (rho[3] * rho[2] - rho[1]) /
                    (1 - rho[3] ^ 2) ^ 1.5  / sqrt(1 - rho[1] ^ 2)
                lnL_rho23 <- (y == 0) * (- Pr123A$dr[, 3] + Pr123B$dr[, 3]) / Numerator  +
                    (y != 0) * ( Pr13$d2 * (resid / sigma + rho[3] / (1 - rho[3] ^ 2) *
                                                (bX3 + rho[3] * resid / sigma) ) /
                                    (1 - rho[3] ^ 2) ^ 0.5 +
                                        Pr13$dr * Drho23) / Pr13$f
                lnL_rho <- cbind(lnL_rho, rho23 = lnL_rho23 * gradientrho[3])
            }
            gradi <- cbind(gradi, lnL_rho)
        }
        
        ##......................................
        ##   4.g Derivative respective to lambda
        ##......................................

        if (dist %in% c("bc", "bc2")){
            lnJ_lambda <- switch(dist,
                                 "bc"  = log2(y) + log(Phi3),
                                 "bc2" = log(Phi3 * y + mu),
                                 "ihs" = - lambda * y ^ 2 * Phi3 ^ 2 / (1 + (lambda * y * Phi3) ^ 2)
                                 )
            Ty_lambda <- (log(Phi3 * y + mu) * (Phi3 * y + mu) ^ lambda - Ty) / lambda
            if (mu == 0) T0_lambda <- ( 1 / lambda ^ 2 * (lambda > 0) + 0 * (lambda < 0))
            else T0_lambda <- (log(mu) * mu ^ lambda - (mu ^ lambda - 1) / lambda) / lambda
            Tymax_lambda <- 0 * (lambda > 0) + (1 / lambda ^ 2) * (lambda < 0)
            PI_lambda <- - sign(lambda) * dnorm( (bX2 + 1 / lambda) / sigma ) / (sigma * lambda ^ 2)
            lnL_lambda <- vector(mode = "numeric", length = length(y))
            lnL_lambda[y == 0] <- ( (PI_lambda - Pr123A$d2 * (- T0_lambda /  sigma) + Pr123B$d2 *
                                         (- Tymax_lambda / sigma)) / Numerator - PI_lambda / PI)[y == 0]
            lnL_lambda[y > 0] <- ( (- resid / sigma ^ 2 +
                                        (Pr13$d1 * rho[1] / sqrt(1 - rho[1] ^ 2) / sigma +
                                             Pr13$d2 * rho[3] / sqrt(1 - rho[3] ^ 2) / sigma) / Pr13$f ) *
                                      Ty_lambda + lnJ_lambda - PI_lambda / PI)[y > 0]
            gradi <- cbind(gradi, tr = lnL_lambda)
        }
        if (dist == "ihs"){
            Ty_lambda <- (y * Phi3) / lambda / sqrt(1 + (lambda * y * Phi3) ^ 2) - Ty / lambda
            lnL_lambda <- vector(mode = "numeric", length = length(y))
            gradi <- cbind(gradi, lnL_lambda)
        }

        ##..................................
        ##   4.h Derivative respective to mu
        ##..................................
        if (dist == "bc2"){
            Ty_mu <- (Phi3 * y + mu) ^ (lambda - 1)
            T0_mu <- mu ^ (lambda - 1)
            lnJ_mu <- (lambda - 1) / (Phi3 * y + mu)
            lnL_mu <- vector(mode = "numeric", length = length(y))
            lnL_mu[y == 0] <- (- Pr123A$d2 * (- T0_mu / sigma) / Numerator)[y == 0]
            lnL_mu[y > 0] <- (( - resid / sigma ^ 2 +
                                (Pr13$d1 * rho[1] / sqrt(1 - rho[1] ^ 2) / sigma +
                                 Pr13$d2 * rho[3] / sqrt(1 - rho[3] ^ 2) / sigma) / Pr13$f
                               ) * Ty_mu + lnJ_mu)[y != 0]
            gradi <- cbind(gradi, mu = lnL_mu * gradientmu)
        }
        if (dist == "ln2"){
            Ty_mu <- 1 / (mu + Phi3 * y)
            T0_mu <- 1 / mu
            lnJ_mu <- - 1 / (mu + Phi3 * y)
            lnL_mu <- vector(mode = "numeric", length = length(y))
            lnL_mu[y == 0] <- ( - Pr123A$d2 * (- T0_mu / sigma ) / Numerator)[y == 0]
            lnL_mu[y != 0] <- (( - resid / sigma ^ 2 +
                                    (Pr13$d1 * rho[1] / sqrt(1 - rho[1] ^ 2) / sigma +
                                         Pr13$d2 * rho[3] / sqrt(1 - rho[3] ^ 2) / sigma) / Pr13$f
            ) * Ty_mu + lnJ_mu)[y != 0]
            gradi <- cbind(gradi, mu = lnL_mu * gradientmu)
        }
        if (any(is.na(gradi))) gradi[is.na(gradi)] <- 0
        # for ln2, when mu = 0, lnL_mu is 0 / 0 = 0 for y = 0, the NaN
        # is correctly replaced by 0
        colnames(gradi) <- names(param)
        attr(lnL, "gradient") <- gradi
    }
    
    #------------------------------------
    # 5. Computation of the fitted values
    #------------------------------------

    if (fitted){
        if (dist %in% c("bc", "bc2")){
            if (length(mzn) == 1) mzn <- rep(mzn, length(y))
            if (length(mzx) == 1) mzx <- rep(mzx, length(y))
            if (length(Phi3) == 1) Phi3 <- rep(Phi3, length(y))
            if (length(Phi1) == 1) Phi1 <- rep(Phi1, length(y))
            
            resid <- function(z, index){
                switch(dist,
                       "ln"  = log2(z) + log(Phi3[index]),
                       "ln2" = log2(z * Phi3[index] + mu),
                       "bc"  = (exp(lambda * log(z * Phi3[index])) - 1) / lambda,
                       "bc2" = (exp(lambda * log(z * Phi3[index] + mu)) - 1) / lambda,
                       "ihs" = log(lambda * z * Phi3[index] + sqrt(1 + (lambda  * z * Phi3[index]) ^ 2)) / lambda,
                       z * Phi3[index]
                       ) - bX2[index]
            }
            lnJ <- function(z, index){
                switch(dist,
                       "ln"  = - log2(z),
                       "ln2" = log(Phi3[index]) - log(mu + Phi3[index] * z),
                       "bc"  = (lambda - 1) * log2(z) + lambda * log(Phi3[index]),
                       "bc2" = (lambda - 1) * log(Phi3[index] * z + mu) + log(Phi3[index]),
                       "ihs" = - 0.5 * log(1 + (lambda * Phi3[index] * z) ^ 2) + log(Phi3[index]),
                       log(Phi3[index])
                       )
            }
            
            arg1 <- function(z, index) if (h1) (bX1[index] + rho[1] * resid(z, index) / sigma) / sqrt(1 - rho[1] ^ 2) else Inf
            arg3 <- function(z, index) if (h3) (bX3[index] + rho[3] * resid(z, index) / sigma) / sqrt(1 - rho[3] ^ 2) else Inf
            Pr13 <- function(z, index) PHI2(arg1(z, index), arg3(z, index),
                                            (rho[2] - rho[1] * rho[3]) / sqrt(1 - rho[1] ^ 2) / sqrt(1 - rho[3] ^ 2))$f
            # Pbnorm plante et pas PHI2        
            PI <- function(index) punorm(mzn[index]) - punorm(mzx[index])
            fplus <- function(z, index){
                x <- z * exp(- log(sigma) + dnorm(resid(z, index) / sigma, log = TRUE) + log(Pr13(z, index)) + lnJ(z, index) - log(PI(index)))
                x
            }        
            if (dist %in% c("bc", "bc2") && lambda < 0) maxint <- max(y) * 3
            else maxint <- +Inf
            E <- sapply(seq_len(length(y)), function(i) integrate(function(x) fplus(x, i), 0, maxint)$value)
        }
        if (dist %in% c("ln", "ln2")){
            if (h1) arg1 <- bX1 + rho[1] * sigma else arg1 <- + Inf
            if (h3) arg3 <- bX3 + rho[3] * sigma else arg3 <- + Inf
            E <- exp(bX2 + 0.5 * sigma ^ 2) / Phi3 * ptnorm(arg1, mz0 + sigma, arg3, rho)
            if (dist == "ln2") E <- E - mu * ptnorm(bX1, mz0, bX3, rho)/ Phi3
        }
        if (dist %in% c("n", "tn")){
            phi2 <- dnorm(bX2 / sigma)
            phi1 <- dnorm(bX1)
            phi3 <- dnorm(bX3)
            Pr13 <- pbnorm((bX1 - rho[1] * bX2 / sigma) / sqrt(1 - rho[1] ^ 2),
                           (bX3 - rho[3] * bX2 / sigma) / sqrt(1 - rho[3] ^ 2),
                           (rho[2] - rho[1] * rho[3]) / sqrt( (1 - rho[1] ^ 2) * (1 - rho[3] ^ 2)))
            Pr23 <- pbnorm((bX2 / sigma - rho[1] * bX1) / sqrt(1 - rho[1] ^ 2),
                           (bX3         - rho[2] * bX1) / sqrt(1 - rho[2] ^ 2),
                           (rho[3] - rho[1] * rho[2]) / sqrt( (1 - rho[1] ^ 2) * (1 - rho[2] ^ 2)))
            
            Pr12 <- pbnorm((bX1         - rho[2] * bX3) / sqrt(1 - rho[2] ^ 2),
                           (bX2 / sigma - rho[3] * bX3) / sqrt(1 - rho[3] ^ 2),
                           (rho[1] - rho[2] * rho[3]) / sqrt(1 - rho[2] ^ 2) / sqrt(1 - rho[3] ^ 2)
                           )
            E <- bX2 / Phi3 + sigma / (Phi3 * Pplus) * (phi2 * Pr13 +
                                                                rho[1] * phi1 * Pr23 / sqrt( (1 - rho[1] ^ 2) * (1 - rho[3] ^ 2)) +
                                                                    rho[3] * phi3 * Pr12 / sqrt( (1 - rho[2] ^ 2) * (1 - rho[3] ^ 2)))

            E <- bX2 / Phi3 + sigma / (Phi3 * Pplus) * (phi2 * Pr13 + rho[1] * phi1 * Pr23 + rho[3] * phi3 * Pr12)

        }
        if (! is.null(attr(y, "geomean"))){
            E <- E * attr(y, "geomean")
        }
        attr(lnL, "fitted") <- cbind(E = E,
                                     Ep = E / Pplus,
                                     p = Pplus)
    }
    lnL
}
