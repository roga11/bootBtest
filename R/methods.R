
# -----------------------------------------------------------------------------
#' @title Compute Bootstrap P-Value
#' 
#' @description This function computes the Bootstrap p-value consistent with a left-tailed test, right-tailed test, two-tailed symmetric test, and two-tailed non-symmetric test. See ch. 4 of Davidson & MacKinnon (2004) or pages 11-12 of Davidson (2013).
#'
#' @param tau test statistic from observed data
#' @param tau_B vector of bootstrap test statistics from bootstrap samples 
#' @param type string specifying type of test. Options are: "geq" for right-tailed test, "leq" for 
#' left-tailed test, "two-tailed-s" for two-tailed symmetric test and "two-tailed-ns" for two-tailed non-symmetric test.
#' 
#' @return Bootstrap p-value
#' 
#' @references Russell Davidson & James G. MacKinnon (2004), Econometric Theory and Methods, New York, Oxford University Press.
#' @references Russell Davidson (2013), The Bootstrap in Econometrics, CEA Conference, May 2013
#' 
#' @export
boot_pval <- function(tau, tau_B, type = 'geq'){
  B <- length(tau_B)
  if (type == "absolute" || type == "geq") {
    pval = sum(tau_B>tau)/B
  }else if (type == "leq") {
    pval = sum(tau_B<tau)/B
  }else if (type == "two-tailed-s") {
    pval = sum(abs(tau_B)>abs(tau))/B
  }else if (type == "two-tailed-ns") {
    pval = 2*min(c(sum(tau_B<=tau)/B, sum(tau_B>tau)/B))
  } else{
    print("type must be one of the following: geq, leq, two-tailed or absolute")
  }
  return(pval)
}

# -----------------------------------------------------------------------------
#' @title Bootstrap Pretest for number of simulations B
#' 
#' @description This function performs the pretest procedure described in Davidson & MacKinnon (2000) to determine the appropriate number of simulations needed to minize loss of power.
#'
#' 
#' @references Russell Davidson & James G. MacKinnon (2000), Bootstrap tests: how many bootstraps?, Econometric Reviews, 19:1, 55-68.
#' 
#' @export
bootBtest <- function(tau, tau_B, type = 'geq', B_max = 12799, alpha = 0.05, beta = 0.001){
  B <- length(tau_B)
  # Compute Bootstrap p-value
  pval <- boot_pval(tau, tau_B, type)
  # check if B is enough
  if (pval<=alpha){
    if (B<(10/alpha)){
      check <- 1 - pbinom(pval*B,B, alpha, lower.tail = FALSE) 
    }else{
      check <- 1 - pnorm(pval*B, B*alpha + 0.5, sqrt(B*alpha*(1-alpha)), lower.tail = FALSE)
    }
  }else if (pval>alpha){
    if (B<(10/alpha)){
      check <- pbinom(pval*B,B,alpha,lower.tail = FALSE)
    }else{
      check <- pnorm(pval*B, B*alpha - 0.5, sqrt(B*alpha*(1-alpha)), lower.tail = FALSE)
    }
  }
  val <- (check<=beta)
  if (val==TRUE){
    # stop 
    B_next = 0
    msg <- paste0("Test is rejected - Sufficient number of simulations")
  }else{
    # continue
    if ((2*B+1)>B_max){
      B_next <- B_max - B
    }else{
      B_next <- B + 1
    }
    msg <- paste0("Fail to reject test - Consider using an additional ", B_next, " number of simulations")
    # if (B>B_max){
    #   val <- FALSE
    #   msg <- paste0("Reached max number of simulations")
    # }
  } 
  test_output <- list()
  test_output$B <- B
  test_output$boot_pval <- pval
  test_output$beta <- beta
  test_output$alpha <- alpha
  test_output$pval <- check
  test_output$val <- val
  test_output$message <- msg
  test_output$B_next <- B_next
  return(test_output)
}



# -----------------------------------------------------------------------------
#' @title Bootstrap Pretest
#' 
#' @description This function performs the pretest procedure described in Davidson & MacKinnon (2000) to determine the appropriate number of simulations needed to minimize loss of power.
#'
#' 
#' @references Russell Davidson & James G. MacKinnon (2000), Bootstrap tests: how many bootstraps?, Econometric Reviews, 19:1, 55-68.
#' 
#' @export
bootandBtest <- function(tau, fn, Bname, type = 'geq', B_min = 99, B_max = 12799, 
                         alpha = 0.05, beta = 0.001, seed = NULL, ...){
  # Do some checks
  jc <- NULL
  if (!is.function(fn) || is.null(fn)) {
    stop("'fn' has to be a R function.")
  }
  # set seed (if specified)
  if (is.null(seed)==FALSE){
    set.seed(seed) 
  }
  # initialize number of bootstraps
  B <- B_min
  B_prime <- B_min
  B_tmp <- B
  # Get random vector tau_B
  tau_B <- matrix(0,0,1)
  stop <- FALSE
  while (stop==FALSE){
    tau_B_tmp <- fn(assign(Bname, B_tmp), ...)
    tau_B <- c(tau_B, tau_B_tmp)
    boottest <- bootBtest(tau, tau_B, type = type, B_max = B_max, alpha = alpha, beta = beta)
    # Compute Bootstrap p-value
    pval <- boottest$boot_pval
    val <- boottest$val
    if (val==TRUE){
      # stop 
      stop <- TRUE
    }else{
      # continue
      B <- 2*B_prime + 1
      B_tmp <- B_prime + 1
      B_prime <- B
      if (B>B_max){
        stop <- TRUE
      }
    } 
  }
  return(c(length(tau_B), pval))
}

# -----------------------------------------------------------------------------
#' @title Bootstrap Pretest
#' 
#' @description This function performs the pretest procedure described in Davidson & MacKinnon (2000) to determine the appropriate number of simulations needed to minimize loss of power.
#'
#' 
#' @references Russell Davidson & James G. MacKinnon (2000), Bootstrap tests: how many bootstraps?, Econometric Reviews, 19:1, 55-68.
#' 
#' @export
boot_and_pretest <- function(tau, gamma_null, std, n, type = 'geq', alpha = 0.05, 
                         B_min = 99, B_max = 12799, beta = 0.001, seed = NULL){
  # initialize number of bootstraps
  B <- B_min
  B_prime <- B_min
  B_tmp <- B
  # Get random vector tau_N
  if (is.null(seed)==FALSE){
    set.seed(seed) 
  }
  tau_B <- matrix(0,0,1)
  stop <- FALSE
  while (stop==FALSE){
    tau_B_tmp <- bootNullDist_ex(gamma_null, std, n, B_tmp)
    tau_B <- c(tau_B, tau_B_tmp)
    boottest <- bootBtest(tau, tau_B, type = type, B_max = B_max, alpha = alpha, beta = beta)
    # Compute Bootstrap p-value
    pval <- boottest$boot_pval
    val <- boottest$val
    if (val==TRUE){
      # stop 
      stop <- TRUE
    }else{
      # continue
      B <- 2*B_prime + 1
      B_tmp <- B_prime + 1
      B_prime <- B
      if (B>B_max){
        stop <- TRUE
      }
    } 
  }
  return(c(length(tau_B), pval))
}


# -----------------------------------------------------------------------------
#' @title Test statistic
#' 
#' @description This function computes test statistics for 'T-test' and 'F-test'. 
#'
#' 
#' @export
compute_stat <- function(mdl, null_vec){
  out     <- summary(mdl)
  coef    <- as.matrix(out$coefficients[,1])
  stderr  <- as.matrix(out$coefficients[,2])
  # run checks
  if (length(coef)!=length(null_vec)){
    stop("the vector 'null_vec' must be of same length as parameter vector. Use vallue 'NA' for arameters not being tested.")
  }
  #std     <- sqrt(sum(mdl$residuals^2)/(n-1)) 
  if (sum(is.na(null_vec)==FALSE)==1){
    theta_hat <- coef[is.na(null_vec)==FALSE,]
    theta_null <- null_vec[is.na(null_vec)==FALSE]
    theta_std_err <- stderr[is.na(null_vec)==FALSE,]
    tau <- (theta_hat - theta_null)/theta_std_err 
  }else{
    stop("T-test only allows one value in 'null_vec' to be different from NA. For test of multiple parameters use 'F-test'.")
  }
  names(tau) <- NULL
  return(tau)
}
