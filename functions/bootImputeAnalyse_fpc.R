bootImputeAnalyse_fpc <- 
function (FPC, imps, analysisfun, ..., quiet = FALSE) 
{
  nBoot <- attributes(imps)$nBoot
  nImp <- attributes(imps)$nImp
  firstResult <- analysisfun(imps[[1]], ...)
  nParms <- length(firstResult)
  ests <- array(0, dim = c(nBoot, nImp, nParms))
  count <- 1
  for (b in 1:nBoot) {
    for (m in 1:nImp) {
      ests[b, m, ] <- analysisfun(imps[[count]], ...)
      count <- count + 1
    }
  }
  est <- array(0, dim = nParms)
  var <- array(0, dim = nParms)
  ci <- array(0, dim = c(nParms, 2))
  df <- array(0, dim = nParms)
## START RA EDIT SECTION
  var_ml <- array(0, dim = nParms)
  var_bd <- array(0, dim = nParms)
## END RA EDIT SECTION
  for (i in 1:nParms) {
    SSW <- sum((ests[, , i] - rowMeans(ests[, , i]))^2)
    SSB <- nImp * sum((rowMeans(ests[, , i]) - mean(ests[, 
                                                         , i]))^2)
    MSW <- SSW/(nBoot * (nImp - 1))
    MSB <- SSB/(nBoot - 1)
    resVar <- MSW
    randIntVar <- (MSB - MSW)/nImp
    if (randIntVar <= 0) {
      warning(paste("Parameter ", i, " has an estimated between bootstrap variance of zero. You should re-run with a larger nBoot value.", 
                    sep = ""))
      randIntVar <- 0
      resVar <- (SSW + SSB)/(nBoot * nImp - 1)
    }
    est[i] <- mean(ests[, , i])
## START RA EDIT SECTION
    var[i] <- (1 + 1/nBoot) * FPC* randIntVar + resVar/(nBoot *nImp)
## END RA EDIT SECTION
    df[i] <- (var[i]^2)/((((nBoot + 1)/(nBoot * nImp))^2 * 
                            MSB^2/(nBoot - 1)) + MSW^2/(nBoot * nImp^2 * (nImp - 
                                                                            1)))
    df[i] <- max(3, df[i])
    ci[i, ] <- c(est[i] - stats::qt(0.975, df[i]) * var[i]^0.5, 
                 est[i] + stats::qt(0.975, df[i]) * var[i]^0.5)
## START RA EDIT SECTION
    var_ml[i] <- randIntVar
    var_bd[i] <- resVar
## END RA EDIT SECTION
  }
  if (quiet == FALSE) {
    resTable <- array(0, dim = c(nParms, 5))
    resTable[, 1] <- est
    resTable[, 2] <- var^0.5
    resTable[, 3] <- ci[, 1]
    resTable[, 4] <- ci[, 2]
    resTable[, 5] <- 2 * stats::pt(abs(est/var^0.5), df = df, 
                                   lower.tail = FALSE)
    colnames(resTable) <- c("Estimate", "Std. error", "95% CI lower", 
                            "95% CI upper", "p")
    rownames(resTable) <- names(firstResult)
    print(resTable)
  }
  list(ests = est, var = var, ci = ci, df = df, var_ml = var_ml, var_bd = var_bd)
}