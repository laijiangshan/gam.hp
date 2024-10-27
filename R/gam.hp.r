#' Hierarchical Partitioning of Adjusted R2 and Explained Deviance for Generalized Additive Models

#' @param  mod  Fitted "gam" model objects.
#' @param  type The type of R-square of gam, either "dev" or "adjR2", in which "dev" is explained deviance and "adjR2" is adjusted R-square, the default is "dev".
#' @param  commonality Logical; If TRUE, the result of commonality analysis (2^N-1 fractions for N predictors) is shown, the default is FALSE. 

#' @details This function conducts hierarchical partitioning to calculate the individual contributions of each predictor towards total adjusted R2 and explained deviance for Generalized Additive Models. The adjusted R2 and explained deviance are is the output of summary.gam()in mgcv package.

#' @return \item{dev}{The R2 for the full model.}
#' @return \item{hierarchical.partitioning}{A matrix containing individual effects and percentage of individual effects towards total adjusted R2 and explained deviance for each predictor.}

#' @author {Jiangshan Lai} \email{lai@njfu.edu.cn}



#' @references
#' \itemize{
#' \item Lai J.,Tang J., Li T., Zhang A.,Mao L.(2024)Evaluating the relative importance of predictors in Generalized Additive Models using the gam.hp R package.Plant Diversity,46(4):542-546<DOI:10.1016/j.pld.2024.06.002>
#' \item Lai J.,Zhu W., Cui D.,Mao L.(2023)Extension of the glmm.hp package to Zero-Inflated generalized linear mixed models and multiple regression.Journal of Plant Ecology,16(6):rtad038<DOI:10.1093/jpe/rtad038>
#' \item Lai J.,Zou Y., Zhang S.,Zhang X.,Mao L.(2022)glmm.hp: an R package for computing individual effect of predictors in generalized linear mixed models.Journal of Plant Ecology,15(6):1302-1307<DOI:10.1093/jpe/rtac096>
#' \item Lai J.,Zou Y., Zhang J.,Peres-Neto P.(2022) Generalizing hierarchical and variation partitioning in multiple regression and canonical analyses using the rdacca.hp R package.Methods in Ecology and Evolution,13(4):782-788<DOI:10.1111/2041-210X.13800>
#' \item Chevan, A. & Sutherland, M. (1991). Hierarchical partitioning. American Statistician, 45, 90-96. doi:10.1080/00031305.1991.10475776
#' \item Nimon, K., Oswald, F.L. & Roberts, J.K. (2013). Yhat: Interpreting regression effects. R package version 2.0.0.
#' }


#'@export
#'@examples
#'library(mgcv)
#'mod1 <- gam(Sepal.Length ~ s(Petal.Length) + s(Petal.Width) + Sepal.Width,data = iris)
#'summary(mod1)
#'gam.hp(mod1)
#'gam.hp(mod1,type="adjR2")
#'gam.hp(mod1,commonality=TRUE)

#'@importFrom(foreach,"%dopar%")

gam.hp <- function(mod, type = "dev", commonality = FALSE, parallel = FALSE, cores = NA) {
  # initial checks
  if (!inherits(mod, c("glm", "lm", "gam"))) stop("gam.hp only supports gam objects at the moment")

  Formu <- strsplit(as.character(mod$call$formula)[3], "")[[1]]
  if ("*" %in% Formu)
    stop("Please put the interaction term as a new variable (i.e. link variables by colon(:)) and  avoid the asterisk (*) and colon(:) in the original model")
  #if GAM，not use orgianl ivname
  #if('gam' %in% class(mod)) ivname <- str_split(str_split(as.character(mod$call$formula)[2] ,'~')[[1]][2] ,'[/+]')[[1]] else ivname <- strsplit(as.character(mod$call$formula)[3],"[+]")[[1]] 
  #try delete call in mod$call$formula, if some question
  #if('gam' %in% class(mod)) ivname <- str_split(as.character(mod$formula)[3] ,'[/+]')[[1]] else ivname <- strsplit(as.character(mod$call$formula)[3],"[/+]")[[1]] 
  ivname <- strsplit(as.character(mod$formula)[3], "[/+]")[[1]]

  iv.name <- ivname
  nvar <- length(iv.name)
  if (nvar < 2)
    stop("Analysis not conducted. Insufficient number of predictors.")
  totalN <- 2^nvar - 1
  binarymx <- matrix(0, nvar, totalN)
  print("Calculating binary matrix.")
  pb <- txtProgressBar(max = 100, style = 3)
  for (i in 1:totalN) {
    setTxtProgressBar(pb, i / totalN * 100)
    binarymx <- creatbin(i, binarymx)
  }
  close(pb)

  if (type == "adjR2") outr2  <- summary(mod)$r.sq
  if (type == "dev") outr2  <- summary(mod)$dev.expl

  r2type <- row.names(outr2)
  nr2type <- length(r2type)
  if (nr2type == 0) {
    nr2type <- 1
    if (commonality) {
      r2type <- "commonality.analysis"
    } else {
      r2type <- "hierarchical.partitioning"
    }
  }

  # Update null model
  if (inherits(mod, c("gam", "glm", "lm"))) {
    dat <- na.omit(eval(mod$model))
    if (!inherits(dat, "data.frame")) {
      stop("Please change the name of data object in the original gam analysis then try again.")
    }
    #判断一下是不是GAM，不是就用原来的to_del
    if ("gam" %in% class(mod)) {
      to_del <- paste(as.character(mod$formula)[2], "~", "1")
    } else {
      to_del <- paste(as.character(mod$call$formula)[2], "~", "1")
    }

    mod_null <-  stats::update(object = mod, formula. = to_del, data = dat)
  }

  outputList  <- list()
  outputList[[1]] <- outr2
  # Update models and return the results
  for (k in 1:nr2type) {
    print("Calculating GAM models.")
    commonM <- matrix(nrow = totalN, ncol = 3)
    if (parallel) {
      if (is.na(cores)) {
        cores <- detectCores() - 1
      }
      cluster <- makeSOCKcluster(cores)
      registerDoSNOW(cluster)
      pb <- txtProgressBar(max = 100, style = 3)
      progress <- function(n) setTxtProgressBar(pb, n  / totalN * 100)
      opts <- list(progress = progress)
      results <- foreach(i = 1:totalN,
                         .options.snow = opts,
                         .packages = c("mgcv", "gam.hp")) %dopar% {
        tmp.name <- iv.name[as.logical(binarymx[, i])]
        to_add <- paste("~", paste(tmp.name, collapse = " + "), sep = " ")
        modnew  <- stats::update(object = mod_null, data = dat, to_add)
        if (type == "dev") {
          out <- summary(modnew)$dev.expl
        }
        if (type == "adjR2") {
          out <- summary(modnew)$r.sq
        }
        rm(tmp.name, to_add, modnew)
        gc()
        out
      }
      close(pb)
      commonM[, 2] <- unlist(results)
      stopCluster(cl = cluster)
    } else {
      pb <- txtProgressBar(max = 100, style = 3)
      for (i in 1:totalN) {
        tmp.name <- iv.name[as.logical(binarymx[, i])]
        to_add <- paste("~", paste(tmp.name, collapse = " + "), sep = " ")
        modnew  <- stats::update(object = mod_null, data = dat, to_add)
        setTxtProgressBar(pb, i / totalN * 100)
        if (type == "dev") commonM[i, 2]  <- summary(modnew)$dev.expl
        if (type == "adjR2") commonM[i, 2]  <- summary(modnew)$r.sq
      }
      close(pb)
    }
  }

  commonlist <- vector("list", totalN)

  seqID <- vector()
  for (i in 1:nvar) {
    seqID[i] = 2^(i-1)
  }

  for (i in 1:totalN) {
    bit <- binarymx[1, i]
    if (bit == 1)
      ivname <- c(0, -seqID[1])
    else ivname <- seqID[1]
    for (j in 2:nvar) {
      bit <- binarymx[j, i]
      if (bit == 1) {
        alist <- ivname
        blist <- genList(ivname, -seqID[j])
        ivname <- c(alist, blist)
      } else {
        ivname <- genList(ivname, seqID[j])
      }
    }
    ivname <- ivname * -1
    commonlist[[i]] <- ivname
  }

  for (i in 1:totalN) {
    r2list <- unlist(commonlist[i])
    numlist  <-  length(r2list)
    ccsum  <-  0
    for (j in 1:numlist) {
      indexs  <-  r2list[[j]]
      indexu  <-  abs(indexs)
      if (indexu != 0) {
        ccvalue  <-  commonM[indexu, 2]
        if (indexs < 0)
          ccvalue  <-  ccvalue * -1
        ccsum  <-  ccsum + ccvalue
      }
    }
    commonM[i, 3]  <-  ccsum
  }

  orderList <- vector("list", totalN)
  index  <-  0
  for (i in 1:nvar) {
    for (j in 1:totalN) {
      nbits  <-  sum(binarymx[, j])
      if (nbits == i) {
        index  <-  index + 1
        commonM[index, 1] <- j
      }
    }
  }

  outputcommonM <- matrix(nrow = totalN + 1, ncol = 2)
  totalRSquare <- sum(commonM[, 3])
  for (i in 1:totalN) {
    outputcommonM[i, 1] <- round(commonM[commonM[i,
                                                 1], 3], digits = 4)
    outputcommonM[i, 2] <- round((commonM[commonM[i,
                                                  1], 3]/totalRSquare) * 100, digits = 2)
  }
  outputcommonM[totalN + 1, 1] <- round(totalRSquare,
                                        digits = 4)
  outputcommonM[totalN + 1, 2] <- round(100, digits = 4)
  rowNames  <-  NULL
  for (i in 1:totalN) {
    ii  <-  commonM[i, 1]
    nbits  <-  sum(binarymx[, ii])
    cbits  <-  0
    if (nbits == 1)
      rowName  <-  "Unique to "
    else rowName  <-  "Common to "
    for (j in 1:nvar) {
      if (binarymx[j, ii] == 1) {
        if (nbits == 1)
          rowName  <-  paste(rowName, iv.name[j], sep = "")
        else {
          cbits  <-  cbits + 1
          if (cbits == nbits) {
            rowName  <-  paste(rowName, "and ", sep = "")
            rowName  <-  paste(rowName, iv.name[j], sep = "")
          }
          else {
            rowName  <-  paste(rowName, iv.name[j], sep = "")
            rowName  <-  paste(rowName, ", ", sep = "")
          }
        }
      }
    }
    rowNames  <-  c(rowNames, rowName)
  }
  rowNames  <-  c(rowNames, "Total")
  rowNames <- format.default(rowNames, justify = "left")
  colNames <- format.default(c("Fractions", " % Total"),
                             justify = "right")
  dimnames(outputcommonM) <- list(rowNames, colNames)

  VariableImportance <- matrix(nrow = nvar, ncol = 4)
  # VariableImportance <- matrix(nrow = nvar, ncol = 2)
  for (i in 1:nvar) {
    VariableImportance[i, 3] <-  round(sum(binarymx[i, ] * (commonM[, 3] / apply(binarymx, 2, sum))), digits = 4)
    #VariableImportance[i, 1] <-  round(sum(binarymx[i, ] * (commonM[,3]/apply(binarymx,2,sum))), digits = 4)
  }

  VariableImportance[, 1] <- outputcommonM[1:nvar, 1]
  VariableImportance[, 2] <- VariableImportance[, 3] - VariableImportance[, 1]

  total <- round(sum(VariableImportance[, 3]), digits = 4)
  #total=round(sum(VariableImportance[,1]),digits = 4)
  VariableImportance[, 4] <- round(100 * VariableImportance[, 3] / total, 2)
  #VariableImportance[, 2] <- round(100*VariableImportance[, 1]/total,2)
  #dimnames(VariableImportance) <- list(iv.name, c("Individual","I.perc(%)"))
  dimnames(VariableImportance) <- list(iv.name, c("Unique", "Average.share", "Individual", "I.perc(%)"))

  if (commonality) {
    return(outputcommonM)
  } else {
    return(VariableImportance)
  }

  outputList  <- c(outr2, do.call(par.result))

  if (type == "adjR2") {
    names(outputList) <- c("adjusted.R2", r2type)
  }
  if (type == "dev") {
    names(outputList) <- c("Explained.deviance", r2type)
  }
  #if(inherits(mod, "lm")&!inherits(mod, "glm")){names(outputList) <- c("Total.R2",r2type)}
  outputList$variables <- iv.name
  if (commonality) {
    outputList$type <- "commonality.analysis"
  }
  if (!commonality) {
    outputList$type <- "hierarchical.partitioning"
  }
  class(outputList) <- "gamhp" # Class definition
  outputList
}