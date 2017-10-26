#' @title Wild bootstrap improvements for inference with small cluser size (OLS)
#'
#' @description The function is used to get accurate and precise inference from small cluster sizes. Through  (restricted) wild clustered bootstrap-t to get p-values for OLS models, the function considerably reduces the probability of over rejection of the null hypotheses.The default distribution used to reconstruct the residuals is Rademacher. When cluster size is smaller than 12, the algorithm will suggest to use "six - point" distribution instead.
#'
#' @param model A linear model estimated using \code{lm}. If the mdoel is  nonlinear, please use \code{swboot} instead.
#' @param cluster A formula to specify the variable of cluster.
#' @param beta A key variable of interest. To denote the key varibale to focus on, one could save a amount of time because the algorithm will only report the outcome of key variable. The defalt value is 'All' which will report all the results of every coefficents.
#'
#' @param type The type of distributions for wild boostrap. Default value is Rademacher. If the number of clusters is less than 12, please use 'six' instead (type='six')
#' @param confid The size of test
#' @param R The number of boostrap replicates. The default number is 200.
#' @param seed A single value passed to \code{set.seed}.
#'
#'
#' @details
#'  linear model
#' impose H0 where H0: beta = 0
#'
#'
#' @author Jiawei Fu \email{jiawei.fu@duke.edu}
#'
#'@examples
#' \dontrun{
#' data(wv6_equ)  # input data
#' wv6_equ <- as.data.frame(wv6_equ)
#' model_equ <- lm(income_equ ~ income + age + gender, data = wv6_equ) # linear model
#' wildboott(wv6_equ, ~country, R = 250) # find original p value is over estimated
#'}
#'
#'@references Cameron, A. Colin, Jonah B. Gelbach, and Douglas L. Miller. 2008. "Bootstrap-Based Improvements for Inference with Clustered Errors" \emph{The Review of Economics and Statistics} 90(3): 414-427.
#'
#'@references Cameron, A. Colin and Douglas L. Miller. 2015. "A Practitioner’s Guide to Cluster-Robust Inference" \emph{The Journal of Human Resources} Vol. 50 No. 2 317-372
#'
#'@references Webb, Matthew D. 2013. "Reworking wild bootstrap based inference for clustered errors" \emph{Queen's Economics Department Working Paper}
#'
#'@importFrom utils txtProgressBar
#'@importFrom utils setTxtProgressBar
#'@importFrom compiler cmpfun
#'@importFrom sandwich estfun
#'@import stats
#'
#'
#'@export

wildboott <- function(model, cluster, beta = "All", type = "Rademacher", confid = 0.05, R = 200, seed = NULL) {

  cl <- match.call()

  if (!is.null(seed))
    set.seed(seed)



  clus_vcov <- function(x, cluster) {
    ef <- estfun(x)
    k <- NCOL(ef)
    n <- NROW(ef)
    rval <- matrix(0, nrow = k, ncol = k, dimnames = list(colnames(ef), colnames(ef)))

    ### try my func begin
    mf <- cmpfun(function(y) {
      colSums(ef[which(cluster[, 1] == y), , drop = FALSE])
    })
    ef <- sapply(sort(unique(cluster[, 1])), mf)
    ef <- t(ef)

    ### try my func end
    g <- length(unique(cluster[[1]]))
    adj <- g/(g - 1L)
    rval <- rval + adj * crossprod(ef)/n
    rval <- (n - 1L)/(n - k) * rval
    # bread
    bread <- summary.lm(x)$cov.unscaled * as.vector(sum(summary.lm(x)$df[1:2]))
    # sandwich
    return(1/n * (bread %*% rval %*% bread))
  }

  super <- cmpfun(clus_vcov)
  myfun <- super  # see super vcovCL.R

  if (!all(class(model) == "lm"))
    stop("wild bootstrap improvement is only avaliable to linear model. Trys scboott for nonlinear models.")



  if (!is.null(model$na.action))
    class(model$na.action) <- "omit"

  ##### use regular expression

  if (!as.character(cl[[3]][[2]]) %in% names(model$coefficients)) {
    cluster_tmp <- expand.model.frame(model, cluster, na.expand = FALSE)
    df <- na.omit(cluster_tmp)
    cluster <- model.frame(cluster, df, na.action = na.pass)
    model$call[[3]] <- quote(df)
    model_use <- eval(model$call)  #refit the model with new df
  } else {
    df <- na.omit(model$model)
    cluster <- model.frame(cluster, df, na.action = na.pass)
    model_use <- model
  }


  # ***check if there is NA in cluster, should I refit original model

  # Factors in our clustering variables can potentially cause problems Blunt fix is to force
  # conversion to characters
  i <- !sapply(cluster, is.numeric)
  cluster[i] <- lapply(cluster[i], as.character)

  # step 1, get wald test in the original model
  ori_clucv <- myfun(model_use, cluster = cluster)

  if (R < 2)
    stop("R should be larger than 1")

  if (length(unique(cluster[, 1])) < 12)
    cat("Considering the cluster size is smaller than 12, using six - point instead of Rademacher distribution",
        "\n")


  indvar <- attr(model_use$terms, "term.labels")  # independent variables

  if (!beta == "All") {
    focus <- beta
  } else {
    focus <- indvar
  }

  dat_loc <- which(names(model_use$call) == "data")
  boot_call <- model_use$call[c(-1, -dat_loc)]
  boot_call$formula <- update.formula(formula(model_use), y_boot ~ .)
  estimator <- cmpfun(eval(model_use$call[[1]]))  #cmpfun!!!!!!!!




  # some preparation
  restrict_df <- df
  output_cont <- rep(NA, 5 * length(focus))  #put the output


  for (v in 1:length(focus)) {
    if (length(focus) == 1) {
      cat("Bootstrap the variable:", beta, "\n")
    } else {
      cat("Bootstrap the variable:", names(coef(model_use))[v + 1], "\n")
    }
    pos <- which(names(coef(model_use)) == focus[v])  # position of the focus beta
    beta_focus <- coef(model_use)[[pos]]
    waldsta_ori <- beta_focus/sqrt(ori_clucv[pos, pos])

    # get restricted estimator and residual

    pos1 <- which(names(restrict_df) == focus[v])
    restrict_df[pos1] <- 0  # as if beta=0
    restrict_mc <- model_use$call  #restrict model call
    restrict_mc[[3]] <- quote(restrict_df)
    restrict_m <- eval(restrict_mc)  #restrict mdoel

    ### ****** deal with the NA

    # step 2, wild bootstrap

    # start for loop

    ### put the wald in each boot

    beta_focus_new <- rep(NA, R)
    new_clucv <- rep(NA, R)
    fit_rm <- fitted(restrict_m)
    res_rm <- residuals(restrict_m)

    pb <- txtProgressBar(min = 0, max = R, initial = 0, style = 3)


    ### add 6- point

    for (i in 1:R) {
      setTxtProgressBar(pb, value = i)
      if (type == "Rademacher") {
        wild <- c(1, -1)[rbinom(unique(cluster[, 1]), size = 1, prob = 0.5) + 1][match(cluster[,
                                                                                               1], unique(cluster[, 1]))]
      } else {
        wild <- sample(c(-sqrt(1.5), -1, -sqrt(0.5), sqrt(0.5), 1, sqrt(1.5)), size = length(unique(cluster[,
                                                                                                            1])), replace = TRUE, prob = rep(1/6, 6))[match(cluster[, 1], unique(cluster[,
                                                                                                                                                                                         1]))]
      }

      y_boot <- fit_rm + res_rm * wild

      model_boot <- estimator(boot_call, data = df)

      beta_focus_new[i] <- coef(model_boot)[[pos]]
      new_clucv[i] <- myfun(model_boot, cluster = cluster)[pos, pos]

    }

    close(pb)

    rec <- beta_focus_new/sqrt(new_clucv)
    pvalue <- 1 - (sum(abs(waldsta_ori) > abs(rec))/R)
    pval <- 2 * min((abs(waldsta_ori) > abs(rec))/R, (abs(waldsta_ori) < abs(rec))/R)
    quant <- (quantile(rec, c(confid/2, 1 - confid/2)))


    output_cont[5 * v - 4] <- beta_focus
    output_cont[5 * v - 3] <- waldsta_ori  # waldsta_ori
    output_cont[5 * v - 2] <- quant[[1]]
    output_cont[5 * v - 1] <- quant[[2]]
    output_cont[5 * v] <- pvalue


  }


  output <- matrix(output_cont, nrow = length(focus), ncol = 5, byrow = TRUE)
  rownames(output) <- focus
  colnames(output) <- c("Estimate", "Ori. Wald", "2.5% btWald", " 97.5% btWald", "P-value")
  cat("The score based wild bootstrap outcomes:", "\n")
  printCoefmat(output, P.values = TRUE, has.Pvalue = TRUE)

  # for extract
  output2 <- list()
  output2$estimate <- output[, 1]
  output2$wald0 <- output[, 2]
  output2$wald1 <- output[, 3]
  output2$wald2 <- output[, 4]
  output2$pvalue <- output[, 5]
  invisible(output2)
}




