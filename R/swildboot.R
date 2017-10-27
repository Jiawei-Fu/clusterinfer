#' @title Score based Wild bootstrap improvements for inference with small cluser size (MLE)
#'
#' @description The function is used to get accurate and precise inference from small cluster sizes. Through  (restricted) wild clustered bootstrap-t to get p-values for MLE models, the function considerably reduces the probability of over rejection of the null hypotheses. The distribution used to reconstruct the residuals is Mammen or Rademacher.

#' @param model A nonlinear model (glm,polr) estimated using \code{lm}. If the mdoel is estimated by OLS, please use \code{wboot} instead.
#' @param cluster A formula to specify the variable of cluster.
#' @param beta A key variable of interest. To denote the key varibale to focus on, one could save a amount of time because the algorithm will only report the outcome of key variable. The defalt value is "All" which will report all the results of every coefficents.
#' @param type The type of distributions for wild boostrap. Default value is Mammen. One could also use Rademacher though the Mammen distribution seems to do much better than Rademacher for small cluster size.
#' @param confid The size of test
#' @param R The number of boostrap replicates. The default number is 400.
#' @param seed A single value passed to \code{set.seed}.
#'
#' @author Jiawei Fu \email{jiawei.fu@duke.edu}
#'
#' @examples
#' \dontrun{
#'
#' ##### glm model #####
#' data(wv6_equ)  # input data
#' # create binomial dependent variable
#' wv6_bi <- as.data.frame(wv6_equ)
#' wv6_bi$income_equ0[wv6_bi$income_equ < 6] <- 0
#' wv6_bi$income_equ0[wv6_bi$income_equ > 5] <- 1
#' glm_equ <- glm(income_equ0 ~ income + age + gender, data = wv6_equ) # glm model
#' swildboott(glm_equ, ~country, R = 500) # find original p values are over estimated
#'
#' ##### polr model #####
#' data(wv6_equ)  # input data
#' require(MASS)
#' wv6_equ <- as.data.frame(wv6_equ)
#' polr_equ <- polr(factor(income_equ) ~ income + age + gender, data = wv6_equ) # polr model
#' swildboott(polr_equ, ~country, R = 500) # find original p value is over estimated
#'}
#'
#'
#' @references Kline, Patrick and Andres Santos. 2012 "A Score Based Approach to Wild Bootstrap Inference" \emph{Journal of Econometric Methods} 2012; 1(1): 23–41
#'
#' @importFrom utils txtProgressBar
#' @importFrom utils setTxtProgressBar
#' @importFrom compiler cmpfun
#' @importFrom sandwich estfun
#' @import stats
#' @import zoo
#' @import MASS
#'
#' @export


swildboott <- function(model, cluster, beta='All',type = 'mammen',confid = 0.05, R = 400, seed = NULL ){


  # check nonlinear model
  # check no c() in the global env

  cl <- match.call()

  if(!is.null(seed)) set.seed(seed)

  ####put in estfun_glm

  estfun_glm <- function(x, y)  #restricted
  {
    xmat <- model.matrix(x)
    ymat <- model.matrix(y) #add ymat including imposed variable
    xmat <- naresid(x$na.action, xmat)
    if(any(alias <- is.na(coef(x)))) xmat <- xmat[, !alias, drop = FALSE]
    wres <- as.vector(residuals(x, "working")) * weights(x, "working")
    dispersion <- if(substr(x$family$family, 1, 17) %in% c("poisson", "binomial", "Negative Binomial")) 1
    else sum(wres^2, na.rm = TRUE)/sum(weights(x, "working"), na.rm = TRUE)
    rval <- wres * ymat / dispersion #use ymat including imposed variable
    attr(rval, "assign") <- NULL
    attr(rval, "contrasts") <- NULL
    res <- residuals(x, type = "pearson")
    if(is.ts(res)) rval <- ts(rval, start = start(res), frequency = frequency(res))
    if(is.zoo(res)) rval <- zoo(rval, index(res), attr(res, "frequency"))
    return(rval)
  }

  ####put in estfun_polr
  estfun_polr <- function(x, yy, pos)   # yy is unrestrict, x is restrict
  {
    ## link processing
    mueta <- x$method
    if(mueta == "logistic") mueta <- "logit"
    mueta <- make.link(mueta)$mu.eta

    ## observations
    ymat <- model.matrix(yy)[, -1L, drop = FALSE]

    n <- nrow(ymat)
    k <- ncol(ymat)   # should diff
    m <- length(x$zeta)

    mf <- model.frame(yy) #maybe y
    y <- as.numeric(model.response(mf)) #maybe y
    w <- model.weights(mf)
    if(is.null(w)) w <- rep(1, n)

    ## estimates
    prob <- x$fitted.values[cbind(1:n, y)]

    coe_un <- yy$coefficients
    coe_re <- x$coefficients

    for(i in 1:length(coe_un)){
      if(i==pos){coe_un[names(coe_un)[i]] <- 0}
      else{coe_un[names(coe_un)[i]] <- coe_re[names(coe_un)[i]]}
    }

    xb <- if(k >= 1L) as.vector(ymat %*% coe_un) else rep(0, n) ## ymat y$coe

    zeta <- x$zeta

    lp <- cbind(0, mueta(matrix(zeta, nrow = n, ncol = m, byrow = TRUE) - xb), 0) #zeta

    ## estimating functions
    rval <- matrix(0, nrow = n, ncol = k + m + 2L) # k
    if(k >= 1L) rval[, 1L:k] <- (-ymat * as.vector(lp[cbind(1:n, y + 1L)] - lp[cbind(1:n, y)]))
    rval[cbind(1:n, k + y)] <- -as.vector(lp[cbind(1:n, y)])
    rval[cbind(1:n, k + y + 1L)] <- as.vector(lp[cbind(1:n, y + 1L)])
    rval <- rval[, -c(k + 1L, k + m + 2L), drop = FALSE]
    rval <- w/prob * rval

    ## dimnames and return
    dimnames(rval) <- list(rownames(ymat), c(colnames(ymat), names(x$zeta)))
    return(rval)
  }


  variables <- all.vars(model$call)  # all variables



  if(!as.character(cl[[3]][[2]]) %in% names(model$coefficients)){
    cluster_tmp <- expand.model.frame(model, cluster, na.expand = FALSE)
    # delete factor in polr model data set
    if(model$call[[1]]=='polr') names(cluster_tmp)[[1]] <- variables[1]
    df <- na.omit(cluster_tmp)
    cluster <- model.frame(cluster, df, na.action = na.pass)
    if(model$call[[1]]=='glm') {
      model$call[[4]] <- quote(df)
    } #glm has family
    else {
      cat('Refitting the polr model costs some seconds.', '\n')
      model$call[[3]] <- quote(df)}
    model_use <- eval(model$call)#refit the model with new df
  } else {
    df <- na.omit(model$model)
    cluster <- model.frame(cluster, df, na.action = na.pass)
    model_use <- model
  }


  # ***check if there is NA in cluster, should I refit original model

  # Factors in our clustering variables can potentially cause problems
  # Blunt fix is to force conversion to characters
  i <- !sapply(cluster, is.numeric)
  cluster[i] <- lapply(cluster[i], as.character)


  ####*****step1:xebar
  cat('Remind: The outcome LM is unscaled', '\n')

  n <- length(unique(cluster[,1])) #total num of clusters

  j <- sapply(sort(unique(cluster[,1])), function(x) sum(cluster==x)) #each cluster size

  "%w/o%" <- function(x, y) x[!x %in% y]

  ### hessian
  j_long <- rep(NA,nrow(df))
  nouse <- Map(function(x,y) {j_long[cluster==x] <<- y}, sort(unique(cluster[,1])), j)


  # original

  score_un_tmp <- apply(estfun(model_use)[,c(1:length(coef(model_use)))], 2L, tapply, cluster[[1]], sum)

  score_un_bar <- (colSums(score_un_tmp/j))/n

  minusbar_un <- t(t(j))%*%score_un_bar  ## matrix
  sigma_tmp_un <- (score_un_tmp-minusbar_un)/j
  sigma_un <- (crossprod(sigma_tmp_un))/(n-1)


  ### prepare to for loop
  indvar <- attr(model_use$terms, "term.labels") # independent variables

  if(!beta=='All') {focus <- beta} else {focus <- indvar}
  if(R<2) stop('R should be larger than 1')

  output_cont <- rep(NA, 5*length(focus)) #put the output

  # start for each var
  for(v in 1:length(focus)){

    if (model$call[[1]]=='glm') {k <- v+1} else {k <- v}
    cat("Bootstrap the variable:", names(coef(model_use))[k],'\n')

    pos <- which(focus[v]==indvar) ### be careful of the pos now is in the IV

    ### get res_formula and est_use
    if (model$call[[1]]=='glm'){
      form_res <- as.formula( paste (variables[1], "~", paste( indvar[1:length(indvar) %w/o% (pos)], collapse= " + " ) ) )
      model_fam <- model$family[[1]]
      model_res <- glm(form_res,family = model_fam, data = df)
      est_use <- estfun_glm(model_res, model_use) ### use glm estfun I modified
      est_ori <- estfun(model_use)
    } else if (model$call[[1]]=='polr'){
      form_res <- as.formula( paste (paste0('factor','(',variables[1],')'), "~", paste( indvar[1:length(indvar) %w/o% (pos)], collapse= " + " ) ) )
      model_res <- polr(form_res, data = df)
      est_use <- estfun_polr(model_res, model_use, pos = pos)[,c(1:length(coef(model_use)))] #use polr estfun I modified
      est_ori <- estfun(model_use)[,c(1:length(coef(model_use)))]
    }

    score_res_tmp <- apply(est_use, 2L, tapply, cluster[[1]], sum)

    hess_tmp <- (t(apply(est_ori,1,function(x) tcrossprod(x))))/j_long
    hess_tmp <- colSums(hess_tmp)/n

    if(model$call[[1]]=='glm') col <- length(indvar)+1
    if(model$call[[1]]=='polr') col <- length(indvar)
    hess <- matrix(hess_tmp, ncol = col, byrow = TRUE)

    r <- t(rep(0,col))  # remember the intercept for glm not polr, row vector
    if(model$call[[1]]=='glm') r[pos+1] <- 1 # remember the intercept
    if(model$call[[1]]=='polr') r[pos] <- 1
    l_tmp <- r %*% solve(hess)


    l_boot <- rep(NA,R)


    for(i in 1:R) {

      if(type == 'mammen') {wild <- sample(c((1 - sqrt(5)) / 2, (sqrt(5) + 1) / 2), size = length(unique(cluster[,1])),replace = TRUE, prob = c((sqrt(5) + 1) / (2 * sqrt(5)), (sqrt(5) - 1) / (2 * sqrt(5))))}
      if(type == 'rademacher') {wild <- c(1, -1)[rbinom(unique(cluster[,1]), size=1, prob=0.5)+ 1]}


      # get modified score

      score_res_bar <- (colSums((score_res_tmp/j)*wild))/n ### t() will be row vector
      minusbar <- t(t(j))%*%score_res_bar  ## matrix
      sigma_tmp <- ((score_res_tmp-minusbar)/j)*wild
      sigma <- (crossprod(sigma_tmp))/(n-1)
      l_boot_bread_2 <-  l_tmp %*% score_res_bar
      l_boot_bread_1 <- t(l_boot_bread_2)
      l_sandwich <- solve(l_tmp %*% sigma %*% solve(hess) %*% t(r))
      l_boot[i] <- l_boot_bread_1 %*% l_sandwich %*% l_boot_bread_2

    }
    # original

    score_un_tmp <- score_res_tmp

    score_un_bar <- (colSums(score_un_tmp/j))/n

    minusbar_un <- t(t(j))%*%score_un_bar  ## matrix
    sigma_tmp_un <- (score_un_tmp-minusbar_un)/j
    sigma_un <- (crossprod(sigma_tmp_un))/(n-1)


    l_bread_2 <-  l_tmp %*% score_un_bar
    l_bread_1 <- t(l_bread_2)
    l_un_sandwich <- solve(l_tmp %*% sigma_un %*% solve(hess) %*% t(r))
    l_un <- l_bread_1 %*% l_un_sandwich %*% l_bread_2

    pvalue <- 1-( sum(abs(l_un[1,1]) > abs(l_boot)) / R)
    quant <- quantile(l_boot,c(confid/2,1-confid/2))

    if(model$call[[1]]=='glm') pos <- pos+1

    output_cont[5*v-4] <- coef(model_use)[pos]
    output_cont[5*v-3] <- l_un
    output_cont[5*v-2] <- quant[[1]]
    output_cont[5*v-1] <- quant[[2]]
    output_cont[5*v] <- pvalue


  }

  output<-matrix(output_cont, nrow=length(focus), ncol=5, byrow = TRUE)

  rownames(output) <- focus
  colnames(output) <- c("Estimate","Ori. LM","2.5% btLM"," 97.5% btLM" , "P-value")

  cat("The score based wild bootstrap outcomes:", "\n")
  printCoefmat(output, P.values= TRUE, has.Pvalue=TRUE)
}
