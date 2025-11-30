# Some of the following code is copied and adapted from the R package BhGLM
#  (https://github.com/nyiuab/BhGLM)
# Reference:
# Yi, N., Tang, Z., Zhang, X., & Guo, B. (2019). BhGLM: Bayesian hierarchical
# GLMs and survival models, with applications to genomics and epidemiology. Bioinformatics, 35(8), 1419-1421.
#
# Including the MIT Licences of BhGLM
#
# THE SOFTWARE IS PROVIDED “AS IS”, WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
# AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
# OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN
# THE SOFTWARE.



#' Title
#'
#' @param y
#' @param lp
#'
#' @returns
#' @export
#'
#' @examples
measure_cox <- function(y, lp)
{
  if (NROW(y)!=NROW(lp))
    stop("y and lp should be of the same length", call. = FALSE)
  ff <- bacoxph(y ~ lp, prior=De(1, 0), verbose=FALSE)
  deviance <- -2 * ff$loglik[2]
  cindex <- summary(ff)$concordance[[1]]
  measures <- list(deviance = deviance, Cindex = cindex)
  round(unlist(measures), digits=3)
}

measure.glm <- function(y, y.fitted, family, dispersion = 1)
{
  if (NROW(y)!=NROW(y.fitted))
    stop("y and y.fitted should be of the same length", call. = FALSE)
  if (is.null(dispersion)) dispersion <- 1

  mu <- y.fitted
  if (substr(family, 1, 6)=="NegBin" | substr(family, 1, 17)=="Negative Binomial"
      | family=="nb")
    family <- "NegBin"
  fam <- c("gaussian", "binomial", "poisson", "quasibinomial", "quasipoisson", "NegBin")
  if (! family %in% fam)
    stop("Measures for this family have not been implemented yet")

  if (family=="gaussian")
    logL <- dnorm(y, mean=mu, sd=sqrt(dispersion), log=TRUE)
  if (family=="binomial" | family=="quasibinomial"){
    if (is.factor(y)) y <- as.numeric(y) - 1
    L <- dbinom(y, size=1, prob=mu, log=FALSE)
    L <- ifelse(L==0, 1e-04, L)
    logL <- log(L)
  }
  if (family=="poisson" | family=="quasipoisson")
    logL <- dpois(y, lambda=mu, log=TRUE)
  if (family == "NegBin")
    logL <- dnbinom(y, size=dispersion, mu=mu, log=TRUE)

  logL <- sum(logL, na.rm=TRUE)
  deviance <- -2 * logL

  mse <- mean((y - mu)^2, na.rm = TRUE)
  mae <- mean(abs(y - mu), na.rm = TRUE)
  measures <- list(deviance=deviance, mse=mse, mae=mae)
  if (family=="gaussian") {
    R2 <- (var(y, na.rm = TRUE) - mse)/var(y, na.rm = TRUE)
    measures <- list(deviance=deviance, R2=R2, mse=mse, mae=mae)
  }
  if (family=="binomial" | family=="quasibinomial") {
    # stop("Not implemented")

    # if (!requireNamespace("pROC")) install.packages("pROC")
    # require(pROC)
    if (length(unique(y)) > 1)
      AUC <- suppressMessages( pROC::auc(y, mu) )
    else AUC <- NULL
    AUC <- as.numeric(AUC)
    misclassification <- mean(abs(y - mu) >= 0.5, na.rm = TRUE)
    measures <- list(deviance=deviance, auc=AUC, mse=mse, mae=mae,
                     misclassification=misclassification)
  }

  round(unlist(measures), digits=3)
}





# only used in simulation
link.vars <- function(group.vars) {
  all.group.vars <- unique(unlist(group.vars))
  n.vars <- length(all.group.vars)
  n.groups <- length(group.vars)
  linked.vars <- vector(mode = "list", length = n.vars)
  names(linked.vars) <- all.group.vars
  for (i in 1:n.vars) {
    for (j in 1:n.groups)
      if (all.group.vars[i] %in% group.vars[[j]])
        linked.vars[[i]] <- unique(c(linked.vars[[i]], group.vars[[j]]))
    d <- which(linked.vars[[i]] %in% all.group.vars[i])
    linked.vars[[i]] <- linked.vars[[i]][-d]
  }
  linked.vars
}
