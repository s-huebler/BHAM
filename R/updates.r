# The following code is copied and adapted from the R package BhGLM
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



# Functions Start Here ---------------------------------------------------

#' Title
#'
#' @param prior
#' @param df
#' @param b0
#' @param ss
#' @param theta
#' @param group.vars
#'
#' @returns
#' @export
#'
#' @examples
update.scale.p.group <- function(prior="mde", df=1, b0, ss, theta, group.vars)
{
  if (prior == "mde"){
    den0 <- (2 * ss[1])^(-1) * exp(-abs(b0)/ss[1]) # de density
    den1 <- (2 * ss[2])^(-1) * exp(-abs(b0)/ss[2])
  }
  if (prior == "mt"){
    den0 <- (ss[1])^(-1) * (1 + b0^2/(df * ss[1]^2))^(-(df + 1)/2) # t density
    den1 <- (ss[2])^(-1) * (1 + b0^2/(df * ss[2]^2))^(-(df + 1)/2)
  }

  p <- rep(NA, length(b0))
  names(p) <- names(b0)

  for(j in 1:length(group.vars)){
    # browser()
    vars <- group.vars[[j]]
    if(length(unique(theta[vars]))!=1)
      stop("Posterior Probability of indicator for a group should be the same")
    tmp_theta <- unique(theta[vars])

    #TODO(boyiguo1): Isolate linear and non-linear using name
    vars.lnr <- grep(pattern = ".null\\d+$", x = vars, value = TRUE)
    vars.non.lnr <- grep(pattern = ".pen\\d+$", x = vars, value = TRUE)

    # TODO(boyiguo1): Treating covariates as lnr terms
    if(length(vars.lnr) + length(vars.non.lnr) == 0)
      vars.lnr <- vars

    if(length(vars.lnr) + length(vars.non.lnr) != length(vars) | length(vars.lnr) < 1)
      stop("Not updating all components when updating p and scale")

    #TODO(boyiguo1): Check how many linear cofficients: 1) one, update p 2) more than 1, stop ("not implemented)
    if(length(vars.lnr) == 1){
      p[vars.lnr] <-  tmp_theta * prod(den1[vars.lnr]) / (tmp_theta * prod(den1[vars.lnr]) + (1 - tmp_theta) * prod(den0[vars.lnr]) + 1e-10)
    } else
      stop("More than one cofficients for linear parts. Not implemented yet.")

    if(length(vars.non.lnr)!=0){ #TODO(boyiguo1): Update p for the remaining part.
      p[vars.non.lnr] <- (tmp_theta^2 * prod(den1[vars.non.lnr])) / (tmp_theta^2 * prod(den1[vars.non.lnr]) + ((1 - tmp_theta^2) * prod(den0[vars.non.lnr])) + 1e-10)
    }
  }

  if(any(is.na(p)))
    stop("Failed to update p. NA values exist among p")

  scale <- 1/((1 - p)/ss[1] + p/ss[2] + 1e-10)

  list(scale=scale, p=p)
}

#' @export
update.prior.sd <- function (prior, beta0, prior.scale, prior.df, sd.x, min.x.sd)
{
  prior.scale <- prior.scale + 1e-04
  J <- length(beta0)
  if (prior == "t" | prior == "mt")
    prior.sd <- sqrt((beta0^2 + prior.df * prior.scale^2)/(1 + prior.df))
  if (prior == "de" | prior == "mde")     # prior.scale = lamda in Exp(1/(2*lamda^2) )
    prior.sd <- sqrt(abs(beta0) * prior.scale)

  prior.sd <- ifelse(prior.sd > 1e+04, 1e+04, prior.sd)
  prior.sd <- ifelse(prior.sd < 1e-04, 1e-04, prior.sd)
  prior.sd <- ifelse(sd.x < min.x.sd, 1e-04, prior.sd)
  if (names(beta0)[1] == "(Intercept)") prior.sd[1] <- 1e+10
  prior.sd
}

#' @export
update.ptheta.group <- function(group.vars, p, w, b) # group-specific probability
{
  f <- function(theta, w, p, bb) { # theta ~ beta(1,b)
    sum(p*log(w*theta) + (1-p)*log(1-w*theta)) + mean((bb-1)*log(1-theta))
  }
  theta <- p
  for (j in 1:length(group.vars)) {
    vars <- group.vars[[j]]
    #    theta[vars] <- mean(p[vars])  # posterior mode with theta~beta(1,1)
    theta[vars] <- stats::optimize(f, interval=c(0, 1),
                                   w=w[vars], p=p[vars], bb=b[j], maximum=T)$maximum
  }
  theta <- ifelse(theta < 0.01, 0.01, theta)
  theta <- ifelse(theta > 0.99, 0.99, theta)
  theta <- w * theta

  theta
}

#' @export
update.ptheta.network <- function(theta, p, w)
{
  phi <- 2
  for (j in 1:length(theta)) {
    mu <- w %*% theta
    m <- mu[j] - w[j,j]*theta[j]
    a <- m*phi
    b <- (1-m)*phi
    theta[j] <- (p[j] + a)/(1 + a + b) # posterior mean
  }
  theta <- ifelse(theta < 0.01, 0.01, theta)
  theta <- ifelse(theta > 0.99, 0.99, theta)

  theta
}

#' @export
update.theta.weights <- function (gvars, theta.weights, inter.hierarchy, inter.parents, p)
{
  if (is.null(inter.parents))
    stop("'inter.parents' should be given")
  if (!is.list(inter.parents))
    stop("'inter.parents' should be a list")
  xnames <- strsplit(gvars, split=":", fixed=T)
  inter <- unlist(lapply(xnames, function(x){length(x)}))
  if (length(inter.parents)!=length(inter[inter==2]))
    stop("interactions are not correctly specified in formula or inter.parents")

  p.main <- p[inter==1]
  if (inter.hierarchy=="strong")
    ww <- lapply(inter.parents,
                 function(x, p.main){ p.main[x[1]] * p.main[x[2]] },
                 p.main)
  if (inter.hierarchy=="weak")
    ww <- lapply(inter.parents,
                 function(x, p.main){ (p.main[x[1]] + p.main[x[2]])/2 },
                 p.main)
  theta.weights[inter==2] <- unlist(ww)
  theta.weights
}

#' @export
auto_scale <- function(x, min.x.sd=1e-04)
{
  scale <- apply(x, 2, sd)
  scale <- ifelse(scale<=min.x.sd, 1, scale)
  two <- which(apply(x, 2, function(u) {length(unique(u))==2}))
  scale[two] <- apply(x[, two, drop=F], 2, function(u){max(u)-min(u)})
  scale
}
