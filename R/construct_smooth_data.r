#' Create a design matrix based on spline configuration
#'
#' @description
#' create_smooth_data creates a design matrix with spline functions.
#' A data frame containing custom spline functions can be supplied or
#' spline arguements can be supplied to
#' in sm_df. The spline design matrix is created using smoothCon from the pacakge mgcv.
#'
#'
#' @param sm_df a data frame that has three columns, _Func_, _Var_, _Args_ (case sensitive). See example
#' @param dat the raw data set
#'
#' @return a list contains
#' \itemize{
#'   \item{data}{the contructed design matrix, the intercept and outcome is not included}
#'   \item{Smooth}{}
#' }
#'
#' @export
#'
#' @importFrom glue glue_data
#' @import mgcv
#' @importFrom  rlang parse_expr .data
#'
#' @examples
#'
#' raw_dat <- sim_Bai(100, 5)$dat %>% data.frame
#'
#' sm_df <- data.frame(
#'  Var = setdiff(names(raw_dat), "y"),
#'  Func = "s",
#'  Args ="bs='cr', k=5"
#' )
#'
#'
#'
#' construct_smooth_data(sm_df, raw_dat)
#'
#'
construct_smooth_data <- function(sm_df, dat){

  # To construct a list of expression
  fml_df <- sm_df |>
    mutate(no_args = is.na(.data$Args)|.data$Args=='',
           arg_str = paste0(', ',.data$Args)) |>
    glue::glue_data("{Func}( {Var}{ifelse(no_args, '', arg_str)})")

  .dat <- data.frame(id = 1:nrow(dat))

  SMs <- list()

  for(i in 1:length(fml_df)){
    raw_sm <- (
      fml_df[i] |>
        rlang::parse_expr() |>
        eval() |>
        smoothCon(object = _, data = dat,
                  scale.penalty = TRUE,
                  absorb.cons = TRUE,
                  null.space.penalty = TRUE,
                  diagonal.penalty = TRUE)
    )[[1]]

    colnames(raw_sm$X) <- create_smooth_name(raw_sm)

    .dat <- cbind(.dat,
                  raw_sm$X)
    SMs[[raw_sm$term]] <- raw_sm
  }


  return(
    list(data = .dat |> select(-id),
         Smooth = SMs)
  )
}


# ************************************************************************************


#' Internal name creation for smooth spline
#'
#' @param smth
#'
#' @returns
#' @keywords internal
#' @examples
create_smooth_name <- function(smth) {

  ret_name <- rep(NA_character_, ncol(smth$X))

  if(length(smth$rank) ==2 ){
    pen.ind <- colSums(smth$S[[1]])!=0
    null.ind <- colSums(smth$S[[2]])!=0
    ret_name[pen.ind] <- paste0(smth$term, ".pen", 1:sum(pen.ind))
    ret_name[null.ind] <- paste0(smth$term, ".null", 1:sum(null.ind))
  } else if(length(smth$rank) == 1 ){
    ret_name <- paste0(smth$term, ".base", 1:ncol(smth$X))
    warning("Abnormal behaviour when creating spline design matrix. No null space.")
  } else{
    stop("Fail to create spline design matrix. Rank length > 2")
  }

  return(ret_name)
}

# ************************************************************************************
