
#' Title
#'  `r lifecycle::badge("experimental")`
#'
#' @param Smooth The smooth object from construct_smooth_data
#' @param dat The testing data to construct the new design matrix.
#'
#' @return a data frame containing the trasnformed desgin matrix for the testing data
#' @export
#'
#' @import purrr
#'
#' @examples
#' raw_dat <- sim_Bai(100, 5)$dat |> data.frame()
#' test_dat <- sim_Bai(100, 5)$dat |> data.frame()
#'
#' sm_df <- data.frame(
#'  Var = setdiff(names(raw_dat), "y"),
#'  Func = "s",
#'  Args ="bs='cr', k=5"
#' )
#'
#' dsn_smooth <- construct_smooth_data(sm_df, raw_dat)$Smooth
#'
#' make_predict_dat(dsn_smooth, test_dat)
#'
make_predict_dat <- function(Smooth, dat){

  map_dfc(Smooth,
          .f= function(sm, .dat){
            ret <- mgcv::PredictMat(sm, data = .dat)
            colnames(ret) <- create_smooth_name(sm)
            ret |> data.frame()
          },
          .dat = dat
  )
}

