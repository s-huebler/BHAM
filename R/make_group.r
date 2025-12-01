

#' Create grouping for spline design matrix
#'  `r lifecycle::badge("experimental")`
#'
#' @param .names spline matrix names, always in the format "var_name.baseX". Directly from Construct_Smooth_Data
#' #param null_group A indicator if the null space are in its own group, i.e. if null space are penalized
#' @param penalize_null do we penalize the null space, i.e. do we fit the null space into group parameter in bglm group parameter.
#' @param shared_null  A indicator if the null space have a shared indicator with the penalized space
#'
#' @return A vector of lists, where each element list contains variables that belong to the same group
#' @export
#'
#' @importFrom unglue unglue_unnest
#' @import dplyr
#' @importFrom rlang .data
#'
#' @examples
#'
#' raw_dat <- sim_Bai(100, 5)$dat |> data.frame()
#'
#' sm_df <- data.frame(
#'  Var = setdiff(names(raw_dat), "y"),
#'  Func = "s",
#'  Args ="bs='cr', k=5"
#' )
#'
#' dsn_mat <- construct_smooth_data(sm_df, raw_dat)$data
#'
#' make_group(names(dsn_mat))
make_group <- function(.names,
                       penalize_null = TRUE,
                       # null_group = TRUE,
                       shared_null = TRUE
){

  # null_group & shared_null should not be set at TRUE at the same time
  data.frame(names = .names, stringsAsFactors = FALSE) |>
    unglue::unglue_unnest(names, "{var}.{part=pen|null}{ind=\\d*}", remove=FALSE) |>
    mutate(ind = as.numeric(.data$ind)) %>%
    {
      if(shared_null){
        # stop("Not Implemented yet")
        # TODO: need have this change reflected in the model and also passing onto the var_selection problem
        # TODO: idea:
        # TODO: 1. set this as an attribute of the output
        # TODO: 1. when constructing variable selection check with this attribute.
        group_by(., .data$var)
      }
      else{
        # stop("Not Implemented yet. Does not support when lienar and non-linear have different theta")
        group_by(., .data$var, .data$part)
      }
    } %>%
    {
      if(!penalize_null)
        dplyr::filter(., .data$part == "pen")
      else
        .
    } %>%
    dplyr::summarize(res =  list(.data$names), .groups = "drop") |>
    dplyr::pull(.data$res)

  # if(null_group & !shared_null) {
  #   concat_list <- data.frame(names = .names, stringsAsFactors = FALSE) |>
  #     unglue::unglue_unnest(names, "{var}.base{ind}", remove=FALSE) |>
  #     mutate(ind = as.numeric(ind)) |>
  #     group_by(var) |>
  #     arrange(ind) %>%
  #     slice(., n()) |>
  #     summarize(res =  list(names), .groups = "drop") |>
  #     pull(res)
  # }

  # ret <- c(ret, concat_list)
}


#' Grouping from BhGLM
#'
#' @description
#'  Grouping function copied from BhGLM
#'
#'
#' @param all.var
#' @param group
#'
#' @returns
#' @export
#'
#' @examples
Grouping <- function(all.var, group)
{
  n.vars <- length(all.var)
  group.vars <- list()

  if (is.matrix(group))
  {
    if (nrow(group)!=ncol(group) | ncol(group)>n.vars)
      stop("wrong dimension for 'group'")
    if (any(rownames(group)!=colnames(group)))
      stop("rownames should be the same as colnames")
    if (any(!colnames(group)%in%all.var))
      stop("variabe names in 'group' not in the model predictors")
    group.vars <- colnames(group)
    group <- abs(group)
    wcol <- rowSums(group) - diag(group)
    group <- group/wcol
  }
  else{
    if (is.list(group)) group.vars <- group
    else
    {
      if (is.numeric(group) & length(group)>1) {
        group <- sort(group)
        if (group[length(group)] > n.vars) stop("wrong grouping")
      }
      if (is.numeric(group) & length(group)==1)
        group <- as.integer(seq(0, n.vars, length.out = n.vars/group + 1))
      if (is.null(group)) group <- c(0, n.vars)
      group <- unique(group)
      for (j in 1:(length(group) - 1))
        group.vars[[j]] <- all.var[(group[j] + 1):group[j + 1]]
    }
  }
  all.group.vars <- unique(unlist(group.vars))

  if (length(all.group.vars) == n.vars) ungroup.vars <- NULL
  else ungroup.vars <- all.var[which(!all.var %in% all.group.vars)]

  group.new <- c(length(ungroup.vars), length(ungroup.vars) + cumsum(lapply(group.vars, length)))
  var.new <- c(ungroup.vars, unlist(group.vars))

  list(group=group, group.vars=group.vars, ungroup.vars=ungroup.vars,
       group.new=group.new, var.new=var.new)
}
