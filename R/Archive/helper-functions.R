


# create_smooth_name <- function(smth) {
#
#   ret_name <- rep(NA_character_, ncol(smth$X))
#
#   if(length(smth$rank) ==2 ){
#     pen.ind <- colSums(smth$S[[1]])!=0
#     null.ind <- colSums(smth$S[[2]])!=0
#     ret_name[pen.ind] <- paste0(smth$term, ".pen", 1:sum(pen.ind))
#     ret_name[null.ind] <- paste0(smth$term, ".null", 1:sum(null.ind))
#   } else if(length(smth$rank) == 1 ){
#     ret_name <- paste0(smth$term, ".base", 1:ncol(smth$X))
#     warning("Abnormal behaviour when creating spline design matrix. No null space.")
#   } else{
#     stop("Fail to create spline design matrix. Rank length > 2")
#   }
#
#   return(ret_name)
# }
#
#
# # coef <- coefficients(bgam_local)
# # cov <- vcov.bh(bgam_local)
# # test_vars <- function(coef, cov){
# #   pos_vec <- data.frame(names = names(coef)) %>%
# #     mutate(pos = 1:n()) %>%
# #     unglue::unglue_unnest(names, "{var}.base{ind}", remove=FALSE) %>%
# #     mutate(ind = as.numeric(ind)) %>%
# #     filter(!is.na(ind)) %>%
# #     group_by(var) %>%
# #     arrange(ind) %>%
# #     summarize(pos =  list(pos), .groups = "drop")
# #
# #
# #     waldtest.bh(coef, cov, pos_vec$pos[[3]])
# #   })
# #
# # }
