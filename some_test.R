# test_tib <- tibble(
#   Method = c("bern1", "bern2",
#              "scKDE_uni", "scKDE_2inf", "scKDE_2infplus",
#              "rdens", "amo"))
# 
# test_tib <- test_tib %>%
#   mutate(col1 = NA,
#          col2 = c(1,2,3,4,5,6,7),
#          col3 = c(0.5,0.25,0.5,0.25,0.5,0.25,0.5),
#          col4 = NA)   
# 
# av_test_tib <- test_tib %>%
#   mutate(
#        av_c1c2c3 = rowMeans(select(., starts_with("col")), na.rm = TRUE),
#        av_c1c2 = rowMeans(cbind(col1, col2), na.rm = TRUE),
#        av_c1c4 = rowMeans(cbind(col1, col4), na.rm = TRUE)
#   )
       
  #     rmse_av = rowMeans(select(., starts_with("rmse_")), na.rm = TRUE),
  #     mae_av = rowMeans(select(., starts_with("mae_")), na.rm = TRUE)
  #   ) %>% select(Method,starts_with(c("mse_av","rmse_av","mae_av")))
  

res[2,2] <- NA
res[2,3] <- NA
res[2,4] <- NA

res[2,8] <- NA
res[2,9] <- NA
res[2,10] <- NA

res[5,5] <- NA
res[5,6] <- NA
res[5,7] <- NA

# count nr of invalid folds per method
na_count <- res %>%
  #group_by(Method) %>%
  select(Method, starts_with("mse")) %>%
  mutate(na_count = rowSums(is.na(.))) %>%
  pull(na_count)

# add column with nr of invalid folds per method
res <- res %>% mutate(nr_na_folds = na_count)

# calculate average error over all folds (excl. invalid folds)
av_tib <- res %>%
  mutate(
         mse_av = rowMeans(select(., starts_with("mse_")), na.rm = TRUE),
         rmse_av = rowMeans(select(., starts_with("rmse_")), na.rm = TRUE),
         mae_av = rowMeans(select(., starts_with("mae_")), na.rm = TRUE)
       ) %>% select(Method,mse_av,rmse_av,mae_av,nr_na_folds)



if (any(av_tib$nr_na_folds > 0)) {
  warning("For at least one method, the average error measures were calculated excluding one or more folds. Check the 'nr_na_folds' column in the 'av_tib' tibble to see how many folds were excluded for which method.")
}
