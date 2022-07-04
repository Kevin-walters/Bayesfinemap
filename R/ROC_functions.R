plot_avg_ROC <- function(causal_snp, data_matrix, num_snps, num_datasets,
                      add_to_plot, max_x, line_col, plot_title = NULL, calc_auc = T,
                      x_labs = NA, y_labs = NA, lty = 1, ...){
  if(requireNamespace("ROCR",quietly = TRUE) == F){
    cat("Error: ROCR package is not loaded.")
    stop()
  }
  causal_ind <- rep(0, num_snps)
  causal_ind[causal_snp] <- 1
  data_for_roc <- vector("list", length = num_datasets)
  causal_snp_list <- vector("list", length = num_datasets)
  auc_ROCR <- vector("numeric", length = num_datasets)
  for(dataset in 1:num_datasets) {
    data_for_roc[[dataset]] <- data_matrix[,dataset]
    causal_snp_list[[dataset]] <- causal_ind
    pred <- ROCR::prediction(data_for_roc[[dataset]], labels = causal_ind)
    auc_by_data <- ROCR::performance(pred, measure = "auc",fpr.stop=0.1)
    #print(auc_by_data)
    auc_ROCR[dataset] <- auc_by_data@y.values[[1]]
  }
  pred <- ROCR::prediction(data_for_roc, labels = causal_snp_list)
  perf <- ROCR::performance(pred, measure="tpr", x.measure="fpr")
  #ROCR::plot(perf,avg="none", xlim = c(0,0.03))
  ROCR::plot(perf,avg="vertical", add = add_to_plot, col = line_col,
             xlim = c(0, max_x), main = plot_title, xlab = x_labs,
             ylab = y_labs, lty = lty, ...)
  #ROCR::plot(perf,avg="vertical",spread.estimate="boxplot",add=TRUE)
  if(calc_auc == T){
    return(auc_ROCR)
  }
}




