#(Internal) missing values from AMBROSE 
#source("/castles/nr/projects/2017/gkoutosg-variant-prediction/Laura/AMBROSE/Analysis_Laura/FunctionsPlots.R")


PlotMissingMagda <- function(x){
  
  R <- is.na(x)
  nmis <- colSums(R)
  # sort columnwise
  R <- matrix(R[, order(nmis)], dim(x))
  pat <- apply(R, 1, function(x) paste(as.numeric(x), collapse = ""))
  # sort rowwise
  sortR <- matrix(R[order(pat), ], dim(x))
  if (nrow(x) == 1) {
    mpat <- is.na(x)
  } else {
    mpat <- sortR[!duplicated(sortR), ]
  }
  
  # update row and column margins
  if (all(!is.na(x))) {
    cat(" /\\     /\\\n{  `---'  }\n{  O   O  }\n==>  V <==")
    cat("  No need for mice. This data set is completely observed.\n")
    cat(" \\  \\|/  /\n  `-----'\n\n")
    mpat <- t(as.matrix(mpat, byrow = TRUE))
    rownames(mpat) <- table(pat)
  } else {
    if (is.null(dim(mpat))) {
      mpat <- t(as.matrix(mpat))
    }
    rownames(mpat) <- table(pat)
  }
  
  r <- cbind(abs(mpat - 1), rowSums(mpat))
  r <- rbind(r, c(nmis[order(nmis)], sum(nmis)))
  
  
  op <- par(mar = rep(0, 4))
  on.exit(par(op))
  plot.new()
  
  R <- t(as.matrix(r[1:nrow(r) - 1, 1:ncol(r) - 1]))
  R <- r[1:nrow(r) - 1, 1:ncol(r) - 1]
  
  adj <- c(0, 0.5)
  srt <- 90
  length_of_longest_colname <- max(nchar(colnames(r))) / 2.6
  
  plot.window(
    xlim = c(-1, ncol(R) + 1), #c(-1, 2),
    ylim = c(-1, nrow(R) + length_of_longest_colname), #c(-1, 100), #,
    asp = 0.5
  )
  
  TextSize <- 0.35
  
  M <- cbind(c(row(R)), c(col(R))) - 1
  shade <- ifelse(R[nrow(R):1, ], mdc(1), mdc(2))
  rect(M[, 2], M[, 1], M[, 2] + 1, M[, 1] + 1, col = shade)
  for (i in 1:ncol(R)) {
    text(i - .5, nrow(R) + .3, colnames(r)[i], adj = adj, srt = srt, cex = TextSize )
    text(i - .5, -2.5, nmis[order(nmis)][i], adj = c(0.5,0.0), srt = 90,  cex = TextSize) #.3
  }
  for (i in 1:nrow(R)) {
    text(ncol(R) + .3, i - .5, r[(nrow(r) - 1):1, ncol(r)][i], adj = 0,  cex = TextSize)
    text(-.3, i - .5, rownames(r)[(nrow(r) - 1):1][i], adj = 1,  cex = TextSize)
  }
  text(ncol(R) + .3, -.3, r[nrow(r), ncol(r)],  cex = TextSize)
  
  return(r)
  
}





or_plot_Lau <- function(data, dependent, explanatory, breaks, table_text_size){
  
  #data <- dfProcessedOneCD
  #dependent <- "clavien_dindo"
  #explanatory <- c("age", "sex")
  #breaks <- c(0.5, 1, 5, 10, 20, 30)
  #breaks <- NULL
  confint_type <- NULL
  random_effect <- NULL
  glmfit <- NULL
  factorlist <- NULL
  #table_text_size <- 3.5
  remove_ref <- FALSE
  column_space <- c(-0.5, 0, 0.5)
  suffix <- ": OR (95% CI, p-value)"
  title_text_size <- 13
  plot_opts <- NULL
  table_opts <- NULL
  prefix <- ""
  
  
  # Generate or format factorlist object
  if(!is.null(factorlist)){
    if(is.null(factorlist$Total)) stop("summary_factorlist function must include total_col=TRUE")
    if(is.null(factorlist$fit_id)) stop("summary_factorlist function must include fit_id=TRUE")
  }
  
  if(is.null(factorlist)){
    factorlist = summary_factorlist(data, dependent, explanatory, total_col=TRUE, fit_id=TRUE)
  }
  
  if(remove_ref){
    factorlist = factorlist %>%  
      dplyr::mutate(label = ifelse(label == "", NA, label)) %>% 
      tidyr::fill(label) %>% 
      dplyr::group_by(label) %>%
      dplyr::filter(dplyr::row_number() != 1 | 
                      dplyr::n() > 2 |
                      levels %in% c("Mean (SD)", "Median (IQR)")
      )%>% 
      rm_duplicate_labels()
  }
  
  if(is.null(breaks)){
    breaks = scales::pretty_breaks()
  }
  
  # Confidence intervals, default to "profile" for glm and "Wald" for glmer
  if(is.null(confint_type) && is.null(random_effect)){
    confint_type = "profile"
  } else if(is.null(confint_type) && (!is.null(random_effect) | inherits(glmfit, "glmerMod"))){
    confint_type = "default"
  }
  
  # Generate or format glm
  if(is.null(glmfit) && is.null(random_effect)){
    glmfit = glmmulti(data, dependent, explanatory)
    glmfit_df_c = fit2df(glmfit, condense = TRUE, estimate_suffix = " (multivariable)",
                         confint_type = confint_type)
  } else if(is.null(glmfit) && !is.null(random_effect)){
    glmfit = glmmixed(.data, dependent, explanatory, random_effect)
    glmfit_df_c = fit2df(glmfit, condense = TRUE, estimate_suffix = " (multilevel)",
                         confint_type = confint_type)
  }
  if(!is.null(glmfit) && is.null(random_effect)){
    glmfit_df_c = fit2df(glmfit, condense = TRUE, estimate_suffix = " (multivariable)",
                         confint_type = confint_type, estimate_name = "OR", exp = TRUE)
  } else if(!is.null(glmfit) && !is.null(random_effect)){
    glmfit_df_c = fit2df(glmfit, condense = TRUE, estimate_suffix = " (multilevel)",
                         confint_type = confint_type, estimate_name = "OR", exp = TRUE)
  }
  
  glmfit_df = fit2df(glmfit, condense = FALSE, confint_type = confint_type,  estimate_name = "OR", exp = TRUE)
  
  # Merge
  df.out = finalfit_merge(factorlist, glmfit_df_c)
  df.out = finalfit_merge(df.out, glmfit_df, ref_symbol = "1.0")
  
  # Remove proportions from total column and make continuous explanatory reflect dataset
  df.out$Total = stringr::str_remove(df.out$Total, " \\(.*\\)") %>% 
    as.numeric()
  df.out$Total[which(df.out$levels %in% c("Mean (SD)", "Median (IQR)"))] = dim(data)[1]
  
  # For continuous variables, remove level label
  df.out$levels[which(df.out$levels %in% c("Mean (SD)", "Median (IQR)"))] = "-"
  
  # Remove unwanted lines, where there are more variables in model than wish to display.
  # These not named in factorlist, creating this problem. Interactions don't show on plot.
  if (any(
    is.na(df.out$label)
  )
  ){
    remove_rows = which(is.na(df.out$label)) # This row doesn't work when is.na == FALSE, hence if()
    df.out = df.out[-remove_rows,]
  } else {
    df.out
  }
  
  # Fix order
  df.out$levels = as.character(df.out$levels)
  df.out$fit_id = factor(df.out$fit_id, levels = df.out$fit_id[order(-df.out$index)])
  
  # Plot
  g1 = ggplot(df.out, aes(x = as.numeric(OR), xmin = as.numeric(L95), xmax  = as.numeric(U95),
                          y = fit_id))+
    geom_errorbarh(height=0.2) +
    geom_vline(xintercept = 1, linetype = "longdash", colour = "black")+
    geom_point(aes(size = Total), shape=22, fill="darkblue")+
    scale_x_continuous(trans="log10", breaks= breaks)+
    xlab("Odds ratio (95% CI, log scale)")+
    theme_classic(11)+
    theme(axis.title.x = element_text(),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.line.y = element_blank(),
          axis.ticks.y = element_blank(),
          legend.position="none")
  
  t1 = ggplot(df.out, aes(x = as.numeric(OR), y = fit_id))+
    annotate("text", x = column_space[1], y = df.out$fit_id, label=df.out[,2], hjust=0, size=table_text_size)+
    annotate("text", x = column_space[2], y = df.out$fit_id, label=df.out[,3], hjust=1, size=table_text_size)+
    annotate("text", x = column_space[3], y = df.out$fit_id, label=df.out[,8], hjust=1, size=table_text_size)+
    theme_classic(14)+
    theme(axis.title.x = element_text(colour = "white"),
          axis.text.x = element_text(colour = "white"),
          axis.title.y = element_blank(),
          axis.text.y = element_blank(),
          axis.ticks.y = element_blank(),
          line = element_blank())
  
  # Add optional arguments
  g1 = g1 + plot_opts
  t1 = t1 + table_opts
  
  # Add dependent name label
  title <- 	paste0(dependent, ": OR (95% CI, p-value)")
  #plot_title(data, dependent, dependent_label = dependent_label, prefix = prefix, suffix = suffix)
  
  pdf("OddRatioTable.pdf", 13, 15)
  gridExtra::grid.arrange(t1, g1, ncol=2, widths = c(3,2),
                          top=grid::textGrob(title, x=0.02, y=0.2,
                                             gp=grid::gpar(fontsize=title_text_size), just="left"))
  dev.off()
  
  gridExtra::grid.arrange(t1, g1, ncol=2, widths = c(3,2),
                          top=grid::textGrob(title, x=0.02, y=0.2,
                                             gp=grid::gpar(fontsize=title_text_size), just="left"))
  
}
