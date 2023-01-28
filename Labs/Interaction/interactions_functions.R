# =========================================================================== #
# functions for interaction lab
# =========================================================================== #

# create QQ plots
plot_qq <- function(data, statistic, df, print_yloc = 5) {
  pvalues <- pchisq(data[, statistic], df = df, lower.tail = F)
  r <- gcontrol2(pvalues, pch=16)
  lambda <- round(median(data[, statistic], na.rm = T) / qchisq(0.5,df), 3)
  lambda <- round(r$lambda,3)
  text(x=1, y=print_yloc, labels=bquote(lambda == .(lambda)), cex=2)
}


# create manhattan plots
plot_manhattan <- function(data, statistic, df) {
  
  # remove NA reults if any
  data <- data[which(!is.na(data[,statistic])), ]
  
  # plot
  neglog.pvalues <- -log10(pchisq(data[,statistic], lower.tail = F, df = df))
  plot(1:nrow(data), neglog.pvalues,
       pch=16, xaxt="n", ylim=c(0, max(neglog.pvalues, 3)),
       ylab="-log(p-value)", xlab="SNPs")
  abline(h=-log10(0.05/nrow(data)), lty=2, lwd=2, col=2)
}


# create table of significant findings
sig_table <- function(data, statistic, df) {
  
  # calculate p-values
  varname = paste0(statistic, "_pval")
  data['lrtdg_pval'] <- pchisq(data[,'lrtdg'], lower.tail = F, df = 1)
  data['lrtgxe_pval'] <- pchisq(data[,'lrtgxe'], lower.tail = F, df = 1)
  data['lrt2df_pval'] <- pchisq(data[,'lrt2df'], lower.tail = F, df = 2)
  data['lrteg_pval'] <- pchisq(data[,'lrteg'], lower.tail = F, df = 1)
  data['lrtctrl_pval'] <- pchisq(data[,'lrtctrl'], lower.tail = F, df = 1)
  data['lrtcase_pval'] <- pchisq(data[,'lrtcase'], lower.tail = F, df = 1)
  data['lrt3df_pval'] <- pchisq(data[,'lrt3df'], lower.tail = F, df = 3)
  
  out <- data[which(data[,varname] < (0.05 / nrow(data))), ]
  
  out[, c("lrtdg_pval", "lrtgxe_pval", "lrt2df_pval", 
           "lrteg_pval", "lrtctrl_pval", "lrtcase_pval", 
           "lrt3df_pval")] <- 
    lapply(out[, c("lrtdg_pval", "lrtgxe_pval", "lrt2df_pval", 
                    "lrteg_pval", "lrtctrl_pval", "lrtcase_pval", 
                    "lrt3df_pval")], 
           function(x) formatC(x, format = 'e', digits = 4))
  
  
  return(out[, c("snp", c("lrtdg_pval", "lrtgxe_pval", "lrt2df_pval", 
                          "lrteg_pval", "lrtctrl_pval", "lrtcase_pval", 
                          "lrt3df_pval"))])
}


# create two-step plots
plot_twostep <- function(data, step1_statistic, step1_df, sizeBin0, alpha) {
  
  # create bins
  m = nrow(output)
  nbins = ceiling(log2(m/sizeBin0 + 1))
  sizeBin = c(sizeBin0 * 2^(0:(nbins-2)), m - sizeBin0 * (2^(nbins-1) - 1) )
  endpointsBin = cumsum(sizeBin)
  rk.pv <- c(1:m)
  grp <- ceiling(log(rk.pv/sizeBin0+1,base=2))
  rep_helper <- c(table(grp))
  alphaBin = alpha * 2 ^ -(1:nbins) / sizeBin
  alphaBin_dat <- rep(alphaBin, rep_helper)
  
  # create pvalues and bin variables
  data$step1p <- pchisq(data[, step1_statistic], df = step1_df, lower.tail = F)
  data$step2p <- pchisq(data[, 'lrtgxe'], df = 1, lower.tail = F)
  data                   <- data[order(-data[,step1_statistic]), ]   # sort by step1_statistic
  data$bin_number        <- as.numeric(grp)
  data$step2p_siglvl     <- as.numeric(alphaBin_dat)
  data$log_step2p_siglvl <- -log10(data$step2p_siglvl)
  data$log_step2p        <- -log10(data$step2p)
  significant_hits <- data[which(data$step2p < data$step2p_siglvl), ]
  # create plot
  # split data into lists (for each bin)
  data_plot <- split(data, f = list(data$bin_number))
  
  # mapinfo just to make sure SNPs are plotted evenly in each bin
  create_mapinfo <- function(x) {
    mapinfo <- seq(unique(x$bin_number) - 1 + 0.1, unique(x$bin_number) - 1 + 0.9, length.out = nrow(x))
    out <- cbind(x, mapinfo)
    return(out)
  }
  data_plot <- lapply(data_plot, create_mapinfo)
  logp_plot_limit = 12
  # PLOT
  binsToPlot = length(data_plot)
  color <- rep(c("#377EB8","#4DAF4A"),100)
  par(mar=c(6, 7, 6, 3))
  bin_to_plot = data_plot[[1]]
  plot(bin_to_plot$mapinfo, bin_to_plot$log_step2p,
       col = ifelse(bin_to_plot$snp %in% significant_hits[, 'snp'], '#E41A1C','#377EB8'),
       pch = ifelse(bin_to_plot$snp %in% significant_hits[, 'snp'], 19, 20),
       cex = ifelse(bin_to_plot$snp %in% significant_hits[, 'snp'], 1.3, 1.7),
       xlab="Bin number for step1 p value",
       ylab="-log10(step2 chiSqGxE p value)",
       xlim=c(0, binsToPlot),
       ylim=c(0, logp_plot_limit),
       axes=F,
       cex.main = 1.7,
       cex.axis = 1.7,
       cex.lab = 1.7,
       cex.sub = 1.7)
  lines(bin_to_plot$mapinfo, bin_to_plot$log_step2p_siglvl, col = "black", lwd=1)
  
  # remaining bins
  for(i in 2:binsToPlot) {
    bin_to_plot = data_plot[[i]]
    points(bin_to_plot$mapinfo, bin_to_plot$log_step2p,
           col = ifelse(bin_to_plot$snp %in% significant_hits$snp, '#E41A1C', color[i]),
           pch = ifelse(bin_to_plot$snp %in% significant_hits$snp, 19, 20),
           cex = ifelse(bin_to_plot$snp %in% significant_hits$snp, 1.3, 1.7),
           cex.main = 1.7,
           cex.axis = 1.7,
           cex.lab = 1.7,
           cex.sub = 1.7)
    lines(bin_to_plot$mapinfo, bin_to_plot$log_step2p_sig, col = "black",lwd = 1)
  }
  axis(1, at = c(-1.5, seq(0.5, binsToPlot-0.2, 1)), label = c(0, seq(1, binsToPlot, 1)), cex.axis = 1.7)
  axis(2, at = c(0:floor(logp_plot_limit)), label = c(0:logp_plot_limit), cex.axis=1.7)
}



# significant findings
# report SNP, step1, step2
sig_twostep <- function(data, step1_statistic, step1_df, sizeBin0, alpha) {
  # create bins
  m = nrow(output)
  nbins = ceiling(log2(m/sizeBin0 + 1))
  sizeBin = c(sizeBin0 * 2^(0:(nbins-2)), m - sizeBin0 * (2^(nbins-1) - 1) )
  endpointsBin = cumsum(sizeBin)
  rk.pv <- c(1:m)
  grp <- ceiling(log(rk.pv/sizeBin0+1,base=2))
  rep_helper <- c(table(grp))
  alphaBin = alpha * 2 ^ -(1:nbins) / sizeBin
  alphaBin_dat <- rep(alphaBin, rep_helper)
  # create pvalues and bin variables
  data$step1p <- pchisq(data[, step1_statistic], df = step1_df, lower.tail = F)
  data$step2p <- pchisq(data[, 'lrtgxe'], df = 1, lower.tail = F)
  data                   <- data[order(-data[,step1_statistic]), ]   # sort by step1_statistic
  data$bin_number        <- as.numeric(grp)
  data$step2p_siglvl     <- as.numeric(alphaBin_dat)
  data$log_step2p_siglvl <- -log10(data$step2p_siglvl)
  data$log_step2p        <- -log10(data$step2p)
  data$step1_pval <- formatC(pchisq(data[, step1_statistic], df = step1_df, lower.tail = F), 
                         format = 'e', digits = 4)
  data$step2_pval <- formatC(pchisq(data[, 'lrtgxe'], df = 1, lower.tail = F),
                         format = 'e', digits = 4)
  significant_hits <- data[which(data$step2p < data$step2p_siglvl), c("snp", "step1_pval", "step2_pval")]
  return(significant_hits)
}


