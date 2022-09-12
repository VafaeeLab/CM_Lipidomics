VolcanoPlot <- function(res, pval_cutoff = 0.05, abs_logFC_cutoff = 1, fdr_cutoff = 0.01, textify = FALSE, title = ""){
  
  require(calibrate)
  
  # Make a basic volcano plot
  with(res, plot(logFC, -log10(adj.P.Val), pch=20, col = "grey", xlim=c(min(logFC),max(logFC)), main = title))
  
  # Add colored points: red if adj.P.Val<0.2, orange of log2FC>1, green if both)
  
  # with(subset(res, P.Value < pval_cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="green"))
  # 
  # with(subset(res, abs(logFC)>abs_logFC_cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="orange"))
  # 
  # with(subset(res, adj.P.Val<fdr_cutoff ), points(logFC, -log10(adj.P.Val), pch=20, col="blue"))
  # 
  # with(subset(res, adj.P.Val<fdr_cutoff &  abs(logFC)>abs_logFC_cutoff), 
  #      points(logFC, -log10(adj.P.Val), pch=20, col="red"))
  
  with(subset(res, logFC>abs_logFC_cutoff & adj.P.Val<fdr_cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="red"))
  with(subset(res, -1*logFC >abs_logFC_cutoff & adj.P.Val<fdr_cutoff), points(logFC, -log10(adj.P.Val), pch=20, col="blue"))
  
  
  # Label points with the textxy function from the calibrate plot
 
  if(textify){
    with(subset(res, adj.P.Val<pval_cutoff & abs(logFC)>abs_logFC_cutoff), textxy(logFC, -log10(adj.P.Val), labs=rownames(res), cex=.8, offset = 0.5))
  }
  abline(h = -log10(fdr_cutoff), v = c(-abs_logFC_cutoff,abs_logFC_cutoff), lty=2) # last parameter make it dashed line
}
