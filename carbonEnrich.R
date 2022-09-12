carbonEnrich <- function(ions, c.chain, lipids){
  n = length(ions)
  N = length(lipids$LipidIon)
  out <- c()
  for(i in 1:nrow(c.chain)){
    k = sum(grepl(c.chain[i,1], ions))-
      sum(grepl(paste(c.chain[i,1], "p", sep = ""), ions))-
      sum(grepl(paste(c.chain[i,1], "e", sep = ""), ions))-
      sum(grepl(paste("d", c.chain[i,1], sep = ""), ions))
    K = sum(grepl(c.chain[i,1], lipids$LipidIon))
    p = phyper(k, K, N-K, n, lower.tail = FALSE) 
    out = rbind(out, c(k,K,k/K,n,p))
    
  }
  rownames(out) <- c.chain[,1]
  colnames(out) <- c("k", "K", "k/K", "n", "p")
  return(out)
}