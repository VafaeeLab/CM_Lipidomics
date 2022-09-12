ud <- rep("None", nrow(res))
ud[which(res[, "adj.P.Val"] < pval_cutoff & res[,"logFC"]>logFC_cutoff)] = "up"
ud[which(res[, "adj.P.Val"] < pval_cutoff & res[,"logFC"]<logFC_cutoff*-1)] = "down"

ud = as.data.frame(ud)
ud$ion <- rownames(res)
type <- sub("\\(.*", "", ud[,2])
ud$Type <- type
ud = ud[,-2]

# tb<- table(ud)
# tb<- t(t(tb)/colSums(tb))*100

ggplot(ud, aes(Type))+geom_bar(aes(fill = ud),stat = "count", position = position_stack(reverse = TRUE)) + 
  scale_fill_manual(values=c("blue","grey","red"))+ geom_text(size = 3, position = position_stack(vjust = 0.5))+
  coord_flip()

