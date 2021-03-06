produce_pie_chart_per_dataset <-
function(data_table,plot_file){
  
  colnames(data_table) <- c("samples", "organs","membership")
  
  # removing colon samoples
  ind_col <- grep("Colorectum|colon|COAD|Colorectal|LARGE_INTESTINE", data_table$organs, ignore.case=TRUE )

  if(length(ind_col) >=1 ) data_table <- data_table[-ind_col, ]
  data_table$organs <- rep("Dataset", length(data_table$organs) )
  data_table$membership <- gsub("\\.score", "",data_table$membership)
  data_table$membership <- gsub("\\.", "\\-",data_table$membership)
  
  # getting the color information for each subtypes
  my_col <- create_subtype_col()
  
  # matching the color information with the 6 subtypes and assigning colours accordingly
 
  sub.col <- as.character( my_col[ my_col[,1] %in% unique(data_table$membership), 2] )
  
 # if(length(w1)==0){
 #   
 #   n <- length( unique(data_table$membership) )
 #   qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
 #   sub.col <- unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))[1:length( unique(data_table$membership))]
#    data_table$membership <- paste0("Subtype-", data_table$membership)
#  }
  
  
  counts <- table(data_table$membership, data_table$organs)
  counts_prop <- prop.table(counts, margin = 2)*100
  
  df.1 <-round (counts_prop , digits = 2 )
  
  df.2 <- melt(df.1)
  
  colnames(df.2) <- c("subtypes","organs", "value")
  g <- ggplot(df.2, aes(x=subtypes, y=value, fill=subtypes, label=value))+scale_fill_manual(values = sub.col)
  g <- g + geom_bar(stat="identity") + facet_wrap(~ organs)+ geom_text() +labs(title="Percentage of samples in each subtypes")+ coord_polar()+scale_colour_Publication()+ theme_Publication_pie()+scale_y_continuous(breaks=NULL)
  
  # now plotting a mosaic plot
  ggsave(g, file=plot_file)
  
}
