#' MyCoV takes a fasta file or multi-fasta file as input
#'
#' @param fasta File name: the file containing the fasta-format sequences to be identified
#' @param temp_dir Folder name: the temporary folder in which to store files (eg. blast databse)
#'
#' @return A data frame summarising the BLAST results from the query

MyCoV = function (fasta=NULL, temp_dir=NULL){

  if (is.null(fasta)) stop("Please provide a fasta file input")
  if (is.null(temp_dir)) temp_dir=tempdir()

  print("Welcome to MyCov")

  print(paste0("Temporary files will be written to this location: ",temp_dir))

  print("Reading query fasta file.")
  tryCatch(expr={fas=Biostrings::readDNAStringSet(fasta)},error=function(e){stop("ERROR: Please provide a valid fasta file.")})
  print(paste0("The query file contains ",length(fas)," sequence(s)."))
  print("Removing spaces from sequence names.")
  names(fas)=gsub(" ","_",names(fas))


  print("Generating BLAST database")
  tryCatch(expr={seqinr::write.fasta(as.list(paste(sequences)),names=as.list(names(sequences)),paste0(temp_dir,"/sequences.fasta"))},
           error=function(e){stop("ERROR: YOU MUST HAVE READ/WRITE PERMISSIONS TO THE PROVIDED TEMPORARY DIRECTORY")})
  system(paste0("makeblastdb -dbtype nucl -input_type fasta -in ",temp_dir,"/sequences.fasta -out ",temp_dir,"/sequences"),ignore.stdout = T,ignore.stderr = T)
  print("Done.")

  print("Writing query data to disk")
  seqinr::write.fasta(as.list(paste(fas)),names=as.list(names(fas)),paste0(temp_dir,"/query.fasta"))
  print("Done.")

  if (length(fas)>100){
    print("Running BLASTN. This might take a while (>100 sequences).")
  } else {
    print("Running BLASTN. This shouldn't take long (<100 sequences).")
  }

  system(paste0("blastn -query ",temp_dir,"/query.fasta -db ",temp_dir,"/sequences -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq' -out ",temp_dir,"/results"))
  print("Done.")

  if(!file.info(paste0(temp_dir,"/results"))$size>0) stop("Your sequences returned no significant BLAST hits, perhaps they are not Coronavirus RdRp sequences from the correct locus?")

  print("Loading BLASTN results.")
  tab=read.table(paste0(temp_dir,"/results"),header=F,stringsAsFactors = F)
  colnames(tab)=c("qaccver", "saccver", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qseq", "sseq")

  tab2=tab %>% dplyr::mutate(pid2=pident*length/387) %>% dplyr::group_by(qaccver) %>% dplyr::summarise(best_hit=saccver[pid2==max(pid2)][1],identity=max(pid2))

  df=data.frame(query=names(fas),best_hit=NA,predicted_subgenus=NA,predicted_genus=NA,posterior_probability=NA,pairwise_identity=NA,stringsAsFactors = F)
  df=df[df$query %in% tab2$qaccver,]
  df$best_hit[match(tab2$qaccver,df$query)]=tab2$best_hit
  df$pairwise_identity[match(tab2$qaccver,df$query)]=tab2$identity
  df$predicted_subgenus=predictions$Assigned_subgenus[match(df$best_hit,predictions$X)]
  df$predicted_genus=predictions$Assigned_genus[match(df$best_hit,predictions$X)]
  df$posterior_probability=sapply(1:dim(df)[1],function(x) eval(parse(text = paste0("predictions$",df$predicted_subgenus[x],"[",match(df$best_hit[x],predictions$X),"]")))*100)
  df$recombination=predictions$recombination[match(df$best_hit,predictions$X)]
  print("Done. Now run tabulate_CoV and/or visualise_CoV on the output.")
  print("Thanks for using MyCoV.")
  print("I hope you found the genus:subgenus combo you were looking for.")
  print("Remember to cite our paper.")

  return(df)

}

#' tabulate_CoV takes the output of MyCoV as its input
#' It then uses 'formattable' to provide a summary of the classifications of the queried sequences.
#'
#' @param df data.frame: the output of the function MyCoV
#'
#' @return A 'formattable' object summarising the data

tabulate_CoV=function(df=NULL){

  write.csv(df,"MyCoV_results.csv")
  formattable::formattable(df,list(posterior_probability=formatter("span", style = function(x) style(display = "inline-block",
                                                                                                     `border-radius` = "0px", `padding-right` = "0px",
                                                                                                     `background-color` = formattable::csscolor(grDevices::colorRampPalette(c("red", "lightgreen"))(100)[round(x)]))),
                                   pairwise_identity=formatter("span", style = function(x) style(display = "inline-block",
                                                                                                 `border-radius` = "0px", `padding-right` = "0px",
                                                                                                 `background-color` = formattable::csscolor(grDevices::colorRampPalette(c("red", "lightblue"))(100)[round(x)]))),
                                   query=formattable::formatter("div", style = function(x) style(display = "inline-block", width = "200px", "align"="left","word-wrap"="break-word")),
                                   best_hit=formattable::formatter("div", style = function(x) style(display = "inline-block",width = "200px", "align"="left","word-wrap"="break-word"))
  ))
}

#' visualise_CoV takes the output of MyCoV as its input
#' It then uses 'ggplot2' to provide a summary of the pairwise distances of the queried sequences to their closes BLAST hits in the context of all pairwise comparisons between Coronaviruses .
#'
#' @param df data.frame: the output of the function MyCoV
#'
#' @return A ggplot object summarising the data

visualise_CoV=function(df=NULL){
  ggplot(CoV_dis)+
    geom_histogram(aes(x=dist,fill=compare),bins=50,colour="black")+
    scale_fill_manual(values = c("darkgreen","red"),name="Pairwise Comparison of Subgenera")+
    geom_vline(data=df,aes(xintercept=1-(pairwise_identity/100)),linetype="dashed",colour="black")+
    facet_wrap(~predicted_genus,scale="free")
}


#' plot_similarity takes a mono-fasta file as an iput
#' It then uses 'ggtree' to provide a summary of how the queried sequence compared to classified sequences in our phylogenetic study.
#' It also provides a visual summary of host and country or origin metadata associated with the classified sequences.
#'
#' Please remember to cite ggtree correctly if you go on to use this representation.
#'
#' @param fasta file: a fasta file containing a single dna sequence in fasta format
#' @param temp_dir Folder name: the temporary folder in which to store files (eg. blast databse)
#'
#' @return Saves the tree plot to "output_plot.pdf" in your working directory.

plot_similarity=function(fasta=NULL, temp_dir=NULL){

  if (is.null(fasta)) stop("Please provide a fasta file input")
  if (is.null(temp_dir)) temp_dir=tempdir()

  print("Welcome to MyCov")

  print(paste0("Temporary files will be written to this location: ",temp_dir))

  print("Reading query fasta file.")
  tryCatch(expr={fas=Biostrings::readDNAStringSet(fasta)},error=function(e){stop("ERROR: Please provide a valid fasta file.")})
  if (length(fas)>1) stop("Please provide only a single sequence for comparison")
  print("Removing spaces from sequence name.")
  names(fas)=gsub(" ","_",names(fas))


  print("Generating BLAST database")
  tryCatch(expr={seqinr::write.fasta(as.list(paste(sequences)),names=as.list(names(sequences)),paste0(temp_dir,"/sequences.fasta"))},
           error=function(e){stop("ERROR: YOU MUST HAVE READ/WRITE PERMISSIONS TO THE PROVIDED TEMPORARY DIRECTORY")})
  system(paste0("makeblastdb -dbtype nucl -input_type fasta -in ",temp_dir,"/sequences.fasta -out ",temp_dir,"/sequences"),ignore.stdout = T,ignore.stderr = T)
  print("Done.")

  print("Writing query data to disk")
  seqinr::write.fasta(as.list(paste(fas)),names=as.list(names(fas)),paste0(temp_dir,"/query.fasta"))
  print("Done.")
  print("Running BLASTN. This shouldn't take long.")

  system(paste0("blastn -query ",temp_dir,"/query.fasta -db ",temp_dir,"/sequences -outfmt '6 qaccver saccver pident length mismatch gapopen qstart qend sstart send evalue bitscore qseq sseq' -out ",temp_dir,"/results"))
  print("Done.")

  if(!file.info(paste0(temp_dir,"/results"))$size>0) stop("Your sequences returned no significant BLAST hits, perhaps they are not Coronavirus RdRp sequences from the correct locus?")

  print("Loading BLASTN results.")
  tab=read.table(paste0(temp_dir,"/results"),header=F,stringsAsFactors = F)
  colnames(tab)=c("qaccver", "saccver", "pident", "length", "mismatch", "gapopen", "qstart", "qend", "sstart", "send", "evalue", "bitscore", "qseq", "sseq")
  tab$pid2=tab$pident*tab$length/387
  tab=subset(tab,tab$pid2>70)
  tab2=data.frame(id=tab$saccver,identity=tab$pid2)

  tr=ape::extract.clade(phy = CoV_tree,node = phytools::findMRCA(tree = CoV_tree,tips = paste(tab2$id),type = "node"))

  d = fortify(tr)
  d = subset(d, isTip)
  tips=with(d, label[order(y, decreasing=T)])

  t=data.frame(id=tips,stringsAsFactors = F)
  t=left_join(t,tab2,by="id")
  t=left_join(t,md,by=c("id"="rep_seq"))

  tr1=ggtree(tr) %<+% t
  tr1$data$x=tr1$data$x/max(tr1$data$x)

  a=tr1 +
    geom_tippoint(aes(colour=identity,shape=identity==max(tab2$identity)))+
    geom_tiplab(aes(label=round(identity,digits = 3)),size=1.5,x=1.01)+
    geom_tiplab(aes(label=Assigned_subgenus),size=1.5,x=1.1)+
    geom_tiplab(aes(label=label),size=1.5,x=1.3)+
    scale_colour_gradientn(colours = c("red","yellow","green"))+
    scale_size_manual(values = c(1,5))+
    scale_shape_manual(values = c(20,18))+
    xlim(0,2.5)+
    ylim(0,length(tr$tip.label)+30)+
    theme_tree()+
    theme(legend.position="none")

  t$host_taxon=factor(t$host_taxon)
  t$continent=factor(t$continent)
  t$id=factor(t$id,levels = rev(tips))
  df3=t %>% group_by(id,host_taxon) %>% summarise(count=n())
  df3t=df3 %>% group_by(host_taxon) %>% summarise(count=n())
  df4=t %>% group_by(id,continent) %>% summarise(count=n())
  df4t=df4 %>% group_by(continent) %>% summarise(count=n())

  b=ggplot(df3)+
    geom_tile(data=df3, aes(x=as.numeric(host_taxon),y=as.numeric(id),fill=ifelse(count>0,"blue",NA))) +
    geom_text(data=df3t,aes(x=as.numeric(host_taxon),label=host_taxon,y=length(tr$tip.label)+1),angle=45,hjust=0,size=3)+
    theme_tree()+
    ylim(0,length(tr$tip.label)+30)+
    xlim(0,max(as.numeric(df3$host_taxon))+7)+
    scale_fill_manual(values="blue")+
    theme(legend.position="none")

  c=ggplot(df4)+
    geom_tile(data=df4, aes(x=as.numeric(continent),y=as.numeric(id),fill=ifelse(count>0,"red",NA))) +
    geom_text(data=df4t,aes(x=as.numeric(continent),label=continent,y=length(tr$tip.label)+1),angle=45,hjust=0,size=3)+
    theme_tree()+
    ylim(0,length(tr$tip.label)+30)+
    xlim(0,max(as.numeric(df4$continent))+3)+
    scale_fill_manual(values="red")+
    theme(legend.position="none")

  pdf(file = "output_plot.pdf",width=20,height=10)
  ggtree::multiplot(a,b,c,ncol=3,widths=c(2,1,1))
  dev.off()

  print("Your plot has been written to 'output_plot.pdf' in your working directory.")
  print("The best hit is highlighted by a diamond.")
  wch=which(t$identity==max(t$identity,na.rm = T))
  print("####THESE ARE THE BEST HIT DETAILS#####")
  print(paste0("The best hit corresponds to :",t$id[wch]))
  print(paste0("host : ",t$host_taxon[wch],"  -  ",t$species[wch]))
  print(paste0("Date : ",t$date[wch]))
  print(paste0("Country of origin : ",t$country[wch]))
  print("Thanks for using MyCoV.")
  print("I hope you found the genus:subgenus combo you were looking for.")
  print("Remember to cite our paper.")
  print("If you use plot_similarity, please also cite ggtree correctly")
}
