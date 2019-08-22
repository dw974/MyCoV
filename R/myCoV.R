# Hello, world!
#
# This is an example function named 'hello'
# which prints 'Hello, world!'.
#
# You can learn more about package authoring with RStudio at:
#
#   http://r-pkgs.had.co.nz/
#
# Some useful keyboard shortcuts for package authoring:
#
#   Build and Reload Package:  'Ctrl + Shift + B'
#   Check Package:             'Ctrl + Shift + E'
#   Test Package:              'Ctrl + Shift + T'

myCoV = function (fasta=NULL, temp_dir=NULL){

  if (is.null(fasta)) stop("Please provide a fasta file input")
  if (is.null(temp_dir)) temp_dir=tempdir()

  print("Welcome to MyCov")
  print("Reading query fasta file.")
  tryCatch(expr={fas=readDNAStringSet(fasta)},error=function(e){stop("ERROR: Please provide a valid fasta file.")})
  print(paste0("The query file contains ",length(fas)," sequences."))

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
  df$best_hit[match(tab2$qaccver,df$query)]=tab2$best_hit
  df$pairwise_identity[match(tab2$qaccver,df$query)]=tab2$identity
  df$predicted_subgenus=predictions$Assigned_subgenus[match(df$best_hit,predictions$X)]
  df$predicted_genus=predictions$Assigned_genus[match(df$best_hit,predictions$X)]
  df$posterior_probability=sapply(1:dim(df)[1],function(x) eval(parse(text = paste0("predictions$",df$predicted_subgenus[x],"[",match(df$best_hit[x],predictions$X),"]")))*100)

  return(df)

}

plot_CoV=function(df=NULL){
  formattable::formattable(df,list(posterior_probability=formatter("span", style = function(x) style(display = "inline-block",
                                                                                                     `border-radius` = "0px", `padding-right` = "0px",
                                                                                                     `background-color` = formattable::csscolor(RColorBrewer::colorRampPalette(c("red", "lightgreen"))(100)[round(x)]))),
                                   pairwise_identity=formatter("span", style = function(x) style(display = "inline-block",
                                                                                                 `border-radius` = "0px", `padding-right` = "0px",
                                                                                                 `background-color` = formattable::csscolor(RColorBrewer::colorRampPalette(c("red", "lightblue"))(100)[round(x)]))),
                                   query=formattable::formatter("div", style = function(x) style(display = "inline-block", width = "200px", "align"="left","word-wrap"="break-word")),
                                   best_hit=formattable::formatter("div", style = function(x) style(display = "inline-block",width = "200px", "align"="left","word-wrap"="break-word"))
  ))
}

