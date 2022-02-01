# normalization of voted output based on reads file- made with prep-guide.pl
#Rscript --vanilla norm_a_voted_circs_df.R /vote_hg19/david_atrt_part2ordered_dcc_approved_by_all_three.csv /reads_per_sample_david_atrt_part2.tsv hg19_normed_out_dcc_appr_voted_circs.csv
# 2 infiles- the voted raw counts table with quants and the file with read counts
args = commandArgs(trailingOnly=TRUE)


# test if there is at least one argument: if not, return an error

if (length(args)<3) {
  stop("you should have given 3 arguments: raw read counts voted circs DF (1.), read counts tsv file correctly ordered (2.) and outfile name (3.)")
}
unused_colnames=c("X","X.1")
# get the df, filter out the info columns and the non-info columns
#args[1]="/home/daric/work_WGGC/circRNA_detection/debug_norm/test_run_BMFZ_2samplesordered_find_circ_approved_by_all_three.csv"
#args[2]="/home/daric/work_WGGC/circRNA_detection/debug_norm/reads_per_sample_test_run_BMFZ_2samples.tsv"
#args[3]="/home/daric/work_WGGC/circRNA_detection/debug_norm/test.tsv"
  
  
heat_dcc=read.csv(file=args[1],check.names = F)

heat_dcc=heat_dcc[,!(colnames(heat_dcc)%in% unused_colnames )]
#heat_dcc=read.csv("/home/daric/work_enclave/hpc_outputs/atrt_david_part2/output_hg19/vote_hg19/david_atrt_part2ordered_dcc_approved_by_all_three.csv",header = T,sep=",")
sapply(heat_dcc,as.numeric)


info_dc <- which(!(sapply(heat_dcc,is.numeric)))# info columns
only_info_columns=sapply(heat_dcc[,info_dc],as.character)


tokeep_dc <- which(sapply(heat_dcc,is.numeric))
only_num_heat_dcc=heat_dcc[ , tokeep_dc,]



clean_only_q=only_num_heat_dcc[,(!(colnames(only_num_heat_dcc) %in% unused_colnames)) ]# filtering out bad colnames- besides the non-numeric ones
# read the reads file
read_counts=read.table(file=args[2], header=T,sep="\t", fill = TRUE,quote = "",check.names = F)

#read_counts=read.table(file="/home/daric/work_enclave/hpc_outputs/atrt_david_part2/output_hg19/reads_per_sample_david_atrt_part2.tsv", header=T,sep="\t", fill = TRUE,quote = "")

if(nrow(read_counts)!=ncol(clean_only_q)){
  print("samplenames cleaned raw counts DF:")
  print(colnames(clean_only_q))
  stop("file read_counts does NOT have exactly the same number of samples as the cleaned DF has, retry with another file ")

}

check_df=cbind(as.character(colnames(clean_only_q)),read_counts)
#head(check_df)
colnames(check_df)=c("samples_raw_df","sample_reads_df","reads")
print("normalizing now with the following numbers: ")

# add fun to check if samplenames are the same everywhere- this lead to failures in the past!


print (check_df)
write.csv(check_df,"normalization_df.csv")


reads_in_mill=as.numeric(as.character(check_df$reads))/1000000
# function to normalize

normalize_only_quants_df <- function(quants_df,n_readsmil_ordered){

  new_df=data.frame(ncol=ncol(quants_df))
  o=0
  for ( col in colnames(quants_df)){
    o=o+1
    col_oi=quants_df[,names(quants_df)==col]
    #col_oi=select(quants_df,col)
    # compare names, check again - 
    name_1=check_df$samples_raw_df[o]
    name_2=check_df$sample_reads_df[o]
    if(!(name_1 == name_2)){
      stop(paste("sample name number",o,"is not the same. please check the normalization_df!","\n names compared:",name_1,"and",name_2))
    }
    
    
    
    num_to_norm_to=as.numeric(as.character(n_readsmil_ordered[o]))
    # 
    print (paste("now normalizing ",col ,"with number ",num_to_norm_to))
    norm_col=col_oi / num_to_norm_to
    new_df=cbind(new_df,norm_col)
  }
  new_df=new_df[,(colnames(new_df)!="ncol")]
  #new_df=select(new_df,-c("ncol"))
  colnames(new_df)=colnames(quants_df)
  return(new_df)
  print ("done.")
}


norm_df_oq=normalize_only_quants_df(clean_only_q,reads_in_mill)

full_norm_df=cbind(only_info_columns,norm_df_oq)
colnames(full_norm_df)=c(as.character(colnames(only_info_columns)),as.character(check_df$sample_reads_df))
# clean DF colnames, again
full_norm_df=full_norm_df[,!(colnames(full_norm_df)%in% unused_colnames )]
write.csv(full_norm_df,file=args[3])
