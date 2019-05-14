##This script will produce only the dates for all of the RS values
.libPaths("~/R/rlib-3.4.3")
#load R modules
library(ape)
library(phangorn)
library(phytools)
library(geiger)
library(readr)

args = commandArgs(trailingOnly=TRUE)
rs_file = args[1]
out_file_num = args[2]


#tree file name
tree_file_loc<-"all_species_no_neanderthal_species.nwk"

#substitution matrix location
sub_mat_loc<-"substitution_matrix.txt"

#folder containing the RS files
rs_file_loc<-"Albers_Paper"

#Input and process the tree file
tree<-read.tree(file=tree_file_loc)
tree<-root(tree, outgroup="dasNov3")
tree<-rotate(tree,33)
tree<-rotate(tree,34)
tree<-rotate(tree,35)
tree<-rotate(tree,36)
tree<-rotate(tree,41)
tree<-rotate(tree,42)
tree<-rotate(tree,45)
tree<-rotate(tree,56)
tree<-rotate(tree,58)
tree<-rotate(tree,59)

#Input and process the substitution matrix
sub_mat<-read_delim(sub_mat_loc, "\t", escape_double = FALSE, trim_ws = TRUE, col_names=TRUE)
sub_mat<-as.matrix(sub_mat)
sub_rows<-sub_mat[,1]

sub_mat<-sub_mat[,-1]
sub_cols<-colnames(sub_mat)
sub_mat<-mapply(sub_mat, FUN=as.numeric)
sub_mat<-matrix(sub_mat, ncol=length(sub_cols))
colnames(sub_mat)<-sub_cols
rownames(sub_mat)<-sub_rows

#Setup the process

sub_spec=c("panTro5","panPan2","gorGor5","ponAbe2","nomLeu3","rheMac8","macFas5","macNem1","cerAty1","papAnu3","chlSab2","manLeu1","nasLar1","colAng1","rhiRox1","rhiBie1","calJac3","saiBol1","cebCap1","aotNan1","tarSyr2","micMur3","proCoq1","eulMac1","eulFla1","otoGar3","mm10","canFam3","dasNov3","Pongo_pygmaeus","hg38")


check_nodes=c("59","58","57","56","45","41","40","35","34","33","32")

check_nodes_p<-c("hg38",check_nodes)
check_nodes_p<-check_nodes_p[-length(check_nodes_p)]

all_results<-data.frame(rs=character(), allele=character(), freq= numeric(), mrc_node = numeric(), other_spec=character(), nspec=numeric(), allele_pres=numeric(), stringsAsFactors = FALSE)

this_rs<-read.phyDat(rs_file, format="fasta", type="USER", contrast=sub_mat)
    this_rs<-subset(this_rs, subset=sub_spec)
    anc.mpr<-ancestral.pars(tree, this_rs, "MPR")
    
    
    
    this_file<-strsplit(rs_file, split=".", fixed=TRUE)[[1]]
    rs_name<-this_file[1]
    allele_name<-this_file[2]
    #allele_freq<-this_file[3]
    
    print(rs_name)
    
    #now i need to figure out how far back the human base goes
    hg38<-anc.mpr$hg38
    hg38<-which(hg38==1)
    ##Should i be keeping those with a length > 1???!!!
    if(ncol(as.character(this_rs))>1){
      print("skipping this rs because length is >1 :")
      print(rs_name)
    }else{
    
    nd<-1
    final_node<-"hg38"
    while(nd<=length(check_nodes)){
      this_node<-check_nodes[nd]
      compare<-anc.mpr[[this_node]]
      if(compare[hg38]>0){
        nd<-nd+1
        final_node<-this_node
      }else{
        #this node no longer contains the human allele which means the previous node is the MRC with this human allele
        nd<-length(check_nodes)+10
      }
    }
    
    
    #print("the MRC witht he human allele is:")
    #print(final_node)
    
    #NOT PLOTTING CURRENTLY
    png.name<-paste(rs_name, ".",allele_name,".png", sep="")
    png(png.name)
    plotAnc(tree,anc.mpr,1,col=c("navy","chartreuse4","magenta4","orange","gray","white"))
    dev.off()
    other_spec=array(character())
    #Get a list of the other species that contain this allele that are not within the clade 
    #don't want to include the decendents of the final_node
    if(final_node != 32){
      if(final_node != "hg38"){
        decendents<-tips(tree,final_node)
        #other_spec=array(data=NA)
      }else{
        decendents<-"hg38"
      }
        d<-1
        while(d<=length(sub_spec)){
          check_spec = sub_spec[d]
          if(check_spec %in% decendents){
            #this species is already in the decedents so don't do anything
          }else{
            check_val<-anc.mpr[[check_spec]]
            if(check_val[hg38]>0){
              #this species has the same allele as the human!
              other_spec<- c(other_spec, check_spec)
            }
          }
      
      
        d<-d+1
        }
    }else{
      other_spec=""
    }
    tot_other = length(other_spec)
    other_spec<-paste(other_spec, collapse=";")
    
    
    #Need to test to see how far back this allele
    # go through the all the check_nodes
    
    nd<-1
    
    last_presence = "32"
    
    
    while(nd<=length(check_nodes_p)){
      this_node<-check_nodes_p[nd]
      if(this_node == "hg38"){
        to_drop = "hg38"
      }else{
        to_drop<-tips(tree,this_node)
      }
      
      to_keep<-sub_spec[! sub_spec %in% to_drop]
      test_this<-subset(this_rs, subset=to_keep)
      test_this<-unique(test_this)
      if(length(test_this)==1){
        #only one character
        test_this<-unlist(test_this)
        if(test_this[1] == 11){
          #outside of this node they are ALL Ns!!!
          last_presence<-this_node
          nd<-length(check_nodes)+10
        }
        
      }
      
      
      nd<-nd+1
    }
    
    
    all_results<-paste(rs_name, allele_name, final_node, other_spec, tot_other, last_presence, sep="\t")
    out_file<-paste(out_file_num,".all_dates.txt",sep="");
    write(all_results,file=out_file,append=TRUE)
    #print(all_results)
    }
 

