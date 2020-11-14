rm(list = ls())

# Parameters

# Initialise the RNG
set.seed(1)

# Set parameters for the analysis
map_key <- 'PRKACA_AGC_HUMAN'
nametag <- '_test'
auc_filter <- 0.60
min_cluster <- 6
min_spec <- 0.30

# Choose the column position that you would like to analyse (from 1 to 15)
h <- 6

# Source Omar's script for the generation of PWMs

source('match-tm.r')

# Iterate across all PWMs and then retrieve the position of interest

PWMs_good <- readRDS('PWMs_type_good.rds')

ST_PKs <- data.frame()

for (i in 1:length(PWMs_good)) {
  
  PWM <- PWMs_good[[i]]
  p1 <- PWM[1:20,h]
  ST_PKs <- rbind(ST_PKs,p1)
  row.names(ST_PKs)[dim(ST_PKs)[1]] <- names(PWMs_good)[i]
}

colnames(ST_PKs) <- c('A', 'R', 'N', 'D', 'C', 'Q', 'E', 'G', 'H', 'I', 'L', 'K', 'M', 'F', 'P', 'S', 'T', 'W', 'Y', 'V')
rownames(ST_PKs) <- toupper(rownames(ST_PKs))

# Correct for the misnaming of PDPK1
rownames(ST_PKs)[grep('PDK1',rownames(ST_PKs))] <- 'PDPK1_AGC_HUMAN'

ST_PKs_orig <- ST_PKs

# Retrieve the yeast position vectors

yeast_PKs <- ST_PKs[rapply(strsplit(rownames(ST_PKs),split="_"), function (x) x[3]) %in% "YEAST",] 

# Use affinity propagation clustering of models to generate specificity clusters --------

library(apcluster)
PK_AP <- apcluster(negDistMat(r=2), ST_PKs)

# Remove clusters with number of members below a chosen threshold

logic_index <- logical(length = length(PK_AP))

for (i in 1:length(PK_AP)) {
  
  if (length(PK_AP[[i]]) < min_cluster ) {
    logic_index[i] <- FALSE
  }
  
  else {
    logic_index[i] <- TRUE
  }
}

PK_AP <- PK_AP[logic_index]

# Remove non-specific clusters (i.e. clusters where there is no strong preference for any residue)

logic_index <- logical(length = length(PK_AP))

for (i in 1:length(PK_AP)) {
  
  if (sum(sort(apply(ST_PKs[PK_AP[[i]],],2,mean),decreasing=TRUE)[c(1,2)]) < min_spec) { 
    logic_index[i] <- FALSE
  }
  
  else {
    logic_index[i] <- TRUE
  }
}

PK_AP <- PK_AP[logic_index]


# Merge similar clusters

# AP clustering has a habit of separating within particular specificities; for example, separating out
# kinases with strong specificity for proline against those with weaker specificities. For now we will 
# re-merge such clusters.

numeric_index <- numeric(length = length(PK_AP))

if (length(PK_AP) > 1) {
  
  for (i in 1:length(PK_AP)) {
    
    numeric_index[i] <- which(apply(ST_PKs[PK_AP[[i]],],2,mean)==max(apply(ST_PKs[PK_AP[[i]],],2,mean)))
    
  }
  
  #Treat S/T as single amino acid; also for R/K; also for I/L  
  
  numeric_index <- replace(numeric_index,numeric_index==17,16)
  numeric_index <- replace(numeric_index,numeric_index==12,2)
  numeric_index <- replace(numeric_index,numeric_index==10,11)
  duplicate_aa <- numeric_index[duplicated(numeric_index)]
  
  if  (length(duplicate_aa) > 0) {
    
    l <- length(duplicate_aa)
    
    for (i in 1:l) {
      
      PK_AP <- c(list(unlist(PK_AP[which(numeric_index == duplicate_aa[i])])),PK_AP[which(numeric_index != duplicate_aa[i])])
      
      numeric_index <- numeric(length = length(PK_AP))
      
      for (j in 1:length(PK_AP)) {
        
        numeric_index[j] <- which(apply(ST_PKs[PK_AP[[j]],],2,mean)==max(apply(ST_PKs[PK_AP[[j]],],2,mean)))
        
      }
    } 
  }
}


# Identify the amino acids in each cluster that define the specificity of that cluster, using an arbitrary threshold ------------------------

spec_aa <- list()

for (i in 1:length(PK_AP)) {
  
  spec_order <- order(apply(ST_PKs[PK_AP[[i]],],2,mean),decreasing=TRUE)
  
  for (j in 1:length(spec_order)) {
    
    if (j > 1) {    
      if (sum(apply(ST_PKs[PK_AP[[i]],spec_order[1:j]],2,mean)) > 0.4) {
        spec_aa <- c(spec_aa,list(colnames(ST_PKs)[spec_order[1:j]]))
        break  
      }  
    }
    
    else  {
      
      if (mean(ST_PKs[PK_AP[[i]],spec_order[1:j]]) > 0.4) {
        
        spec_aa <- c(spec_aa,list(colnames(ST_PKs)[spec_order[1:j]]))
        break 
      }
    }
  }
}


# Next I will retrieve co-specific protein kinases that were missed from the clustering procedure:

for (i in 1:length(spec_aa)) {
  
  fav_aa <- spec_aa[[i]][1]  
  lq_freq <- unname(quantile(ST_PKs[PK_AP[[i]],fav_aa],probs=seq(0,1,0.1))[5])
  
  cospec_kin <- rownames(ST_PKs[ST_PKs[,fav_aa] > lq_freq,])
  cluster_kin <- rownames(ST_PKs[PK_AP[[i]],])
  new_kin <- union(cluster_kin,cospec_kin)
  PK_AP[[i]] <- new_kin
  
  jvec <- NULL
  
  # Here we remove potential 'false positives' cluster members by excluding preferred amino acids not found in 'spec_aa';
  # usually these mis-classifications are due to noisy PWMs generated from a small number of substrates
  
  for (j in 1:length(PK_AP[[i]])) {
    
     
    maxres <- colnames(ST_PKs)[which(ST_PKs[PK_AP[[i]][j],] %in% max(ST_PKs[PK_AP[[i]][j],]))]
    diffres <- setdiff(maxres,spec_aa[[i]])
    
    if (length(diffres) > 0){
      jvec <- c(jvec,j)
    }
    
    
  } 
  
  if (length(jvec) > 0) {
    PK_AP[[i]] <- PK_AP[[i]][-(jvec)]
  }
  
}
  
  # Sort out sequences on basis of specificity and name; general processing ------------------------------
  
  # Retrieve kinase names from pre-existing alignment
  
  library(seqinr)
  Human_and_GRID <- read.fasta('report_kinases_al.fasta',seqtype='AA')
  
  # Retrieve the sequence annotations
  Annot <- getAnnot(Human_and_GRID)
  Index <- regexpr("GN",unlist(Annot))
  
  Annot_string <- character()
  Species_string <- character()
  
  for (i in 1:length(Index)) {
    
    # Retrieve gene name
    sub_Annot <- paste(unlist(strsplit(unlist(Annot)[i],""))[Index[i]:length(unlist(strsplit(unlist(Annot)[i],"")))],collapse="")
    
    # Retrieve species
    sub_Species <- unlist(strsplit(unlist(strsplit(unlist(Annot[i]),"_"))[2]," "))[1]
    
    Annot_string <- c(Annot_string,sub_Annot)
    Species_string <- c(Species_string,sub_Species)
  }
  
  # The purpose of this block of code is to match the name of the kinases in the FASTA files
  # with the names used to label the kinase PWMs
  
  yeast_index <- grep("YEAST",Species_string)[1]
  GN_string <- matrix(unlist(strsplit(Annot_string," ")),byrow=T,ncol=3)[,1]
  Gene_name <- matrix(unlist(strsplit(GN_string,"=")),byrow=T,ncol=2)[,2]
  Gene_name1 <- Gene_name[1:(yeast_index-1)]
  New_index <- match(Gene_name1,unlist(lapply(strsplit(rownames(ST_PKs_orig[1:(yeast_index-1),]),"_"), `[[`, 1)))
  Gene_name2 <- Gene_name[yeast_index:length(Gene_name)]
  Kinase_name1 <- paste(Gene_name1,"_",rapply(strsplit(rownames(ST_PKs[1:length(Gene_name1),]),split="_"), function(x) x[2]),"_",Species_string[1:(yeast_index-1)],sep="")
  Kinase_name2 <- paste(Gene_name2,"_",rapply(strsplit(rownames(yeast_PKs),split="_"), function(x) x[2]),"_",Species_string[yeast_index:length(Species_string)],sep="")
  Kinase_name <- c(Kinase_name1,Kinase_name2)
  famindex2 <- which(rapply(strsplit(Kinase_name,split="_"), function (x) x[2]) != "")
  Kinase_name <- Kinase_name[famindex2]
  kinase_gn <- rapply(strsplit(Kinase_name,split="_"),function(x) x[c(1)])
  kinase_sp <- rapply(strsplit(Kinase_name,split="_"),function(x) x[c(3)])
  kinase <- paste(kinase_gn,"_",kinase_sp,sep="")
  id <- paste(Gene_name,"_",Species_string,sep="")
  index <- id %in% kinase
  Human_and_GRID <- Human_and_GRID[index]
  
  # Write out segregated fasta files
  
  Kinase_name <- toupper(Kinase_name)
  
  int_DF <-list()
  spec_name_vec <- character()
  
  command_p1 <- "python /home/david/Documents/Work/PhD/First_year/Python_Scripts/Alignment_map_current/Groupsim_prep.py"
  
  # Iterate through each cluster and partition the alignment into those sequences that match the cluster and
  # those that do not
  
  for (i in 1:length(spec_aa)) {
    
    # Sequences matching the cluster
    spec_seq <- Human_and_GRID[match(rownames(ST_PKs[PK_AP[[i]],]),Kinase_name)]
    
    # Sequences that do not match the cluster
    nonspec_seq <- Human_and_GRID[match(setdiff(rownames(ST_PKs),rownames(ST_PKs[PK_AP[[i]],])),Kinase_name)]
    order_seq <- c(spec_seq,nonspec_seq)
    
    # Length of cluster and non-cluster
    spec_num <- length(rownames(ST_PKs[PK_AP[[i]],]))
    nonspec_num <- length(rownames(ST_PKs)) - spec_num  
    
    # Prepare names of sequence files
    spec_name <- paste(spec_aa[[i]],collapse="")
    spec_name <- paste(spec_name,h,sep="")
    spec_file_name <- paste("spec","_",spec_name,'.fasta',sep="")
    non_spec_file_name <- paste("nonspec","_",spec_name,'.fasta',sep="")
    order_file_name <- paste("order","_",spec_name,'.fasta',sep="")
    
    write.fasta(spec_seq,names = Kinase_name[match(rownames(ST_PKs[PK_AP[[i]],]),Kinase_name)], file.out=spec_file_name,open = "w")
    write.fasta(nonspec_seq,names = Kinase_name[match(setdiff(rownames(ST_PKs),rownames(ST_PKs[PK_AP[[i]],])),Kinase_name)], file.out=non_spec_file_name,open = "w")
    #write.fasta(order_seq,names = Kinase_name[match(c(rownames(ST_PKs[PK_AP[[i]],]),setdiff(rownames(ST_PKs),rownames(ST_PKs[PK_AP[[i]],]))),Kinase_name)], file.out=order_file_name,open = "w")
    
    # The purpose of these commands is to assign tags to 'cluster' and 'non-cluster' kinases
    # so that Groupsim can differentiate them 
    
    arg1 <- spec_file_name
    arg2 <- paste("spec","_",spec_name,'_name.fasta',sep="")
    command <- paste(command_p1,arg1,'g1',arg2)
    system(command)
    
    arg3 <- non_spec_file_name
    arg4 <- paste("nonspec","_",spec_name,'_name.fasta',sep="")
    command <- paste(command_p1,arg3,'g2',arg4)
    system(command)
    
    arg5 <- paste(spec_name,'_name.fasta',sep="")
    com_del <- paste('>',arg5,sep="")
    system(com_del)
    command <- paste('cat',arg2,arg4,'>>',arg5)
    system(command)
    
    ## Remove all columns with gaps in more than 20% of the sequences using trimAl
    
    # trimal
    #setwd("~/Documents/Software/trimAl/source")
    #Sys.setenv(PATH=paste(Sys.getenv("PATH"),"/home/david/Documents/Software/trimAl/source",sep=":"))
    trim_com1 <- 'trimal -in'
    trim_arg1 <- arg5
    trim_arg2 <- paste('-out ',spec_name,nametag,'_name_trim.fasta',sep="")
    trim_arg3 <- paste('-gt 0.8 -colnumbering > ',spec_name,'_colnum.txt',sep="")
    trim_command <- paste(trim_com1,trim_arg1,trim_arg2,trim_arg3)
    system(trim_command)
    
    # Remove blank lines:
    
    colnum <- paste(spec_name,'_colnum.txt',sep="")
    colnum_new <- paste(spec_name,'_colnum_short.txt',sep="")
    noblank <- paste("grep '[^[:blank:]]' <",colnum,">",colnum_new)
    system(noblank)
    
    # This script maps between the full alignment positions and the trimmed alignment positions
    
    DF_arg1 <- 'python /home/david/Documents/Work/PhD/First_year/Python_Scripts/Alignment_map_current/corresponding_numbering.py'
    DF_arg3 <- paste(spec_name,'_colnum_DF.txt',sep="")
    DF_command <- paste(DF_arg1,colnum_new,DF_arg3)
    system(DF_command)
    
    # If the cluster members are predominantly CMGC kinases, then we use a structural mapping
    # from a CDK2 structure. Otherwise, use a structural mapping from PKA
    
    spec_group_table <- table(rapply(strsplit(PK_AP[[i]],split="_"), function (x) x[2]))
    
    if (is.na(spec_group_table['CMGC']) == TRUE) {
      struc <- 'PRKACA_AGC_HUMAN'
      struc_tag <- 'PRKACA'
      structure <- '1atp'
    } else if (spec_group_table['CMGC'] > length(PK_AP[[i]])/2) {
      struc <- 'CDK2_CMGC_HUMAN'
      struc_tag <- 'CDK2'
      structure <- '2cci'
    } else {
      struc <- 'PRKACA_AGC_HUMAN'
      struc_tag <- 'PRKACA'
      structure <- '1atp'
    }
    
    # These are tables of neighbouring residues for each structure
    # This is needed as an input for Multi-Relief 3D
    
    ## 6 A-squared threshold : http://ligin.weizmann.ac.il/cma/
    
    # neighbours_in: raw_input
    # neighbours_out: processed output
    
    neighbours_in <- paste(struc_tag,'_neighbours.txt',sep="")
    neighbours_out <- paste(struc_tag,'_neighbours_new.txt',sep="")
 
    ### PDB mapping to UniProt mapping (parse the XML from PDBe)
    
    library(XML)
    library(RCurl)
    
    xml_file <- paste(structure,'_map.xml',sep="")
    SIFTS_XML <- xmlParse(xml_file)
    SIFTS_XML_root <- xmlRoot(SIFTS_XML)
    nchains <- xmlSize(SIFTS_XML_root) - 2 
    
    siftsmap <- character()
    
    for (x in 1:nchains){
      nseg <- xmlSize(SIFTS_XML_root[[2+x]][[1]])
      for (w in 1:nseg){
        nres <- xmlSize(SIFTS_XML_root[[2+x]][[w]][[1]])
        if (is.null(SIFTS_XML_root[[2+x]][[w]][[1]][[1]][[1]]) == FALSE) {
          chain <- as.character(xmlAttrs(SIFTS_XML_root[[2+x]][[w]][[1]][[1]][[1]])[6])
          for (y in 1:nres){
            PDB_resno <- xmlAttrs(SIFTS_XML_root[[2+x]][[w]][[1]][[y]][[1]])[4]
            PDB_resname <- xmlAttrs(SIFTS_XML_root[[2+x]][[w]][[1]][[y]][[1]])[5]
            if (length(SIFTS_XML_root[[2+x]][[w]][[1]][[y]][[2]]) > 0){
              UniP_resno <- xmlAttrs(SIFTS_XML_root[[2+x]][[w]][[1]][[y]][[2]])[4]
            } else {UniP_resno <- "-"}
            res_pair <- c(chain,PDB_resno,UniP_resno,PDB_resname)
            
            siftsmap <- rbind(siftsmap,res_pair)
          }
        }
      }
    }
    
    if (structure == "2cci") {
      siftsmap <- siftsmap[siftsmap[,1] == "A",]
      replace <- "A"
    } else {
      siftsmap <- siftsmap[siftsmap[,1] == "E",]
      replace <- "E"
    }
    
    neighbours <- read.table(neighbours_in,sep="",skip=2, stringsAsFactors = FALSE)
    neighbours <- neighbours[,c(2,4)]
    neighbours <- apply(neighbours,2, function(x) gsub(replace,"",x))
    neighbours <- apply(neighbours,2, function(x) as.numeric(x))
    neighbours <- unique(neighbours)
  
    
    neighbours <- neighbours
    
    for (a in 1:nrow(neighbours)) {
      for (b in 1:ncol(neighbours)) {
        
        index <- which(as.numeric(siftsmap[,2]) == neighbours[a,b])
        
        neighbours[a,b] <- as.numeric(siftsmap[index,3])
        
      }
    }  
    
    write.table(neighbours,neighbours_out,quote=FALSE,row.names=FALSE,col.names=FALSE)
    
    # This script maps between the untrimmed and trimmed alignment for neighbouring residues of CDK2 or PRKACA
    
    WI_com1 <- 'python /home/david/Documents/Work/PhD/First_year/Python_Scripts/WHATIF_parsing/WHATIF_interchain_parser_neighbours_almap_notrim.py'
    WI_com2 <- neighbours_out
    WI_com3 <- DF_arg3
    WI_com4 <- arg5
    WI_com5 <- struc
    WI_com6 <- paste(struc_tag,'_neighbours_listed_almap_',spec_name,'_paralogue.txt',sep="")
    WI_command <- paste(WI_com1,WI_com2,WI_com3,WI_com4,WI_com5,WI_com6)
    system(WI_command)
    
    
    trunc_com1 <- "python /home/david/Documents/Work/PhD/First_year/Python_Scripts/Alignment_map_current/name_truncate.py"
    trunc_com2 <- paste(spec_name,nametag,'_name_trim.fasta',sep="")
    trunc_com3 <- paste(spec_name,nametag,'_name_trim_trunc.fasta',sep="")
    trunc_com <- paste(trunc_com1,trunc_com2,trunc_com3)
    system(trunc_com)
    
    
    # Multi-RELIEF-3D (custom code based on algorithm description) ---------------------------------------------------------
    
    
    library(seqinr)
    seq_com <- paste(spec_name,nametag,'_name_trim.fasta',sep="")
    Alignment <- read.alignment(seq_com,format="fasta")
    Seq_dis_mat<- as.matrix(dist.alignment(Alignment,matrix="identity",gap=1))
    Seq_dis_mat1 <- Seq_dis_mat[grepl("g1",rownames(Seq_dis_mat)),]
    Seq_dis_mat2 <- Seq_dis_mat[grepl("g2",rownames(Seq_dis_mat)),]
    
    set.seed(1)
    Al_length <- nchar(Alignment[[3]][[1]])
    WeightsM <- numeric(Al_length) 
    iter <- 200
    
    if (length(PK_AP[[i]]) >= 10) {
      sample_size <- 10 
    } else {
      sample_size = length(PK_AP[[i]])
    }
    
    for (i in 1:iter){
      
      Weights <- numeric(Al_length)   
      
      rand_index1 <- sample(nrow(Seq_dis_mat1),sample_size,replace=FALSE)
      rand_index2 <- sample(nrow(Seq_dis_mat2),sample_size,replace=FALSE)
      
      Seq_dis_mat1_rand <- Seq_dis_mat1[rand_index1,]
      Seq_dis_mat2_rand <- Seq_dis_mat2[rand_index2,]
      
      Seq_dis_mat_rand <- rbind(Seq_dis_mat1_rand,Seq_dis_mat2_rand)
      Seq_dis_mat_rand <- Seq_dis_mat_rand[,c(rand_index1,nrow(Seq_dis_mat1)+rand_index2)]
      
      Seq_dis_mat1_rand <- Seq_dis_mat_rand[grepl("g1",rownames(Seq_dis_mat_rand)),]
      Seq_dis_mat2_rand <- Seq_dis_mat_rand[grepl("g2",rownames(Seq_dis_mat_rand)),]
      
      for (j in 1:ncol(Seq_dis_mat_rand)) {
        
        if (grepl("g1",colnames(Seq_dis_mat_rand)[j]) == TRUE) {
          
          nnin <- names(Seq_dis_mat1_rand[,j])[which(Seq_dis_mat1_rand[,j] == sort(Seq_dis_mat1_rand[,j])[2])]
          nnin <- sample(nnin,1)
          nnout <- names(Seq_dis_mat2_rand[,j])[which(Seq_dis_mat2_rand[,j] == sort(Seq_dis_mat2_rand[,j])[1])]
          nnout <- sample(nnout,1)
        }  else  {  
          nnin <- names(Seq_dis_mat2_rand[,j])[which(Seq_dis_mat2_rand[,j] == sort(Seq_dis_mat2_rand[,j])[2])]
          nnin <- sample(nnin,1)
          nnout <- names(Seq_dis_mat1_rand[,j])[which(Seq_dis_mat1_rand[,j] == sort(Seq_dis_mat1_rand[,j])[1])]
          nnout <- sample(nnout,1)
        }
        
        indexin <- which(rownames(Seq_dis_mat_rand) == nnin)
        indexout <- which(rownames(Seq_dis_mat_rand) == nnout)
        Alignment_rand <- Alignment[[3]][c(rand_index1,nrow(Seq_dis_mat1)+rand_index2)] 
        
        for(k in 1:Al_length){
          
          if (unlist(strsplit(Alignment_rand[[j]],""))[k] == unlist(strsplit(Alignment_rand[[indexin]],""))[k]) {
            Weights[k] <- Weights[k] + 1/(sample_size*2)
          } else {
            
            Weights[k] <- Weights[k]
          }
          
          if (unlist(strsplit(Alignment_rand[[j]],""))[k] == unlist(strsplit(Alignment_rand[[indexout]],""))[k]) {
            Weights[k] <- Weights[k] - 1/(sample_size*2)
          } else {
            
            Weights[k] <- Weights[k]
          }  
          
        }
      }
      
      WeightsM <- rbind(WeightsM,Weights)
      
    }
    
    WeightsM <- WeightsM[-1,]
    WeightsMM <- numeric(ncol(WeightsM))
    
    # Take sign-sensitive means 
    
    for (i in 1:ncol(WeightsM)) {
      
      weight.logical <- which(WeightsM[,i] > 0)
      
      if (length(weight.logical) == 0) {
        
        weight2.logical <- which(WeightsM[,i] != 0)
        
        if (length(weight2.logical) == 0) {
          
          WeightsMM[i] <- 0
          
        }  else { WeightsMM[i] <- mean(WeightsM[which(WeightsM[,i] < 0 ),i])}
        
      } else { WeightsMM[i] <- mean(WeightsM[which(WeightsM[,i] > 0 ),i])}
      
    }
    
    order(WeightsMM)
    sort(WeightsMM)
    
    # 3D-neighbour weighting (relies upon the previous mapping between the trimmed and untrimmed alignment)
    
    Al_contacts <- read.table(WI_com6,na.strings="9999")
    Al_contacts <- na.omit(Al_contacts)
    write.table(Al_contacts,paste(struc_tag,'_neighbours_new_al.txt',sep=""),quote=FALSE,row.names=FALSE,col.names=FALSE)
    
    head(Al_contacts)
    partner_list <- list()
    
    for (i in 1:length(unique(Al_contacts[,1]))) {
      index <- which(Al_contacts[,1] == unique(Al_contacts[,1])[i])
      partners <- list(Al_contacts[index,2])
      partner_list <- c(partner_list,partners)
    }
    
    names(partner_list) <- paste("Site_",unique(Al_contacts[,1]),sep="")
    
    WeightsMM3D <- WeightsMM
    
    for (i in 1:length(unique(Al_contacts[,1]))) {
      partner_indexes <- unlist(partner_list[i])
      seq_index <- unique(Al_contacts[,1])[i]
      partner_weight <- WeightsMM[partner_indexes]
      partner_weight_mean <- mean(partner_weight)
      WeightsMM3D[seq_index] <- WeightsMM[seq_index] + partner_weight_mean
    }
    
    MRoutput <- paste(spec_name,nametag,'_name_trim_mrelief.out',sep="")
    write.table(WeightsMM3D,MRoutput,quote=FALSE,row.names=FALSE,col.names=FALSE,sep=" ")
    
    # Multi-RELIEF-mapping (map from the trimmed alignment to the untrimmed alignment and to a reference protein e.g. PRKACA)
    
    MRmap_com1 <- 'python /home/david/Documents/Work/PhD/First_year/Python_Scripts/Alignment_map_current/mrelief_mapping_new_ensembl_notrim.py'
    MRmap_com2 <- MRoutput
    MRmap_com3 <- DF_arg3
    MRmap_com4 <- arg5
    MRmap_com5 <- map_key
    MRmap_com6 <- paste(spec_name,nametag,'_mrelief_mapped.txt',sep="")
    MRmap_command <- paste(MRmap_com1,MRmap_com2,MRmap_com3,MRmap_com4,MRmap_com5,MRmap_com6)
    system(MRmap_command)
    
    
    # SPEER --------------------------------------------------------------------
    
    #SPEER
    
    #setwd("~/Documents/Software/Speer")
    #Sys.setenv(PATH=paste(Sys.getenv("PATH"),"/home/david/Documents/Software/Speer/src",sep=":"))
    
    ss_com1 <- 'SPEER -i'
    ss_com2 <- paste(spec_name,nametag,'_name_trim_trunc.fasta',sep="")
    ss_com3 <- paste("-o ",spec_name,"_SPEER_pre.out",sep="")
    ss_com4 <- paste('-wERate 0','-g 1','-f 1')
    ss_com5 <- paste(spec_num,nonspec_num,sep=" ")
    ss_com <- paste(ss_com1,ss_com2,ss_com3,ss_com4,ss_com5)
    
    command <- paste('head -n -3 ',spec_name,"_SPEER_pre.out > ",spec_name,"_SPEER.out",sep="")
    
    system(ss_com)
    
    system(command)
    
    
    #SPEER mapping (map from the trimmed alignment to the untrimmed alignment and to a reference protein e.g. PRKACA)
    
    ssmap_com1 <- 'python /home/david/Documents/Work/PhD/First_year/Python_Scripts/Alignment_map_current/SPEER_mapping_final_notrim.py'
    ssmap_com2 <- paste(spec_name,"_SPEER.out",sep="")
    ssmap_com3 <- DF_arg3
    ssmap_com4 <- arg5
    ssmap_com5 <- map_key
    ssmap_com6 <- paste(spec_name,nametag,'_SPEER_mapped.txt',sep="")
    ssmap_command <- paste(ssmap_com1,ssmap_com2,ssmap_com3,ssmap_com4,ssmap_com5,ssmap_com6)
    system(ssmap_command)
    
    
    # Gsim --------------------------------------------------------------------
    
    gsim_com1 <- 'python /home/david/Documents/Software/group_sim_sdp/group_sim_sdp.py -m blosum62.bla -w 3 -l 0.50 -o'
    #gsim_com1 <- 'python /home/david/Documents/Software/group_sim_sdp/group_sim_sdp.py -w 3 -l 0.50 -o'
    arg6 <- paste(spec_name,'_gsim.out',sep="")
    arg7 <- paste('\\',spec_name,nametag,'_name_trim_trunc.fasta',sep="")
    command <- paste(gsim_com1,arg6,arg7,'g1','g2')
    system(command)
    
    
    # GroupSim mapping
    
    GSmap_com1 <- 'python /home/david/Documents/Work/PhD/First_year/Python_Scripts/Alignment_map_current/groupsim_mapping_notrim.py'
    GSmap_com2 <- arg6
    GSmap_com3 <- DF_arg3
    GSmap_com4 <- arg5
    GSmap_com5 <- map_key
    GSmap_com6 <- paste(spec_name,nametag,'_gsim_mapped.txt',sep="")
    GSmap_command <- paste(GSmap_com1,GSmap_com2,GSmap_com3,GSmap_com4,GSmap_com5,GSmap_com6)
    system(GSmap_command)
    
    
    # Gsim SDRs
    
    gsim_table <- read.table(GSmap_com6,header=TRUE)
    gsim_table_order <- gsim_table[order(gsim_table[,4],decreasing=TRUE),]
    gsim_table_order <- gsim_table_order[gsim_table_order[,4]!="None",]
    gsim_table_order <- gsim_table_order[1:15,c(3,4)]
    colnames(gsim_table_order)[2] <- 'GS_score'
    
    
    # SS SDRs
    
    ss_table <- read.table(ssmap_com6,header=TRUE,stringsAsFactors = FALSE)
    ss_table_order <- ss_table[order(ss_table[,7],decreasing=FALSE),c(3,7)]
    ss_table_order <- ss_table_order[1:15,]
    
    
    # MRelief SDRs:
    
    mr_table <- read.table(MRmap_com6,header=TRUE)
    mr_table_order <- mr_table[order(mr_table[,4],decreasing=TRUE),]
    mr_table_order <- mr_table_order[1:15,c(3,4)]
    
    
    int <- intersect(intersect(gsim_table_order[,1],round(as.numeric(ss_table_order[,1]))),mr_table_order[,1])
    
    int_DF <- c(int_DF,list(int))
    spec_name_vec <- c(spec_name_vec,spec_name)
    
    
  }
  
  # Output --------------------------------------------------------------------
  
  End_time <- Sys.time()
  print(End_time-Start_time)
  
  names(int_DF) <- spec_name_vec
  
  
  for (i in 1:length(PK_AP)) {
    
    number <- length(PK_AP[[i]])
    aa <- spec_aa[i]
    plot_title <- paste("PWM column ",h,", number = ",number,sep="")
    
    barplot(apply(ST_PKs[PK_AP[[i]],],2,mean),xlab="Amino acids",ylab='Within-cluster PWM score average',main=plot_title,col=c('green','red','purple','blue','yellow','purple','blue','green','red','green','green','red','yellow','pink','orange','brown','brown','pink','pink','green'),ylim=c(0.05,1.0), yaxt="n")
    axis(2, at = seq(0, 1, 0.1), las = 1)
  }
  
  print(int_DF)
  
  
# This code maps between the PKA full length sequence and the Pfam protein kinase domain

col1 <- sort(ss_table[as.numeric(ss_table[,3]) %in% int_DF[[1]],1])
col2 <- sort(as.numeric(int_DF[[1]])-1)

# Pfam domain map
library(seqinr)
prkaca_fasta <- read.fasta('prkaca_hits.fasta',seqtype='AA')
prkaca_fasta <- prkaca_fasta[[2]]

start <- as.numeric(unlist(strsplit(unlist(strsplit(attr(prkaca_fasta,'name'),split="/"))[2],split="-"))[1]) - 1

sequence <- getSequence(prkaca_fasta)

count_prkaca <- start
count_domain <- 0

prkaca_vec <- NULL
domain_vec <- NULL

for (i in 1:length(sequence)) {
  
  if (sequence[i] == "-") {
    
    count_prkaca <- count_prkaca
    count_domain <- count_domain + 1
    
    
    prkaca_vec <- c(prkaca_vec,"-")
    domain_vec <- c(domain_vec,count_domain)
    
  } else if (sequence[i] == tolower(sequence[i])) {
    
    count_prkaca <- count_prkaca + 1
    count_domain <- count_domain 
    
    prkaca_vec <- c(prkaca_vec,count_prkaca)
    domain_vec <- c(domain_vec,"-")
    
    
  } else {
    
    count_prkaca <- count_prkaca + 1
    count_domain <- count_domain + 1
    
    prkaca_vec <- c(prkaca_vec,count_prkaca)
    domain_vec <- c(domain_vec,count_domain)
    
  }
}

pfam_map <- cbind(as.numeric(prkaca_vec),as.numeric(domain_vec))

col3 <- pfam_map[pfam_map[,1] %in% int_DF[[2]],2]

SDRs <- cbind(col1,col2,col3)
