### cell type classification

####### 1.construct weighted graph #########
## function to construct weighted graph based on blocks file from biclustering.
## Input is .blocks file, output is three-column weighted graph, with nodes being cells and edge being weight
## this is an intermediate function, the only purpose is prepare weighted graph for following clustering 

GRAPH <-function(blocks){
  F <-readLines(blocks)
  TEMP <-grep('Conds',F,value=T) ## extract condition lines in each BC
  BC <-sapply(strsplit(TEMP,':',2),'[',2) # only keep cell names
  
  CONDS <-as.character()   # store the conditions 
  label_C <-as.numeric()   # store the occurence of one condistions
  
  for (j in 1:length(BC)){
    BCcond <-unlist(strsplit(BC[j], split = " "))
    BCcond <-BCcond[BCcond!=""]  # exclude the blank string
    CONDS <-c(BCcond,CONDS)
    label_C <-c(label_C,rep(j,length(BCcond)))
  }
  
  df_C <-data.frame(conds=CONDS,label=label_C)
  uniq_C <-df_C$conds[!duplicated(df_C$conds)]   # unique conditions
  Node <-t(combn(uniq_C,2))
  
  Wt <-rep(-1,dim(Node)[1])
  for (k in 1:dim(Node)[1]){
    member1 <-df_C[which(df_C$conds %in% Node[k,1]),]   # identify which BC the k th Node appear
    member2 <-df_C[which(df_C$conds %in% Node[k,2]),]
    Wt[k] <-length(intersect(member1[,2],member2[,2])) # the weight between two node
  }
  GRAPH <-data.frame(Node[,1],Node[,2],Wt)
  names(GRAPH) <-c('Node1','Node2','Weigth')
  if (dim(GRAPH)[1]!=0)	{
    write.csv(subset(GRAPH,Wt!=0),paste(basename(blocks),'_graph.csv',sep=''),row.names=FALSE)
  }
}

####### 2. cell type prediction  #######
## cell type prediction based on weighted graph

library(igraph)
library(mclust)
library(MCL)
library(clues)
library(anocva)

## clustering function 
MCL <-function(Raw,graph){   # Raw is the original expression matrix
  RAW <-read.table(Raw,header=T,sep='\t')
  CellNum <-dim(RAW)[2]-1  # the number of cells
  Graph <-read.csv(graph,header=T,sep=",") 
  G <-graph.data.frame(Graph,directed = FALSE)  # convert file into graph
  A <- as_adjacency_matrix(G,type="both",attr="Weight",names=TRUE,sparse=FALSE)  # convert graph into adjacency matrix
  V_name <-rownames(A)   # the vertix
  Covered <-length(V_name)  # the #of covered cells
  
  CLUST <-list()
  for (i in 1:100){
    CLUST[[i]] <-mcl(A,addLoops = FALSE,inflation =i,max.iter=200)
  }
  KK <- as.data.frame(do.call(rbind,lapply(CLUST,'[[',1)))  # extract the number of clusters
  CAN_I <-c(which(as.numeric(as.character(KK$V1))>=2)) 	# results that has more than 5 clusters
  tt <-as.numeric(as.character(KK$V1))
  tt <-sort(table(tt),decreasing=T)[1]
  Final_K <-as.numeric(names(tt))
  
  if (length(CAN_I)!=0){
    MATRIX <-rep(0,Covered)%o%rep(0,Covered)
    for (k in 1:length(CAN_I)){	
      MCL_label <-CLUST[[CAN_I[k]]]$Cluster  # record the label
      ClusterNum <-unique(MCL_label)   # record the number of clusters
      TEMP <-rep(0,Covered)%o%rep(0,Covered)
      temp <-rep(0,Covered) %o% rep(0,length(ClusterNum))
      for (n in 1:length(ClusterNum)){
        index <-which(MCL_label==ClusterNum[n])
        temp[index,n] <-1
        TEMP <-TEMP+temp[,n]%o%temp[,n] 
      }
      MATRIX <-MATRIX+TEMP
    }
    MATRIX <-MATRIX/length(CAN_I)
    rownames(MATRIX) <-colnames(MATRIX) <-rownames(A)
    hc <-hclust(dist(MATRIX))
    memb <-cutree(hc,k=Final_K)
    if (length(rownames(A)) ==CellNum){
      label <-memb
    }else{
      LEFT <-setdiff(names(RAW)[-1],V_name)
      LEFT_Cluster <-rep(Final_K+1,length(LEFT))
      df_cell_label <-data.frame(cell=c(names(memb),LEFT),cluster=c(memb,LEFT_Cluster),K=rep(Final_K+1,CellNum))				
      label <-df_cell_label$cluster
    }	
  }
  return(label)
}


SC <-function(Raw,graph,K){
  RAW <-read.table(Raw,header=T,sep='\t')  # expression data
  CellNum <-dim(RAW)[2]-1  # the number of cells 
  Graph <-read.csv(graph,header=T,sep=",") 
  G <-graph.data.frame(Graph,directed = FALSE)  # convert file into graph
  A <- as_adjacency_matrix(G,type="both",attr="Weight",names=TRUE,sparse=FALSE)  # convert graph into adjacency matrix
  V_name <-rownames(A)   # the vertix
  Covered <-length(V_name)  # the #of covered cells
  
  sc <-spectralClustering(A,k=K)
  names(sc) <-rownames(A)
  if (length(rownames(A)) ==CellNum){
    label <-sc
  }else{
    LEFT <-setdiff(names(RAW)[-1],V_name)
    LEFT_Cluster <-rep(K+1,length(LEFT))
    df_cell_label <-data.frame(cell=c(names(sc),LEFT),cluster=c(sc,LEFT_Cluster),K=rep(K+1,CellNum))				
    label <-df_cell_label$cluster
  }
  return(label)
}

## this should be the function user really use, output cell label (always) and ARI score (only when reference label is provided)
## Raw is the path to the original expression matrix
## graph is the weighted graph from GRAPH function
## method should be either 'MCL' or 'SC', and if 'SC', user also need to specify k
## ref is optional input, which is reference label. Must have two columns, one is 'Cell_type', denoting the name of cells, the other is 'Cluster', denoting the membership of each cell
CLUSTERING <- function(Raw,graph,method='MCL',K=null,ref='Yan_cell_label.csv'){
  RST <-list()
  if (method=='MCL'){
    RST[['Label']] <-MCL(Raw,graph)
  }else if (method =='SC'){
    RST[['Label']] <-SC(Raw,graph,k)
  }
  
  ### I want to add a judgment here, if user provide reference label, proceed to calculate ARI, FM and JI
  ### if not, just return Label	
  if (!is.null(ref)){
    target <-read.table(ref,header=T,sep=',')
    aa <-names(RST[['Label']])
    bb <-target$Cell_type
    # judge if the cell names are consistent
    # if consistent, continue to calculate ARI ect
    if (identical(sort(aa),sort(as.character(bb)))=='TRUE'){
      sorted <-RST[['Label']][match(target$Cell_type,names(RST[['Label']]))] # sort the predicted label
      ARI <-adjustedRandIndex(sorted,target$Cluster)  
      RI <-adjustedRand(sorted,target$Cluster,randMethod='Rand')
      FM <-adjustedRand(sorted,target$Cluster,randMethod='FM')
      JI <-adjustedRand(sorted,target$Cluster,randMethod='Jaccard')
      df <-data.frame(ARI=ARI, RandIndex=RI,FolkesMallow=FM, Jaccard=JI)
      
      RST[['df']] <-df
    }else{
      warning(paste0('Cell names in the reference label and predicted label are inconsistent, please double check !'))
    }
  }
  RST
}

#' @export
## final function
## i is the input, K is an optional parameter, used only when method=='SC'
final <- function(i, method = 'MCL', K, ref, R = FALSE, F = FALSE, d = FALSE, f = 0.85, k = 13, c = 0.90, o = 5000){
  qubic(i, R = R, F = F, d = d, f = f, k = k, c = c, o = o)
  if (R) {
    GRAPH(paste0(i,'.chars.blocks'))
    CLUSTERING(i, paste0(i,'.chars.blocks_graph.csv'), method, K = K, ref)    # not sure how to deal with that K 
  } else {
    GRAPH(paste0(i,'.blocks'))
    CLUSTERING(i, paste0(i,'.blocks_graph.csv'), method, K = K, ref)    # not sure how to deal with that K 
  }
}
