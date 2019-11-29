setwd('/workdir')
setwd('/Users/wk/Downloads/')

input <- "./input/"
output <- "./output/"

library(glmnet)
load( paste0( input,"marker_all.Rdata" ) ) 
# -1为aa,1为AA,0为miss

marker <- marker_all
marker <- t(marker)
n <- dim(marker)[1]
marker_list <- list()
name <- colnames(marker)

cat("Constructing the items for each differential equation... ...\n")
cat("equation:\n")

for (col in 1:length(marker[1,])) {
  #tim <- proc.time()
  cat(col,"\t")
  if(col%%10==0){cat("\n")}
  if(col == length(marker[1,])){cat("\n")}
  m <- marker[,col]
  M <- marker[,-col]
  vec <- rep(NA,length(M[1,]))
  for (i in 1:length(M[1,])) {
    vec[i] <- cor(m,M[,i])
  }
  x_matrix <- M[,which( vec %in% -sort(-vec)[1:(n/log(n))] )]
  ridge1_cv <- cv.glmnet(x = x_matrix, y = m,type.measure = "mse",nfold = 10,alpha = 0)
  best_ridge_coef <- as.numeric(coef(ridge1_cv, s = ridge1_cv$lambda.min))[-1]
  fit_res <- cv.glmnet(x = x_matrix, y = m,type.measure = "mse",nfold = 10,alpha = 1,penalty.factor = 1 / abs(best_ridge_coef),keep = TRUE)
  best_alasso_coef1 <- coef(fit_res, s = fit_res$lambda.min)
  marker_list_one <- list()
  marker_list_one[[1]] <- name[col]#第一个列表是直接qtl的名字
  marker_list_one[[2]] <- best_alasso_coef1@Dimnames[[1]][best_alasso_coef1@i[-1]+1]#第二个列表是间接qtl的名字
  marker_list_one[[3]] <- best_alasso_coef1@x[-1]#第三个列表是变量选择系数
  marker_list[[col]] <- marker_list_one
  #proc.time() - tim
}


load( paste0( input,"effect_par_all.Rdata" ) )
# 前5个元素为aa基因型的平均参数
# 后5个元素为AA基因型的平均参数
# 5号,10号元素为最佳曲线方程编号


get_LOPm <- function(X){
  len = length(X)
  LOP <- function(r){
    t <- seq(-1,1,2/(len-1))
    temp <- rep(0,len)
    for (m in 0:as.integer(r/2)) {
      temp <- temp  + (-1)^m*gamma(2*r - 2*m + 1)/(2^r*gamma(m+1)*gamma(r-m+1)*gamma(r-2*m + 1)) * t^(r-2*m)
    }
    return(temp)
  }
  LOPm <- cbind(LOP(0),LOP(1),LOP(2),LOP(3),LOP(4),LOP(5),LOP(6) )
  return(LOPm[,1:ORDER])
}
ORDER <- 6
library(mvtnorm)
f1 <- function(x,t){
  y <- t[1]*exp(-exp(t[2]*exp(1)*(t[3]-x)/t[1] + 1))
  return(y)
}
f2 <- function(x,t){
  y <- t[1]/(1+exp(4*t[2]*(t[3]-x)/t[1] + 2))
  return(y)
}
f3 <- function(x,t){
  y <- t[1]*(1+t[4]*exp(1+t[4])*exp((t[2]/t[1])*(1+t[4])^(1+1/t[4])*(t[3]-x)))^(-1/t[4])
  return(y)
}
func <- list(f1,f2,f3)
fy <- function(t,X){
  func[[t[5]]](X,t[1:4][which(!is.na(t[1:4]))] )-func[[t[10]]](X,t[6:9][which(!is.na(t[6:9]))] )
}
get_origin <- function(dy,X,y0){
  y0 <- c(y0)
  for (i in 2:(length(X)-1)) {
    slope <- dy[i-1]
    y_before <- y0[length(y0)]
    add <- y_before + slope*(X[2]-X[1])
    y0 <- c(y0,add)
  }
  return(y0)
}
fl_new <- function(t,X,dep,ind,dep_per,ind_per,LOPm){
  ydep <- fy(c(dep_per),X)
  d <- 2*(ydep[-1] - ydep[-length(ydep)])/(X[2]-X[1])
  tm <- matrix(t,ncol=ORDER,byrow = T)
  temp1 <- LOPm[-1,]%*%t(tm) # mp * m个线
  temp1 <- temp1*matrix(rep(ydep[-1],length(ind)+1),ncol = length(ind)+1,byrow = F)
  for (i in 1:length(ind)) {
    yind <- fy( c(ind_per[i,]) ,X)
    temp1[,i+1] <- temp1[,i+1]*yind[-1]
  }
  x0 <- X[-length(X)]
  temp0 <- LOPm[-length(X),]%*%t(tm)
  temp0 <- temp0*matrix(rep(ydep[-length(X)],length(ind)+1),ncol = length(ind)+1,byrow = F)
  for (i in 1:length(ind)) {
    yind <- fy(c(ind_per[i,]),X)
    temp0[,i+1] <- temp0[,i+1]*yind[-length(X)]
  }
  #----------------
  d_mat <- LOPm%*%t(tm)
  ydep <- fy(c(dep_per),X) # effect
  for (i in 1:length(tm[,1])) {
    if (i == 1) {
      d_mat[,i] <- d_mat[,i]*ydep
    }else{
      d_mat[,i] <- d_mat[,i]*fy(c(ind_per[i-1,]),X)
    }
  }
  # 将d_mat转化为o_mat
  o_mat <- c()
  for (i in 1:length(d_mat[1,]) ) {
    if(i == 1){
      o_mat <- cbind(o_mat,get_origin(d_mat[,i],X,ydep[1]))
    }else{
      o_mat <- cbind(o_mat,get_origin(d_mat[,i],X,0))
    }
  }
  o_mat <- as.data.frame(o_mat)
  y <- ydep[-length(ydep)]
  e <- colSums(t(o_mat))
  ssr <- sum((y-e)^2)
  sst <- sum((y-mean(y))^2)
  r <- 1-(ssr/sst)
  #----------------
  return(sum((d - colSums(t(temp1 + temp0)) )^2) + abs(1-r))
}

for (col in 1:length(marker_list) ) {
  cat(col,'\n')
  dep <- marker_list[[col]][[1]]
  ind <- marker_list[[col]][[2]]
  if( length(ind) == 0 ){
    next
  }
  dep_per <- matrix(effect_par_all[[dep]],byrow = T)
  ind_per <- matrix(NA, ncol = 10,nrow = length(ind))
  for (i in 1:length(ind)) {
    ind_per[i,] <- effect_par_all[[ind[i]]]
  }
  
  X <- seq(0,50,50/((length(ind)+1)*ORDER*4))
  
  # 等式的右边有两部分,1~n/0~n-1
  t0 <- rep(0.001,(length(ind)+1)*ORDER)
  
  # ptm <- proc.time()
  itimes <- 1
  repeat{
    s1 <- optim(t0,fl_new,method = 'Nelder-Mead',X = X,dep = dep,ind = ind, 
                dep_per = dep_per, ind_per = ind_per, LOPm = get_LOPm(X))
    r1 <- s1$par
    s2 <- optim(r1,fl_new,method = 'Nelder-Mead',X = X,dep = dep,ind = ind, 
                dep_per = dep_per, ind_per = ind_per, LOPm = get_LOPm(X))
    cat(col,'-',itimes,s2$value,'\n')
    itimes <- itimes + 1
    if(all( abs(r1-s2$par) == 0 )||itimes == 10){ #*** itimes越高精度越高,计算速度越慢,有条件部署在集群时,应该尽可能大与1000 ***#
      break
    }else{
      t0 <- s2$par
    }
  }
  # roc.time() - ptm
  
  marker_list[[col]][[4]] <- matrix(s2$par,ncol=ORDER,byrow=TRUE)
  tm <-  matrix(s2$par,ncol=ORDER,byrow=TRUE)
  
  ydep <- fy(c(dep_per),X) # effect
  d_mat <- get_LOPm(X)%*%t(tm)
  
  for (i in 1:length(tm[,1])) {
    if (i == 1) {
      d_mat[,i] <- d_mat[,i]*ydep
    }else{
      d_mat[,i] <- d_mat[,i]*fy(c(ind_per[i-1,]),X)
    }
  }
  
  # 将d_mat转化为o_mat
  o_mat <- c()
  for (i in 1:length(d_mat[1,]) ) {
    if(i == 1){
      o_mat <- cbind(o_mat,get_origin(d_mat[,i],X,ydep[1]))
    }else{
      o_mat <- cbind(o_mat,get_origin(d_mat[,i],X,0))
    }
  }
  
  o_mat <- as.data.frame(o_mat)

  if( dim(o_mat)[2] <= 2 ){
    marker_list[[col]][[5]] <- sum(o_mat[,length(o_mat[1,])]*(X[2]-X[1]))
  }else{
    marker_list[[col]][[5]] <- colSums(o_mat[,2:(length(ind)+1)]*(X[2]-X[1]))
  }
}

# 将marker_list转化为after table
after <- c()
for (i in 1:length(marker_list)){
  dep <- marker_list[[i]][[1]]
  ind <- marker_list[[i]][[2]]
  effect <- marker_list[[i]][[5]]
  one <- c()
  for (j in 1:length(ind)) {
    if(effect[j] >= 0){
      type <- '+'
    }else{
      type <- '-'
    }
    one <- rbind(one,c(ind[j],dep,abs(effect[j]),type))
  }
  after <- rbind(after,one)
}

# 网络
origin <- after
nodes <- data.frame( names=names(table(c(origin[,1],origin[,2]))) )  
nodes <- cbind( id=paste0( 's',1:length( nodes[,1] ) ),nodes )

nodes$id <- as.character(nodes$id)
nodes$names <- as.character(nodes$names)

for (i in 1:length(origin[,1]) ){
  origin[i,1] <- nodes$id[ which( nodes$names==origin[i,1] ) ]
  origin[i,2] <- nodes$id[ which( nodes$names==origin[i,2] ) ]
}

origin[,c(3,4)] <- origin[,c(4,3)]

colnames(origin)<-c("from","to","type","weight")
origin <- as.data.frame(origin)
origin$from <- as.character(origin$from)
origin$to <- as.character(origin$to)
origin$type <- as.character(origin$type)
origin$weight <- abs(as.numeric(as.character( origin$weight )))

links <- origin
links$edge.color <- NA
for (i in 1:length( links[,1] )) {
  if(links$type[i] == "-"){
    links$edge.color[i]<-grDevices::adjustcolor("#377eb8", alpha=0.4)
  }
  if(links$type[i] == "+"){
    links$edge.color[i]<-grDevices::adjustcolor("#e41a1c", alpha=0.4)
  }
}

nodes$node.color <- "#fdae61"

for (i in 1:length(links[,1])) {
  links$weight[i] = (links$weight[i]-min(links$weight))/(max(links$weight)-min(links$weight))
}

library("igraph")
net <- graph_from_data_frame( d=links,vertices = nodes,directed = T )
set.seed(1)
coor_w1<-layout.auto(net)

#得到标准的坐标
stand_coor <- data.frame( names=nodes$names,coor_w1 )
stand_coor$names <- as.character( stand_coor$names )
stand_coor$X1 <- as.numeric( stand_coor$X1 )
stand_coor$X2 <- as.numeric( stand_coor$X2 )

setwd(output)
pdf(file = 'net.pdf',width = 10,height = 10)
plot(net,edge.arrow.size=1,
     vertex.label=V(net)$names,
     vertex.size=5,
     vertex.label.cex=0.55,
     edge.curved=0,
     edge.color=E(net)$edge.color,
     edge.frame.color=E(net)$edge.color,
     edge.width=E(net)$weight*10,
     vertex.color=V(net)$node.color,
     layout=coor_w1
)
dev.off()



