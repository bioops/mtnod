# set work directory
setwd("~/github/mtnod/anova")

# load packages
library(dplyr)
library(tidyr)
library(igraph)
library(glasso)
library(networkD3)
library(htmlwidgets)
library(topGO)
#library(preprocessCore)

# 1. read gene expression files
mtgea.expr<-NULL
files<-dir(path="./expr/all/",full.names=T)
for (file in files){
  tmp<-read.table(file,header=T,sep="\t",stringsAsFactors = F)
  rownames(tmp)<-tmp[,1]
  tmp<-tmp[,-1]
  tmp<-log(t(tmp))
  tmp.names<-rownames(tmp)
  tmp.names<-gsub("_\\d$","",tmp.names, perl=T)
  tmp<-data.frame(condition=tmp.names,tmp)
  mtgea.expr<-rbind(mtgea.expr,tmp)
}
########

# 2. quantile normalize, not necessary
#mtgea.temp<-t(data.frame(normalize.quantiles(as.matrix(t(mtgea.expr[,-1])))))
#mtgea.temp<-data.frame(condition=mtgea.expr[,1],mtgea.temp)
#colnames(mtgea.temp)<-colnames(mtgea.expr)
#colnames(mtgea.temp)<-colnames(mtgea.expr)
#mtgea.expr<-mtgea.temp
########

# 3. 15 pre-selected probesets
start<-read.table("./input/start.blast.edit.tab",header=F,stringsAsFactors = F)
########

# 4. anova model (very slow, use pre calculated data instead via "load("fish.RData")")
fish<-NULL
for (gene in start$V2){
  pvals.1<-NULL
  pvals.2<-NULL
  pvals.3<-NULL
  for (i in 2:ncol(mtgea.expr)){
    long<-cbind(mtgea.expr[,c(1,i)],mtgea.expr[,gene]) %>% gather(gene,expr, -c(1))
    long.sum<-summary(aov(expr~condition*gene,data=long))
    pvals.1<-c(pvals.1,unlist(long.sum)[17])
    pvals.2<-c(pvals.2,unlist(long.sum)[18])
    pvals.3<-c(pvals.3,unlist(long.sum)[19])
  }
  fish<-c(fish,(which(pvals.1<0.01/length(pvals.1) & pvals.3>0.1))+1)
}
fish<-unique(fish)
########


# 5. map probesets to Mt4.0 genes
mtgea.select<-mtgea.expr[,fish]
mtgea.cor<-cor(mtgea.select)
probe.select<-colnames(mtgea.select)
map<-read.csv("./input/mapping4.csv",stringsAsFactors = F)
map.select<-map[map$probe_set %in% probe.select,]
map.select$gene<-gsub(".\\d$","",map.select$gene_id, perl=T)
map.agg<-aggregate(map.select, by=list(map.select$probe_set),min)
rownames(map.agg)<-map.agg$probe_set
mtgea.agg<-mtgea.select[,probe.select %in% map.agg$probe_set]
mtgea.agg<-t(mtgea.agg)
mtgea.final.t<-aggregate(mtgea.agg, by=list(map.agg[rownames(mtgea.agg),]$gene_id),mean)
rownames(mtgea.final.t)<-mtgea.final.t[,1]
mtgea.final.t<-mtgea.final.t[,-1]
mtgea.final<-t(mtgea.final.t)
########

# 5. glasso, network reconstruction
mtgea.final.cor<-cor(mtgea.final)
genes<-colnames(mtgea.final)
p<-ncol(mtgea.final)
n<-nrow(mtgea.final)
s<-mtgea.final.cor

# tune parameter
rholist=seq(from=0.4,to=0.9,by=0.02)

PLMs<-NULL
for(i in 1:(length(rholist))){
  g<-glasso(s,rholist[i])
  mat <- g$wi
  diag(mat) <- 0
  mat[which(mat != 0)] <- 1
  degree <- colSums(mat)
  power.lm.data<-data.frame(x=log(as.numeric(names(table(degree+1)))),y=log(as.numeric(table(degree+1))))
  PLM<-summary(lm(y~x,data=power.lm.data))$r.squared
  PLMs<-c(PLMs,PLM)
  print(rholist[i])
}

# figure 1
plot(rholist,PLMs, xlab=expression(rho),ylab=expression(R^2))
abline(h=0.8,v=0.82,col="gray")


# optimal glasso model
rho<-0.82
g<-glasso(s,rho)
mat <- g$wi
colnames(mat)<-genes
rownames(mat)<-genes
diag(mat) <- 0
mat[which(mat != 0)] <- 1
degree <- colSums(mat)
power.lm.data<-data.frame(x=log(as.numeric(names(table(degree+1)))),y=log(as.numeric(table(degree+1))))
expr.graph <- graph.adjacency(mat, mode="undirected")

# figure 3
f2lm<-lm(y~x,data=power.lm.data)
plot(power.lm.data,xlab="log(Degree+1)",ylab="Log(Node Number)")
abline(f2lm, col="gray")
legend("topright",expression(R^2==0.8417),box.col="gray")


# plot network
plot(expr.graph,vertex.label=NA,vertex.size=1)
########

# 6. module detection
g<-expr.graph
seed<-1
fc <- edge.betweenness.community(g)
mem<-fc$membership
mem[mem %in% c(5:9,11:91)]<-0
modules<-as.numeric(names(table(mem)))
V(g)$color <- mem
g$layout <- layout.fruchterman.reingold
palette(c("lightskyblue","red","green3","black","cyan",
          "palevioletred","yellow","gray","yellow","yellow","yellow"))

# figure 2
plot(g, vertex.label=NA,vertex.size=2)
########

# 7. plot expression pattern
# figure 4
par(oma=c(2,2,2,2),mar=c(6,5,5,5))
for(i in modules){
  expr<-scale(mtgea.final[,which(mem==i)])
  expr.wide1<-matrix(expr[,1],ncol=3,byrow =T)
  expr.mean1<-apply(expr.wide1,1,median)
  plot(expr.mean1,ylab="Relative Expression",xlab="",
       type="l",ylim=c(-2,2), xaxt = "n")
  axis(1, at=1:9, labels=FALSE)
  text(1:9, par("usr")[3] - 0.03*(par()$usr[4]-par()$usr[3]),
       labels = c("ART_00dpi", "ART_00mer", "Nod_03dpi",
                  "Nod_04dpi", "Nod_06dpi", "Nod_10dpi",
                  "Nod_14dpi", "Nod_20dpi","Nod_28dpi"),
       srt = 45, adj = 1, xpd = TRUE)
  title(xlab="Time Points",mgp=c(5.2,1,0))
  for (i in 2:ncol(expr)){
    expr.wide<-matrix(expr[,i],ncol=3,byrow =T)
    expr.mean<-apply(expr.wide,1,median)
    lines(expr.mean)
  }
}

########

# 8. output D3.js html
MyClickScript <- 
  '      d3.select(this).select("circle").transition()
        .duration(750)
.attr("r", 30)'

Links<-data.frame(get.edgelist(g,names=F))-1
colnames(Links)<-c("source","target")
Nodes<-data.frame(name=V(g)$name,group=mem)
forceNetwork(Links = Links , Nodes = Nodes, Source = "source", 
             Target = "target", NodeID = "name", 
             Group = "group", charge=-50, linkDistance=5, colourScale = JS("d3.scale.category10()")) %>%
  saveNetwork(file = 'medicago_net.html')
########


# 9. output all gene and annotation

gene.mod<-data.frame(gene=igraph::V(g)$name,module=mem, stringsAsFactors=F)
gene.mod$degree<-igraph::degree(g)

annot<-read.table("./input/Mt4.0v1_GenesProteinSeq_20130731_1800.annot.tab",
                  header=F,sep="\t",stringsAsFactors=F,quote="")
rownames(annot)<-annot[,1]

gene.mod$description<-annot[gene.mod$gene,2]

GO<-read.table("./input/Mt4.0v1_IPR2GO_20131219_1700.tsv",
               header=F,sep="\t",stringsAsFactors=F,quote ="")
rownames(GO)<-GO[,1]
gene.mod$GO<-GO[gene.mod$gene,2]

# Mt4.0 genes on Uniprot
mtgene.all<-read.table("./input/Mt4.0v1_GenesProteinSeq_20130731_1800.uniprot.select.tab",
                       header=F,sep="\t",stringsAsFactors=F,quote="")
# gene description
uniprot.all<-read.table("./input/uniprot.head",
                        header=F,sep="\t",stringsAsFactors=F,quote="")
# gene symbol
uniprot.gn<-read.table("./input/uniprot.head2",
                       header=F,sep="\t",stringsAsFactors=F,quote="")

rownames(uniprot.all)<-uniprot.all[,1]
mtgene.all$description<-uniprot.all[mtgene.all$V2,2]
rownames(uniprot.gn)<-uniprot.gn[,1]
mtgene.all$gene_name<-uniprot.gn[mtgene.all$V2,2]

rownames(mtgene.all)<-mtgene.all[,1]
gene.mod$uniprot<-mtgene.all[gene.mod$gene,2]
gene.mod$gene_name<-mtgene.all[gene.mod$gene,]$gene_name
gene.mod$description<-mtgene.all[gene.mod$gene,]$description

gene.mod$description_Mt40<-annot[gene.mod$gene,2]

mtgene<-data.frame(gene_id=V(g)$name)
map.net<-map.agg[map.agg$gene_id %in% V(g)$name,c(2,3)]
map.start<-data.frame(uniprot=start$V1,probe_set=start$V2)
map.merge<-merge(map.net,map.start,by="probe_set",all.x=T)
map.merge2<-merge(mtgene,map.merge,by="gene_id",all.x=T)
map.merge2$gene_id<-as.vector(map.merge2$gene_id)
map.merge2$uniprot<-as.vector(map.merge2$uniprot)
map.merge.agg<-aggregate(map.merge2, by=list(map.merge2$gene_id),min)
gene.mod$preselected<-map.merge.agg$uniprot

write.csv(gene.mod,"./output/gene.nod.all.csv")

## negative regulated genes in module 4
negative<-NULL
expr<-scale(mtgea.final[,which(mem==4)])
for (i in 1:ncol(expr)){
  expr.wide<-matrix(expr[,i],ncol=3,byrow =T)
  expr.mean<-apply(expr.wide,1,median)
  negative<-cbind(negative, expr.mean)
}
colnames(negative)<-colnames(expr)
negative.gene<-names(which(negative[2,]<0))
neg.gene.annot<-gene.mod[gene.mod$gene %in% negative.gene,]

########

# 10. GO enrichment test (F, P, C)
geneID2GO <- readMappings(file = "./input/gene.all.goC")
geneNames <- names(geneID2GO)

gene.go<-data.frame(Term=NULL,module=NULL)
for (i in modules[-1]){
  sig.gene<-gene.mod$uniprot[gene.mod$module==i]
  sig.gene<-gsub("^\\w{2}\\|","",sig.gene, perl=T)
  sig.gene<-gsub("\\|.*$","",sig.gene, perl=T)
  
  geneList <- factor(as.integer(geneNames %in% sig.gene))
  names(geneList) <- geneNames
  
  # ontology="MF", "BP", or "CC"
  GOdata <- new("topGOdata", ontology = "CC", allGenes = geneList,
                annot = annFUN.gene2GO, gene2GO = geneID2GO)
  
  resultFisher <- runTest(GOdata, algorithm = "weight01", statistic = "fisher")
  resultFisher
  pvalFis <- score(resultFisher)
  sig.term<-Term(names(pvalFis[pvalFis<0.01]))
  if(length(sig.term)>0){
    gene.go<-rbind(gene.go,data.frame(module=i,term=names(sig.term),pvalue=pvalFis[names(sig.term)],description=sig.term))
  }
}

gene.go
write.csv(gene.go,"./output/gene.nod.goC.csv")
########


# all done!
########




##########################################################
##                                                      ##
##  The following code is redundant, but maybe useful.  ##
##                                                      ##
##########################################################
# remove singleton
g<-delete.vertices(expr.graph, which(igraph::degree(expr.graph) < 1))
g<-delete.vertices(g,c("Medtr2g064310.1","Medtr4g060460.1"))
fc <- edge.betweenness.community(g)
com<-community.to.membership(g, fc$merges, steps= which.max(fc$modularity)-1)
mem<-com$membership+1
mem[mem>=9]<-0
V(g)$color <- mem
g$layout <- layout.fruchterman.reingold
plot(g, vertex.label=NA,vertex.size=5)

############


# plot network with 11 pre-selected gene (only 9 left)
V(g)$name<-map.merge.agg$uniprot
plot(g,vertex.size=5, vertex.label.color = "black")


pdf("crop1.pdf",width=10,height=10)
E(g)$weight=5
E(g)$color="black"
plot(g,vertex.size=3, vertex.label.color = "black")
dev.off()

colors<-c(palette(),"white")
gene.mod<-data.frame(gene=NULL,module=NULL)
for (i in 1:max(fc$membership)){
  tmp.mod<-data.frame(gene=fc$names[fc$membership==i],module=i)
  gene.mod<-rbind(gene.mod,tmp.mod)
}
gene.mod$gene<-as.vector(gene.mod$gene)

annot<-read.table("./input/Mt4.0v1_GenesProteinSeq_20130731_1800.annot.tab",
                  header=F,sep="\t",stringsAsFactors=F,quote="")
rownames(annot)<-annot[,1]
gene.mod$description<-annot[gene.mod$gene,2]

GO<-read.table("./input/Mt4.0v1_IPR2GO_20131219_1700.tsv",
               header=F,sep="\t",stringsAsFactors=F,quote ="")
rownames(GO)<-GO[,1]
gene.mod$GO<-GO[gene.mod$gene,2]


write.csv(gene.mod,"./output/gene.mod.csv")
#write.graph(expr.graph)

gene.all<-read.csv("./input/gene.mod.all.csv",stringsAsFactors = F)
mtgene.all<-read.table("Mt4.0v1_GenesProteinSeq_20130731_1800.uniprot.select.tab",
                       header=F,sep="\t",stringsAsFactors=F,quote="")
uniprot.all<-read.table("uniprot.head",
                       header=F,sep="\t",stringsAsFactors=F,quote="")
uniprot.gn<-read.table("uniprot.head2",
                       header=F,sep="\t",stringsAsFactors=F,quote="")

rownames(uniprot.all)<-uniprot.all[,1]
mtgene.all$description<-uniprot.all[mtgene.all$V2,2]
rownames(uniprot.gn)<-uniprot.gn[,1]
mtgene.all$gene_name<-uniprot.gn[mtgene.all$V2,2]

rownames(mtgene.all)<-mtgene.all[,1]
gene.all$uniprot<-mtgene.all[gene.all$gene,2]
gene.all$gene_name<-mtgene.all[gene.all$gene,]$gene_name
gene.all$description<-mtgene.all[gene.all$gene,]$description

gene.all$description_Mt40<-annot[gene.all$gene,2]
gene.all$GO_Mt40<-GO[gene.all$gene,2]

write.csv(gene.all,"./output/gene.nod.all.csv")


# go enrichment


test<-mtgea.expr[,1:20] %>% gather(gene,expr, -c(1))

test.sum<-summary(aov(expr~condition*gene,data=test))
unlist(test.sum)[17]


RT_0dpi<-read.table("expr/RT_0dpi.txt",header=T,sep="\t",stringsAsFactors = F)

start<-read.table("start.blast.edit.tab",header=F,stringsAsFactors = F)

# figure 1
par(mfrow=c(4,4))
for (i in 1:15){
  long<-data.frame(condition=mtgea.expr[,1],expr=mtgea.expr[,start$V2[i]])
  boxplot(expr~condition,data=long,main=paste(start$V1[i],"-",start$V2[i]),xaxt="n")
  axis(1, at=1:9, labels=c("RM",0,3,4,6,10,14,20,28))
}


par(mfrow=c(4,4),mar=c(2,2,2,2))
for (i in 1:16){
  long<-data.frame(condition=mtgea.expr[,1],expr=mtgea.select[,i])
  boxplot(expr~condition,data=long,main=start$V2[i])
}

long<-data.frame(condition=mtgea.expr[,1], mtgea.expr[,start$V2[3:4]]) %>% gather(gene,expr, -c(1))
summary(aov(expr~condition*gene,data=long))

tmp.long<-tmp %>% gather(gene,expr, -c(1)) 

tmp2<-read.table("expr/Nod_4dpi.txt",header=T,sep="\t",stringsAsFactors = F)
rownames(tmp2)<-tmp2[,1]
tmp2<-tmp2[,-1]
tmp2<-t(tmp2)
tmp.names<-rownames(tmp2)
tmp.names<-gsub("_\\d$","",tmp.names, perl=T)
tmp2<-data.frame(condition=tmp.names,tmp2)
tmp2.long<-tmp2 %>% gather(gene,expr, -c(1)) 

