#307686155@qq.com
#18666983305
#更多课程请关注生信碱移
#af

library("clusterProfiler")
library("org.Hs.eg.db")
library("enrichplot")
library("ggplot2")
library("pathview")
library("ggnewscale")
library("DOSE")
library(stringr)

pvalueFilter=0.05     #过滤p值    
qvalueFilter=0.05     #过滤矫正后的p值，可改为1，其他富集同理    
shownum=7            #展示前多少条通路
highcol="gray"      #高颜色(P值更大)
lowcol="red3"          #低颜色

rt=read.table("Ferroptosis signal.txt",sep="\t",check.names=F,header=F)      
genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)  
entrezIDs <- as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)
colnames(rt)=c("symbol","entrezID") 
rt=rt[is.na(rt[,"entrezID"])==F,]                        
gene=rt$entrezID
gene=unique(gene)

kk=enrichGO(gene = gene,OrgDb = org.Hs.eg.db, pvalueCutoff =1, qvalueCutoff = 1, ont="all", readable =T)
GO=as.data.frame(kk)
GO=GO[(GO$pvalue<pvalueFilter & GO$qvalue<qvalueFilter),]
write.table(GO,file="GO.xls",sep="\t",quote=F,row.names = F)
GO=GO[c(which(GO$ONTOLOGY=="BP")[1:shownum],which(GO$ONTOLOGY=="CC")[1:shownum],which(GO$ONTOLOGY=="MF")[1:shownum]),]
names(GO)
rt=GO[,c(1,2,3,13,4,9,11)]       
names(rt)=c("Ontology","ID","Term","Count","Ratio","pvalue","qvalue")
#添加编号
rt$ID=gsub(":","",rt$ID)
for (i in 1:nrow(rt)) {
  rt[i,3]=paste0(rt[i,2],":",rt[i,3])
}
#分列
split_b<-str_split(rt$Ratio,"/")
b<-sapply(split_b,"[",1)
c<-sapply(split_b,"[",2)
rt$Ratio=as.numeric(rt$Count)/as.numeric(c[1])
#按照Ratio对Term排序
labels=rt[order(rt$Ratio),"Term"]
rt$Term = factor(rt$Term,levels=labels)
#绘制
p = ggplot(rt,aes(Ratio, Term)) + 
  geom_point(aes(size=Count, color=qvalue))
p1 = p + 
  scale_colour_gradient(high=highcol, low = lowcol) + #midpoint为颜色梯度中间节段点，low、mid、high分别为对应三颜色截断
  labs(color="qvalue",size="Count",x="Gene ratio",y="Term")+     
  theme(axis.text.x=element_text(color="black", size=10),axis.text.y=element_text(color="black", size=10)) + 
  scale_size_continuous(range=c(4,9))+      #控制点的大小范围，即4到9
  theme_bw()+
  facet_grid( Ontology~. ,scales="free")
ggsave("GO_bubble.pdf", width=7, height=6)      #保存

rt=GO[,c(1,2,3,13,4,9,11)]       
names(rt)=c("Ontology","ID","Term","Count","Ratio","pvalue","qvalue")
#添加编号
rt$ID=gsub(":","",rt$ID)
for (i in 1:nrow(rt)) {
  rt[i,3]=paste0(rt[i,2],":",rt[i,3])
}
#按FDR排序
labels=rt[order(rt$qvalue,decreasing =T),"Term"]
rt$Term = factor(rt$Term,levels=labels)
#绘制,midpoint为颜色梯度中间节段点
p=ggplot(data=rt)+geom_bar(aes(x=Term, y=Count, fill=qvalue), stat='identity')+
  coord_flip() + scale_fill_gradient(high=highcol, low = lowcol) +     #midpoint为颜色梯度中间节段点，low、mid、high分别为对应三颜色截断
  xlab("Term") + ylab("Gene count") +          #x、y轴名称
  theme(axis.text.x=element_text(color="black", size=10),axis.text.y=element_text(color="black", size=10)) + 
  scale_y_continuous(expand=c(0, 0)) + scale_x_discrete(expand=c(0,0))+
  theme_bw()+
  facet_grid( Ontology~. ,scales="free")
print(p)
ggsave("GO_barplot.pdf", width=7, height=6)        #保存图片

save.image(file="Go.Ddata")
#307686155@qq.com
#18666983305
#更多课程请关注生信碱移
#af
