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
showNum=20            #展示前多少条通路
highcol="grey"      #高颜色(P值更小)
lowcol="red3"          #低颜色

#输入文件
rt=read.table("Ferroptosis signal.txt",sep="\t",check.names=F,header=F)      
genes=as.vector(rt[,1])
entrezIDs <- mget(genes, org.Hs.egSYMBOL2EG, ifnotfound=NA)  
entrezIDs <- as.character(entrezIDs)
rt=cbind(rt,entrezID=entrezIDs)
colnames(rt)=c("symbol","entrezID") 
rt=rt[is.na(rt[,"entrezID"])==F,]                        
gene=rt$entrezID
gene=unique(gene)
R.utils::setOption( "clusterProfiler.download.method",'auto' )

kk <- enrichKEGG(gene = gene, organism = "hsa", pvalueCutoff =1, qvalueCutoff =1)
KEGG=as.data.frame(kk)
KEGG$geneID=as.character(sapply(KEGG$geneID,function(x)paste(rt$symbol[match(strsplit(x,"/")[[1]],as.character(rt$entrezID))],collapse="/")))
KEGG=KEGG[(KEGG$pvalue<pvalueFilter & KEGG$qvalue<qvalueFilter),]
write.table(KEGG,file="KEGG.xls",sep="\t",quote=F,row.names = F)

rt=KEGG[1:showNum,c(3,2,14,5,10,11)]       
names(rt)=c("ID","Term","Count","Ratio","pvalue","qvalue")
#添加编号
for (i in 1:nrow(rt)) {
  rt[i,2]=paste0(rt[i,1],":",rt[i,2])
}
#分列
split_b<-str_split(rt$Ratio,"/")
b<-sapply(split_b,"[",1)
c<-sapply(split_b,"[",2)
rt$Ratio=as.numeric(rt$Count)/as.numeric(c[1])
#按照Ratio对Term排序
labels=rt[order(rt$Ratio),"Term"]
duplicated(labels)  # 返回逻辑向量，显示哪些元素重复
any(duplicated(labels))  # 返回TRUE表示存在重复
rt$Term = factor(rt$Term,levels=labels)
#绘制
p = ggplot(rt,aes(Ratio, Term)) + 
  geom_point(aes(size=Count, color=qvalue))
p1 = p + 
  scale_colour_gradient(high=highcol, low = lowcol) + #midpoint为颜色梯度中间节段点，low、mid、high分别为对应三颜色截断
  labs(color="qvalue",size="Count",x="Gene ratio",y="Term")+     
  theme(axis.text.x=element_text(color="black", size=10),axis.text.y=element_text(color="black", size=10)) + 
  scale_size_continuous(range=c(4,9))+      #控制点的大小范围，即4到9
  theme_bw()
ggsave("KEGG_bubble.pdf", width=7, height=6)      #保存

rt=KEGG[1:showNum,c(3,2,14,5,10,11)]          
names(rt)=c("ID","Term","Count","Ratio","pvalue","qvalue")
#添加编号
for (i in 1:nrow(rt)) {
  rt[i,2]=paste0(rt[i,1],":",rt[i,2])
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
  theme_bw()
print(p)
ggsave("KEGG_barplot.pdf", width=7, height=6)        #保存图片

save.image(file= "kegg.Rdata" )
#307686155@qq.com
#18666983305
#更多课程请关注生信碱移
#af
