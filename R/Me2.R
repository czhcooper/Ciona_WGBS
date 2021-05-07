library(dplyr)
library(tidyverse)
# prase HT.gtf file
read.table("~/bio/ciona_i/genome/HT/Genemodels/corrected_HT.gtf",sep = "\t",header = F,quote = "",skip = 1)->HT.gtf
read.table("~/bio/ciona_i/genome/HT/Genemodels/HT.exon.gtf",sep = "\t",quote = "",header = F)->HT.exon.gtf

#select the longest transcript from each gene
HT.gtf %>% separate(V9,into = c("transcript_id","gene_id"),sep = ";")->HT.gtf
sub("transcript_id","",HT.gtf$transcript_id)->HT.gtf$transcript_id
sub(" ","",HT.gtf$transcript_id)->HT.gtf$transcript_id
sub(" ","",HT.gtf$transcript_id)->HT.gtf$transcript_id
sub("\"", "",HT.gtf$transcript_id)->HT.gtf$transcript_id
sub("\"", "",HT.gtf$transcript_id)->HT.gtf$transcript_id

sub("gene_id","",HT.gtf$gene_id)->HT.gtf$gene_id
sub(" ","",HT.gtf$gene_id)->HT.gtf$gene_id
sub(" ","",HT.gtf$gene_id)->HT.gtf$gene_id
sub("\"", "",HT.gtf$gene_id)->HT.gtf$gene_id
sub("\"", "",HT.gtf$gene_id)->HT.gtf$gene_id

Longest.mRNA<-HT.gtf %>% filter(V3=="mRNA") %>% group_by(gene_id) %>% summarise(Len=V5-V4,transcript_id=transcript_id) %>% group_by(gene_id) %>% filter(Len==max(Len))
# some gene might have the mRNA wiht same length, so only select one of them.
Longest.mRNA<-Longest.mRNA[!duplicated(Longest.mRNA$gene_id),]


#prase HT.exon.gtf
sub("gene_id","",HT.exon.gtf$V9)->HT.exon.gtf$V9
sub("transcript_id","",HT.exon.gtf$V9)->HT.exon.gtf$V9
sub(" ","",HT.exon.gtf$V9)->HT.exon.gtf$V9
sub(" ","",HT.exon.gtf$V9)->HT.exon.gtf$V9
sub("\"","",HT.exon.gtf$V9)->HT.exon.gtf$V9
sub("\"","",HT.exon.gtf$V9)->HT.exon.gtf$V9
sub(";","",HT.exon.gtf$V9)->HT.exon.gtf$V9
sub("KY2019:","",HT.exon.gtf$V9)->HT.exon.gtf$V9
#select the exon of each longest mRNA
HT.gene2exon<-HT.exon.gtf %>% filter(V9 %in% Longest.mRNA$transcript_id)

HT.gene2exon$geneID<-Longest.mRNA$gene_id[match(HT.gene2exon$V9,Longest.mRNA$transcript_id)]
HT.gene2exon %>% select(V1,V3,V4,V5,V7,V8,geneID)->HT.gene2exon
colnames(HT.gene2exon)<-c("Chr","type","Start","End","Strand","Score","geneID")
HT.gene2exon.list<-split(HT.gene2exon,HT.gene2exon$geneID)

#annotate intron
HT.gene2exon2intron.list<-lapply(HT.gene2exon.list,function(x){
  y<-dim(x)[1]
  x<-x %>% arrange(Start)
  if(y==1){
    return(x)
  } else {for (i in 1:(y-1) ){
    x[y+i,"Chr"]<-x[1,"Chr"]
    x[y+i,"type"]<-"intron"
    x[y+i,"Start"]<-x[i,"End"]+1
    x[y+i,"End"]<-x[i+1,"Start"]-1
    x[y+i,"Strand"]<-x[1,"Strand"]
    x[y+i,"Score"]<-x[1,"Score"]
    x[y+i,"geneID"]<-x[1,"geneID"]
  }
    return(x)}

})

#annotate promoter and downstream
HT.gene2exon2intron2promoter2downstream.list<-lapply(HT.gene2exon2intron.list,function(x){
  y<-dim(x)[1]
  if(x[1,"Strand"]=="+"){
    exon<-x %>% filter(type=="exon") %>% arrange(Start)
    exon1<-exon[1,]
    if(exon1[1,"Start"]>=1000){
      #promoter 1kb
      x[y+1,"Chr"]<-x[1,"Chr"]
      x[y+1,"type"]<-"promoter1"
      x[y+1,"Start"]<-exon1[1,"Start"]-1000+1
      x[y+1,"End"]<-exon1[1,"Start"]-1
      x[y+1,"Strand"]<-x[1,"Strand"]
      x[y+1,"Score"]<-x[1,"Score"]
      x[y+1,"geneID"]<-x[1,"geneID"]

      #promoter 2kb
      if (exon1[1,"Start"]>=2000){
        x[y+2,"Chr"]<-x[1,"Chr"]
        x[y+2,"type"]<-"promoter2"
        x[y+2,"Start"]<-exon1[1,"Start"]-2000+1
        x[y+2,"End"]<-exon1[1,"Start"]-1000-1
        x[y+2,"Strand"]<-x[1,"Strand"]
        x[y+2,"Score"]<-x[1,"Score"]
        x[y+2,"geneID"]<-x[1,"geneID"]
      } else {}

    } else{ }
   # downstream 1kb
    exon_1<-exon[dim(exon)[1],]
    if((HT.size$V2[HT.size$V1==exon_1[1,"Chr"]]- exon_1[1,"End"]) >= 1000){
      x[y+3,"Chr"]<-x[1,"Chr"]
      x[y+3,"type"]<-"downstream1"
      x[y+3,"Start"]<-exon_1[1,"End"]+1
      x[y+3,"End"]<-exon_1[1,"End"]+1000
      x[y+3,"Strand"]<-x[1,"Strand"]
      x[y+3,"Score"]<-x[1,"Score"]
      x[y+3,"geneID"]<-x[1,"geneID"]
      #downstream 2kb
      if((HT.size$V2[HT.size$V1==exon_1[1,"Chr"]]- exon_1[1,"End"]) >= 2000){
        x[y+4,"Chr"]<-x[1,"Chr"]
        x[y+4,"type"]<-"downstream2"
        x[y+4,"Start"]<-exon_1[1,"End"]+1+1000
        x[y+4,"End"]<-exon_1[1,"End"]+2000
        x[y+4,"Strand"]<-x[1,"Strand"]
        x[y+4,"Score"]<-x[1,"Score"]
        x[y+4,"geneID"]<-x[1,"geneID"]
      } else{}

    } else{}

    return(x)
  } else if(x[1,"Strand"]=="-"){

    exon<-x %>% filter(type=="exon") %>% arrange(Start)
    y1<-dim(exon)[1]
    exon1<-exon[y1, ]
    if((HT.size$V2[HT.size$V1==exon1[1,"Chr"]] - exon1[1,"End"]) >=1000){
      #promoter 1kb
      x[y+1,"Chr"]<-x[1,"Chr"]
      x[y+1,"type"]<-"promoter1"
      x[y+1,"Start"]<-exon1[1,"End"]+1
      x[y+1,"End"]<-exon1[1,"End"]+1000
      x[y+1,"Strand"]<-x[1,"Strand"]
      x[y+1,"Score"]<-x[1,"Score"]
      x[y+1,"geneID"]<-x[1,"geneID"]

      #promoter 2kb
      if((HT.size$V2[HT.size$V1==exon1[1,"Chr"]] - exon1[1,"End"]) >=2000){
        x[y+2,"Chr"]<-x[1,"Chr"]
        x[y+2,"type"]<-"promoter2"
        x[y+2,"Start"]<-exon1[1,"End"]+1000+1
        x[y+2,"End"]<-exon1[1,"End"]+2000
        x[y+2,"Strand"]<-x[1,"Strand"]
        x[y+2,"Score"]<-x[1,"Score"]
        x[y+2,"geneID"]<-x[1,"geneID"]

      } else{}
      #downstream 1kb
      exon_1<-exon[1,]
      if(exon_1[1,"Start"]>=1000){
        x[y+3,"Chr"]<-x[1,"Chr"]
        x[y+3,"type"]<-"downstream1"
        x[y+3,"Start"]<-exon_1[1,"Start"]-1000
        x[y+3,"End"]<-exon_1[1,"Start"]-1
        x[y+3,"Strand"]<-x[1,"Strand"]
        x[y+3,"Score"]<-x[1,"Score"]
        x[y+3,"geneID"]<-x[1,"geneID"]
        #downstream 2kb
        if(exon_1[1,"Start"]>=2000){
          x[y+4,"Chr"]<-x[1,"Chr"]
          x[y+4,"type"]<-"downstream2"
          x[y+4,"Start"]<-exon_1[1,"Start"]-2000
          x[y+4,"End"]<-exon_1[1,"Start"]-1-1000
          x[y+4,"Strand"]<-x[1,"Strand"]
          x[y+4,"Score"]<-x[1,"Score"]
          x[y+4,"geneID"]<-x[1,"geneID"]
        } else{}

      } else {}

    } else {}
    return(x)
  }
})



#calculate the number of exon  and among all genes
ee<-lapply(HT.gene2exon2intron2promoter2downstream.list,function(x){
  exon<-x %>% filter(type=="exon")
  en<-dim(exon)[1]
  intron<-x %>% filter(type=="intron")
  inN<-dim(intron)[1]
  ee<-data.frame(en=en,inN=inN)
  return(ee)
})

ee<-bind_rows(ee,.id = "geneID")


summary(ee$en)

#extract upstream and downstreamgene within each operon



#calculate median length of each feature:exon1,exon2,exon3,etc..
featureLen<-lapply(gene2exon2intron2promoter2downstream,function(x){
  exon<-x %>% filter(type=="exon")
  el<-dim(exon)[1]
  intron <- x %>% filter(type=="intron")
  il<-dim(intron)[1]

  if(exon[1,"Strand"]=="+"){
    exon1<-exon[1,]
    e1_len<-exon1[["End"]]-exon1[["Start"]]
    if(el>1){
      exon2<- exon[2,]
      e2_len<-exon2[["End"]]-exon2[["Start"]]
    } else {e2_len<-NA}

    if(el>2){
      exon3<- exon[3,]
      e3_len<-exon3[["End"]]-exon3[["Start"]]
    } else {e3_len<-NA}

    exon_1<- exon[el,]
    e_1_len<-exon_1[["End"]]-exon_1[["Start"]]
    if(el>1){
      exon_2<-exon[el-1,]
      e_2_len<-exon_2[["End"]]-exon_2[["Start"]]
    } else{e_2_len<-NA}

    if(el>2){
      exon_3<-exon[el-1,]
      e_3_len<-exon_3[["End"]]-exon_3[["Start"]]
    } else{e_3_len<-NA}


    if(dim(intron)[1]!=0){
      #
      intron1<-intron[1,]
      in1_len<-intron1[["End"]]-intron1[["Start"]]
      if(il>1){
        intron2<-intron[2,]
        in2_len<-intron2[["End"]]-intron2[["Start"]]
      } else { in2_len<-NA}

      if(il>2){
        intron3<-intron[3,]
        in3_len<-intron3[["End"]]-intron3[["Start"]]
      } else { in3_len<-NA}
      #
      intron_1<-intron[il,]
      in_1_len<-intron_1[["End"]]-intron_1[["Start"]]
      if(il>1){
        intron_2<-intron[il-1,]
        in_2_len<-intron_2[["End"]]-intron_2[["Start"]]
      } else { in_2_len<-NA}
      if(il>2){
        intron_3<-intron[il-2,]
        in_3_len<-intron_3[["End"]]-intron_3[["Start"]]
      } else { in_3_len<-NA}
    } else {}


    #if strand is "-"
  } else if(exon[1,"Strand"]=="-"){
    exon1<-exon[dim(exon)[1],]
    e1_len<-exon1[["End"]]-exon1[["Start"]]
    if(el>1){
      exon2<- exon[el-1,]
      e2_len<-exon2[["End"]]-exon2[["Start"]]
    } else {e2_len<-NA}

    if(el>2){
      exon3<- exon[el-2,]
      e3_len<-exon3[["End"]]-exon3[["Start"]]
    } else {e3_len<-NA}
    exon_1<-exon[1,]
    e_1_len<-exon_1[["End"]]-exon_1[["Start"]]
    if(el>1){
      exon_2<- exon[2,]
      e_2_len<-exon_2[["End"]]-exon_2[["Start"]]
    } else {e_2_len<-NA}

    if(el>2){
      exon3<- exon[3,]
      e_3_len<-exon_3[["End"]]-exon_3[["Start"]]
    } else {e_3_len<-NA}

    if(il !=0){
      intron1<-intron[il,]
      in1_len<-intron1[["End"]]-intron1[["Start"]]
      if(il>1){
        intron2<-intron[il-1,]
        in2_len<-intron2[["End"]]-intron2[["Start"]]
      } else{in2_len<-NA}
      if(il>2){
        intron3<-intron[il-2,]
        in3_len<-intron3[["End"]]-intron3[["Start"]]
      } else{in3_len<-NA}

      intron_1<-intron[1,]
      in_1_len<-intron_1[["End"]]-intron_1[["Start"]]
      if(il>1){
        intron_2<-intron[2,]
        in_2_len<-intron_2[["End"]]-intron_2[["Start"]]
      } else{in_2_len<-NA}
      if(il>2){
        intron_3<-intron[3,]
        in_3_len<-intron_3[["End"]]-intron_3[["Start"]]
      } else{in_3_len<-NA}
    } else {}


  }

  fl<-rbind(e1_len,e2_len,e3_len,e_1_len,e_2_len,e_3_len,in1_len,in2_len,in3_len,in_1_len,in_2_len,in_3_len)
  fl<-data.frame(type=rownames(fl),len=fl)
  return(fl)
  })

featureLength<-bind_rows(aa,.id = "geneID")
featureLen<-na.omit(featureLength) %>% filter() %>%group_by(type) %>% summarise(mediaLen=median(len))


get_Mbin<-function(x){
  # get each feature of gene
  if(x[1,"Strand"]=="+"){
    exon <- x %>% filter(type=="exon")
    exon1<- exon[1,]
    el<-dim(exon)[1]
    if(el>1){ exon2<- exon[2,]} else {exon2<-NA}
    if(el>2){ exon3<- exon[3,]} else {exon3<-NA}
    exon_1<- exon[el,]
    if(el>1){ exon_2<-exon[el-1,] } else{exon_2<-NA}
    if(el>2){ exon_3<-exon[el-2,] } else{exon_3<-NA}

    intron <- x %>% filter(type=="intron")
    il<-dim(intron)[1]
    if(il!=0){intron1<-intron[1,]} else {intron1<-NA}
    if(il>1){intron2<-intron[2,]} else { intron2<-NA}
    if(il>2){intron3<-intron[3,]} else { intron3<-NA}

    if(il!=0){intron_1<-intron[il,]} else {intron_1<-NA}
    if(il>1){ intron_2<-intron[il-1,] } else {intron_2<-NA}
    if(il>2){ intron_3<-intron[il-2,] } else {intron_3<-NA}

    promoter1<- x %>% filter(type=="promoter1")
    promoter2<- x %>% filter(type=="promoter2")
    downstream1<- x %>% filter(type=="downstream1")
    downstream2<- x %>% filter(type=="downstream2")

  } else if(x[1,"Strand"]=="-"){
    exon <- x %>% filter(type=="exon")
    el<-dim(exon)[1]
    exon1<- exon[el,]
    if(el>1){ exon2<- exon[el-1,]} else {exon2<-NA}
    if(el>2){ exon3<- exon[el-2,]} else {exon3<-NA}
    exon_1<- exon[1,]
    if(el>1){ exon_2<-exon[2,] } else{exon_2<-NA}
    if(el>2){ exon_3<-exon[3,] } else{exon_3<-NA}

    intron <- x %>% filter(type=="intron")
    il<-dim(intron)[1]
    if(il!=0){intron1<-intron[il,]} else {intron1<-NA}
    if(il>1){intron2<-intron[il-1,]} else { intron2<-NA}
    if(il>2){intron3<-intron[il-2,]} else { intron3<-NA}

    if(il!=0){intron_1<-intron[1,]} else {intron_1<-NA}
    if(il>1){ intron_2<-intron[2,] } else {intron_2<-NA}
    if(il>2){ intron_3<-intron[3,] } else {intron_3<-NA}

    promoter1<- x %>% filter(type=="promoter1")
    promoter2<- x %>% filter(type=="promoter2")
    downstream1<- x %>% filter(type=="downstream1")
    downstream2<- x %>% filter(type=="downstream2")
  }

  feature<-rbind(exon1,exon2,exon3,exon_1,exon_2,exon_3,intron1,intron2,intron3,intron_1,
                 intron_2,intron_3,promoter1,promoter2,downstream1,downstream2)
  feature$type<-c("exon1","exon2","exon3","exon_1","exon_2","exon_3","intron1","intron2",
                  "intron3","intron_1","intron_2","intron_3","promoter1","promoter2","downstream1","downstream2")

  smin<-min(na.omit(feature$Start))
  smax<-max(na.omit(feature$End))

  site<-character()
  for (i in 1:(smax-smin)){
    site[i]<-paste0(feature[["Chr"]],":",smin+i-1)
  }

  site_me<-wgbs.data[wgbs.data$Site %in% site,]
  if(dim(site_me)[1] !=0){
    as.numeric(site_me$Start)->site_me$Start
    feature<-na.omit(feature)
    site_me<-site_me %>% group_by(Site) %>% summarise(Site=Site,Start=Start,Methyl= Methyl, UnMethyl=UnMethyl,group=c("a","b","c","d","e","f"),.groups = "drop")
    site_me[,'type']<-NA
    datalist=list()
    for (i in 1:dim(feature)[1]){ #since some feature might have exon less than 3, then the exon1 is also a exon_1. therefore site_me should be counted by twice.

      aa<-feature[i,]
      site_me$type[site_me$Start>= aa$Start & site_me$Start <= aa$End]<-aa$type
      datalist[[i]]<-site_me
    }
    site_me<-do.call(rbind,datalist)
    site_me<-unique(na.omit(site_me))
    site_me<-site_me %>% separate(Site,into = c("Chr","Start"),sep = ":",remove = F)
    as.integer(site_me$Start)->site_me$Start

    for (i in featureLen$type){
      if(i %in% feature$type){
        bin_nu<-ceiling(featureLen$mediaLen[featureLen$type==i]/50)
        f_len<-feature$End[feature$type==i]-feature$Start[feature$type==i]+1
        bin_len<-ceiling(f_len/bin_nu)

        for(bin in 1:bin_nu){
          if(feature[1,"Strand"]=="+"){
            s<-feature$Start[feature$type==i]+bin_len*(bin-1)
            e<-feature$Start[feature$type==i]+bin_len*bin-1
            if(s<=feature$End[feature$type==i] ){
              s<-s
            } else{}

            if(feature$End[feature$type==i] >=e){
              e<-e
            } else {e<-feature$End[feature$type==i]}

            site_me[site_me$Start %in% (s:e),"bin"]<-bin # some bins might have no methylation detected!

          } else if(feature[1,"Strand"]=="-"){
            s<-feature$Start[feature$type==i]+bin_len*(bin-1)
            e<-feature$Start[feature$type==i]+bin_len*bin-1

            if(feature$End[feature$type==i] >=e){
              e<-e
            } else {e<-feature$End[feature$type==i]}

            site_me[site_me$Start %in% (s:e),"bin"]<-bin_nu-bin+1 # some bins might have no methylation detected!
          }

        }

      } else{}

    }
    #summarise the result
    site_me<-site_me %>% group_by(group,type,bin) %>% summarise(Fre=length(which(Methyl>0))/length(bin),.groups = "drop")
    site_me<-site_me %>% group_by(type,bin) %>% summarise(Fre=median(Fre),.groups = "drop") # .groups="drop" : Don't show messages on console.

    return(site_me)

  } else {return(NA)}

}

get_bin<-function(x){  # if unexpected error happen, this function will skip the error!
  return( tryCatch(get_Mbin(x),error=function(e) NULL))
}

ciMethyl.list<-lapply(HT.gene2exon2intron2promoter2downstream.list,get_bin)

ciMethyl.list->aa
aa[!is.na(aa)]->cc
bb<-bind_rows(cc,.id = "geneID")

for (i in somethingWrong$type){
 bb<- bb %>% filter( !(type==i & bin> somethingWrong$len[somethingWrong$type==i]) )
}




bb->ciMethyl.data

ciMethyl.data<-na.omit(ciMethyl.data)


aa<-ciMethyl.data %>% group_by(type,bin) %>% summarise(Fre=median(Fre))
aa<-ciMethyl.data %>% group_by(type,bin) %>% summarise(Fre=mean(Fre))


aa1 <-aa %>% filter(type =="promoter2")
aa2 <-aa %>% filter(type =="promoter1")
aa3 <-aa %>% filter(type =="exon1")
aa4 <-aa %>% filter(type =="intron1")
aa5 <-aa %>% filter(type =="exon2")
aa6 <-aa %>% filter(type =="intron2")
aa7 <-aa %>% filter(type =="exon3")
aa8<-data.frame(type=" ",bin=c(1:10),Fre=0)
aa9 <-aa %>% filter(type =="exon_3")
aa10 <-aa %>% filter(type =="intron_2")
aa11 <-aa %>% filter(type =="exon_2")
aa12 <-aa %>% filter(type =="intron_1")
aa13 <-aa %>% filter(type =="exon_1")
aa14 <-aa %>% filter(type =="downstream1")
aa15 <-aa %>% filter(type =="downstream2")

aa<-rbind(aa1,aa2,aa3,aa4,aa5,aa6,aa7,aa8,aa9,aa10,aa11,aa12,aa13,aa14,aa15)
aa$seq<-seq(1,dim(aa)[1])
aa$seq<-factor(aa$seq,levels = c(1:dim(aa)[1]))

ggplot(aa,aes(x=seq,y=Fre,fill=type))+geom_col(width = 1)+theme_classic()+ylab("Methylated positions per unit \n genomic length (arb. units) ")+
  geom_line(aes(x=seq,y=Fre*0.8),group=1,color="#d1495b" )+xlab("")+
  theme(axis.ticks.length = unit(0.25,"cm"),axis.line.x = element_blank(),axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(colour = "black",size = 1),
        axis.text.x = element_blank(),legend.position = "none",
        axis.line.y = element_line(color="black",lineend = "butt" ,linetype = "solid",size = 0.8),
        axis.text.y = element_text(margin = margin(-10,unit = "cm"),size = 13),
        axis.title.y = element_text(size = 20))+
  scale_fill_manual(values = rep("grey60",length(aa$type)))+
  #promoter
  geom_segment(aes(x=0.5,y=-0.02,xend=20.5,yend=-0.02),size=4,colour="#E08B00")+ #promoter 2kb
  geom_segment(aes(x=20.5,y=-0.02,xend=40.5,yend=-0.02),size=4,colour="#E08B00")+ #promoter 1kb
  #
  geom_segment(aes(x=40.5,y=-0.02,xend=45.5,yend=-0.02),size=8,colour="#00A5FF")+ #exon1
  geom_segment(aes(x=45.5,y=-0.02,xend=54.5,yend=-0.02),size=4,colour="#00B8E5")+ #intron1
  geom_segment(aes(x=54.5,y=-0.02,xend=57.5,yend=-0.02),size=8,colour="#00A5FF")+ #exon2
  geom_segment(aes(x=57.5,y=-0.02,xend=66.5,yend=-0.02),size=4,colour="#00B8E5")+ #intron2
  geom_segment(aes(x=66.5,y=-0.02,xend=69.5,yend=-0.02),size=8,colour="#00A5FF")+ #exon3
  geom_segment(aes(x=69.5,y=-0.02,xend=79.5,yend=-0.02),size=4,colour="#E08B00")+ #nothing


  geom_segment(aes(x=79.5,y=-0.02,xend=82.5,yend=-0.02),size=8,colour="#00A5FF")+ #exon_3
  geom_segment(aes(x=82.5,y=-0.02,xend=91.5,yend=-0.02),size=4,colour="#00B8E5")+ #intron_2
  geom_segment(aes(x=91.5,y=-0.02,xend=94.5,yend=-0.02),size=8,colour="#00A5FF")+ #exon_2
  geom_segment(aes(x=94.5,y=-0.02,xend=102.5,yend=-0.02),size=4,colour="#00B8E5")+ #intron_1
  geom_segment(aes(x=102.5,y=-0.02,xend=113.5,yend=-0.02),size=8,colour="#00A5FF")+ #exon_1
  geom_segment(aes(x=113.5,y=-0.02,xend=133.5,yend=-0.02),size=4,colour="#E08B00")+ #downstream 1kb
  geom_segment(aes(x=133.5,y=-0.02,xend=153.5,yend=-0.02),size=4,colour="#E08B00")+ # downstream 2kb

  annotate("text",x=10,,y=-0.03,label="2kb \npromoter ",angle=45,size=4)+ #promoter 2kb
  annotate("text",x=30,,y=-0.03,label="1kb \npromoter ",angle=45,size=4)+#promoter 1kb
  annotate("text",x=42.5,,y=-0.03,label="exon 1",angle=45,size=4)+#exon1
  annotate("text",x=50.5,,y=-0.03,label="intron 1 ",angle=45,size=4)+#intron1
  annotate("text",x=56.5,,y=-0.03,label="exon 2 ",angle=45,size=4)+ #exon2
  annotate("text",x=62,,y=-0.03,label="intron 2 ",angle=45,size=4)+ #intron2
  annotate("text",x=68,,y=-0.03,label="exon 3 ",angle=45,size=4)+ #exon3
  annotate("text",x=81.5,,y=-0.03,label="exon -3 ",angle=45,size=4)+ #exon -3
  annotate("text",x=87.5,y=-0.03,label="intron -2 ",angle=45,size=4)+  #intron -2
  annotate("text",x=94,y=-0.03,label="exon -2 ",angle=45,size=4)+ #exon -2
  annotate("text",x=99,y=-0.03,label="intron -1 ",angle=45,size=4)+ #intron -1
  annotate("text",x=108,y=-0.03,label="exon -1 ",angle=45,size=4)+
  annotate("text",x=128,,y=-0.03,label="1kb \n downstream",angle=45,size=4)+
  annotate("text",x=142,,y=-0.03,label="2kb \n downstream",angle=45,size=4)+
  geom_segment(aes(x=40.5,xend=40.5,y=-0.0175,yend=0.12),size=0.7,colour="black") + #exon1

  geom_segment(aes(x=45.5,xend=45.5,y=-0.0175,yend=0.12),size=0.7,colour="black",linetype="dashed")+
  geom_segment(aes(x=54.5,xend=54.5,y=-0.0175,yend=0.12),size=0.7,colour="black",linetype="dashed")+
  geom_segment(aes(x=57.5,xend=57.5,y=-0.0175,yend=0.12),size=0.7,colour="black",linetype="dashed")+
  geom_segment(aes(x=66.5,xend=66.5,y=-0.0175,yend=0.12),size=0.7,colour="black",linetype="dashed")+
  geom_segment(aes(x=69.5,xend=69.5,y=-0.0175,yend=0.12),size=0.7,colour="black",linetype="dashed")+
  geom_segment(aes(x=79.5,xend=79.5,y=-0.0175,yend=0.12),size=0.7,colour="black",linetype="dashed")+
  geom_segment(aes(x=82.5,xend=82.5,y=-0.0175,yend=0.12),size=0.7,colour="black",linetype="dashed")+
  geom_segment(aes(x=91.5,xend=91.5,y=-0.0175,yend=0.12),size=0.7,colour="black",linetype="dashed")+
  geom_segment(aes(x=94.5,xend=94.5,y=-0.0175,yend=0.12),size=0.7,colour="black",linetype="dashed")+
  geom_segment(aes(x=102.5,xend=102.5,y=-0.0175,yend=0.12),size=0.7,colour="black",linetype="dashed")+
  geom_segment(aes(x=113.5,xend=113.5,y=-0.0175,yend=0.12),size=0.7,colour="black")+
  geom_segment(aes(x=-1,xend=-1,y=0,yend=0.12),size=0.7,colour="black",linetype="dashed")




#nonoperon
ciMethyl.nonoperon.data<-ciMethyl.data %>% filter(! geneID %in% ci.HT.operon$Geneid)# select gene that not in operon

aa<-ciMethyl.nonoperon.data %>% group_by(type,bin) %>% summarise(Fre=median(Fre))
aa<-ciMethyl.nonoperon.data %>% group_by(type,bin) %>% summarise(Fre=mean(Fre))


aa1 <-aa %>% filter(type =="promoter2")
aa2 <-aa %>% filter(type =="promoter1")
aa3 <-aa %>% filter(type =="exon1")
aa4 <-aa %>% filter(type =="intron1")
aa5 <-aa %>% filter(type =="exon2")
aa6 <-aa %>% filter(type =="intron2")
aa7 <-aa %>% filter(type =="exon3")
aa8<-data.frame(type=" ",bin=c(1:10),Fre=0)
aa9 <-aa %>% filter(type =="exon_3")
aa10 <-aa %>% filter(type =="intron_2")
aa11 <-aa %>% filter(type =="exon_2")
aa12 <-aa %>% filter(type =="intron_1")
aa13 <-aa %>% filter(type =="exon_1")
aa14 <-aa %>% filter(type =="downstream1")
aa15 <-aa %>% filter(type =="downstream2")

aa<-rbind(aa1,aa2,aa3,aa4,aa5,aa6,aa7,aa8,aa9,aa10,aa11,aa12,aa13,aa14,aa15)
aa$seq<-seq(1,dim(aa)[1])
aa$seq<-factor(aa$seq,levels = c(1:dim(aa)[1]))

ggplot(aa,aes(x=seq,y=Fre,fill=type))+geom_col(width = 1)+theme_classic()+ylab("Methylated positions per unit \n genomic length (arb. units) ")+
  geom_line(aes(x=seq,y=Fre*0.8),group=1,color="#d1495b" )+xlab("")+
  theme(axis.ticks.length = unit(0.25,"cm"),axis.line.x = element_blank(),axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(colour = "black",size = 1),
        axis.text.x = element_blank(),legend.position = "none",
        axis.line.y = element_line(color="black",lineend = "butt" ,linetype = "solid",size = 0.8),
        axis.text.y = element_text(margin = margin(-10,unit = "cm"),size = 13),
        axis.title.y = element_text(size = 20))+
  scale_fill_manual(values = rep("grey60",length(aa$type)))+
  #promoter
  geom_segment(aes(x=0.5,y=-0.02,xend=20.5,yend=-0.02),size=4,colour="#E08B00")+ #promoter 2kb
  geom_segment(aes(x=20.5,y=-0.02,xend=40.5,yend=-0.02),size=4,colour="#E08B00")+ #promoter 1kb
  #
  geom_segment(aes(x=40.5,y=-0.02,xend=45.5,yend=-0.02),size=8,colour="#00A5FF")+ #exon1
  geom_segment(aes(x=45.5,y=-0.02,xend=54.5,yend=-0.02),size=4,colour="#00B8E5")+ #intron1
  geom_segment(aes(x=54.5,y=-0.02,xend=57.5,yend=-0.02),size=8,colour="#00A5FF")+ #exon2
  geom_segment(aes(x=57.5,y=-0.02,xend=66.5,yend=-0.02),size=4,colour="#00B8E5")+ #intron2
  geom_segment(aes(x=66.5,y=-0.02,xend=69.5,yend=-0.02),size=8,colour="#00A5FF")+ #exon3
  geom_segment(aes(x=69.5,y=-0.02,xend=79.5,yend=-0.02),size=4,colour="#E08B00")+ #nothing


  geom_segment(aes(x=79.5,y=-0.02,xend=82.5,yend=-0.02),size=8,colour="#00A5FF")+ #exon_3
  geom_segment(aes(x=82.5,y=-0.02,xend=91.5,yend=-0.02),size=4,colour="#00B8E5")+ #intron_2
  geom_segment(aes(x=91.5,y=-0.02,xend=94.5,yend=-0.02),size=8,colour="#00A5FF")+ #exon_2
  geom_segment(aes(x=94.5,y=-0.02,xend=102.5,yend=-0.02),size=4,colour="#00B8E5")+ #intron_1
  geom_segment(aes(x=102.5,y=-0.02,xend=113.5,yend=-0.02),size=8,colour="#00A5FF")+ #exon_1
  geom_segment(aes(x=113.5,y=-0.02,xend=133.5,yend=-0.02),size=4,colour="#E08B00")+ #downstream 1kb
  geom_segment(aes(x=133.5,y=-0.02,xend=153.5,yend=-0.02),size=4,colour="#E08B00")+ # downstream 2kb

  annotate("text",x=10,,y=-0.03,label="2kb \npromoter ",angle=45,size=4)+ #promoter 2kb
  annotate("text",x=30,,y=-0.03,label="1kb \npromoter ",angle=45,size=4)+#promoter 1kb
  annotate("text",x=42.5,,y=-0.03,label="exon 1",angle=45,size=4)+#exon1
  annotate("text",x=50.5,,y=-0.03,label="intron 1 ",angle=45,size=4)+#intron1
  annotate("text",x=56.5,,y=-0.03,label="exon 2 ",angle=45,size=4)+ #exon2
  annotate("text",x=62,,y=-0.03,label="intron 2 ",angle=45,size=4)+ #intron2
  annotate("text",x=68,,y=-0.03,label="exon 3 ",angle=45,size=4)+ #exon3
  annotate("text",x=81.5,,y=-0.03,label="exon -3 ",angle=45,size=4)+ #exon -3
  annotate("text",x=87.5,y=-0.03,label="intron -2 ",angle=45,size=4)+  #intron -2
  annotate("text",x=94,y=-0.03,label="exon -2 ",angle=45,size=4)+ #exon -2
  annotate("text",x=99,y=-0.03,label="intron -1 ",angle=45,size=4)+ #intron -1
  annotate("text",x=108,y=-0.03,label="exon -1 ",angle=45,size=4)+
  annotate("text",x=128,,y=-0.03,label="1kb \n downstream",angle=45,size=4)+
  annotate("text",x=142,,y=-0.03,label="2kb \n downstream",angle=45,size=4)+
  geom_segment(aes(x=40.5,xend=40.5,y=-0.0175,yend=0.12),size=0.7,colour="black") + #exon1

  geom_segment(aes(x=45.5,xend=45.5,y=-0.0175,yend=0.12),size=0.7,colour="black",linetype="dashed")+
  geom_segment(aes(x=54.5,xend=54.5,y=-0.0175,yend=0.12),size=0.7,colour="black",linetype="dashed")+
  geom_segment(aes(x=57.5,xend=57.5,y=-0.0175,yend=0.12),size=0.7,colour="black",linetype="dashed")+
  geom_segment(aes(x=66.5,xend=66.5,y=-0.0175,yend=0.12),size=0.7,colour="black",linetype="dashed")+
  geom_segment(aes(x=69.5,xend=69.5,y=-0.0175,yend=0.12),size=0.7,colour="black",linetype="dashed")+
  geom_segment(aes(x=79.5,xend=79.5,y=-0.0175,yend=0.12),size=0.7,colour="black",linetype="dashed")+
  geom_segment(aes(x=82.5,xend=82.5,y=-0.0175,yend=0.12),size=0.7,colour="black",linetype="dashed")+
  geom_segment(aes(x=91.5,xend=91.5,y=-0.0175,yend=0.12),size=0.7,colour="black",linetype="dashed")+
  geom_segment(aes(x=94.5,xend=94.5,y=-0.0175,yend=0.12),size=0.7,colour="black",linetype="dashed")+
  geom_segment(aes(x=102.5,xend=102.5,y=-0.0175,yend=0.12),size=0.7,colour="black",linetype="dashed")+
  geom_segment(aes(x=113.5,xend=113.5,y=-0.0175,yend=0.12),size=0.7,colour="black")+
  geom_segment(aes(x=-1,xend=-1,y=0,yend=0.12),size=0.7,colour="black",linetype="dashed")


#operon
ciMethyl.operon.data<-ciMethyl.data %>% filter(geneID %in% ci.HT.operon$Geneid)# select gene that in operon

aa<-ciMethyl.operon.data %>% group_by(type,bin) %>% summarise(Fre=median(Fre))
aa<-ciMethyl.operon.data %>% group_by(type,bin) %>% summarise(Fre=mean(Fre))


aa1 <-aa %>% filter(type =="promoter2")
aa2 <-aa %>% filter(type =="promoter1")
aa3 <-aa %>% filter(type =="exon1")
aa4 <-aa %>% filter(type =="intron1")
aa5 <-aa %>% filter(type =="exon2")
aa6 <-aa %>% filter(type =="intron2")
aa7 <-aa %>% filter(type =="exon3")
aa8<-data.frame(type=" ",bin=c(1:10),Fre=0)
aa9 <-aa %>% filter(type =="exon_3")
aa10 <-aa %>% filter(type =="intron_2")
aa11 <-aa %>% filter(type =="exon_2")
aa12 <-aa %>% filter(type =="intron_1")
aa13 <-aa %>% filter(type =="exon_1")
aa14 <-aa %>% filter(type =="downstream1")
aa15 <-aa %>% filter(type =="downstream2")

aa<-rbind(aa1,aa2,aa3,aa4,aa5,aa6,aa7,aa8,aa9,aa10,aa11,aa12,aa13,aa14,aa15)
aa$seq<-seq(1,dim(aa)[1])
aa$seq<-factor(aa$seq,levels = c(1:dim(aa)[1]))

ggplot(aa,aes(x=seq,y=Fre,fill=type))+geom_col(width = 1)+theme_classic()+ylab("Methylated positions per unit \n genomic length (arb. units) ")+
  geom_line(aes(x=seq,y=Fre*0.8),group=1,color="#d1495b" )+xlab("")+
  theme(axis.ticks.length = unit(0.25,"cm"),axis.line.x = element_blank(),axis.ticks.x = element_blank(),
        axis.ticks.y = element_line(colour = "black",size = 1),
        axis.text.x = element_blank(),legend.position = "none",
        axis.line.y = element_line(color="black",lineend = "butt" ,linetype = "solid",size = 0.8),
        axis.text.y = element_text(margin = margin(-10,unit = "cm"),size = 13),
        axis.title.y = element_text(size = 20))+
  scale_fill_manual(values = rep("grey60",length(aa$type)))+
  #promoter
  geom_segment(aes(x=0.5,y=-0.02,xend=20.5,yend=-0.02),size=4,colour="#E08B00")+ #promoter 2kb
  geom_segment(aes(x=20.5,y=-0.02,xend=40.5,yend=-0.02),size=4,colour="#E08B00")+ #promoter 1kb
  #
  geom_segment(aes(x=40.5,y=-0.02,xend=45.5,yend=-0.02),size=8,colour="#00A5FF")+ #exon1
  geom_segment(aes(x=45.5,y=-0.02,xend=54.5,yend=-0.02),size=4,colour="#00B8E5")+ #intron1
  geom_segment(aes(x=54.5,y=-0.02,xend=57.5,yend=-0.02),size=8,colour="#00A5FF")+ #exon2
  geom_segment(aes(x=57.5,y=-0.02,xend=66.5,yend=-0.02),size=4,colour="#00B8E5")+ #intron2
  geom_segment(aes(x=66.5,y=-0.02,xend=69.5,yend=-0.02),size=8,colour="#00A5FF")+ #exon3
  geom_segment(aes(x=69.5,y=-0.02,xend=79.5,yend=-0.02),size=4,colour="#E08B00")+ #nothing


  geom_segment(aes(x=79.5,y=-0.02,xend=82.5,yend=-0.02),size=8,colour="#00A5FF")+ #exon_3
  geom_segment(aes(x=82.5,y=-0.02,xend=91.5,yend=-0.02),size=4,colour="#00B8E5")+ #intron_2
  geom_segment(aes(x=91.5,y=-0.02,xend=94.5,yend=-0.02),size=8,colour="#00A5FF")+ #exon_2
  geom_segment(aes(x=94.5,y=-0.02,xend=102.5,yend=-0.02),size=4,colour="#00B8E5")+ #intron_1
  geom_segment(aes(x=102.5,y=-0.02,xend=113.5,yend=-0.02),size=8,colour="#00A5FF")+ #exon_1
  geom_segment(aes(x=113.5,y=-0.02,xend=133.5,yend=-0.02),size=4,colour="#E08B00")+ #downstream 1kb
  geom_segment(aes(x=133.5,y=-0.02,xend=153.5,yend=-0.02),size=4,colour="#E08B00")+ # downstream 2kb

  annotate("text",x=10,,y=-0.03,label="2kb \npromoter ",angle=45,size=4)+ #promoter 2kb
  annotate("text",x=30,,y=-0.03,label="1kb \npromoter ",angle=45,size=4)+#promoter 1kb
  annotate("text",x=42.5,,y=-0.03,label="exon 1",angle=45,size=4)+#exon1
  annotate("text",x=50.5,,y=-0.03,label="intron 1 ",angle=45,size=4)+#intron1
  annotate("text",x=56.5,,y=-0.03,label="exon 2 ",angle=45,size=4)+ #exon2
  annotate("text",x=62,,y=-0.03,label="intron 2 ",angle=45,size=4)+ #intron2
  annotate("text",x=68,,y=-0.03,label="exon 3 ",angle=45,size=4)+ #exon3
  annotate("text",x=81.5,,y=-0.03,label="exon -3 ",angle=45,size=4)+ #exon -3
  annotate("text",x=87.5,y=-0.03,label="intron -2 ",angle=45,size=4)+  #intron -2
  annotate("text",x=94,y=-0.03,label="exon -2 ",angle=45,size=4)+ #exon -2
  annotate("text",x=99,y=-0.03,label="intron -1 ",angle=45,size=4)+ #intron -1
  annotate("text",x=108,y=-0.03,label="exon -1 ",angle=45,size=4)+
  annotate("text",x=128,,y=-0.03,label="1kb \n downstream",angle=45,size=4)+
  annotate("text",x=142,,y=-0.03,label="2kb \n downstream",angle=45,size=4)+
  geom_segment(aes(x=40.5,xend=40.5,y=-0.0175,yend=0.12),size=0.7,colour="black") + #exon1

  geom_segment(aes(x=45.5,xend=45.5,y=-0.0175,yend=0.12),size=0.7,colour="black",linetype="dashed")+
  geom_segment(aes(x=54.5,xend=54.5,y=-0.0175,yend=0.12),size=0.7,colour="black",linetype="dashed")+
  geom_segment(aes(x=57.5,xend=57.5,y=-0.0175,yend=0.12),size=0.7,colour="black",linetype="dashed")+
  geom_segment(aes(x=66.5,xend=66.5,y=-0.0175,yend=0.12),size=0.7,colour="black",linetype="dashed")+
  geom_segment(aes(x=69.5,xend=69.5,y=-0.0175,yend=0.12),size=0.7,colour="black",linetype="dashed")+
  geom_segment(aes(x=79.5,xend=79.5,y=-0.0175,yend=0.12),size=0.7,colour="black",linetype="dashed")+
  geom_segment(aes(x=82.5,xend=82.5,y=-0.0175,yend=0.12),size=0.7,colour="black",linetype="dashed")+
  geom_segment(aes(x=91.5,xend=91.5,y=-0.0175,yend=0.12),size=0.7,colour="black",linetype="dashed")+
  geom_segment(aes(x=94.5,xend=94.5,y=-0.0175,yend=0.12),size=0.7,colour="black",linetype="dashed")+
  geom_segment(aes(x=102.5,xend=102.5,y=-0.0175,yend=0.12),size=0.7,colour="black",linetype="dashed")+
  geom_segment(aes(x=113.5,xend=113.5,y=-0.0175,yend=0.12),size=0.7,colour="black")+
  geom_segment(aes(x=-1,xend=-1,y=0,yend=0.12),size=0.7,colour="black",linetype="dashed")

#wow
