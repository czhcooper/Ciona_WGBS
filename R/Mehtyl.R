

# read samples files
setwd("/Volumes/data_cooper/ciona_i/WGBS")
wgbs.files<-list.files("/Volumes/data_cooper/ciona_i/WGBS")
aa<-read.table(wgbs.files[1],sep ="\t",quote = "",header = F )
bb<-read.table(wgbs.files[2],sep ="\t",quote = "",header = F )
cc<-read.table(wgbs.files[3],sep ="\t",quote = "",header = F )
dd<-read.table(wgbs.files[4],sep ="\t",quote = "",header = F )
ee<-read.table(wgbs.files[5],sep ="\t",quote = "",header = F )
ff<-read.table(wgbs.files[6],sep ="\t",quote = "",header = F )

aa$sample<-"aa"
bb$sample<-"bb"
cc$sample<-"cc"
dd$sample<-"dd"
ee$sample<-"ee"
ff$sample<-"ff"

d<-1
rbind(aa,bb,cc,dd,ee,ff)->gg




data <- do.call("rbind", lapply(wgbs.files, function(fn)
  data.frame(Filename=fn, read.table(fn,sep ="\t",quote = "",header = F)
  )))




colnames(aa)<-c("Chr","Start","End","Percentage","Methy","UnMethy")


#annotate downstream region
gene2exon2intron2promoter2downstream<-lapply(gene2exon2intron2promoter1,function(x){

  fl<-dim(x)[1]
  if(x[1,"Strand"]=="+"){

    exon<-x %>% filter(type=="exon") %>% arrange(Start)
    last_exon<- exon[dim(exon)[1],]
    #downstream 1kb
    if(HT.size$V2[HT.size$V1==x[1,"Chr"]]-last_exon[1,"End"]>=1000){
      x[fl+1,"Chr"]<-x[1,"Chr"]
      x[fl+1,"type"]<-"downstream1"
      x[fl+1,"Start"]<-last_exon[1,"End"]+1
      x[fl+1,"End"]<-x[fl+1,"Start"]+1000
      x[fl+1,"geneID"]<-x[1,"geneID"]
      x[fl+1,"Strand"]<-"+"
      x[fl+1,"Score"]<-"."
    } else {}

    #donwstream 2kb

    if(HT.size$V2[HT.size$V1==x[1,"Chr"]]-last_exon[1,"End"]>=2000){
      x[fl+2,"Chr"]<-x[1,"Chr"]
      x[fl+2,"type"]<-"downstream2"
      x[fl+2,"Start"]<-x[fl+1,"End"]+1
      x[fl+2,"End"]<-x[fl+2,"Start"]+1000-1
      x[fl+2,"geneID"]<-x[1,"geneID"]
      x[fl+2,"Strand"]<-"+"
      x[fl+2,"Score"]<-"."
    } else { }

  } else {
    exon<-x %>% filter(type=="exon") %>% arrange(Start)
    # reverse strand
    colnames(exon)[c(3,4)]<-c("End","Start")
    last_exon<- exon[1,]
    colnames(last_exon)[c(3,4)]<-c("Start","End")
    # downstream 1kb
    if(last_exon[1,"Start"]>=1000){
      x[fl+1,"Chr"]<-x[1,"Chr"]
      x[fl+1,"type"]<-"downstream1"
      x[fl+1,"Start"]<-last_exon[1,"Start"]-1000
      x[fl+1,"End"]<-x[fl+1,"Start"]+1000-1
      x[fl+1,"geneID"]<-x[1,"geneID"]
      x[fl+1,"Strand"]<-"-"
      x[fl+1,"Score"]<-"."
    } else {}

    #donwstream 2kb

    if(last_exon[1,"Start"]>=2000){
      x[fl+2,"Chr"]<-x[1,"Chr"]
      x[fl+2,"type"]<-"downstream2"
      x[fl+2,"Start"]<-last_exon[1,"Start"]-2000
      x[fl+2,"End"]<-x[fl+2,"Start"]+1000-1
      x[fl+2,"geneID"]<-x[1,"geneID"]
      x[fl+2,"Strand"]<-"-"
      x[fl+2,"Score"]<-"."
    } else { }
    }
   x<- na.omit(x)
   return(x)
})




#

cal_Mbin<-function(x){
  if(x[1,5]=="+"){
    exon1<-x[1,]
    if(dim(exon1)[1]==1){
      ex_site<-character()
      for ( i in 0:(exon1[["End"]]-exon1[["Start"]]) ) {
        ex_site[i+1]<- paste0(exon1[["Chr"]],":",exon1[["Start"]]+i)
      }

      ex_me<-wgbs.data[wgbs.data$Site %in% ex_site,]

      if(dim(ex_me)[1]!=0){
        exon_len<-exon1[["End"]]-exon1[["Start"]]+1
        bin_len<-round(exon_len/bin_nu)
        ex_me %>% separate(Site,into = c("Chr","Start"),sep = ":",remove = F)->ex_me

        for (bin in 1:bin_nu ){
          s<-bin_len*(bin-1)+exon1[["Start"]]
          e<-bin_len*bin+exon1[["Start"]]-1

          ex_me[ex_me$Start %in% (s:e),"bin"]<-bin # some bins might have no methylation detected!

        }
        na.omit(ex_me)->ex_me
        ex_me %>% group_by( Site) %>% summarise(Methyl= Methyl, UnMethyl=UnMethyl,group=c("a","b","c","d","e","f"),bin=bin,.groups = "drop")->ex_me
        ex_bin_me <-ex_me %>% group_by( group,bin) %>% summarise(Fre=length(Methyl>0)/bin_len,.groups = "drop")
        ex_bin_me <-ex_bin_me  %>% group_by(bin) %>% summarise(Fre=median(Fre))

        return(ex_bin_me)
      } else {
        return(NA)
      }

    } else { return(NA) }
  } else if(x[1,5]=="-"){

    exon<-x[x$type=="exon",]
    exon_last<-exon[dim(exon)[1],]
    exon1<-exon_last
    if(dim(exon1)[1]==1){
      ex_site<-character()
      for ( i in 0:(exon1[["End"]]-exon1[["Start"]]) ) {
        ex_site[i+1]<- paste0(exon1[["Chr"]],":",exon1[["Start"]]+i)
      }

      ex_me<-wgbs.data[wgbs.data$Site %in% ex_site,]
      if(dim(ex_me)[1]!=0){
        exon1_len<-exon1[["End"]]-exon1[["Start"]]+1
        bin_len<-round(exon_len/bin_nu)
        ex_me %>% separate(Site,into = c("Chr","Start"),sep = ":",remove = F)->ex_me

        for(bin in 1:bin_nu){
          s<-exon1[["End"]]-bin_len*(bin)+1
          e<-exon1[["End"]]-bin_len*(bin-1)
          ex_me[ex_me$Start %in% (s:e),"bin"]<-bin # some bins might have no methylation detected!
        }
        na.omit(ex_me)->ex_me
        ex_me %>% group_by( Site) %>% summarise(Methyl= Methyl, UnMethyl=UnMethyl,group=c("a","b","c","d","e","f"),bin=bin,.groups = "drop")->ex_me
        ex_bin_me <-ex_me %>% group_by( group,bin) %>% summarise(Fre=length(Methyl>0)/bin_len,.groups = "drop")
        ex_bin_me <-ex_bin_me  %>% group_by(bin) %>% summarise(Fre=median(Fre))

        return(ex_bin_me)


      } else{ return(NA)}


    } else{return(NA)}

  }


}


get_Mbin<-function(x){  # if error happen, this function will skip the error!
  return( tryCatch(cal_Mbin(x),error=function(e) NULL))
}

exon1_Methyl[!is.na(exon1_Methyl)]->aa

aa<-bind_rows(aa,.id = "geneID")

aa %>% group_by(bin) %>% summarise(Fre=median(Fre))->bb


