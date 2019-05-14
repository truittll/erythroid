source("~/Documents/code/Utilities.R")
setwd("/Volumes/dirgroup/TSCBB/LAB/Lauren/1. Projects/12. Fanny- nrbc")

jd76=read.combined("JD76_combined_20180928.txt")
zk22=read.combined("ZK22_combined_20180823.txt")
zh19=read.combined("ZH19_combined_20180823.txt")
zh33=read.combined("ZH33_combined_20180806.txt")
rq859=read.combined("RQ859_combined_20190305.txt")
rq3600=read.combined("RQ3600_combined_20190411.txt")

jd76=barcodetrackR::threshold(jd76)
zk22=barcodetrackR::threshold(zk22)
zh19=barcodetrackR::threshold(zh19)
zh33=barcodetrackR::threshold(zh33)
rq85=barcodetrackR::threshold(rq859)
rq3600=barcodetrackR::threshold(rq3600)

jd76=apply(jd76,2,function(x){x/sum(x)})
zk22=apply(zk22,2,function(x){x/sum(x)})
zh19=apply(zh19,2,function(x){x/sum(x)})
zh33=apply(zh33,2,function(x){x/sum(x)})
rq859=apply(rq859,2,function(x){x/sum(x)})
rq3600=apply(rq3600,2,function(x){x/sum(x)})

jd76.lbm.samples=c("JD76_3hm_8n9n2017_LBM_CD34p_sampled_LT2927_FX42_i505_S5_L008_R1_001.fastq",
                   "JD76_3hm_8n9n2017_LBM_T_CELL_sampled_LT2927_FX42_i506_S6_L008_R1_001.fastq",
                   "JD76_3hm_8n9n2017_LBM_B_CELL_sampled_LT2927_FX42_i507_S7_L008_R1_001.fastq",
                   "JD76_3hm_8n9n2017_LBM_MONO_CELL_sampled_LT2927_FX42_i510_S8_L008_R1_001.fastq",
                   "JD76_3hm_8n9n2017_LBM_Grans_sampled_LT2927_FX42_i511_S9_L008_R1_001.fastq",
                   "JD76_3hm_8n9n2017_LBM_NRBC_sampled_LT2927_FX42_i514_S11_L008_R1_001.fastq")

jd76.lbm.data=jd76[,jd76.lbm.samples]

jd76.rbm.samples=c("JD76_3hmm_8_9_2017_RBM_CD34p_sampled_FX39_i511_S9_L003_R1_001.fastq",
                   "JD76_3hmm_8_9_2017_RBM_T_CELL_sampled_FX39_i512_S10_L003_R1_001.fastq",
                   "JD76_3hmm_8_9_2017_RBM_B_CELL_sampled_FX39_i514_S11_L003_R1_001.fastq",
                   "JD76_3hmm_8_9_2017_RBM_MONO_CELL_sampled_FX39_i515_S12_L003_R1_001.fastq",
                   "JD76_3hmm_8_9_2017_RBM_Grans_sampled_FX39_i516_S13_L003_R1_001.fastq",
                   "JD76_3hmm_8_9_2017_RBM_NRBC_sampled_FX39_i519_S15_L003_R1_001.fastq")

jd76.rbm.data=jd76[,jd76.rbm.samples]

zk22.lbm.samples=c("zk22_11hm_LBM_CD34p_20170123_sampled_DE2242_CD143_N7015_S13_L008_R1_001.fastq",
                   "zk22_11hm_LBM_T_20170123_sampled_FX25_N703_S3_L008_R1_001.fastq",
                   "zk22_11hm_LBM_B_20170123_sampled_FX25_N704_S4_L008_R1_001.fastq",
                   "zk22_11hm_LBM_Mono_20170123_sampled_FX25_N705_S5_L008_R1_001.fastq",
                   "zk22_11hm_LBM_Gr_20170123_sampled_DE2242_CD143_N7011_S9_L008_R1_001.fastq",
                   "zk22_11hm_LBM_Gr_20170123_sampled_DE2242_CD143_N7011_S9_L008_R1_001.fastq")

zk22.lbm.data=zk22[,zk22.lbm.samples]

zh19.rbm.samples=c("zh19_45hm_20170302_RBM_CD34p_sampled_FX31_N701_S41_L007_R1_001.fastq",
                   "zh19_45hm_RBM_T_20170302_trimmed_sampled_FX33_R4_S48_L004_R1_001.fastq",
                   "zh19_45hm_RBM_B_20170302_trimmed_sampled_FX33_R5_S49_L004_R1_001.fastq",
                   "zh19_45hm_RBM_B_20170302_trimmed_sampled_FX33_R5_S49_L004_R1_001.fastq",
                   "zh19_45hm_20170302_RBM_Gr_sampled_DE2382_FX30_N705_S15_L004_R1_001.fastq",
                   "zh19_45hm_RBM_nRBC_20170302_trimmed_sampled_FX33_R12_S51_L004_R1_001.fastq")

zh19.rbm.data=zh19[,zh19.rbm.samples]

zh33.lbm.samples=c("zh33_46m_LBM_CD34p_20160325_pp_50trim_CD125_R24.fastq",
                    "sampled_trimmed_ZH33_46m_03_25_16_BM_T_CELL_R1_FX18.fastq",
                   "sampled_trimmed_ZH33_46m_03_25_16__BM_B_CELL_R2_FX18.fastq",
                    "sampled_trimmed_ZH33_46m_03_25_16__BM_MONO_CELL_R3_FX18.fastq",
                   "zh33_46m_LBM_Gr_1_pp_50trim_CD123_R24.fastq",
                   "sampled_trimmed_ZH33_46m_03_25_16__BM_NRBC_R4_FX18.fastq")

zh33.lbm.data=zh33[,zh33.lbm.samples]

zh33.rbm.samples=c("zh33_46m_RBM_CD34p_20160325_pp_50trim_CD125_R28.fastq",
                   "sampled_trimmed_ZH33_46m_R_BM_03_25_16_T_CELL_R1_FX21.fastq",
                   "sampled_trimmed_ZH33_46m_R_BM_03_25_16_B_CELL_R2_FX21.fastq",
                   "sampled_trimmed_ZH33_46m_R_BM_03_25_16_MONO_CELL_R4_FX21.fastq",
                   "zh33_46m_RBM_Gr_1_pp_50trim_CD123_R28.fastq",
                   "sampled_trimmed_ZH33_46m_R_BM_03_25_16_NRBC_CELL_R5_FX21.fastq")

zh33.rbm.data=zh33[,zh33.rbm.samples]

rq859.samples=c("RQ859 40m 20181107 CD34p CD41n BM  CD164",
                "RQ859 40m 20181107 T BM  CD164",
                "RQ859 40m 20181107 B BM  CD164",
                "RQ859 40m 20181107 Grans BM  CD164",
                "RQ859 40m 20181107 nRBC BM  CD164")
rq859.samples=filename(rq859.samples,read.delim("RQ859_key_20190305.txt",header=TRUE,row.names=NULL))

rq859.data=rq859[,rq859.samples]

rq3600.samples=c("RQ3600_48hm_20190312_BM_CD34p_CD41n____sampled_FX53_i510_S32_L006_R1_001.fastq",
                 "RQ3600_48hm_20190312_BM_T____sampled_FX53_i502_S26_L006_R1_001.fastq",
                 "RQ3600_48hm_20190312_BM_B____sampled_FX53_i503_S27_L006_R1_001.fastq",
                 "RQ3600_48hm_20190312_BM_Grans____sampled_FX53_i501_S25_L006_R1_001.fastq",
                 "RQ3600_48hm_20190314_BM_MKP____sampled_FX53_i507_S31_L006_R1_001.fastq")

rq3600.data=rq3600[,rq3600.samples]

assign.bias=function(data,index=ncol(data)){
  other=data[,-index]
  selected=data[,index]
  
  max.other=apply(other,1,max)
  
  bias=selected/max.other
  
  bias=bias[!is.na(bias)]
  bias=as.data.frame(bias)
  
  barcodes.g10x=rownames(bias)[bias>=10]
  barcodes.5xg10x=rownames(bias)[bias<10&bias>=5]
  barcodes.2xg5x=rownames(bias)[bias<5&bias>=2]
  barcodes.1xg2x=rownames(bias)[bias<2&bias>=1]
  barcodes.1xg2xA=rownames(bias)[bias<1&bias>=1/2]
  barcodes.2xg5xA=rownames(bias)[bias<1/2&bias>=1/5]
  barcodes.5xg10xA=rownames(bias)[bias<1/5&bias>=1/10]
  barcodes.g10xA=rownames(bias)[bias<1/10]
  
  return(list(barcodes.g10x=barcodes.g10x,barcodes.5xg10x=barcodes.5xg10x,barcodes.2xg5x=barcodes.2xg5x,
              barcodes.1xg2x=barcodes.1xg2x,barcodes.1xg2xA=barcodes.1xg2xA,barcodes.2xg5xA=barcodes.2xg5xA,
              barcodes.5xg10xA=barcodes.5xg10xA,barcodes.g10xA=barcodes.g10xA))
}

bias.plot=function(data,index=ncol(data)){
  bias=assign.bias(data,index)
  
  barcodes.g10x=bias$barcodes.g10x
  barcodes.5xg10x=bias$barcodes.5xg10x
  barcodes.2xg5x=bias$barcodes.2xg5x
  barcodes.1xg2x=bias$barcodes.1xg2x
  barcodes.1xg2xA=bias$barcodes.1xg2xA
  barcodes.2xg5xA=bias$barcodes.2xg5xA
  barcodes.5xg10xA=bias$barcodes.5xg10xA
  barcodes.g10xA=bias$barcodes.g10xA
  
  df=NULL
  for(i in 1:8){
    barcodes=bias[[i]]
    values=data[barcodes,index]
    cat=c(">10x","5-10x","2-5x","1-2x","1-2x A","2-5x A","5-10x A",">10x A")[i]
    if(length(barcodes)<1){
      barcodes=NA
      values=NA
    }
    temp.df=data.frame(barcodes=barcodes,values=values,cat=cat)
    df=rbind(df,temp.df)
  }
  


  return(df)
  if(sum(df$values,na.rm=TRUE)==1){
    print(ggplot(df,aes(x=cat,y=values,category=barcodes))+
            geom_bar(colour="black",stat="identity",size=.3)+
            scale_fill_manual(values=rep("yellow",times=nrow(df)))+
            scale_y_continuous("",limits=c(0,1),breaks=seq(0,1,by=.2),labels=c("0%","20%","40%","60%","80%","100%"))+
            scale_x_discrete(""))
  }
  

  
}

bias.plot(zh19.rbm.data,index=5)
barcodetrackR::barcode_ggheatmap(zh19.rbm.data)



jpeg(file="jd76_lbm_nrbc.jpg")
bias.plot(jd76.lbm.data)
dev.off()
jpeg(file="jd76_rbm_nrbc.jpg")
bias.plot(jd76.rbm.data)
dev.off()
jpeg(file="zk22_lbm_nrbc.jpg")
bias.plot(zk22.lbm.data)
dev.off()
jpeg(file="zh19_rbm_nrbc.jpg")
bias.plot(zh19.rbm.data)
dev.off()
jpeg(file="zh33_lbm_nrbc.jpg")
bias.plot(zh33.lbm.data)
dev.off()
jpeg(file="zh33_rbm_nrbc.jpg")
bias.plot(zh33.rbm.data)
dev.off()


jpeg(file="jd76_lbm_cd34.jpg")
bias.plot(jd76.lbm.data,index=1)
dev.off()
jpeg(file="jd76_rbm_cd34.jpg")
bias.plot(jd76.rbm.data,index=1)
dev.off()
jpeg(file="zk22_lbm_cd34.jpg")
bias.plot(zk22.lbm.data,index=1)
dev.off()
jpeg(file="zh19_rbm_cd34.jpg")
bias.plot(zh19.rbm.data,index=1)
dev.off()
jpeg(file="zh33_lbm_cd34.jpg")
bias.plot(zh33.lbm.data,index=1)
dev.off()
jpeg(file="zh33_rbm_cd34.jpg")
bias.plot(zh33.rbm.data,index=1)
dev.off()

jpeg(file="jd76_lbm_t.jpg")
bias.plot(jd76.lbm.data,index=2)
dev.off()
jpeg(file="jd76_rbm_t.jpg")
bias.plot(jd76.rbm.data,index=2)
dev.off()
jpeg(file="zk22_lbm_t.jpg")
bias.plot(zk22.lbm.data,index=2)
dev.off()
jpeg(file="zh19_rbm_t.jpg")
bias.plot(zh19.rbm.data,index=2)
dev.off()
jpeg(file="zh33_lbm_t.jpg")
bias.plot(zh33.lbm.data,index=2)
dev.off()
jpeg(file="zh33_rbm_t.jpg")
bias.plot(zh33.rbm.data,index=2)
dev.off()

jpeg(file="jd76_lbm_b.jpg")
bias.plot(jd76.lbm.data,index=3)
dev.off()
jpeg(file="jd76_rbm_b.jpg")
bias.plot(jd76.rbm.data,index=3)
dev.off()
jpeg(file="zk22_lbm_b.jpg")
bias.plot(zk22.lbm.data,index=3)
dev.off()
jpeg(file="zh19_rbm_b.jpg")
bias.plot(zh19.rbm.data,index=3)
dev.off()
jpeg(file="zh33_lbm_b.jpg")
bias.plot(zh33.lbm.data,index=3)
dev.off()
jpeg(file="zh33_rbm_b.jpg")
bias.plot(zh33.rbm.data,index=3)
dev.off()

jpeg(file="jd76_lbm_mono.jpg")
bias.plot(jd76.lbm.data,index=4)
dev.off()
jpeg(file="jd76_rbm_mono.jpg")
bias.plot(jd76.rbm.data,index=4)
dev.off()
jpeg(file="zk22_lbm_mono.jpg")
bias.plot(zk22.lbm.data,index=4)
dev.off()
jpeg(file="zh19_rbm_mono.jpg")
bias.plot(zh19.rbm.data,index=4)
dev.off()
jpeg(file="zh33_lbm_mono.jpg")
bias.plot(zh33.lbm.data,index=4)
dev.off()
jpeg(file="zh33_rbm_mono.jpg")
bias.plot(zh33.rbm.data,index=4)
dev.off()


jpeg(file="jd76_lbm_grans.jpg")
bias.plot(jd76.lbm.data,index=5)
dev.off()
jpeg(file="jd76_rbm_grans.jpg")
bias.plot(jd76.rbm.data,index=5)
dev.off()
jpeg(file="zk22_lbm_grans.jpg")
bias.plot(zk22.lbm.data,index=5)
dev.off()
jpeg(file="zh19_rbm_grans.jpg")
bias.plot(zh19.rbm.data,index=5)
dev.off()
jpeg(file="zh33_lbm_grans.jpg")
bias.plot(zh33.lbm.data,index=5)
dev.off()
jpeg(file="zh33_rbm_grans.jpg")
bias.plot(zh33.rbm.data,index=5)
dev.off()

jpeg(file="rq859_cd34.jpg")
bias.plot(rq859.data,index=1)
dev.off()
jpeg(file="rq859_t.jpg")
bias.plot(rq859.data,index=2)
dev.off()
jpeg(file="rq859_b.jpg")
bias.plot(rq859.data,index=3)
dev.off()
jpeg(file="rq859_grans.jpg")
bias.plot(rq859.data,index=4)
dev.off()
jpeg(file="rq859_nrbc.jpg")
bias.plot(rq859.data,index=5)
dev.off()


jpeg(file="rq3600_cd34.jpg")
bias.plot(rq3600.data,index=1)
dev.off()
jpeg(file="rq3600_t.jpg")
bias.plot(rq3600.data,index=2)
dev.off()
jpeg(file="rq3600_b.jpg")
bias.plot(rq3600.data,index=3)
dev.off()
jpeg(file="r3600_grans.jpg")
bias.plot(rq3600.data,index=4)
dev.off()
jpeg(file="r3600_nrbc.jpg")
bias.plot(rq3600.data,index=5)
dev.off()

barcodetrackR::barcode_ggheatmap(rq859.data,cellnote_option = "percents",cellnote_size=0,label_size = 5,n_clones=100)
temp=temp[temp$cat==">10x",]

tempbar=barcodetrackR::barcode_ggheatmap(rq859.data,cellnote_option = "percents",label_size = 5,printtable=TRUE,n_clones=100)[rownames(temp),]
tempbar=tempbar[,colSums(tempbar)>0]
barcodetrackR::barcode_ggheatmap(tempbar,cellnote_option = "percents",label_size = 5)
