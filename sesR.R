# SES normalization 
# Normalization, bias correction, and peak calling for ChIP-seq
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC3342857/

# a package using this method
# Differential peak calling of ChIP-seq signals with replicates with THOR
# https://www.ncbi.nlm.nih.gov/pmc/articles/PMC5175345/#B42

# specific for histone modifications


library(csaw)
library(readxl)

standard.chr <- paste0("chr", c("I","II","III","IV","V","VI","VII","VIII","IX","X","XI","XII","XIII","XIV","XV","XVI"))

# this is for PE ChIP-seq, not for RNA-seq or TL-seq

param <- readParam(minq=30, restrict=standard.chr,
                   max.frag=1000, pe="both", dedup=TRUE, BPPARAM=MulticoreParam(workers = 8))
# -----


# -----
# functions
# -----
calculate_factor <- function(dat, name){
  dat = dat[order(dat[,1]),]
  dat$rank = seq(1:nrow(dat))
  dat$bin_pct = dat$rank/nrow(dat)
  dat$dat1_tag_pct = cumsum(dat[,1])/sum(dat[,1])
  dat$dat2_tag_pct = cumsum(dat[,2])/sum(dat[,2])
  dat$max = abs(dat$dat1_tag_pct-dat$dat2_tag_pct)
  maxpq = dat[which.max(dat$max),]
  png(paste(name,'.png'), width=800, height=600, unit='px')
  plot(dat$bin_pct, dat$dat1_tag_pct, type='line',xlim=c(0,1), col='red')
  points(dat$bin_pct, dat$dat2_tag_pct, type='line', col='black')
  abline(v=maxpq$bin_pct,lty=2)
  dev.off()
  k = maxpq$rank
  factor = sum(dat[,1][1:k])/sum(dat[,2][1:k])
  total_scale = sum(dat[,1])/sum(dat[,2])
  return(list('ses'=factor, 'sds'=total_scale))
  #
  
}

#only do count1#
count1 = function(bamfiles, width, name){
  print(bamfiles)
  count_data = windowCounts(bamfiles, width=width, bin=TRUE, param=param, filter=0)
  save(count_data, file=paste0(name, '_count.RData'))
}

main = function(count.marker, count.h3, width){
  l = list()
  for(i in 1:dim(count.marker)[2]){
    dat = data.frame(count.marker[,i], count.h3[,i])
    factor = calculate_factor(dat,width)
    l[[i]] = factor
  }
  return(l)
}
# -----

# -----
# input bam
# -----

#label=read_excel('/Volumes/Data1/ChIP_seq_11_2018/sequencing\ libraries\ info.xlsx')

#bamdir='/Volumes/Data1/ChIP_seq_11_2018/bam/remove_duplication'
#setwd(bamdir)
#bam.files=list.files(bamdir)[grepl('.bam.gz$', list.files(bamdir))]

dir="/mnt/pci-0000:05:00.0-sas-phy1-lun-0/ruofan/hdaproject/2018novchip/bam/neededbam"
setwd(dir)
f=list.files()
bams = f[grepl('.bam$',f)]
wth3.bams=bams[grepl('wt_h3_',bams)]
count1(wth3.bams,1000,"wt_h3_bin1000")

hdah3.bams=bams[grepl('hda_h3_',bams)]
count1(hdah3.bams,1000,"hda_h3_bin1000")

wt18.bams=bams[grepl('wt_k18_',bams)]
count1(wt18.bams,1000,"wt_k18_bin1000")

hda18.bams=bams[grepl('hda_k18_',bams)]
count1(hda18.bams,1000,"hda_k18_bin1000")


# -----


# ------
# count
# ------

l = list()
width=1000
for(width in c(1000)){
  load(paste0('wt_h3_bin',width,'_count.RData'))
  count.h3 = count_data
  load(paste0('wt_k18_bin',width,'_count.RData'))
  count.marker = count_data
  count = data.frame(assay(count.marker))
  count.sum=data.frame(count=rowSums(count))
  counth3 = data.frame(assay(count.h3))
  count.h3.sum=data.frame(count=rowSums(counth3))
  dat = data.frame(count.sum$count, count.h3.sum$count)
  factor = calculate_factor(dat,width)
  }
factor

#step2 tmm use RDATA
# TMM normalization 
# 1. count in bin = 2 kb for both histone markers and H3 control
# 2. use SES method to norm. subtraction: histone markers (x) = x - x*sesFactor 
# 3. use TMM from edgeR to get TMM factor
# 4. calculate norm.lib.size and final factor
# the final factor should be used after H3 total norm. 

library(csaw)
library(readxl)
library(edgeR)

file1='/mnt/pci-0000:05:00.0-sas-phy1-lun-0/ruofan/hdaproject/2018novchip/bam/neededbam/hda_k18_bin1000_count.RData'
file2='/mnt/pci-0000:05:00.0-sas-phy1-lun-0/ruofan/hdaproject/2018novchip/bam/neededbam/hda_h3_bin1000_count.RData'
file3='/mnt/pci-0000:05:00.0-sas-phy1-lun-0/ruofan/hdaproject/2018novchip/bam/neededbam/wt_k18_bin1000_count.RData'
file4='/mnt/pci-0000:05:00.0-sas-phy1-lun-0/ruofan/hdaproject/2018novchip/bam/neededbam/wt_h3_bin1000_count.RData'


load(file1)
counts.hdamarker = count_data
load(file2)
counts.hdah3 = count_data
load(file3)
counts.wtmarker = count_data
load(file4)
counts.wth3 = count_data

counts.hdamarker.df = data.frame(hdamarker=rowSums(data.frame(assay(counts.hdamarker))))
counts.hdah3.df = data.frame(hdah3=rowSums(data.frame(assay(counts.hdah3))))
counts.wtmarker.df = data.frame(wtmarker=rowSums(data.frame(assay(counts.wtmarker))))
counts.wth3.df = data.frame(wth3=rowSums(data.frame(assay(counts.wth3))))

#input scale:factor1 is treatment,2 is control
factor1 = 0.4788616
factor2 = 0.4489195 
# subtract or ratio with SES factor

region = data.frame(rowRanges(counts.h3))[,1:3]

i = 1
new_sub = counts.marker.df[,i] - (factor*counts.h3.df[,i])



# TMM factor

new_sub = data.frame(kda=(counts.hdamarker.df[,i] - (factor1*counts.hdah3.df[,i])),
                     wt=counts.wtmarker.df[,i] - (factor2*counts.wth3.df[,i]))

new_sub[new_sub<0] = 0
keep=rowSums(new_sub>0)>0

# use subtraction for following analysis
new_sub = new_sub[keep,]

y = DGEList(new_sub)
y <- calcNormFactors(y)

norm.lib.size = y$samples$lib.size*y$samples$norm.factors
tmm.factor = norm.lib.size/mean(norm.lib.size)

y












