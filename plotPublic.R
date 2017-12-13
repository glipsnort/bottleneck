{
source("plotSim.R")
scale = 78.4
wrong_anc = .023
#scale = scale * 1.25  # 25% higher mutation rate
colors = c("#00008F", "#005AFF", "#23FFDC", "#ECFF13", "#FF4A00", "#800000")

# 1kg data
d = read.delim("freq_1kG_hist.txt")
afr = d$count[d$freq > 0.01 & d$freq < 0.99]
freq = d$freq[d$freq > 0.01 & d$freq < 0.99]
afrSeg = sum(d$count[d$freq > 0.6 & d$freq < 0.7])

plot(freq,afr,col="black",main="1000 Genomes data",pch=16,xlab="Derived allele frequency", 
      ylab="Number of variants")

# 100kya 4k
d = read.delim("100kya_4k.txt")
plotSim(d, "Bottleneck 100 kya (pop size = 4k)", "100 kya simulation")

# 250kya 8k
#d = read.delim("250kya_8k.txt")
#plotSim(d, "Bottleneck 250 kya (pop size = 8k)", "250 kya simulation")

}
