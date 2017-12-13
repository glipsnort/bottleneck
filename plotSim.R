plotSim = function(d, m, l) {
  print(m)
  raw = scale * d$raw_counts[d$derived_freq > 0.01 & d$derived_freq < 0.99]
  raw = (1-wrong_anc) * raw + wrong_anc * rev(raw)
  pre = scale * d$pre_bneck[d$derived_freq > 0.01 & d$derived_freq < 0.99]
  pre = (1-wrong_anc) * pre + wrong_anc * rev(pre)
  post = raw - pre
  thisSeg = sum(raw[d$derived_freq > 0.6 & d$derived_freq < 0.7])
  thisSeg_pre = sum(pre[d$derived_freq > 0.6 & d$derived_freq < 0.7])
  bscale = 1 + (afrSeg  - thisSeg) / thisSeg_pre
  print(paste("ancestral multiplier", bscale))
  scaled = pre * bscale + post
  plot(freq,afr,col="black",main=m,xlab="Derived allele frequency", 
      ylab="Number of variants", lwd=2, type="l")
  points(freq, scaled ,col="red", type="l", lwd=2)
  points(c(0,1), c(0,0),type="l")
  legend("topright", c("1000 Genomes data", l), col=c("black", "red"), lwd=2)

}