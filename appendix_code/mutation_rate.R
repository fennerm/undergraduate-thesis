library(boot)
library(perm)

# Simple mutation rate calculation.
mut_rate <- function(muts, gen, nuc) {
  summed = sum(muts)
  summed/(gen*nuc)
}

# Indexed version of mutation rate calculation for bootstrapping confidence intervals
mut_rate_indexed <- function(data, indices, gen, nuc) {
  d = data[indices]
  summed = sum(data[indices])
  return((summed)/(gen*nuc))
}

# Bootstrap mutation rate estimate (10000 replicates)
bootstrap_CI <- function(data, gen=gen, nuc) {
  bt = boot(data, statistic=mut_rate_indexed, R=10000, gen=gen,nuc=nuc)
  boot.ci(bt)$basic
}

# Compare mutation rates using permutation
permute_diffs <- function(data1, data2, gen1, gen2, nuc1, nuc2) {
  grouped.dat = data.frame(af=c(data1, data2), group=as.numeric(c(rep(1, each=length(data1)), rep(2, each=length(data2)))))

  # Actual difference in mutation rate estimates
  t0 = mut_rate(grouped.dat[which(grouped.dat$group==1), "af"], nuc=nuc1, gen=gen1) - mut_rate(grouped.dat[which(grouped.dat$group==2), "af"], nuc=nuc2, gen=gen2)

  b5 = numeric(10000)
  for (i in 1:10000) {
    # Permute group labels
    grouped.dat$group = sample(grouped.dat$group, length(grouped.dat$group))
    # Recalculate Mutation Rate    
    b5[i] = mut_rate(grouped.dat[which(grouped.dat$group==1), "af"],     nuc=nuc1, gen=gen1) - mut_rate(grouped.dat[which(grouped.dat$group==2), "af"], nuc=nuc2, gen=gen2) 
  }
  mean(abs(b5) > abs(t0))
}

# Compare indel:bs mutation rate ratios by permuting variants
permute_bs_ind <- function(data1, data2, type1, type2, gen1, gen2, nuc1, nuc2) {
  grouped.dat = data.frame(af=c(data1, data2), type=c(type1, type2), group=as.numeric(c(rep(1, each=length(data1)), rep(2, each=length(data2)))))
  # Compute base substitution mutation rates
  bs1 = mut_rate(grouped.dat[which((grouped.dat$group==1) & (grouped.dat$type=="bs")), "af"], nuc=nuc1, gen=gen1) 
  bs2 = mut_rate(grouped.dat[which((grouped.dat$group==2) & (grouped.dat$type=="bs")), "af"], nuc=nuc2, gen=gen2)
  # Compute indel mutation rates
  ind1 = mut_rate(grouped.dat[which((grouped.dat$group==1) & (grouped.dat$type=="ind")), "af"], nuc=nuc1, gen=gen1) 
  ind2 = mut_rate(grouped.dat[which((grouped.dat$group==2) & (grouped.dat$type=="ind")), "af"], nuc=nuc2, gen=gen2
  # Calculate observed difference in indel:bs
  t0 = (bs1/ind1) - (bs2/ind2)
  tp = numeric(10000)

  for (i in 1:10000) {
    grouped.dat$group = sample(grouped.dat$group, length(grouped.dat$group))
    bs1 = mut_rate(grouped.dat[which((grouped.dat$group==1) & (grouped.dat$type=="bs")), "af"], nuc=nuc1, gen=gen1) 
    bs2 = mut_rate(grouped.dat[which((grouped.dat$group==2) & (grouped.dat$type=="bs")), "af"], nuc=nuc2, gen=gen2)
    ind1 = mut_rate(grouped.dat[which((grouped.dat$group==1) & (grouped.dat$type=="ind")), "af"], nuc=nuc1, gen=gen1) 
    ind2 = mut_rate(grouped.dat[which((grouped.dat$group==2) & (grouped.dat$type=="ind")), "af"], nuc=nuc2, gen=gen2)
    tp[i] = (bs1/ind1) - (bs2/ind2)
  }
  mean(abs(tp) > abs(t0))
}
