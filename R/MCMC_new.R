
source("/groups/reich/jp203/projects/mtdna/run_all/funcs_new.R")

a = commandArgs(T)
index1 = a[1]
index2 = a[2]
cprior = as.numeric(a[3])
fpropose = as.numeric(a[4])
tpropose = as.numeric(a[5])
d = read.table(a[6], as.is= T, head = T)
outfile = a[7]

ca = as.numeric(d[index1,3:ncol(d)])
cb = as.numeric(d[index2,3:ncol(d)])
p = ca+cb
ca = ca[p>0]
cb = cb[p>0]


init_theta = list( rep(1/length(ca), length(ca)), log(1), rep(1/length(ca), length(ca))) 
dprior = rep(0.5, length(ca))
cpropose = 0.4
tpropose = 0.00001
N = 50000
Nburn = 2000
Nsamp = 10


#outfile = paste(d[index1,1], "_", d[index2,1], ".out", sep = "")
samps = run_mcmc(ca, cb, init_theta, dprior, cprior, cpropose, tpropose, fpropose, N, Nburn, Nsamp, outfile)
c = chisq.test(rbind(ca, cb))
c = c$statistic
f = mean(samps[[1]] >=c)
s = sum(ca >0)
fs = mean(samps[[2]] <= s)
sink(file = outfile, append = T)
cat("truec", c, "\n")
cat("p-value", f, "\n")
cat("trues", s, "\n")
cat("p-value_s", fs, "\n")
if (f ==0){ f = 1/( (N-Nburn)/Nsamp)}
if (fs ==0){ fs = 1/( (N-Nburn)/Nsamp)}
cp = 1-pchisq( -2*( log(f) + log(fs)), df = 4)
cat("combined_p", cp, "\n") 
sink()


