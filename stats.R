getwd()
data<-read.table("2L_bos", header=T)	#### user edit "file_for_R.txt"

avg<-mean(data$OneRatio)
std<-sd(data$OneRatio)
N<-nrow(data)
contig_col<-data$Contig

p<-numeric()
for (x in data$OneRatio) p<-append(p,pnorm(x, avg, std))	#calculate pnorm (pval)
plist<-data.frame(binID=data[,4],Position=data[,5],OneRatio=data[,6],pval=p )	#add pnorm (pval) to dataframe

sorted<-plist[with(plist, order(pval)),]	#sort dataframe by pval
#sorted: binID	Position	OneRatio	pval
fdrcalc<-function(i,N) {	#calculate FDR
	i*(0.05/N)	###### FALSE DISCOVERY RATE HERE #####
}
fdr<-numeric()
for (i in 1:N) fdr<-append(fdr,fdrcalc(i,N))	#associate sorted pvals with equivalent FDR
combined<-cbind(sorted,fdr)

#combined: binID	Position	OneRatio	pval[4]	FDR[5]
sig<-character()
dif<-numeric()
for (i in 1:N) dif<-append(dif,(combined[i,5] - combined[i,4]))	#check if FDR>pval or FDR<pval via subtraction
for (x in dif) if (x < 0) sig<-append(sig,"n") else sig<-append(sig,"y")	#if FDR-pval is negative, not significant / if FDR-pval is positive, significant
combined<-cbind(combined,sig)

combined<-cbind(contig_col, combined)

final<-combined[with(combined, order(Position)),]

write.table(final, file="2L.test", sep="\t", col.names=TRUE, row.names=FALSE)	#### user edit "file.txt"