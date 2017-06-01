data = read.table('/Users/ericthered/BC15/results/pipeline_v3/2L.test', header=T)
pdf(file="/Users/ericthered/BC15/results/pipeline_v3/2L.pdf", height = 5, width = 12)
plot(data$OneRatio~data$Position, ylim=c(0,1),
	xlab="",type="l",frame.plot=T, lwd=2, yaxt="n", ylab="", main="")
axis(2,at=c(0.2,0.4,0.6,0.8))
#abline(a=0.5,b=0, lty=2)
#segments(13.35,.81,14.3,.81)
#segments(15.65,.81,16.2,.81)
#segments(19.85,.81,19.98,.81)
#mtext("*",side=3,at=13.8, line=-3.5)
#mtext("*",side=3,at=15.93, line=-3.5)
#mtext("*",side=3,at=19.87, line=-3.5)
dev.off()





mtext("*",side=1,at=12, line=-0.7)














