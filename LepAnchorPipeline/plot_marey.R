library(scales)

map=read.table(gzfile("marey.data.gz")) 

NN=unique(map$V3)
MM=unique(map$V4)

cl=rainbow(length(MM))


for(lg in NN) {
	print(lg)
	png(paste0("marey", lg, ".png"))
	ymax = 0
	for (m in MM) {
		map_tmp = map[map$V3==lg & map$V4==m,]
		ymax = max(ymax, max(map_tmp$V5))
	}

	index = 1
	for (m in MM) {
		map_tmp = map[map$V3==lg & map$V4==m,]
		if (glm(map_tmp$V2 ~ map_tmp$V5)$coefficients[2] < 0)
			map_tmp$V5 = max(map_tmp$V5) - map_tmp$V5
		if (index == 1)
			plot(map_tmp$V2, map_tmp$V5, xlab="Position (Mb)",ylab="Recombination Distance",xaxt="n", main=paste0("LG", lg), col=alpha(cl[index],0.2), pch=20,cex=1.5, ylim=c(0,ymax))
		else
			points(map_tmp$V2, map_tmp$V5, col=alpha(cl[index],0.2), pch=20, cex=1.5)
		index = index + 1
	}

	agp <- read.table(paste0("chr", lg, ".agp"))
	agp=agp[agp$V1==paste0("LG",lg) & agp$V5=="W",c(1:3)]

	segments(agp$V2,2*c(1:2),agp$V3,2*c(1:2),col=c("red","blue"),lwd=2,lend=1)
	axis(1,seq(0,50000000,5000000),seq(0,50,5))
	segments(agp$V2,0,agp$V2,ymax,lwd=0.5,col="darkgray")
	dev.off()
}

