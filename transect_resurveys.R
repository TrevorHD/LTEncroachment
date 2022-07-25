resurv <- read.csv("https://github.com/TrevorHD/LTEncroachment/raw/master/Data/creosote_transect_resurvey.csv")

par(mfrow=c(1,2),mar=c(4,4,1,1))
plot(resurv$Quad[resurv$Year==2001 & resurv$Transect=="East"],
     resurv$LATR_percent[resurv$Year==2001 & resurv$Transect=="East"],
     type="l",lwd=2,col="gray",ylim=c(0,50),
     xlab="Meter",ylab="Percent cover")
points(resurv$Quad[resurv$Year==2013 & resurv$Transect=="East"]+5,
       resurv$LATR_percent[resurv$Year==2013 & resurv$Transect=="East"],
       type="l",lwd=2,col="black")
title("A",adj=0,font=3)
plot(resurv$Quad[resurv$Year==2001 & resurv$Transect=="West"],
     resurv$LATR_percent[resurv$Year==2001 & resurv$Transect=="West"],
     type="l",lwd=2,col="gray",ylim=c(0,50),
     xlab="Meter",ylab="Percent cover")
points(resurv$Quad[resurv$Year==2013 & resurv$Transect=="West"]+5,
       resurv$LATR_percent[resurv$Year==2013 & resurv$Transect=="West"],
       type="l",lwd=2,col="black")
title("B",adj=0,font=3)
legend("topright",legend=c("2001","2013"),bty="n",lwd=2,col=c("gray","black"))

##how many quads had zero creosote in 2001 and non-zero in 2013
resurv$trans_quad <- interaction(resurv$Transect,resurv$Quad)
zero2001 <- resurv[which(resurv$LATR_percent[resurv$Year==2001]==0),]$trans_quad
nonzero2013 <- resurv[resurv$trans_quad%in%zero2001 & resurv$Year==2013,]$LATR_percent>0
zero_to_nonzero <- sum(nonzero2013,na.rm=T)/length(na.omit(nonzero2013))

##now the reverse: what fraction of quads with nonzero in 2001 had zero in 2013?
nonzero2001 <- resurv[which(resurv$LATR_percent[resurv$Year==2001]>0),]$trans_quad
zero2013 <- resurv[resurv$trans_quad%in%nonzero2001 & resurv$Year==2013,]$LATR_percent==0
nonzero_to_zero <- sum(zero2013,na.rm=T)/length(na.omit(zero2013))

## what was the overall mean cover in both years?
mean_cover <- resurv %>% 
  group_by(Transect,as.factor(Year)) %>% 
  summarise(mean(LATR_percent,na.rm=T))
