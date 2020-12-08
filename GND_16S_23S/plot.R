require(ggplot2);require(scales)

require(chngpt)
a= read.csv('combined_data_qdist75_weighted.txt',sep=" ")

fit=chngptm(formula.1=qdist~1, formula.2=~mean_gnd, data=a[a$both_unresolved<130,] , type="stegmented",
            family="gaussian")
summary(fit,show.slope.post.threshold=TRUE,verbose = TRUE)
summary(lm(qdist~mean_gnd,data=a[a$both_unresolved<130&a$mean_gnd<fit[[2]][5],]))
summary(lm(qdist~mean_gnd,data=a[a$both_unresolved<130&a$mean_gnd>fit[[2]][5],]))

print(paste(fit[[1]]$coefficients[2],fit[[1]]$coefficients[1],sep=" x+ "))
print(paste(fit[[1]]$coefficients[4]+fit[[1]]$coefficients[2],fit[[1]]$coefficients[4]*-fit[[2]][5] +fit[[1]]$coefficients[1]+fit[[1]]$coefficients[3],sep=" x+ "))

ggplot(aes(x=mean_gnd/100,y=qdist), data=a[a$both_unresolved<130,])+
  geom_point(alpha=1/3,color="blue")+
  theme_classic()+xlab("median pairwise GND")+ylab("quartet distance (>0.75 BS)")+
  geom_smooth(aes(group = mean_gnd<fit[[2]][5]),method="lm",se=F,color="red")+
  annotate(geom="text", label=paste(round(fit[[1]]$coefficients[2],5),round(fit[[1]]$coefficients[1],5),sep=" x+ "),x = 0.081,y=0.5)+
  annotate(geom="text", label=paste(round(fit[[1]]$coefficients[4]+fit[[1]]$coefficients[2],5),round(fit[[1]]$coefficients[4]*-fit[[2]][5] +fit[[1]]$coefficients[1]+fit[[1]]$coefficients[3],5),sep=" x+ "),y = 0.05,x=0.2)+
  scale_x_continuous(labels=percent) #+ scale_y_log10()+coord_cartesian(ylim=c(0.02,1.1))

ggsave("combined_data.pdf",width=7,height = 4)
