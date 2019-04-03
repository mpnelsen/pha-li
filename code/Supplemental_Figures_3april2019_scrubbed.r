###########################################
###########################################
###########################################
#Figure S1
#Sparse Ascomycota w Boots
ml.sparse<-read.tree(file="RAxML_bipartitions.result")
ml.sparse.lad<-ladderize(ml.sparse,FALSE)
plot(ml.sparse.lad,cex=0.35)
node.support(ml.sparse.lad$node.label,cutoff=70,cex=0.4)
#so used this instead, which is ugly
#ml.dense.lad$node.label[ml.dense.lad$node.label<70]<-""
#nodelabels(ml.dense.lad$node.label,adj=c(0.5,0.5),frame="none",bg="white",cex=0.75)
add.scale.bar(x=0.4,y=5,cex=0.5)

###########################################
###########################################
###########################################
#Figure S2
#Sparse Ascomycota w HPD
beast.tree<-read.beast(file="fungi_redo_silog_nonuni_wstart_mcc_medhts.tre")
plot(ladderize(beast.tree,FALSE),cex=0.7,edge.width=0.25,x.lim=c(-100, 800))
HPDbars(ladderize(beast.tree,FALSE),col="darkolivegreen3",lwd=5)
axisPhylo(cex=0.75)

###########################################
###########################################
###########################################
#Figure S3
#Detailed ML Algae w boots
require(phyloch)
require(diversitree)
alg<-read.tree(file="RAxML_bipartitions.result")
alg$tip.label[alg$tip.label %in% "Dictyochloropsis_irregularis"]<-"Symbiochloris_irregularis"
alg$tip.label[alg$tip.label %in% "Dictyochloropsis_sp_MPN124"]<-"Symbiochloris_pauciautosporica_MPN124"
alg$tip.label[alg$tip.label %in% "Dictyochloropsis_sp_MPN168"]<-"Symbiochloris_sp_MPN168"
alg$tip.label[alg$tip.label %in% "Dictyochloropsis_sp_MPN181"]<-"Symbiochloris_sp_MPN181"
alg$tip.label[alg$tip.label %in% "Dictyochloropsis_splendida"]<-"Symbiochloris_pauciautosporica"
alg$tip.label[alg$tip.label %in% "Dictyochloropsis_symbiontica_CAUP"]<-"Symbiochloris_splendida_CAUP"
alg$tip.label[alg$tip.label %in% "Dictyochloropsis_symbiontica_SAG"]<-"Symbiochloris_splendida_SAG"
alg<-ladderize(alg,FALSE)
plot(alg,cex=0.35)
node.support(alg$node.label,cutoff=70,cex=0.4)
add.scale.bar(x=0.4,y=5,cex=0.5)

###########################################
###########################################
###########################################
#Figure S4
require(phyloch)
require(diversitree)
alg<-read.beast(file="mcc.med.hts.tre")
alg$tip.label[alg$tip.label %in% "Dictyochloropsis_irregularis"]<-"Symbiochloris_irregularis"
alg$tip.label[alg$tip.label %in% "Dictyochloropsis_sp_MPN124"]<-"Symbiochloris_pauciautosporica_MPN124"
alg$tip.label[alg$tip.label %in% "Dictyochloropsis_sp_MPN168"]<-"Symbiochloris_sp_MPN168"
alg$tip.label[alg$tip.label %in% "Dictyochloropsis_sp_MPN181"]<-"Symbiochloris_sp_MPN181"
alg$tip.label[alg$tip.label %in% "Dictyochloropsis_splendida"]<-"Symbiochloris_pauciautosporica"
alg$tip.label[alg$tip.label %in% "Dictyochloropsis_symbiontica_CAUP"]<-"Symbiochloris_splendida_CAUP"
alg$tip.label[alg$tip.label %in% "Dictyochloropsis_symbiontica_SAG"]<-"Symbiochloris_splendida_SAG"
plot(ladderize(alg,FALSE),cex=0.7,edge.width=0.25,x.lim=c(-600, 3000))
HPDbars(alg,col="darkolivegreen3",lwd=5)
axisPhylo(cex=0.75)

###########################################
###########################################
###########################################
#Figure S5
require(diversitree)
require(circlize)
alg<-read.nexus(file="mcc.med.hts.tre")
aa<-read.csv(file="all4_organized_for_plotting.updated.csv",stringsAsFactors=FALSE,na.strings="")

aa$TaxToUse[aa$Species %in% "Dictyochloropsis_irregularis"]<-"Symbiochloris"
aa$TaxToUse[aa$Species %in% "Dictyochloropsis_sp_MPN124"]<-"Symbiochloris"
aa$TaxToUse[aa$Species %in% "Dictyochloropsis_sp_MPN168"]<-"Symbiochloris"
aa$TaxToUse[aa$Species %in% "Dictyochloropsis_sp_MPN181"]<-"Symbiochloris"
aa$TaxToUse[aa$Species %in% "Dictyochloropsis_splendida"]<-"Symbiochloris"
aa$TaxToUse[aa$Species %in% "Dictyochloropsis_symbiontica_CAUP"]<-"Symbiochloris"
aa$TaxToUse[aa$Species %in% "Dictyochloropsis_symbiontica_SAG"]<-"Symbiochloris"
alg$tip.label[alg$tip.label %in% "Dictyochloropsis_irregularis"]<-"Symbiochloris_irregularis"
aa$Genus[aa$Species %in% "Dictyochloropsis_irregularis"]<-"Symbiochloris"
aa$Species[aa$Species %in% "Dictyochloropsis_irregularis"]<-"Symbiochloris_irregularis"
alg$tip.label[alg$tip.label %in% "Dictyochloropsis_sp_MPN124"]<-"Symbiochloris_pauciautosporica_MPN124"
aa$Genus[aa$Species %in% "Dictyochloropsis_sp_MPN124"]<-"Symbiochloris"
aa$Species[aa$Species %in% "Dictyochloropsis_sp_MPN124"]<-"Symbiochloris_pauciautosporica_MPN124"
alg$tip.label[alg$tip.label %in% "Dictyochloropsis_sp_MPN168"]<-"Symbiochloris_sp_MPN168"
aa$Genus[aa$Species %in% "Dictyochloropsis_sp_MPN168"]<-"Symbiochloris"
aa$Species[aa$Species %in% "Dictyochloropsis_sp_MPN168"]<-"Symbiochloris_sp_MPN168"
alg$tip.label[alg$tip.label %in% "Dictyochloropsis_sp_MPN181"]<-"Symbiochloris_sp_MPN181"
aa$Genus[aa$Species %in% "Dictyochloropsis_sp_MPN181"]<-"Symbiochloris"
aa$Species[aa$Species %in% "Dictyochloropsis_sp_MPN181"]<-"Symbiochloris_sp_MPN181"
alg$tip.label[alg$tip.label %in% "Dictyochloropsis_splendida"]<-"Symbiochloris_pauciautosporica"
aa$Genus[aa$Species %in% "Dictyochloropsis_splendida"]<-"Symbiochloris"
aa$Species[aa$Species %in% "Dictyochloropsis_splendida"]<-"Symbiochloris_pauciautosporica"
alg$tip.label[alg$tip.label %in% "Dictyochloropsis_symbiontica_CAUP"]<-"Symbiochloris_splendida_CAUP"
aa$Genus[aa$Species %in% "Dictyochloropsis_symbiontica_CAUP"]<-"Symbiochloris"
aa$Species[aa$Species %in% "Dictyochloropsis_symbiontica_CAUP"]<-"Symbiochloris_splendida_CAUP"
alg$tip.label[alg$tip.label %in% "Dictyochloropsis_symbiontica_SAG"]<-"Symbiochloris_splendida_SAG"
aa$Genus[aa$Species %in% "Dictyochloropsis_symbiontica_SAG"]<-"Symbiochloris"
aa$Species[aa$Species %in% "Dictyochloropsis_symbiontica_SAG"]<-"Symbiochloris_splendida_SAG"
aa$Species[!aa$Species %in% alg$tip.label]

rownames(aa)<-aa$Species
tip.states<-as.data.frame(aa$lichen.clade,row.names=aa$Species)
colnames(tip.states)<-"lichen"
#reduce tip labels to higher-level clades
class.high<-as.data.frame(aa$Class.or.Higher.Level,row.names=aa$Species,stringsAsFactors=FALSE)
colnames(class.high)<-"Class.or.Higher.Level"
class.low<-as.data.frame(aa$Genus,row.names=aa$Species,stringsAsFactors=FALSE)
colnames(class.low)<-"Genus"
class.plot<-as.data.frame(aa$TaxToUse,row.names=aa$Species,stringsAsFactors=FALSE)
colnames(class.plot)<-"TaxToUse"
tiplabs<-alg$tip.label
tax2.high<-class.high[match(tiplabs,rownames(class.high)),]
tax2.low<-class.low[match(tiplabs,rownames(class.low)),]
tax2.plot<-class.plot[match(tiplabs,rownames(class.plot)),]
class.high2<-aa[,"Class.or.Higher.Level"]
names(class.high2)<-aa
class.low2<-aa[,"Genus"]

pdf("Fig_3_algae.ages.geoscale.23nov2018.pdf",width=11.5,height=12.5)

trait.plot(ladderize(alg,FALSE),tip.states,cols=list(lichen=c("white","green4")),legend=FALSE,edge.width=1,cex.lab=1,class=tax2.plot,w=1/100)
ht<-max(nodeHeights(alg))

timescale<-read.csv(file="timescale_ics2015_modified.csv",stringsAsFactors=FALSE)

for(x in 1:nrow(timescale)){
	timescale$RGB[x]<-rgb(timescale$Col_R[x]/255,timescale$Col_G[x]/255,timescale$Col_B[x]/255,alpha=0.2)
}
require(phytools)
timey<-timescale[timescale$Type %in% "Period",]
timey<-timey[timey$End<max(nodeHeights(alg)),]


circos.par(gap.degree=0,start.degree = 90,canvas.xlim=c(-1.25,1.25),canvas.ylim=c(-1.25,1.25))
circos.initialize("a", xlim = c(0, 1),sector.width=1) # 'a` just means there is one sector 

#for(x in 1:nrow(timey)){
#	draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(ht-timey$End[x])/ht,col=timey$RGB[x],border=NA)
#}

ht

for(x in 1:nrow(timey)){
	draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(ht-timey$End[x])/ht,rou2=(ht-timey$Start[x])/ht,col=timey$RGB[x],border=NA)
	draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(ht-200)/ht,col=NA,border="gray77",lwd=0.25,lty=3)
	draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(ht-400)/ht,col=NA,border="gray77",lwd=0.25,lty=3)
	draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(ht-600)/ht,col=NA,border="gray77",lwd=0.25,lty=3)
	draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(ht-800)/ht,col=NA,border="gray77",lwd=0.25,lty=3)
	draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(ht-1000)/ht,col=NA,border="gray77",lwd=0.25,lty=3)
	draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(ht-1200)/ht,col=NA,border="gray77",lwd=0.25,lty=3)
	draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(ht-1400)/ht,col=NA,border="gray77",lwd=0.25,lty=3)
}

circos.clear()
par(new = T)
trait.plot(ladderize(alg,FALSE),tip.states,cols=list(lichen=c("white","green4")),legend=FALSE,edge.width=1,cex.lab=1,class=tax2.plot,w=1/100)

#legend("bottomright",title="Node States",c("Non-Lichenized, Stable","Non-Lichenized","Lichenized","Lichenized, Stable"),col=cols[c(1,3,4,2)],pch=19,bty="n")
legend("bottomright",title=expression(bold("Tip States")),"Lichen-Associated Clade",col="green4",cex=1.5,pch=19,bty="n",pt.cex=2)
segments(x0=850,y0=-1950,x1=1925,y1=-1950)

#axis(1,pos=-0.008*ht,lwd=0,at=c(41.457,241.457,441.457,641.4579,841.457,1041.457,1241.457),labels=c("1.4 Ga","1.2 Ga","1 Ga","800 Ma","600 Ma","400 Ma","200 Ma"),lwd.ticks=0,padj=-3.5,cex.axis=0.7)

for(x in 3:17){
	#print(x)
	text(x=ht-rev(timey$Midpoint[x]),y=seq(from=0,to=-7,by=-7/16)[x],labels=timey$Abbrev[x],cex=0.75)
}

dev.off()


###########################################
###########################################
###########################################
#Figure S6
#Detailed ML Ascomycota w boots
ml.dense<-read.tree(file="RAxML_bipartitions.bestw500boots")
ml.dense.lad<-ladderize(ml.dense,FALSE)
plot(ml.dense.lad,cex=0.35)
node.support(ml.dense.lad$node.label,cutoff=70,cex=0.5)
#so used this instead, which is ugly
#ml.dense.lad$node.label[ml.dense.lad$node.label<70]<-""
#nodelabels(ml.dense.lad$node.label,adj=c(0.5,0.5),frame="none",bg="white",cex=0.75)
add.scale.bar(x=0.4,y=2,cex=0.75)


###########################################
###########################################
###########################################
#Figure S7

col2rgb("darkgreen")
trans.dk.green<-rgb(0,100,0,max=255,alpha=125,names="transdkgreen")

timescale<-read.csv(file="timescale_ics2015_modified.csv",stringsAsFactors=FALSE)

for(x in 1:nrow(timescale)){
	timescale$RGB[x]<-rgb(timescale$Col_R[x]/255,timescale$Col_G[x]/255,timescale$Col_B[x]/255,alpha=0.2)
}
require(phytools)
timey<-timescale[timescale$Type %in% "Period",]
timey<-timey[timey$End<500,]
timey<-timey[timey$Start>200,]
timey$End[1]<-200
timey$Start[nrow(timey)]<-500
timey$Midpoint[nrow(timey)]<-timey$End[nrow(timey)]+(timey$Start[nrow(timey)]-timey$End[nrow(timey)])/2


bootswmed<-read.csv(file="corhmm_summary_root_really_fixed.csv",stringsAsFactors=FALSE)


par(mfrow=c(2,2))
dens<-density(bootswmed$MaxAgeStartLich*-1)

#hist AND density - nice. CROWN
hist(bootswmed$MaxAgeEndLich*-1,breaks=seq(-500,-200,by=10),freq=TRUE,col="light grey",xlab="Crown Age (Ma) Earliest LFF",main=NA,ylim=c(0,150),cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
for(x in 1:nrow(timey)){
	rect(xleft=-timey$Start[x],ybottom=-50,xright=-timey$End[x],ytop=150,col=timey$RGB[x],border=NA)
}
#rect(-385,-0.00095,-200,0.025,col="green",border=NA)
hist(bootswmed$MaxAgeEndLich*-1,breaks=seq(-500,-200,by=10),freq=TRUE,col="light grey",xlab="Crown Age (Ma) Earliest LFF",main=NA,ylim=c(0,150),cex.lab=1.5,cex.axis=1.5,cex.main=1.5,add=T)
#polygon(density(bootswmed$MaxAgeEndLich*-1),col=trans.dk.green)
arrows(-245,140,-245,113,length=0.1)
text(x=-245,y=145,cex=1.5,labels="ML Tree")
text(x=-500+(300*0.05),y=143,cex=2,labels="A.)")
for(x in 2:nrow(timey)){
	text(x=-timey$Midpoint[x],y=-3.5,cex=1,col="grey35",labels=timey$Abbrev[x])
}

#hist AND density - nice. STEM
hist(bootswmed$MaxAgeStartLich*-1,breaks=seq(-500,-200,by=10),freq=TRUE,col="light grey",xlab="Stem Age (Ma) Earliest LFF",main=NA,ylim=c(0,150),cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
for(x in 1:nrow(timey)){
	rect(xleft=-timey$Start[x],ybottom=-50,xright=-timey$End[x],ytop=150,col=timey$RGB[x],border=NA)
}
hist(bootswmed$MaxAgeStartLich*-1,breaks=seq(-500,-200,by=10),freq=TRUE,col="light grey",xlab="Stem Age (Ma) Earliest LFF",main=NA,ylim=c(0,150),cex.lab=1.5,cex.axis=1.5,cex.main=1.5,add=T)
#polygon(density(bootswmed$MaxAgeStartLich*-1),col=trans.dk.green)
arrows(-295,140,-295,90,length=0.1)
text(x=-295,y=145,cex=1.5,labels="ML Tree")
text(x=-500+(300*0.05),y=143,cex=2,labels="B.)")
for(x in 2:nrow(timey)){
	text(x=-timey$Midpoint[x],y=-3.5,cex=1,col="grey35",labels=timey$Abbrev[x])
}


bootswsample<-read.csv(file="boots_corhmm_summary_root_really_fixed_combined_correct.csv",stringsAsFactors=FALSE)
timey<-timescale[timescale$Type %in% "Period",]
timey<-timey[timey$End<600,]
timey<-timey[timey$Start>150,]
timey$End[1]<-150
timey$Start[nrow(timey)]<-600
timey$Midpoint[nrow(timey)]<-timey$End[nrow(timey)]+(timey$Start[nrow(timey)]-timey$End[nrow(timey)])/2


#dens<-density(bootswsample$MaxAgeStartLich*-1)

#hist AND density - nice. CROWN
hist(bootswsample$MaxAgeEndLich*-1,breaks=seq(-600,-150,by=10),freq=TRUE,col="light grey",xlab="Crown Age (Ma) Earliest LFF",main=NA,ylim=c(0,625),cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
for(x in 1:nrow(timey)){
	rect(xleft=-timey$Start[x],ybottom=-50,xright=-timey$End[x],ytop=625,col=timey$RGB[x],border=NA)
}
#rect(-385,-0.00095,-200,0.025,col="green",border=NA)
hist(bootswsample$MaxAgeEndLich*-1,breaks=seq(-600,-150,by=10),freq=TRUE,col="light grey",,xlab="Crown Age (Ma) Earliest LFF",main=NA,ylim=c(0,625),cex.lab=1.5,cex.axis=1.5,cex.main=1.5,add=T)
#polygon(density(bootswsample$MaxAgeEndLich*-1),col=trans.dk.green)
arrows(-255,590,-255,450,length=0.1)
text(x=-255,y=605,cex=1.5,labels="ML Tree (Mean)")
text(x=-590+(300*0.05),y=600,cex=2,labels="C.)")
for(x in 2:nrow(timey)){
	text(x=-timey$Midpoint[x],y=-14,cex=1,col="grey35",labels=timey$Abbrev[x])
}

#hist AND density - nice. STEM
hist(bootswsample$MaxAgeStartLich*-1,breaks=seq(-600,-150,by=10),freq=TRUE,col="light grey",xlab="Stem Age (Ma) Earliest LFF",main=NA,ylim=c(0,625),cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
for(x in 1:nrow(timey)){
	rect(xleft=-timey$Start[x],ybottom=-50,xright=-timey$End[x],ytop=625,col=timey$RGB[x],border=NA)
}
hist(bootswsample$MaxAgeStartLich*-1,breaks=seq(-600,-150,by=10),freq=TRUE,col="light grey",xlab="Stem Age (Ma) Earliest LFF",main=NA,ylim=c(0,625),cex.lab=1.5,cex.axis=1.5,cex.main=1.5,add=T)
#polygon(density(bootswsample$MaxAgeStartLich*-1),col=trans.dk.green)
arrows(-305,590,-305,575,length=0.1)
text(x=-305,y=605,cex=1.5,labels="ML Tree (Mean)")
text(x=-590+(300*0.05),y=600,cex=2,labels="D.)")
for(x in 2:nrow(timey)){
	text(x=-timey$Midpoint[x],y=-14,cex=1,col="grey35",labels=timey$Abbrev[x])
}



#Figure S8
#SI FIGURE
####Make SI w losses
col2rgb("darkgreen")
#trans.dk.green<-rgb(0,100,0,max=255,alpha=125,names="transdkgreen")

bootswmed<-read.csv(file="corhmm_summary_root_really_fixed.csv",stringsAsFactors=FALSE)
bootswsample<-read.csv(file="boots_corhmm_summary_root_really_fixed_combined_correct.csv",stringsAsFactors=FALSE)

par(mfrow=c(2,2))
#hist AND density - nice.
hist(bootswmed$NoGains,breaks=seq(0,25,by=1),ylim=c(0,100),freq=TRUE,col="light grey",xlab="Transitions To LFF",main=NA,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
#polygon(density(bootswmed$NoGains),col=trans.dk.green)
arrows(16.5,97,16.5,45,length=0.1)
text(x=16.5,y=100,cex=1.5,labels="ML Tree")
text(x=20*0.05,y=98,cex=2,labels="A.)")

#hist AND density - nice.
hist(bootswmed$NoLosses,breaks=seq(0,40,by=1),ylim=c(0,100),freq=TRUE,col="light grey",xlab="Transitions From LFF",main=NA,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
#polygon(density(bootswmed$NoLosses),col=trans.dk.green)
arrows(4.5,97,4.5,55,length=0.1)
text(x=7.5,y=100,cex=1.5,labels="ML Tree")
text(x=20*0.05,y=98,cex=2,labels="B.)")

#hist AND density - nice.
hist(bootswsample$NoGains,breaks=seq(0,25,by=1),ylim=c(0,600),freq=TRUE,col="light grey",xlab="Transitions To LFF",main=NA,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
#polygon(density(bootswsample$NoGains),col=trans.dk.green)
arrows(16.5,590,16.5,435,length=0.1)
text(x=16.5,y=605,cex=1.5,labels="ML Tree (Mean)")
text(x=20*0.05,y=598,cex=2,labels="C.)")

#hist AND density - nice.
hist(bootswsample$NoLosses,breaks=seq(0,50,by=1),ylim=c(0,600),freq=TRUE,col="light grey",xlab="Transitions From LFF",main=NA,cex.lab=1.5,cex.axis=1.5,cex.main=1.5)
#polygon(density(bootswsample$NoLosses),col=trans.dk.green)
arrows(5,590,5,530,length=0.1)
text(x=14.5,y=605,cex=1.5,labels="ML Tree (Mean)")
text(x=20*0.05,y=598,cex=2,labels="D.)")



