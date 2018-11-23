#Figure 1

require(diversitree)
require(circlize)
require(phytools)
tip.states<-read.csv(file="~/Desktop/deep_lichens/asco_phylo_all_medhts/bad_idea_lichenization_for_corhmm.csv",stringsAsFactors=FALSE)
tip.states<-as.data.frame(tip.states$lichen,row.names=tip.states$edited_name)
colnames(tip.states)<-"lichen"
sum.tab<-read.csv(file="~/Desktop/deep_lichens/asco_phylo_all_medhts/MLTREE_2Rate_rootreallyfixed_corhmm_states.csv",stringsAsFactors=FALSE)
colnames(sum.tab)<-c("Node","0R1","1R1","0R2","1R2","Age")
sum.tab.mat<-as.matrix(sum.tab)
cols<-c("gray75","green4","gray95","darkolivegreen1")
tr<-read.tree(file="~/Desktop/deep_lichens/asco_phylo_all_medhts/MLTREE_2Rate_rootreallyfixed_corhmm_states.tre")
#reduce tip labels to higher-level clades
class<-read.csv(file="~/Desktop/deep_lichens/asco_phylo_all_medhts/bad_idea_lichenization_for_corhmm_w_higher.csv",stringsAsFactors=FALSE)
rownames(class)<-class$edited_name
tiplabs<-tr$tip.label
tax2<-class[match(tiplabs,rownames(class)),]
class2<-class[,"Class"]

pdf("~/Desktop/deep_lichens/Fig_1_lichenization.MLTree.2rate.tip.cols.lich.rootreallyfixed.geoscale_21nov2018.pdf",width=11,height=11)

trait.plot(ladderize(tr,FALSE),tip.states,cols=list(Lichen=c("gray75","green4")),legend=FALSE,edge.width=0.25,cex.lab=0.75,class=tax2$Class)
#nodelabels(frame="circle",col=cols[as.numeric(tr$node.label)],bg=cols[as.numeric(tr$node.label)],cex=0.15)
nodelabels(pie=sum.tab.mat[,2:5],piecol=cols,bg=cols,cex=0.25)
legend("topleft",title="Character State",c("Non-Lichenized, Stable","Non-Lichenized","Lichenized","Lichenized, Stable"),fill=cols[c(1,3,4,2)])

timescale<-read.csv(file="/Volumes/welwitschia/li'l rascal iv/Users/matthewnelsen/Documents/projects/carboniferous_buildup/carboniferous_revisions/25nov2015_data_pull_sap/timescale_ics2015_modified.csv",stringsAsFactors=FALSE)

for(x in 1:nrow(timescale)){
	timescale$RGB[x]<-rgb(timescale$Col_R[x]/255,timescale$Col_G[x]/255,timescale$Col_B[x]/255,alpha=0.2)
}
require(phytools)
timey<-timescale[timescale$Type %in% "Period",]
timey<-timey[timey$End<max(nodeHeights(tr)),]


#add age rings
trait.plot(ladderize(tr,FALSE),tip.states,cols=list(Lichen=c("gray75","green4")),legend=FALSE,edge.width=0.25,cex.lab=0.75,class=tax2$Class)
ht<-max(nodeHeights(tr))
circos.par(gap.degree=0,start.degree = 90,canvas.xlim=c(-1.25,1.25),canvas.ylim=c(-1.25,1.25))
circos.initialize("a", xlim = c(0, 1),sector.width=1) # 'a` just means there is one sector 

ht

for(x in 1:nrow(timey)){
	draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(ht-timey$End[x])/ht,rou2=(ht-timey$Start[x])/ht,col=timey$RGB[x],border=NA)
	draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(ht-100)/ht,col=NA,border="gray77",lwd=0.25,lty=3)
	draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(ht-200)/ht,col=NA,border="gray77",lwd=0.25,lty=3)
	draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(ht-300)/ht,col=NA,border="gray77",lwd=0.25,lty=3)
	draw.sector(center=c(0,0),start.degree=0,end.degree=360,rou1=(ht-400)/ht,col=NA,border="gray77",lwd=0.25,lty=3)
}

circos.clear()
par(new = T)
trait.plot(ladderize(tr,FALSE),tip.states,cols=list(Lichen=c("gray75","green4")),legend=FALSE,edge.width=0.25,cex.lab=0.75,class=tax2$Class)
nodelabels(pie=sum.tab.mat[,2:5],piecol=cols,bg=cols,cex=0.25,border=FALSE)
legend(x=350,y=-525,title=expression(bold("           Tip States")),c("Non-Lichenized","Lichenized"),col=cols[c(1,2)],pch=19,bty="n",cex=1,pt.cex=2)
segments(x0=350,y0=-550,x1=600,y1=-550)
legend(x=350,y=-375,title=expression(bold("Node States")),c("Non-Lichenized, Stable","Non-Lichenized","Lichenized","Lichenized, Stable"),col=cols[c(1,3,4,2)],pch=19,bty="n",cex=1,pt.cex=2)
segments(x0=350,y0=-405,x1=600,y1=-405)

for(x in 2:10){
	#print(x)
	text(x=ht-rev(timey$Midpoint[x]),y=seq(from=0,to=-7,by=-7/9)[x],labels=timey$Abbrev[x],cex=0.75)
}

dev.off()

#Figure 2

col2rgb("darkgreen")
trans.dk.green<-rgb(0,100,0,max=255,alpha=125,names="transdkgreen")

timescale<-read.csv(file="/Volumes/welwitschia/li'l rascal iv/Users/matthewnelsen/Documents/projects/carboniferous_buildup/carboniferous_revisions/25nov2015_data_pull_sap/timescale_ics2015_modified.csv",stringsAsFactors=FALSE)

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


bootswmed<-read.csv(file="~/Desktop/deep_lichens/asco_phylo_all_medhts/individual_bs_trees/corhmm_summary_root_really_fixed.csv",stringsAsFactors=FALSE)


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


bootswsample<-read.csv(file="~/Desktop/deep_lichens/asco_phylo_ml_100medhts/asco_phylo_boot_25medhts/boots_corhmm_summary_root_really_fixed_combined_correct.csv",stringsAsFactors=FALSE)
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






#Figure 3
require(diversitree)
require(circlize)
alg<-read.nexus(file="~/Desktop/deep_lichens/algae_inflog_unun_noi_emp/mcc.med.hts.tre")
aa<-read.csv(file="~/Documents/all4_organized_for_plotting.updated.csv",stringsAsFactors=FALSE,na.strings="")

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

pdf("~/Desktop/deep_lichens/Fig_3_algae.ages.geoscale.23nov2018.pdf",width=11.5,height=12.5)

trait.plot(ladderize(alg,FALSE),tip.states,cols=list(lichen=c("white","green4")),legend=FALSE,edge.width=1,cex.lab=1,class=tax2.plot,w=1/100)
ht<-max(nodeHeights(alg))



timescale<-read.csv(file="/Volumes/welwitschia/li'l rascal iv/Users/matthewnelsen/Documents/projects/carboniferous_buildup/carboniferous_revisions/25nov2015_data_pull_sap/timescale_ics2015_modified.csv",stringsAsFactors=FALSE)

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