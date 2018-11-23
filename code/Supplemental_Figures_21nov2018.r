#Figure S1
#Detailed ML Algae w boots
require(phyloch)
require(diversitree)
alg<-read.tree(file="~/Desktop/deep_lichens/deep_lichens_to_folks/stuff/algae/raxml/RAxML_bipartitions.result")
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

#Figure S2
require(phyloch)
require(diversitree)
alg<-read.beast(file="~/Desktop/deep_lichens/algae_inflog_unun_noi_emp/mcc.med.hts.tre")
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

#Figure S3
#Sparse Ascomycota w Boots
ml.sparse<-read.tree(file="~/Desktop/deep_lichens/deep_lichens_to_folks/stuff/fungi_restricted/new_analyses_wo_bad_dactylina_seq/raxml_FA/RAxML_bipartitions.result")
ml.sparse.lad<-ladderize(ml.sparse,FALSE)
plot(ml.sparse.lad,cex=0.35)
node.support(ml.sparse.lad$node.label,cutoff=70,cex=0.4)
#so used this instead, which is ugly
#ml.dense.lad$node.label[ml.dense.lad$node.label<70]<-""
#nodelabels(ml.dense.lad$node.label,adj=c(0.5,0.5),frame="none",bg="white",cex=0.75)
add.scale.bar(x=0.4,y=5,cex=0.5)

#Figure S4
#Detailed ML Ascomycota w boots
ml.dense<-read.tree(file="~/Desktop/deep_lichens/asco_phylo_all_medhts/RAxML_w_500_boots_on_ML_tree/RAxML_bipartitions.bestw500boots")
ml.dense.lad<-ladderize(ml.dense,FALSE)
plot(ml.dense.lad,cex=0.35)
node.support(ml.dense.lad$node.label,cutoff=70,cex=0.5)
#so used this instead, which is ugly
#ml.dense.lad$node.label[ml.dense.lad$node.label<70]<-""
#nodelabels(ml.dense.lad$node.label,adj=c(0.5,0.5),frame="none",bg="white",cex=0.75)
add.scale.bar(x=0.4,y=2,cex=0.75)

#Figure S5
#Sparse Ascomycota w HPD
beast.tree<-read.beast(file="~/Desktop/deep_lichens/asco_phylo_all_medhts/fungi_redo_silog_nonuni_wstart_mcc_medhts.tre")
plot(ladderize(beast.tree,FALSE),cex=0.7,edge.width=0.25,x.lim=c(-100, 800))
HPDbars(ladderize(beast.tree,FALSE),col="darkolivegreen3",lwd=5)
axisPhylo(cex=0.75)

#Figure S6
#SI FIGURE
####Make SI w losses
col2rgb("darkgreen")
#trans.dk.green<-rgb(0,100,0,max=255,alpha=125,names="transdkgreen")

bootswmed<-read.csv(file="~/Desktop/deep_lichens/asco_phylo_all_medhts/individual_bs_trees/corhmm_summary_root_really_fixed.csv",stringsAsFactors=FALSE)
bootswsample<-read.csv(file="~/Desktop/deep_lichens/asco_phylo_ml_100medhts/asco_phylo_boot_25medhts/boots_corhmm_summary_root_really_fixed_combined_correct.csv",stringsAsFactors=FALSE)

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



#Figure S7
require(diversitree)
require(circlize)
require(phytools)
tip.states<-read.csv(file="~/Desktop/deep_lichens/asco_phylo_all_medhts/bad_idea_lichenization_for_corhmm.csv",stringsAsFactors=FALSE)
tip.states<-as.data.frame(tip.states$lichen,row.names=tip.states$edited_name)
colnames(tip.states)<-"lichen"
sum.tab<-read.csv(file="~/Desktop/deep_lichens/asco_phylo_all_medhts/MLTREE_2Rate_rootreallyfixed_corhmm_states.csv",stringsAsFactors=FALSE)
colnames(sum.tab)<-c("Node","0R1","1R1","0R2","1R2","Age")
sum.tab.mat<-as.matrix(sum.tab)
#cols<-c("gray55","green4","gray93","greenyellow","black","blue")
cols<-c("gray75","green4","gray95","darkolivegreen1")
tr<-read.tree(file="~/Desktop/deep_lichens/asco_phylo_all_medhts/MLTREE_2Rate_rootreallyfixed_corhmm_states.tre")
#reduce tip labels to higher-level clades
class<-read.csv(file="~/Desktop/deep_lichens/asco_phylo_all_medhts/bad_idea_lichenization_for_corhmm_w_higher.csv",stringsAsFactors=FALSE)
rownames(class)<-class$edited_name
tiplabs<-tr$tip.label
tax2<-class[match(tiplabs,rownames(class)),]
class2<-class[,"Class"]
mult<-read.csv(file="~/Desktop/deep_lichens/bad_idea_lichenization_for_multicorhmm.csv",stringsAsFactors=FALSE)
mult.tip.states<-as.data.frame(mult$multistate,row.names=mult$edited_name)
tip.states$mult<-NA
for(x in 1:nrow(tip.states)){
	tip.states$mult[x]<-mult$multistate[mult$edited_name %in% rownames(tip.states)[x]]
}
multcols<-c("white","green4","darkorange","cadetblue1")

pdf("~/Desktop/deep_lichens/Fig_S7_mults_21nov2018.pdf",width=11,height=11)

trait.plot(ladderize(tr,FALSE),mult.tip.states,cols=list(multistate=c("grey75","green4","darkorange","cadetblue1")),legend=FALSE,edge.width=0.25,cex.lab=0.75,class=tax2$Class)
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
trait.plot(ladderize(tr,FALSE),mult.tip.states,cols=list(multistate=c("grey75","green4","darkorange","cadetblue1")),legend=FALSE,edge.width=0.25,cex.lab=0.75,class=tax2$Class)
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
trait.plot(ladderize(tr,FALSE),mult.tip.states,cols=list(multistate=c("grey75","green4","darkorange","cadetblue1")),legend=FALSE,edge.width=0.25,cex.lab=0.75,class=tax2$Class)
nodelabels(pie=sum.tab.mat[,2:5],piecol=cols,bg=cols,cex=0.25,border=FALSE)
legend(x=325,y=-490,title=expression(bold("Tip States")),c("Non-Lichenized","Lichenized-Trebouxiophyceae","Lichenized-Trentepohliales","Lichenized-Cyanobacteria"),col=multcols,pch=19,bty="n",cex=1,pt.cex=2)
segments(x0=325,y0=-515,x1=605,y1=-515)
legend(x=325,y=-340,title=expression(bold("        Node States")),c("Non-Lichenized, Stable","Non-Lichenized","Lichenized","Lichenized, Stable"),col=cols[c(1,3,4,2)],pch=19,bty="n",cex=1,pt.cex=2)
segments(x0=325,y0=-370,x1=605,y1=-370)
for(x in 2:10){
	#print(x)
	text(x=ht-rev(timey$Midpoint[x]),y=seq(from=0,to=-7,by=-7/9)[x],labels=timey$Abbrev[x],cex=0.75)
}

dev.off()