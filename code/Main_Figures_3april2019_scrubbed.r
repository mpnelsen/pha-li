###########################################
###########################################
###########################################
#Figure 1 (get ages for lichen fungal algal clades)
#funcion to get age of a clade when providing a single tree
get.age<-function(tree,clade,position){
require(phytools)
node.depth.edgelength(tree)->node.ages
max(node.ages)->depth
depth-node.ages->new.node.ages
if(length(clade)==1){
	tip.number<-which(tree$tip.label==clade)
	parent.node<-Ancestors(tree,tip.number,type="parent")
	age<-new.node.ages[parent.node]
	} else if(length(clade)>1){
		crown.node<-findMRCA(tree,clade)
		if(position=="crown"){
			age<-new.node.ages[crown.node]
			} else if(position=="stem"){
				stem.node<-Ancestors(tree,crown.node,type="parent")
				if(stem.node==0){
					age<-NA
					} else if(stem.node!=0){
						age<-new.node.ages[stem.node]
			}
		}
	}
return(age)
}

#implement above function for multiple trees
get.age.sample<-function(trees,clade,position){
	ages<-NULL
	for(i in 1:length(trees)){
		ages[i]<-get.age(trees[[i]],clade,position)
	}
return(ages)
}

#read algal trees
require(phyloch)
alg<-read.nexus(file="~resamp40k.tre")

#rename clades
aa<-read.csv(file="all4_organized_for_plotting.updated.csv",stringsAsFactors=FALSE,na.strings="")
aa$TaxToUse[aa$Species %in% "Dictyochloropsis_irregularis"]<-"Symbiochloris"
aa$TaxToUse[aa$Species %in% "Dictyochloropsis_sp_MPN124"]<-"Symbiochloris"
aa$TaxToUse[aa$Species %in% "Dictyochloropsis_sp_MPN168"]<-"Symbiochloris"
aa$TaxToUse[aa$Species %in% "Dictyochloropsis_sp_MPN181"]<-"Symbiochloris"
aa$TaxToUse[aa$Species %in% "Dictyochloropsis_splendida"]<-"Symbiochloris"
aa$TaxToUse[aa$Species %in% "Dictyochloropsis_symbiontica_CAUP"]<-"Symbiochloris"
aa$TaxToUse[aa$Species %in% "Dictyochloropsis_symbiontica_SAG"]<-"Symbiochloris"
aa$Genus[aa$Species %in% "Dictyochloropsis_irregularis"]<-"Symbiochloris"
aa$Genus[aa$Species %in% "Dictyochloropsis_sp_MPN124"]<-"Symbiochloris"
aa$Genus[aa$Species %in% "Dictyochloropsis_sp_MPN168"]<-"Symbiochloris"
aa$Genus[aa$Species %in% "Dictyochloropsis_sp_MPN181"]<-"Symbiochloris"
aa$Genus[aa$Species %in% "Dictyochloropsis_splendida"]<-"Symbiochloris"
aa$Genus[aa$Species %in% "Dictyochloropsis_symbiontica_CAUP"]<-"Symbiochloris"
aa$Genus[aa$Species %in% "Dictyochloropsis_symbiontica_SAG"]<-"Symbiochloris"
alg[[1]]$tip.label[!alg[[1]]$tip.label %in% aa$Species]

#get taxa in clades
Trentepohliales<-aa$Species[aa$Order %in% "Trentepohliales"]
names(Trentepohliales)[1:length(Trentepohliales)]<-"Trentepohliales"
Diplosphaera<-aa$Species[aa$Genus %in% "Diplosphaera"]
names(Diplosphaera)[1:length(Diplosphaera)]<-"Diplosphaera"
Symbiochloris<-aa$Species[aa$Genus %in% "Symbiochloris"]
Symbiochloris<-Symbiochloris[!Symbiochloris %in% c("Dictyochloropsis_symbiontica_CAUP","Dictyochloropsis_symbiontica_SAG")]
#Symbiochloris<-c("Dictyochloropsis_irregularis","Dictyochloropsis_sp_MPN124","Dictyochloropsis_sp_MPN168","Dictyochloropsis_sp_MPN181","Dictyochloropsis_splendida")
names(Symbiochloris)[1:length(Symbiochloris)]<-"Symbiochloris"
Asterochloris<-aa$Species[aa$Genus %in% "Asterochloris"]
names(Asterochloris)[1:length(Asterochloris)]<-"Asterochloris"
Asterochloris_Vulcanochloris<-aa$Species[aa$Genus %in% c("Asterochloris","Vulcanochloris")]
names(Asterochloris_Vulcanochloris)[1:length(Asterochloris_Vulcanochloris)]<-"Asterochloris_Vulcanochloris"
Trebouxia<-aa$Species[aa$Genus %in% "Trebouxia"]
names(Trebouxia)[1:length(Trebouxia)]<-"Trebouxia"
Trebouxiales<-aa$Species[aa$Genus %in% c("Trebouxia","Myrmecia","Asterochloris","Vulcanochloris")]
names(Trebouxiales)[1:length(Trebouxiales)]<-"Trebouxiales"

aspecs<-c(Trentepohliales,Diplosphaera,Symbiochloris,Asterochloris,Asterochloris_Vulcanochloris,Trebouxia,Trebouxiales)

#make empty data frame and get ages for multiple clades
algal.clades<-c("Trentepohliales","Diplosphaera","Symbiochloris","Asterochloris","Asterochloris_Vulcanochloris","Trebouxia","Trebouxiales")
alg.df<-as.data.frame(matrix(data=NA,nrow=length(alg),ncol=length(algal.clades)))
colnames(alg.df)<-algal.clades
for(x in algal.clades){
	alg.df[,x]<-get.age.sample(tree=alg,clade=c(aa$Species[aa$Species %in% aspecs[names(aspecs) %in% x]]),position="crown")
}

require(reshape2)
alg.df.lf<-melt(alg.df)

#get fungal trees and clade info
fun<-read.nexus(file="combined_trees.tre")
ff<-read.csv(file="SI_Table_4_broad_asco_7b_4_what_was_in_128test.csv",stringsAsFactors=FALSE,na.strings="")
fun[[1]]$tip.label[!fun[[1]]$tip.label %in% ff$edited_name]

#get taxa in clades
Trypetheliaceae<-ff$edited_name[ff$Family %in% "Trypetheliaceae"]
names(Trypetheliaceae)[1:length(Trypetheliaceae)]<-"Trypetheliaceae"
Arthoniales<-ff$edited_name[ff$Order %in% "Arthoniales"]
names(Arthoniales)[1:length(Arthoniales)]<-"Arthoniales"
Arthoniales<-Arthoniales[!Arthoniales %in% c("Chrysothrix_caesia","Chrysothrix_chrysophthalma")]
Ostropales<-ff$edited_name[ff$Order %in% "Ostropales"]
names(Ostropales)[1:length(Ostropales)]<-"Ostropales"
Ostropales<-Ostropales[!Ostropales %in% c("Odontotrema_phacidioides","Stictis_radiata")]
Pyrenulaceae<-ff$edited_name[ff$Family %in% "Pyrenulaceae"]
names(Pyrenulaceae)[1:length(Pyrenulaceae)]<-"Pyrenulaceae"
Strigulaceae<-ff$edited_name[ff$Family %in% "Strigulaceae"]
names(Strigulaceae)[1:length(Strigulaceae)]<-"Strigulaceae"
Monoblastiaceae<-ff$edited_name[ff$Family %in% "Monoblastiaceae"]
names(Monoblastiaceae)[1:length(Monoblastiaceae)]<-"Monoblastiaceae"
Monoblastiaceae<-Monoblastiaceae[!Monoblastiaceae %in% "Heleiosa_barbatula"]
Verrucariaceae<-ff$edited_name[ff$Family %in% "Verrucariaceae"]
names(Verrucariaceae)[1:length(Verrucariaceae)]<-"Verrucariaceae"
Phlyctidaceae<-ff$edited_name[ff$Family %in% "Phlyctidaceae"]
names(Phlyctidaceae)[1:length(Phlyctidaceae)]<-"Phlyctidaceae"
Teloschistales<-ff$edited_name[ff$Order %in% "Teloschistales"]
names(Teloschistales)[1:length(Teloschistales)]<-"Teloschistales"
Cladoniineae<-ff$edited_name[ff$Family %in% c("Cladoniaceae","Stereocaulaceae","Squamarinaceae")]
names(Cladoniineae)[1:length(Cladoniineae)]<-"Cladoniineae"
Caliciales<-ff$edited_name[ff$Order %in% "Caliciales"]
names(Caliciales)[1:length(Caliciales)]<-"Caliciales"
Chrysothrichaceae<-ff$edited_name[ff$Family %in% "Chrysothrichaceae"]
names(Chrysothrichaceae)[1:length(Chrysothrichaceae)]<-"Chrysothrichaceae"
Teloschistaceae<-ff$edited_name[ff$Family %in% "Teloschistaceae"]
names(Teloschistaceae)[1:length(Teloschistaceae)]<-"Teloschistaceae"
Lecanoromycetes<-ff$edited_name[ff$Class %in% "Lecanoromycetes"]
names(Lecanoromycetes)[1:length(Lecanoromycetes)]<-"Lecanoromycetes"

fspecs<-c(Trypetheliaceae,Arthoniales,Ostropales,Pyrenulaceae,Strigulaceae,Monoblastiaceae,Verrucariaceae,Phlyctidaceae,Teloschistales,Cladoniineae,Caliciales,Chrysothrichaceae,Teloschistaceae,Lecanoromycetes)

#make empty data frame and get ages for clades
fungal.clades<-c("Trypetheliaceae","Arthoniales","Ostropales","Pyrenulaceae","Strigulaceae","Monoblastiaceae","Verrucariaceae","Phlyctidaceae","Teloschistales","Cladoniineae","Caliciales","Chrysothrichaceae","Teloschistaceae","Lecanoromycetes")
fun.df<-as.data.frame(matrix(data=NA,nrow=length(fun),ncol=length(fungal.clades)))
colnames(fun.df)<-fungal.clades
for(x in fungal.clades){
	fun.df[,x]<-get.age.sample(tree=fun,clade=c(ff$edited_name[ff$edited_name %in% fspecs[names(fspecs) %in% x]]),position="crown")
}

require(reshape2)
fun.df.lf<-melt(fun.df)

#order them
xvf<-rbind(alg.df.lf,fun.df.lf)
colnames(xvf)<-c("Clade","Age")
xvf$NodeNo<-NA
xvf$NodeNo[xvf$Clade %in% "Trentepohliales"]<-1
xvf$NodeNo[xvf$Clade %in% "Ostropales"]<-2
xvf$NodeNo[xvf$Clade %in% "Arthoniales"]<-3
xvf$NodeNo[xvf$Clade %in% "Trypetheliaceae"]<-4
xvf$NodeNo[xvf$Clade %in% "Pyrenulaceae"]<-5
xvf$NodeNo[xvf$Clade %in% "Strigulaceae"]<-6
xvf$NodeNo[xvf$Clade %in% "Monoblastiaceae"]<-7
xvf$NodeNo[xvf$Clade %in% "Trebouxia"]<-8
xvf$NodeNo[xvf$Clade %in% "Caliciales"]<-9
xvf$NodeNo[xvf$Clade %in% "Chrysothrichaceae"]<-10
xvf$NodeNo[xvf$Clade %in% "Teloschistaceae"]<-11
xvf$NodeNo[xvf$Clade %in% "Diplosphaera"]<-12
xvf$NodeNo[xvf$Clade %in% "Verrucariaceae"]<-13
xvf$NodeNo[xvf$Clade %in% "Symbiochloris"]<-14
xvf$NodeNo[xvf$Clade %in% "Teloschistales"]<-15
xvf$NodeNo[xvf$Clade %in% "Phlyctidaceae"]<-16
xvf$NodeNo[xvf$Clade %in% "Asterochloris"]<-17
xvf$NodeNo[xvf$Clade %in% "Asterochloris_Vulcanochloris"]<-18
xvf$NodeNo[xvf$Clade %in% "Cladoniineae"]<-19
xvf$NodeNo[xvf$Clade %in% "Trebouxiales"]<-20
xvf$NodeNo[xvf$Clade %in% "Lecanoromycetes"]<-21

#assign color scheme
xvf["Color"]<-NA
xvf$Color[xvf$Clade %in% algal.clades]<-"Algae"
xvf$Color[xvf$Clade %in% fungal.clades]<-"Fungi"

#add min/max hpd
min.hpd<-c(258.7,123.2,102.7,99.1,75.3,59.9,41.3,65.5,50.3,43.6,35.4,64.9,50.8,74.9,76.9,40.4,13,34.9,36.8,173.8,199.7)
max.hpd<-c(489.7,206.4,196.4,183.2,160.4,139.2,107.4,188.3,167.1,163.7,92.8,191.7,113.1,177.9,151.8,138.0,51,121.1,85.8,403.3,303.0)
ordered.clades<-c("Trentepohliales","Ostropales","Arthoniales","Trypetheliaceae","Pyrenulaceae","Strigulaceae","Monoblastiaceae","Trebouxia","Caliciales","Chrysothrichaceae","Teloschistaceae","Diplosphaera","Verrucariaceae","Symbiochloris","Teloschistales","Phlyctidaceae","Asterochloris","Asterochloris_Vulcanochloris","Cladoniineae","Trebouxiales","Lecanoromycetes")
xvf$minHPD<-NA
xvf$maxHPD<-NA
for(x in 1:length(ordered.clades)){
	xvf$minHPD[xvf$Clade==ordered.clades[x]]<-min.hpd[x]
	xvf$maxHPD[xvf$Clade==ordered.clades[x]]<-max.hpd[x]
}

#organize
combo<-xvf[order(xvf$NodeNo),]

#rename Asterochloris_Vulcanochloris
levels(combo$Clade)[levels(combo$Clade)=="Asterochloris_Vulcanochloris"]<-"Asterochloris (Stem)"

#wrap clade labels
require(stringr)
combo$Clade<-str_wrap(combo$Clade,width=20)

#make factors
combo$Node2<-factor(combo$Clade,levels=unique(combo$Clade))
combo$Age<-as.numeric(combo$Age)

#multiply by -1
combo[,c("Age","minHPD","maxHPD")]<-combo[,c("Age","minHPD","maxHPD")]*(-1)

#assign plain or italics
face.plot<-c("plain","plain","plain","plain","plain","plain","plain","italic","plain","plain","plain","italic","plain","italic","plain","plain","italic","italic","plain","plain","plain")

#get timescale
timescale<-read.csv(file="timescale_ics2015_modified.csv",stringsAsFactors=FALSE)
for(x in 1:nrow(timescale)){
	timescale$RGB[x]<-rgb(timescale$Col_R[x]/255,timescale$Col_G[x]/255,timescale$Col_B[x]/255)
}
timey<-timescale[timescale$Type %in% "Period",]
timey<-timey[timey$End<675,]
timey$Name<-c(NA,"Ng","Pg","K","J","Tr","P","C","D","S","O","Cm","Ed","Cr")
timey$Start[nrow(timey)]<-675

#multiply by -1
timey[,c("Start","Midpoint","End")]<-timey[,c("Start","Midpoint","End")]*(-1)

require(ggplot2)

#create plot
p<-ggplot(combo,aes(x=Node2,y=Age,ymin=minHPD,ymax=maxHPD))+scale_y_continuous(limits=c(-675,0),expand=c(0,0))+geom_violin(aes(x=Node2,fill=Color,color=Color))+scale_color_manual(values=c("white","white"))+scale_fill_manual(values=c("white","white"))+annotate("rect", ymin=timey$End[1],ymax=timey$Start[1],xmin=-Inf,xmax=Inf,fill=timey$RGB[1],color=NA,alpha=0.5)+annotate("rect", ymin=timey$End[2],ymax=timey$Start[2],xmin=-Inf,xmax=Inf,fill=timey$RGB[2],color=NA,alpha=0.5)+annotate("rect", ymin=timey$End[3],ymax=timey$Start[3],xmin=-Inf,xmax=Inf,fill=timey$RGB[3],color=NA,alpha=0.5)+annotate("rect", ymin=timey$End[4],ymax=timey$Start[4],xmin=-Inf,xmax=Inf,fill=timey$RGB[4],color=NA,alpha=0.5)+annotate("rect", ymin=timey$End[5],ymax=timey$Start[5],xmin=-Inf,xmax=Inf,fill=timey$RGB[5],color=NA,alpha=0.5)+annotate("rect", ymin=timey$End[6],ymax=timey$Start[6],xmin=-Inf,xmax=Inf,fill=timey$RGB[6],color=NA,alpha=0.5)+annotate("rect", ymin=timey$End[7],ymax=timey$Start[7],xmin=-Inf,xmax=Inf,fill=timey$RGB[7],color=NA,alpha=0.5)+annotate("rect",ymin=timey$End[8],ymax=timey$Start[8],xmin=-Inf,xmax=Inf,fill=timey$RGB[8],color=NA,alpha=0.5)+annotate("rect",ymin=timey$End[9],ymax=timey$Start[9],xmin=-Inf,xmax=Inf,fill=timey$RGB[9],color=NA,alpha=0.5)+annotate("rect",ymin=timey$End[10],ymax=timey$Start[10],xmin=-Inf,xmax=Inf,fill=timey$RGB[10],color=NA,alpha=0.5)+annotate("rect",ymin=timey$End[11],ymax=timey$Start[11],xmin=-Inf,xmax=Inf,fill=timey$RGB[11],color=NA,alpha=0.5)+annotate("rect",ymin=timey$End[12],ymax=timey$Start[12],xmin=-Inf,xmax=Inf,fill=timey$RGB[12],color=NA,alpha=0.5)+annotate("rect",ymin=timey$End[13],ymax=timey$Start[13],xmin=-Inf,xmax=Inf,fill=timey$RGB[13],color=NA,alpha=0.5)+annotate("rect",ymin=timey$End[14],ymax=timey$Start[14],xmin=-Inf,xmax=Inf,fill=timey$RGB[14],color=NA,alpha=0.5)+geom_violin(aes(x=Node2,fill=Color,color=Color))+scale_color_manual(values=c("palegreen4","navajowhite4"))+scale_fill_manual(values=c("palegreen4","navajowhite4"))+geom_errorbar(aes(x=Node2),size=0.05,width=0.25,color="white")+stat_summary(aes(group=Node2),fun.y=median,geom="point",size=1,color="white")+theme(axis.text.y=element_text(size=12),axis.text.x  = element_text(angle=90, vjust=0.5, hjust=1,size=12,face=face.plot),plot.title=element_text(size=16,hjust=0.5),axis.line = element_line(color = "black"),panel.grid.major = 'element_blank'(),panel.grid.minor = 'element_blank'(),panel.border = 'element_blank'(),panel.background = 'element_blank'(),legend.position="right", legend.title='element_blank'(),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14),legend.text=element_text(size=14))+ggtitle("Ages of Interacting Clades")+ylab("Age (Ma)")+xlab("Clade")+geom_vline(xintercept=c(7.5,11.5,13.5,16.5,19.5),linetype="dotted",color="grey50")+geom_vline(xintercept=0.5,color="black")+annotate("text",y=-654,x=0,label=timey$Name[14],size=3,colour="grey35")+annotate("text",y=timey$Midpoint[13],x=0,label=timey$Name[13],size=3,colour="grey35")+annotate("text",y=timey$Midpoint[12],x=0,label=timey$Name[12],size=3,colour="grey35")+annotate("text",y=timey$Midpoint[11],x=0,label=timey$Name[11],size=3,colour="grey35")+annotate("text",y=timey$Midpoint[10],x=0,label=timey$Name[10],size=3,colour="grey35")+annotate("text",y=timey$Midpoint[9],x=0,label=timey$Name[9],size=3,colour="grey35")+annotate("text",y=timey$Midpoint[8],x=0,label=timey$Name[8],size=3,colour="grey35")+annotate("text",y=timey$Midpoint[7],x=0,label=timey$Name[7],size=3,colour="grey35")+annotate("text",y=timey$Midpoint[6],x=0,label=timey$Name[6],size=3,colour="grey35")+annotate("text",y=timey$Midpoint[5],x=0,label=timey$Name[5],size=3,colour="grey35")+annotate("text",y=timey$Midpoint[4],x=0,label=timey$Name[4],size=3,colour="grey35")+annotate("text",y=timey$Midpoint[3],x=0,label=timey$Name[3],size=3,colour="grey35")+annotate("text",y=timey$Midpoint[2],x=0,label=timey$Name[2],size=3,colour="grey35")+expand_limits(x=-0.5)

#plot
p

###########################################
###########################################
###########################################
#Figure 2
require(diversitree)
require(circlize)
require(phytools)
tip.states<-read.csv(file="bad_idea_lichenization_for_corhmm.csv",stringsAsFactors=FALSE)
tip.states<-as.data.frame(tip.states$lichen,row.names=tip.states$edited_name)
colnames(tip.states)<-"lichen"
sum.tab<-read.csv(file="MLTREE_2Rate_rootreallyfixed_corhmm_states.csv",stringsAsFactors=FALSE)
colnames(sum.tab)<-c("Node","0R1","1R1","0R2","1R2","Age")
sum.tab.mat<-as.matrix(sum.tab)
#cols<-c("gray55","green4","gray93","greenyellow","black","blue")
cols<-c("gray75","green4","gray95","darkolivegreen1")
tr<-read.tree(file="MLTREE_2Rate_rootreallyfixed_corhmm_states.tre")
#reduce tip labels to higher-level clades
class<-read.csv(file="bad_idea_lichenization_for_corhmm_w_higher.csv",stringsAsFactors=FALSE)
rownames(class)<-class$edited_name
tiplabs<-tr$tip.label
tax2<-class[match(tiplabs,rownames(class)),]
class2<-class[,"Class"]
mult<-read.csv(file="bad_idea_lichenization_for_multicorhmm.csv",stringsAsFactors=FALSE)
mult.tip.states<-as.data.frame(mult$multistate,row.names=mult$edited_name)
tip.states$mult<-NA
for(x in 1:nrow(tip.states)){
	tip.states$mult[x]<-mult$multistate[mult$edited_name %in% rownames(tip.states)[x]]
}
multcols<-c("grey75","green4","darkorange","cadetblue1")

pdf("Fig_2_mults_3april2019.pdf",width=11,height=11)

trait.plot(ladderize(tr,FALSE),mult.tip.states,cols=list(multistate=c("grey75","green4","darkorange","cadetblue1")),legend=FALSE,edge.width=0.25,cex.lab=0.8,class=tax2$Class)
#nodelabels(frame="circle",col=cols[as.numeric(tr$node.label)],bg=cols[as.numeric(tr$node.label)],cex=0.15)
nodelabels(pie=sum.tab.mat[,2:5],piecol=cols,bg=cols,cex=0.25)
legend("topleft",title="Character State",c("Non-Lichenized, Stable","Non-Lichenized","Lichenized","Lichenized, Stable"),fill=cols[c(1,3,4,2)])

timescale<-read.csv(file="timescale_ics2015_modified.csv",stringsAsFactors=FALSE)

for(x in 1:nrow(timescale)){
	timescale$RGB[x]<-rgb(timescale$Col_R[x]/255,timescale$Col_G[x]/255,timescale$Col_B[x]/255,alpha=0.2)
}
require(phytools)
timey<-timescale[timescale$Type %in% "Period",]
timey<-timey[timey$End<max(nodeHeights(tr)),]

#add age rings
trait.plot(ladderize(tr,FALSE),mult.tip.states,cols=list(multistate=c("grey75","green4","darkorange","cadetblue1")),legend=FALSE,edge.width=0.25,cex.lab=0.8,class=tax2$Class)
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
trait.plot(ladderize(tr,FALSE),mult.tip.states,cols=list(multistate=c("grey75","green4","darkorange","cadetblue1")),legend=FALSE,edge.width=0.25,cex.lab=0.8,class=tax2$Class)
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


###########################################
###########################################
###########################################
#Figure 3 (deep_fungal_algal_comparison)
#funcion to get age of a clade when providing a single tree
get.age<-function(tree,clade,position){
require(phytools)
node.depth.edgelength(tree)->node.ages
max(node.ages)->depth
depth-node.ages->new.node.ages
if(length(clade)==1){
	tip.number<-which(tree$tip.label==clade)
	parent.node<-Ancestors(tree,tip.number,type="parent")
	age<-new.node.ages[parent.node]
	} else if(length(clade)>1){
		crown.node<-findMRCA(tree,clade)
		if(position=="crown"){
			age<-new.node.ages[crown.node]
			} else if(position=="stem"){
				stem.node<-Ancestors(tree,crown.node,type="parent")
				if(stem.node==0){
					age<-NA
					} else if(stem.node!=0){
						age<-new.node.ages[stem.node]
			}
		}
	}
return(age)
}

#implement above function for multiple trees
get.age.sample<-function(trees,clade,position){
	ages<-NULL
	for(i in 1:length(trees)){
		ages[i]<-get.age(trees[[i]],clade,position)
	}
return(ages)
}

#read algal trees
require(phyloch)
alg<-read.nexus(file="resamp40k.tre")

#rename
aa<-read.csv(file="all4_organized_for_plotting.updated.csv",stringsAsFactors=FALSE,na.strings="")
aa$TaxToUse[aa$Species %in% "Dictyochloropsis_irregularis"]<-"Symbiochloris"
aa$TaxToUse[aa$Species %in% "Dictyochloropsis_sp_MPN124"]<-"Symbiochloris"
aa$TaxToUse[aa$Species %in% "Dictyochloropsis_sp_MPN168"]<-"Symbiochloris"
aa$TaxToUse[aa$Species %in% "Dictyochloropsis_sp_MPN181"]<-"Symbiochloris"
aa$TaxToUse[aa$Species %in% "Dictyochloropsis_splendida"]<-"Symbiochloris"
aa$TaxToUse[aa$Species %in% "Dictyochloropsis_symbiontica_CAUP"]<-"Symbiochloris"
aa$TaxToUse[aa$Species %in% "Dictyochloropsis_symbiontica_SAG"]<-"Symbiochloris"
aa$Genus[aa$Species %in% "Dictyochloropsis_irregularis"]<-"Symbiochloris"
aa$Genus[aa$Species %in% "Dictyochloropsis_sp_MPN124"]<-"Symbiochloris"
aa$Genus[aa$Species %in% "Dictyochloropsis_sp_MPN168"]<-"Symbiochloris"
aa$Genus[aa$Species %in% "Dictyochloropsis_sp_MPN181"]<-"Symbiochloris"
aa$Genus[aa$Species %in% "Dictyochloropsis_splendida"]<-"Symbiochloris"
aa$Genus[aa$Species %in% "Dictyochloropsis_symbiontica_CAUP"]<-"Symbiochloris"
aa$Genus[aa$Species %in% "Dictyochloropsis_symbiontica_SAG"]<-"Symbiochloris"
alg[[1]]$tip.label[!alg[[1]]$tip.label %in% aa$Species]

#get taxa in clades
Trentepohliales<-aa$Species[aa$Order %in% "Trentepohliales"]
names(Trentepohliales)[1:length(Trentepohliales)]<-"Trentepohliales"
Diplosphaera_Prasiola<-aa$Species[aa$Genus %in% c("Diplosphaera","Stichococcus")]
names(Diplosphaera_Prasiola)[1:length(Diplosphaera_Prasiola)]<-"Diplosphaera_Prasiola"
Symbiochloris<-aa$Species[aa$Genus %in% "Symbiochloris"]
names(Symbiochloris)[1:length(Symbiochloris)]<-"Symbiochloris"
Trebouxiales<-aa$Species[aa$Genus %in% c("Trebouxia","Myrmecia","Asterochloris","Vulcanochloris")]
names(Trebouxiales)[1:length(Trebouxiales)]<-"Trebouxiales"
Heveochlorella_Heterochlorella<-aa$Species[aa$Genus %in% c("Heterochlorella","Heveochlorella")]
names(Heveochlorella_Heterochlorella)[1:length(Heveochlorella_Heterochlorella)]<-"Heveochlorella_Heterochlorella"
Coccomyxa_Elliptochloris<-aa$Species[aa$Genus %in% c("Coccomyxa","Elliptochloris")]
names(Coccomyxa_Elliptochloris)[1:length(Coccomyxa_Elliptochloris)]<-"Coccomyxa_Elliptochloris"

aspecs<-c(Trentepohliales,Coccomyxa_Elliptochloris,Trebouxiales,Symbiochloris,Diplosphaera_Prasiola,Heveochlorella_Heterochlorella)

#make empty data frame and get ages of clades
algal.clades<-c("Trentepohliales","Coccomyxa_Elliptochloris","Trebouxiales","Symbiochloris","Diplosphaera_Prasiola","Heveochlorella_Heterochlorella")
alg.df<-as.data.frame(matrix(data=NA,nrow=length(alg),ncol=length(algal.clades)))
colnames(alg.df)<-algal.clades
for(x in algal.clades){
	alg.df[,x]<-get.age.sample(tree=alg,clade=c(aa$Species[aa$Species %in% aspecs[names(aspecs) %in% x]]),position="crown")
}

require(reshape2)
alg.df.lf<-melt(alg.df)

#read summaries of gains of lichenization in fungi
bootswmed<-read.csv(file="corhmm_summary_root_really_fixed.csv",stringsAsFactors=FALSE)
bootswsample<-read.csv(file="boots_corhmm_summary_root_really_fixed_combined_correct.csv",stringsAsFactors=FALSE)
mlwsample<-read.csv(file="corhmm_summary_root_really_fixed.csv")
mlwmed<-248.6
fun.df.lf<-as.data.frame(matrix(data=NA,nrow=(nrow(bootswmed)+nrow(bootswsample)+nrow(mlwsample)+1),ncol=2))
colnames(fun.df.lf)<-c("variable","value")
fun.df.lf$value[1:nrow(bootswmed)]<-bootswmed$MaxAgeEndLich
fun.df.lf$variable[1:nrow(bootswmed)]<-"LFF_BM"
fun.df.lf$value[(nrow(bootswmed)+1):(nrow(bootswmed)+nrow(bootswsample))]<-bootswsample$MaxAgeEndLich
fun.df.lf$variable[(nrow(bootswmed)+1):(nrow(bootswmed)+nrow(bootswsample))]<-"LFF_BS"
fun.df.lf$value[(nrow(bootswmed)+nrow(bootswsample)+1):(nrow(bootswmed)+nrow(bootswsample)+nrow(mlwsample))]<-mlwsample$MaxAgeEndLich
fun.df.lf$variable[(nrow(bootswmed)+nrow(bootswsample)+1):(nrow(bootswmed)+nrow(bootswsample)+nrow(mlwsample))]<-"LFF_MS"
fun.df.lf$value[(nrow(bootswmed)+nrow(bootswsample)+nrow(mlwsample)+1)]<-mlwmed
fun.df.lf$variable[(nrow(bootswmed)+nrow(bootswsample)+nrow(mlwsample)+1)]<-"LFF_MM"

fun.df.lf$variable<-as.factor(fun.df.lf$variable)

fungal.clades<-c("LFF_BM","LFF_BS","LFF_MS","LFF_MM")

#order them
xvf<-rbind(alg.df.lf,fun.df.lf)
colnames(xvf)<-c("Clade","Age")
xvf$NodeNo<-NA
xvf$NodeNo[xvf$Clade %in% "Trentepohliales"]<-5
xvf$NodeNo[xvf$Clade %in% "Coccomyxa_Elliptochloris"]<-6
xvf$NodeNo[xvf$Clade %in% "Trebouxiales"]<-7
xvf$NodeNo[xvf$Clade %in% "Symbiochloris"]<-8
xvf$NodeNo[xvf$Clade %in% "Diplosphaera_Prasiola"]<-9
xvf$NodeNo[xvf$Clade %in% "Heveochlorella_Heterochlorella"]<-10
xvf$NodeNo[xvf$Clade %in% "LFF_MM"]<-1
xvf$NodeNo[xvf$Clade %in% "LFF_BM"]<-2
xvf$NodeNo[xvf$Clade %in% "LFF_MS"]<-3
xvf$NodeNo[xvf$Clade %in% "LFF_BS"]<-4

#assign color scheme
xvf["Color"]<-NA
xvf$Color[xvf$Clade %in% algal.clades]<-"Algae"
xvf$Color[xvf$Clade %in% fungal.clades]<-"Fungi"

#add min/max hpd...also add median node age estimates for algae, and ml estimate for LFF w median, and mean estimate for LFF w sample
min.hpd<-c(248.6,248.6,210.5,232.6,258.7,196.2,173.8,165.9,91.8,75.2)
max.hpd<-c(248.6,350.2,332.7,396.9,489.7,447.5,403.3,342.1,267.2,220.3)
median.est<-c(248.6,289.7,257.4,291.5,368.2,316.9,279.8,247.6,169.7,140.4)
ordered.clades<-c("LFF_MM","LFF_BM","LFF_MS","LFF_BS","Trentepohliales","Coccomyxa_Elliptochloris","Trebouxiales","Symbiochloris","Diplosphaera_Prasiola","Heveochlorella_Heterochlorella")
xvf$minHPD<-NA
xvf$maxHPD<-NA
xvf$median_est<-NA
for(x in 1:length(ordered.clades)){
	xvf$minHPD[xvf$Clade==ordered.clades[x]]<-min.hpd[x]
	xvf$maxHPD[xvf$Clade==ordered.clades[x]]<-max.hpd[x]
	xvf$median_est[xvf$Clade==ordered.clades[x]]<-median.est[x]
}


#organize
combo<-xvf[order(xvf$NodeNo),]

#rename some clades
levels(combo$Clade)[levels(combo$Clade)=="Coccomyxa_Elliptochloris"]<-"Coccomyxa + Elliptochloris"
levels(combo$Clade)[levels(combo$Clade)=="Diplosphaera_Prasiola"]<-"Diplosphaera + Prasiola"
levels(combo$Clade)[levels(combo$Clade)=="Heveochlorella_Heterochlorella"]<-"Heveochlorella + Heterochlorella"
levels(combo$Clade)[levels(combo$Clade)=="LFF_MM"]<-"LFF (ML-M)"
levels(combo$Clade)[levels(combo$Clade)=="LFF_BM"]<-"LFF (BS-M)"
levels(combo$Clade)[levels(combo$Clade)=="LFF_MS"]<-"LFF (ML-S)"
levels(combo$Clade)[levels(combo$Clade)=="LFF_BS"]<-"LFF (BS-S)"

#wrap clade labels
require(stringr)
combo$Clade<-str_wrap(combo$Clade,width=16)

#make factors
combo$Node2<-factor(combo$Clade,levels=unique(combo$Clade))
combo$Age<-as.numeric(combo$Age)

#multiply by -1
combo[,c("Age","minHPD","maxHPD","median_est")]<-combo[,c("Age","minHPD","maxHPD","median_est")]*(-1)

#assign plain or italics
face.plot<-c("plain","plain","plain","plain","plain","italic","plain","italic","italic","italic")

#get timescale
timescale<-read.csv(file="timescale_ics2015_modified.csv",stringsAsFactors=FALSE)
for(x in 1:nrow(timescale)){
	timescale$RGB[x]<-rgb(timescale$Col_R[x]/255,timescale$Col_G[x]/255,timescale$Col_B[x]/255)
}
timey<-timescale[timescale$Type %in% "Period",]
timey<-timey[timey$End<675,]
timey$Name<-c(NA,"Ng","Pg","K","J","Tr","P","C","D","S","O","Cm","Ed","Cr")
timey$Start[nrow(timey)]<-675

#multiply by -1
timey[,c("Start","Midpoint","End")]<-timey[,c("Start","Midpoint","End")]*(-1)

require(ggplot2)

#create plot
p<-ggplot(combo,aes(x=Node2,y=Age,ymin=minHPD,ymax=maxHPD))+scale_y_continuous(limits=c(-675,0),expand=c(0,0))+geom_violin(aes(x=Node2,fill=Color,color=Color))+scale_color_manual(values=c("white","white"))+scale_fill_manual(values=c("white","white"))+annotate("rect", ymin=timey$End[1],ymax=timey$Start[1],xmin=-Inf,xmax=Inf,fill=timey$RGB[1],color=NA,alpha=0.5)+annotate("rect", ymin=timey$End[2],ymax=timey$Start[2],xmin=-Inf,xmax=Inf,fill=timey$RGB[2],color=NA,alpha=0.5)+annotate("rect", ymin=timey$End[3],ymax=timey$Start[3],xmin=-Inf,xmax=Inf,fill=timey$RGB[3],color=NA,alpha=0.5)+annotate("rect", ymin=timey$End[4],ymax=timey$Start[4],xmin=-Inf,xmax=Inf,fill=timey$RGB[4],color=NA,alpha=0.5)+annotate("rect", ymin=timey$End[5],ymax=timey$Start[5],xmin=-Inf,xmax=Inf,fill=timey$RGB[5],color=NA,alpha=0.5)+annotate("rect", ymin=timey$End[6],ymax=timey$Start[6],xmin=-Inf,xmax=Inf,fill=timey$RGB[6],color=NA,alpha=0.5)+annotate("rect", ymin=timey$End[7],ymax=timey$Start[7],xmin=-Inf,xmax=Inf,fill=timey$RGB[7],color=NA,alpha=0.5)+annotate("rect",ymin=timey$End[8],ymax=timey$Start[8],xmin=-Inf,xmax=Inf,fill=timey$RGB[8],color=NA,alpha=0.5)+annotate("rect",ymin=timey$End[9],ymax=timey$Start[9],xmin=-Inf,xmax=Inf,fill=timey$RGB[9],color=NA,alpha=0.5)+annotate("rect",ymin=timey$End[10],ymax=timey$Start[10],xmin=-Inf,xmax=Inf,fill=timey$RGB[10],color=NA,alpha=0.5)+annotate("rect",ymin=timey$End[11],ymax=timey$Start[11],xmin=-Inf,xmax=Inf,fill=timey$RGB[11],color=NA,alpha=0.5)+annotate("rect",ymin=timey$End[12],ymax=timey$Start[12],xmin=-Inf,xmax=Inf,fill=timey$RGB[12],color=NA,alpha=0.5)+annotate("rect",ymin=timey$End[13],ymax=timey$Start[13],xmin=-Inf,xmax=Inf,fill=timey$RGB[13],color=NA,alpha=0.5)+annotate("rect",ymin=timey$End[14],ymax=timey$Start[14],xmin=-Inf,xmax=Inf,fill=timey$RGB[14],color=NA,alpha=0.5)+geom_violin(aes(x=Node2,fill=Color,color=Color))+scale_color_manual(values=c("palegreen4","navajowhite4"))+scale_fill_manual(values=c("palegreen4","navajowhite4"))+geom_errorbar(aes(x=Node2),size=0.05,width=0.25,color="white")+geom_point(aes(y=median_est),size=1,color="white")+theme(axis.text.y=element_text(size=12),axis.text.x  = element_text(angle=90, vjust=0.5, hjust=1,size=12,face=face.plot),plot.title=element_text(size=16,hjust=0.5),axis.line = element_line(color = "black"),panel.grid.major = 'element_blank'(),panel.grid.minor = 'element_blank'(),panel.border = 'element_blank'(),panel.background = 'element_blank'(),legend.position="right", legend.title='element_blank'(),axis.title.x=element_text(size=14),axis.title.y=element_text(size=14),legend.text=element_text(size=14))+ggtitle("Ages of Lichen-Forming Clades")+ylab("Age (Ma)")+xlab("Clade")+geom_vline(xintercept=0.5,color="black")+geom_hline(yintercept=-425,linetype="dotted",color="grey35")+annotate("text",y=-654,x=0,label=timey$Name[14],size=3,colour="grey35")+annotate("text",y=timey$Midpoint[13],x=0,label=timey$Name[13],size=3,colour="grey35")+annotate("text",y=timey$Midpoint[12],x=0,label=timey$Name[12],size=3,colour="grey35")+annotate("text",y=timey$Midpoint[11],x=0,label=timey$Name[11],size=3,colour="grey35")+annotate("text",y=timey$Midpoint[10],x=0,label=timey$Name[10],size=3,colour="grey35")+annotate("text",y=timey$Midpoint[9],x=0,label=timey$Name[9],size=3,colour="grey35")+annotate("text",y=timey$Midpoint[8],x=0,label=timey$Name[8],size=3,colour="grey35")+annotate("text",y=timey$Midpoint[7],x=0,label=timey$Name[7],size=3,colour="grey35")+annotate("text",y=timey$Midpoint[6],x=0,label=timey$Name[6],size=3,colour="grey35")+annotate("text",y=timey$Midpoint[5],x=0,label=timey$Name[5],size=3,colour="grey35")+annotate("text",y=timey$Midpoint[4],x=0,label=timey$Name[4],size=3,colour="grey35")+annotate("text",y=timey$Midpoint[3],x=0,label=timey$Name[3],size=3,colour="grey35")+annotate("text",y=timey$Midpoint[2],x=0,label=timey$Name[2],size=3,colour="grey35")+expand_limits(x=-0.5)+annotate("rect",xmin=1.55,xmax=3.45,ymin=-522,ymax=-484,fill="white",color="white",linetype="solid")+annotate("text",x=2.5,y=-500,label="Tracheophytes",color="grey35",size=3.5)+geom_segment(x=2.5,xend=2.5,y=-482,yend=-425,color="grey35",arrow=arrow(length=unit(0.2,"cm"),type="closed"))

#plot
p