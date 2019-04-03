require(phyloch)
alg<-read.nexus(file="resamp40k.tre")
aa<-read.csv(file="all4_organized_for_plotting.updated.csv",stringsAsFactors=FALSE,na.strings="")

ml.age<-get.age.sample(tree=alg,clade=c(aa$Species[aa$TaxToUse %in% "Trentepohliales"]),position="crown")

#95% HPD for Trentepohliales Crown is 
hpd.min<-258.6591681342161
hpd.max<-489.7007404370182

#retain ages in 95% HPD
hpd.age<-ml.age[hpd.max>=ml.age & ml.age>=hpd.min]

length(hpd.age[hpd.age>=425])/length(hpd.age)
length(hpd.age[hpd.age<425])/length(hpd.age)
length(ml.age[ml.age>=425])/length(ml.age)
length(ml.age[ml.age<425])/length(ml.age)



require(phyloch)
alg<-read.nexus(file="resamp40k.tre")
aa<-read.csv(file="all4_organized_for_plotting.updated.csv",stringsAsFactors=FALSE,na.strings="")
ml.age<-get.age.sample(tree=alg,clade=c(aa$Species[aa$TaxToUse %in% c("Coccomyxa","Elliptochloris")]),position="crown")

#95% HPD for Coccomyxa+Elliptochloris Crown is 
hpd.min<-196.2403
hpd.max<-447.5478

#retain ages in 95% HPD
hpd.age<-ml.age[hpd.max>=ml.age & ml.age>=hpd.min]

length(hpd.age[hpd.age>=425])/length(hpd.age)
length(hpd.age[hpd.age<425])/length(hpd.age)
length(ml.age[ml.age>=425])/length(ml.age)
length(ml.age[ml.age<425])/length(ml.age)

#when a single tree
get.age<-function(tree,clade,position){
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

#when multiple trees
get.age.sample<-function(trees,clade,position){
	ages<-NULL
	for(i in 1:length(trees)){
		ages[i]<-get.age(trees[[i]],clade,position)
	}
return(ages)
}
