
##################################################################################################################################
#############Record changes in ml tree timescaled with median node ages (mltree_timescaled_w_median_node_heights)
##################################################################################################################################
require(ape)
require(phytools)
require(phangorn)

tree.file.path<-"~/Desktop/asco_phylo_all_medhts/"
tip.states<-read.csv(file="~/Desktop/asco_phylo_all_medhts/bad_idea_lichenization_for_corhmm.csv",stringsAsFactors=FALSE)

#read csv file
mat<-read.csv(file="~/Desktop/asco_phylo_all_medhts/MLTREE_rootreallyfixed_corhmm_results.csv",stringsAsFactors=FALSE,row.names=1)

#fsummary of best models
table(mat$Best)

#1Rate 2Rate 
#   54   458 

sum.tab<-data.frame(matrix(nrow=nrow(mat),ncol=207),stringsAsFactors=FALSE)

cols<-c(paste("GainNodeBegin",1:25,sep="_"),paste("GainNodeEnd",1:25,sep="_"),paste("GainBegin",1:25,sep="_"),paste("GainEnd",1:25,sep="_"),paste("LossNodeBegin",1:25,sep="_"),paste("LossNodeEnd",1:25,sep="_"),paste("LossBegin",1:25,sep="_"),paste("LossEnd",1:25,sep="_"))
colnames(sum.tab)<-c("TreeNo","NoRates","RootState","NoGains","NoLosses","MaxAgeStartLich","MaxAgeEndLich",cols)


#if root is 0,2

for(x in 1:nrow(mat)){
	#open the tree w. best
	tr<-read.tree(file=paste(tree.file.path,"MLTREE_",mat$Best[x],"_rootreallyfixed_corhmm_states.tre",sep=""))
	#tr.re<-tr
	#tr.re$node.labels<-(length(tr$tip.label)+1):(nrow(tr.tab.mat)+length(tr$tip.label))
	#open attached states file
	#sts<-read.csv(file=paste(tree.file.path,"bs_",x,"_",mat$Best[x],"_corhmm_states.csv",sep=""),stringsAsFactors=FALSE,row.names=1)
	tr.tab.df<-data.frame(branching.times(tr),stringsAsFactors=FALSE)
	rownames(tr.tab.df)<-(length(tr$tip.label)+1):(nrow(tr.tab.df)+length(tr$tip.label))
	tr.edge.df<-data.frame(tr$edge,stringsAsFactors=FALSE)
	colnames(tr.edge.df)<-c("From","To")
	tr.edge.df[,c("FromAge","ToAge","FromState","ToState","Change")]<-NA
	for(z in 1:nrow(tr.edge.df)){
		tr.edge.df$FromAge[z]<-tr.tab.df[rownames(tr.tab.df) %in% tr.edge.df$From[z],]
		tr.edge.df$FromState[z]<-tr$node.label[tr.edge.df$From[z]-length(tr$node.label)-1]
		if(tr.edge.df$To[z]<=length(tr$tip.label)){
			tr.edge.df$ToAge[z]<-0
			tr.edge.df$ToState[z]<-tip.states$lichen[tip.states$edited_name %in% tr$tip.label[tr.edge.df$To[z]]]
			if(tr.edge.df$ToState[z] %in% "0"){
				tr.edge.df$ToState[z]<-"Non"
			}
			if(tr.edge.df$ToState[z] %in% "1"){
				tr.edge.df$ToState[z]<-"Lich"
			}
			if(tr.edge.df$FromState[z] %in% c(1,3) && tr.edge.df$ToState[z] %in% "Non"){
				tr.edge.df$Change[z]<-"NoChange"
			}
			if(tr.edge.df$FromState[z] %in% c(1,3) && tr.edge.df$ToState[z] %in% "Lich"){
				tr.edge.df$Change[z]<-"Gain"
			}
			if(tr.edge.df$FromState[z] %in% c(2,4) && tr.edge.df$ToState[z] %in% "Non"){
				tr.edge.df$Change[z]<-"Loss"
			}
			if(tr.edge.df$FromState[z] %in% c(2,4) && tr.edge.df$ToState[z] %in% "Lich"){
				tr.edge.df$Change[z]<-"NoChange"
			}
		}
		if(tr.edge.df$To[z]>length(tr$tip.label)){
			tr.edge.df$ToAge[z]<-tr.tab.df[rownames(tr.tab.df) %in% tr.edge.df$To[z],]
			tr.edge.df$ToState[z]<-tr$node.label[tr.edge.df$To[z]-length(tr$node.label)-1]
			if(tr.edge.df$FromState[z] %in% c(1,3) && tr.edge.df$ToState[z] %in% c(1,3)){
				tr.edge.df$Change[z]<-"NoChange"
			}
			if(tr.edge.df$FromState[z] %in% c(1,3) && tr.edge.df$ToState[z] %in% c(2,4)){
				tr.edge.df$Change[z]<-"Gain"
			}
			if(tr.edge.df$FromState[z] %in% c(2,4) && tr.edge.df$ToState[z] %in% c(1,3)){
				tr.edge.df$Change[z]<-"Loss"
			}
			if(tr.edge.df$FromState[z] %in% c(2,4) && tr.edge.df$ToState[z] %in% c(2,4)){
				tr.edge.df$Change[z]<-"NoChange"
			}
		}
	}	
	sum.tab$TreeNo[x]<-x
	sum.tab$NoRates[x]<-mat$Best[x]
	if(sum.tab$NoRates[x] %in% "1Rate"){
		if(tr.edge.df$FromState[1] %in% "1"){
			sum.tab$RootState[x]<-"Non"
		}
		if(tr.edge.df$FromState[1] %in% "2"){
			sum.tab$RootState[x]<-"Lich"
		}
	}
	if(sum.tab$NoRates[x] %in% "2Rate"){
		if(tr.edge.df$FromState[1] %in% c(1,3)){
			sum.tab$RootState[x]<-"Non"
		}
		if(tr.edge.df$FromState[1] %in% c(2,4)){
			sum.tab$RootState[x]<-"Lich"
		}
	}
	if(sum.tab$NoRates[x] %in% "3Rate"){
		if(tr.edge.df$FromState[1] %in% c(1,3,5)){
			sum.tab$RootState[x]<-"Non"
		}
		if(tr.edge.df$FromState[1] %in% c(2,4,6)){
			sum.tab$RootState[x]<-"Lich"
		}
	}
	if(sum.tab$RootState[x] %in% "Non"){
		gains<-tr.edge.df[tr.edge.df$Change=="Gain",]
		sum.tab$NoGains[x]<-nrow(gains)
		sum.tab$MaxAgeStartLich[x]<-max(gains$FromAge)
		sum.tab$MaxAgeEndLich[x]<-max(gains$ToAge)
		for(p in 1:nrow(gains)){
			sum.tab[x,paste("GainNodeBegin",p,sep="_")]<-gains$From[p]
			sum.tab[x,paste("GainNodeEnd",p,sep="_")]<-gains$To[p]
			sum.tab[x,paste("GainBegin",p,sep="_")]<-gains$FromAge[p]
			sum.tab[x,paste("GainEnd",p,sep="_")]<-gains$ToAge[p]
		}	
	}
	if(sum.tab$RootState[x] %in% "Lich"){
		if(nrow(tr.edge.df[tr.edge.df$Change=="Gain",])==0){
			sum.tab$NoGains[x]<-0
			sum.tab$MaxAgeStartLich[x]<-tr.edge.df$FromAge[1]
			sum.tab$MaxAgeEndLich[x]<-tr.edge.df$FromAge[1]
			sum.tab[x,paste("GainNodeBegin",p,sep="_")]<-tr.edge.df$From[1]
			sum.tab[x,paste("GainNodeEnd",p,sep="_")]<-tr.edge.df$From[1]
			sum.tab[x,paste("GainBegin",p,sep="_")]<-tr.edge.df$FromAge[1]
			sum.tab[x,paste("GainEnd",p,sep="_")]<-tr.edge.df$FromAge[1]		
		}
		if(nrow(tr.edge.df[tr.edge.df$Change=="Gain",])>0){
			gains<-tr.edge.df[tr.edge.df$Change=="Gain",]
			sum.tab$MaxAgeStartLich[x]<-tr.edge.df$FromAge[1]
			sum.tab$MaxAgeEndLich[x]<-tr.edge.df$FromAge[1]
			sum.tab[x,paste("GainNodeBegin",p,sep="_")]<-tr.edge.df$From[1]
			sum.tab[x,paste("GainNodeEnd",p,sep="_")]<-tr.edge.df$From[1]
			sum.tab[x,paste("GainBegin",p,sep="_")]<-tr.edge.df$FromAge[1]
			sum.tab[x,paste("GainEnd",p,sep="_")]<-tr.edge.df$FromAge[1]		
			sum.tab$NoGains[x]<-nrow(gains)
			for(p in 1:nrow(gains)){
				sum.tab[x,paste("GainNodeBegin",p,sep="_")]<-gains$From[p]
				sum.tab[x,paste("GainNodeEnd",p,sep="_")]<-gains$To[p]
				sum.tab[x,paste("GainBegin",p,sep="_")]<-gains$FromAge[p]
				sum.tab[x,paste("GainEnd",p,sep="_")]<-gains$ToAge[p]
			}		
		}
	}
	if(nrow(tr.edge.df[tr.edge.df$Change=="Loss",])>0){
		losses<-tr.edge.df[tr.edge.df$Change=="Loss",]
		sum.tab$NoLosses[x]<-nrow(losses)
		for(p in 1:nrow(losses)){
			sum.tab[x,paste("LossNodeBegin",p,sep="_")]<-losses$From[p]
			sum.tab[x,paste("LossNodeEnd",p,sep="_")]<-losses$To[p]
			sum.tab[x,paste("LossBegin",p,sep="_")]<-losses$FromAge[p]
			sum.tab[x,paste("LossEnd",p,sep="_")]<-losses$ToAge[p]
		}
	}
	if(nrow(tr.edge.df[tr.edge.df$Change=="Loss",])==0){
		sum.tab$NoLosses[x]<-0
		sum.tab$MaxAgeStartLoss[x]<-0
		sum.tab$MaxAgeEndLoss[x]<-0
	}
}

sum(sum.tab$MaxAgeStartLich>443)
write.csv(sum.tab,file="~/Desktop/asco_phylo_all_medhts/MLTree_corhmm_summary_root_really_fixed.csv",row.names=FALSE)

##################################################################################################################################
#############Record changes in ml tree timescaled with sample of trees (mltree_timescaled_w_sample)
##################################################################################################################################
require(ape)
require(phytools)
require(phangorn)

tree.file.path<-"~/Desktop/deep_lichens/asco_phylo_ml_100medhts/asco_phylo_ml_100medhts/"
tip.states<-read.csv(file="~/Desktop/deep_lichens/asco_phylo_ml_100medhts/asco_phylo_ml_100medhts/bad_idea_lichenization_for_corhmm.csv",stringsAsFactors=FALSE)

#read csv file
mat<-read.csv(file="~/Desktop/deep_lichens/asco_phylo_ml_100medhts/asco_phylo_ml_100medhts/rootreallyfixed_corhmm_results.csv",stringsAsFactors=FALSE)
cn<-c("Tree",colnames(mat)[2:ncol(mat)])
colnames(mat)<-cn
mat.clean<-mat[!is.na(mat$Best),]
write.csv(mat.clean,file="~/Desktop/deep_lichens/asco_phylo_ml_100medhts/asco_phylo_ml_100medhts/rootreallyfixed_corhmm_results_CLEAN.csv",row.names=FALSE)

mat<-read.csv(file="~/Desktop/deep_lichens/asco_phylo_ml_100medhts/asco_phylo_ml_100medhts/rootreallyfixed_corhmm_results_CLEAN.csv",stringsAsFactors=FALSE)


#fsummary of best models
table(mat$Best)

#1Rate 2Rate 
#   54   458 

sum.tab<-data.frame(matrix(nrow=nrow(mat),ncol=207),stringsAsFactors=FALSE)

cols<-c(paste("GainNodeBegin",1:25,sep="_"),paste("GainNodeEnd",1:25,sep="_"),paste("GainBegin",1:25,sep="_"),paste("GainEnd",1:25,sep="_"),paste("LossNodeBegin",1:25,sep="_"),paste("LossNodeEnd",1:25,sep="_"),paste("LossBegin",1:25,sep="_"),paste("LossEnd",1:25,sep="_"))
colnames(sum.tab)<-c("TreeNo","NoRates","RootState","NoGains","NoLosses","MaxAgeStartLich","MaxAgeEndLich",cols)

notrees<-mat$Tree

#if root is 0,2

for(x in 1:nrow(mat)){
	#open the tree w. best
	tr<-read.tree(file=paste(tree.file.path,"ML_w_",mat$Tree[x],"_",mat$Best[x],"_ages_rootreallyfixed_corhmm_states.tre",sep=""))
	#tr.re<-tr
	#tr.re$node.labels<-(length(tr$tip.label)+1):(nrow(tr.tab.mat)+length(tr$tip.label))
	#open attached states file
	#sts<-read.csv(file=paste(tree.file.path,"bs_",x,"_",mat$Best[x],"_corhmm_states.csv",sep=""),stringsAsFactors=FALSE,row.names=1)
	tr.tab.df<-data.frame(branching.times(tr),stringsAsFactors=FALSE)
	rownames(tr.tab.df)<-(length(tr$tip.label)+1):(nrow(tr.tab.df)+length(tr$tip.label))
	tr.edge.df<-data.frame(tr$edge,stringsAsFactors=FALSE)
	colnames(tr.edge.df)<-c("From","To")
	tr.edge.df[,c("FromAge","ToAge","FromState","ToState","Change")]<-NA
	for(z in 1:nrow(tr.edge.df)){
		tr.edge.df$FromAge[z]<-tr.tab.df[rownames(tr.tab.df) %in% tr.edge.df$From[z],]
		tr.edge.df$FromState[z]<-tr$node.label[tr.edge.df$From[z]-length(tr$node.label)-1]
		if(tr.edge.df$To[z]<=length(tr$tip.label)){
			tr.edge.df$ToAge[z]<-0
			tr.edge.df$ToState[z]<-tip.states$lichen[tip.states$edited_name %in% tr$tip.label[tr.edge.df$To[z]]]
			if(tr.edge.df$ToState[z] %in% "0"){
				tr.edge.df$ToState[z]<-"Non"
			}
			if(tr.edge.df$ToState[z] %in% "1"){
				tr.edge.df$ToState[z]<-"Lich"
			}
			if(tr.edge.df$FromState[z] %in% c(1,3) && tr.edge.df$ToState[z] %in% "Non"){
				tr.edge.df$Change[z]<-"NoChange"
			}
			if(tr.edge.df$FromState[z] %in% c(1,3) && tr.edge.df$ToState[z] %in% "Lich"){
				tr.edge.df$Change[z]<-"Gain"
			}
			if(tr.edge.df$FromState[z] %in% c(2,4) && tr.edge.df$ToState[z] %in% "Non"){
				tr.edge.df$Change[z]<-"Loss"
			}
			if(tr.edge.df$FromState[z] %in% c(2,4) && tr.edge.df$ToState[z] %in% "Lich"){
				tr.edge.df$Change[z]<-"NoChange"
			}
		}
		if(tr.edge.df$To[z]>length(tr$tip.label)){
			tr.edge.df$ToAge[z]<-tr.tab.df[rownames(tr.tab.df) %in% tr.edge.df$To[z],]
			tr.edge.df$ToState[z]<-tr$node.label[tr.edge.df$To[z]-length(tr$node.label)-1]
			if(tr.edge.df$FromState[z] %in% c(1,3) && tr.edge.df$ToState[z] %in% c(1,3)){
				tr.edge.df$Change[z]<-"NoChange"
			}
			if(tr.edge.df$FromState[z] %in% c(1,3) && tr.edge.df$ToState[z] %in% c(2,4)){
				tr.edge.df$Change[z]<-"Gain"
			}
			if(tr.edge.df$FromState[z] %in% c(2,4) && tr.edge.df$ToState[z] %in% c(1,3)){
				tr.edge.df$Change[z]<-"Loss"
			}
			if(tr.edge.df$FromState[z] %in% c(2,4) && tr.edge.df$ToState[z] %in% c(2,4)){
				tr.edge.df$Change[z]<-"NoChange"
			}
		}
	}	
	sum.tab$TreeNo[x]<-mat$Tree[x]
	sum.tab$NoRates[x]<-mat$Best[x]
	if(sum.tab$NoRates[x] %in% "1Rate"){
		if(tr.edge.df$FromState[1] %in% "1"){
			sum.tab$RootState[x]<-"Non"
		}
		if(tr.edge.df$FromState[1] %in% "2"){
			sum.tab$RootState[x]<-"Lich"
		}
	}
	if(sum.tab$NoRates[x] %in% "2Rate"){
		if(tr.edge.df$FromState[1] %in% c(1,3)){
			sum.tab$RootState[x]<-"Non"
		}
		if(tr.edge.df$FromState[1] %in% c(2,4)){
			sum.tab$RootState[x]<-"Lich"
		}
	}
	if(sum.tab$RootState[x] %in% "Non"){
		gains<-tr.edge.df[tr.edge.df$Change=="Gain",]
		sum.tab$NoGains[x]<-nrow(gains)
		sum.tab$MaxAgeStartLich[x]<-max(gains$FromAge)
		sum.tab$MaxAgeEndLich[x]<-max(gains$ToAge)
		for(p in 1:nrow(gains)){
			sum.tab[x,paste("GainNodeBegin",p,sep="_")]<-gains$From[p]
			sum.tab[x,paste("GainNodeEnd",p,sep="_")]<-gains$To[p]
			sum.tab[x,paste("GainBegin",p,sep="_")]<-gains$FromAge[p]
			sum.tab[x,paste("GainEnd",p,sep="_")]<-gains$ToAge[p]
		}	
	}
	if(sum.tab$RootState[x] %in% "Lich"){
		if(nrow(tr.edge.df[tr.edge.df$Change=="Gain",])==0){
			sum.tab$NoGains[x]<-0
			sum.tab$MaxAgeStartLich[x]<-tr.edge.df$FromAge[1]
			sum.tab$MaxAgeEndLich[x]<-tr.edge.df$FromAge[1]
			sum.tab[x,paste("GainNodeBegin",p,sep="_")]<-tr.edge.df$From[1]
			sum.tab[x,paste("GainNodeEnd",p,sep="_")]<-tr.edge.df$From[1]
			sum.tab[x,paste("GainBegin",p,sep="_")]<-tr.edge.df$FromAge[1]
			sum.tab[x,paste("GainEnd",p,sep="_")]<-tr.edge.df$FromAge[1]		
		}
		if(nrow(tr.edge.df[tr.edge.df$Change=="Gain",])>0){
			gains<-tr.edge.df[tr.edge.df$Change=="Gain",]
			sum.tab$MaxAgeStartLich[x]<-tr.edge.df$FromAge[1]
			sum.tab$MaxAgeEndLich[x]<-tr.edge.df$FromAge[1]
			sum.tab[x,paste("GainNodeBegin",p,sep="_")]<-tr.edge.df$From[1]
			sum.tab[x,paste("GainNodeEnd",p,sep="_")]<-tr.edge.df$From[1]
			sum.tab[x,paste("GainBegin",p,sep="_")]<-tr.edge.df$FromAge[1]
			sum.tab[x,paste("GainEnd",p,sep="_")]<-tr.edge.df$FromAge[1]		
			sum.tab$NoGains[x]<-nrow(gains)
			for(p in 1:nrow(gains)){
				sum.tab[x,paste("GainNodeBegin",p,sep="_")]<-gains$From[p]
				sum.tab[x,paste("GainNodeEnd",p,sep="_")]<-gains$To[p]
				sum.tab[x,paste("GainBegin",p,sep="_")]<-gains$FromAge[p]
				sum.tab[x,paste("GainEnd",p,sep="_")]<-gains$ToAge[p]
			}		
		}
	}
	if(nrow(tr.edge.df[tr.edge.df$Change=="Loss",])>0){
		losses<-tr.edge.df[tr.edge.df$Change=="Loss",]
		sum.tab$NoLosses[x]<-nrow(losses)
		for(p in 1:nrow(losses)){
			sum.tab[x,paste("LossNodeBegin",p,sep="_")]<-losses$From[p]
			sum.tab[x,paste("LossNodeEnd",p,sep="_")]<-losses$To[p]
			sum.tab[x,paste("LossBegin",p,sep="_")]<-losses$FromAge[p]
			sum.tab[x,paste("LossEnd",p,sep="_")]<-losses$ToAge[p]
		}
	}
	if(nrow(tr.edge.df[tr.edge.df$Change=="Loss",])==0){
		sum.tab$NoLosses[x]<-0
		sum.tab$MaxAgeStartLoss[x]<-0
		sum.tab$MaxAgeEndLoss[x]<-0
	}
}

sum(sum.tab$MaxAgeStartLich>425)
write.csv(sum.tab,file="~/Desktop/deep_lichens/asco_phylo_ml_100medhts/asco_phylo_ml_100medhts/corhmm_summary_root_really_fixed.csv",row.names=FALSE)


##################################################################################################################################
#############Record changes in bootstrap-derived topologies timescaled with median node ages (boots_timescaled_w_median)
##################################################################################################################################
require(ape)
require(phytools)
require(phangorn)

tree.file.path<-"~/Desktop/deep_lichens/asco_phylo/individual_bs_trees/"
tip.states<-read.csv(file="~/Desktop/deep_lichens/asco_phylo/individual_bs_trees/bad_idea_lichenization_for_corhmm.csv",stringsAsFactors=FALSE)

#read csv file
mat<-read.csv(file="~/Desktop/deep_lichens/asco_phylo/individual_bs_trees/rootreallyfixed_corhmm_results.csv",stringsAsFactors=FALSE,row.names=1)

#fsummary of best models
table(mat$Best)

#1Rate 2Rate 
#   54   458 

sum.tab<-data.frame(matrix(nrow=nrow(mat),ncol=207),stringsAsFactors=FALSE)

cols<-c(paste("GainNodeBegin",1:25,sep="_"),paste("GainNodeEnd",1:25,sep="_"),paste("GainBegin",1:25,sep="_"),paste("GainEnd",1:25,sep="_"),paste("LossNodeBegin",1:25,sep="_"),paste("LossNodeEnd",1:25,sep="_"),paste("LossBegin",1:25,sep="_"),paste("LossEnd",1:25,sep="_"))
colnames(sum.tab)<-c("TreeNo","NoRates","RootState","NoGains","NoLosses","MaxAgeStartLich","MaxAgeEndLich",cols)


#if root is 0,2

for(x in 1:nrow(mat)){
	#open the tree w. best
	tr<-read.tree(file=paste(tree.file.path,"bs_",x,"_",mat$Best[x],"_rootreallyfixed_corhmm_states.tre",sep=""))
	#tr.re<-tr
	#tr.re$node.labels<-(length(tr$tip.label)+1):(nrow(tr.tab.mat)+length(tr$tip.label))
	#open attached states file
	#sts<-read.csv(file=paste(tree.file.path,"bs_",x,"_",mat$Best[x],"_corhmm_states.csv",sep=""),stringsAsFactors=FALSE,row.names=1)
	tr.tab.df<-data.frame(branching.times(tr),stringsAsFactors=FALSE)
	rownames(tr.tab.df)<-(length(tr$tip.label)+1):(nrow(tr.tab.df)+length(tr$tip.label))
	tr.edge.df<-data.frame(tr$edge,stringsAsFactors=FALSE)
	colnames(tr.edge.df)<-c("From","To")
	tr.edge.df[,c("FromAge","ToAge","FromState","ToState","Change")]<-NA
	for(z in 1:nrow(tr.edge.df)){
		tr.edge.df$FromAge[z]<-tr.tab.df[rownames(tr.tab.df) %in% tr.edge.df$From[z],]
		tr.edge.df$FromState[z]<-tr$node.label[tr.edge.df$From[z]-length(tr$node.label)-1]
		if(tr.edge.df$To[z]<=length(tr$tip.label)){
			tr.edge.df$ToAge[z]<-0
			tr.edge.df$ToState[z]<-tip.states$lichen[tip.states$edited_name %in% tr$tip.label[tr.edge.df$To[z]]]
			if(tr.edge.df$ToState[z] %in% "0"){
				tr.edge.df$ToState[z]<-"Non"
			}
			if(tr.edge.df$ToState[z] %in% "1"){
				tr.edge.df$ToState[z]<-"Lich"
			}
			if(tr.edge.df$FromState[z] %in% c(1,3) && tr.edge.df$ToState[z] %in% "Non"){
				tr.edge.df$Change[z]<-"NoChange"
			}
			if(tr.edge.df$FromState[z] %in% c(1,3) && tr.edge.df$ToState[z] %in% "Lich"){
				tr.edge.df$Change[z]<-"Gain"
			}
			if(tr.edge.df$FromState[z] %in% c(2,4) && tr.edge.df$ToState[z] %in% "Non"){
				tr.edge.df$Change[z]<-"Loss"
			}
			if(tr.edge.df$FromState[z] %in% c(2,4) && tr.edge.df$ToState[z] %in% "Lich"){
				tr.edge.df$Change[z]<-"NoChange"
			}
		}
		if(tr.edge.df$To[z]>length(tr$tip.label)){
			tr.edge.df$ToAge[z]<-tr.tab.df[rownames(tr.tab.df) %in% tr.edge.df$To[z],]
			tr.edge.df$ToState[z]<-tr$node.label[tr.edge.df$To[z]-length(tr$node.label)-1]
			if(tr.edge.df$FromState[z] %in% c(1,3) && tr.edge.df$ToState[z] %in% c(1,3)){
				tr.edge.df$Change[z]<-"NoChange"
			}
			if(tr.edge.df$FromState[z] %in% c(1,3) && tr.edge.df$ToState[z] %in% c(2,4)){
				tr.edge.df$Change[z]<-"Gain"
			}
			if(tr.edge.df$FromState[z] %in% c(2,4) && tr.edge.df$ToState[z] %in% c(1,3)){
				tr.edge.df$Change[z]<-"Loss"
			}
			if(tr.edge.df$FromState[z] %in% c(2,4) && tr.edge.df$ToState[z] %in% c(2,4)){
				tr.edge.df$Change[z]<-"NoChange"
			}
		}
	}	
	sum.tab$TreeNo[x]<-x
	sum.tab$NoRates[x]<-mat$Best[x]
	if(sum.tab$NoRates[x] %in% "1Rate"){
		if(tr.edge.df$FromState[1] %in% "1"){
			sum.tab$RootState[x]<-"Non"
		}
		if(tr.edge.df$FromState[1] %in% "2"){
			sum.tab$RootState[x]<-"Lich"
		}
	}
	if(sum.tab$NoRates[x] %in% "2Rate"){
		if(tr.edge.df$FromState[1] %in% c(1,3)){
			sum.tab$RootState[x]<-"Non"
		}
		if(tr.edge.df$FromState[1] %in% c(2,4)){
			sum.tab$RootState[x]<-"Lich"
		}
	}
	if(sum.tab$RootState[x] %in% "Non"){
		gains<-tr.edge.df[tr.edge.df$Change=="Gain",]
		sum.tab$NoGains[x]<-nrow(gains)
		sum.tab$MaxAgeStartLich[x]<-max(gains$FromAge)
		sum.tab$MaxAgeEndLich[x]<-max(gains$ToAge)
		for(p in 1:nrow(gains)){
			sum.tab[x,paste("GainNodeBegin",p,sep="_")]<-gains$From[p]
			sum.tab[x,paste("GainNodeEnd",p,sep="_")]<-gains$To[p]
			sum.tab[x,paste("GainBegin",p,sep="_")]<-gains$FromAge[p]
			sum.tab[x,paste("GainEnd",p,sep="_")]<-gains$ToAge[p]
		}	
	}
	if(sum.tab$RootState[x] %in% "Lich"){
		if(nrow(tr.edge.df[tr.edge.df$Change=="Gain",])==0){
			sum.tab$NoGains[x]<-0
			sum.tab$MaxAgeStartLich[x]<-tr.edge.df$FromAge[1]
			sum.tab$MaxAgeEndLich[x]<-tr.edge.df$FromAge[1]
			sum.tab[x,paste("GainNodeBegin",p,sep="_")]<-tr.edge.df$From[1]
			sum.tab[x,paste("GainNodeEnd",p,sep="_")]<-tr.edge.df$From[1]
			sum.tab[x,paste("GainBegin",p,sep="_")]<-tr.edge.df$FromAge[1]
			sum.tab[x,paste("GainEnd",p,sep="_")]<-tr.edge.df$FromAge[1]		
		}
		if(nrow(tr.edge.df[tr.edge.df$Change=="Gain",])>0){
			gains<-tr.edge.df[tr.edge.df$Change=="Gain",]
			sum.tab$MaxAgeStartLich[x]<-tr.edge.df$FromAge[1]
			sum.tab$MaxAgeEndLich[x]<-tr.edge.df$FromAge[1]
			sum.tab[x,paste("GainNodeBegin",p,sep="_")]<-tr.edge.df$From[1]
			sum.tab[x,paste("GainNodeEnd",p,sep="_")]<-tr.edge.df$From[1]
			sum.tab[x,paste("GainBegin",p,sep="_")]<-tr.edge.df$FromAge[1]
			sum.tab[x,paste("GainEnd",p,sep="_")]<-tr.edge.df$FromAge[1]		
			sum.tab$NoGains[x]<-nrow(gains)
			for(p in 1:nrow(gains)){
				sum.tab[x,paste("GainNodeBegin",p,sep="_")]<-gains$From[p]
				sum.tab[x,paste("GainNodeEnd",p,sep="_")]<-gains$To[p]
				sum.tab[x,paste("GainBegin",p,sep="_")]<-gains$FromAge[p]
				sum.tab[x,paste("GainEnd",p,sep="_")]<-gains$ToAge[p]
			}		
		}
	}
	if(nrow(tr.edge.df[tr.edge.df$Change=="Loss",])>0){
		losses<-tr.edge.df[tr.edge.df$Change=="Loss",]
		sum.tab$NoLosses[x]<-nrow(losses)
		for(p in 1:nrow(losses)){
			sum.tab[x,paste("LossNodeBegin",p,sep="_")]<-losses$From[p]
			sum.tab[x,paste("LossNodeEnd",p,sep="_")]<-losses$To[p]
			sum.tab[x,paste("LossBegin",p,sep="_")]<-losses$FromAge[p]
			sum.tab[x,paste("LossEnd",p,sep="_")]<-losses$ToAge[p]
		}
	}
	if(nrow(tr.edge.df[tr.edge.df$Change=="Loss",])==0){
		sum.tab$NoLosses[x]<-0
		sum.tab$MaxAgeStartLoss[x]<-0
		sum.tab$MaxAgeEndLoss[x]<-0
	}
}

sum(sum.tab$MaxAgeStartLich>443)
write.csv(sum.tab,file="~/Desktop/deep_lichens/asco_phylo/individual_bs_trees/corhmm_summary_root_really_fixed.csv",row.names=FALSE)


##################################################################################################################################
#############Record changes in 100 bootstrap-derived topologies timescaled with posterior sample (50) of trees (boots_timescaled_w_sample)
##################################################################################################################################

tree.file.path<-"~/Desktop/deep_lichens/asco_phylo_ml_100medhts/asco_phylo_boot_25medhts/"
tip.states<-read.csv(file="~/Desktop/deep_lichens/asco_phylo_ml_100medhts/asco_phylo_boot_25medhts/bad_idea_lichenization_for_corhmm.csv",stringsAsFactors=FALSE)

#notrees<-c(153, 244, 1165, 1193, 1423, 1770, 2206, 2451, 2635, 3493, 3602, 3699, 3715, 4229, 4491, 4853, 4976, 4986, 6445, 6586, 6776, 6943, 7081, 7172, 8159)
notrees<-sort(c(153, 244, 1165, 1193, 1423, 1770, 2206, 2451, 2635, 3493, 3602, 3699, 3715, 4229, 4491, 4853, 4976, 4986, 6445, 6586, 6776, 6943, 7081, 7172, 8159, 228, 364, 1260, 1599, 2030, 2701, 2791, 3637, 3914, 4199, 4384, 5042, 5097, 5143, 5220, 5434, 5488, 5568, 6492, 7158, 7434, 7438, 7830, 7928, 8222))



for(b in 1:100){
	#read csv file
	mat.second<-read.csv(file=paste("~/Desktop/deep_lichens/asco_phylo_ml_100medhts/asco_phylo_boot_25medhts/","bs_",b,"_rootreallyfixed_corhmm_results_secondbatch.csv",sep=""),stringsAsFactors=FALSE,row.names=1)
	mat.first<-read.csv(file=paste("~/Desktop/deep_lichens/asco_phylo_ml_100medhts/asco_phylo_boot_25medhts/","bs_",b,"_rootreallyfixed_corhmm_results.csv",sep=""),stringsAsFactors=FALSE,row.names=1)
	mat<-rbind(mat.first,mat.second)
	#then order by Tree
	mat<-mat[order(mat$Tree),]
	table(mat$Best)
	sum.tab<-data.frame(matrix(nrow=nrow(mat),ncol=208),stringsAsFactors=FALSE)
	cols<-c(paste("GainNodeBegin",1:25,sep="_"),paste("GainNodeEnd",1:25,sep="_"),paste("GainBegin",1:25,sep="_"),paste("GainEnd",1:25,sep="_"),paste("LossNodeBegin",1:25,sep="_"),paste("LossNodeEnd",1:25,sep="_"),paste("LossBegin",1:25,sep="_"),paste("LossEnd",1:25,sep="_"))
colnames(sum.tab)<-c("BootNo","CalTreeNo","NoRates","RootState","NoGains","NoLosses","MaxAgeStartLich","MaxAgeEndLich",cols)
	for(x in 1:nrow(mat)){
		#open the tree w. best
		tr<-read.tree(file=paste(tree.file.path,"bs_",b,"_w_",mat$Tree[x],"_",mat$Best[x],"_ages_rootreallyfixed_corhmm_states.tre",sep=""))
		tr.tab.df<-data.frame(branching.times(tr),stringsAsFactors=FALSE)
		rownames(tr.tab.df)<-(length(tr$tip.label)+1):(nrow(tr.tab.df)+length(tr$tip.label))
		tr.edge.df<-data.frame(tr$edge,stringsAsFactors=FALSE)
		colnames(tr.edge.df)<-c("From","To")
		tr.edge.df[,c("FromAge","ToAge","FromState","ToState","Change")]<-NA
		for(z in 1:nrow(tr.edge.df)){
			tr.edge.df$FromAge[z]<-tr.tab.df[rownames(tr.tab.df) %in% tr.edge.df$From[z],]
			tr.edge.df$FromState[z]<-tr$node.label[tr.edge.df$From[z]-length(tr$node.label)-1]
			if(tr.edge.df$To[z]<=length(tr$tip.label)){
				tr.edge.df$ToAge[z]<-0
				tr.edge.df$ToState[z]<-tip.states$lichen[tip.states$edited_name %in% tr$tip.label[tr.edge.df$To[z]]]
				if(tr.edge.df$ToState[z] %in% "0"){
					tr.edge.df$ToState[z]<-"Non"
				}
				if(tr.edge.df$ToState[z] %in% "1"){
					tr.edge.df$ToState[z]<-"Lich"
				}
				if(tr.edge.df$FromState[z] %in% c(1,3) && tr.edge.df$ToState[z] %in% "Non"){
					tr.edge.df$Change[z]<-"NoChange"
				}
				if(tr.edge.df$FromState[z] %in% c(1,3) && tr.edge.df$ToState[z] %in% "Lich"){
					tr.edge.df$Change[z]<-"Gain"
				}
				if(tr.edge.df$FromState[z] %in% c(2,4) && tr.edge.df$ToState[z] %in% "Non"){
					tr.edge.df$Change[z]<-"Loss"
				}
				if(tr.edge.df$FromState[z] %in% c(2,4) && tr.edge.df$ToState[z] %in% "Lich"){
					tr.edge.df$Change[z]<-"NoChange"
				}
			}
			if(tr.edge.df$To[z]>length(tr$tip.label)){
				tr.edge.df$ToAge[z]<-tr.tab.df[rownames(tr.tab.df) %in% tr.edge.df$To[z],]
				tr.edge.df$ToState[z]<-tr$node.label[tr.edge.df$To[z]-length(tr$node.label)-1]
				if(tr.edge.df$FromState[z] %in% c(1,3) && tr.edge.df$ToState[z] %in% c(1,3)){
					tr.edge.df$Change[z]<-"NoChange"
				}
				if(tr.edge.df$FromState[z] %in% c(1,3) && tr.edge.df$ToState[z] %in% c(2,4)){
					tr.edge.df$Change[z]<-"Gain"
				}
				if(tr.edge.df$FromState[z] %in% c(2,4) && tr.edge.df$ToState[z] %in% c(1,3)){
					tr.edge.df$Change[z]<-"Loss"
				}
				if(tr.edge.df$FromState[z] %in% c(2,4) && tr.edge.df$ToState[z] %in% c(2,4)){
					tr.edge.df$Change[z]<-"NoChange"
				}
			}
		}	
		sum.tab$BootNo[x]<-b
		sum.tab$CalTreeNo[x]<-mat$Tree[x]
		sum.tab$NoRates[x]<-mat$Best[x]
		if(sum.tab$NoRates[x] %in% "1Rate"){
			if(tr.edge.df$FromState[1] %in% "1"){
				sum.tab$RootState[x]<-"Non"
			}
			if(tr.edge.df$FromState[1] %in% "2"){
				sum.tab$RootState[x]<-"Lich"
			}
		}
		if(sum.tab$NoRates[x] %in% "2Rate"){
			if(tr.edge.df$FromState[1] %in% c(1,3)){
				sum.tab$RootState[x]<-"Non"
			}
			if(tr.edge.df$FromState[1] %in% c(2,4)){
				sum.tab$RootState[x]<-"Lich"
			}
		}
		if(sum.tab$NoRates[x] %in% "3Rate"){
			if(tr.edge.df$FromState[1] %in% c(1,3,5)){
				sum.tab$RootState[x]<-"Non"
			}
			if(tr.edge.df$FromState[1] %in% c(2,4,6)){
				sum.tab$RootState[x]<-"Lich"
			}
		}
		if(sum.tab$RootState[x] %in% "Non"){
			gains<-tr.edge.df[tr.edge.df$Change=="Gain",]
			sum.tab$NoGains[x]<-nrow(gains)
			sum.tab$MaxAgeStartLich[x]<-max(gains$FromAge)
			sum.tab$MaxAgeEndLich[x]<-max(gains$ToAge)
			for(p in 1:nrow(gains)){
				sum.tab[x,paste("GainNodeBegin",p,sep="_")]<-gains$From[p]
				sum.tab[x,paste("GainNodeEnd",p,sep="_")]<-gains$To[p]
				sum.tab[x,paste("GainBegin",p,sep="_")]<-gains$FromAge[p]
				sum.tab[x,paste("GainEnd",p,sep="_")]<-gains$ToAge[p]
			}	
		}
		if(sum.tab$RootState[x] %in% "Lich"){
			if(nrow(tr.edge.df[tr.edge.df$Change=="Gain",])==0){
				sum.tab$NoGains[x]<-0
				sum.tab$MaxAgeStartLich[x]<-tr.edge.df$FromAge[1]
				sum.tab$MaxAgeEndLich[x]<-tr.edge.df$FromAge[1]
				sum.tab[x,paste("GainNodeBegin",p,sep="_")]<-tr.edge.df$From[1]
				sum.tab[x,paste("GainNodeEnd",p,sep="_")]<-tr.edge.df$From[1]
				sum.tab[x,paste("GainBegin",p,sep="_")]<-tr.edge.df$FromAge[1]
				sum.tab[x,paste("GainEnd",p,sep="_")]<-tr.edge.df$FromAge[1]		
			}
			if(nrow(tr.edge.df[tr.edge.df$Change=="Gain",])>0){
				gains<-tr.edge.df[tr.edge.df$Change=="Gain",]
				sum.tab$MaxAgeStartLich[x]<-tr.edge.df$FromAge[1]
				sum.tab$MaxAgeEndLich[x]<-tr.edge.df$FromAge[1]
				sum.tab[x,paste("GainNodeBegin",p,sep="_")]<-tr.edge.df$From[1]
				sum.tab[x,paste("GainNodeEnd",p,sep="_")]<-tr.edge.df$From[1]
				sum.tab[x,paste("GainBegin",p,sep="_")]<-tr.edge.df$FromAge[1]
				sum.tab[x,paste("GainEnd",p,sep="_")]<-tr.edge.df$FromAge[1]		
				sum.tab$NoGains[x]<-nrow(gains)
				for(p in 1:nrow(gains)){
					sum.tab[x,paste("GainNodeBegin",p,sep="_")]<-gains$From[p]
					sum.tab[x,paste("GainNodeEnd",p,sep="_")]<-gains$To[p]
					sum.tab[x,paste("GainBegin",p,sep="_")]<-gains$FromAge[p]
					sum.tab[x,paste("GainEnd",p,sep="_")]<-gains$ToAge[p]
				}		
			}
		}
		if(nrow(tr.edge.df[tr.edge.df$Change=="Loss",])>0){
			losses<-tr.edge.df[tr.edge.df$Change=="Loss",]
			sum.tab$NoLosses[x]<-nrow(losses)
			for(p in 1:nrow(losses)){
				sum.tab[x,paste("LossNodeBegin",p,sep="_")]<-losses$From[p]
				sum.tab[x,paste("LossNodeEnd",p,sep="_")]<-losses$To[p]
				sum.tab[x,paste("LossBegin",p,sep="_")]<-losses$FromAge[p]
				sum.tab[x,paste("LossEnd",p,sep="_")]<-losses$ToAge[p]
			}
		}
		if(nrow(tr.edge.df[tr.edge.df$Change=="Loss",])==0){
			sum.tab$NoLosses[x]<-0
			sum.tab$MaxAgeStartLoss[x]<-0
			sum.tab$MaxAgeEndLoss[x]<-0
		}
	}
sum(sum.tab$MaxAgeStartLich>443)
write.csv(sum.tab,file=paste("~/Desktop/deep_lichens/asco_phylo_ml_100medhts/asco_phylo_boot_25medhts/","bs_",b,"_corhmm_summary_root_really_fixed_combined_correct.csv",sep=""),row.names=FALSE)
}


#NOW SUM BS Reps
sum.tab<-data.frame(matrix(nrow=0,ncol=208),stringsAsFactors=FALSE)
cols<-c(paste("GainNodeBegin",1:25,sep="_"),paste("GainNodeEnd",1:25,sep="_"),paste("GainBegin",1:25,sep="_"),paste("GainEnd",1:25,sep="_"),paste("LossNodeBegin",1:25,sep="_"),paste("LossNodeEnd",1:25,sep="_"),paste("LossBegin",1:25,sep="_"),paste("LossEnd",1:25,sep="_"))
colnames(sum.tab)<-c("BootNo","CalTreeNo","NoRates","RootState","NoGains","NoLosses","MaxAgeStartLich","MaxAgeEndLich",cols)

for(b in 1:100){
	require(plyr)
	#read csv file
	dat<-read.csv(file=paste("~/Desktop/deep_lichens/asco_phylo_ml_100medhts/asco_phylo_boot_25medhts/","bs_",b,"_corhmm_summary_root_really_fixed_combined_correct.csv",sep=""),stringsAsFactors=FALSE)
	sum.tab<-rbind.fill(sum.tab,dat)
}

write.csv(sum.tab,file="~/Desktop/deep_lichens/asco_phylo_ml_100medhts/asco_phylo_boot_25medhts/boots_corhmm_summary_root_really_fixed_combined_correct.csv",row.names=FALSE)


