##This will summarize MEDUSA results from a multitree run. You'll need a table with "clean" results i.e. columns for Tree No., Shift order, Div rate, Extinction Rate and Clade/group.

data<-read.table("/Users/SAI/Documents/Molecular/MedusaMANYresults/Medusa_Many_clean", header = TRUE)

maxTree<-max(data[,1])

sets<-character(length=maxTree)
for(i in 1:maxTree) {
	ok<-data[,1]==i
	sets[i]<-paste(sort(data[ok,5]), collapse="")
	}
	
sets<-as.factor(sets)

table(sets)	

for(i in levels(sets)) {
	cat("Set: ", i, "Trees: \n")
	print(which(sets==i))
	}

for(i in levels(sets)[-1]) {
	trees<-which(sets==i)
	ok<-data[,1]%in%trees
	write.csv(file="Data_in_sets.txt", data[ok,], append=T)
	}

# This part will give you statistics about everything

for(i in levels(sets)[-1]) {
	trees<-which(sets==i)
	ok<-data[,1]%in%trees
	dd<-data[ok,]
	mm<-aggregate(dd[,3:4], by=list(dd[,5]), FUN=mean)
	vv<-aggregate(dd[,3:4], by=list(dd[,5]), FUN=var)
	nn<-aggregate(dd[,3:4], by=list(dd[,5]), FUN=length)
	res<-cbind(mm, vv[,-1], nn[,2])
	colnames(res)[2:6]<-c("diver_mean", "extinc_mean", "diver_var", "extinc_var", "n")
	print(res)

	}