#' perform compensation and logicle transformation
#'
#' This function performs logicle transformation on all fcs files in a folder and outputs transforemd data in text format
#'
#' @param input.directory Directory of the folder contains all fcs files
#' @param output.directory Directory of the folder to store all transformed data in text format
#' @param marker.index Indices of markers that need to be logicle transformed
#' @param marker.names Names of markers in the same order as the indices
#' @return NULL
#' @export
preprocess = function(input.directory, output.directory, marker.index, marker.names){
	dir.create(output.directory)
	infiles = list.files(input.directory)
	#----------------------------------------------------gating----------------------------------------------------------------------------------
	for (i in 1:length(infiles)) {
		frame = flowCore::read.FCS(paste0(input.directory,infiles[i]))
		name = sub(".fcs","",infiles[i])
		frame = flowCore::compensate(frame,frame@description$SPILL)
		lgcl = flowCore::estimateLogicle(frame, channels = colnames(frame)[marker.index])
		frame = transform(frame,lgcl)
		print(name)
		outdata = exprs(frame)[,marker.index]
		colnames(outdata) = marker.names
		write.table(outdata,paste0(output.directory,name,".txt"),append=F,row.names=F,sep="\t")
	}
}

#' perform CytoBinning on all files in a folder
#'
#' This function performs CytoBinning on all data files in a folder and output a matrix where each row is the binning result for a data file
#'
#' @param bin.number A vector of selected number of bins
#' @param markers A vector contains names of markers to be used in binning
#' @param input.dir Directory of the folder that contains pre-processed text data files
#' @param output.dir Directory of the folder that contains binning results, each file stores one marker pair of a given bin number
#' @return NULL
#' @export
CytoBinning = function(bin.number, markers, input.dir, output.dir) {
	dir.create(output.dir)
	infiles = list.files(input.dir)
	for (i_file in 1:length(infiles)) {
		indata = read.table(paste0(input.dir,infiles[i_file]),header=T,colClasses="numeric",sep="\t")
		name = sub(".txt","",infiles[i_file])
		for (i_marker in 1:(length(markers)-1)) {
			for (j_marker in (i_marker+1):length(markers)) {
				marker_pair = c(markers[i_marker],markers[j_marker])
				ind = which(colnames(indata)%in%marker_pair)
				to_use = indata[,ind]
				for (i_bin in 1:length(bin.number)) {
					numBin = bin.number[i_bin]
					gates = apply(to_use,2,quantile,c(0:numBin)/numBin)
					temp = c()
					for (i_binn in 1:numBin) {
						for (j_bin in 1:numBin) {
							clause1 <- to_use[,1]>gates[i_binn,1] & to_use[,1]<gates[(i_binn+1),1]
							clause2 <- to_use[,2]>gates[j_bin,2] & to_use[,2]<gates[(j_bin+1),2]
							together <- clause1 & clause2
							percent <- 100*sum(together)/nrow(to_use)
							temp <- c(temp,percent)
						}
					}
					write(file=paste0(output.dir,markers[i_marker],"_",markers[j_marker],"_",numBin,".txt"),c(name,as.numeric(temp)),ncol=(numBin^3+1),sep="\t",append=TRUE)
				}
			}
		}
		print(paste("file",i_file,"done!!!"))
	}
}

#' performs SVM classification and testing
#'
#' This function applies SVM in kernlab package to classify two groups of data points stored in two matrices, 
#' applies its classification boundary to testing dataset and outputs corresponding accuracy.
#'
#' @param controls Training matrix of healthy cells to be classified
#' @param exper Training matrix of diseased cells to be classified
#' @param feature_index A vector contain index of measurements used in classification, should be a vector of integers, length must be larger than 1
#' @param test.control Testing matrix of healthy cells
#' @param test.exp Testing matrix of diseased cells
#' @return accuracy Classificaiton accuracy of the training dataset
#' @return test Classification accuracy of the test dataset
#' @return weight Normalized weights of the features
#' @export
SVM_cv <- function(controls,exper,feature_index,test.control,test.exp) {
	library(kernlab)
	library(ks)
	AvsB <- matrix(1,nrow(controls),1)
	controls <- cbind(controls,AvsB)
	AvsB <- matrix(-1,nrow(exper),1)
	exper <- cbind(exper,AvsB)
	myTable <- rbind(controls,exper)
	AvsB <- factor(myTable[,ncol(myTable)])
	top_measures <- feature_index

	n_meas <- length(top_measures)
	x <- as.matrix(myTable[,top_measures])
	means <- apply(x,2,mean)
	sds <- apply(x,2,sd)
	x <- scale(x)
	# x[is.na(x)] <- 0
	test <- rbind(test.control,test.exp)
	y <- as.matrix(test[,top_measures])
	y <- scale(y,center=means,scale=sds)
	# y[is.na(y)] <- 0
	myTable2 <- data.frame(x,AvsB)
	myModel <- ksvm(AvsB ~ ., data=myTable2,type="C-svc", kernel="vanilladot", C = 10, prob.model=TRUE)

	SVM_coeff <- alpha(myModel)[[1]]
	SVM_coeff2 <- coef(myModel)[[1]]
	SVM_index <- alphaindex(myModel)[[1]]
	weight <- rep(0,each=ncol(x))
	for (i in 1:length(SVM_coeff)) {
		weight <- weight + SVM_coeff2[i]*x[SVM_index[i],]
	}
	norm <- sqrt(sum(weight**2))
	weightnorm <- weight/norm
	weigh <- cbind(as.matrix(weightnorm**2),as.matrix(weightnorm))
	to_weight <- weightnorm**2
	label <- sign(SVM_coeff/SVM_coeff2)
	SVM_b <- 0
	for (i in 1:length(SVM_coeff)) {	
		SVM_b <- SVM_b + (sum(weight*x[SVM_index[i],])-label[i])
	}
	SVM_b <- SVM_b/length(SVM_coeff)
	SVM_bn <- SVM_b/norm

	dist <- rep(0,each=nrow(x))
	for (i in 1:nrow(x)) {
		dist[i] <- sum(weightnorm*x[i,])-SVM_bn
	}
	accuracy <- 100*sum(dist*as.numeric(as.character(AvsB))>0)/nrow(x)
	truth <- c(rep(1,nrow(test.control)),rep(-1,nrow(test.exp)))
	test.dist <- y%*%weigh[,2] - SVM_bn
	CV <- 100*sum(test.dist*truth>0)/nrow(y)
	# gap size
	n_pat_class_1 <- sum(AvsB == levels(AvsB)[2])
	gap_size <- min(dist[1:n_pat_class_1])-max(dist[(n_pat_class_1+1):nrow(x)])
	to_delete <- which(feature_index%in%feature_index[which(to_weight==min(to_weight))])
	results <- list(accuracy, CV, weightnorm)
	names(results) <- c("accuracy","test","weight")
	return(results)
}