#######################results.txt
argv <- commandArgs(TRUE)
file.in <- argv[1]
file.out <- argv[2]
filetype <- argv[3]
ci.in.mutation <- read.table(file=file.in, header=T, comment.char = "", quote = "", colClasses = "character", sep="\t",check=F)
p_value <- c()
CI <- c()
if (dim(ci.in.mutation)[1]>0){
	if(grepl("_",ci.in.mutation$Cov_max[1])){
		for(i in 1:dim(ci.in.mutation)[1]){
			covs <- strsplit(ci.in.mutation$Cov_max[i], "_")
			tumor_cov <- as.numeric(covs[[1]][1])
			var_covs <- strsplit(ci.in.mutation$Var_Cov_max[i], "_")
			tumor_varcov <- as.numeric(var_covs[[1]][1])
			var_freqs <- strsplit(ci.in.mutation$VarFreq_precent_max[i], "_")
			tumor_varfreq <- (as.numeric(var_freqs[[1]][1]))/100
			if((tumor_cov==-1)||(tumor_varcov==-1)){
				p_tumor <- "na"
			}else{
				p_tumor <- binom.test(c(tumor_varcov,tumor_cov-tumor_varcov), p = 0.001, alternative = c("greater"), conf.level = 0.95)
				p_tumor <- formatC(p_tumor$p.value, format="e", digits=2)
			}
			p_value <- c(p_value, p_tumor)
			if (tumor_cov==0){
				ci_tumor_lower <- 0.00
				ci_tumor_uper <- 0.00
			} else if((tumor_cov==-1)||(tumor_varcov==-1)){
				ci_tumor_lower <- "na"
				ci_tumor_uper <- "na"
			}else{
				ci_tumor <- sqrt(tumor_varfreq*(1-tumor_varfreq)/tumor_cov)
				ci_tumor_lower <- round(((tumor_varfreq-1.96*ci_tumor-0.5/tumor_cov)*100),2)
				ci_tumor_lower <- formatC(ci_tumor_lower, format="f", digits=2)
				ci_tumor_uper <- round(((tumor_varfreq+1.96*ci_tumor+0.5/tumor_cov)*100),2)
				ci_tumor_uper <- formatC(ci_tumor_uper, format="f", digits=2)
			}
			ci_value <- paste("[",ci_tumor_lower,",",ci_tumor_uper,"]",sep="")
			CI <- c(CI,ci_value)
		}
	}else{
		for(i in 1:dim(ci.in.mutation)[1]){
			tumor_cov <- as.numeric(ci.in.mutation$Cov_max[i])
			tumor_varcov <- as.numeric(ci.in.mutation$Var_Cov_max[i])
			tumor_varfreq <- (as.numeric(ci.in.mutation$VarFreq_precent_max[i]))/100
			software_source <- ci.in.mutation$Source[i]
			if((tumor_cov==-1)||(tumor_varcov==-1)){
				p_tumor <- "na"
			}else if (software_source=="Mpileup"){
				if(tumor_cov>=5){
					p_tumor <- binom.test(c(tumor_varcov,tumor_cov-tumor_varcov), p = 0.001, alternative = c("greater"), conf.level = 0.95)
					p_tumor <- formatC(p_tumor$p.value, format="e", digits=2)
				}else{
					p_tumor <- 1.00
				}
			}else{
				#print(tumor_cov)
				#print(tumor_varcov)
				p_tumor <- binom.test(c(tumor_varcov,tumor_cov-tumor_varcov), p = 0.001, alternative = c("greater"), conf.level = 0.95)
				p_tumor <- formatC(p_tumor$p.value, format="e", digits=2)
			}
			p_value <- c(p_value, p_tumor)
			if (tumor_cov==0){
				ci_tumor_lower <- 0.00
				ci_tumor_uper <- 0.00
				p_tumor <- 0.00
			}else if ((tumor_cov==-1)||(tumor_varcov==-1)){
				ci_tumor_lower <- "na"
				ci_tumor_uper <- "na"
			}else{
				ci_tumor <- sqrt(tumor_varfreq*(1-tumor_varfreq)/tumor_cov)
				ci_tumor_lower <- round(((tumor_varfreq-1.96*ci_tumor-0.5/tumor_cov)*100),2)
				ci_tumor_lower <- formatC(ci_tumor_lower, format="f", digits=2)
				ci_tumor_uper <- round(((tumor_varfreq+1.96*ci_tumor+0.5/tumor_cov)*100),2)
				ci_tumor_uper <- formatC(ci_tumor_uper, format="f", digits=2)
			}
			ci_value <- paste("[",ci_tumor_lower,",",ci_tumor_uper,"]",sep="")
			CI <- c(CI,ci_value)
		}
	}
	ci.in.mutation.p <- cbind(ci.in.mutation, CI, p_value)
	if ((filetype!=0)&&(filetype!=4)){
		ci.in.mutation.sig <- ci.in.mutation.p[(as.numeric(as.vector(ci.in.mutation.p$p_value))<0.001), ]
	}
	if ((filetype==0)||(filetype==4)){
		write.table(file=file.out, ci.in.mutation.p, sep="\t", quote = F, row.names=F)
	}else{
		write.table(file=file.out, ci.in.mutation.sig, sep="\t", quote = F, row.names=F)
	}
}else{
	header<-c(names(ci.in.mutation),"CI")
	header<-c(header,"p_value")
	write.table(file=file.out, t(header), sep="\t", quote = F, row.names=F,col.names=F)
}
