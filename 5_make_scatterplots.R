library(dplyr)
library(ggplot2)
library(stringr)


for (dir in c("human/bigwigs","human/bigwigs_subtract", "human/bigwigs_nopseudo","mouse/bigwigs","mouse/bigwigs_subtract", "mouse/bigwigs_nopseudo")){
	
	dir.create(paste0(dir, "_plots/scatterplots/"), showWarnings = FALSE)
	dir.create(paste0(dir, "_plots/scatterplots_GFP/"), showWarnings = FALSE)

	means <- read.table(paste0(dir, "_plots/average_signal.txt"), sep="\t", header=T)

	if (str_detect(dir, "human")) {
		species = "hs"
	}else{
		species = "mm"
	}
	samples <- unique(paste0(means$Sample, "_", means$Method))
	samples2 <- unique(paste0(means$Sample, "_", means$Method, "-", means$Tag))

	#- sample_plus versus sample_minus
	#- Do this for any sample that has both plus and minus data
	for (sample in samples){
		if ((paste0(sample, "-", "PlusTag") %in% samples2) && (paste0(sample, "-", "MinusTag") %in% samples2)){
			for (type in c("exon","intron")){

				fileplus <- read.table(paste0(dir, "/", str_split_fixed(sample, "_", n=2)[1], "_PlusTag_", str_split_fixed(sample, "_", n=2)[2], ".merged.mapped.",species,".inputnormalised_",type,".txt"), sep="\t", header=F)
				fileminus <- read.table(paste0(dir, "/", str_split_fixed(sample, "_", n=2)[1], "_MinusTag_", str_split_fixed(sample, "_", n=2)[2], ".merged.mapped.",species,".inputnormalised_",type,".txt"), sep="\t", header=F)
				colnames(fileplus) <- c("name","size","covered_p","sum_p","mean0_p","mean_p", "type")
				colnames(fileminus) <- c("name","size","covered_m","sum_m","mean0_m","mean_m", "type")

				file <- full_join(fileplus, fileminus, by=c("name", "size", "type"))

				if(dim(file)[1] != dim(fileplus)[1]){
					print("ERROR: Files do not have the same number of rows!")
					print(paste0(dir, "/", str_split_fixed(sample, "_", n=2)[1], "_PlusTag_", str_split_fixed(sample, "_", n=2)[2], ".merged.mapped.",species,".inputnormalised_",type,".txt"))
					print(paste0(dir, "/", str_split_fixed(sample, "_", n=2)[1], "_MinusTag_", str_split_fixed(sample, "_", n=2)[2], ".merged.mapped.",species,".inputnormalised_",type,".txt"))
					print(dim(file))
					print(dim(fileplus))
					print(dim(fileminus))
					stop("ERROR: Files do not have the same number of rows!")
				}

				#- make a scatterplot 
				png(paste0(dir, "_plots/scatterplots/scatterplot_", sample, "_",type,".png"))
				range_vals <- range(c(file$mean_m , file$mean_p))
				p <- ggplot(file, aes(x=mean_m, y=mean_p)) +
					geom_point() +
					xlim (range_vals) +
					ylim (range_vals) + 
 					geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "red") +  # Add diagonal line
					ggtitle(paste0(sample, " ", type)) +
					xlab("MinusTag") +
					ylab("PlusTag") +
					theme_minimal() 
				print(p)
				dev.off()
 
			}
		}
		#- sample_plus versus GFP_plus
		#– Do this for any plus sample that is not GFP
		if ((str_split_fixed(sample, "_", n=2)[1] != "GFP") && 
			(paste0(sample, "-", "PlusTag") %in% samples2) &&
			(paste0("GFP_", str_split_fixed(sample, "_", n=2)[2], "-PlusTag") %in% samples2)){
			for (type in c("exon","intron")){

				fileplus <- read.table(paste0(dir, "/", str_split_fixed(sample, "_", n=2)[1], "_PlusTag_", str_split_fixed(sample, "_", n=2)[2], ".merged.mapped.",species,".inputnormalised_",type,".txt"), sep="\t", header=F)
				gfpplus <- read.table(paste0(dir, "/GFP", "_PlusTag_", str_split_fixed(sample, "_", n=2)[2], ".merged.mapped.",species,".inputnormalised_",type,".txt"), sep="\t", header=F)
				colnames(fileplus) <- c("name","size","covered_p","sum_p","mean0_p","mean_p", "type")
				colnames(gfpplus) <- c("name","size","covered_m","sum_m","mean0_m","mean_m", "type")

				file <- full_join(fileplus, gfpplus, by=c("name", "size", "type"))

				if(dim(file)[1] != dim(fileplus)[1]){
					print("ERROR: Files do not have the same number of rows!")
					print(paste0(dir, "/", str_split_fixed(sample, "_", n=2)[1], "_PlusTag_", str_split_fixed(sample, "_", n=2)[2], ".merged.mapped.",species,".inputnormalised_",type,".txt"))
					print(paste0(dir, "/GFP", "_PlusTag_", str_split_fixed(sample, "_", n=2)[2], ".merged.mapped.",species,".inputnormalised_",type,".txt"))
					print(dim(file))
					print(dim(fileplus))
					print(dim(gfpplus))
					stop("ERROR: Files do not have the same number of rows!")
				}

				#- make a scatterplot
				png(paste0(dir, "_plots/scatterplots_GFP/scatterplot_", sample, "_",type,"_vsGFP.png"))
				range_vals <- range(c(file$mean_m, file$mean_p))
				p <- ggplot(file, aes(x=mean_m, y=mean_p)) +
					geom_point() +
					xlim (range_vals) +
					ylim (range_vals) + 
 					geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "red") +  # Add diagonal line
					ggtitle(paste0(sample, " ", type)) +
					xlab("GFP PlusTag") +
					ylab(paste0(str_split_fixed(sample, "_", n=2)[1], " PlusTag")) +
					theme_minimal() 
				print(p)
				dev.off()
 
			}
		}
	}
}


