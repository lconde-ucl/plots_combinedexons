library(dplyr)
library(ggplot2)
library(stringr)


for (dir in c("human/bigwigs","mouse/bigwigs")){
	
	dir.create(paste0(dir, "_plots/scatterplots_logs/"), showWarnings = FALSE)

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
			
			type = "exonintron"

			#-exon
			fileplusexon <- read.table(paste0(dir, "/", str_split_fixed(sample, "_", n=2)[1], "_PlusTag_", str_split_fixed(sample, "_", n=2)[2], ".merged.mapped.",species,".inputnormalised_exon.txt"), sep="\t", header=F)
			fileminusexon <- read.table(paste0(dir, "/", str_split_fixed(sample, "_", n=2)[1], "_MinusTag_", str_split_fixed(sample, "_", n=2)[2], ".merged.mapped.",species,".inputnormalised_exon.txt"), sep="\t", header=F)
			colnames(fileplusexon) <- c("name","size","covered_p","sum_p","mean0_p","mean_p", "type")
			colnames(fileminusexon) <- c("name","size","covered_m","sum_m","mean0_m","mean_m", "type")

			#-intron
			fileplusintron <- read.table(paste0(dir, "/", str_split_fixed(sample, "_", n=2)[1], "_PlusTag_", str_split_fixed(sample, "_", n=2)[2], ".merged.mapped.",species,".inputnormalised_intron.txt"), sep="\t", header=F)
			fileminusintron <- read.table(paste0(dir, "/", str_split_fixed(sample, "_", n=2)[1], "_MinusTag_", str_split_fixed(sample, "_", n=2)[2], ".merged.mapped.",species,".inputnormalised_intron.txt"), sep="\t", header=F)
			colnames(fileplusintron) <- c("name","size","covered_p","sum_p","mean0_p","mean_p", "type")
			colnames(fileminusintron) <- c("name","size","covered_m","sum_m","mean0_m","mean_m", "type")

			fileexon <- full_join(fileplusexon, fileminusexon, by=c("name", "size", "type"))
			fileintron <- full_join(fileplusintron, fileminusintron, by=c("name", "size", "type"))

			if(dim(fileexon)[1] != dim(fileplusexon)[1]){
				print("ERROR: Files do not have the same number of rows!")
				print(paste0(dir, "/", str_split_fixed(sample, "_", n=2)[1], "_PlusTag_", str_split_fixed(sample, "_", n=2)[2], ".merged.mapped.",species,".inputnormalised_exon.txt"))
				print(paste0(dir, "/", str_split_fixed(sample, "_", n=2)[1], "_MinusTag_", str_split_fixed(sample, "_", n=2)[2], ".merged.mapped.",species,".inputnormalised_exon.txt"))
				print(dim(fileexon))
				print(dim(fileplusexon))
				print(dim(fileminusexon))
				stop("ERROR: Files do not have the same number of rows!")
			}
			if(dim(fileintron)[1] != dim(fileplusintron)[1]){
				print("ERROR: Files do not have the same number of rows!")
				print(paste0(dir, "/", str_split_fixed(sample, "_", n=2)[1], "_PlusTag_", str_split_fixed(sample, "_", n=2)[2], ".merged.mapped.",species,".inputnormalised_intron.txt"))
				print(paste0(dir, "/", str_split_fixed(sample, "_", n=2)[1], "_MinusTag_", str_split_fixed(sample, "_", n=2)[2], ".merged.mapped.",species,".inputnormalised_intron.txt"))
				print(dim(fileintron))
				print(dim(fileplusintron))
				print(dim(fileminusintron))
				stop("ERROR: Files do not have the same number of rows!")
			}

			file <- rbind(fileexon, fileintron)

			#- make a scatterplot of the log
			png(paste0(dir, "_plots/scatterplots_logs/scatterplot_", sample, "_",type,".png"))
			#	geom_bin_2d(bins = 50) +
			p <- ggplot(file, aes(x = mean_m, y = mean_p)) +
			    geom_point() +
			    scale_x_continuous(trans = scales::pseudo_log_trans(base = 10)) +
			    scale_y_continuous(trans = scales::pseudo_log_trans(base = 10)) +
			    geom_abline(intercept = 0, slope = 1, linetype = "dotted", color = "red") +  # Add diagonal line
			    ggtitle(paste0(sample, " ", type)) +
			    xlab("-tag (pseudo-log scale)") +
			    ylab(paste0("+tag (pseudo-log scale)")) +
			    theme_minimal(base_family = "Arial") +   # set font family
			    theme(
			        axis.title = element_text(size = 32, family = "Arial"),   # axis labels
			        axis.text  = element_text(size = 30, family = "Arial")    # axis tick values
			    )
			print(p)
			dev.off()

		}
	}
}


