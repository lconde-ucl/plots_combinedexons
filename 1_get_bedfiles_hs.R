library(dplyr)
library(stringr)


#- Human genome
canonical_transcript="knownCanonical.txt"
exons_file="wgEncodeGencodeBasicV40.tsv"
#exons_file="exons.txt" #- for testing
genome="human/hs"

#- Mouse genome
#canonical_transcript="knownCanonical.Mm.txt"
#exons_file="wgEncodeGencodeBasicVM23.tsv"
#genome="mouse/mm"


dir.create(paste0(genome, "_bedfiles"), recursive=T)

#- get 1 transcript per gene only (remove non canonical chromsomes)
canonical<-read.table(canonical_transcript, sep="\t", header=F)
colnames(canonical)<-c("chr","start","end","id","transcript_id","gene_id")
dim(canonical)
canonical<-canonical[!str_detect(canonical$chr, "_"),]
canonical <- filter(canonical, chr != "chrM")
dim(canonical)


#- get exons in anonical transcripts only
#- some transcripts appear in both chrX and chrY in the exons file (eg  ENST00000627721.3)
#- so search for transcript-chr pairs
transcripts<-paste0(canonical$transcript_id, "-", canonical$chr)
exons<-read.table(exons_file, header=T, sep="\t")
dim(exons)
exons<-filter(exons, paste0(name, "-",chrom)  %in% transcripts)
dim(exons)
dim(unique(exons))


exon <- data.frame()
intron <- data.frame()

for (i in c(1:length(exons$name))){

  d<-exons[i,]

  start_positions <- as.numeric(unlist(strsplit(d$exonStarts, ",")))
  end_positions <- as.numeric(unlist(strsplit(d$exonEnds, ",")))

  for (j in c(1:length(start_positions))){
    exon <- rbind(exon, cbind(d$chrom, sprintf("%.0f", start_positions[j]), sprintf("%.0f", end_positions[j]), paste0(d$name, "-", j)))
  }
  if(d$exonCount > 1) {
    for (k in c(1:(length(start_positions)-1))){
      intron <- rbind(intron, cbind(d$chrom, sprintf("%.0f", end_positions[k]+1), sprintf("%.0f", start_positions[k+1]-1), paste0(d$name, "-", k)))
    }
  }

}

write.table(exon, paste0(genome, "_bedfiles/exon.bed"), sep="\t", row.names=F, col.names=F, quote=F)
write.table(intron, paste0(genome, "_bedfiles/intron.bed"), sep="\t", row.names=F, col.names=F, quote=F)
