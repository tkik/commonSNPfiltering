library(RnBeads)
library(dplyr)

library(RnBeads.hg19)
anno<-rnb.get.annotation(type = "probesEPIC", assembly = "hg19")

chromosomes<-paste("chr", c(1:22, "X", "Y"), sep="")
#in.folder <- "E:/Reka/WES/database/"
in.folder <- "P:/old_dkfz_project_codes/00-random_scripts/common_snp_filtering/results_hg19_155/"
out.folder <- "P:/old_dkfz_project_codes/00-random_scripts/common_snp_filtering/results_hg19_155/filters"
#out.folder <- "Z:/Reka/00-random_scripts/common_snp_filtering/results/"
all_overlap<-list()
all_overlap_3<-list()

for (chr in chromosomes) {

  filename <- file.path(in.folder, paste0("common_all_chr_", chr, ".vcf"))
  common_all<- read.delim(filename, header=FALSE, stringsAsFactors = F)
  colnames(common_all)[1:3]<-c("chr", "start", "rsid")
  common_all$chr[!grepl("chr", common_all$chr)]<-paste0("chr", common_all$chr[!grepl("chr", common_all$chr)])

  common_all <- common_all %>%
    mutate(end=start+nchar(V4)) %>%
    mutate(dbGaP_PopFreq=ifelse(grepl("dbGaP_PopFreq", V8), gsub("(RS=.*)(dbGaP_PopFreq:)([[:digit:]]*\\.*[[:digit:]]*),(.*)",  "\\3", V8), NA)) %>%
    mutate(thousandGenomes=ifelse(grepl("1000Genomes", V8), gsub("(RS=.*)(1000Genomes:)([[:digit:]]*\\.[[:digit:]]*),(.*)",  "\\3", V8), NA)) %>%
    mutate(GnomAD=ifelse(grepl("GnomAD:", V8), gsub("(RS=.*)(GnomAD:)([[:digit:]]*\\.[[:digit:]]*),(.*)",  "\\3", V8), NA)) %>%
    filter(dbGaP_PopFreq<=0.99 | thousandGenomes<=0.99 | GnomAD<=0.99) %>%
    makeGRangesFromDataFrame( keep.extra.columns = T, ignore.strand = T)



anno_chr<-anno[[chr]]
start(anno_chr[strand(anno_chr)=="-"])<-start(anno_chr[strand(anno_chr)=="-"])+1
end(anno_chr[strand(anno_chr)=="+"])<-end(anno_chr[strand(anno_chr)=="+"])-1

overlapping <- subsetByOverlaps(anno_chr, common_all, type = "any")
#snps<- subsetByOverlaps(common_all, anno_chr, type = "start")
overlaps<-findOverlaps(overlapping, common_all, type = "any", select="all")
if (length(overlapping)!=0){
mcols(overlapping)$SNP<-NA
for ( i in 1:overlaps@nLnode) {
  mcols(overlapping)$SNP[i]<-paste(mcols(common_all)$rsid[overlaps@to[overlaps@from==i]], collapse = ", ")
  mcols(overlapping)$dbGaP_PopFreq[i]<-paste(mcols(common_all)$dbGaP_PopFreq[overlaps@to[overlaps@from==i]], collapse = ", ")
  mcols(overlapping)$thousandGenomes[i]<-paste(mcols(common_all)$thousandGenomes[overlaps@to[overlaps@from==i]], collapse = ", ")
  mcols(overlapping)$GnomAD[i]<-paste(mcols(common_all)$GnomAD[overlaps@to[overlaps@from==i]], collapse = ", ")
}} else {
  overlapping<-NA
}
all_overlap<-c(all_overlap, overlapping)
start(anno_chr)<-start(anno_chr)-3
end(anno_chr)<-end(anno_chr)+3
rm(overlapping, overlaps)
overlapping <- subsetByOverlaps(anno_chr, common_all, type = "any")
#snps<- subsetByOverlaps(common_all, anno_chr, type = "start")
overlaps<-findOverlaps(overlapping, common_all, type = "any", select="all")
if (length(overlapping)!=0){
mcols(overlapping)$SNP_3<-NA
for ( i in 1:overlaps@nLnode) {
  mcols(overlapping)$SNP_3[i]<-paste(mcols(common_all)$rsid[overlaps@to[overlaps@from==i]], collapse = ", ")
  mcols(overlapping)$dbGaP_PopFreq[i]<-paste(mcols(common_all)$dbGaP_PopFreq[overlaps@to[overlaps@from==i]], collapse = ", ")
  mcols(overlapping)$thousandGenomes[i]<-paste(mcols(common_all)$thousandGenomes[overlaps@to[overlaps@from==i]], collapse = ", ")
  mcols(overlapping)$GnomAD[i]<-paste(mcols(common_all)$GnomAD[overlaps@to[overlaps@from==i]], collapse = ", ")
}} else {
  overlapping<-NA

}
all_overlap_3<-c(all_overlap_3, overlapping)
rm(overlapping, overlaps, anno_chr, common_all)
}

names(all_overlap)<-chromosomes
names(all_overlap_3)<-chromosomes

combined_all_overlap<-GRanges()
combined_all_overlap_3<-GRanges()

for (i in 1:24){
  if (!is.null(all_overlap[[i]]))
  combined_all_overlap_3<-append(combined_all_overlap_3, all_overlap_3[[i]])
}
combined_all_overlap_3<-as.data.frame(combined_all_overlap_3)

for (i in 1:24){
  if (!all(is.na(all_overlap[[i]])))
  combined_all_overlap<-append(combined_all_overlap, all_overlap[[i]])
}
combined_all_overlap<-as.data.frame(combined_all_overlap)

saveRDS(combined_all_overlap, file = file.path(out.folder, "cg_o_SNP_155_EPIC.RDS"))

saveRDS(combined_all_overlap_3, file = file.path(out.folder, "cg_o3_SNP_155_EPIC.RDS"))

cgs<-data.frame(rownames(combined_all_overlap))
write.table(cgs, file = file.path(out.folder, "blacklist_overlapping_SNP_155_EPIC.txt"), sep="\t", row.names=F, col.names=F, quote=F)

cgs_3<-data.frame(rownames(combined_all_overlap_3))
write.table(cgs_3, file = file.path(out.folder, "blacklist_overlapping_SNP_3_155_EPIC.txt"), sep="\t", row.names=F, col.names=F, quote=F)



################some additional examination

for (i in c(1:22, "X", "Y")){
  anno_3<-anno
  anno_3<-anno[[paste0("chr", i)]][which(anno[[paste0("chr", i)]]$"SNPs 3"!=0),]
  print(sum(!(names(anno_3) %in% rownames(combined_all_overlap_3[combined_all_overlap_3$seqnames==paste0("chr", i),]) )))
 print(sum(!(rownames(combined_all_overlap_3[combined_all_overlap_3$seqnames==paste0("chr", i),]) %in%  names(anno_3))))
}

