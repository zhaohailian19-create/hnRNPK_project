library(DiffBind)
samples= read.table("tmp.txt",header=T,row.names=1)
tamoxifen=dba(sampleSheet=samples)
tamoxifen

pdf("1.pdf")
plot(tamoxifen)
dev.off()

#tamoxifen <- dba.contrast(tamoxifen, reorderMeta=list(Condition="LM2")


tamoxifen1 <- dba.count(tamoxifen, bUseSummarizeOverlaps=TRUE)
pdf("2.pdf")
dba.plotPCA(tamoxifen1,  attributes=DBA_FACTOR, label=DBA_ID)
dev.off()


tamoxifen2 <- dba.contrast(tamoxifen1, categories=DBA_FACTOR,minMembers = 2)
tamoxifen2
tamoxifen3 <- dba.analyze(tamoxifen2, method=DBA_ALL_METHODS)
tamoxifen3

dba.show(tamoxifen3, bContrasts=T)
pdf("3.pdf")
dba.plotVenn(tamoxifen3,contrast=1,method=DBA_ALL_METHODS)
dev.off()

comp1.edgeR <- dba.report(tamoxifen3, method=DBA_EDGER, contrast = 1, th=1)
comp1.deseq <- dba.report(tamoxifen3, method=DBA_DESEQ2, contrast = 1, th=1)

# EdgeR
out <- as.data.frame(comp1.edgeR)
write.table(out, file="out1.edgeR.txt", sep="\t", quote=F, col.names = NA)
# DESeq2
out <- as.data.frame(comp1.deseq)
write.table(out, file="out1.deseq2.txt", sep="\t", quote=F, col.names = NA)


