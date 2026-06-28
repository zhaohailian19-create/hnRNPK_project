### clean
trim_galore -q 10 --paired --phred33 -e 0.1 --fastqc --length 36 --stringency 3 -o ./ ${id}_R1.fq.gz ${id}_R2.fq.gz

### mapping
bwa mem -t 10 -M GRCh37.p13.genome.clean.fa ${id}_R1_val_1.fq.gz /${id}_R2_val_2.fq.gz |samtools sort -@ 10 -o ${id}.to_hg19.sort.bam
java -jar picard.jar MarkDuplicates I=${id}.to_hg19.sort.bam O=${id}.to_hg19.rmdup.bam M=${id}.to_hg19.md.txt REMOVE_DUPLICATES=true
samtools view -@ 8 -h -q 20 -b ${id}.to_hg19.rmdup.bam > ${id}.to_hg19.uniq.rmdup.bam
samtools index ${id}.to_hg19.uniq.rmdup.bam
bedtools bamtobed -i ${id}.to_hg19.uniq.rmdup.bam > ${id}.bed

### call peaks
macs2 callpeak -t hnRNPK_cuttag_merge.bam -p 1e-20 -f BAMPE -g hs -n out0.hnRNPK_cuttag
macs2 callpeak -t hnRNPK_cuttag_siK_merge.bam -p 1e-20 -f BAMPE -g hs -n out0.hnRNPK_cuttag_siK
macs2 callpeak -t NC_cuttag_merge.bam -p 1e-20 -f BAMPE -g hs -n out0.NC_cuttag


## diffbind by DiffBind
Rscript diffbind.r
