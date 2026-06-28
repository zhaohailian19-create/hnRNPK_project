### hnRNPK CUTTag and pol2 ChIP-seq data pre-process
### clean
trim_galore -q 10 --paired --phred33 -e 0.1 --fastqc --length 36 --stringency 3 -o ./ ${id}_R1.fq.gz ${id}_R2.fq.gz

### mapping
bwa mem -t 10 -M GRCh37.p13.genome.clean.fa ${id}_R1_val_1.fq.gz ${id}_R2_val_2.fq.gz |samtools sort -@ 10 -o ${id}.to_hg19.sort.bam
java -jar picard.jar MarkDuplicates I=${id}.to_hg19.sort.bam O=${id}.to_hg19.rmdup.bam M=${id}.to_hg19.md.txt REMOVE_DUPLICATES=true
samtools view -@ 8 -h -q 20 -b ${id}.to_hg19.rmdup.bam > ${id}.to_hg19.uniq.rmdup.bam
samtools index ${id}.to_hg19.uniq.rmdup.bam
bedtools bamtobed -i ${id}.to_hg19.uniq.rmdup.bam > ${id}.bed

### call peaks
macs2 callpeak -t hnRNPK_cuttag_merge.bam -p 1e-20 -f BAMPE -g hs -n out0.hnRNPK_cuttag
macs2 callpeak -t hnRNPK_cuttag_siK_merge.bam -p 1e-20 -f BAMPE -g hs -n out0.hnRNPK_cuttag_siK
macs2 callpeak -t NC_cuttag_merge.bam -p 1e-20 -f BAMPE -g hs -n out0.NC_cuttag

### identify high-confidence hnRNPK binding peaks compared with control (siK and NC group)
perl calculate_enrichRatio.pl hnRNPK_cuttag_peaks.narrowPeak hnRNPK_cuttag_merge.bed hnRNPK_cuttag_siK_merge.bed > hnRNPK_cuttag.with_siK.count.bed
awk -F '\t' 'BEGIN{OFS="\t"}{print $0,($11+1)/($12+1)}' hnRNPK_cuttag.with_siK.count.bed > hnRNPK_cuttag.with_siK.fold.bed
awk -F '\t' 'BEGIN{OFS="\t"}{if(int($13)>=int(5)) print $0 }' hnRNPK_cuttag.with_siK.fold.bed > hnRNPK_cuttag.with_siK.fold5.bed
perl calculate_enrichRatio.pl hnRNPK_cuttag.with_siK.fold5.bed hnRNPK_cuttag_merge.bed NC_cuttag_merge.bed > hnRNPK_cuttag.with_siK.fold5.with_NC.count.bed
awk -F '\t' 'BEGIN{OFS="\t"}{print $0,($14+1)/($15+1)}' hnRNPK_cuttag.with_siK.fold5.with_NC.count.bed > hnRNPK_cuttag.with_siK.fold5.with_NC.fold.bed
awk -F '\t' 'BEGIN{OFS="\t"}{if(int($16)>=int("'5'")) print $0 }' hnRNPK_cuttag.with_siK.fold5.with_NC.fold.bed > hnRNPK_cuttag.with_siK.fold5.with_NC.fold5.bed
intersectBed -v -wa -a hnRNPK_cuttag.with_siK.fold5.with_NC.fold5.bed -b hnRNPK_cuttag_siK_peaks.narrowPeak > tmp.bed
intersectBed -v -wa -a tmp.bed -b NC_cuttag_peaks.narrowPeak > hnRNPK_cuttag.with_siK.fold5.with_NC.fold5.rm_siK_NC_peaks.bed

### identify pol2 ChIP-seq peaks 
macs2 callpeak -t pol2_ChIP_seq_siNC_merge.bam -c pol2_ChIP_seq_siNC-Input_merge.bam -p 0.01 -f BAMPE -g hs -n pol2ChIP_siNC
macs2 callpeak -t pol2_ChIP_seq_siK_merge.bam -c pol2_ChIP_seq_siK-Input_merge.bam -p 0.01 -f BAMPE -g hs -n pol2ChIP_siK
 
### identify the different peaks between siNC and siK
Rscript diffbind.r
awk -F '\t' 'BEGIN{OFS="\t"}{if(($2!="seqnames")&&($10<(-0.584963))&&($11<0.05)) print $2,$3,$4,$10,$11 }' out1.deseq2.txt > pol2_siNC_vs_siK_diffbind_deseq2.up.bed
awk -F '\t' 'BEGIN{OFS="\t"}{if(($2!="seqnames")&&($10>0.584963)&&($11<0.05)) print $2,$3,$4,$10,$11 }' out1.deseq2.txt > pol2_siNC_vs_siK_diffbind_deseq2.down.bed
intersectBed -wa -F 0.1 -a gencode_TSS_2Kb.v19.bed -b pol2_siNC_vs_siK_diffbind_deseq2.up.bed |sort |uniq > pol2_siNC_vs_siK_diffbind_deseq2.up.genes.bed
intersectBed -wa -F 0.1 -a gencode_TSS_2Kb.v19.bed -b pol2_siNC_vs_siK_diffbind_deseq2.down.bed |sort |uniq > pol2_siNC_vs_siK_diffbind_deseq2.down.genes.bed

### HiChIP data process
trim_galore -q 10 --paired --phred33 -e 0.1 --fastqc --length 36 --stringency 3 -o ./ ${id}_R1.fq.gz ${id}_R2.fq.gz
HiC-Pro --input input --output output --conf ./config-hicpro.txt

