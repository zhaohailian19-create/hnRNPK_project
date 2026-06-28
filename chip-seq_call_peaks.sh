### call peaks
macs2 callpeak -t hela_953-Pol2ChIP_merge.bam -c hela_953-Pol2ChIP-Input_merge.bam -p 0.01 -f BAMPE -g hs -n out0.hela_953-Pol2ChIP_WithInput_p0.01
macs2 callpeak -t hela_WT-Pol2ChIP_merge.bam -c hela_WT-Pol2ChIP-input_merge.bam -p 0.01 -f BAMPE -g hs -n out0.hela_WT-Pol2ChIP_WithInput_p0.01
macs2 callpeak -t hela_siNC-Pol2ChIP_merge.bam -c hela_siNC-Pol2ChIP-Input_merge.bam -p 0.01 -f BAMPE -g hs -n out0.hela_siNC-Pol2ChIP_p0.01
macs2 callpeak -t hela_siK-Pol2ChIP_merge.bam -c hela_siK-Pol2ChIP-Input_merge.bam -p 0.01 -f BAMPE -g hs -n out0.hela_siK-Pol2ChIP_p0.01
macs2 callpeak -t ttf_953-Pol2ChIP_merge.bam -c ttf_953-Pol2ChIP-Input_merge.bam -p 0.01 -f BAMPE -g mm -n out0.ttf_953-Pol2ChIP_p0.01
macs2 callpeak -t ttf_WT-Pol2ChIP_merge.bam -c ttf_WT-Pol2ChIP-Input_merge.bam -p 0.01 -f BAMPE -g mm -n out0.ttf_WT-Pol2ChIP_p0.01

## diffbind by DiffBind
Rscript diffbind.r
