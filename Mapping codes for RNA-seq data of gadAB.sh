#!/bin/bash
#SBATCH -J WT
#SBATCH -p cpu
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=8
#SBATCH --mem-per-cpu=4G
#SBATCH -o WT.out


#mkdir -p 01_clean 02_bam 03_counts
#index= ~/Reference/Escherichia_coli_str_k_12_substr_mg1655_gca_000005845.ASM584v2.dna.toplevel.fa.gz
3gtf= ~/Reference/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.46.gtf

#for sample in Tn-1_R2.part_001 Tn-1_R2.part_002 Tn-1_R2.part_003; do
for sample in Tn-1.1 Tn-1.2 Tn-1.3 Tn-2 Tn-3 Tn-4 Tn-5 Tn-6 Tn-7 Tn-8 Tn-9 Tn-10; do
    # 1. 质控
    #fastp -i in.R1.fastq.gz -I in.R2.fastq.gz -o out.R1.fastq.gz -O out.R2.fastq.gz    
    #fastp -i 01.RawData/${sample}/${sample}_1.fq.gz -I 01.RawData/${sample}/${sample}_2.fq.gz \
     #     -o 01_clean/${sample}_R1.fq.gz -O 01_clean/${sample}_R2.fq.gz \
      #    -g -q 5 -u 50 -n 15 -l 140 \  # -l 150---140!!!!!
       #   --overlap_diff_limit 1 --overlap_diff_percent_limit 10 \
        #  --detect_adapter_for_pe \
         # -j 01_clean/${sample}.json -h 01_clean/${sample}.html -w 8

    # 2. 比对 + 排序
    zcat 01_clean/${sample}_R2.fq.gz | seqkit grep -s -i -p acaggttggatgataagtccccggtctctagaa > ${sample}.1.fastq
    cutadapt -u -120 -o ${sample}.2.fastq ${sample}.1.fastq
    cutadapt -a acaggttggatgataa -o ${sample}.3.fastq ${sample}.2.fastq
    cat ${sample}.3.fastq | seqkit grep -s -r -i -p TA$ > ${sample}.4.fastq
    cat ${sample}.4.fastq | seqkit seq -m 14 -g > ${sample}.result.fastq
#rm ${i}.1.fastq ${i}.2.fastq ${i}.3.fastq ${i}.4.fastq
    bwa aln ~/Reference/ST-GeneBack-GCA_000022165.1_ASM2216v1_genomic.fna.gz ${sample}.result.fastq > ${sample}.sai
    bwa samse -f ${sample}.sam ~/Reference/ST-GeneBack-GCA_000022165.1_ASM2216v1_genomic.fna.gz ${sample}.sai ${sample}.result.fastq
    samtools sort -@ 8 -o 02_bam/${sample}.sorted.bam ${sample}.sam
    samtools index 02_bam/${sample}.sorted.bam

    #samtools view -Sb ${i}.sam > ${i}.bam
    #featureCounts -T 6 -t exon -g gene_id -a /home/yanxiaod/Reference/Escherichia_coli_str_k_12_substr_mg1655.ASM584v2.46.gtf -o ${i}.txt  ${i}.bam
    rm ${sample}.sai #${i}.result.fastq ${i}.txt.summary

    gawk '(and(16, $2))' ${sample}.sam > ${sample}.reverse_mapped_reads.sam
    gawk '(!and(16, $2))' ${sample}.sam > ${sample}.forward_mapped_reads.sam
    awk '{print $4}' ${sample}.reverse_mapped_reads.sam > ${sample}.reverse.start.txt
    awk '{print $4}' ${sample}.forward_mapped_reads.sam > ${sample}.forward.start.txt

    sort ${sample}.forward.start.txt | uniq -c > ${sample}f.txt
    sort ${sample}.reverse.start.txt | uniq -c > ${sample}r.txt
done

# 3. 一次性定量
featureCounts -T 8  -t CDS -g gene_id -a ~/Reference/ST-GeneBack-GCA_000022165.1_ASM2216v1_genomic.gtf.gz -o 03_counts/all_counts-ST.txt 02_bam/*.sorted.bam

