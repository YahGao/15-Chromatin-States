#=================#
# Epi Genome data #
#=================#
# Step 1. Convert aligned reads to bed format
bedtools bamtobed -i example.bam > $example.bed

# Step 2. Create binned and binarized tracks
java -mx20000M -jar ChromHMM.jar BinarizeBed CHROMSIZES/bosTau8.txt bedfiles_folder/ cellmarkfiletable output_folder1/

# Step 3. Train the model and segment the genome
java -mx20000M -jar ChromHMM.jar LearnModel -p 20 output_folder1/ output_folder2/ 15 bosTau8



#==============#
# RNA-Seq data #
#==============#
# Step 1. Build index
STAR --runMode genomeGenerate \
--genomeDir . \
--genomeFastaFiles Bos_taurus.UMD3.1.dna.toplevel.fa \
--sjdbGTFfile Bos_taurus.UMD3.1.94.gtf \
--runThreadN 5

# Step 2. Quality control
java -jar trimmomatic_folder/trimmomatic-0.39.jar PE \
-phred33 sample_1.fq.gz sample_2.fq.gz sample_1.clean.fq.gz sample_1_unpaired.fq.gz sample_2.clean.fq.gz sample_2_unpaired.fq.gz \
-threads 5 ILLUMINACLIP:trimmomatic_folder/adapters/TruSeq3-PE.fa:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36

# Step 3. Mapping to genome
STAR --runThreadN 10 \
--genomeDir Bovine_genome \
--sjdbGTFfile Bos_taurus.UMD3.1.94.gtf \
--quantMode TranscriptomeSAM \
--outSAMtype BAM SortedByCoordinate \
--outSAMunmapped Within \
--readFilesCommand zcat \
--outFilterMismatchNmax 3 \
--readFilesIn sample_1.clean.fq.gz sample_2.clean.fq.gz \
--outFileNamePrefix sample-STAR

# Step 4. FPKM quantification
stringtie -p 5 -e -B -G Bos_taurus.UMD3.1.94.gtf -o sample.gtf -A sample.tsv sample-STARAligned.sortedByCoord.out.bam

# Step 5. Run featureCounts
featureCounts -T 10 -p -t exon -g gene_id -a Bos_taurus.UMD3.1.94.gtf -o sample.featureCounts.txt sample-STARAligned.sortedByCoord.out.bam

# Step 6. Merge gtf files
stringtie --merge -p 8 -G Bos_taurus.UMD3.1.94.gtf -o all_merge.gtf all_merge_gtf.txt
# all_merge_gtf.txt includes all .gtf files generated from last step.

# Step 7. DEG analysis
cuffdiff -o all_diff -b Bos_taurus.UMD3.1.dna.toplevel.fa -p 20 -L BW,AW -u all_merge.gtf BW1-STARAligned.sortedByCoord.out.bam,BW2-STARAligned.sortedByCoord.out.bam,BW3-STARAligned.sortedByCoord.out.bam AW1-STARAligned.sortedByCoord.out.bam,AW2-STARAligned.sortedByCoord.out.bam,AW3-STARAligned.sortedByCoord.out.bam



#==============================#
# WGBS data (BW as an example) #
#==============================#
# Step 1.
bismark_folder/bismark_genome_preparation --bowtie2 genome_folder/

# Step 2.
cd fq_folder/
cat *_1.fq.gz > BW_1.fq.gz
cat *_2.fq.gz > BW_2.fq.gz

trim_galore --paired --fastqc --max_n 15 -o trim_reads BW_1.fq.gz BW_2.fq.gz  

bismark --multicore 2 --bowtie2 --gzip -p 4 -N 0 -o bamfile_folder/ genome_folder/ -1 BW_1_val_1.fq.gz -2 BW_2_val_2.fq.gz

cd bamfile_folder/ 
deduplicate_bismark -p --bam BW_1_val_1_bismark_bt2_pe.bam
bismark_methylation_extractor -p --gzip --ignore_r2 6 --multicore 8 --bedgraph -o methylation_folder/ --cytosine_report --genome_folder genome_folder/ *.deduplicated.bam

# Step 3.
cd methylation_folder

awk -F "[_.\t]" '{print $i"_"$2}' BW_1_val_1_bismark_bt2_pe.deduplicated.CpG_report.txt.CpG_report.txt > BW_1_val_1_bismark_bt2_pe.deduplicated.CpG_report.txt.CpG_report1.txt
awk 'OFS="\t"{if (($4+$5)>0) print $1,$2,$3,"CpG",$4/($4+$5),$4+$5}'  BW_1_val_1_bismark_bt2_pe.deduplicated.CpG_report.txt.CpG_report1.txt | sort -k1,1 -k2,2n > BW.meth
symmetric-cpgs BW.meth > BW.cpg.meth

# Step 4.
sed -i 's/chr//g' BW_15_segments.bed
sort -k1,1 -k2,2n BW_15_segments.bed > BW_15_segments1.bed
roimethstat -o BW_final.meth BW_15_segments1.bed BW.cpg.meth
