# QC
fastqc <folder_name>/* -o fastqc_results -t 64
multiqc --outdir multiqc_results --filename fastqc_summary.html --title "QC Report" fastqc_results

# To check the ribosomal contamination
bbduk.sh in=<filename_R1.fastq> in2=<filename_R2.fastq> out=<out_filename_R1.fastq> out2=<out_filename_R2.fastq> ref=ref/human_ribosomal_ru.fa threads=40

# To align the reads
STAR --runThreadN 64 --outSAMunmapped Within --twopassMode Basic --sjdbGTFfile ref/gencode.v49.primary_assembly.annotation.gtf --readFilesCommand zcat --genomeDir ref/STAR_index --readFilesIn <filename_R1.fastq.gz> <filename_R2.fastq.gz> --outSAMtype BAM SortedByCoordinate --outFileNamePrefix <out_filename>

# To sort bam files
samtools sort -N <filename_Aligned.sortedByCoord.out.bam> -O bam -T temp_sort -o<out_filename_Aligned.sortedByReadNames.out.bam>

# To use check the quality of alignment
qualimap bamqc -bam <filename_Aligned.sortedByCoord.out.bam> -gff ref/gencode.v49.primary_assembly.annotation.gtf -gd HUMAN --java-mem-size=200G -nt 40 -outdir <dir_name>
qualimap rnaseq -bam <filename_Aligned.sortedByReadNames.out.bam> --sorted --paired -gtf ref/gencode.v49.primary_assembly.annotation.gtf --java-mem-size=400G -outdir <dir_name>

# To get the raw counts
featureCounts -p --countReadPairs -B -P -C -T 64 -a ref/gencode.v49.primary_assembly.annotation.gtf -o <filename.txt> <dir_name>/*ByReadNames*.bam 

# To use get TPM
salmon quant -i ref/salmon_decoyed_index -l A -1 <filename_R1.fastq> -2 <filename_R2.fastq> --gcBias -p 64 --validateMappings --numBootstraps 100 -o <dir_name> 
