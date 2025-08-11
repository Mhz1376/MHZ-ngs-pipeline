#!/bin/bash



# Function to process each SRR identifier
process_srr() {
  SRR=$1
  # Define common directories
  
  BASE_DIR="/path/to/base/directory/SRA/$SRR"
  RESULTS_DIR="/path/to/results/directory/$SRR"
  SRA_DATA_DIR="/path/to/sra/data/directory/sra/$SRR.sra"
  
  echo "Processing $SRR..."

  # Create individual SRR directory and base subdirectories
  mkdir -p $BASE_DIR/QC $BASE_DIR/ALIGNMENT $BASE_DIR/POST_ALIGNMENT $BASE_DIR/BQSR $BASE_DIR/VCF $BASE_DIR/VQSR $BASE_DIR/ANNOTATION

  # Define log files
  LOG_FILE=$BASE_DIR/ANNOTATION/${SRR}_process.log
  ERR_FILE=$BASE_DIR/ANNOTATION/${SRR}_process.err

  {
    echo "Start processing $SRR at $(date)"

    # Copy SRR files
    echo "Copying SRR files..."
    cp $SRA_DATA_DIR $BASE_DIR/QC

    # Navigate to QC directory
    cd $BASE_DIR/QC

    # Perform fastq dump
    echo "Running fasterq-dump..."
    fasterq-dump $SRR.sra --split-files --skip-technical
    
    # Remove intermediate file
    echo "Cleaning up intermediate files..."
    rm $SRR.sra

	# Step 1: Run fastp for initial processing
	echo "Running fastp..."
	fastp --in1 ${SRR}_1.fastq --in2 ${SRR}_2.fastq --out1 ${SRR}_1.trimmed.fastq --out2 ${SRR}_2.trimmed.fastq --adapter_sequence AGATCGGAAGAGC --adapter_sequence_r2 AGATCGGAAGAGC --cut_window_size 4 --cut_mean_quality 20 --length_required 36 -g -h ${SRR}.html -t 12 &> ${SRR}.log

	# Step 2: Run BBDuk for further cleaning and adapter removal
	echo "Running BBDuk..."
	bbduk.sh in1=${SRR}_1.trimmed.fastq in2=${SRR}_2.trimmed.fastq out1=${SRR}_1.bbduk.cleaned.fastq out2=${SRR}_2.bbduk.cleaned.fastq ref=/path/to/BBTools/adapters.fa ktrim=r k=23 mink=11 hdist=1 qtrim=rl trimq=10 stats=${SRR}_bbduk.stats threads=12

	# Step 3: Quality control with IlluQC (optional)
	echo "Running IlluQC..."
	IlluQC.pl -pe ${SRR}_1.bbduk.cleaned.fastq ${SRR}_2.bbduk.cleaned.fastq N A -l 70 -s 20 -p 12 -t 1 -o . -z g

	# Step 4: Run FastQC on the cleaned fastq files
	echo "Running FastQC..."
	fastqc -t 12 ${SRR}_1.bbduk.cleaned.fastq_filtered.gz ${SRR}_2.bbduk.cleaned.fastq_filtered.gz -o .

	# Optional: Summarize results with MultiQC
	echo "Running MultiQC..."
	multiqc .

    # Clean up intermediate files (optional)
    echo "Cleaning up intermediate files..."
    rm ${SRR}_1.fastq ${SRR}_2.fastq ${SRR}_1.trimmed.fastq ${SRR}_2.trimmed.fastq ${SRR}_1.bbduk.cleaned.fastq ${SRR}_2.bbduk.cleaned.fastq


    # Navigate to ALIGNMENT directory
    cd $BASE_DIR/ALIGNMENT

    # Perform alignment with BWA
    echo "Running BWA MEM..."
    bwa mem -t 12 -M /path/to/reference/genome/human_g1k_v37.fasta ../QC/${SRR}_1.bbduk.cleaned.fastq_filtered.gz ../QC/${SRR}_2.bbduk.cleaned.fastq_filtered.gz > ALIGNMENT.sam
    
    # Navigate back to QC directory to clean up
    cd $BASE_DIR/QC

    # Remove intermediate files
    echo "Cleaning up intermediate files..."
    rm  ${SRR}_1.bbduk.cleaned.fastq_filtered.gz ${SRR}_2.bbduk.cleaned.fastq_filtered.gz
    
    # Navigate to ALIGNMENT directory
    cd $BASE_DIR/ALIGNMENT
    
    # Convert SAM to BAM
    echo "Running Samtools view..."
    samtools view -@ 12 -b ALIGNMENT.sam > ALIGNMENT2.bam


    # Navigate to POST_ALIGNMENT directory
    cd $BASE_DIR/POST_ALIGNMENT

    # Sort BAM file
    echo "Running Picard SortSam..."
    picard SortSam -I ../ALIGNMENT/ALIGNMENT2.bam -O POST_ALIGNMENT_Sorted.bam -SO coordinate

	# Define the input BAM file and output files
	echo "Calculating the depth..."
	bam_file="POST_ALIGNMENT_Sorted.bam"
	depth_file="file_name.depth.txt"
	average_depth_file="$BASE_DIR/QC/average_depth.txt"
	tmp_dir="tmp_depths"

	# Create a temporary directory for parallel processing
	mkdir -p "$tmp_dir"

	# Calculate depth in parallel using samtools
	samtools depth "$bam_file" | split -l 1000000 - "$tmp_dir/depth_chunk_"

	# Process each chunk to calculate average depth
	echo "Processing depth files..."
	find "$tmp_dir" -name "depth_chunk_*" | xargs -P 12 -I {} awk '{sum+=$3; count++} END {if (count > 0) print sum/count; else print "No data"}' {} > "$average_depth_file"

	# Clean up temporary files
	rm -r "$tmp_dir"

	echo "Average depth of coverage has been saved to $average_depth_file"
    
    # Mark duplicates
    echo "Running Picard MarkDuplicates..."
    picard MarkDuplicates -I POST_ALIGNMENT_Sorted.bam -O POST_ALIGNMENT_Sorted_remove-dup.bam -M POST_ALIGNMENT_Sorted_remove-dup.txt -REMOVE_DUPLICATES true -AS true

    # Add or replace read groups
    echo "Running Picard AddOrReplaceReadGroups..."
    picard AddOrReplaceReadGroups -I POST_ALIGNMENT_Sorted_remove-dup.bam -O POST_ALIGNMENT_Sorted_remove-dup_RG.bam -LB PRJNA824495 -PL HiSeq2500 -SM $SRR -PU $SRR

    # Build BAM index
    echo "Running Picard BuildBamIndex..."
    picard BuildBamIndex -I POST_ALIGNMENT_Sorted_remove-dup_RG.bam
    
    # Clean up and move results
    cd $BASE_DIR
    rm -rf $BASE_DIR/ALIGNMENT

    # Navigate to BQSR directory
    cd $BASE_DIR/BQSR

    # Base recalibration
    echo "Running GATK BaseRecalibrator..."
    gatk BaseRecalibrator  -R /path/to/reference/genome/human_g1k_v37.fasta -I ../POST_ALIGNMENT/POST_ALIGNMENT_Sorted_remove-dup_RG.bam --known-sites /path/to/gatk/resources/dbsnp_138.b37.vcf.gz --known-sites /path/to/gatk/resources/Mills_and_1000G_gold_standard.indels.b37.vcf.gz -O before-recal.table

    # Apply BQSR
    echo "Running GATK ApplyBQSR..."
    gatk ApplyBQSR  -R /path/to/reference/genome/human_g1k_v37.fasta -I ../POST_ALIGNMENT/POST_ALIGNMENT_Sorted_remove-dup_RG.bam --bqsr-recal-file before-recal.table -O recal.bam
    
    # Clean up and move results
    cd $BASE_DIR
    rm -rf $BASE_DIR/POST_ALIGNMENT

    # Navigate to VCF directory
    cd $BASE_DIR/VCF

    # Variant calling
    echo "Running GATK HaplotypeCaller..."
    gatk HaplotypeCaller  -R /path/to/reference/genome/human_g1k_v37.fasta -I ../BQSR/recal.bam -O var.vcf.gz -D /path/to/gatk/resources/dbsnp_138.b37.vcf.gz 
    
    # Make index file
    echo "Running GATK IndexFeatureFile..."
    gatk IndexFeatureFile -I var.vcf.gz
    
    # Clean up and move results
    cd $BASE_DIR
    rm -rf $BASE_DIR/BQSR

    # Navigate to VQSR directory
    cd $BASE_DIR/VQSR

    # SNP recalibration
    echo "Running GATK VariantRecalibrator for SNPs..."
    gatk VariantRecalibrator  -R /path/to/reference/genome/human_g1k_v37.fasta -mode SNP -V ../VCF/var.vcf.gz --resource:hapmap,known=false,training=true,truth=true,prior=15.0 /path/to/gatk/resources/hapmap_3.3.b37.vcf.gz --resource:omni,known=false,training=true,truth=true,prior=12.0 /path/to/gatk/resources/1000G_omni2.5.b37.vcf.gz -resource:1000G,known=false,training=true,truth=false,prior=10.0 /path/to/gatk/resources/1000G_phase1.snps.high_confidence.b37.vcf.gz --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /path/to/gatk/resources/dbsnp_138.b37.vcf.gz -tranche 99.5 -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -O recal_snp.recal --tranches-file recal_snp.tranches 

    # Apply SNP recalibration
    echo "Running GATK ApplyVQSR for SNPs..."
    gatk ApplyVQSR  -R /path/to/reference/genome/human_g1k_v37.fasta -mode SNP -V ../VCF/var.vcf.gz -ts-filter-level 99.0 --recal-file recal_snp.recal --tranches-file recal_snp.tranches -O recal_snp_raw_indel.vcf.gz

    # INDEL recalibration
    echo "Running GATK VariantRecalibrator for INDELs..."
    gatk VariantRecalibrator  -R /path/to/reference/genome/human_g1k_v37.fasta -mode INDEL -V recal_snp_raw_indel.vcf.gz --resource:mills,known=false,training=true,truth=true,prior=12.0 /path/to/gatk/resources/Mills_and_1000G_gold_standard.indels.b37.vcf.gz  --resource:dbsnp,known=true,training=false,truth=false,prior=2.0 /path/to/gatk/resources/dbsnp_138.b37.vcf.gz -tranche 99.5 -an QD -an MQ -an MQRankSum -an ReadPosRankSum -an FS -an SOR -O recal_indel.recal --tranches-file recal_indel.tranches 

    # Apply INDEL recalibration
    echo "Running GATK ApplyVQSR for INDELs..."
    gatk ApplyVQSR  -R /path/to/reference/genome/human_g1k_v37.fasta -mode INDEL -V recal_snp_raw_indel.vcf.gz -ts-filter-level 99.0 --recal-file recal_indel.recal --tranches-file recal_indel.tranches -O recal_snp_recal_indel.vcf.gz
    
    # Clean up and move results
    cd $BASE_DIR
    rm -rf $BASE_DIR/VCF
    
    # Navigate to ANNOTATION directory
    cd $BASE_DIR/ANNOTATION

    # Annotation with ANNOVAR
    echo "Running ANNOVAR..."
    table_annovar.pl ../VQSR/recal_snp_recal_indel.vcf.gz /path/to/Annovar/database/ -buildver hg19 -out $SRR -otherinfo -nastring . -vcfinput -polish -protocol refGeneWithVer,clinvar_20221231,dbnsfp47a,intervar_20180118,icgc28,popfreq_all_20150413,snp142,dbscsnv11,cosmic99_coding,cosmic99_noncoding,cytoBand -operation g,f,f,f,f,f,f,f,f,f,r -thread 12 -remove

	echo "Filtering variants..."

	# Define the PGL file path
	PGL_FILE="/path/to/gene/panel/file"  # Updated path to the PGL file

	# Check if PGL file exists
	if [[ ! -f "$PGL_FILE" ]]; then
	  echo "PGL file not found: $PGL_FILE"
	  exit 1
	fi

	# Read the list of genes into an associative array for faster lookups
	declare -A PGL_GENES_MAP
	while IFS= read -r gene; do
	  PGL_GENES_MAP["$gene"]=1
	  echo "Loaded gene: $gene"  # Debugging line to ensure genes are loaded
	done < "$PGL_FILE"

	# Define output file names based on the current working directory
	matched_pgl="matched_pgl.vcf"
	AD_GT_matched_pgl="AD_GT_matched_pgl.vcf"
	filtered_without_Y_and_benign_AD_GT_matched_pgl="filtered_without_Y_and_benign_AD_GT_matched_pgl.vcf"
	passed_filtered_without_Y_and_benign_AD_GT_matched_pgl="passed_filtered_without_Y_and_benign_AD_GT_matched_pgl.vcf"

	echo "Processing files in the current directory..."

	# Empty the output file if it exists and write VCF header
	echo -e "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE" > "$matched_pgl"

	# Function to process the VCF file
	process_vcf_file() {
	  VCF_FILE="$1"

	  echo "Processing $VCF_FILE..."

	  # Check if VCF file exists
	  if [[ ! -f "$VCF_FILE" ]]; then
		echo "VCF file not found: $VCF_FILE"
		exit 1
	  fi

	  # Use awk to filter based on gene list and write directly to the output file
	  awk -v genes="${!PGL_GENES_MAP[*]}" '
		BEGIN {
		  split(genes, geneArray, " ");
		  for (i in geneArray) geneMap[geneArray[i]] = 1;
		}
		/^#/ { print; next }
		{
		  split($8, infoFields, ";");
		  for (i in infoFields) {
			if (infoFields[i] ~ /^Gene.refGeneWithVer=/) {
			  split(infoFields[i], geneField, "=");
			  split(geneField[2], geneNames, ",");
			  for (j in geneNames) {
				if (geneNames[j] in geneMap) {
				  print;  # Output the matched line
				  next;
				}
			  }
			}
		  }
		}
	  ' "$VCF_FILE" >> "$matched_pgl"
	}

	# Identify the VCF file in the directory (assuming there's only one)
	VCF_FILE=$(ls ./*.vcf 2>/dev/null)

	if [[ -f "$VCF_FILE" ]]; then
	  process_vcf_file "$VCF_FILE"
	else
	  echo "No VCF files found in the current directory."
	  exit 1
	fi

	# Step 2: Process the output from step 1 to filter based on AD and GT values
	if [[ -f "$matched_pgl" ]]; then
	  echo "Filtering based on AD and GT from $matched_pgl to $AD_GT_matched_pgl..."

	  awk '
	  /^#/ { print; next }
	  {
	  # 1) Find indices of GT, AD, and DP in the FORMAT field··
	  split($9, fmt, ":")··
	  gtIx=adIx=dpIx=0··
	  for (i in fmt) {  
		  if (fmt[i]=="GT") gtIx=i  
		  if (fmt[i]=="AD") adIx=i  
		  if (fmt[i]=="DP") dpIx=i  
	  }  
	  # 2) If all necessary fields are present, parse the SAMPLE column··
	  if (gtIx && adIx && dpIx) {  
		  split($10, vals, ":")  
		  genotype = vals[gtIx]  
		  split(vals[adIx], ad, ",")  
		  refCount = ad[1]  
		  altCount = ad[2]  
		  dpValue  = vals[dpIx]  

		  # 3) Apply filtering by category:··
		  # a) Mosaic variants (genotype 1/2) → skip··
		  if (genotype == "1/2") {  
		  next  
		  }
		  # b) Heterozygous variants (0/1) → refCount ≥ 8, altCount ≥ 8, altCount/DP ≥ 0.25··
		  else if (genotype == "0/1") {  
		  if (refCount >= 8 && altCount >= 8 && (altCount / dpValue) >= 0.25) {  
			  print  
		  }  
		  }
		  # c) Homozygous alternate variants (1/1) → altCount ≥ 8··
		  else if (genotype == "1/1") {  
		  if (altCount >= 8) {  
		  	  print  
		  }  
		  }
	  }  
	  }
	  ' "$matched_pgl" > "$AD_GT_matched_pgl"



	  echo "Filtered: $matched_pgl -> $AD_GT_matched_pgl"
	else
	  echo "File not found: $matched_pgl"
	fi

	# Step 3: Filter based on Y chromosome and benign status
	if [[ -f "$AD_GT_matched_pgl" ]]; then
	  echo "Filtering based on chromosome and clinical significance from $AD_GT_matched_pgl to $filtered_without_Y_and_benign_AD_GT_matched_pgl..."

	  awk '
		BEGIN { FS="\t"; OFS="\t" }
		/^#/ { print; next }
		$1 != "Y" && !($8 ~ /CLNSIG=Benign(;|$)/) { print }
	  ' "$AD_GT_matched_pgl" > "$filtered_without_Y_and_benign_AD_GT_matched_pgl"

	  echo "Final filtering done: $filtered_without_Y_and_benign_AD_GT_matched_pgl"
	else
	  echo "Input file $AD_GT_matched_pgl does not exist."
	fi

	# Step 4: Filter for variants with "PASS" in the FILTER column
    echo "Running final PASS filter..."

    # Keep all header lines and only variants with FILTER field equal to PASS
    awk '$7=="PASS" || /^#/' "$filtered_without_Y_and_benign_AD_GT_matched_pgl" > "$passed_filtered_without_Y_and_benign_AD_GT_matched_pgl"

    echo "PASS-filtered VCF saved as $passed_filtered_without_Y_and_benign_AD_GT_matched_pgl"


    echo "Finished processing $SRR at $(date)"
  } >$LOG_FILE 2>$ERR_FILE
  
  # Clean up and move results
  cd $BASE_DIR
  rm -rvf  $BASE_DIR/VQSR
  mkdir -p $RESULTS_DIR
  mv -v $BASE_DIR/QC $BASE_DIR/ANNOTATION $RESULTS_DIR

  # Remove SRR directory
  rm -rvf $BASE_DIR

  echo "Completed $SRR with logs saved in $RESULTS_DIR/$SRR/ANNOTATION."
}

export -f process_srr

# Read SRR identifiers from file and process each one in parallel
cat /path/to/SRR_Acc_List.txt | parallel -j 1 process_srr
