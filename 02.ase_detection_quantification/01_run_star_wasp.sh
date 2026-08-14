#!/bin/bash

set -euo pipefail

if [[ $# -ne 9 ]]; then
  echo "Usage: $0 <sample_id> <read1.fq.gz> <read2.fq.gz> <het.vcf.gz> <genome_dir> <annotation.gtf> <threads> <output_dir> <stats_tsv>" >&2
  exit 1
fi

sample_id="$1"
read1="$2"
read2="$3"
het_vcf="$4"
genome_dir="$5"
annotation_gtf="$6"
threads="$7"
output_dir="$8"
stats_tsv="$9"

mkdir -p "${output_dir}"

vcf_for_star="${output_dir}/${sample_id}.het.biallelic.snv.vcf"
prefix="${output_dir}/${sample_id}_"
raw_bam="${prefix}Aligned.sortedByCoord.out.bam"
filtered_bam="${output_dir}/${sample_id}.WASP_filtered.bam"
log_final="${prefix}Log.final.out"

bcftools view -m2 -M2 -v snps -Ov "${het_vcf}" > "${vcf_for_star}"

STAR \
  --runThreadN "${threads}" \
  --genomeDir "${genome_dir}" \
  --sjdbGTFfile "${annotation_gtf}" \
  --twopassMode Basic \
  --outFilterMismatchNmax 999 \
  --readFilesIn "${read1}" "${read2}" \
  --readFilesCommand zcat \
  --outFileNamePrefix "${prefix}" \
  --outSAMtype BAM SortedByCoordinate \
  --limitBAMsortRAM 80000000000 \
  --varVCFfile "${vcf_for_star}" \
  --waspOutputMode SAMtag \
  --outSAMattributes NH HI AS nM NM MD jM jI MC vA vG vW

samtools index -@ "${threads}" "${raw_bam}"

samtools view -h -@ "${threads}" "${raw_bam}" | \
awk 'BEGIN{OFS="\t"}
     /^@/ {print; next}
     {
       tagged=0; pass=0;
       for(i=12;i<=NF;i++){
         if($i ~ /^vW:i:/){
           tagged=1;
           if($i=="vW:i:1"){pass=1}
           break
         }
       }
       if(tagged==0 || pass==1){print}
     }' | \
samtools view -@ "${threads}" -b -o "${filtered_bam}" -

samtools index -@ "${threads}" "${filtered_bam}"

vcf_n=$(grep -vc '^#' "${vcf_for_star}" || true)
raw_primary=$(samtools view -@ "${threads}" -c -F 2308 "${raw_bam}")
filtered_primary=$(samtools view -@ "${threads}" -c -F 2308 "${filtered_bam}")

get_star_value() {
  local key="$1"
  awk -F'|' -v k="$key" '
    $1 ~ k {
      gsub(/^[ \t]+|[ \t]+$/, "", $2)
      print $2
      exit
    }' "${log_final}"
}

{
  echo -e "sample_id\tvcf_biallelic_snvs\tstar_input_reads\tuniquely_mapped_reads\tuniquely_mapped_rate\tmulti_mapped_reads\tmulti_mapped_rate\tprimary_mapped_before_wasp\tprimary_mapped_after_wasp\tfiltered_bam"
  echo -e "${sample_id}\t${vcf_n}\t$(get_star_value "Number of input reads")\t$(get_star_value "Uniquely mapped reads number")\t$(get_star_value "Uniquely mapped reads %")\t$(get_star_value "Number of reads mapped to multiple loci")\t$(get_star_value "% of reads mapped to multiple loci")\t${raw_primary}\t${filtered_primary}\t${filtered_bam}"
} > "${stats_tsv}"
