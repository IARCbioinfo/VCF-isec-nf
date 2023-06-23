#! /usr/bin/env nextflow
nextflow.enable.dsl = 2

// Copyright (C) 2023 IARC/WHO

// This program is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.

// This program is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.

// You should have received a copy of the GNU General Public License
// along with this program.  If not, see <http://www.gnu.org/licenses/>.

params.help = null

log.info ""
log.info "--------------------------------------------------------"
log.info "  Nextflow Pipeline: Intersection of VCF Files          "
log.info "--------------------------------------------------------"
log.info "Copyright (C) 2023 IARC/WHO"
log.info "This program comes with ABSOLUTELY NO WARRANTY; for details see LICENSE"
log.info "This is free software, and you are welcome to redistribute it"
log.info "under certain conditions; see LICENSE for details."
log.info "--------------------------------------------------------"
log.info ""

if (params.help) {
    log.info "--------------------------------------------------------"
    log.info "  USAGE                                                 "
    log.info "--------------------------------------------------------"
    log.info ""
    log.info "nextflow run pipeline.nf [-with-docker] [OPTIONS]"
    log.info ""
    log.info "Mandatory arguments:"
    log.info "--input_folder1   <DIR>   Folder containing VCF files"
    log.info "--input_folder2   <DIR>   Folder containing other VCF files"
    log.info "--output_folder      <DIR>   Output folder"
    log.info ""
    log.info "Optional arguments:"
    log.info "--vcfSuffix1_snvs   <STRING>  Suffix (without the extension) of files in"
    log.info "                              input_folder1 containing SNVs"
    log.info "--vcfSuffix1_indels <STRING>  Suffix (without the extension) of files in"
    log.info "                              input_folder1 containing indels and MNVs"
    log.info "--vcfSuffix2_snvs   <STRING>  Suffix (without the extension) of files in"
    log.info "                              input_folder2 containing SNVs"
    log.info "--vcfSuffix2_indels <STRING>  Suffix (without the extension) of files in"
    log.info "                              input_folder2 containing indels and MNVs"
    log.info "--ext               <STRING>  Extension of variant calling files"

// Cluster options
params.mem = 8
params.cpu = 2
    exit 0
} else {
    /* Software information */
    log.info "help:               ${params.help}"
    log.info "input_folder1:    ${params.input_folder1}"
    log.info "input_folder2:   ${params.input_folder2}"
    log.info "output_folder:       ${params.output_folder}"
}


// Define the pipeline parameters
params.input_folder1 = '/path/to/VCF1/vcf/files'
params.input_folder2 = '/path/to/VCF2/vcf/files'

// Output folder
params.output_folder = 'output/'

// File suffixes
params.vcfSuffix1_snvs   = '_filtered_PASS_norm.vcf.hg38_multianno.vcf'
params.vcfSuffix1_indels = '_filtered_PASS_norm.vcf.hg38_multianno.vcf'
params.vcfSuffix2_snvs   = '.somatic.snvs_norm.vcf'
params.vcfSuffix2_indels = '.somatic.indels_norm.vcf'
params.ext = '.gz'

// Cluster options
params.mem = 8
params.cpu = 2

process concatsnvsindels {
  memory params.mem+'GB'
  cpus params.cpu

  input:
  tuple val(sampleId), path(SNVs), path(SNVstbi), path(indels), path(indelstbi)

  output:
  tuple val(sampleId), path("${sampleId}.mnps_indels.vcf.gz") , path("${sampleId}.mnps_indels.vcf.gz.tbi") 

  script:
  """
  # Concatenate SNV and indels/MNPs files
  mkdir sort_tmp
  bcftools concat $SNVs $indels -a -Ou | bcftools view -v indels,mnps -Ou | bcftools sort -m ${params.mem}G -T sort_tmp/ -Oz -o ${sampleId}.mnps_indels.vcf.gz

  # Index the concatenated VCF file
  bcftools index -t ${sampleId}.mnps_indels.vcf.gz
  """
}

process extractSNVs {
  memory params.mem+'GB'
  cpus params.cpu

  input:
  tuple val(sampleId), path(VCF1File), path(VCF1Filetbi), path(VCF2File), path(VCF2Filetbi)

  output:
  tuple val(sampleId), path("snvs_${sampleId}.vcf.gz")

  script:
  """
  # Extract SNVs from VCF1 VCF
  bcftools view -v snps $VCF1File -Oz -o snvs_${sampleId}.vcf.gz
  """
}

process intersectIndelsMNPs {
  memory params.mem+'GB'
  cpus params.cpu

  input:
  tuple val(sampleId), path(VCF1File), path(VCF1Filetbi), path(VCF2File), path(VCF2Filetbi)

  output:
  tuple val(sampleId), path("indels_mnps_${sampleId}.vcf.gz") 

  script:
  """
  # Intersect VCF1 and VCF2 VCF files for indels and MNPs
  bcftools isec -Oz -n =2 $VCF1File $VCF2File -p intersect_${sampleId}
  mv intersect_${sampleId}/0000.vcf.gz indels_mnps_${sampleId}.vcf.gz
  """
}


process concatAndIndexVCF {
  memory params.mem+'GB'
  cpus params.cpu

  input:
  tuple val(sampleId), path(SNVs), path(indels_MNPs) 

  output:
  path("${sampleId}.snvs_plus_shared_indels.vcf.hg38_multianno.vcf.gz*") 

  publishDir "${params.output_folder}", mode: 'move'

  script:
  """
  # Index input files
  bcftools index -t $SNVs
  bcftools index -t $indels_MNPs

  # Concatenate SNV and indels/MNPs files
  mkdir sort_tmp
  bcftools concat $SNVs $indels_MNPs -a -Ou | bcftools sort -m ${params.mem}G -T sort_tmp/ -Oz -o ${sampleId}.snvs_plus_shared_indels.vcf.hg38_multianno.vcf.gz

  # Index the concatenated VCF file
  bcftools index -t ${sampleId}.snvs_plus_shared_indels.vcf.hg38_multianno.vcf.gz
  """
}

workflow {

// Create separate channels for VCF1 and VCF2 files
// Channel with SNVs
  VCF1Files_snvs = Channel
  .fromPath("$params.input_folder1/*$params.vcfSuffix1_snvs$params.ext")
  .map { file ->  [file.baseName.replaceAll(params.vcfSuffix1_snvs, ""), file,file+".tbi"]}
  
if(params.vcfSuffix1_snvs != params.vcfSuffix1_indels ){
  println "Separate SNV and indel files for VCF 1, concatenating"
  // Channel with indels
  VCF1Files_indels = Channel
  .fromPath("$params.input_folder1/*$params.vcfSuffix1_indels$params.ext")
  .map { file ->  [file.baseName.replaceAll(params.vcfSuffix1_indels, ""), file,file+".tbi"]}
  // Channel with both files
  VCF1Files2concat = VCF1Files_snvs.join(VCF1Files_indels)
  // Channel with concatenated files
  VCF1Files = concatsnvsindels(VCF1Files2concat)
}else{
  VCF1Files = VCF1Files_snvs
}
//VCF1Files = VCF1Files.view()

// same for second set of VCFs
VCF2Files_snvs = Channel
  .fromPath("$params.input_folder2/*$params.vcfSuffix2_snvs$params.ext")
  .map { file ->  [file.baseName.replaceAll(params.vcfSuffix2_snvs, ""), file,file+".tbi"]}

if(params.vcfSuffix2_snvs != params.vcfSuffix2_indels ){
  println "Separate SNV and indel files for VCF 2, concatenating"
  // Channel with indels
  VCF2Files_indels = Channel
  .fromPath("$params.input_folder2/*$params.vcfSuffix2_indels$params.ext")
  .map { file ->  [file.baseName.replaceAll(params.vcfSuffix2_indels, ""), file,file+".tbi"]}
  // Channel with both files
  VCF2Files2concat = VCF2Files_snvs.join(VCF2Files_indels).view()
  // Channel with concatenated files
  VCF2Files = concatsnvsindels(VCF2Files2concat)
}else{
  VCF2Files = VCF2Files_snvs
}
VCF2Files = VCF2Files

  // Join VCF1 and VCF2 channels using sample ID as key
intersectChannel = VCF1Files
  .join(VCF2Files).view()

  // Run the processes in parallel
  extractSNVs(intersectChannel)
  intersectIndelsMNPs(intersectChannel)

  // Join outputs
  intersectChannel2 = extractSNVs.out
                                 .join(intersectIndelsMNPs.out)

  // Run the final process after the previous ones complete
  concatAndIndexVCF(intersectChannel2)
}