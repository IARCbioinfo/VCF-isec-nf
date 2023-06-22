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
    log.info "--vcfFolder1   <DIR>   Folder containing VCF1 VCF files"
    log.info "--vcfFolder2  <DIR>   Folder containing VCF2 VCF files"
    log.info "--outputFolder      <DIR>   Output folder"
    log.info ""
    log.info "Optional arguments:"
    log.info ""
    log.info "Flags:"
    log.info ""
    exit 0
} else {
    /* Software information */
    log.info "help:               ${params.help}"
    log.info "vcfFolder1:    ${params.vcfFolder1}"
    log.info "vcfFolder2:   ${params.vcfFolder2}"
    log.info "outputFolder:       ${params.outputFolder}"
}


// Define the pipeline parameters
params.vcfFolder1 = '/path/to/VCF1/vcf/files'
params.vcfFolder2 = '/path/to/VCF2/vcf/files'

// Output folder
params.outputFolder = 'output/'

// File suffixes
params.vcfSuffix1 = '_filtered_PASS_norm.vcf.hg38_multianno.vcf'
params.vcfSuffix2 = '.somatic.indels_norm.vcf.hg38_multianno.vcf'
params.ext = '.gz'

params.mem = 8

    // Separate channels for VCF1 and VCF2 files
VCF1Files = Channel
  .fromPath("$params.vcfFolder1/*$params.vcfSuffix1$params.ext")
  .map { file ->  [file.baseName.replaceAll(params.vcfSuffix1, ""), file,file+".tbi"]}

VCF2Files = Channel
  .fromPath("$params.vcfFolder2/*$params.vcfSuffix2$params.ext")
  .map { file ->  [file.baseName.replaceAll(params.vcfSuffix2, ""), file,file+".tbi"]}


  // Join VCF1 and VCF2 channels using sample ID as key
intersectChannel = VCF1Files
  .join(VCF2Files).view()

process extractSNVs {
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

  input:
  tuple val(sampleId), path(SNVs), path(indels_MNPs) 

  output:
  path("${sampleId}.snvs_plus_shared_indels.vcf.hg38_multianno.vcf.gz*") 

  publishDir "${params.outputFolder}", mode: 'move'

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
  // Run the processes in parallel
  extractSNVs(intersectChannel)
  intersectIndelsMNPs(intersectChannel)

  // Join outputs
  intersectChannel2 = extractSNVs.out
                                 .join(intersectIndelsMNPs.out)

  // Run the final process after the previous ones complete
  concatAndIndexVCF(intersectChannel2)
}