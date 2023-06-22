#! /usr/bin/env nextflow

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
    log.info "--vcfFolderMutect   <DIR>   Folder containing Mutect VCF files"
    log.info "--vcfFolderStrelka  <DIR>   Folder containing Strelka VCF files"
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
    log.info "vcfFolderMutect:    ${params.vcfFolderMutect}"
    log.info "vcfFolderStrelka:   ${params.vcfFolderStrelka}"
    log.info "outputFolder:       ${params.outputFolder}"
}

// Define the pipeline name and description
manifest {
  name = 'vcf-intersection-pipeline'
  description = 'Nextflow pipeline for intersecting VCF files using bcftools'
}

// Define the pipeline parameters
params {
  // Input folders
  vcfFolderMutect = '/path/to/mutect/vcf/files'
  vcfFolderStrelka = '/path/to/strelka/vcf/files'

  // Output folder
  outputFolder = '/path/to/output/folder'

  // File suffixes
  vcfSuffixMutect = '_filtered_PASS_norm.vcf.hg38_multianno.vcf.gz'
  vcfSuffixStrelka = 'somatic.snvs_norm.vcf.hg38_multianno.vcf.gz'
}

// Separate channels for Mutect and Strelka files
mutectFiles = Channel
  .fromPath("$params.vcfFolderMutect/*$params.vcfSuffixMutect")
  .map { file -> file.baseName =~ /S(\d+)_/; [sampleId: $1, file: file] }

strelkaFiles = Channel
  .fromPath("$params.vcfFolderStrelka/*$params.vcfSuffixStrelka")
  .map { file -> file.baseName =~ /S(\d+)\./; [sampleId: $1, file: file] }

// Join Mutect and Strelka channels using sample ID as key
intersectChannel = mutectFiles
  .keyBy { it.sampleId }
  .join(strelkaFiles.keyBy { it.sampleId })

process extractSNVs {
  input:
  tuple val(sampleId), file(mutectFile), file(strelkaFile) from intersectChannel

  output:
  file("snvs_${sampleId}.vcf.gz") into snvFiles

  script:
  """
  # Extract SNVs from Mutect VCF
  bcftools view -f 'TYPE="snp"' $mutectFile -Oz -o snvs_${sampleId}.vcf.gz
  """
}

process intersectIndelsMNPs {
  input:
  tuple val(sampleId), file(mutectFile), file(strelkaFile) from intersectChannel

  output:
  file("indels_mnps_${sampleId}.vcf.gz") into indelMnpFiles

  script:
  """
  # Intersect Mutect and Strelka VCF files for indels and MNPs
  bcftools isec -Oz -n =2 $mutectFile $strelkaFile -p intersect_${sampleId}
  mv intersect_${sampleId}/0002.vcf.gz indels_mnps_${sampleId}.vcf.gz
  rm -r intersect_${sampleId}
  """
}

process concatAndIndexVCF {
  input:
  file snvFile from snvFiles
  file indelMnpFile from indelMnpFiles

  output:
  file("${params.outputFolder}/concatenated.vcf.gz") into finalVCF

  script:
  """
  # Concatenate SNV and indels/MNPs files
  bcftools concat $snvFile $indelMnpFile -Oz -o ${params.outputFolder}/concatenated.vcf.gz

  # Index the concatenated VCF file
  tabix -p vcf ${
