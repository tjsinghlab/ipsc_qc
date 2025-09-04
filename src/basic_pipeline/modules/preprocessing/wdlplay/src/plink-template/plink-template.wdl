version 1.0

workflow plink_workflow {
  input {
    String input_dir
    String output_dir=input_dir
    String input_prefix
    String docker
  }

  String inter_dir = output_dir+'chr/'

  call Split_PLINK_By_Chromosome {
    input:
      input_dir = input_dir,
      output_dir = inter_dir,
      input_prefix = input_prefix,
      docker = docker
  }

  scatter (bim_file in Split_PLINK_By_Chromosome.bim_files) {
    call PLINK_QC {
      input:
        input_dir = inter_dir,
        output_dir = output_dir+'clean/',
        input_prefix = basename(bim_file, '.bim'),
        docker = docker
    }
  }

  output {
    Array[File] clean_bed_files = PLINK_QC.clean_dataset_bed
    Array[File] clean_bim_files = PLINK_QC.clean_dataset_bim
    Array[File] clean_fam_files = PLINK_QC.clean_dataset_fam
  }
}

task Split_PLINK_By_Chromosome {
  input {
    String input_dir
    String output_dir
    String input_prefix
    Int start_chr = 1
    Int end_chr = 22
    String docker
  }

  command {
    set -e
    mkdir -p ~{output_dir}
    for chr in $(seq ~{start_chr} ~{end_chr})
    do
      plink --bfile ~{input_dir}~{input_prefix} --chr $chr --make-bed --threads 8 --out ~{output_dir}~{input_prefix}_chr$chr
    done
  }

  output {
    Array[File] bed_files = glob("~{output_dir}*_chr*.bed")
    Array[File] bim_files = glob("~{output_dir}*_chr*.bim")
    Array[File] fam_files = glob("~{output_dir}*_chr*.fam")
  }

  runtime {
    cpus: 8
    cpu: 8
    mem: "16G"
    memory: "16G"
    singularity: docker
  }
}

task PLINK_QC {
  input {
    String input_dir
    String output_dir
    String input_prefix
    Float geno_threshold = 0.05
    Float maf_threshold = 0.01
    Float mind_threshold = 0.02
    String docker
  }

  command {
    mkdir -p ~{output_dir}
    # Variant filtering based on missingness and MAF
    plink --bfile ~{input_dir}${input_prefix} --geno ${geno_threshold} --maf ${maf_threshold} --make-bed --out intermediate_dataset

    # Sample filtering based on missing genotype data
    plink --bfile intermediate_dataset --mind ${mind_threshold} --make-bed --out ~{output_dir}${input_prefix}_clean
  }

  output {
    File clean_dataset_bed = "~{output_dir}${input_prefix}_clean.bed"
    File clean_dataset_bim = "~{output_dir}${input_prefix}_clean.bim"
    File clean_dataset_fam = "~{output_dir}${input_prefix}_clean.fam"
  }

  runtime {
    cpus: 8
    cpu: 8
    mem: "16G"
    memory: "16G"
    singularity: docker
  }
}