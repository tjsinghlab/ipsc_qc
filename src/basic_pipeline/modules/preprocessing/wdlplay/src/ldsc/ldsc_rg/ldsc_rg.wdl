version 1.0
workflow LDSC {
    input {
        Array[File] raw_sumstats_gz
        String merge_alleles
        Array[Int] sample_size
        Int? chunksize
        String munge_outdir
        String rg_outdir
        String panel="1KG"
        Array[String] phenotypes
    }

    scatter (index in range(length(raw_sumstats_gz))) {
        call munge_sumstats {
            input:
            raw_sumstats=raw_sumstats_gz[index],
            merge_alleles=merge_alleles,
            sample_size=sample_size[index],
            chunksize=chunksize,
            panel=panel,
            munge_outdir=munge_outdir
        }
    }

    scatter (index in range(length(munge_sumstats.munged_sumstats))) {
        scatter (index2 in range(length(munge_sumstats.munged_sumstats))) {
            if (index < index2) {
                call ldsc_rg_job {
                    input:
                    outdir=rg_outdir,
                    infiles=[munge_sumstats.munged_sumstats[index], munge_sumstats.munged_sumstats[index2]],
                    phenotypes=[phenotypes[index], phenotypes[index2]]
                }
            }
        }
    }    
        
    output {
        Array[Array[File?]] rg_job_log = ldsc_rg_job.rg_job_log
    }
}


task munge_sumstats {
    input {
        String merge_alleles
        Int sample_size
        Int chunksize=500000
        String munge_outroot = basename(raw_sumstats,".tsv.gz")+".munged"
        File raw_sumstats
        String panel
        String munge_outdir
    }

    command {
        if [ ! -f "~{munge_outdir}/~{munge_outroot}.~{panel}.sumstats.gz" ]; then
            mkdir -p ~{munge_outdir}
            if [ "~{panel}" == "1KG" ]; then
                python /ldsc/munge_sumstats.py \
                --merge-alleles ~{merge_alleles} \
                --N ~{sample_size} \
                --chunksize ~{chunksize} \
                --ignore OR,Z,v,rsid,AF,palin,flip_palin,flip \
                --sumstats ~{raw_sumstats} \
                --out ~{munge_outdir}"/"~{munge_outroot}"."~{panel}
            elif [ "~{panel}" == "UKBB" ]; then
                python /ldsc/munge_sumstats.py \
                --out ~{munge_outdir}"/"~{munge_outroot} \
                --N ~{sample_size} \
                --chunksize ~{chunksize} \
                --ignore SNP,OR,Z,rsid,AF,palin,flip_palin,flip \
                --snp v \
                --sumstats ~{raw_sumstats}
            else
                echo "Panel not recognized. Please use 1KG or UKBB."
            fi
        fi
    }

    runtime {
        cpu: 1
        memory: "16 GB"
        #singularity: "docker://gcr.io/singh-comp-d-271c/ldsc-base:latest"
        docker: "gcr.io/singh-comp-d-271c/ldsc-base:latest"
        singularity: "/gpfs/commons/groups/singh_lab/software/docker/ldsc-base/ldsc-base_latest.sif"
    }

    output {
       # File munged_sumstats_log = "~{munge_outdir}/~{munge_outroot}.~{panel}.log"
        File munged_sumstats = "~{munge_outdir}/~{munge_outroot}.~{panel}.sumstats.gz"
    }
}

task ldsc_rg_job {
    input {
        String outdir
        Array[File] infiles
        Array[String] phenotypes
        String rgout="~{outdir}/~{phenotypes[0]}_~{phenotypes[1]}"
    }

    command {
        mkdir -p ~{outdir}
        python /ldsc/ldsc.py \
        --rg ~{infiles[0]},~{infiles[1]} \
        --ref-ld-chr /ldsc/eur_w_ld_chr/ \
        --w-ld-chr /ldsc/eur_w_ld_chr/ \
        --out rged > ~{rgout}".log"
    }

    runtime {
        cpu:1
        memory: "4GB"
        # singularity: "docker://gcr.io/singh-comp-d-271c/ldsc-base:latest"
        docker: "gcr.io/singh-comp-d-271c/ldsc-base:latest"
        singularity: "/gpfs/commons/groups/singh_lab/software/docker/ldsc-base/ldsc-base_latest.sif"
    }

    output {
        File rg_job_log = "~{rgout}.log"
    }
}
