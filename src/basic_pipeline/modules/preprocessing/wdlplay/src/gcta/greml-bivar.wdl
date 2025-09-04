version 1.0

workflow GREML_biv {
    input {
        File pheno
        File mgrm
        String outdir
        Array[String] phenotypes
        File covar
        File qcovar
        Int mem
    }

    scatter(k in range(length(phenotypes))) {
        call calc_covar {
            input:
            index=k,
            pheno=pheno,
            pheno1=phenotypes[k],
            npheno=length(phenotypes),
            mgrm=mgrm,
            outdir=outdir,
            phenotypes=phenotypes,
            covar=covar,
            qcovar=qcovar,
            mem=mem
        }
    }

    output {
        Array[File?] hsq = calc_covar.hsq
    }
}

task calc_covar {
    input {
        Int index
        File pheno
        String pheno1
        File mgrm
        String outdir
        Array[String] phenotypes
        Int npheno
        File covar
        File qcovar
        Int mem
    }

    command <<<
        pheno_array=('~{sep="' '" phenotypes}')
        END=~{npheno}
        for ((i=~{index};i<$((END-1));i++)); do
        gcta --mgrm "~{mgrm}" --pheno "~{pheno}" \
        --reml-bivar ~{index+1} $((i+2)) --qcovar "~{qcovar}" \
        --covar "~{covar}" --thread-num 8 \
        --out "~{outdir}/~{pheno1}_"${pheno_array[(($i+1))]}
        done
    >>>

    runtime {
        cpu: 4
        memory: "~{mem} GB"
        # singularity: "docker://gcr.io/singh-comp-d-271c/genetics:0.2.18"
        singularity: "/gpfs/commons/groups/singh_lab/software/docker/genetics_0.2.18.sif"
    }

    output {
        File hsq = "~{outdir}/~{pheno1}_~{phenotypes[(npheno-1)]}.hsq"
    }
}