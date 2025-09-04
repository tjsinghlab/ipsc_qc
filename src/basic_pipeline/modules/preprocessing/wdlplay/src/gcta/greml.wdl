version 1.0

workflow GREML {
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
        call calc_h2 {
            input:
            index=k,
            pheno=pheno,
            mgrm=mgrm,
            outdir=outdir,
            phenotype=phenotypes[k],
            covar=covar,
            qcovar=qcovar,
            mem=mem
        }
    }

    output {
        Array[File?] hsq = calc_h2.hsq
    }
}

task calc_h2 {
    input {
        Int index
        File pheno
        File mgrm
        String outdir
        String phenotype
        File covar
        File qcovar
        Int mem
    }

    command <<<
        gcta \
        --mgrm "~{mgrm}" --pheno "~{pheno}" --reml \
        --mpheno ~{index+1} --qcovar "~{qcovar}" \
        --covar "~{covar}" --thread-num 8 \
        --reml-lrt 1 2 3 4 \
        --reml-maxit 1000 \
        --out "~{outdir}/~{phenotype}"
    >>>

    runtime {
        cpu: 4
        memory: "~{mem} GB"
        # singularity: "docker://gcr.io/singh-comp-d-271c/genetics:0.2.18"
        singularity: "/gpfs/commons/groups/singh_lab/software/docker/genetics_0.2.18.sif"
    }

    output {
        File hsq = "~{outdir}/~{phenotype}.hsq"
    }
}