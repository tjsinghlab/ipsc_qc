version 1.0

workflow LDSC {
    input {
        File raw_sumstats_gz
        Int sample_size
        Int? chunksize
        String ldsc_h2_outdir
        String phenotype
        String liability="observed"
        String ref="default"
        String samprev="NA"
        String poprev="NA" 
        String ancestry="EUR"
    }

    # call unzip_sumstats {
    #     input:
    #     raw_sumstats_gz=raw_sumstats_gz
    # }

    call munge_sumstats {
        input:
        raw_sumstats=raw_sumstats_gz,
        sample_size=sample_size,
        chunksize=chunksize,
        ldsc_h2_outdir=ldsc_h2_outdir,
        ref=ref
    }

    call calculate_ldsc_h2_job {
        input:
            ldsc_h2_outdir=ldsc_h2_outdir,
            infile=munge_sumstats.munged_sumstats,
            phenotype=phenotype,
            liability=liability,
            ref=ref,
            samprev=samprev,
            poprev=poprev,
            ancestry=ancestry
    }

    output {
        File munged_sumstats = munge_sumstats.munged_sumstats
        File h2_log = calculate_ldsc_h2_job.h2_log
    }
}

task unzip_sumstats {
    input {
        File raw_sumstats_gz
        String raw_sumstats_path = basename(raw_sumstats_gz,".gz")
    }

    command {
        gunzip -cd ~{raw_sumstats_gz} > ~{raw_sumstats_path}
        echo ~{raw_sumstats_path}
    }

    runtime {
        cpu: 1
        memory: "1 GB"
    }

    output {
        File raw_sumstats = "~{raw_sumstats_path}"
    }
}

task munge_sumstats {
    input {
        Int sample_size
        Int chunksize=500000
        String ldsc_h2_outdir
        String munge_outroot = ldsc_h2_outdir+"/"+basename(raw_sumstats,".tsv.gz")+".munged"
        File raw_sumstats
        String ref
        Boolean ukb=(ref=="UKB")
    }

    command {
        python /ldsc/munge_sumstats.py \
        ~{true="" false="--merge-alleles /ldsc/w_hm3.snplist " ukb}\
        --N ~{sample_size} \
        --chunksize ~{chunksize} \
        --ignore ~{if ukb then "SNP,OR,Z,rsid,AF,palin,flip_palin,flip" else "OR,Z,v,rsid,AF,palin,flip_palin,flip"} \
        ~{true="--snp v " false="" ukb}\
        --sumstats ~{raw_sumstats} \
        --out ~{munge_outroot} > ~{munge_outroot}".log"
    }

    runtime {
        cpu: 1
        memory: "16 GB"
        # singularity: "docker://gcr.io/singh-comp-d-271c/ldsc-base:latest"
        singularity: "/gpfs/commons/groups/singh_lab/software/docker/ldsc-base/ldsc-base_latest.sif"
    }

    output {
        File munged_sumstats_log = "~{munge_outroot}.log"
        File munged_sumstats = "~{munge_outroot}.sumstats.gz"
    }
}

task calculate_ldsc_h2_job {
    input {
        String ldsc_h2_outdir
        File infile
        String phenotype
        String liability="observed"
        String ref="default"
        String samprev="NA"
        String poprev="NA"
        String ancestry="EUR"
        Boolean lib_true=(liability=="liability")
    }

    command {
        if [ "~{ref}" == "default" ]; then
            if [ "~{liability}" == "liability" ]; then
                echo "Running with 1kg as ref (default) and on liability scale"
            else
                echo "Running with 1kg as ref (default) and on observed scale (default)"
            fi
            python /ldsc/ldsc.py \
            --out heritability \
            --ref-ld-chr /ldsc/1000G_Phase3_ldscores/LDscore. \
            --w-ld-chr /ldsc/1000G_Phase3_ldscores/LDscore. \
            --frqfile-chr /ldsc/1000G_Phase3_frq/1000G.EUR.QC. \
            --maf 0.05 \
            --chisq-max 9999 \
            --h2 ~{infile} \
            ~{if lib_true then "--samp-prev " + samprev + " --pop-prev " + poprev else ""} > "~{ldsc_h2_outdir}/~{phenotype}_h2.log"
        elif [ "~{ref}" == "hapmap" ]; then
            if [ "~{liability}" == "liability" ]; then
                echo "Running with hapmap as ref and on liability scale"
            else
                echo "Running with hapmap as ref and on observed scale (default)"
            fi
            python /ldsc/ldsc.py \
            --out heritability \
            --ref-ld /ldsc/eur_w_ld_chr/ \
            --w-ld /ldsc/eur_w_ld_chr/ \
            --h2 "~{infile}" \
            ~{if lib_true then "--samp-prev " + samprev + " --pop-prev " + poprev else ""} > "~{ldsc_h2_outdir}/~{phenotype}_h2.log"                
        elif [ "~{ref}" == "UKB" ]; then
            if [ "~{liability}" == "liability" ]; then
                echo "Running with UKB."~{ancestry}" as ref and on liability scale"
            else
                echo "Running with UKB."~{ancestry}" as ref and on observed scale (default)"
            fi
            python /ldsc/ldsc.py \
            --out heritability \
            --ref-ld /ldsc/UKBB.ALL.ldscore/UKBB."~{ancestry}" \
            --w-ld /ldsc/UKBB.ALL.ldscore/UKBB."~{ancestry}" \
            --h2 "~{infile}" \
            ~{if lib_true then "--samp-prev " + samprev + " --pop-prev " + poprev else ""} > "~{ldsc_h2_outdir}/~{phenotype}_h2.log"                
        else
            echo "Invalid ref argument. Please use 'default', 'hapmap', or 'UKB'."
        fi
    }
    
    runtime {
        cpu: 1
        memory: "7 GB"
        # singularity: "docker://gcr.io/singh-comp-d-271c/ldsc-annot:0.2.19"
        singularity: "/gpfs/commons/groups/singh_lab/software/docker/ldsc-annot/ldsc-annot_0.2.19.sif"
    }

    output {
        File h2_log = "~{ldsc_h2_outdir}/~{phenotype}_h2.log"
    }
}

# task run_sldsc_job {
#     input {
#         String outdir
#         File sumstats
#         String phe_label
#         String geneset_zip?
#         String geneset_prefix?
#     }

#     command {
#         String outroot = $outdir + "/" + $phe_labelo
#         String ref_ld_chr = '/ldsc/1000G_Phase3_baselineLD_v2.2_ldscores/baselineLD.'

#         if (defined(geneset_zip) && defined(geneset_prefix)) {
#             String ref_ld_chr = "~{geneset_prefix},~{ref_ld_chr}"
#             String outroot = outroot + "_~{geneset_prefix}"
#             tar -xvf ~{geneset_zip}
#         }
#         python ldsc.py 
#             --h2 ~{sumstats} \
#             --ref-ld-chr ~{ref_ld_chr} \
#             --w-ld-chr /ldsc/1000G_Phase3_weights_hm3_no_MHC/weights.hm3_noMHC. \
#             --frqfile-chr /ldsc/1000G_Phase3_frq/1000G.EUR.QC. \
#             --overlap-annot \
#             --maf 0.05 \
#             --chisq-max 9999 \
#             --print-coefficients \
#             --out ~{outroot} + ".log"
#     }

#     runtime {
#         cpu:2
#         memory: "8 GB"
#     }

#     output {
#         File outlog = ~{outroot} + ".log"
#     }
# }