version 1.0

workflow GRM {
    input {
        String bfile
        Boolean ldms=false
        Int ld_window=200
        String outdir
    }

    if (ldms) {
        call segment_ld {
            input:
                bfile=bfile,
                ld_window=ld_window,
                outdir=outdir
        }
        
        call stratify_ld {
            input:
                ldscore=segment_ld.ldscore,
                outdir=outdir
        }
        
        scatter(snp in stratify_ld.snp_groups) {
            call ldms_grm {
                input:
                    bfile=bfile,
                    snp=snp,
                    outdir=outdir
            }
        }
        # output {
        #     Array[File] grm = ldms_grm.grm
        # }
    }
    if (!ldms){
        call make_grm {
            input:
                bfile=bfile,
                outdir=outdir
        }
        # output {
        #     File grm = make_grm.grm
        # }
    }

    output {
        Array[File]? grm_list = ldms_grm.grm
        File? grm = make_grm.grm
    }
}

task segment_ld {
    input {
        String bfile
        Int ld_window     
        String outdir
        String out=outdir+basename(bfile)
    }

    command {
        if [ -e "~{out}.score.ld" ]; then
            echo "LD score file already exists"
        else
            gcta --bfile ~{bfile} --autosome --ld-score-region ~{ld_window} \
            --out ~{out} --thread-num 4
        fi
    }

    runtime {
        cpu: 4
        memory: "50 GB"
        #singularity: "docker://gcr.io/singh-comp-d-271c/genetics:0.2.18"
        singularity: "/gpfs/commons/groups/singh_lab/software/docker/genetics_0.2.18.sif"
    }

    output {
        File ldscore = "~{out}.score.ld"
    }
}

task stratify_ld {
    input {
        File ldscore
        String outdir
    }

    command {
        R --no-save <<CODE
        if(all(file.exists(c("~{outdir}/snp_group1.txt", "~{outdir}/snp_group2.txt",
        "~{outdir}/snp_group3.txt", "~{outdir}/snp_group4.txt")))
        ) quit()
        lds_seg = read.table("~{ldscore}",header = T,
        colClasses = c("character", rep("numeric", 8)))
        quartiles = summary(lds_seg[, "ldscore_SNP"])
    
        lb1 = which(lds_seg[, "ldscore_SNP"] <= quartiles[2])
        lb2 = which(lds_seg[, "ldscore_SNP"] > quartiles[2] & lds_seg[, "ldscore_SNP"] <= quartiles[3])
        lb3 = which(lds_seg[, "ldscore_SNP"] > quartiles[3] & lds_seg[, "ldscore_SNP"] <= quartiles[5])
        lb4 = which(lds_seg[, "ldscore_SNP"] > quartiles[5])
        
        lb1_snp = lds_seg[lb1, "SNP"]
        lb2_snp = lds_seg[lb2, "SNP"]
        lb3_snp = lds_seg[lb3, "SNP"]
        lb4_snp = lds_seg[lb4, "SNP"]
        
        write.table(lb1_snp, file.path("~{outdir}", "snp_group1.txt"),
        row.names = F,quote = F, col.names = F)
        write.table(lb2_snp, file.path("~{outdir}", "snp_group2.txt"),
        row.names = F,quote = F,col.names = F)
        write.table(lb3_snp, file.path("~{outdir}", "snp_group3.txt"),
        row.names = F, quote = F, col.names = F)
        write.table(lb4_snp, file.path("~{outdir}", "snp_group4.txt"),
        row.names = F, quote = F, col.names = F)
        CODE
        }

    output {
        Array[File] snp_groups=[
            "~{outdir}/snp_group1.txt",
            "~{outdir}/snp_group2.txt",
            "~{outdir}/snp_group3.txt",
            "~{outdir}/snp_group4.txt"
        ]
    }

    runtime {
        cpu: 1
        memory: "30 GB"
        singularity: "/gpfs/commons/groups/singh_lab/software/docker/r-base_latest.sif"
    }

}

task ldms_grm {
    input {
        String bfile
        File snp
        String outdir
    }

    command {
        if [ -e "~{outdir}/~{basename(bfile)}_~{basename(snp)}.grm.bin" ]; then
            echo "GRM file already exists"
        else
            gcta --bfile ~{bfile} --extract ~{snp} --autosome --make-grm \
            --out ~{outdir}/~{basename(bfile)}_~{basename(snp)} \
            --thread-num 4
        fi
    }

    runtime {
        cpu: 4
        memory: "30 GB"
        #singularity: "docker://gcr.io/singh-comp-d-271c/genetics:0.2.18"
        singularity: "/gpfs/commons/groups/singh_lab/software/docker/genetics_0.2.18.sif"
    }

    output {
        File grm = "~{outdir}/~{basename(bfile)}_~{basename(snp)}.grm.bin"
    }
}
task make_grm {
    input {
        String bfile
        String outdir
    }

    command {
       if [ -e "~{outdir}/~{basename(bfile)}.grm.bin" ]; then
            echo "GRM file already exists"
        else 
            gcta --bfile ~{bfile} --autosome --make-grm --out ~{outdir}/~{basename(bfile)} \
            --thread-num 4
        fi
    }

    runtime {
        cpu: 4
        memory: "30 GB"
        #singularity: "docker://gcr.io/singh-comp-d-271c/genetics:0.2.18"
        singularity: "/gpfs/commons/groups/singh_lab/software/docker/genetics_0.2.18.sif"
    }

    output {
        File grm = "~{outdir}/~{basename(bfile)}.grm.bin"
    }
}