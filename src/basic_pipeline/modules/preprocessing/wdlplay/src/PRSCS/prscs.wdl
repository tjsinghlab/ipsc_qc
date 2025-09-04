version 1.0

workflow PRSCS {
    input {
        File raw_sumstats
        String tag
        String bim
        Int sample_size
        String ld_pop="eur" # eur or afr
        String outdir
        Int? n_cpu
        Int? n_iter
        Int? n_burnin
        Int? thinning
        Int export_mem=8
    }

    call export_prscs_input {
        input: 
            raw_sumstats=raw_sumstats,
            tag=tag,
            mem=export_mem
    }

    call clean_tsv {
        input:
            input_sumstats=export_prscs_input.input_sumstats,
            tag=tag
    }

    scatter (chr in range(22)) {
        call run_prscs {
            input: 
                bim=bim,
                sumstats=clean_tsv.cleaned_sumstats,
                chromosome=chr+1,
                sample_size=sample_size,
                ld_pop=ld_pop,
                tag=tag,
                n_cpu=n_cpu,
                n_iter=n_iter,
                n_burnin=n_burnin,
                thinning=thinning
        }
    }

    call combine_outputs {
        input:
            prscs_outputs=run_prscs.prscs_output,
            tag=tag,
            outdir=outdir
    }

    output {
        File input_sumstats = export_prscs_input.input_sumstats
        Array[File] prscs_outputs = run_prscs.prscs_output
        File combined_output = combine_outputs.combined_output
    }
}

task export_prscs_input {
    input {
        File raw_sumstats
        String tag
        Int mem
    }

    command {
        python3 -c "
        import pandas as pd
        data = pd.read_table('~{raw_sumstats}', compression='gzip', low_memory=False)
        data = data[['SNP', 'A1', 'A2', 'BETA', 'P']]
        data.to_csv('~{tag}_prscs.tsv', sep='\t', header=True, index=False)
        "
    }

    runtime {
        cpu: 1
        memory: "~{mem} GB"
        singularity: "/gpfs/commons/groups/singh_lab/resources/images/hailgenetics-python-dill-pandas-batchutils_latest.sif"
    }

    output {
        File input_sumstats = "~{tag}_prscs.tsv"
    }
}

task clean_tsv {
    input {
        File input_sumstats
        String tag
    }

    command {
        grep -v -e '^[[:space:]]*$' ~{input_sumstats} > ~{tag}_prscs_cleaned.tsv
    }

    runtime {
        cpu:1
        memory: "8GB"
    }

    output {
        File cleaned_sumstats = "~{tag}_prscs_cleaned.tsv"
    }
}

task run_prscs {
    input {
        String bim
        File sumstats
        Int chromosome
        Int sample_size
        String ld_pop
        String tag
        Int n_cpu = 1
        Int n_iter = 1000
        Int n_burnin = 500
        Int thinning = 5
    }

    command {
        python /PRScs/PRScs.py \
            --ref_dir=/PRScs/ldblk_ukbb_~{ld_pop} \
            --bim_prefix=~{sub(bim,"\.bim$","")} \
            --sst_file=~{sumstats} \
            --n_gwas=~{sample_size} \
            --n_iter=~{n_iter} \
            --n_burnin=~{n_burnin} \
            --thin=~{thinning} \
            --chrom=~{chromosome} \
            --seed=827 \
            --out_dir=output 
    }

    runtime {
        cpu: n_cpu
        memory: "16 GB"
        singularity: "/gpfs/commons/groups/singh_lab/resources/images/prscs_latest.sif"
    }

    output {
        File prscs_output = "output_pst_eff_a1_b0.5_phiauto_chr~{chromosome}.txt"
    }
}

task combine_outputs {
    input {
        Array[File] prscs_outputs
        String tag
        String outdir
    }

    command {
        # cat ~{sep=" " prscs_outputs} > prscs_outputs_~{tag}.txt
        mkdir -p ~{outdir}
        for file in ~{sep=" " prscs_outputs}; do
            echo "$file"
            cat "$file" >> ~{outdir}prscs_outputs_~{tag}.txt
        done
    }

    runtime {
        cpu:1
        memory: "8 GB"
    }

    output {
        File combined_output = "~{outdir}prscs_outputs_~{tag}.txt"
    }
}