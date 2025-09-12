version 1.0

task star {

    input {
        File fastq1
        File? fastq2
        String prefix
        String star_index

        # STAR options
        Int? outFilterMultimapNmax
        Int? alignSJoverhangMin
        Int? alignSJDBoverhangMin
        Int? outFilterMismatchNmax
        Float? outFilterMismatchNoverLmax
        Int? alignIntronMin
        Int? alignIntronMax
        Int? alignMatesGapMax
        String? outFilterType
        Float? outFilterScoreMinOverLread
        Float? outFilterMatchNminOverLread
        Int? limitSjdbInsertNsj
        String? outSAMstrandField
        String? outFilterIntronMotifs
        String? alignSoftClipAtReferenceEnds
        String? quantMode
        String? outSAMattrRGline
        String? outSAMattributes
        File? varVCFfile
        String? waspOutputMode
        Int? chimSegmentMin
        Int? chimJunctionOverhangMin
        String? chimOutType
        Int? chimMainSegmentMultNmax
        Int? chimOutJunctionFormat
        File? sjdbFileChrStartEnd

        String outdir

        Int memory=64
        Int disk_space=200
        Int num_threads=8
        Int num_preempt=0
    }

    command {
        set -euo pipefail

        if [[ ${fastq1} == *".tar" || ${fastq1} == *".tar.gz" ]]; then
            tar -xvvf ${fastq1}
            fastq1_abs=$(for f in *_1.fastq*; do echo "$(pwd)/$f"; done | paste -s -d ',')
            fastq2_abs=$(for f in *_2.fastq*; do echo "$(pwd)/$f"; done | paste -s -d ',')
            if [[ $fastq1_abs == *"*_1.fastq*" ]]; then  # no paired-end FASTQs found; check for single-end FASTQ
                fastq1_abs=$(for f in *.fastq*; do echo "$(pwd)/$f"; done | paste -s -d ',')
                fastq2_abs=''
            fi
        else
            # make sure paths are absolute
            fastq1_abs=${fastq1}
            fastq2_abs=${fastq2}
            if [[ $fastq1_abs != /* ]]; then
                fastq1_abs=$PWD/$fastq1_abs
                fastq2_abs=$PWD/$fastq2_abs
            fi
        fi

        echo "FASTQs:"
        echo $fastq1_abs
        echo $fastq2_abs

        # extract index
        # echo $(date +"[%b %d %H:%M:%S] Extracting STAR index")
        # mkdir star_index
        # tar -xvvf ${star_index} -C star_index --strip-components=1

        mkdir -p ${outdir + "star_out"}
        # placeholders for optional outputs
        touch ${outdir + "star_out"}/${prefix}.Aligned.toTranscriptome.out.bam
        touch ${outdir + "star_out"}/${prefix}.Chimeric.out.sorted.bam
        touch ${outdir + "star_out"}/${prefix}.Chimeric.out.sorted.bam.bai
        touch ${outdir + "star_out"}/${prefix}.ReadsPerGene.out.tab  # run_STAR.py will gzip

        /src/run_STAR.py \
            ${star_index} $fastq1_abs $fastq2_abs ${prefix} \
            --output_dir ${outdir + "star_out"} \
            ${"--outFilterMultimapNmax " + outFilterMultimapNmax} \
            ${"--alignSJoverhangMin " + alignSJoverhangMin} \
            ${"--alignSJDBoverhangMin " + alignSJDBoverhangMin} \
            ${"--outFilterMismatchNmax " + outFilterMismatchNmax} \
            ${"--outFilterMismatchNoverLmax " + outFilterMismatchNoverLmax} \
            ${"--alignIntronMin " + alignIntronMin} \
            ${"--alignIntronMax " + alignIntronMax} \
            ${"--alignMatesGapMax " + alignMatesGapMax} \
            ${"--outFilterType " + outFilterType} \
            ${"--outFilterScoreMinOverLread " + outFilterScoreMinOverLread} \
            ${"--outFilterMatchNminOverLread " + outFilterMatchNminOverLread} \
            ${"--limitSjdbInsertNsj " + limitSjdbInsertNsj} \
            ${"--outSAMstrandField " + outSAMstrandField} \
            ${"--outFilterIntronMotifs " + outFilterIntronMotifs} \
            ${"--alignSoftClipAtReferenceEnds " + alignSoftClipAtReferenceEnds} \
            ${"--quantMode " + quantMode} \
            ${"--outSAMattrRGline " + outSAMattrRGline} \
            ${"--outSAMattributes " + outSAMattributes} \
            ${"--varVCFfile " + varVCFfile} \
            ${"--waspOutputMode " + waspOutputMode} \
            ${"--chimSegmentMin " + chimSegmentMin} \
            ${"--chimJunctionOverhangMin " + chimJunctionOverhangMin} \
            ${"--chimOutType " + chimOutType} \
            ${"--chimMainSegmentMultNmax " + chimMainSegmentMultNmax} \
            ${"--chimOutJunctionFormat " + chimOutJunctionFormat} \
            ${"--sjdbFileChrStartEnd " + sjdbFileChrStartEnd} \
            --threads ${num_threads}
    }

    output {
        File bam_file = "${outdir}/star_out/${prefix}.Aligned.sortedByCoord.out.bam"
        File bam_index = "${outdir}/star_out/${prefix}.Aligned.sortedByCoord.out.bam.bai"
        File transcriptome_bam = "${outdir}/star_out/${prefix}.Aligned.toTranscriptome.out.bam"
        File chimeric_junctions = "${outdir}/star_out/${prefix}.Chimeric.out.junction.gz"
        File chimeric_bam_file = "${outdir}/star_out/${prefix}.Chimeric.out.sorted.bam"
        File chimeric_bam_index = "${outdir}/star_out/${prefix}.Chimeric.out.sorted.bam.bai"
        File read_counts = "${outdir}/star_out/${prefix}.ReadsPerGene.out.tab.gz"
        File junctions = "${outdir}/star_out/${prefix}.SJ.out.tab.gz"
        File junctions_pass1 = "${outdir}/star_out/${prefix}._STARpass1/${prefix}.SJ.pass1.out.tab.gz"
        Array[File] logs = ["${outdir}/star_out/${prefix}.Log.final.out", "${outdir}/star_out/${prefix}.Log.out", "${outdir}/star_out/${prefix}.Log.progress.out"]
    }

    runtime {
        singularity: "/ref/gtex_rnaseq_V10.sif"
        memory: "${memory}GB"
        disks: "local-disk ${disk_space} HDD"
        cpu: "${num_threads}"
        preemptible: "${num_preempt}"
    }

    meta {
        author: "Francois Aguet"
    }
}