
rule input_reads:
    input:
        reads=get_input_reads,
    output:
        reads=temp("temp/{sm}/{sm}.{read_type}.{idx}.reads.fa.gz"),
    threads: 8
    resources:
        mem_mb=8 * 1024,
        runtime=60 * 4,
    conda:
        "../envs/env.yml"
    shell:
        """
        # make a fasta if bam cram or sam
        if [[ {input.reads} =~ .*\\.(fa|fa.gz|fq|fq.gz|fasta|fasta.gz|fastq|fastq.gz) ]]; then
            echo "linking {input.reads} to {output.reads}"
            ln -s $(realpath {input.reads}) {output.reads}
        elif [[ {input.reads} =~ .*\\.(bam|sam|cram) ]]; then
            echo "converting {input.reads} to {output.reads}"
            samtools fasta -@ {threads} {input.reads} | bgzip -@ {threads} > {output.reads}
        fi
        """


rule merge_input_reads:
    input:
        reads=get_inputs_to_merge,
    output:
        reads=temp("temp/{sm}/{sm}.{read_type}.reads.fa.gz"),
    threads: 8
    resources:
        mem_mb=8 * 1024,
        runtime=60 * 4,
    conda:
        "../envs/env.yml"
    shell:
        """
        cat {input} > {output.reads}
        """


rule yak:
    input:
        parental_reads=get_parental_reads,
    output:
        yak=temp("temp/{sm}/{sm}.{parental}.yak"),
    threads: 16
    resources:
        mem_mb=100 * 1024,
        runtime=60 * 4,
    conda:
        "../envs/env.yml"
    shell:
        """
        yak count -k31 -b37 -t {threads} -o {output.yak} {input.parental_reads} {input.parental_reads}
        """


rule hifiasm:
    input:
        unpack(asm_inputs),
    output:
        hap1="results/{sm}/{sm}.{asm_type}.hap1.p_ctg.gfa",
        lowQ1="results/{sm}/{sm}.{asm_type}.hap1.p_ctg.lowQ.bed",
        noseq1=temp("results/{sm}/{sm}.{asm_type}.hap1.p_ctg.noseq.gfa"),
        hap2="results/{sm}/{sm}.{asm_type}.hap2.p_ctg.gfa",
        lowQ2="results/{sm}/{sm}.{asm_type}.hap2.p_ctg.lowQ.bed",
        noseq2=temp("results/{sm}/{sm}.{asm_type}.hap2.p_ctg.noseq.gfa"),
        utg="results/{sm}/{sm}.{asm_type}.p_utg.gfa",
        lowQutg="results/{sm}/{sm}.{asm_type}.p_utg.lowQ.bed",
        nosequtg=temp("results/{sm}/{sm}.{asm_type}.p_utg.noseq.gfa"),
        r_utg="results/{sm}/{sm}.{asm_type}.r_utg.gfa",
        lowQr_utg="results/{sm}/{sm}.{asm_type}.r_utg.lowQ.bed",
        noseqr_utg=temp("results/{sm}/{sm}.{asm_type}.r_utg.noseq.gfa"),
        #ec_bin="results/{sm}/{sm}.ec.bin",
    threads: config.get("asm-threads", 48)
    resources:
        mem_mb=asm_mem_mb,
        runtime=60 * 24,
    params:
        extra=extra_asm_options,
    conda:
        "../envs/env.yml"
    shell:
        """
        out_dir=results/{wildcards.sm}/{wildcards.sm}
        mkdir -p $out_dir
        hifiasm -o $out_dir -t {threads} {input.reads} {params.extra}
        """


rule gfa_to_fa:
    """Convert GFA to FASTA with minimal PanSN-spec headers (no metadata yet)."""
    input:
        gfa="results/{sm}/{sm}.{asm_type}.{hap}.p_ctg.gfa",
    output:
        fa=temp("temp/{sm}/{sm}.{asm_type}.{hap}.raw.fa.gz"),
        fai=temp("temp/{sm}/{sm}.{asm_type}.{hap}.raw.fa.gz.fai"),
    threads: 4
    resources:
        mem_mb=8 * 1024,
        runtime=60 * 4,
    params:
        hap_num=get_haplotype_number,
    conda:
        "../envs/env.yml"
    shell:
        """
        # Convert GFA to FASTA with basic PanSN-spec headers
        # Format: sample#haplotype#contig (metadata added after orientation)
        gfatools gfa2fa {input.gfa} \
            | awk -v sample="{wildcards.sm}" \
                  -v hap="{params.hap_num}" \
                  'BEGIN {{{{OFS=""}}}}
                   /^>/ {{{{
                       contig = substr($1, 2)
                       print ">" sample "#" hap "#" contig
                       next
                   }}}}
                   {{{{print}}}}' \
            | bgzip -@ {threads} > {output.fa}
        samtools faidx {output.fa}
        """


rule finalize_assembly_no_ref:
    """Finalize assembly with metadata when no reference is provided."""
    input:
        fa=rules.gfa_to_fa.output.fa,
    output:
        fa="results/assemblies/{sm}.{asm_type}.{hap}.fa.gz",
        fai="results/assemblies/{sm}.{asm_type}.{hap}.fa.gz.fai",
    threads: 4
    resources:
        mem_mb=8 * 1024,
        runtime=60,
    params:
        hap_num=get_haplotype_number,
        phasing=get_phasing_description,
        ultralong=get_ultralong_used,
        hifiasm_version=HIFIASM_VERSION,
    conda:
        "../envs/env.yml"
    shell:
        """
        # Add metadata to headers without orientation info
        zcat {input.fa} \
            | awk -v version="{params.hifiasm_version}" \
                  -v phasing="{params.phasing}" \
                  -v ultralong="{params.ultralong}" \
                  '/^>/ {{{{
                       # Add metadata to existing PanSN header
                       print $0 " assembler=hifiasm version=" version " phasing=" phasing " ultralong=" ultralong
                       next
                   }}}}
                   {{{{print}}}}' \
            | bgzip -@ {threads} > {output.fa}
        samtools faidx {output.fa}
        """


# Quick alignment with mashmap for orientation detection
rule quick_align:
    """Fast whole-genome alignment using mashmap for orientation detection."""
    input:
        fa=rules.gfa_to_fa.output.fa,
        ref=get_ref,
    output:
        paf=temp("temp/{sm}/{ref}/{sm}.{asm_type}.{hap}.mashmap.paf"),
    threads: 8
    resources:
        mem_mb=16 * 1024,
        runtime=60,
    conda:
        "../envs/env.yml"
    shell:
        """
        mashmap -r {input.ref} -q {input.fa} \
            -t {threads} --pi 95 -s 50000 \
            -o {output.paf}
        """


rule analyze_orientation:
    """Analyze mashmap PAF to determine contig orientation relative to reference."""
    input:
        paf=rules.quick_align.output.paf,
    output:
        tsv=temp("temp/{sm}/{ref}/{sm}.{asm_type}.{hap}.orientation.tsv"),
    threads: 1
    resources:
        mem_mb=4 * 1024,
        runtime=60,
    conda:
        "../envs/env.yml"
    shell:
        """
        python {workflow.basedir}/scripts/orient_contigs.py {input.paf} {output.tsv}
        """


rule orient_assembly:
    """Orient contigs, add all header metadata, and flip sequences as needed."""
    input:
        fa=rules.gfa_to_fa.output.fa,
        orientation=rules.analyze_orientation.output.tsv,
    output:
        fa="results/assemblies/{ref}/{sm}.{asm_type}.{hap}.fa.gz",
        fai="results/assemblies/{ref}/{sm}.{asm_type}.{hap}.fa.gz.fai",
    threads: 4
    resources:
        mem_mb=8 * 1024,
        runtime=60 * 2,
    params:
        phasing=get_phasing_description,
        ultralong=get_ultralong_used,
        hifiasm_version=HIFIASM_VERSION,
    conda:
        "../envs/env.yml"
    shell:
        """
        python {workflow.basedir}/scripts/apply_orientation.py \
            {input.fa} {input.orientation} /dev/stdout \
            --assembler hifiasm \
            --version "{params.hifiasm_version}" \
            --phasing {params.phasing} \
            --ultralong {params.ultralong} \
            | bgzip -@ {threads} > {output.fa}
        samtools faidx {output.fa}
        """


# Full alignment with minimap2 on oriented assemblies
rule align:
    """Align oriented assembly to reference using minimap2."""
    input:
        fa=rules.orient_assembly.output.fa,
        ref=get_ref,
    output:
        bam="results/alignments/{ref}/{sm}.{asm_type}.{hap}.bam",
        index="results/alignments/{ref}/{sm}.{asm_type}.{hap}.bam.csi",
    threads: 16
    resources:
        mem_mb=64 * 1024,
        runtime=60 * 4,
    conda:
        "../envs/env.yml"
    params:
        mm2_opts=config.get("mm2_opts", "-x asm20 --secondary=no -s 25000 -K 8G"),
    shell:
        """
        minimap2 --MD --cs --eqx -a {params.mm2_opts} \
             {input.ref} {input.fa} \
            | samtools view -F 4 -u -@ {threads} \
            | samtools sort -m 2G -@ {threads} \
                --write-index -o {output.bam}
        """


rule bam_to_paf:
    """Convert BAM alignment to PAF format."""
    input:
        bam=rules.align.output.bam,
    output:
        paf="results/alignments/{ref}/{sm}.{asm_type}.{hap}.paf",
    threads: 4
    resources:
        mem_mb=16 * 1024,
        runtime=60 * 4,
    conda:
        "../envs/env.yml"
    shell:
        """
        samtools view -h -@ {threads} {input.bam} \
            | paftools.js sam2paf - \
            > {output.paf}
        """
