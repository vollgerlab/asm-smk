
def get_ref_url(wc):
    """Get URL for reference download."""
    return DEFAULT_REFS.get(wc.ref_name, {}).get("url", "")


rule download_reference:
    """Download and decompress reference genome from URL."""
    output:
        fa="resources/references/{ref_name}.fa",
        fai="resources/references/{ref_name}.fa.fai",
    params:
        url=get_ref_url,
    threads: 4
    resources:
        mem_mb=8 * 1024,
        runtime=60 * 2,
    conda:
        "../envs/env.yml"
    shell:
        """
        mkdir -p resources/references
        curl -L {params.url} | gunzip -c > {output.fa}
        samtools faidx {output.fa}
        """


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
    """Run hifiasm assembly. Uses shadow to auto-cleanup intermediate files."""
    input:
        unpack(asm_inputs),
    output:
        hap1="results/{sm}/{sm}.{asm_type}.hap1.p_ctg.gfa",
        hap2="results/{sm}/{sm}.{asm_type}.hap2.p_ctg.gfa",
    shadow:
        "minimal"
    threads: config.get("asm-threads", 64)
    resources:
        mem_mb=asm_mem_mb,
        runtime=60 * 24,
    params:
        extra=extra_asm_options,
        prefix=lambda wc, output: output.hap1.replace(".hap1.p_ctg.gfa", ""),
    conda:
        "../envs/env.yml"
    shell:
        """
        hifiasm -o {params.prefix} -t {threads} {input.reads} {params.extra}
        """


rule gfa_to_fa:
    """Convert GFA to FASTA (raw contig names, headers added in orient_assembly)."""
    input:
        gfa="results/{sm}/{sm}.{asm_type}.{hap}.p_ctg.gfa",
    output:
        fa=temp("temp/{sm}/{sm}.{asm_type}.{hap}.raw.fa.gz"),
        fai=temp("temp/{sm}/{sm}.{asm_type}.{hap}.raw.fa.gz.fai"),
    threads: 4
    resources:
        mem_mb=8 * 1024,
        runtime=60 * 4,
    conda:
        "../envs/env.yml"
    shell:
        """
        gfatools gfa2fa {input.gfa} | bgzip -@ {threads} > {output.fa}
        samtools faidx {output.fa}
        """


# Quick alignment with minimap2 for orientation detection (no base-level alignment)
rule quick_align:
    """Fast whole-genome alignment using minimap2 for orientation detection."""
    input:
        fa=rules.gfa_to_fa.output.fa,
        ref=get_ref,
    output:
        paf=temp("temp/{sm}/{ref}/{sm}.{asm_type}.{hap}.orient.paf"),
    threads: 8
    resources:
        mem_mb=16 * 1024,
        runtime=60,
    conda:
        "../envs/env.yml"
    shell:
        """
        minimap2 -s 25000 -x asm20 -t {threads} --secondary=no \
            {input.ref} {input.fa} > {output.paf}
        """


rule analyze_orientation:
    """Analyze PAF to determine contig orientation relative to reference."""
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
        hap_num=get_haplotype_number,
        phasing=get_phasing_description,
        ultralong=get_ultralong_used,
        hifiasm_version=HIFIASM_VERSION,
    conda:
        "../envs/env.yml"
    shell:
        """
        python {workflow.basedir}/scripts/apply_orientation.py \
            {input.fa} {input.orientation} /dev/stdout \
            --sample {wildcards.sm} \
            --haplotype {params.hap_num} \
            --assembler hifiasm \
            --version "{params.hifiasm_version}" \
            --phasing {params.phasing} \
            --ultralong {params.ultralong} \
            --orient-ref {wildcards.ref} \
            | bgzip -@ {threads} > {output.fa}
        samtools faidx {output.fa}
        """
