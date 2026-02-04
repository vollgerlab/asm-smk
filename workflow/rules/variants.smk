# Variant calling from assembly alignments


rule call_small_variants:
    """Call SNPs and small indels from assembly alignment using paftools.js."""
    input:
        bam=rules.align.output.bam,
        ref=get_ref,
    output:
        vcf=temp("temp/{sm}/{ref}/{sm}.{asm_type}.{hap}.small.vcf.gz"),
        tbi=temp("temp/{sm}/{ref}/{sm}.{asm_type}.{hap}.small.vcf.gz.tbi"),
    threads: 4
    resources:
        mem_mb=8 * 1024,
        runtime=60 * 2,
    log:
        "results/logs/call_small_variants/{ref}/{sm}.{asm_type}.{hap}.log",
    conda:
        "../envs/env.yml"
    shell:
        "(samtools view -h -@ {threads} {input.bam} "
        "    | paftools.js call -f {input.ref} -L50000 -l10000 - "
        "    | bgzip -@ {threads} > {output.vcf} && "
        "tabix -p vcf {output.vcf}) &> {log}"


rule phase_small_variants:
    """Convert haplotype-specific variant calls to phased diploid VCF."""
    input:
        hap1=lambda wc: f"temp/{wc.sm}/{wc.ref}/{wc.sm}.{wc.asm_type}.hap1.small.vcf.gz",
        hap2=lambda wc: f"temp/{wc.sm}/{wc.ref}/{wc.sm}.{wc.asm_type}.hap2.small.vcf.gz",
        ref=get_ref,
    output:
        vcf=temp("temp/{sm}/{ref}/{sm}.{asm_type}.small.phased.vcf.gz"),
        tbi=temp("temp/{sm}/{ref}/{sm}.{asm_type}.small.phased.vcf.gz.tbi"),
    threads: 4
    resources:
        mem_mb=8 * 1024,
        runtime=60 * 2,
    log:
        "results/logs/phase_small_variants/{ref}/{sm}.{asm_type}.log",
    conda:
        "../envs/env.yml"
    shell:
        "(python {workflow.basedir}/scripts/phase_haplotype_vcfs.py "
        "    --hap1 {input.hap1} "
        "    --hap2 {input.hap2} "
        "    --sample {wildcards.sm} "
        "    --ref {input.ref} "
        "    --output /dev/stdout "
        "    | bgzip -@ {threads} > {output.vcf} && "
        "tabix -p vcf {output.vcf}) &> {log}"


rule filter_small_variants:
    """Filter small variants to exclude those >= min_sv_size (handled by SVIM-asm)."""
    input:
        vcf=lambda wc: f"temp/{wc.sm}/{wc.ref}/{wc.sm}.{wc.asm_type}.small.phased.vcf.gz",
    output:
        vcf=temp("temp/{sm}/{ref}/{sm}.{asm_type}.small.filtered.vcf.gz"),
        tbi=temp("temp/{sm}/{ref}/{sm}.{asm_type}.small.filtered.vcf.gz.tbi"),
    params:
        min_sv_size=config.get("min_sv_size", 50),
    threads: 4
    resources:
        mem_mb=4 * 1024,
        runtime=60,
    log:
        "results/logs/filter_small_variants/{ref}/{sm}.{asm_type}.log",
    conda:
        "../envs/env.yml"
    shell:
        "(bcftools view -i 'abs(ILEN) < {params.min_sv_size}' "
        "    -Oz -o {output.vcf} {input.vcf} && "
        "tabix -p vcf {output.vcf}) &> {log}"


rule call_svs:
    """Call structural variants from diploid assembly using SVIM-asm."""
    input:
        hap1=lambda wc: rules.align.output.bam.format(
            ref=wc.ref, sm=wc.sm, asm_type=wc.asm_type, hap="hap1"
        ),
        hap1_idx=lambda wc: rules.align.output.index.format(
            ref=wc.ref, sm=wc.sm, asm_type=wc.asm_type, hap="hap1"
        ),
        hap2=lambda wc: rules.align.output.bam.format(
            ref=wc.ref, sm=wc.sm, asm_type=wc.asm_type, hap="hap2"
        ),
        hap2_idx=lambda wc: rules.align.output.index.format(
            ref=wc.ref, sm=wc.sm, asm_type=wc.asm_type, hap="hap2"
        ),
        ref=rules.prepare_uncompressed_ref.output.fa,
        ref_idx=rules.prepare_uncompressed_ref.output.fai,
    output:
        vcf=temp("temp/{sm}/{ref}/{sm}.{asm_type}.svim-asm.vcf.gz"),
        tbi=temp("temp/{sm}/{ref}/{sm}.{asm_type}.svim-asm.vcf.gz.tbi"),
    params:
        min_sv_size=config.get("min_sv_size", 50),
    shadow:
        "minimal"
    threads: 4
    resources:
        mem_mb=16 * 1024,
        runtime=60 * 4,
    log:
        "results/logs/call_svs/{ref}/{sm}.{asm_type}.log",
    conda:
        "../envs/env.yml"
    shell:
        "(svim-asm diploid svim_workdir {input.hap1} {input.hap2} {input.ref} "
        "    --sample {wildcards.sm} "
        "    --min_sv_size {params.min_sv_size} && "
        "bgzip -@ {threads} -c svim_workdir/variants.vcf > {output.vcf} && "
        "tabix -p vcf {output.vcf}) &> {log}"


rule merge_variants:
    """Merge small variants and SVs into a single phased VCF."""
    input:
        small=lambda wc: f"temp/{wc.sm}/{wc.ref}/{wc.sm}.{wc.asm_type}.small.filtered.vcf.gz",
        small_tbi=lambda wc: f"temp/{wc.sm}/{wc.ref}/{wc.sm}.{wc.asm_type}.small.filtered.vcf.gz.tbi",
        svs=lambda wc: f"temp/{wc.sm}/{wc.ref}/{wc.sm}.{wc.asm_type}.svim-asm.vcf.gz",
        svs_tbi=lambda wc: f"temp/{wc.sm}/{wc.ref}/{wc.sm}.{wc.asm_type}.svim-asm.vcf.gz.tbi",
        ref=get_ref,
    output:
        vcf="results/alignments/{ref}/{sm}.{asm_type}.vcf.gz",
        tbi="results/alignments/{ref}/{sm}.{asm_type}.vcf.gz.tbi",
    threads: 4
    resources:
        mem_mb=8 * 1024,
        runtime=60 * 2,
    log:
        "results/logs/merge_variants/{ref}/{sm}.{asm_type}.log",
    conda:
        "../envs/env.yml"
    shell:
        "(bcftools concat -a -D {input.small} {input.svs} "
        "    | bcftools sort -m 4G "
        "    | bcftools norm -f {input.ref} -c s -m -any "
        "    | bgzip -@ {threads} > {output.vcf} && "
        "tabix -p vcf {output.vcf}) &> {log}"
