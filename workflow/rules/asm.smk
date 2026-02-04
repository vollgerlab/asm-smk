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
    log:
        "results/logs/download_reference/{ref_name}.log",
    conda:
        "../envs/env.yml"
    shell:
        "(curl -L {params.url} | gunzip -c > {output.fa} && "
        "samtools faidx {output.fa}) &> {log}"


rule prepare_uncompressed_ref:
    """Prepare uncompressed reference for tools that cannot read gzipped FASTA (e.g., svim-asm)."""
    input:
        ref=get_ref,
    output:
        fa=temp("temp/references/{ref}.fa"),
        fai=temp("temp/references/{ref}.fa.fai"),
    threads: 4
    resources:
        mem_mb=8 * 1024,
        runtime=60,
    log:
        "results/logs/prepare_uncompressed_ref/{ref}.log",
    conda:
        "../envs/env.yml"
    shell:
        "(if [[ {input.ref} == *.gz ]]; then "
        "    gunzip -c {input.ref} > {output.fa}; "
        "else "
        "    ln -s $(realpath {input.ref}) {output.fa}; "
        "fi && "
        "samtools faidx {output.fa}) &> {log}"


rule input_reads:
    input:
        reads=get_input_reads,
    output:
        reads=temp("temp/{sm}/{sm}.{read_type}.{idx}.reads.fa.gz"),
    threads: 8
    resources:
        mem_mb=8 * 1024,
        runtime=60 * 4,
    log:
        "results/logs/input_reads/{sm}.{read_type}.{idx}.log",
    conda:
        "../envs/env.yml"
    shell:
        "(if [[ {input.reads} =~ .*\\.(fa|fa.gz|fq|fq.gz|fasta|fasta.gz|fastq|fastq.gz) ]]; then "
        "    echo 'linking {input.reads} to {output.reads}' && "
        "    ln -s $(realpath {input.reads}) {output.reads}; "
        "elif [[ {input.reads} =~ .*\\.(bam|sam|cram) ]]; then "
        "    echo 'converting {input.reads} to {output.reads}' && "
        "    samtools fasta -@ {threads} {input.reads} | bgzip -@ {threads} > {output.reads}; "
        "fi) &> {log}"


rule merge_input_reads:
    input:
        reads=get_inputs_to_merge,
    output:
        reads=temp("temp/{sm}/{sm}.{read_type}.reads.fa.gz"),
    threads: 8
    resources:
        mem_mb=8 * 1024,
        runtime=60 * 4,
    log:
        "results/logs/merge_input_reads/{sm}.{read_type}.log",
    conda:
        "../envs/env.yml"
    shell:
        "cat {input} > {output.reads} 2> {log}"


rule yak:
    input:
        parental_reads=get_parental_reads,
    output:
        yak=temp("temp/{sm}/{sm}.{parental}.yak"),
    threads: 16
    resources:
        mem_mb=100 * 1024,
        runtime=60 * 4,
    log:
        "results/logs/yak/{sm}.{parental}.log",
    conda:
        "../envs/env.yml"
    shell:
        "yak count -k31 -b37 -t {threads} -o {output.yak} "
        "{input.parental_reads} {input.parental_reads} &> {log}"


rule hifiasm:
    """Run hifiasm assembly. Uses shadow to auto-cleanup intermediate files."""
    input:
        unpack(asm_inputs),
    output:
        hap1=temp("temp/{sm}/{sm}.{asm_type}.hap1.p_ctg.gfa"),
        hap2=temp("temp/{sm}/{sm}.{asm_type}.hap2.p_ctg.gfa"),
    shadow:
        "minimal"
    threads: config.get("asm-threads", 64)
    resources:
        mem_mb=asm_mem_mb,
        runtime=60 * 24,
    log:
        "results/logs/hifiasm/{sm}.{asm_type}.log",
    params:
        extra=extra_asm_options,
        # hifiasm adds asm_type (.bp/.dip/.hic) to prefix automatically
        prefix=lambda wc, output: output.hap1.replace(
            f".{wc.asm_type}.hap1.p_ctg.gfa", ""
        ),
    conda:
        "../envs/env.yml"
    shell:
        "hifiasm -o {params.prefix} -t {threads} {input.reads} {params.extra} &> {log}"


rule gfa_to_fa:
    """Convert GFA to FASTA (raw contig names, headers added in orient_assembly)."""
    input:
        gfa="temp/{sm}/{sm}.{asm_type}.{hap}.p_ctg.gfa",
    output:
        fa=temp("temp/{sm}/{sm}.{asm_type}.{hap}.raw.fa.gz"),
        fai=temp("temp/{sm}/{sm}.{asm_type}.{hap}.raw.fa.gz.fai"),
    threads: 4
    resources:
        mem_mb=8 * 1024,
        runtime=60 * 4,
    log:
        "results/logs/gfa_to_fa/{sm}.{asm_type}.{hap}.log",
    conda:
        "../envs/env.yml"
    shell:
        "(gfatools gfa2fa {input.gfa} | bgzip -@ {threads} > {output.fa} && "
        "samtools faidx {output.fa}) &> {log}"


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
    log:
        "results/logs/quick_align/{sm}.{ref}.{asm_type}.{hap}.log",
    conda:
        "../envs/env.yml"
    shell:
        "minimap2 -s 25000 -x asm20 -t {threads} --secondary=no "
        "{input.ref} {input.fa} > {output.paf} 2> {log}"


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
    log:
        "results/logs/analyze_orientation/{sm}.{ref}.{asm_type}.{hap}.log",
    conda:
        "../envs/env.yml"
    shell:
        "python {workflow.basedir}/scripts/orient_contigs.py "
        "{input.paf} {output.tsv} &> {log}"


rule orient_assembly:
    """Orient contigs, add all header metadata, and flip sequences as needed."""
    input:
        fa=rules.gfa_to_fa.output.fa,
        orientation=lambda wc: f"temp/{wc.sm}/{ORIENT_REF}/{wc.sm}.{wc.asm_type}.{wc.hap}.orientation.tsv",
    output:
        fa="results/assemblies/{sm}.{asm_type}.{hap}.fa.gz",
        fai="results/assemblies/{sm}.{asm_type}.{hap}.fa.gz.fai",
    threads: 4
    resources:
        mem_mb=8 * 1024,
        runtime=60 * 2,
    log:
        "results/logs/orient_assembly/{sm}.{asm_type}.{hap}.log",
    params:
        hap_num=get_haplotype_number,
        phasing=get_phasing_description,
        ultralong=get_ultralong_used,
        hifiasm_version=HIFIASM_VERSION,
        orient_ref=ORIENT_REF,
    conda:
        "../envs/env.yml"
    shell:
        "(python {workflow.basedir}/scripts/apply_orientation.py "
        "    {input.fa} {input.orientation} /dev/stdout "
        "    --sample {wildcards.sm} "
        "    --haplotype {params.hap_num} "
        "    --assembler hifiasm "
        '    --version "{params.hifiasm_version}" '
        "    --phasing {params.phasing} "
        "    --ultralong {params.ultralong} "
        "    --orient-ref {params.orient_ref} "
        "    | bgzip -@ {threads} > {output.fa} && "
        "samtools faidx {output.fa}) &> {log}"


rule merge_haplotype_fasta:
    """Merge hap1 and hap2 FASTA files into a single diploid assembly."""
    input:
        hap1=lambda wc: rules.orient_assembly.output.fa.format(
            sm=wc.sm, asm_type=wc.asm_type, hap="hap1"
        ),
        hap2=lambda wc: rules.orient_assembly.output.fa.format(
            sm=wc.sm, asm_type=wc.asm_type, hap="hap2"
        ),
    output:
        fa="results/assemblies/{sm}.{asm_type}.dip.fa.gz",
        fai="results/assemblies/{sm}.{asm_type}.dip.fa.gz.fai",
    threads: 4
    resources:
        mem_mb=8 * 1024,
        runtime=60,
    log:
        "results/logs/merge_haplotype_fasta/{sm}.{asm_type}.log",
    conda:
        "../envs/env.yml"
    shell:
        "(cat {input.hap1} {input.hap2} > {output.fa} && "
        "samtools faidx {output.fa}) &> {log}"
