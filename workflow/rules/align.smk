
# Full alignment with minimap2 on oriented assemblies
rule align:
    """Align oriented assembly to reference using minimap2."""
    input:
        fa=lambda wc: rules.orient_assembly.output.fa.format(sm=wc.sm, asm_type=wc.asm_type, hap=wc.hap),
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
        rg=lambda wc: f"@RG\\tID:{wc.sm}\\tSM:{wc.sm}",
    shell:
        """
        minimap2 --MD --cs --eqx -a {params.mm2_opts} \
            -R '{params.rg}' \
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


rule merge_haplotype_fasta:
    """Merge hap1 and hap2 FASTA files into a single diploid assembly."""
    input:
        hap1=lambda wc: rules.orient_assembly.output.fa.format(sm=wc.sm, asm_type=wc.asm_type, hap="hap1"),
        hap2=lambda wc: rules.orient_assembly.output.fa.format(sm=wc.sm, asm_type=wc.asm_type, hap="hap2"),
    output:
        fa="results/assemblies/{sm}.{asm_type}.dip.fa.gz",
        fai="results/assemblies/{sm}.{asm_type}.dip.fa.gz.fai",
    threads: 4
    resources:
        mem_mb=8 * 1024,
        runtime=60,
    conda:
        "../envs/env.yml"
    shell:
        """
        cat {input.hap1} {input.hap2} > {output.fa}
        samtools faidx {output.fa}
        """


rule merge_haplotype_bam:
    """Merge hap1 and hap2 BAM files into a single diploid alignment."""
    input:
        hap1=lambda wc: rules.align.output.bam.format(ref=wc.ref, sm=wc.sm, asm_type=wc.asm_type, hap="hap1"),
        hap2=lambda wc: rules.align.output.bam.format(ref=wc.ref, sm=wc.sm, asm_type=wc.asm_type, hap="hap2"),
    output:
        bam="results/alignments/{ref}/{sm}.{asm_type}.dip.bam",
        index="results/alignments/{ref}/{sm}.{asm_type}.dip.bam.csi",
    threads: 8
    resources:
        mem_mb=16 * 1024,
        runtime=60 * 2,
    conda:
        "../envs/env.yml"
    shell:
        """
        samtools merge -@ {threads} --write-index -o {output.bam} {input.hap1} {input.hap2}
        """
