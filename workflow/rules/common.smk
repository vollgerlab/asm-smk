# Parse tool versions from env.yml
import yaml
import re

with open(workflow.source_path("../envs/env.yml")) as f:
    _env_config = yaml.safe_load(f)
    _versions = {}
    for dep in _env_config.get("dependencies", []):
        if isinstance(dep, str):
            # Handle version specifiers: ==, >=, <=, >, <, =
            match = re.match(r"^([a-zA-Z0-9_-]+)(?:[><=]+)(.+)$", dep)
            if match:
                name, version = match.groups()
                _versions[name] = version
    HIFIASM_VERSION = _versions.get("hifiasm", "unknown")


def get_mem_mb(wildcards, attempt):
    if attempt < 3:
        return attempt * 1024 * 8
    return attempt * 1024 * 16


# get the memory for the assembly
def asm_mem_mb(wc, attempt):
    # 500GB for the first attempt, 1TB for the second, etc.
    return 500 * 1024 * attempt


def asm_inputs(wc):
    rtn = {}
    rtn["reads"] = rules.merge_input_reads.output.reads.format(
        sm=wc.sm,
        read_type="hifi",
    )
    # Trio phasing inputs
    if HAS_PARENTAL[wc.sm]:
        rtn["pat"] = expand(rules.yak.output.yak, parental="pat", allow_missing=True)
        rtn["mat"] = expand(rules.yak.output.yak, parental="mat", allow_missing=True)
    # Hi-C phasing inputs
    if HAS_HIC[wc.sm]:
        rtn["hic_r1"] = rules.merge_input_reads.output.reads.format(
            sm=wc.sm,
            read_type="hic_r1",
        )
        rtn["hic_r2"] = rules.merge_input_reads.output.reads.format(
            sm=wc.sm,
            read_type="hic_r2",
        )
    # Ultra-long reads
    if HAS_ULTRALONG[wc.sm]:
        rtn["ultralong"] = rules.merge_input_reads.output.reads.format(
            sm=wc.sm,
            read_type="ultralong",
        )
    return rtn


def get_input_reads(wc):
    idx = int(wc.idx)
    if wc.read_type == "hifi":
        return tbl.loc[wc.sm, "hifi"][idx]
    elif wc.read_type == "mat":
        return tbl.loc[wc.sm, "maternal"][idx]
    elif wc.read_type == "pat":
        return tbl.loc[wc.sm, "paternal"][idx]
    elif wc.read_type == "hic_r1":
        return tbl.loc[wc.sm, "hic_r1"][idx]
    elif wc.read_type == "hic_r2":
        return tbl.loc[wc.sm, "hic_r2"][idx]
    elif wc.read_type == "ultralong":
        return tbl.loc[wc.sm, "ultralong"][idx]
    else:
        raise ValueError(f"Unknown read type: {wc.read_type}")


def get_inputs_to_merge(wc):
    # Map short names to column names
    column_map = {
        "pat": "paternal",
        "mat": "maternal",
        "hifi": "hifi",
        "hic_r1": "hic_r1",
        "hic_r2": "hic_r2",
        "ultralong": "ultralong",
    }
    column_name = column_map.get(wc.read_type, wc.read_type)
    files = tbl.loc[wc.sm, column_name]
    if files is None or (isinstance(files, float) and pd.isna(files)):
        return []
    n_files = len(files)
    return expand(
        rules.input_reads.output.reads,
        sm=wc.sm,
        read_type=wc.read_type,
        idx=[i for i in range(n_files)],
    )


def get_parental_reads(wc):
    if wc.parental == "pat":
        return rules.merge_input_reads.output.reads.format(
            sm=wc.sm,
            read_type="pat",
        )
    elif wc.parental == "mat":
        return rules.merge_input_reads.output.reads.format(
            sm=wc.sm,
            read_type="mat",
        )
    else:
        raise ValueError(f"Unknown parental type: {wc.parental}")


def extra_asm_options(wc):
    options = []

    # Phasing options (mutually exclusive - validated in Snakefile)
    if wc.asm_type == "dip":
        pat_yak = rules.yak.output.yak.format(parental="pat", asm_type="dip", sm=wc.sm)
        mat_yak = rules.yak.output.yak.format(parental="mat", asm_type="dip", sm=wc.sm)
        options.append(f"-1 {pat_yak} -2 {mat_yak}")
    elif wc.asm_type == "hic":
        hic_r1 = rules.merge_input_reads.output.reads.format(sm=wc.sm, read_type="hic_r1")
        hic_r2 = rules.merge_input_reads.output.reads.format(sm=wc.sm, read_type="hic_r2")
        options.append(f"--h1 {hic_r1} --h2 {hic_r2}")
    # bp type has no phasing options

    # Ultra-long reads (can combine with any phasing mode)
    if HAS_ULTRALONG[wc.sm]:
        ul_reads = rules.merge_input_reads.output.reads.format(sm=wc.sm, read_type="ultralong")
        options.append(f"--ul {ul_reads}")

    return " ".join(options)


def get_ref(wc):
    return REFS[wc.ref]


def get_haplotype_number(wc):
    """Convert hap1/hap2 to 1/2 for PanSN-spec naming."""
    return "1" if wc.hap == "hap1" else "2"


def get_phasing_description(wc):
    """Get human-readable phasing description for header comments."""
    if wc.asm_type == "dip":
        return "trio"
    elif wc.asm_type == "hic":
        return "hic"
    else:
        return "none"


def get_ultralong_used(wc):
    """Return whether ultra-long ONT reads were used."""
    return "yes" if HAS_ULTRALONG[wc.sm] else "no"


# Default reference URLs for download (ref_name -> url)
DEFAULT_REFS = {
    "T2T-CHM13v2.0": {
        "url": "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz",
    },
    "GRCh38": {
        "url": "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz",
    },
}


def set_references():
    """Set up references for alignment. Returns dict of ref_name -> path."""
    references = config.get("references", {})
    if not references:
        # Check for local paths first, otherwise use download paths
        t2t_local = "/mmfs1/gscratch/stergachislab/assemblies/T2Tv2.0_maskedY.fa"
        hg38_local = "/mmfs1/gscratch/stergachislab/assemblies/simple-names/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna"

        if os.path.exists(t2t_local):
            references["T2T-CHM13v2.0"] = t2t_local
        else:
            references["T2T-CHM13v2.0"] = "resources/references/T2T-CHM13v2.0.fa"

        if os.path.exists(hg38_local):
            references["GRCh38"] = hg38_local
        else:
            references["GRCh38"] = "resources/references/GRCh38.fa"
    return references


def get_orient_ref():
    """Get the single reference used for orientation (first reference, T2T preferred)."""
    refs = set_references()
    if refs:
        # Prefer T2T for orientation
        if "T2T-CHM13v2.0" in refs:
            return "T2T-CHM13v2.0"
        return list(refs.keys())[0]
    return None


# The single reference used for orientation
ORIENT_REF = get_orient_ref()
