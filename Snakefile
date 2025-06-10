###############################################################################
# CONFIG
###############################################################################
import os, glob, re, sys
from snakemake.io import glob_wildcards

configfile: "config.yaml"

###############################################################################
# FUNCTIONS
###############################################################################

def derive_samples(reads_dir, suffix, remove_index):
    """
    Discover samples in reads_dir that match *suffix* (which contains one '*').
    If remove_index is True, strip '_S<digits>' plus the suffix.
    Otherwise strip only the literal tail after the '*'.
    """
    # fixed part after the *
    paths = glob.glob(os.path.join(reads_dir, f"*{suffix}"))
    names = [os.path.basename(p) for p in paths]
    if remove_index:
        # remove '_S<digits>' plus the fixed suffix
        pattern = r"_S\d+" + re.escape(suffix) + r"$"
        samples = [re.sub(pattern, "", n) for n in names]
    else:
        # keep the _S<digits> part, drop only the fixed suffix
        pattern = re.escape(suffix) + r"$"
        samples = [re.sub(pattern, "", n) for n in names]
    return samples

def _get_cores_from_argv():
    for flag in ("--cores", "-j"):
        if flag in sys.argv:
            i = sys.argv.index(flag)
            if i+1 < len(sys.argv):
                try:
                    return int(sys.argv[i+1])
                except ValueError:
                    pass
    return os.cpu_count()

def all_inputs(wildcards):
    # get the mapping file from the checkpoint
    map_file = checkpoints.parse_profile.get().output[0]
    # read it
    with open(map_file) as fh:
        samples = [line.split("\t",1)[0] for line in fh]
    # build per-sample consensus + nextclade targets
    consensus = [
        os.path.join(OUTDIR, s, f"{s}.consensus.fasta")
        for s in samples
    ]
    consensus_ambiguities = [
        os.path.join(OUTDIR, s, f"{s}.IUPAC_consensus.fasta")
        for s in samples
    ]
    reports   = [
        os.path.join(OUTDIR, s, f"{s}.nextclade_report.tsv")
        for s in samples
    ]
    # and your final QC outputs
    others = [
        os.path.join(OUTDIR, "overall_samples_qc_summary.csv"),
        os.path.join(OUTDIR, "all_nextclade_merged.tsv")
    ]
    # plus the global nextclade‐ready flags
    flags = [
        os.path.join(OUTDIR, f".nextclade_{st}_ready.txt")
        for st in REFS.keys()
    ]
    return [map_file] + flags + consensus + consensus_ambiguities + reports + others


def get_mapped_samples():
    import pandas as pd
    df = pd.read_csv(MAP, sep="\t", header=None, names=["sample","subtype"])
    return df["sample"].tolist()

###############################################################################
# PARAMS
###############################################################################

THREADS           = _get_cores_from_argv()
# READS_DIR         = config["reads_dir"]
OUTDIR            = config["outdir"]
PRIMERS           = config["primers_bed"]
MIN_DEPTH         = config["min_depth"]
REFS              = config["refs"]      # dict: subtype -> fasta
SYLDB             = config["syldb"]
SNP_THRESH        = config["snp_threshold"]
CONSENSUS_THRESH  = config["consensus_threshold"]
SUFFIX1           = config["suffix1"]
SUFFIX2           = config["suffix2"]

# file produced by the parser
MAP = os.path.join(OUTDIR, "selected_sample_subtypes.tsv")

###############################################################################
# DISCOVER SAMPLES
###############################################################################
# SAMPLES = derive_samples(READS_DIR, SUFFIX1, config["remove_sample_index"])

REMOVE_INDEX = config["remove_sample_index"]       # True / False

###############################################################################
# RULE: all
###############################################################################
rule all:
    input:
        all_inputs
###############################################################################
# RULE: sylph_profile (one global job)
###############################################################################
rule sylph_profile:
    params:
        suffix1 = config["suffix1"],
        suffix2 = config["suffix2"],
        READS_DIR = config['reads_dir']
    output:
        custom=os.path.join(OUTDIR, "sylph_subtype_profiling.tsv")
    threads: THREADS
    shell:
        """
        mkdir -p {OUTDIR}
        sylph profile {SYLDB} -c 100 \
            -1 {params.READS_DIR}/*{params.suffix1} \
            -2 {params.READS_DIR}/*{params.suffix2} \
            -t {threads} \
        > {output.custom} 2> /dev/null
        """

###############################################################################
# RULE: parse_profile -> selected_sample_subtypes
###############################################################################
checkpoint parse_profile:
    input:
        profile = rules.sylph_profile.output.custom
    output:
        MAP
    run:
        import os, re, pandas as pd

        remove_index = config["remove_sample_index"]
        tail = SUFFIX1  # literal suffix

        # choose correct regex
        patt = re.compile(
            (r"_S\d+" if remove_index else "") + re.escape(tail) + r"$"
        )

        df = pd.read_csv(input.profile, sep="\t", comment="#")
        df["sample"] = df["Sample_file"].apply(
            lambda p: patt.sub("", os.path.basename(p))
        )
        idx = df.groupby("sample")["Taxonomic_abundance"].idxmax()
        mapping = df.loc[idx, ["sample", "Contig_name"]]
        os.makedirs(os.path.dirname(output[0]), exist_ok=True)
        mapping.to_csv(output[0], sep="\t", header=False, index=False)

###############################################################################
# RULE: choose_ref (per-sample)  -> ref.txt  (holds real fasta path)
###############################################################################
rule choose_ref:
    input:  MAP
    output: ref_txt=os.path.join(OUTDIR, "{sample}", "ref.txt")
    run:
        os.makedirs(os.path.dirname(output.ref_txt), exist_ok=True)
        subtype = None
        with open(input[0]) as fh:
            for line in fh:
                s, st = line.strip().split("\t")
                if s == wildcards.sample:
                    subtype = st
                    break
        if subtype is None:
            raise ValueError(f"{wildcards.sample} not found in {input[0]}")
        ref_path = REFS[subtype]
        with open(output.ref_txt, "w") as fh:
            fh.write(ref_path + "\n")

###############################################################################
# RULES: per-sample pipeline (all take ref.txt as *input*, read it in shell)
###############################################################################
rule fastp:
    input:
        MAP,
        r1 = lambda wc: glob.glob(os.path.join(
                 config["reads_dir"],
                 f"{wc.sample}_S*{config["suffix1"]}" if REMOVE_INDEX
                 else f"{wc.sample}{config["suffix1"]}"
             ))[0],
        r2 = lambda wc: glob.glob(os.path.join(
                 config["reads_dir"],
                 f"{wc.sample}_S*{config["suffix2"]}" if REMOVE_INDEX
                 else f"{wc.sample}{config["suffix2"]}"
             ))[0]
    output:
        r1   = temp(os.path.join(OUTDIR, "{sample}", "{sample}_R1.trimmed.fastq.gz")) \
                   if not config["keep_trimmed_reads"] \
                   else os.path.join(OUTDIR, "{sample}", "{sample}_R1.trimmed.fastq.gz"),
        r2   = temp(os.path.join(OUTDIR, "{sample}", "{sample}_R2.trimmed.fastq.gz")) \
                   if not config["keep_trimmed_reads"] \
                   else os.path.join(OUTDIR, "{sample}", "{sample}_R2.trimmed.fastq.gz"),
        html=os.path.join(OUTDIR, "{sample}", "{sample}_fastp.html"),
        json=os.path.join(OUTDIR, "{sample}", "{sample}_fastp.json")
    threads: THREADS
    shell:
        """
        mkdir -p {OUTDIR}/{wildcards.sample}
        fastp -i {input.r1} -I {input.r2} \
              -o {output.r1} -O {output.r2} \
              -w {threads} \
              --detect_adapter_for_pe \
              --cut_front \
              --cut_tail \
              --cut_mean_quality 20 \
              --correction \
              --length_required 50 \
              -h {output.html} \
              -j {output.json} 2> /dev/null
        """

rule index_refs:
    input:
        ready = expand(os.path.join(OUTDIR, ".nextclade_{subtype}_ready.txt"),
                subtype=REFS.keys()),
    params:
        refs=list(REFS.values()),
    output:
        idx=[
            f"{ref}{ext}"
            for ref in REFS.values()
            for ext in (".bwt", ".pac", ".ann", ".amb", ".sa", ".fai")
        ]
    threads: 1
    shell:
        """
        for REF in {params.refs}; do
            # index only if missing – avoids pointless re-runs
            [ -f "$REF.bwt" ] || bwa index "$REF" 2> /dev/null
            [ -f "$REF.fai" ] || samtools faidx "$REF" 2> /dev/null
        done
        """

rule bwa_mem:
    input:
        idx_rule = rules.index_refs.output,          # ensure global index exists
        ref_txt  = rules.choose_ref.output.ref_txt,  # per-sample path file
        r1       = rules.fastp.output.r1,
        r2       = rules.fastp.output.r2
    output:
        bam = temp(os.path.join(OUTDIR, "{sample}", "{sample}.sorted.bam")),
        csi = temp(os.path.join(OUTDIR, "{sample}", "{sample}.sorted.bam.csi"))
    threads: 4
    shell:
        """
        REF=$(cat {input.ref_txt})
        bwa mem -t {threads} "$REF" {input.r1} {input.r2} 2> /dev/null \
        | samtools sort -@ {threads} | \
        samtools view -@ {threads} -F 4 --write-index -o {output.bam} -
        """

rule choose_bed:
    input:  MAP
    output: bed_txt=os.path.join(OUTDIR, "{sample}", "bed.txt")
    run:
        import os
        os.makedirs(os.path.dirname(output.bed_txt), exist_ok=True)
        subtype = None
        with open(input[0]) as fh:
            for line in fh:
                s, st = line.rstrip().split("\t")
                if s == wildcards.sample:
                    subtype = st
                    break
        if subtype is None:
            raise ValueError(f"{wildcards.sample} not found in {input[0]}")
        bed_path = config["primers_bed"][subtype]
        with open(output.bed_txt, "w") as f:
            f.write(bed_path + "\n")

rule clip:
    input:
        bed_txt = rules.choose_bed.output.bed_txt,
        bam     = rules.bwa_mem.output.bam
    output:
        bam = temp(os.path.join(OUTDIR, "{sample}", "{sample}.clipped.bam")),
        bai = temp(os.path.join(OUTDIR, "{sample}", "{sample}.clipped.bam.bai")),
        log = os.path.join(OUTDIR, "{sample}", "{sample}.ampliconclip.log")
    threads: THREADS
    shell:
        """
        BED=$(cat {input.bed_txt})
        samtools ampliconclip -b $BED -@ {threads} \
        --strand --both-ends  -o - {input.bam} 2> {output.log} \
        | samtools sort -@ {threads} -o {output.bam} - 2> /dev/null
        samtools index {output.bam}
        """

rule indelqual:
    input:
        ref_txt=rules.choose_ref.output.ref_txt,
        bam    =rules.clip.output.bam
    output:
        bam=os.path.join(OUTDIR, "{sample}", "{sample}.indelq.bam"),
        bai=os.path.join(OUTDIR, "{sample}", "{sample}.indelq.bam.bai")
    threads: THREADS
    shell:
        """
        REF=$(cat {input.ref_txt})
        lofreq indelqual --dindel -f $REF -o {output.bam} {input.bam} 2> /dev/null
        samtools index {output.bam}
        """

rule call_variants:
    input:
        ref_txt=rules.choose_ref.output.ref_txt,
        bam    =rules.indelqual.output.bam
    output:
        vcf=os.path.join(OUTDIR, "{sample}", "{sample}.raw.vcf.gz"),
        tbi=os.path.join(OUTDIR, "{sample}", "{sample}.raw.vcf.gz.tbi")
    threads: 4
    shell:
        """
        REF=$(cat {input.ref_txt})
        lofreq call-parallel --no-baq --call-indels --pp-threads {threads} \
            -f $REF -o {wildcards.sample}.raw.vcf {input.bam} 2> /dev/null
        bgzip -f {wildcards.sample}.raw.vcf
        tabix -p vcf {wildcards.sample}.raw.vcf.gz
        mv {wildcards.sample}.raw.vcf.gz     {output.vcf}
        mv {wildcards.sample}.raw.vcf.gz.tbi {output.tbi}
        """

rule filter_variants:
    input: vcf=rules.call_variants.output.vcf
    params: 
        snp_thresh=SNP_THRESH,
        min_depth=MIN_DEPTH,
        min_qual=20,
        indel_thresh=config['indel_threshold'],
    output:
        filt=os.path.join(OUTDIR, "{sample}", "{sample}.filt.vcf.gz"),
        filt_tbi=os.path.join(OUTDIR, "{sample}", "{sample}.filt.vcf.gz.tbi")
    threads: THREADS
    shell:
        """
        bcftools +fill-tags {input.vcf} -Ou -- -t "TYPE" | \
        bcftools norm -Ou -a -m -  2> /dev/null | \
        bcftools view -f 'PASS' \
            -i '(
                  (TYPE="indel" && INFO/AF >= {params.indel_thresh}  && INFO/DP >= {params.min_depth} && QUAL >= {params.min_qual})
                  ||
                  (TYPE!="indel"   && INFO/AF >= {params.snp_thresh}    && INFO/DP >= {params.min_depth} && QUAL >= {params.min_qual})
                )' \
            -Oz -o {output.filt}
        tabix -f -p vcf {output.filt}
        """

rule add_fake_gt:
    input: filt=rules.filter_variants.output.filt
    output:
        gt=os.path.join(OUTDIR, "{sample}", "{sample}.gt.vcf.gz"),
        gt_tbi=os.path.join(OUTDIR, "{sample}", "{sample}.gt.vcf.gz.tbi")
    shell:
        """
        bash scripts/add_fake_genotype.sh \
          -i {input.filt} -g 1/1 -o {output.gt}
        tabix -f -p vcf {output.gt}
        """

rule set_vcf_genotype:
    input:
        vcf_file=rules.add_fake_gt.output.gt,
    output:
        temp_vcf=temp(os.path.join(OUTDIR, "{sample}", "{sample}.temp.vcf.gz")),
        vcf_file=os.path.join(OUTDIR, "{sample}", "{sample}.final.vcf.gz"),
    log:
        os.path.join(OUTDIR, "{sample}", "{sample}.bcftools_setGT.log"),
    params:
        snv_freq=config["snp_threshold"],
        con_freq=config["consensus_threshold"],
        indel_freq=config["indel_threshold"],
    message:
        "setting conditional GT for {wildcards.sample}"
    shell:
        """
        cp {input.vcf_file} {output.temp_vcf}
        bcftools index -f {output.temp_vcf}
        bcftools +setGT {output.temp_vcf} -- -t q -i 'GT="1/1" && INFO/AF < {params.con_freq}' -n 'c:0/1' 2>> {log} | \
        bcftools +setGT -- -t q -i 'TYPE="indel" && INFO/AF < {params.indel_freq}' -n . 2>> {log} | \
        bcftools +setGT -o {output.vcf_file} -- -t q -i 'GT="1/1" && INFO/AF >= {params.con_freq}' -n 'c:1/1' 2>> {log}
        bcftools index -f {output.vcf_file}
        """

rule variants_bed:
    input: vcf=rules.set_vcf_genotype.output.vcf_file
    params: min_depth=MIN_DEPTH
    output: variants_bed=temp(os.path.join(OUTDIR, "{sample}", "{sample}.variants.bed"))
    shell:
        """
        bcftools query -f'%CHROM\t%POS0\t%END\n' {input.vcf} > {output.variants_bed}
        """

rule depth_mask:
    input: 
        bam=rules.indelqual.output.bam,
        variants_bed=rules.variants_bed.output.variants_bed,
    params: 
        min_depth=MIN_DEPTH
    output:
        mask=os.path.join(OUTDIR, "{sample}", "{sample}.lowdepth.bed")
    shell:
        """
        bedtools genomecov -bga -ibam {input.bam} \
        | awk '$4 < {params.min_depth}' \
        | bedtools subtract -a - -b {input.variants_bed} > {output.mask}
        """

rule consensus:
    input:
        ref_txt = rules.choose_ref.output.ref_txt,
        mask    = rules.depth_mask.output.mask,
        vcf     = rules.set_vcf_genotype.output.vcf_file
    params:
        prefix="{sample}",
        consensus_threshold = CONSENSUS_THRESH,
        min_depth=MIN_DEPTH,
    output:
        fasta=os.path.join(OUTDIR, "{sample}", "{sample}.consensus.fasta")
    log:
        os.path.join(OUTDIR, "{sample}", "{sample}.bcftools_consensus.log")
    shell:
        """
        REF=$(cat {input.ref_txt})
        bcftools consensus -p "{params.prefix} " \
            -f $REF \
            --mark-del '-' \
            -m {input.mask} \
            -i 'INFO/DP >= {params.min_depth} & INFO/AF >= {params.consensus_threshold} & GT!="mis"' {input.vcf} -o {output.fasta} 2 > {log}
        """

rule consensus_ambiguities:
    input:
        ref_txt = rules.choose_ref.output.ref_txt,
        mask    = rules.depth_mask.output.mask,
        vcf     = rules.set_vcf_genotype.output.vcf_file
    params:
        prefix="{sample}",
        consensus_threshold = CONSENSUS_THRESH,
        min_depth=MIN_DEPTH,
    output:
        fasta=os.path.join(OUTDIR, "{sample}", "{sample}.IUPAC_consensus.fasta")
    log:
        os.path.join(OUTDIR, "{sample}", "{sample}.bcftools_iupac_consensus.log")
    shell:
        """
        REF=$(cat {input.ref_txt})
        bcftools consensus -p "{params.prefix} IUPAC " \
            -f $REF \
            --mark-del '-' \
            -m {input.mask} \
            -H I \
            -i 'INFO/DP >= {params.min_depth} & GT!="mis"' {input.vcf} -o {output.fasta} 2 > {log}
        """


rule update_nextclade:
    wildcard_constraints:
        subtype="RSV[AB]"
    output:
        ready=os.path.join(OUTDIR, ".nextclade_{subtype}_ready.txt")
    params:
        dir   = lambda wc: config["nextclade_dataset"][wc.subtype],
        dname = lambda wc: config["nextclade_name"][wc.subtype]
    threads: 1
    shell:
        """
        mkdir -p {params.dir}
        nextclade dataset get --name {params.dname} -o {params.dir} 2> /dev/null
        echo "dataset ready" > {output.ready}
        """

rule choose_dataset:
    input:  MAP
    output: ds_txt=os.path.join(OUTDIR, "{sample}", "dataset.txt")
    run:
        import os
        # look up subtype for this sample
        sub = None
        for line in open(input[0]):
            s, st = line.rstrip().split("\t")
            if s == wildcards.sample:
                sub = st
                break
        if sub is None:
            raise ValueError(f"{wildcards.sample} not found in {input[0]}")
        ds_path = config["nextclade_dataset"][sub]
        os.makedirs(os.path.dirname(output.ds_txt), exist_ok=True)
        with open(output.ds_txt, "w") as fh:
            fh.write(ds_path + "\n")

rule nextclade:
    input:
        dataset_txt = rules.choose_dataset.output.ds_txt,
        fasta       = os.path.join(OUTDIR, "{sample}", "{sample}.consensus.fasta"),
        # depend on the global ready flag for each subtype
        ready_rsva  = os.path.join(OUTDIR, ".nextclade_RSVA_ready.txt"),
        ready_rsvb  = os.path.join(OUTDIR, ".nextclade_RSVB_ready.txt")
    output:
        report = os.path.join(OUTDIR, "{sample}", "{sample}.nextclade_report.tsv")
    threads: 1
    shell:
        """
        DATASET=$(cat {input.dataset_txt})
        nextclade run \
          --input-dataset=$DATASET \
          --output-tsv={output.report} \
          --jobs 1 {input.fasta} 2> /dev/null
        """

###############################################################################
# RULE: qc_samples → results/overall_samples_qc_summary.csv
###############################################################################
rule qc_samples:
    input:
        profiling = rules.sylph_profile.output.custom,
        reports   = lambda wc: [
            os.path.join(OUTDIR, s, f"{s}.nextclade_report.tsv")
            for s in get_mapped_samples()
        ]
    output:
        csv        = os.path.join(OUTDIR, "overall_samples_qc_summary.csv"),
        merged_nc  = os.path.join(OUTDIR, "all_nextclade_merged.tsv")
    run:
        import os, re, pandas as pd

        # ------------------------------------------------------------------ #
        remove_index = config["remove_sample_index"]
        tail_r1      = SUFFIX1
        mixed_thresh = float(config["mixed_threshold"])
        # ------------------------------------------------------------------ #

        # ---------------- 1 - SYLPH summary --------------------------------
        sylph = pd.read_csv(input.profiling, sep="\t", comment="#")

        patt = re.compile(
            (r"_S\d+" if remove_index else "") + re.escape(tail_r1) + r"$"
        )
        sylph["sample"] = sylph["Sample_file"].apply(
            lambda p: patt.sub("", os.path.basename(p))
        )

        idx = sylph.groupby("sample")["Taxonomic_abundance"].idxmax()
        sylph = sylph.loc[idx, [
            "sample", "Contig_name",
            "Taxonomic_abundance", "Mean_cov_geq1", "Median_cov"
        ]].rename(columns={
            "Contig_name":           "subtype",
            "Taxonomic_abundance":   "subtype_relative_abundance",
            "Mean_cov_geq1":         "mean_depth",
            "Median_cov":            "median_depth"
        })

        # ---------------- 2 - Nextclade QC ----------------------------------
        qc_rows   = []
        nc_frames = []

        for rpt in input.reports:
            df = pd.read_csv(rpt, sep="\t")
            if df.empty:
                continue
            nc_frames.append(df)                        # keep full report

            row       = df.iloc[0]
            sample_id = str(row["seqName"]).split(" ", 1)[0]

            qc_rows.append({
                "sample":                       sample_id,
                "clade":                        row["clade"],
                "coverage":                     round(float(row["coverage"])*100, 2),
                "qc.overallStatus":             row["qc.overallStatus"],
                "totalSubstitutions":           row["totalSubstitutions"],
                "totalDeletions":               row["totalDeletions"],
                "totalInsertions":              row["totalInsertions"],
                "totalMissing":                 row["totalMissing"],
                "totalAminoacidSubstitutions":  row["totalAminoacidSubstitutions"],
                "totalAminoacidDeletions":      row["totalAminoacidDeletions"],
                "totalAminoacidInsertions":     row["totalAminoacidInsertions"]
            })

        nextclade_qc = pd.DataFrame(qc_rows)

        # ---------------- 3 - Merge & mixed flag ---------------------------
        merged = pd.merge(sylph, nextclade_qc, on="sample", how="inner")
        cutoff = (1.0 - mixed_thresh) * 100.0
        merged["mixed_sample"] = merged["subtype_relative_abundance"] < cutoff

        merged = merged[[
            "sample", "subtype", "clade",
            "subtype_relative_abundance", "mixed_sample", "coverage",
            "mean_depth", "median_depth",
            "qc.overallStatus",
            "totalSubstitutions", "totalDeletions", "totalInsertions",
            "totalMissing",
            "totalAminoacidSubstitutions",
            "totalAminoacidDeletions",
            "totalAminoacidInsertions"
        ]]

        # write outputs
        os.makedirs(os.path.dirname(output.csv), exist_ok=True)
        merged.to_csv(output.csv, index=False)

        if nc_frames:
            pd.concat(nc_frames, ignore_index=True) \
              .drop(columns=["index"], errors="ignore") \
              .to_csv(output.merged_nc, sep="\t", index=False)
        else:
            # ensure file exists even if no reports
            pd.DataFrame().to_csv(output.merged_nc, sep="\t", index=False)
