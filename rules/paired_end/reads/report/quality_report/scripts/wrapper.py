import os
import tempfile

from snakemake.shell import shell


with tempfile.TemporaryDirectory() as tmpdir:
    shell(
        "fastqc "
        " --outdir {tmpdir}"
        " --extract"
        " --threads {snakemake.threads}"
        " --dir {tmpdir}"
        " {snakemake.input.reads}"
        " >  {snakemake.log.out}"
        " 2> {snakemake.log.err}"
    )

    PRESUMED_SUFFIX = ".fastq.gz"
    if not snakemake.input.reads.endswith(PRESUMED_SUFFIX):
        raise ValueError(f"{snakemake.input.reads} does not ends with {PRESUMED_SUFFIX}")

    base_name = os.path.basename(snakemake.input.reads).replace(".fastq.gz", "")
    print("basename",base_name)
    html_path = os.path.join(tmpdir, f"{base_name}_fastqc.html")
    zip_path = os.path.join(tmpdir, f"{base_name}_fastqc.zip")

    fastqc_datapath = os.path.join(tmpdir, f"{base_name}_fastqc", "fastqc_data.txt")
    summary_path = os.path.join(tmpdir, f"{base_name}_fastqc", "summary.txt")

    shell("ls {tmpdir}")
    shell("mv {html_path} {snakemake.output.html}")
    shell("mv {zip_path} {snakemake.output.zip}")
    shell("mv {fastqc_datapath} {snakemake.output.data}")
    shell("mv {summary_path} {snakemake.output.txt}")

