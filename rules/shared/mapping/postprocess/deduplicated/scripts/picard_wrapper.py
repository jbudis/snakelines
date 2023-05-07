import tempfile
from snakemake.shell import shell


sys.stderr = open(snakemake.log.err, "w")
sys.stdout = open(snakemake.log.out, "w")


with tempfile.TemporaryDirectory() as tmpdir:
    shell(
        "picard MarkDuplicates"
        " I={snakemake.input.bam}"
        " O={snakemake.output.bam}"
        " M={snakemake.log.stat}"
        " TMP_DIR={tmpdir}"
        " VALIDATION_STRINGENCY=SILENT"
        " 1>> {snakemake.log.out} 2>> {snakemake.log.err}"
    )
