import tempfile
from snakemake.shell import shell


sys.stderr = open(snakemake.log.err, "w")
sys.stdout = open(snakemake.log.out, "w")


with tempfile.TemporaryDirectory() as tmpdir:
    shell(
        "gatk SelectVariants"
        " --variant {snakemake.input.vcf}"
        " --reference {snakemake.input.fasta}"
        " --output {snakemake.output.vcf}"
        " --tmp-dir {tmpdir}"
        " 1>> {snakemake.log.out} 2>> {snakemake.log.err}"
    )
