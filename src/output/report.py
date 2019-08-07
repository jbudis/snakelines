import subprocess
import sys


def snakemake_report(outdir: str, snakefile: str, configfile: str):
    report_file = outdir + '/_execution/report.html'
    args = [
        'snakemake',
        '--snakefile',
        snakefile,
        '--configfile',
        configfile,
        '--report',
        report_file
    ]
    subprocess.run(args)
