import subprocess
from itertools import chain


def runtime_report(outdir: str, snakefile: str, configfile: str):
    report_file = outdir + '/runtime_report.html'
    args = [
        'snakemake',
        '--snakefile',
        snakefile,
        '--configfile',
        configfile,
        '--report',
        report_file,
        '--nolock'
    ]
    subprocess.run(args)


def write_output_filenames(out_filename: str, workflow):
    out_list = []
    for job in workflow.persistence.dag.jobs:
        # noinspection PyProtectedMember
        out_list.append(job.output._plainstrings())
    
    with open(out_filename, 'w') as out_file:
        out_file.write('\n'.join(chain(*out_list)))


def multiqc_report(outdir: str, out_list_filename: str):
    basename = 'multiqc_report'
    args = [
        'multiqc',
        '--file-list',
        out_list_filename,
        '--outdir',
        outdir,
        '--filename',
        '%s.html' % basename,
        '--data-format',
        'tsv',
        '--force',
        '--verbose'
    ]
    subprocess.run(args, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
