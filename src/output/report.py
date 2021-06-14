import json
import os
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


def multiqc_report(out_dir: str, out_list_filename: str):
    basename = 'multiqc_report'
    args = [
        'multiqc',
        '--file-list',
        out_list_filename,
        '--outdir',
        out_dir,
        '--filename',
        '%s.html' % basename,
        '--data-format',
        'tsv',
        '--force',
        '--verbose'
    ]
    subprocess.run(args, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def dig_value(data: dict, keys: list):
    value = data
    for key in keys:
        if isinstance(value, dict):
            if key in value:
                value = value[key]
        elif isinstance(value, list):
            # merge listed dicts
            merged = {}
            for list_value in value:
                if isinstance(list_value, dict) and key in list_value:
                    merged.update(list_value[key])
            value = merged
        else:
            break

    if isinstance(value, dict):
        # the value was not found
        value = None

    return value


def validate_multiqc_summary_config(validator: list, value, default='fail'):
    result = default
    passed = False
    for validation_class in validator:
        class_label = next(iter(validation_class))
        range_list = validation_class[class_label]
        for from_to in range_list:
            from_value = from_to.get('from')
            to_value = from_to.get('to')

            if from_value is not None and to_value is not None:
                assert float(from_value) <= float(to_value)

            from_passed = True
            if from_value is not None:
                from_passed = float(value) >= float(from_value)

            to_passed = True
            if to_value is not None:
                to_passed = float(value) <= float(to_value)

            if from_passed and to_passed:
                passed = True
                break

        if passed:
            result = class_label
            break

    return result


def format_value(value, format_str: str, humanize: str) -> str:
    formatted_value = float(value)
    if humanize == 'G':
        formatted_value = formatted_value / 1e9
    if humanize == 'M':
        formatted_value = formatted_value / 1e6
    elif humanize == 'K':
        formatted_value = formatted_value / 1e3

    if format_str:
        formatted_value = format_str.format(formatted_value)

    return str(formatted_value)


def multiqc_summary_report(multiqc_dir: str, sample_dir: str, samples: list, metrics: list):
    data_filepath = os.path.join(multiqc_dir, 'multiqc_report_data', 'multiqc_data.json')
    with open(data_filepath, 'rt') as data_file:
        data = json.load(data_file)

        for sample in samples:
            summary_filepath = os.path.join(sample_dir, sample, 'summary.tsv')

            with open(summary_filepath, 'wt') as sample_file:

                for metric in metrics:
                    assert 'key' in metric
                    key = metric['key'].split('.')
                    name = metric.get('name', key[-1])
                    validator = metric.get('validator')
                    format_str = metric.get('format')
                    unit = metric.get('humanize')

                    sample_key = [sample if segment == '{sample}' else segment for segment in key]
                    value = dig_value(data, sample_key)

                    if value is not None:
                        sample_file.write(f'{name}-value\t{value}\n')
                        sample_file.write(f'{name}-text\t{format_value(value, format_str, unit)}\n')

                        if validator is not None:
                            result = validate_multiqc_summary_config(validator, value)
                            sample_file.write(f'{name}-validation\t{result}\n')
                        else:
                            sample_file.write(f'{name}-validation\tNOT_APPLICABLE\n')

                    else:  # missing value
                        sample_file.write(f'{name}-value\t\n')
                        sample_file.write(f'{name}-text\t\n')
                        sample_file.write(f'{name}-validation\tNOT_APPLICABLE\n')
