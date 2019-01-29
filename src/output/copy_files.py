import shutil
import errno
import os


# TODO this should be renamed to matching or deleted - unused
def copy_input_files_with_consistent_output_names(input_files, output_files):
    """
    Compare names of the inputs and outputs of the rule. All items with the same name in inputs and outputs are copied,
    therefore generate output files.

    :param input_files: input files of the rule
    :param output_files: output files of the rule
    """

    # Iterate through all defined outputs in the rule
    for item_name in output_files.__dict__:

        # Skip internal snakemake attributes
        if item_name.startswith('_'):
            continue

        # Skip output items do not have matching partner in input
        if not hasattr(input_files, item_name):
            continue

        # Check if items with the same name have also same data type, both should be list or str
        input_item, output_item = getattr(input_files, item_name), getattr(output_files, item_name)
        assert type(input_item) == type(output_item), 'Output and input items with name {} for the rule have different data types {} and {}'.format(
            item_name, type(input_item), type(output_item))

        # Check if items with the same name have also same length
        assert len(input_item) == len(output_item), 'Output and input items with name {} \
                for the rule have different lengths {} and {}'.format(item_name, len(input_item), len(output_item))

        # Single file items are directly copied
        if type(input_item) == str:
            shutil.copy(input_item, output_item)

        # List items are copied item by item
        for in_item, out_item in zip(input_item, output_item):
            shutil.copy(in_item, out_item)


def copy_config(report_dir, workflow):
    """
    Copy configuration file (snakemake --configfile argument) to the report directory
    :param report_dir - output directory for reported files
    :param workflow: str - global snakemake variable
    :return: None
    """
    config_file = workflow.overwrite_configfile
    if config_file:
        report_file = '{}/{}'.format(report_dir, os.path.basename(config_file))
        copy_with_makedirs(config_file, report_file)


def copy(src, dest):
    """
    Copy a file or directory
    :param src: str - source file or directory
    :param dest: str - destination file or directory
    :return: None
    :raise IOError, PermissionError: if cannot copy, e.g. the directory structure does not exists
    """
    if os.path.isdir(src):
        if os.path.exists(dest):
            shutil.rmtree(dest)
        shutil.copytree(src, dest)
    else:
        shutil.copy(src, dest)


def copy_with_makedirs(src, dest):
    """
    Copy a file or directory and create directory structure if necessary.
    :param src: str - source file or directory
    :param dest: str - destination file or directory
    :return: None
    :raise IOError: if cannot copy
    """
    try:
        copy(src, dest)
    except IOError as e:
        # ENOENT(2): file does not exist, raised also on missing dest parent dir
        if e.errno != errno.ENOENT:
            raise
        # try creating parent directories
        os.makedirs(os.path.dirname(dest))
        copy(src, dest)


def copy_input_files_with_consisent_names_dict(input_dict, output_dict):
    """
    Compare names of the inputs and outputs of the two dictionaries. All items with the same name in inputs and outputs are copied,
    therefore generate output files.
    :param input_dict: dict - input file names and locations
    :param output_dict: dict - output file names and locations
    :return: int - number of files copied
    """
    files_copied = 0

    for item_name, output_item in output_dict.items():
        if item_name in input_dict:
            input_item = input_dict[item_name]

            # Check if items are str of lists
            assert type(input_item) in [str, list], 'Input item {} \
                    is not str or list({})'.format(input_item, type(input_item))
            assert type(output_item) in [str, list], 'Output item {} \
                    is not str or list ({})'.format(output_item, type(output_item))

            # Check if items with the same name have also same data type
            assert type(input_item) == type(output_item), 'Output and input items with name {} \
                    for the rule have different data types {} and {}'.format(item_name, type(input_item), type(output_item))

            # Check if items with the same name have also same length
            if type(input_item) is list and type(output_item) is list:
                assert len(input_item) == len(output_item), 'Output and input items with name {} \
                        for the rule have different lengths {} and {}'.format(item_name, len(input_item), len(output_item))

            # Single file items are directly copied
            if type(input_item) == str:
                copy_with_makedirs(input_item, output_item)
                files_copied += 1
            else:
                # List items are copied item by item
                for in_item, out_item in zip(input_item, output_item):
                    copy_with_makedirs(in_item, out_item)
                    files_copied += 1

    return files_copied