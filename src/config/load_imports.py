from collections import OrderedDict


def import_files_in_config(subconfig, item_path):
    """
    Recursively traverse configuration file and identify possible snake imports annotated by "method:" node

    :param subconfig: part of configuration file to traverse
    :param item_path: already traversed path in the configuration file
    """

    # Directories are searched recursively
    if type(subconfig) == dict or type(subconfig) == OrderedDict:

        # Search all elements of the current item
        for next_item_name in subconfig:
            if next_item_name == 'method':
                yield '{}/{}.snake'.format(item_path, subconfig[next_item_name]), subconfig
                continue

            next_item_path = '{}/{}'.format(item_path, next_item_name)

            yield from import_files_in_config(subconfig[next_item_name], next_item_path)

    # Items at the end of the configuration tree are checked as well
    # if type(subconfig) == str:
    #    yield from import_files_in_config(None, '{}/{}'.format(item_path, subconfig))
