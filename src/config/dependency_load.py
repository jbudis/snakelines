import oyaml
import sys

from src.config import load_imports
from src.config import values_retrieval


def load_yaml(file_name):
    """
    Loads yaml file
    :param dependency_file:
    :return: dict - YAML dict
    """
    with open(file_name) as f:
        dependency_dict = oyaml.load(f)

    return dependency_dict


def load_config(config):
    """
    Loads the YAML config "config.yaml" if not provided as --configfile.
    :param config: OrderedDict - config dictionary so far
    :return: OrderedDict - loaded config dictionary
    """
    # configfile: srcdir('config.yaml') - do not use it like this, since configfile uses dict instead of OrderedDict
    # try if we have good config:
    try:
        config['samples']
    except (NameError, KeyError):
        print('WARNING: "--configfile" not specified, using "config.yaml".')
        # load default config.yaml
        try:
            cfg = load_yaml('config.yaml')
            cfg.update(config)
            config = cfg
        except FileNotFoundError:
            print('ERROR: "config.yaml" not found in the current dir.')

    # return loaded config
    return config


def add_temp(text, add_temp=True):
    """
    Add temporary 'temp(<text>)' flag to text.
    :param text: str - text
    :param add_temp: bool - if we should add the temporary flag
    :return: str - text with temporary flag
    """
    if add_temp:
        return 'temp(' + text + ')'
    return text


def check_dependency(dependency_file, config_dict):
    """
    Checks the dependencies of all things from config file.
    :param dependency_file: str - filename for dependency file
    :param config_dict: dict - dict for config file
    :return: str - empty if no warning is raised
    """

    # empty inputs/outputs:
    inputs = {}
    outputs = {}

    # load dependency file
    dependency_dict = load_yaml(dependency_file)

    # for every import file check the dependency and outputs
    import_files = load_imports.import_files_in_config(config_dict, "")
    for import_file, _ in import_files:
        # strip the '/method.snake' and check if exists in dependency file
        path = '.'.join(import_file.split('.')[:-1]).strip('/')
        path = '/'.join(path.split('/')[:-1])
        subd = values_retrieval.extract_from_config(dependency_dict, path.split('/'))

        assert subd is not None, "ERROR: missing values in dependency config file: {}".format(path)

        # resolve dependency
        if 'depends' in subd:
            for dependency in subd['depends']:
                subc = values_retrieval.extract_from_config(config_dict, dependency.split('/'))

                if subc is None:
                    print("WARNING: dependency of {} not met: {}".format(path, dependency), file=sys.stderr)
                else:
                    print("Dependency of {} met: {}".format(path, dependency), file=sys.stderr)

        # resolve output
        assert 'output' in subd, "ERROR: missing 'output' value in {} (define it at least with emtpy list as 'output: []')".format(path)

        if 'output' in subd:
            for output_file in subd['output']:
                path_underscored = path.replace('/', '_') + '_'
                subc = values_retrieval.extract_from_config(config_dict, path.split('/'))
                temporary = 'temporary' in subc and subc['temporary']
                if 'from' in subd['output'][output_file] and 'to' in subd['output'][output_file]:
                    inputs[path_underscored + output_file] = add_temp(subd['output'][output_file]['from'], temporary)
                    outputs[path_underscored + output_file] = add_temp(subd['output'][output_file]['to'], temporary)
                elif type(subd['output'][output_file]) is str:
                    inputs[path_underscored + output_file] = add_temp(subd['output'][output_file], temporary)
                else:
                    print("WARNING: weird things in dependency file {}".format(subd['output'][output_file]), file=sys.stderr)

    # return required outputs and inputs
    return inputs, outputs
