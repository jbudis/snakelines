def extract_from_config(config, items, default_value = None):
    """
    Load items from the configuration file based on series of keys. If some key is not found in hierarchy,
    return default value specified as the parameter.

    :param config: dict: configuration dictionary
    :param items: list str: series of keys to hierarchically find item in the configuration file
    :param default_value: value to return if the series of keys are not present in the configuration
    :return: dict/value - subconfig of the configuration
    """
    subconfig = config
    for item in items:
        if item not in subconfig:
            return default_value
        subconfig = subconfig[item]
    return subconfig
