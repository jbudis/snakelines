import oyaml


def impute_mapping_types(pipeline, workflow, config):
    """
    Impute map types into the mapping post-process steps
    :param pipeline: Pipeline - pipeline to modify
    :param workflow: ?
    :param config: dict - configuration dictionary
    :return: None
    """
    ORIGINAL_MAPPING = 'original'
    setattr(pipeline, 'processed_map_type', ORIGINAL_MAPPING)

    for config_file in workflow.configfiles:

        # Load configuration with dictionary that keeps order of the preprocess steps
        file_config = oyaml.load(open(config_file))

        # Find configuration file with reads->preprocess attributes
        if 'postprocess' not in file_config.get('mapping', {}):
            continue

        # Add new read_type attribute to the main config, so rules for preprocess steps may deduce input read files
        previous_type = ORIGINAL_MAPPING
        for map_type in file_config['mapping']['postprocess'].keys():
            config['mapping']['postprocess'][map_type]['input_map_type'] = previous_type
            previous_type = map_type

        setattr(pipeline, 'postprocessed_map_type', previous_type)


def impute_read_types(pipeline, workflow, config):
    """
    Impute read types into the mapping post-process steps
    :param pipeline: Pipeline - pipeline to modify
    :param workflow: ?
    :param config: dict - configuration dictionary
    :return: None
    """
    # Impute read types into the read preprocess steps
    ORIGINAL_READS = 'original'

    # TODO preprocessed should be called processed to unite with mapping
    setattr(pipeline, 'preprocessed_read_type', ORIGINAL_READS)

    for config_file in workflow.configfiles:

        # TODO avoid using oyaml, return to list yaml representation
        # Load configuration with dictionary that keeps order of the preprocess steps
        file_config = oyaml.load(open(config_file))

        # Find configuration file with reads->preprocess attributes
        if 'preprocess' not in file_config.get('reads', {}):
            continue

        # Add new read_type attribute to the main config, so rules for preprocess steps may deduce input read files
        previous_type = ORIGINAL_READS
        for read_type in file_config['reads']['preprocess'].keys():
            config['reads']['preprocess'][read_type]['input_read_type'] = previous_type
            previous_type = read_type

        setattr(pipeline, 'preprocessed_read_type', previous_type)