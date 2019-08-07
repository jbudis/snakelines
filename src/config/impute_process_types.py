def impute_read_types(pipeline, config):
    """
    Impute read types into the mapping post-process steps
    :param pipeline: Pipeline - pipeline to modify
    :param workflow: ?
    :param config: dict - configuration dictionary
    :return: None
    """
    # Impute read types into the read preprocess steps
    original_type = 'original'
    
    # TODO preprocessed should be called processed to unite with mapping
    setattr(pipeline, 'preprocessed_read_type', original_type)
    
    # Find configuration file with reads->preprocess attributes
    if 'preprocess' not in config.get('reads', {}):
        return
    
    # Add new read_type attribute to the main config, so rules for preprocess steps may deduce input read files
    previous_type = original_type
    for this_type in config['reads']['preprocess'].keys():
        config['reads']['preprocess'][this_type]['input_read_type'] = previous_type
        previous_type = this_type
    
    setattr(pipeline, 'preprocessed_read_type', previous_type)


def impute_mapping_types(pipeline, config):
    """
    Impute map types into the mapping post-process steps
    :param pipeline: Pipeline - pipeline to modify
    :param workflow: ?
    :param config: dict - configuration dictionary
    :return: None
    """
    original_type = 'original'
    setattr(pipeline, 'postprocessed_map_type', original_type)
    
    # Find configuration file with reads->preprocess attributes
    if 'postprocess' not in config.get('mapping', {}):
        return
    
    # Add new read_type attribute to the main config, so rules for preprocess steps may deduce input read files
    previous_type = original_type
    for this_type in config['mapping']['postprocess'].keys():
        config['mapping']['postprocess'][this_type]['input_map_type'] = previous_type
        previous_type = this_type
    
    setattr(pipeline, 'postprocessed_map_type', previous_type)


def impute_variant_types(pipeline, config):
    """
    Impute variant types into the variant post-process steps
    :param pipeline: Pipeline - pipeline to modify
    :param workflow: ?
    :param config: dict - configuration dictionary
    :return: None
    """
    original_type = 'original'
    setattr(pipeline, 'postprocessed_variant_type', original_type)
    
    # Find configuration file with reads->preprocess attributes
    if 'postprocess' not in config.get('variant', {}):
        return
    
    # Add new read_type attribute to the main config, so rules for preprocess steps may deduce input read files
    previous_type = original_type
    if config['variant']['postprocess'] is not None:
        for this_type in config['variant']['postprocess'].keys():
            config['variant']['postprocess'][this_type]['input_variant_type'] = previous_type
            previous_type = this_type
    
    setattr(pipeline, 'postprocessed_variant_type', previous_type)
