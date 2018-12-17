import tempfile
if not 'tmp_dir' in config:
    config['tmp_dir'] = tempfile.gettempdir()