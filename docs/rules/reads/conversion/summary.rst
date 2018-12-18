Pigz - Unzip File
---------------------

set +euo pipefail
pigz \
    --decompress \
    --keep \
    --processes {threads} \
    {input.gzipped}

**Location**

- *Filepath:* <SnakeLines_dir>/rules/reads/conversion/unzip.snake
- *Rule name:* pigz__unzip_file

**Input(s):**


