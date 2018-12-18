Picard - Mark Duplicates
----------------------------

TEMP_DIR=$(mktemp -d {params.tmp_dir}/XXXXXXXXXXX)

picard MarkDuplicates \
    I={input.bam} \
    O={output.bam} \
    M={log.stat} \
    TMP_DIR=$TEMP_DIR \
    VALIDATION_STRINGENCY=SILENT \
>  {log.out} \
2> {log.err}

rm -r $TEMP_DIR

**Location**

- *Filepath:* <SnakeLines_dir>/rules/mapping/postprocess/deduplicated/picard.snake
- *Rule name:* picard__mark_duplicates

**Input(s):**


