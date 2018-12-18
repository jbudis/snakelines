Qualimap - Mapping Quality Report Accross Reference
-------------------------------------------------------

qualimap bamqc \
    --java-mem-size=100G \
    -bam {input.bam} \
    --paint-chromosome-limits \
    -outdir {params.outdir} \
    -outformat PDF:HTML \
    -nt {threads} \
>  {log.out} \
2> {log.err}

**Location**

- *Filepath:* <SnakeLines_dir>/rules/mapping/report/quality_report/qualimap.snake
- *Rule name:* qualimap__mapping_quality_report_accross_reference

**Input(s):**


Qualimap - Mapping Quality Report Accross Panel
---------------------------------------------------

qualimap bamqc \
    --java-mem-size=100G \
    -bam {input.bam} \
    --feature-file {input.bed} \
    -outdir {params.outdir} \
    -outformat PDF:HTML \
    -nt {threads} \
>  {log.out} \
2> {log.err}

**Location**

- *Filepath:* <SnakeLines_dir>/rules/mapping/report/quality_report/qualimap.snake
- *Rule name:* qualimap__mapping_quality_report_accross_panel

**Input(s):**


Qualimap - Summarize Quality Reports
----------------------------------------

# Prepare sample list file for the qualimap
## Column1: Sample name
## Column2: Path to precomputed bamqc report

for REPORT in {input.reports}; do
    REPORT_DIR=`dirname $REPORT`
    SAMPLE=`basename $REPORT_DIR`
    echo -e "$SAMPLE\t$REPORT_DIR" >> {params.sample_list}
done

# Run qualimap on samples specified in sample list
qualimap multi-bamqc \
    -d {params.sample_list} \
    -outdir {params.out_dir} \
    -outformat PDF:HTML \
>  {log.out} \
2> {log.err}

**Location**

- *Filepath:* <SnakeLines_dir>/rules/mapping/report/quality_report/qualimap.snake
- *Rule name:* qualimap__summarize_quality_reports

**Input(s):**


