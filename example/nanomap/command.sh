echo "[INFO] --includingAllContigs enabled"

python /home/adog/Clair/clair/../clair.py callVarBam --chkpnt_fn "/home/adog/Clair/models/ont/model/model" --ref_fn "/home/adog/snakelines/example/nanomap/reference/sars_cov_2/sars_cov_2.fa" --bam_fn "/home/adog/snakelines/example/nanomap/mapping/sars_cov_2/original/example_1_sorted.bam" --threshold "0.2" --minCoverage "4" --pypy "pypy3" --samtools "samtools" --delay "10" --threads "4" --sampleName "SAMPLE" --ctgName "NC_045512.2" --ctgStart "0" --ctgEnd "29903" --call_fn "variant/tmp/a.NC_045512.2_0_29903.vcf"
