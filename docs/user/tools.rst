Additional tools
================

Although Snakelines installs required tools automatically, some may offer installation or configuration of addtional components, for example plugins.

Variant effect predictor (VEP)
------------------------------
Installation is based on `<https://www.ensembl.org/info/docs/tools/vep/script/vep_download.html#installer>`_
Localize the conda installation directory, for example /docker/tools/miniconda3/envs/ensembl-vep-101.0/bin/, and use the following commands to install VEP cache.

::

    vep_install -a c --species homo_sapiens --assembly GRCh38 --cachedir /data/genome/human/grch38/annotation/vep
    vep_convert_cache --dir /data/genome/human/grch38/annotation/vep/ --species all --version all --compress "gzip -dc"

Use the following commands as reference to install the dbNSFP plugin.

::

    vep_install -a p --PLUGINS dbNSFP
    wget ftp://dbnsfp:dbnsfp@dbnsfp.softgenetics.com/dbNSFP4.1c.zip
    unzip dbNSFP4.1c.zip
    gzip -dc dbNSFP4.1c_variant.chr1.gz | head -n 1 | bgzip > header.gz
    cat dbNSFP4.1c_variant.chr{1..22}.gz dbNSFP4.1c_variant.chrX.gz dbNSFP4.1c_variant.chrY.gz dbNSFP4.1c_variant.chrM.gz | zgrep -v '#chr' | bgzip -@ 8 > dbNSFP4.1c.noheader.gz
    cat header.gz dbNSFP4.1c.noheader.gz > dbNSFP4.1c.gz
    tabix -s 1 -b 2 -e 2 dbNSFP4.1c.gz

