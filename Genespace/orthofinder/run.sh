# workflow adapted from https://gitlab.com/salk-tm/snake_orthofinder

### run orthofinder
nohup snakemake -p -j 1 --use-conda --tibanna \
    --snakefile Snakefile \
    --default-remote-prefix path/to/orthofinder/dir/ \
    --default-resources mem_mb=30000 disk_mb=1000000 \
    --tibanna-config root_ebs_size=16 log_bucket=log_dir \
    --config \
        proteins=myproteins.tar.gz \
        outbase=orthofinder \
        threads=32 \
    > orthofinder.log 2>&1 &
