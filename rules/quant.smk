def get_fq(wildcards):
    if config["trimming"]["skip"]:
        # no trimming, use raw reads
        print(units.loc[(wildcards.sample, wildcards.unit), "dir" + ["fq1","fq2"]].dropna())
        return units.loc[(wildcards.sample, wildcards.unit), "dir" + ["fq1","fq2"]].dropna()
    else:
        # yes trimming, use trimmed data
        if not is_single_end(**wildcards):
            # paired-end sample
            return expand("trimmed/{sample}-{unit}.{group}.fastq.gz",
                          group=[1, 2], **wildcards)
        # single end sample
        return "trimmed/{sample}-{unit}.fastq.gz".format(**wildcards)
            

rule quant:
    input:
        sample=get_fq
    output:
        counts="kallisto/{sample}-{unit}/abundance.tsv",
        dir="kallisto/{sample}-{unit}/"
        
    log:
        "logs/kallisto/{sample}-{unit}.log"
    params:
        index=config["ref"]["index"],
        # optional parameters
        extra="--threads 1 "
    threads: 1
    shell: "module load kallisto ; kallisto quant -i {params.index} -o {output} {input} &> {log}"
