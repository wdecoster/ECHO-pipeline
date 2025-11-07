# rules/log_generation.smk



rule create_log:
    output:
        LOGFILE
    shell:
        r"""
        set -euo pipefail

        echo "================ Pipeline Execution Log ================" > {output}
        echo "Unique Run ID: {LOGTIMESTAMP}" >> {output}
        echo "Date & Time: $(date)" >> {output}
        echo "Executed on Server: $(hostname)" >> {output}
        echo "" >> {output}

        echo "--------------- CONFIGURATION PARAMETERS ---------------" >> {output}
        echo "SAMPLES: {SAMPLES}" >> {output}
        echo "START_FROM: {START_FROM}" >> {output}
        echo "OUTPUT_DIR: {OUTPUT_DIR}" >> {output}
        echo "INPUT_DIR: {INPUT_DIR}" >> {output}
        echo "REFERENCE: {REFERENCE}" >> {output}
        echo "REFERENCE_TE: {REFERENCE_TE}" >> {output}
        echo "TR_CATALOG: {TR_CATALOG}" >> {output}
        echo "TE_CATALOG: {TE_CATALOG}" >> {output}
        echo "FLANKING_LENGTH_BP: {FLANKING_LENGTH_BP}" >> {output}
        echo "CONSENSUS_EXTENSION: {CONSENSUS_EXTENSION}" >> {output}
        echo "HAPLOID_CHRS: {HAPLOID_CHRS}" >> {output}
        echo "TYPE_OF_TE: {TYPE_OF_TE}" >> {output}
        echo "TYPE_OF_TR: {TYPE_OF_TR}" >> {output}

        echo "" >> {output}
        echo "------------------- SOFTWARE VERSION -------------------" >> {output}
        GIT_TAG=$(git describe --tags --abbrev=0 2>/dev/null || echo "NA")
        GIT_COMMIT_HASH=$(git rev-parse HEAD 2>/dev/null || echo "NA")
        GIT_BRANCH=$(git rev-parse --abbrev-ref HEAD 2>/dev/null || echo "NA")
        GIT_REMOTE_URL=$(git config --get remote.origin.url 2>/dev/null || echo "NA")
        
        echo "snakemake version:" $(snakemake --version) >> {output}
        echo "python version:" $(python --version) >> {output}
        echo "ECHO pipeline" >> {output}
        echo "Git tag:          $GIT_TAG" >> {output}
        echo "Git commit hash:  $GIT_COMMIT_HASH" >> {output}
        echo "Git branch:       $GIT_BRANCH" >> {output}
        echo "Git remote URL:   $GIT_REMOTE_URL" >> {output}
        echo "=========================================================" >> {output}
        echo "" >> {output} 

        echo "------------------- CLUSTER VERSION -------------------" >> {output}
        echo "Cluster config file: {CLUSTER_CONFIG}" >> {output}
       """


