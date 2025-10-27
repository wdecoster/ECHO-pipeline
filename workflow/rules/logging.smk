# rules/logging.smk
rule create_log:
    output:
        LOGFILE
    shell:
        """
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

        echo "------------------- RESOURCE SUMMARY -------------------" >> {output}
        echo "CPUs on this node: $(nproc)" >> {output}
        echo "RAM for this node $(free -h | grep Mem | awk '{{print "total: " $2, "free: " $7}}')" >> {output}
        echo "=========================================================" >> {output}
        """

