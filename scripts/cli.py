# cli.py
import argparse

def parse_args():
    p = argparse.ArgumentParser(
        description="Create a Snakemake config.yaml for the ECHO repeatome pipeline."
    )

    sub = p.add_subparsers(dest="cmd", required=True)

    init = sub.add_parser("init", help="Create a new config.yaml")

    # required outputs
    init.add_argument("--output", required=True, help="Path to write config.yaml")

    # samples
    init.add_argument("--samples", nargs="+", required=True,
                      help="Sample IDs (space-separated)")

    # pipeline control
    init.add_argument("--start-from",
                      choices=["pod5", "ubam", "fastq", "bam"],
                      required=True)

    init.add_argument("--input-dir", required=True)
    init.add_argument("--output-dir", required=True)

    # reference genome
    init.add_argument("--reference", required=True,
                      help="Reference genome FASTA")
    init.add_argument("--reference-name",
                      choices=["GRCh38", "chm13v2"],
                      default="GRCh38")

    # catalog selection
    init.add_argument("--use-bundled-db", action="store_true",
                      help="Use ECHO-provided catalogs is resources/echoDB_v1")

    init.add_argument("--tr-catalog",
                      help="Custom TR catalog BED (disables ECHO defaults)")
    init.add_argument("--te-catalog",
                      help="Custom TE catalog BED (disables ECHO defaults)")

    init.add_argument("--tr-type",
                      choices=["genome-wide", "genome-wide-hipstr", "genome-wide-str-cpg", "genome-wide-vntr-cpg", "pathogenic", "forensic"])
    init.add_argument("--te-type",
                      choices=["all", "LINE", "SINE", "LTR", "DNA", "helitron", "retroposon"])

    # CpG filtering
    init.add_argument("--cpg-filter", action="store_true")
    init.add_argument("--no-cpg-filter", dest="cpg_filter", action="store_false")
    init.set_defaults(cpg_filter=None)

    init.add_argument("--cpg-str-bed",
                      help="CpG STR BED (required for custom catalogs + cpg-filter)")

    # read filtering
    init.add_argument("--min-read-quality", type=int, default=7)
    init.add_argument("--min-read-length", type=int, default=500)

    # analysis params
    init.add_argument("--flanking-length-bp", type=int, default=250)
    init.add_argument("--extension-repeat-consensus", type=int, default=1000)

    return p.parse_args()

