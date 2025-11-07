# rules/log_generation.smk

rule create_log:
    output:
        LOGFILE
    run:
        import subprocess
        from datetime import datetime
        import socket
        import sys
        import snakemake
        import yaml
        import os
        from pathlib import Path


        def get_cluster_config_from_args():
            cluster_config = {}
            profile_path = None
            
            # First, try to get profile from workflow object
            if hasattr(workflow, 'overwrite_configfiles'):
                print(f"DEBUG: workflow.overwrite_configfiles = {workflow.overwrite_configfiles}")
            
            # Try workflow.default_remote_prefix (sometimes contains profile info)
            if hasattr(workflow, 'configfiles'):
                print(f"DEBUG: workflow.configfiles = {workflow.configfiles}")
            
            # Debug: print sys.argv
            print(f"DEBUG: sys.argv = {sys.argv}")
            
            # Parse command line arguments to find --profile
            args = sys.argv
            for i, arg in enumerate(args):
                if arg == "--profile" and i + 1 < len(args):
                    profile_path = Path(args[i + 1])
                    break
                elif arg.startswith("--profile="):
                    profile_path = Path(arg.split("=", 1)[1])
                    break
            
            print(f"DEBUG: profile_path from sys.argv = {profile_path}")
            
            # Alternative: try to find profile in common locations if not found
            if not profile_path:
                # Check if there's a profile mentioned in workflow attributes
                for attr in dir(workflow):
                    if 'profile' in attr.lower():
                        print(f"DEBUG: workflow.{attr} = {getattr(workflow, attr, 'N/A')}")
                
                # Fallback: try common profile locations
                potential_profiles = [
                    Path("profiles/brando_profile"),
                    Path(".snakemake/profiles/default"),
                ]
                for p in potential_profiles:
                    if (p / "config.yaml").exists():
                        profile_path = p
                        print(f"DEBUG: Found profile at fallback location: {profile_path}")
                        break
            
            # If profile path found, read config.yaml
            if profile_path:
                config_file = profile_path / "config.yaml"
                print(f"DEBUG: Looking for config file at: {config_file}")
                print(f"DEBUG: Config file exists: {config_file.exists()}")
                if config_file.exists():
                    with open(config_file) as f:
                        cluster_config = yaml.safe_load(f)
                    print(f"DEBUG: Successfully loaded cluster config")
            
            return cluster_config, str(profile_path) if profile_path else "Not specified"


        def get_git_info(cmd, default="NA"):
            """Helper to run git commands safely."""
            try:
                return subprocess.check_output(cmd, shell=True, text=True).strip()
            except subprocess.CalledProcessError:
                return default

        # Collect system info
        run_id = str(LOGTIMESTAMP)
        hostname = socket.gethostname()
        datetime_now = datetime.now().strftime("%Y-%m-%d %H:%M:%S")

        # Collect Git info
        git_version = get_git_info("git describe --tags --abbrev=0")
        git_commit = get_git_info("git rev-parse HEAD")
        git_branch = get_git_info("git rev-parse --abbrev-ref HEAD")
        git_remote = get_git_info("git config --get remote.origin.url")

        # Collect environment info
        python_version = sys.version.split()[0]  # e.g. "3.10.14"
        snakemake_version = snakemake.__version__

        # Get cluster config from command line profile
        cluster_config, profile_path = get_cluster_config_from_args()        

        with open(output[0], "w") as f:
            f.write("================ Pipeline Execution Log ================\n")
            f.write(f"Unique Run ID: {run_id}\n")
            f.write(f"Date & Time: {datetime_now}\n")
            f.write(f"Executed on Server: {hostname}\n\n")
            
            f.write("--------------- SOFTWARE ENVIRONMENT -------------------\n")
            f.write(f"Snakemake version: {snakemake_version}\n")
            f.write(f"Python version: {python_version}\n\n")

            f.write("--------------- GIT REPOSITORY INFO --------------------\n")
            f.write(f"ECHO version: {git_version}\n")
            f.write(f"ECHO git commit Hash: {git_commit}\n")
            f.write(f"ECHO git Branch: {git_branch}\n")
            f.write(f"ECHO github remote URL: {git_remote}\n\n")

            f.write("--------------- CONFIGURATION PARAMETERS ---------------\n")
            f.write(f"SAMPLES: {SAMPLES}\n")
            f.write(f"START_FROM: {START_FROM}\n")
            f.write(f"OUTPUT_DIR: {OUTPUT_DIR}\n")
            f.write(f"INPUT_DIR: {INPUT_DIR}\n")
            f.write(f"REFERENCE: {REFERENCE}\n")
            f.write(f"REFERENCE_TE: {REFERENCE_TE}\n")
            f.write(f"TR_CATALOG: {TR_CATALOG}\n")
            f.write(f"TE_CATALOG: {TE_CATALOG}\n")
            f.write(f"FLANKING_LENGTH_BP: {FLANKING_LENGTH_BP}\n")
            f.write(f"CONSENSUS_EXTENSION: {CONSENSUS_EXTENSION}\n")
            f.write(f"HAPLOID_CHRS: {HAPLOID_CHRS}\n")
            f.write(f"TYPE_OF_TE: {TYPE_OF_TE}\n")
            f.write(f"TYPE_OF_TR: {TYPE_OF_TR}\n")

            # Write cluster config section
            f.write("--------------- CLUSTER CONFIGURATION ------------------\n")
            if cluster_config:
                f.write(yaml.dump(cluster_config, default_flow_style=False, sort_keys=False))
            else:
                f.write("No cluster configuration found\n")
            f.write("\n")
 

    
