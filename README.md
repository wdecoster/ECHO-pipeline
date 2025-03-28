To run the pipeline, input files must be in one of the following formats:

- `.pod5`
- `.ubam`
- `.bam`

---

## 🗂 Project Directory Structure

If you start the pipeline from `pod5`, ensure that the `pod5` files follow this directory structure:

```
/ifs/data/research/unique/projects/{project_name}/00_raw_data/pod5/{sample_name}
```

If you start the pipeline from `ubam`, the `ubam` file should be in this one:

```
/ifs/data/research/unique/projects/{project_name}/00_raw_data/basecalled/ubam/{sample_name}
```

If you start the pipeline from `bam`, the `bam` and the index `bai` files should be in this one:

```
/ifs/data/research/unique/projects/{project_name}/01_alignment/{sample_name}
```



Replace `{project_name}` and `{sample_name}` with your actual project and sample identifiers.

---

## ⚙️ Configuration File

Before running the pipeline, you need to create a `config.yaml` file that includes the following:

- Sample ID
- Input format (`.pod5`, `.ubam`, or `.bam`)
- Input directory
- Reference genome path
- TE catalog path
- TR catalog path
- Length of flanking regions for TE and TR analysis

> An example `config.yaml` file is provided. Copy it and customize it for your own analysis.

---

## 🚀 Running the Pipeline

To run the pipeline, use the following command:

```bash
snakemake --use-conda --use-singularity --jobs 4 --configfile config.yaml --cluster-config cluster-config.yaml --cluster "sbatch --partition={cluster.partition} --mem={cluster.mem} --cpus-per-task={cluster.cpus} --time={cluster.time}" --conda-frontend conda --singularity-args "-B /ifs/data/research/unique"
```

---

## 🔁 Restarting After Interruption

If the pipeline is interrupted, you can resume from where it left off using the `--rerun-incomplete` flag:

```bash
snakemake --use-conda --use-singularity --jobs 4 --configfile config.yaml --cluster-config cluster-config.yaml --cluster "sbatch --partition={cluster.partition} --mem={cluster.mem} --cpus-per-task={cluster.cpus} --time={cluster.time}" --conda-frontend conda --singularity-args "-B /ifs/data/research/unique" --rerun-incomplete
```

