version 1.0

import "https://raw.githubusercontent.com/broadinstitute/hailrunner/refs/heads/main/hailrunner/wdl/hailrunner_run.wdl" as hailrunner

workflow batch_e {
    input {
        # ---- batch_e analysis inputs (required) ----
        String input_path    # VCF glob or MatrixTable path (auto-detected)
        String ancestry_tsv
        String comparison_tsv
        String comparison_col
        String output_dir

        # ---- batch_e analysis inputs (optional) ----
        Array[String] intervals = [
            "ACMG59=https://raw.githubusercontent.com/broadinstitute/batch-e/refs/heads/main/intervals/acmg59_allofus_19dec2019.GRC38.wGenes.NEW.bed.gz",
            "Low_Mappability=https://raw.githubusercontent.com/broadinstitute/batch-e/refs/heads/main/intervals/GRCh38_lowmappabilityall.bed.gz",
            "GC_gt_85=https://raw.githubusercontent.com/broadinstitute/batch-e/refs/heads/main/intervals/GRCh38_gc85_slop50.bed.gz",
            "GC_lt_25=https://raw.githubusercontent.com/broadinstitute/batch-e/refs/heads/main/intervals/GRCh38_gclt25_merged.bed.gz",
            "HighConf_Genome=https://raw.githubusercontent.com/broadinstitute/batch-e/refs/heads/main/intervals/giab_highconf_wgs_calling_regions_hg38_intersection.downsampled.bed.gz"
        ]
        String? ancestry_col
        Array[String]? ancestries
        String? comparison_name
        Array[String]? comparison_values
        String? sample_id_col
        String? samples_per_group
        Boolean no_cache = false
        Boolean force_reimport = false

        # ---- Script URLs ----
        String batch_e_script = "https://raw.githubusercontent.com/broadinstitute/batch-e/refs/heads/main/batch_e/batch_e.py"
        String reporter_script = "https://raw.githubusercontent.com/broadinstitute/batch-e/refs/heads/main/batch_e/batch_e_reporter.py"

        # ---- hailrunner cluster config ----
        String staging_bucket
        String? project
        String region = "us-central1"
        String subnet = "subnetwork"
        Int workers = 32
        Int preemptibles = 0
        String worker_type = "n1-highmem-8"
        String driver_type = "n1-highmem-32"
        Int worker_disk_gb = 300
        Int driver_disk_gb = 500
        String disk_type = "pd-standard"
        Int max_idle = 60
        Int max_age = 1440
        String? cluster_name
        String? service_account

        # ---- hailrunner spark config ----
        Int executor_cores = 4
        String executor_memory = "26g"
        Int driver_cores = 4
        String driver_memory = "26g"

        # ---- hailrunner WDL runtime ----
        String hailrunner_memory = "4GB"
        Int hailrunner_disk_size_gb = 50
        Int hailrunner_cpu = 2
        String hailrunner_docker = "us-docker.pkg.dev/broad-dsde-methods/hailrunner/hailrunner:0.1.0"

        # ---- Reporter inputs ----
        String report_title = "Batch Effect Report"
        Boolean report_no_sample_stats = false
        Float report_effect_threshold = 0.5
        String reporter_docker = "python:3.11-slim"
        String reporter_memory = "8GB"
        Int reporter_disk_gb = 50
    }

    # ================================================================
    # Stage intervals (download HTTPS → GCS if needed, GCS passthrough)
    # ================================================================
    scatter (spec in intervals) {
        call stage_interval {
            input:
                interval_spec = spec,
                staging_bucket = staging_bucket
        }
    }

    # ================================================================
    # Build script_args array from typed inputs
    # ================================================================

    # Required args -- always present
    Array[String] required_args = [
        "--input-path", input_path,
        "--ancestry-tsv", ancestry_tsv,
        "--comparison-tsv", comparison_tsv,
        "--comparison-col", comparison_col,
        "--output-dir", output_dir
    ]

    # Interval args from staged specs
    scatter (staged in stage_interval.staged_spec) {
        Array[String] one_interval = ["--interval", staged]
    }
    Array[String] interval_args = flatten(one_interval)

    # Optional scalar string args
    Array[String] ancestry_col_args = if defined(ancestry_col) then ["--ancestry-col", select_first([ancestry_col])] else []
    Array[String] comparison_name_args = if defined(comparison_name) then ["--comparison-name", select_first([comparison_name])] else []
    Array[String] sample_id_col_args = if defined(sample_id_col) then ["--sample-id-col", select_first([sample_id_col])] else []
    Array[String] samples_per_group_args = if defined(samples_per_group) then ["--samples-per-group", select_first([samples_per_group])] else []

    # Optional array args (nargs='+')
    Array[String] ancestries_args = if defined(ancestries) then flatten([["--ancestries"], select_first([ancestries])]) else []
    Array[String] comparison_values_args = if defined(comparison_values) then flatten([["--comparison-values"], select_first([comparison_values])]) else []

    # Boolean flags
    Array[String] no_cache_args = if no_cache then ["--no-cache"] else []
    Array[String] force_reimport_args = if force_reimport then ["--force-reimport"] else []

    # Combine all into final script_args
    Array[String] all_script_args = flatten([
        required_args,
        interval_args,
        ancestry_col_args,
        comparison_name_args,
        sample_id_col_args,
        samples_per_group_args,
        ancestries_args,
        comparison_values_args,
        no_cache_args,
        force_reimport_args
    ])

    # ================================================================
    # Build output_specs for hailrunner (delocalize results from GCS)
    # ================================================================
    Array[String] output_specs = [
        output_dir + "/pairwise_comparisons.tsv:pairwise_comparisons.tsv",
        output_dir + "/group_summaries.tsv:group_summaries.tsv",
        output_dir + "/sample_stats.tsv:sample_stats.tsv",
        output_dir + "/config.json:config.json",
        output_dir + "/timing_metrics.json:timing_metrics.json"
    ]

    # ================================================================
    # Phase 1: Run batch_e analysis on Dataproc via hailrunner
    # ================================================================
    call hailrunner.hailrunner_run as analysis {
        input:
            project          = project,
            staging_bucket   = staging_bucket,
            script           = batch_e_script,
            script_args      = all_script_args,
            output_specs     = output_specs,
            region           = region,
            subnet           = subnet,
            workers          = workers,
            preemptibles     = preemptibles,
            worker_type      = worker_type,
            driver_type      = driver_type,
            worker_disk_gb   = worker_disk_gb,
            driver_disk_gb   = driver_disk_gb,
            disk_type        = disk_type,
            max_idle         = max_idle,
            max_age          = max_age,
            cluster_name     = cluster_name,
            service_account  = service_account,
            executor_cores   = executor_cores,
            executor_memory  = executor_memory,
            driver_cores     = driver_cores,
            driver_memory    = driver_memory,
            total_memory     = hailrunner_memory,
            disk_size_gb     = hailrunner_disk_size_gb,
            cpu              = hailrunner_cpu,
            docker_image     = hailrunner_docker
    }

    # ================================================================
    # Phase 2: Generate HTML report (lightweight, no Dataproc)
    # ================================================================
    call generate_report {
        input:
            analysis_outputs    = analysis.output_files,
            results_dir         = output_dir,
            reporter_script_url = reporter_script,
            title               = report_title,
            no_sample_stats     = report_no_sample_stats,
            effect_threshold    = report_effect_threshold,
            docker              = reporter_docker,
            memory              = reporter_memory,
            disk_gb             = reporter_disk_gb
    }

    output {
        Array[File] analysis_files = analysis.output_files
        File report_html = generate_report.report_html
    }
}

# ================================================================
# Stage interval files: download HTTPS → GCS, pass through GCS paths
# ================================================================
task stage_interval {
    input {
        String interval_spec   # "NAME=PATH" or "NAME=https://..."
        String staging_bucket
        String docker = "google/cloud-sdk:slim"
    }

    command <<<
        set -euo pipefail
        spec="~{interval_spec}"
        name="${spec%%=*}"
        path="${spec#*=}"

        if [[ "$path" == http://* ]] || [[ "$path" == https://* ]]; then
            gcs_dest="~{staging_bucket}/staged_intervals/${name}.bed.gz"
            curl -fsSL "$path" -o interval_file
            gsutil cp interval_file "$gcs_dest"
            echo "${name}=${gcs_dest}" > result.txt
        else
            echo "$spec" > result.txt
        fi
    >>>

    output {
        String staged_spec = read_string("result.txt")
    }

    runtime {
        docker: docker
        memory: "2GB"
        disks: "local-disk 10 SSD"
        cpu: 1
    }
}

task generate_report {
    input {
        Array[File] analysis_outputs
        String results_dir
        String reporter_script_url
        String title = "Batch Effect Report"
        Boolean no_sample_stats = false
        Float effect_threshold = 0.5
        String docker = "python:3.11-slim"
        String memory = "8GB"
        Int disk_gb = 50
    }

    command <<<
        set -euo pipefail

        pip install --quiet pandas numpy matplotlib seaborn gcsfs scikit-learn

        curl -fsSL "~{reporter_script_url}" -o batch_e_reporter.py

        python3 batch_e_reporter.py \
            "~{results_dir}" \
            -o report.html \
            --title "~{title}" \
            ~{if no_sample_stats then "--no-sample-stats" else ""} \
            --effect-threshold ~{effect_threshold}
    >>>

    output {
        File report_html = "report.html"
    }

    runtime {
        docker: docker
        memory: memory
        disks: "local-disk ~{disk_gb} SSD"
        cpu: 2
    }
}
