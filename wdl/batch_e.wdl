import "https://raw.githubusercontent.com/broadinstitute/hailrunner/refs/heads/main/hailrunner/wdl/hailrunner_run.wdl" as hailrunner

struct IntervalSpec {
    String name
    String path
}

workflow batch_e {
    input {
        # ---- batch_e analysis inputs (required) ----
        String ancestry_tsv
        String comparison_tsv
        String comparison_col
        Array[IntervalSpec] intervals
        String output_dir

        # ---- batch_e analysis inputs (optional) ----
        String data_source = "mt"
        String? vcf_path
        String? mt_path
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
    # Build script_args array from typed inputs
    # ================================================================

    # Required args -- always present
    Array[String] required_args = [
        "--ancestry-tsv", ancestry_tsv,
        "--comparison-tsv", comparison_tsv,
        "--comparison-col", comparison_col,
        "--data-source", data_source,
        "--output-dir", output_dir
    ]

    # Interval args (scatter over struct array, flatten pairs)
    scatter (spec in intervals) {
        Array[String] one_interval = ["--interval", spec.name + "=" + spec.path]
    }
    Array[String] interval_args = flatten(one_interval)

    # Optional scalar string args
    Array[String] mt_path_args = if defined(mt_path) then ["--mt-path", select_first([mt_path])] else []
    Array[String] vcf_path_args = if defined(vcf_path) then ["--vcf-path", select_first([vcf_path])] else []
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
        mt_path_args,
        vcf_path_args,
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
