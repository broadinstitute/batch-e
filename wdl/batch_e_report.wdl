version 1.0

import "https://raw.githubusercontent.com/broadinstitute/batch-e/refs/heads/main/wdl/batch_e.wdl" as be

workflow batch_e_report {
    input {
        String results_dir
        String title = "Batch Effect Report"
        Boolean no_sample_stats = false
        Float effect_threshold = 0.5
        String batch_e_docker = "us-docker.pkg.dev/broad-dsde-methods/batch-e/batch-e:latest"
        String memory = "8GB"
        Int disk_gb = 50
    }

    call be.generate_report {
        input:
            analysis_outputs = [],
            results_dir      = results_dir,
            title            = title,
            no_sample_stats  = no_sample_stats,
            effect_threshold = effect_threshold,
            docker           = batch_e_docker,
            memory           = memory,
            disk_gb          = disk_gb
    }

    output {
        File report_html = generate_report.report_html
    }
}
