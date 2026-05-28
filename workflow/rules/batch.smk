if config["batch_correction"]["enabled"]:

    rule combatseq:
        input:
            counts  = "results/counts/raw_counts.tsv",
            samples = config["samples"],
        output:
            corrected = "results/batch/corrected_counts.tsv",
        params:
            batch_col = config["batch_correction"]["batch_column"],
        conda: "../../envs/r.yaml"
        log: "logs/batch/combatseq.log"
        script:
            "../../scripts/run_combatseq.R"
