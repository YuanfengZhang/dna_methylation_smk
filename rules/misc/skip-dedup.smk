rule fake_dedup:
    input:
        bam = "result/{BaseName}/{DedupParentDir}/{BaseName}.bam",
        bai = "result/{BaseName}/{DedupParentDir}/{BaseName}.bam.bai"
    output:
        bam = "result/{BaseName}/{DedupParentDir}/no-dedup/{BaseName}.bam",
        bai = "result/{BaseName}/{DedupParentDir}/no-dedup/{BaseName}.bam.bai"
    threads: 1
    run:
        from pathlib import Path


        input_ls = [input.bam, input.bai]
        output_ls = [output.bam, output.bai]

        for original_path_str, target_path_str in zip(input_ls, output_ls):
            original_path = Path(original_path_str)
            target_path = Path(target_path_str)

            if not original_path.exists(follow_symlinks=True):
                raise FileNotFoundError(f"{original_path} does not exist")

            real_path = original_path.resolve()
            target_path.parent.mkdir(parents=True, exist_ok=True)
            target_path.symlink_to(real_path)
