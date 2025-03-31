rule fake_pre_calibrate:
    input:
        r1 = "result/{BaseName}/{DoubleSkipParentDir}/{BaseName}.R1.fq.gz",
        r2 = "result/{BaseName}/{DoubleSkipParentDir}/{BaseName}.R2.fq.gz"
    output:
        r1 = "result/{BaseName}/{DoubleSkipParentDir}/no-pre-dedup/no-pre-calibration/{BaseName}.R1.fq.gz",
        r2 = "result/{BaseName}/{DoubleSkipParentDir}/no-pre-dedup/no-pre-calibration/{BaseName}.R2.fq.gz"
    threads: 1
    run:
        from pathlib import Path


        input_ls = [input.r1, input.r2]
        output_ls = [output.r1, output.r2]

        for original_path_str, target_path_str in zip(input_ls, output_ls):
            original_path = Path(original_path_str)
            target_path = Path(target_path_str)

            if not original_path.exists(follow_symlinks=True):
                raise FileNotFoundError(f"{original_path} does not exist")

            real_path = original_path.resolve()
            target_path.parent.mkdir(parents=True, exist_ok=True)
            target_path.symlink_to(real_path)
