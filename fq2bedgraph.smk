configfile: "config/runtime_config.yaml"
wildcard_constraints:
    BaseName=r"[A-Z]{2}\d_[A-Z]{2}_\d"


from utils.sample_sheet_parser import (fq2bedgraph_file_ls,
                                       fq2bedgraph_tool_ls)

to_generate_files = fq2bedgraph_file_ls(config['sample_sheet'])

rule all:
    input:
        to_generate_files

tool_dict = fq2bedgraph_tool_ls(config['sample_sheet'])

if tool_dict["TRIMMER"]:
    for tool in tool_dict["TRIMMER"]:
        include: f"rules/{tool}/trim.smk"

if tool_dict["FASTQCER"]:
    for tool in tool_dict["FASTQCER"]:
        include: f"rules/{tool}/qc.smk"

if tool_dict["ALIGNER"]:
    for tool in tool_dict["ALIGNER"]:
        if 'bismark' in tool:
            include: f"rules/bismark/align.smk"
        elif 'msuite2' in tool:
            include: f"rules/msuite2/align.smk"
        else:
            include: f"rules/{tool}/align.smk"

if tool_dict["DEDUPER"]:
    for tool in tool_dict["DEDUPER"]:
        if tool != 'no-dedup':
            if tool == 'gatk-dedup':
                include: "rules/gatk/dedup.smk"
            else:
                include: f"rules/{tool}/dedup.smk"

if tool_dict["COUNTER"]:
    for tool in tool_dict["COUNTER"]:
        include: f"rules/{tool}/count.smk"

if tool_dict["CALIBRATOR"]:
    for tool in tool_dict["CALIBRATOR"]:
        if tool != 'no-calibrate':
            if tool == 'gatk-cali':
                include: "rules/gatk/calibrate.smk"
            else:
                include: f"rules/{tool}/calibrate.smk"

if any(tool_dict["BAMStatist"]):
    for tool in ['samtools', 'bismark', 'methyldackel',
                 'biscuit', 'qualimap']:
        include: f"rules/{tool}/stats.smk"

if any('no-pre-dedup' in file for file in to_generate_files):
    include: "rules/misc/skip-pre-dedup.smk"
if any('no-pre-calibrate' in file for file in to_generate_files):
    include: "rules/misc/skip-pre-calibrate.smk"
if any('no-dedup' in file for file in to_generate_files):
    include: "rules/misc/skip-dedup.smk"
if any('no-calibrate' in file for file in to_generate_files):
    include: "rules/misc/skip-calibrate.smk"
