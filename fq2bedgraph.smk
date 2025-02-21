configfile: "config/runtime_config.yaml"

from utils.sample_sheet_parser import generate_file_ls, generate_tool_ls

rule all:
    input:
        generate_file_ls(config['sample_sheet'])

tool_dict = generate_tool_ls(config['sample_sheet'])

if tool_dict["TRIMMER"]:
    for tool in tool_dict["TRIMMER"]:
        include: f"rules/{tool}/trim.smk"

if tool_dict["QC_REPORTER"]:
    for tool in tool_dict["QC_REPORTER"]:
        include: f"rules/{tool}/qc.smk"

if tool_dict["ALIGNER"]:
    for tool in tool_dict["ALIGNER"]:
        if 'bismark' in tool:
            include: f"rules/bismark/align.smk"
        else:
            include: f"rules/{tool}/align.smk"

if tool_dict["DEDUPER"]:
    for tool in tool_dict["DEDUPER"]:
        if tool != 'no_dup':
            include: f"rules/{tool}/dedup.smk"
        else:
            include: f"rules/misc/skip_dup.smk"

if tool_dict["COUNTER"]:
    for tool in tool_dict["COUNTER"]:
        include: f"rules/{tool}/count.smk"

if any(tool_dict["RECALIBRATE"]):
    include: "rules/gatk/bqsr.smk"

if any(tool_dict["STATS"]):
    for tool in ['astair', 'samtools', 'bismark', 'methyldackel']:
        include: f"rules/{tool}/stats.smk"

    if any(tool_dict["RECALIBRATE"]) and not all(tool_dict["RECALIBRATE"]):
        include: "rules/qualimap/qualimap_bqsr.smk"
        include: "rules/qualimap/qualimap.smk"

    if not any(tool_dict["RECALIBRATE"]):
        include: "rules/qualimap/qualimap.smk"
    elif all(tool_dict["RECALIBRATE"]):
        include: "rules/qualimap/qualimap_bqsr.smk"
    else:
        include: "rules/qualimap/qualimap_bqsr.smk"
        include: "rules/qualimap/qualimap.smk"
