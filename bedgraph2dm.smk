configfile: "config/runtime_config.yaml"

from utils.sample_sheet_parser import (bedgraph2dm_file_ls,
                                       bedgraph2dm_tool_ls)

rule all:
    input:
        bedgraph2dm_file_ls(config['sample_sheet'])


for tool in bedgraph2dm_tool_ls(config['sample_sheet']):
    include: f"rules/{tool}/dm.smk"

include: "fq2bedgraph.smk"
