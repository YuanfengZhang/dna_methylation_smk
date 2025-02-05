from utils.sample_sheet_parser import read_sample_sheet, get_required_tools

configfile: "config/runtime_config.yaml"

rule all:
    input:
        read_sample_sheet(configfile['sample_sheet'])
