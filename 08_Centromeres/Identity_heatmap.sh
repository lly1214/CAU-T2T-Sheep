#!/bin/bash
set -ex

snakemake   --snakefile  StainedGlass/workflow/Snakefile   \
	--configfile  config.yaml  \
	 --keep-going    --cores 10   make_figures   --rerun-incomplete  