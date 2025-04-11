# SPDX-FileCopyrightText: Contributors to technology-data <https://github.com/pypsa/technology-data>
#
# SPDX-License-Identifier: GPL-3.0-only

# coding: utf-8

.ONESHELL:

.PHONY: test unit-test

# Run default tests
test:
	set -e
	snakemake --cores all -f compile_cost_assumptions --configfile config.yaml
	snakemake --cores all -f compile_cost_assumptions_usa --configfile config.yaml
	echo "All tests completed successfully."

unit-test:
	pytest test
