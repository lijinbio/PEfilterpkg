#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

set -v
tmpfile=$(mktemp --suffix=.bam)
../src/pefiltertrad/pefiltertrad LC1_chr_1k.bam "$tmpfile"
