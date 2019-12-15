#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

set -v
tmpdir=$(mktemp -d)
../src/pefiltertag/pefiltertag -i LC1_chr_1k.bam -o "$tmpdir/outfile.bam" -t 4
tree "$tmpdir"
