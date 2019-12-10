#!/usr/bin/env bash
# vim: set noexpandtab tabstop=2:

set -v
../src/pefilter/pefilter -i LC1_chr_1k.bam -s -t 4
tmpdir=$(mktemp -d)
../src/pefilter/pefilter -i LC1_chr_1k.bam -o "$tmpdir/outfile.bam" -t 4
tree "$tmpdir"
