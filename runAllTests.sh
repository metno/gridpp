#!/bin/sh
# Test suite

# make test

set -x  # Output every command
for program in $(dirname $0)/testing/*.exe; do
   "${program}" || exit $?
done
