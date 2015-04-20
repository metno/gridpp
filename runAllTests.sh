#!/bin/sh
# Test suite

# make test

set -x  # Output every command
for program in testing/*.exe; do
   "${program}" || exit $?
done
