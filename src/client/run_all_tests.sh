#!/bin/bash

set -x  # Output every command
for program in $(dirname $0)/src/client/*.exe; do
   "${program}" || exit $?
done
