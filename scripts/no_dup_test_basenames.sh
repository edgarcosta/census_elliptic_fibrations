#!/bin/sh
# Make's check/%.m and debug/%.m targets share a namespace across directories.
set -eu
dups=$(find Tests examples -type f -name '*.m' -exec basename {} \; | sort | uniq -d)
if [ -n "$dups" ]; then
    echo "Duplicate test basenames across Tests/ and examples/:"
    echo "$dups"
    exit 1
fi
