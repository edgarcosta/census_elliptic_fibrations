#!/bin/sh
# Arguments: absolute script, working directory, spec, result file, log file.
set -eu
[ "$#" -eq 5 ] || exit 2
magma_command=$(command -v "${MAGMA:-magma}")
magma_command=$(cd "$(dirname "$magma_command")" && pwd)/$(basename "$magma_command")
script=$1
workdir=$2
specpath=$3
resultpath=$4
logpath=$5
cd "$workdir"
exec "$magma_command" -b "specpath:=$specpath" "resultpath:=$resultpath" "$script" > "$logpath" 2>&1
