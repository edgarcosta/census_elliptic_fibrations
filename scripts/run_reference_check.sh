#!/bin/sh
# Arguments: check path from the repository root, absolute result and log paths.
set -eu
[ "$#" -eq 3 ] || exit 2
magma_command=$(command -v "${MAGMA:-magma}")
magma_command=$(cd "$(dirname "$magma_command")" && pwd)/$(basename "$magma_command")
root=$(CDPATH= cd -- "$(dirname -- "$0")/.." && pwd)
check=$1
resultpath=$2
logpath=$3
cd "$root"
# The launcher enables error exits before the standalone check loads its inputs.
mkdir -p "$root/artifacts"
launcher=$(mktemp "$root/artifacts/reference-launch-XXXXXXXX")
trap 'rm -f "$launcher"' EXIT HUP INT TERM
python3 - "$launcher" "$check" <<'PYTHON'
import json
import sys
from pathlib import Path
Path(sys.argv[1]).write_text("SetQuitOnError(true);\nload " + json.dumps(sys.argv[2]) + ";\n")
PYTHON
"$magma_command" -b "checkpath:=$check" "resultpath:=$resultpath" "$launcher" > "$logpath" 2>&1 || exit 1
printf 'REFERENCE_CHECK_OK\n' | cmp -s - "$resultpath"
