#!/usr/bin/env python3
import argparse
import fcntl
import json
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
from pipeline.sha2_evidence import DEFAULT, Sha2Evidence, digest, encoded
from pipeline.sha2 import require


def immutable(path, data):
    if path.exists():
        require(path.read_bytes() == data, f'changed immutable output: {path.name}')
    else:
        temporary = path.with_suffix(path.suffix + '.new')
        temporary.write_bytes(data)
        temporary.replace(path)


def main():
    parser = argparse.ArgumentParser(description='Import and verify saved F_2 Sha[2] refinements.')
    parser.add_argument('command', choices=['import-sha2', 'verify-sha2'])
    parser.add_argument('--output', type=Path, required=True)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    with (args.output / '.lock').open('w') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX | fcntl.LOCK_NB)
        expected_manifest = (DEFAULT / 'manifest.json').read_bytes()
        if args.command == 'import-sha2':
            archive = (DEFAULT / 'evidence.tar.gz').read_bytes()
        else:
            require((args.output / 'source-manifest.json').read_bytes() == expected_manifest,
                    'changed source manifest')
            archive = (args.output / 'evidence.tar.gz').read_bytes()
        source = Sha2Evidence(archive, expected_manifest)
        records = source.records()
        output = b''.join(encoded(r) for r in records)
        summary = source.summary(records)
        summary['records_sha256'] = digest(output)
        files = {'evidence.tar.gz': archive, 'source-manifest.json': expected_manifest,
                 'import-producer.json': encoded(source.producer), 'refinements.jsonl': output,
                 'import.json': encoded(summary)}
        for name, data in files.items():
            path = args.output / name
            if args.command == 'import-sha2':
                immutable(path, data)
            else:
                require(path.read_bytes() == data, f'changed import output: {name}')
        print(json.dumps(summary, sort_keys=True))


if __name__ == '__main__':
    try:
        main()
    except (ValueError, KeyError, OSError) as error:
        sys.exit(str(error))
