#!/usr/bin/env python3
"""Verify an extracted Supplementary File S6 against its SHA256 manifest."""
from pathlib import Path
import csv
import hashlib
import sys

root=Path(__file__).resolve().parent
manifest=root/'archive_manifest.tsv'
if not manifest.exists():
    raise SystemExit('Run the verifier from the extracted S6 archive, beside archive_manifest.tsv.')
with manifest.open(encoding='utf-8') as stream:
    records=list(csv.DictReader(stream,delimiter='\t'))
failed=[]
for r in records:
    file=root/r['archive_path']
    if not file.is_file():failed.append((r['archive_path'],'missing'));continue
    h=hashlib.sha256()
    with file.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    if file.stat().st_size!=int(r['bytes']) or h.hexdigest()!=r['sha256']:
        failed.append((r['archive_path'],'size or SHA256 differs'))
if failed:
    for name,reason in failed:print(name+': '+reason)
    sys.exit(1)
print(f'All {len(records)} files match the S6 manifest.')
