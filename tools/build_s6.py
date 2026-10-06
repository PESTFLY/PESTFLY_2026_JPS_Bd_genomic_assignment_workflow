#!/usr/bin/env python3
"""Build Supplementary File S6 from the checked repository using Python's standard library."""
from pathlib import Path
import argparse
import collections
import csv
import hashlib
import json
import zipfile

ROOT=Path(__file__).resolve().parents[1]
GUIDES=ROOT/'supplementary/S6'

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()

def table(path,rows,headers):
    with path.open('w',newline='',encoding='utf-8') as f:
        w=csv.DictWriter(f,fieldnames=headers,delimiter='\t',lineterminator='\n')
        w.writeheader();w.writerows(rows)

def inventory():
    candidates=[]
    allowed_docs={'CONCEPTS.md','INSTALL.md','METADATA.md','REPORTING_CRITERIA.md',
                  'S6_SOURCE_MAP.md','cli_parameters.tsv','runtime_version_evidence.tsv',
                  'reporting_criteria_diagnostic.tsv','location_index_changes.tsv',
                  'textual_output_changes.tsv','NATIVE_R_VERIFICATION.md',
                  'licensing_decision.json'}
    for p in sorted(ROOT.rglob('*')):
        if not p.is_file():continue
        rel=p.relative_to(ROOT)
        if rel.parts[0] in {'steps','results','data','environment'}:
            candidates.append((p,rel.as_posix()))
        elif rel.parts[0]=='docs' and (p.name in allowed_docs or rel.parts[1]=='native_R_verification'):
            candidates.append((p,rel.as_posix()))
        elif rel.parts[0]=='tools' and p.name=='check_implementation.R':
            candidates.append((p,rel.as_posix()))
        elif len(rel.parts)==1 and p.name in {'run_public_pipeline.R','run_validation_pipeline.R','verify_release.R','CITATION.cff','LICENSE','LICENSE_CONTENT.txt','LICENSING.md','ALIGNMENT_ACCESS.md'}:
            candidates.append((p,rel.as_posix()))
    # Only return_bundle export copies are candidates for removal. Matching is
    # confined to the same supplementary analysis and requires the same SHA256.
    canonical=collections.defaultdict(list)
    hashed={}
    for p,rel in candidates:
        h=sha(p);hashed[rel]=h
        parts=Path(rel).parts
        if len(parts)>3 and parts[:2]==('results','validation') and 'return_bundle' not in parts:
            canonical[(parts[2],h)].append(rel)
    aliases=[];retained=[]
    for p,rel in candidates:
        parts=Path(rel).parts
        matches=canonical.get((parts[2],hashed[rel]),[]) if len(parts)>3 and parts[:2]==('results','validation') and 'return_bundle' in parts else []
        if matches:
            target=min(matches,key=lambda x:(len(Path(x).parts),len(x),x))
            aliases.append({'omitted_export':rel,'retained_file':target,'bytes':p.stat().st_size,'sha256':hashed[rel]})
        else:retained.append((p,rel))
    for name in ['README.md','S6_reader_guide.pdf','source_map.tsv','verify_archive.py']:
        p=GUIDES/name
        if not p.exists():raise FileNotFoundError(f'Missing S6 guide resource: {p}')
        retained.append((p,name))
    return retained,aliases

def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--output',type=Path,default=ROOT.parent/'PESTFLY_Supplementary_File_S6.zip')
    args=ap.parse_args()
    retained,aliases=inventory()
    table(GUIDES/'duplicate_exports.tsv',aliases,['omitted_export','retained_file','bytes','sha256'])
    retained.append((GUIDES/'duplicate_exports.tsv','duplicate_exports.tsv'))
    counts=collections.Counter()
    for _,rel in retained:
        parts=Path(rel).parts
        if len(parts)>2 and parts[:2]==('results','validation'):counts[parts[2]]+=1
    info={'archive_role':'Supplementary File S6 with seven supplementary analyses and baseline resources',
          'assembly_date':'2026-10-05','publication_checkpoint':7,
          'copyright_holders':['Massimiliano Virgilio','Wannes Dermauw'],
          'code_license':'MIT','project_content_license':'CC-BY-4.0',
          'alignment_access':'distributed separately upon request; excluded from public licence grants','omitted_identical_export_copies':len(aliases),
          'omitted_bytes':sum(x['bytes'] for x in aliases),'analysis_result_files':dict(counts),
          'statistical_outputs_recomputed':False,
          'native_R_execution_during_assembly':'User supplied R 4.5.1 Windows execution on 2026-10-05; isolated implementation checks, reporting checks and formatter passed. Native exports independently verified.',
          'historical_run_records':'Preserved. Machine paths and early internal labels describe the original runs.',
          'authority_report':'Verified native R exports. The original checker stopped on Step 4 versus Step 04 descriptive labels; original logs and completed independent comparisons are retained in docs/native_R_verification/. Earlier reconstruction retained separately.',
          'manifest_scope':'All archive files except archive_manifest.tsv itself; omitted identical copies are mapped in duplicate_exports.tsv'}
    (GUIDES/'assembly_info.json').write_text(json.dumps(info,indent=2)+'\n')
    retained.append((GUIDES/'assembly_info.json','assembly_info.json'))
    rows=[{'archive_path':rel,'repository_source':p.relative_to(ROOT).as_posix(),
           'bytes':p.stat().st_size,'sha256':sha(p)} for p,rel in retained]
    table(GUIDES/'archive_manifest.tsv',rows,['archive_path','repository_source','bytes','sha256'])
    retained.append((GUIDES/'archive_manifest.tsv','archive_manifest.tsv'))
    args.output.parent.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(args.output,'w',zipfile.ZIP_DEFLATED,compresslevel=6) as z:
        for path,rel in retained:z.write(path,'PESTFLY_S6/'+rel)
    with zipfile.ZipFile(args.output) as z:
        bad=z.testzip()
        if bad:raise ValueError(f'ZIP integrity check failed for {bad}')
        for row in rows:
            content=z.read('PESTFLY_S6/'+row['archive_path'])
            if len(content)!=row['bytes'] or hashlib.sha256(content).hexdigest()!=row['sha256']:
                raise ValueError(f'Archived file mismatch: {row["archive_path"]}')
    print(json.dumps({'archive':str(args.output),'files':len(retained),
                      'size_MiB':round(args.output.stat().st_size/1024/1024,2),
                      'sha256':sha(args.output),'omitted_identical_export_copies':len(aliases),
                      'complete_retained_file_hash_check':True},indent=2))

if __name__=='__main__':main()
