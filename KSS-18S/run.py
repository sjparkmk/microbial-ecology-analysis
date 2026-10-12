#!/usr/bin/env python3
"""Check the frozen input package or explicitly execute ONE existing stage.

The default is read-only validation, never an analysis. Commands invoke the
frozen numerical definitions; they do not recreate the publication layout.
"""
from pathlib import Path
import argparse
import csv
import hashlib
import json
import os
import shutil
import subprocess
import sys

PACKAGE = Path(__file__).resolve().parent
STAGES = ('prepare', 'baseline', 'composition', 'community', 'daa', 'daa-retry',
          'aldex', 'rf-grouped', 'rf-lopo', 'breadth-depth', 'depth-statistics', 'a1')

def contained(root, relative):
    p=(root/relative).resolve()
    if root.resolve() not in p.parents: raise ValueError('Not a contained input: '+str(relative))
    return p

def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as f:
        for b in iter(lambda:f.read(1024*1024),b''): h.update(b)
    return h.hexdigest()

def records(path):
    with path.open(encoding='utf-8-sig',newline='') as f:
        return list(csv.DictReader(f,delimiter='\t'))

def count_schema(path, expected_samples):
    ids=set(); totals=[0]*len(expected_samples)
    with path.open(encoding='utf-8-sig',newline='') as f:
        reader=csv.reader(f,delimiter='\t'); head=next(reader)
        assert head[1:]==expected_samples, str(path)+': sample order'
        for row in reader:
            assert len(row)==len(head) and row[0] not in ids
            ids.add(row[0]); values=[int(x) for x in row[1:]]
            assert all(x>=0 for x in values)
            totals=[a+b for a,b in zip(totals,values)]
    assert ids and all(x>0 for x in totals)
    return ids,totals

def validate(root):
    manifest=json.loads((root/'input_fingerprints.json').read_text(encoding='utf-8'))
    assert len({r['path'] for r in manifest})==len(manifest)
    for r in manifest:
        p=contained(root,r['path'])
        assert sha(p)==r['sha256'], 'Input hash mismatch: '+r['path']
    summary={}
    for dataset,n in [('kss',15),('d2024',75),('eukbank',50)]:
        if dataset=='kss':
            md=records(root/'kss/inputs/metadata_analysis.tsv'); key='SampleID'
            count=root/'kss/inputs/otutab_raw_corrected.tsv'
            taxonomy=root/'kss/inputs/taxonomy_pr2v511.tsv'
        else:
            md=records(root/dataset/'subsets/metadata_samples.tsv')
            key='SampleID' if dataset=='d2024' else 'sample'
            count=root/dataset/'subsets/ALL_EUK/asv_table.tsv'
            taxonomy=root/dataset/'subsets/ALL_EUK/taxonomy.tsv'
        samples=[r[key] for r in md]
        assert len(samples)==n and len(set(samples))==n
        # External metadata may store a different row order: equality of the
        # sample set is required; count order is preserved, not changed.
        with count.open(encoding='utf-8-sig') as f: columns=f.readline().rstrip('\r\n').split('\t')[1:]
        assert set(columns)==set(samples)
        features,totals=count_schema(count,columns)
        with taxonomy.open(encoding='utf-8-sig') as f:
            reader=csv.reader(f,delimiter='\t'); next(reader)
            taxids=[r[0] for r in reader]
        assert len(taxids)==len(set(taxids)) and features<=set(taxids)
        summary[dataset]={'samples':n,'ASVs':len(features)}
        if dataset=='kss':
            rare_ids,rare_totals=count_schema(root/'kss/inputs/otutab_rare_20000.tsv',samples)
            assert rare_ids==features and rare_totals==[20000]*15
    md=records(root/'comparative/inputs/metadata.tsv')
    samples=[r['sample'] for r in md]
    assert len(samples)==140 and len(set(samples))==140
    assert {d:sum(r['dataset']==d for r in md) for d in ('ours','d2024','eukbank')}=={'ours':15,'d2024':75,'eukbank':50}
    assert len({r['project'] for r in md})==34
    for rank in ('Genus','Family'): count_schema(root/f'comparative/inputs/counts_{rank}.tsv',samples)
    print(json.dumps({'status':'PASS','files_hashed':len(manifest),'inputs':summary,
                      'merged_samples':140,'source_projects':34,'analysis_executed':False},indent=2))

def command(stage, inp, out, rscript):
    s=PACKAGE/'scripts'; comp=inp/'comparative/inputs'; py=sys.executable
    destinations={'prepare':'prepared','baseline':'baseline','composition':'composition','community':'community',
        'daa':'daa','daa-retry':'daa_ours_family_retry_v1','aldex':'aldex','rf-grouped':'rf_grouped',
        'rf-lopo':'rf_lopo','breadth-depth':'breadth_rf_depth','depth-statistics':'depth_statistics','a1':'a1_sensitivity'}
    dest=out/destinations[stage]
    if dest.exists(): raise FileExistsError('Refusing existing stage output: '+str(dest))
    if stage in ('community','daa','aldex','daa-retry','depth-statistics','a1'):
        if not (Path(rscript).is_file() or shutil.which(rscript)):
            raise FileNotFoundError('Rscript not found; provide --rscript PATH. Nothing is installed automatically.')
    if stage=='prepare': cmd=[py,str(s/'prepare_comparison.py'),'--root',str(inp),'--out',str(dest)]
    elif stage=='baseline': cmd=[py,str(s/'build_baseline_support.py'),'--out',str(dest)]
    elif stage in ('composition','rf-lopo','breadth-depth'):
        name={'composition':'composition_numeric.py','rf-lopo':'rf_ecological_context.py','breadth-depth':'breadth_rf_depth.py'}[stage]
        cmd=[py,str(s/name),'--run']
    elif stage in ('community','daa','aldex','daa-retry'):
        name='retry_ours_family_numeric.R' if stage=='daa-retry' else 'comparison_statistics.R'
        cmd=[rscript,str(s/name),'--input',str(comp),'--out',str(dest),'--mode','daa' if stage=='daa-retry' else stage]
    elif stage=='rf-grouped':
        cmd=[py,str(s/'comparison_python.py'),'--input',str(comp),'--out',str(dest),'--mode','rf']
    elif stage=='a1': cmd=[rscript,str(s/'a1_restricted_sensitivity_v1.R'),'--run']
    else:
        # This R script writes next to distance files. Copy only its six inputs
        # into a new output stage; never pass the read-only input directory.
        source=out/'breadth_rf_depth/tables'
        needed=[source/f'Depth_distance_{rank}_{rep}.tsv' for rank in ('Genus','Family') for rep in ('raw_relative','mean_14000','mean_10000')]
        if not all(p.is_file() for p in needed): raise FileNotFoundError('Run breadth-depth first in the same output root.')
        dest.mkdir(parents=True); (dest/'tables').mkdir()
        for p in needed: shutil.copy2(p,dest/'tables'/p.name)
        cmd=[rscript,str(s/'depth_statistics_v1.R'),str(dest/'tables'),str(comp)]
    return cmd

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--input-root',type=Path,default=Path('/input/ur/data'))
    p.add_argument('--output-root',type=Path,default=Path('/output/ur/data'))
    mode=p.add_mutually_exclusive_group(); mode.add_argument('--check',action='store_true'); mode.add_argument('--run',choices=STAGES)
    p.add_argument('--rscript',default='Rscript')
    a=p.parse_args(); inp=a.input_root.resolve(); out=a.output_root.resolve()
    if inp==out or inp in out.parents or out in inp.parents:
        p.error('Input and output roots must be separate, non-nested directories.')
    validate(inp)
    if not a.run: return
    env=os.environ.copy(); env.update(KSS_INPUT_ROOT=str(inp),KSS_OUTPUT_ROOT=str(out),KSS_SCRIPT_ROOT=str(PACKAGE/'scripts'))
    cmd=command(a.run,inp,out,a.rscript)
    out.mkdir(parents=True,exist_ok=True)
    print('Explicit analysis execution requested:',a.run,flush=True)
    subprocess.run(cmd,env=env,check=True)

if __name__=='__main__': main()
