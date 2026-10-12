"""Prepare provenance-checked Run A + selected external count tables (no inference)."""
from __future__ import annotations
import argparse
import hashlib
import json
from pathlib import Path
import numpy as np
import pandas as pd

VERSION = '1.0.0'
RANKS = ['Domain','Supergroup','Division','Subdivision','Class','Order','Family','Genus']

def read(path, **kw):
    return pd.read_csv(path, sep='\t', keep_default_na=False, **kw)

def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def save(df, path, index=False):
    path.parent.mkdir(parents=True, exist_ok=True)
    data = df.to_csv(sep='\t',index=index,lineterminator='\n').encode('utf-8')
    if path.exists():
        if path.read_bytes() != data:
            raise FileExistsError(f'Refusing changed output: {path}; use a new output version')
    else:
        path.write_bytes(data)

def load_source(root, name):
    base = root/'kss'
    run = root/name/'subsets'
    if name == 'ours':
        paths = [base/'inputs/otutab_raw_corrected.tsv',base/'inputs/taxonomy_pr2v511.tsv',base/'inputs/metadata_analysis.tsv']
        count,tax,md = read(paths[0],index_col=0),read(paths[1],index_col=0),read(paths[2])
        md=md.loc[md.IncludeInGroupTests.eq('yes')].copy().set_index('SampleID')
        ids=md.index
        meta=pd.DataFrame({'source_sample':ids,'dataset':name,
            'sample_type':np.where(md.sample_type.eq('soil'),'land_soil','land_sediment'),
            'substrate':md.sample_type,'exposure':md.AnalysisGroup,'site':md.site_label,
            'region':md.RegionGroup,'season':'unavailable','cluster':name+':'+md.index,
            'project':name,'latitude':np.nan,'longitude':np.nan,'pH':pd.to_numeric(md.pH_mid)},index=ids)
    else:
        paths=[run/'ALL_EUK/asv_table.tsv',run/'ALL_EUK/taxonomy.tsv',run/'metadata_samples.tsv']
        count,tax,md=read(paths[0],index_col=0),read(paths[1],index_col=0),read(paths[2])
        key='SampleID' if name=='d2024' else 'sample'
        assert not md[key].duplicated().any()
        md=md.set_index(key).loc[count.columns]
        ids=md.index
        if name=='d2024':
            meta=pd.DataFrame({'source_sample':ids,'dataset':name,
                'sample_type':md.type.map({'water':'marine_water','soil':'marine_sediment'}),
                'substrate':md.type.map({'water':'water','soil':'sediment'}),'exposure':'unavailable',
                'site':md.site,'region':md.Group,'season':md.time,'cluster':name+':'+md.site,
                'project':name,'latitude':np.nan,'longitude':np.nan,'pH':pd.to_numeric(md.pH,errors='coerce')},index=ids)
        else:
            meta=pd.DataFrame({'source_sample':ids,'dataset':name,'sample_type':md.envplot,
                'substrate':md.envplot.str.split('_').str[-1],'exposure':'unavailable',
                'site':md.biosample,'region':'unavailable','season':'unavailable',
                'cluster':name+':'+md.project,'project':md.project,
                'latitude':pd.to_numeric(md.latitude,errors='coerce'),
                'longitude':pd.to_numeric(md.longitude,errors='coerce'),'pH':np.nan},index=ids)
    assert count.index.is_unique and tax.index.is_unique and count.columns.is_unique
    count=count.loc[:,ids].apply(pd.to_numeric)
    assert np.isfinite(count.to_numpy()).all() and (count.to_numpy()>=0).all()
    assert (count.to_numpy()==np.floor(count.to_numpy())).all()
    assert (count.sum(axis=0)>0).all()
    assert set(count.index)<=set(tax.index)
    assert not meta.sample_type.isna().any()
    count.columns = name+':'+pd.Index(ids)
    meta.index=count.columns
    meta.index.name='sample'
    meta['raw_eukaryote_reads']=count.sum()
    return count.astype('int64'),tax.loc[count.index],meta,paths

def aggregate(count,tax,rank):
    boot=rank+('_bootstrap' if rank+'_bootstrap' in tax else '_boot')
    labels=tax[rank].astype(str)
    assigned=~labels.str.strip().isin(['','NA','nan','Unassigned_'+rank])
    b=pd.to_numeric(tax[boot],errors='coerce')
    if (assigned & (b<50)).any():
        raise ValueError('Named rank below bootstrap threshold')
    assigned &= b.ge(50)
    prefix=RANKS[:RANKS.index(rank)+1]
    lineage=tax[prefix].astype(str).replace('', '?').agg('|'.join,axis=1)
    named=count.loc[assigned].groupby(lineage.loc[assigned],sort=True).sum()
    named.index.name='feature'
    unresolved=count.loc[~assigned].sum()
    assert np.array_equal(named.sum().to_numpy()+unresolved.to_numpy(),count.sum().to_numpy())
    manifest=pd.DataFrame({'feature':named.index,'rank':rank,'label':named.index.str.split('|').str[-1]})
    return named,unresolved,manifest

def run(root,out,check_only=False):
    loaded={n:load_source(root,n) for n in ['ours','d2024','eukbank']}
    meta=pd.concat([v[2] for v in loaded.values()])
    assert meta.index.is_unique and meta.groupby('dataset').size().to_dict()=={'d2024':75,'eukbank':50,'ours':15}
    matrices={}; manifests={}; assignment=[]; audits=[]
    for rank in ['Genus','Family']:
        portions=[]; maps=[]
        for name,(count,tax,md,paths) in loaded.items():
            tab,unknown,mapping=aggregate(count,tax,rank)
            assert (tab.sum()>0).all()
            portions.append(tab);maps.append(mapping)
            assignment.append(pd.DataFrame({'sample':md.index,'dataset':name,'rank':rank,
                'total_reads':count.sum().to_numpy(),'assigned_reads':tab.sum().to_numpy(),
                'unassigned_reads':unknown.to_numpy(),'assigned_fraction':(tab.sum()/count.sum()).to_numpy()}))
        matrices[rank]=pd.concat(portions,axis=1).fillna(0).astype('int64').sort_index().loc[:,meta.index]
        manifests[rank]=pd.concat(maps).drop_duplicates('feature')
        assert manifests[rank].feature.nunique()==len(matrices[rank])
    for name,(count,tax,md,paths) in loaded.items():
        audits.append({'dataset':name,'samples':count.shape[1],'asvs':count.shape[0],
            'reads':int(count.to_numpy().sum()),'min_reads':int(count.sum().min()),
            'zero_fraction':float((count.to_numpy()==0).mean())})
    hashes=[{'dataset':n,'path':str(p),'sha256':digest(p),'bytes':p.stat().st_size}
            for n,v in loaded.items() for p in v[3]]
    print(json.dumps({'audit':audits,'ranks':{r:len(t) for r,t in matrices.items()},'check_only':check_only},indent=2))
    if check_only:return
    for rank,tab in matrices.items():
        save(tab,out/f'counts_{rank}.tsv',True)
        save(tab.div(tab.sum(axis=0),axis=1),out/f'relabund_{rank}.tsv',True)
        save(manifests[rank],out/f'features_{rank}.tsv')
    save(meta,out/'metadata.tsv',True)
    save(pd.concat(assignment),out/'assignment_coverage.tsv')
    save(pd.DataFrame(audits),out/'input_audit.tsv')
    save(pd.DataFrame(hashes),out/'source_manifest.tsv')
    save(pd.crosstab(meta.dataset,meta.sample_type),out/'dataset_by_habitat.tsv',True)
    ours=loaded['ours']
    save(ours[0].loc[ours[0].sum(axis=1)>0],out/'ours_ASV_counts.tsv',True)
    save(ours[2],out/'ours_metadata.tsv',True)
    save(ours[1],out/'ours_taxonomy.tsv',True)

def self_test():
    counts=pd.DataFrame({'s1':[1,2,3],'s2':[0,4,2]},index=['a','b','c'])
    tax=pd.DataFrame({r:['D','D','D'] for r in RANKS},index=counts.index)
    tax['Genus']=['same','same',''];tax['Family']=['f1','f2','f3'];tax['Genus_boot']=[50,99,1]
    named,unknown,man=aggregate(counts,tax,'Genus')
    assert len(named)==2 and unknown.tolist()==[3,2]
    assert man.label.tolist()==['same','same']
    print('PASS: lineage homonyms separate; unresolved excluded; counts conserved')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--version',action='version',version=VERSION)
    p.add_argument('--root',type=Path)
    p.add_argument('--out',type=Path)
    p.add_argument('--check',action='store_true')
    p.add_argument('--self-test',action='store_true')
    a=p.parse_args()
    if a.self_test:self_test()
    elif a.root and a.out:run(a.root.resolve(),a.out.resolve(),a.check)
    else:p.error('--root and --out required unless --self-test')
