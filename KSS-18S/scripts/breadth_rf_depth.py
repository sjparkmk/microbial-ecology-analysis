"""Frozen-scope v4 additions; no original input/output modifications or DAA fits."""
from pathlib import Path
import importlib.util, hashlib, json, platform, datetime, argparse
import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import spearmanr
from sklearn.ensemble import RandomForestClassifier
from sklearn.dummy import DummyClassifier
from sklearn.metrics import balanced_accuracy_score, accuracy_score, confusion_matrix
import sklearn, scipy

from paths import INPUT_ROOT, OUTPUT_ROOT, require_run
SCRIPT_DIR=Path(__file__).resolve().parent
OUT=OUTPUT_ROOT/'breadth_rf_depth'
FROOT=INPUT_ROOT
COMP=INPUT_ROOT/'comparative/inputs'
ALIASES={'ours':'KSS-18S','d2024':'Coastal reference','eukbank':'EukBank subset'}
SEED=20261001
RANKS=['Domain','Supergroup','Division','Subdivision','Class','Order','Family','Genus']

def sha(p):return hashlib.sha256(p.read_bytes()).hexdigest()
def save(df,name,index=False):
    p=OUT/'tables'/name
    p.parent.mkdir(parents=True,exist_ok=True)
    assert not p.exists(), f'Immutable run output already exists: {p}'
    df.to_csv(p,sep='\t',index=index,lineterminator='\n',float_format='%.12g')

def mapping(tax,rank):
    boot=rank+('_bootstrap' if rank+'_bootstrap' in tax else '_boot')
    assigned=~tax[rank].astype(str).str.strip().isin(['','NA','nan','Unassigned_'+rank])
    assigned &= pd.to_numeric(tax[boot],errors='coerce').ge(50)
    lin=tax[RANKS[:RANKS.index(rank)+1]].astype(str).replace('','?').agg('|'.join,axis=1)
    return lin.where(assigned,'')

def relative_class(count,tax):
    lab=mapping(tax,'Class'); tab=count.loc[lab.ne('')].groupby(lab[lab.ne('')]).sum()
    return tab.div(count.sum(axis=0),axis=1)

def breadth(rel,meta,dataset,stratum):
    rows=[]
    groups=sorted(meta[stratum].astype(str).unique())
    assert len(groups)>=2
    means=pd.DataFrame({g:rel.loc[:,meta.index[meta[stratum].astype(str).eq(g)]].mean(axis=1) for g in groups})
    for feature,r in means.iterrows():
        freq=int((rel.loc[feature]>0).sum())
        if freq<2 or r.sum()<=0:continue
        p=r.to_numpy()/r.sum();b=1/np.sum(p*p);ba=(b-1)/(len(groups)-1)
        assert -1e-12<=ba<=1+1e-12
        rows.append({'dataset_id':dataset,'dataset':ALIASES[dataset],'stratum':stratum,'feature':feature,'Class':feature.split('|')[-1],
                     'k_strata':len(groups),'libraries':len(meta),'prevalence':freq,'mean_read_fraction':float(rel.loc[feature].mean()),
                     'Levins_B':float(b),'Levins_BA':float(np.clip(ba,0,1)),
                     'stratum_mean_read_fractions':json.dumps({g:float(r[g]) for g in groups},sort_keys=True),
                     'status':'descriptive_only_no_null_or_specialist_label'})
    return rows,means

def main(draws):
    import matplotlib
    for row in json.loads((INPUT_ROOT/'input_fingerprints.json').read_text(encoding='utf-8')):
        assert sha(INPUT_ROOT/row['path'])==row['sha256'], f"Changed input: {row['path']}"
    assert draws==100,'Prospective specification fixes 100 draws'
    spec=importlib.util.spec_from_file_location('source_prepare',SCRIPT_DIR/'prepare_comparison.py')
    prep=importlib.util.module_from_spec(spec);spec.loader.exec_module(prep)
    loaded={d:prep.load_source(FROOT,d) for d in ALIASES}
    md=pd.read_csv(COMP/'metadata.tsv',sep='\t',index_col=0,keep_default_na=False)
    assert len(md)==140 and md.index.is_unique
    assert all(np.array_equal(c.columns,m.index) for c,t,m,paths in loaded.values())
    for d,(c,t,m,paths) in loaded.items():
        assert m.index.tolist()==md.index[md.dataset.eq(d)].tolist()
        assert np.array_equal(c.sum().to_numpy(),pd.to_numeric(md.loc[m.index,'raw_eukaryote_reads']).to_numpy())
    pd.DataFrame({'library':md.index,'dataset':[ALIASES[d] for d in md.dataset],'sample_type':md.sample_type,'raw_eukaryote_reads':md.raw_eukaryote_reads}).to_csv(OUT/'input_library_register.tsv',sep='\t',index=False)
    # Prospectively recorded classification; never replace canonical metadata.
    oursmd=pd.read_csv(INPUT_ROOT/'kss/metadata_analysis_current.tsv',sep='\t',keep_default_na=False).set_index('SampleID')
    oursmd=oursmd.loc[loaded['ours'][2].source_sample]
    oursmd.index=loaded['ours'][2].index
    rare=pd.read_csv(INPUT_ROOT/'kss/inputs/otutab_rare_20000.tsv',sep='\t',index_col=0)
    rare=rare.loc[:,loaded['ours'][2].source_sample];rare.columns=oursmd.index
    assert np.array_equal(rare.sum(),np.repeat(20000,15))
    primary_tax=loaded['ours'][1].reindex(rare.index)
    assert primary_tax['Domain'].eq('Eukaryota').all()
    classrels={d:relative_class(c,t) for d,(c,t,m,paths) in loaded.items()}
    classrels['ours']=relative_class(rare,primary_tax)
    breadth_rows=[];coverage=[];means_record=[];phrows=[];phprofiles=[]
    for d,rel in classrels.items():
        m=loaded[d][2].copy()
        if d=='ours':
            for col in ['sample_type','pH_class','AnalysisGroup','RegionGroup']:
                m[col]=oursmd[col]
            strata=['sample_type','pH_class','AnalysisGroup','RegionGroup']
        else:strata=['sample_type']
        save(rel,f'Class_relative_all_Eukaryota_{d}.tsv',True)
        coverage.extend([{'sample':s,'dataset':ALIASES[d],'Class_assigned_fraction':float(rel[s].sum())} for s in m.index])
        for stratum in strata:
            rows,means=breadth(rel,m,d,stratum);breadth_rows.extend(rows)
            for feature,r in means.iterrows():
                for group,value in r.items():means_record.append({'dataset':ALIASES[d],'stratum':stratum,'feature':feature,'Class':feature.split('|')[-1],'group':group,'mean_read_fraction':float(value),'n_libraries':int(m[stratum].astype(str).eq(group).sum())})
        if d=='eukbank':continue
        if d=='ours':
            units=rel.T.copy();pH=pd.to_numeric(oursmd.pH_mid)
            for sample in m.index:
                phprofiles.append({'dataset':ALIASES[d],'unit':sample,'pH':float(pH[sample]),'record_type':'library-associated existing point/range midpoint','n_libraries':1,'site':m.loc[sample,'site'],'pH_raw':oursmd.loc[sample,'pH_raw']})
        else:
            units=rel.T.groupby(m.site).mean()
            pH=pd.to_numeric(m.pH).groupby(m.site).first()
            assert m.groupby('site').pH.nunique().eq(1).all() and len(units)==7
            for site in units.index:
                phprofiles.append({'dataset':ALIASES[d],'unit':site,'pH':float(pH[site]),'record_type':'site-level record reused across libraries','n_libraries':int(m.site.eq(site).sum()),'site':site,'pH_raw':str(pH[site])})
        for feature in units.columns:
            w=units[feature].to_numpy()
            if (w>0).sum()<2 or w.sum()==0:continue
            p=w/w.sum();x=pH.loc[units.index].to_numpy(dtype=float);mu=float(p@x);sd=float(np.sqrt(p@((x-mu)**2)))
            phrows.append({'dataset':ALIASES[d],'feature':feature,'Class':feature.split('|')[-1],'independent_record_units':len(units),'occupied_units':int((w>0).sum()),'weighted_pH_mean':mu,'weighted_pH_SD':sd,'min_recorded_pH':float(x.min()),'max_recorded_pH':float(x.max()),'status':'descriptive_sampled_pH_distribution_not_tolerance'})
    save(pd.DataFrame(breadth_rows),'R9_breadth_all_eligible_Classes.tsv')
    save(pd.DataFrame(coverage),'R9_Class_assignment_coverage.tsv')
    save(pd.DataFrame(means_record),'R9_stratum_means_and_sizes.tsv')
    save(pd.DataFrame(phrows),'R9_pH_weighted_profile.tsv')
    save(pd.DataFrame(phprofiles),'R9_pH_independent_record_register.tsv')
    common=set.intersection(*(set(r.index) for r in classrels.values()))
    equal_mean=pd.DataFrame({ALIASES[d]:r.mean(axis=1) for d,r in classrels.items()}).fillna(0)
    common_rank=equal_mean.loc[sorted(common)].mean(axis=1).sort_values(ascending=False)
    save(common_rank.head(20).rename('equal_dataset_mean_read_fraction').reset_index().rename(columns={'index':'feature'}),'R9_plot_taxon_selection.tsv')
    eligibility=[]
    for d,(c,t,m,paths) in loaded.items():
        for target in ['exposure','geography','pH']:
            col={'exposure':'exposure','geography':'region','pH':'pH'}[target]
            valid=m[col].astype(str).str.strip().ne('unavailable') & m[col].notna() & m[col].astype(str).ne('')
            if d=='ours':status='available_within_one_project_only'
            elif d=='d2024' and target=='pH':status='descriptive_7_site_records_not_75_independent_pH'
            elif d=='d2024' and target=='geography':status='source_specific_coastal_label_not_shared_geographic_target'
            else:status='not_estimable_missing_harmonized_target'
            eligibility.append({'dataset':ALIASES[d],'target':target,'libraries':len(m),'nonmissing_local_labels':int(valid.sum()),'n_projects':m.project.nunique(),'status':status})
    save(pd.DataFrame(eligibility),'R10_target_eligibility.tsv')
    # New within-project exploratory prediction, no favorable-rank or target selection.
    scores=[];predictions=[];folds=[];confusions=[]
    for rank in ['Genus','Family']:
        X=pd.read_csv(COMP/f'counts_{rank}.tsv',sep='\t',index_col=0).loc[:,oursmd.index].T
        X=X.div(X.sum(axis=1),axis=0)
        for target,col in [('exposure','AnalysisGroup'),('geography','RegionGroup')]:
            y=oursmd[col].astype(str).to_numpy();classes=np.unique(y)
            assert len(classes)==2
            yp=[];yd=[]
            for i in range(len(X)):
                train=np.delete(np.arange(len(X)),i);test=np.array([i])
                assert not set(train)&set(test) and set(y[train])==set(classes)
                keep=(X.iloc[train].gt(0).mean(axis=0)>=.05)&(X.iloc[train].var(axis=0)>0)
                assert keep.any()
                model=RandomForestClassifier(n_estimators=100,max_features='sqrt',class_weight='balanced',random_state=SEED,n_jobs=1)
                model.fit(X.iloc[train].loc[:,keep],y[train]);p=model.predict(X.iloc[test].loc[:,keep])[0]
                dummy=DummyClassifier(strategy='most_frequent');dummy.fit(X.iloc[train],y[train]);q=dummy.predict(X.iloc[test])[0]
                yp.append(p);yd.append(q)
                predictions.append({'rank':rank,'target':target,'sample':X.index[i],'truth':y[i],'RF_prediction':p,'dummy_prediction':q,'fold':i,'n_train':len(train),'n_features':int(keep.sum())})
                folds.append({'rank':rank,'target':target,'fold':i,'test_sample':X.index[i],'training_samples':json.dumps(X.index[train].tolist()),'class_counts_train':json.dumps(pd.Series(y[train]).value_counts().to_dict()),'no_overlap':True,'training_only_feature_filter':True})
            for model,pr in [('RF',yp),('dummy',yd)]:
                scores.append({'dataset':'KSS-18S','rank':rank,'target':target,'model':model,'n':len(y),'accuracy':accuracy_score(y,pr),'balanced_accuracy':balanced_accuracy_score(y,pr),'class_sizes':json.dumps(pd.Series(y).value_counts().to_dict()),'validation':'leave_one_library_out_within_one_project','status':'exploratory_not_external_transfer'})
                cm=confusion_matrix(y,pr,labels=classes)
                for a,truth in enumerate(classes):
                    for b,pred in enumerate(classes):confusions.append({'rank':rank,'target':target,'model':model,'true_class':truth,'predicted_class':pred,'n':int(cm[a,b]),'true_class_n':int(cm[a].sum())})
            print(f'RF {rank} {target} complete',flush=True)
    save(pd.DataFrame(scores),'R10_RF_summary.tsv');save(pd.DataFrame(predictions),'R10_RF_predictions.tsv')
    save(pd.DataFrame(folds),'R10_RF_fold_checks.tsv');save(pd.DataFrame(confusions),'R10_RF_confusion.tsv')
    # Matched per-ASV hypergeometric draws, rank aggregation after subsampling.
    combined={};codes={};positive={}
    for rank in ['Genus','Family']:
        original=pd.read_csv(COMP/f'counts_{rank}.tsv',sep='\t',index_col=0)
        rel=original.div(original.sum(axis=0),axis=1)
        combined[rank]=original
        save(pd.DataFrame(squareform(pdist(rel.T,'braycurtis')),index=md.index,columns=md.index),f'Depth_distance_{rank}_raw_relative.tsv',True)
        for d,(c,t,m,paths) in loaded.items():
            labels=mapping(t,rank);codes[(d,rank)]=pd.Index(original.index).get_indexer(labels)
    for d,(c,t,m,paths) in loaded.items():
        a=c.to_numpy(dtype=np.int64)
        positive[d]=[(np.flatnonzero(a[:,j]),a[a[:,j]>0,j]) for j in range(a.shape[1])]
    depth_summary=[];draw_diagnostics=[];read_assignment=[]
    for depth in [14000,10000]:
        meanD={r:np.zeros((140,140)) for r in combined};draw_D_var={r:np.zeros((140,140)) for r in combined}
        for iteration in range(draws):
            rng=np.random.default_rng(np.random.SeedSequence([SEED,depth,iteration]))
            agg={r:np.zeros((len(c),140),dtype=np.int64) for r,c in combined.items()}
            offset=0
            for d,(c,t,m,paths) in loaded.items():
                for j,(ix,colors) in enumerate(positive[d]):
                    draw=rng.multivariate_hypergeometric(colors,depth,method='marginals')
                    assert draw.sum()==depth and (draw<=colors).all() and (draw>=0).all()
                    for rank in combined:
                        code=codes[(d,rank)][ix];good=code>=0
                        agg[rank][:,offset+j]=np.bincount(code[good],weights=draw[good],minlength=len(combined[rank])).astype(np.int64)
                        assigned=int(draw[good].sum());assert assigned>0
                        read_assignment.append({'depth':depth,'iteration':iteration,'rank':rank,'sample':m.index[j],'dataset':ALIASES[d],'assigned_reads':assigned,'unassigned_reads':depth-assigned,'assigned_fraction':assigned/depth})
                offset+=len(m)
            assert offset==140
            for rank,a in agg.items():
                rel=a/a.sum(axis=0)
                D=squareform(pdist(rel.T,'braycurtis'))
                assert np.isfinite(D).all() and np.allclose(D,D.T) and np.allclose(np.diag(D),0)
                meanD[rank]+=D/draws;draw_D_var[rank]+=D**2/draws
                raw=pd.read_csv(OUT/'tables'/f'Depth_distance_{rank}_raw_relative.tsv',sep='\t',index_col=0).to_numpy()
                upper=np.triu_indices(140,1)
                draw_diagnostics.append({'depth':depth,'iteration':iteration,'rank':rank,'distance_spearman_vs_raw':spearmanr(D[upper],raw[upper]).statistic,'mean_abs_distance_delta':np.mean(np.abs(D[upper]-raw[upper]))})
            if iteration%20==0 or iteration==draws-1:print(f'Depth {depth}: draw {iteration+1}/{draws}',flush=True)
        for rank,D in meanD.items():
            save(pd.DataFrame(D,index=md.index,columns=md.index),f'Depth_distance_{rank}_mean_{depth}.tsv',True)
            raw=pd.read_csv(OUT/'tables'/f'Depth_distance_{rank}_raw_relative.tsv',sep='\t',index_col=0).to_numpy()
            upper=np.triu_indices(140,1);sd=np.sqrt(np.maximum(draw_D_var[rank]-D**2,0))
            depth_summary.append({'rank':rank,'depth':depth,'draws':draws,'n_libraries':140,'minimum_input_reads':int(pd.to_numeric(md.raw_eukaryote_reads).min()),'mean_distance_spearman_vs_raw':spearmanr(D[upper],raw[upper]).statistic,'mean_abs_distance_delta':float(np.mean(np.abs(D[upper]-raw[upper]))),'max_abs_distance_delta':float(np.max(np.abs(D[upper]-raw[upper]))),'mean_draw_distance_SD':float(sd[upper].mean()),'status':'sensitivity_not_primary_result_replacement'})
    save(pd.DataFrame(depth_summary),'Depth_geometry_summary.tsv')
    save(pd.DataFrame(draw_diagnostics),'Depth_all_draw_geometry.tsv')
    save(pd.DataFrame(read_assignment),'Depth_all_draw_assignment.tsv')
    environment={'utc':datetime.datetime.now(datetime.timezone.utc).isoformat(),'python':platform.python_version(),'numpy':np.__version__,'pandas':pd.__version__,'scipy':scipy.__version__,'sklearn':sklearn.__version__,'matplotlib':matplotlib.__version__,'seed':SEED,'draws':draws,'depths':[14000,10000],'spec_sha256':sha(SCRIPT_DIR.parent/'provenance/ANALYSIS_CONTRACT.md'),'script_sha256':sha(Path(__file__)),'unchanged_originals':True}
    (OUT/'environment_python_v1.json').write_text(json.dumps(environment,indent=2),encoding='utf-8')
    print('ALL PYTHON ADDITIONS COMPLETE',flush=True)

if __name__=='__main__':
    require_run()
    OUT.mkdir(parents=True,exist_ok=False)
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--draws',type=int,default=100);args=p.parse_args();main(args.draws)
