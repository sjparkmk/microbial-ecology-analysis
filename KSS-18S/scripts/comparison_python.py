"""Grouped RF, compositional ordination and pre-specified network diagnostics."""
from __future__ import annotations
import argparse, itertools, json
from pathlib import Path
import numpy as np
import pandas as pd
import networkx as nx
from scipy import stats
from scipy.special import expit
from sklearn.decomposition import PCA
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import adjusted_rand_score,confusion_matrix,accuracy_score,balanced_accuracy_score
from sklearn.model_selection import GroupKFold
from prepare_comparison import save

VERSION='1.0.0';SEED=20260923

def read(path):return pd.read_csv(path,sep='\t',index_col=0)
def clr(x):
    y=np.log(np.asarray(x,float)+1.)
    return y-y.mean(axis=1,keepdims=True)
def ba(y,p,labels):
    return float(np.mean([np.mean(p[y==k]==k) for k in labels if np.any(y==k)]))
def shuffle_within(y,groups,rng):
    z=y.copy()
    for g in np.unique(groups):
        ix=np.flatnonzero(groups==g);z[ix]=rng.permutation(z[ix])
    return z

def rf_run(inp,out,nperm):
    md=read(inp/'metadata.tsv').fillna('unavailable')
    results=[];predictions=[];imps=[];nulls=[]
    groups=md.cluster.to_numpy()
    folds=list(GroupKFold(n_splits=5).split(md,groups=groups))
    save(pd.DataFrame([{'sample':md.index[i],'fold':fold} for fold,(_,te) in enumerate(folds) for i in te]),out/'cv_folds.tsv')
    for rank in ['Genus','Family']:
        raw=read(inp/f'counts_{rank}.tsv').T.loc[md.index]
        rel=raw.div(raw.sum(axis=1),axis=0)
        pc=PCA(n_components=2,svd_solver='full')
        ll=np.log(rel.to_numpy()+1e-6);ll-=ll.mean(axis=1,keepdims=True)
        coords=pc.fit_transform(ll)
        save(pd.DataFrame({'sample':md.index,'PC1':coords[:,0],'PC2':coords[:,1]}),out/f'{rank}_clr_pca.tsv')
        save(pd.DataFrame({'axis':[1,2],'variance_fraction':pc.explained_variance_ratio_,'relative_pseudocount':1e-6}),out/f'{rank}_pca_variance.tsv')
        prepared=[]
        for tr,te in folds:
            keep=(np.mean(raw.iloc[tr].to_numpy()>0,axis=0)>=.05)&(np.var(rel.iloc[tr].to_numpy(),axis=0)>0)
            prepared.append((tr,te,rel.to_numpy()[tr][:,keep],rel.to_numpy()[te][:,keep],np.flatnonzero(keep)))
        for target in ['sample_type','dataset']:
            y=md[target].to_numpy(dtype=str);labels=np.unique(y)
            def predict(yy,importance=False):
                pp=np.empty_like(yy);imp=np.zeros(raw.shape[1])
                for k,(tr,te,xt,xv,ii) in enumerate(prepared):
                    rf=RandomForestClassifier(n_estimators=100,class_weight='balanced',max_features='sqrt',random_state=SEED+k,n_jobs=1)
                    rf.fit(xt,yy[tr]);pp[te]=rf.predict(xv)
                    if importance:imp[ii]+=rf.feature_importances_/len(folds)
                return pp,imp
            p,importance=predict(y,True)
            observed=ba(y,p,labels)
            rng=np.random.default_rng(SEED)
            samples=[]
            for b in range(nperm):
                yp=shuffle_within(y,md.dataset.to_numpy(),rng) if target=='sample_type' else rng.permutation(y)
                pp,_=predict(yp);samples.append(ba(yp,pp,labels))
                if (b+1)%25==0:print('RF',rank,target,'permutation',b+1,flush=True)
            perm_p=(1+np.count_nonzero(np.asarray(samples)>=observed-1e-12))/(nperm+1)
            # OOB is a secondary random-library estimate and is deliberately separated.
            rr=RandomForestClassifier(n_estimators=500,class_weight='balanced',oob_score=True,random_state=SEED,n_jobs=1)
            rr.fit(rel,y)
            rng=np.random.default_rng(SEED);ug=np.unique(groups);accboot=[]
            for _ in range(1999):
                selected=rng.choice(ug,len(ug),replace=True)
                ix=np.concatenate([np.flatnonzero(groups==g) for g in selected])
                accboot.append(np.mean(y[ix]==p[ix]))
            ci=np.quantile(accboot,[.025,.975])
            results.append({'rank':rank,'target':target,'n':len(y),'clusters':len(ug),'cv_accuracy':accuracy_score(y,p),
                'cv_balanced_accuracy':observed,'accuracy_cluster_bootstrap_lower':ci[0],'accuracy_cluster_bootstrap_upper':ci[1],
                'oob_accuracy_secondary':rr.oob_score_,'permutation_p':perm_p,'permutations':nperm,
                'null':'labels permuted within dataset' if target=='sample_type' else 'labels permuted across libraries; diagnostic only',
                'validation':'group-blocked internal cross-validation; not an independent cohort'})
            predictions.append(pd.DataFrame({'sample':md.index,'rank':rank,'target':target,'truth':y,'predicted':p,'cluster':groups}))
            imps.append(pd.DataFrame({'rank':rank,'target':target,'feature':raw.columns,'mean_training_importance':importance}))
            nulls.append(pd.DataFrame({'rank':rank,'target':target,'draw':np.arange(nperm),'balanced_accuracy':samples}))
            cm=confusion_matrix(y,p,labels=labels)
            save(pd.DataFrame(cm,index=labels,columns=labels),out/f'{rank}_{target}_confusion.tsv',True)
            print('RF COMPLETE',rank,target,observed,perm_p,flush=True)
    save(pd.DataFrame(results),out/'rf_summary.tsv')
    save(pd.concat(predictions),out/'rf_oof_predictions.tsv')
    save(pd.concat(imps),out/'rf_training_importance.tsv')
    save(pd.concat(nulls),out/'rf_permutation_null.tsv')

def partition(g,seed):
    if g.number_of_edges()==0:return [{v} for v in g]
    return nx.community.louvain_communities(g,weight='weight',seed=seed)
def labels_from(communities,nodes):
    mm={v:i for i,c in enumerate(communities) for v in c}
    return np.asarray([mm[v] for v in nodes])
def chi_stat(y,env):
    yy=np.unique(y);ee=np.unique(env)
    tab=np.asarray([[np.count_nonzero((y==a)&(env==b)) for b in ee] for a in yy])
    if min(tab.shape)<2:return 0.
    return float(stats.chi2_contingency(tab,correction=False)[0])
def majority_cv(y,env,folds):
    pred=np.empty_like(y)
    for tr,te in folds:
        classes,n=np.unique(y[tr],return_counts=True);default=classes[n.argmax()]
        for level in np.unique(env[te]):
            jj=tr[env[tr]==level]
            if len(jj):v,n=np.unique(y[jj],return_counts=True);value=v[n.argmax()]
            else:value=default
            pred[te[env[te]==level]]=value
    return pred

def network_run(inp,out):
    md=read(inp/'metadata.tsv').fillna('unavailable');raw=read(inp/'counts_Genus.tsv').T.loc[md.index]
    ek=md.dataset.eq('eukbank').to_numpy();total=md.raw_eukaryote_reads.to_numpy()
    occurrence=np.mean(raw.to_numpy()[ek]/total[ek,None]>.002,axis=0)
    keep=occurrence>=.10
    x=raw.to_numpy()[:,keep];nodes=raw.columns[keep].tolist();z=clr(x)
    n=len(nodes)
    if n<5:raise ValueError('Fewer than 5 retained network genera')
    c=np.corrcoef(z[ek],rowvar=False);ii,jj=np.triu_indices(n,1)
    rr=c[ii,jj];tv=rr*np.sqrt(48/np.maximum(1-rr**2,1e-15));pv=2*stats.t.sf(abs(tv),48)
    from statsmodels.stats.multitest import multipletests
    q=multipletests(pv,method='fdr_bh')[1]
    edges=pd.DataFrame({'from':np.asarray(nodes)[ii],'to':np.asarray(nodes)[jj],'r':rr,'p':pv,'q':q})
    save(edges,out/'all_tested_edges.tsv')
    graphs={};parts={};stability=[]
    for qt in [.01,.05,.10]:
        g=nx.Graph();g.add_nodes_from(nodes)
        for a,b,w in edges.loc[(edges.r>0)&(edges.q<=qt),['from','to','r']].itertuples(index=False,name=None):g.add_edge(a,b,weight=w)
        graphs[qt]=g;parts[qt]=partition(g,SEED)
    for a,b in itertools.combinations(graphs,2):
        stability.append({'q_a':a,'q_b':b,'ARI':adjusted_rand_score(labels_from(parts[a],nodes),labels_from(parts[b],nodes))})
    save(pd.DataFrame(stability),out/'threshold_stability.tsv')
    g=graphs[.05];communities=parts[.05]
    primary=labels_from(communities,nodes)
    obs=nx.community.modularity(g,communities,weight='weight') if g.number_of_edges() else 0.
    rng=np.random.default_rng(SEED);weights=[d['weight'] for _,_,d in g.edges(data=True)];null=[]
    for b in range(999):
        h=g.copy()
        if h.number_of_edges()>=2:
            try:nx.double_edge_swap(h,nswap=5*h.number_of_edges(),max_tries=100*h.number_of_edges(),seed=SEED+b)
            except (nx.NetworkXAlgorithmError,nx.NetworkXError) as exc:raise RuntimeError('Degree-preserving null failed') from exc
            shuffled=rng.permutation(weights)
            for (u,v),w in zip(h.edges,shuffled):h[u][v]['weight']=float(w)
            assert sorted(dict(h.degree()).values())==sorted(dict(g.degree()).values())
            pc=partition(h,SEED+b)
            null.append(nx.community.modularity(h,pc,weight='weight'))
        else:null.append(0.)
        if (b+1)%100==0:print('C5 modularity null',b+1,flush=True)
    k1p=(1+np.count_nonzero(np.asarray(null)>=obs-1e-12))/1000
    save(pd.DataFrame({'draw':np.arange(999),'modularity_Q':null}),out/'modularity_null.tsv')
    mapping=pd.DataFrame({'feature':nodes,'module':['module'+str(i+1) for i in primary],'occurrence_fraction':occurrence[keep]})
    save(mapping,out/'module_membership.tsv')
    center=z[ek].mean(axis=0);zz=z-center
    scores=np.full((len(md),len(communities)),np.nan)
    nodeidx={v:i for i,v in enumerate(nodes)}
    for m,members in enumerate(communities):
        ee=[(nodeidx[u],nodeidx[v],d['weight']) for u,v,d in g.subgraph(members).edges(data=True)]
        if not ee:continue
        u,v,w=np.array(ee).T;u=u.astype(int);v=v.astype(int)
        scores[:,m]=(np.minimum(expit(zz[:,u]),expit(zz[:,v]))@w)/sum(w)
    has=np.any(np.isfinite(scores),axis=1)
    typ=np.where(has,np.argmax(np.nan_to_num(scores,nan=-1),axis=1),-1)
    coverage=np.mean(x>0,axis=1)
    proj=pd.DataFrame({'sample':md.index,'dataset':md.dataset,'sample_type':md.sample_type,
        'community_type':['undefined' if i<0 else 'module'+str(i+1) for i in typ],
        'observed_network_fraction':coverage,'max_edge_satisfaction':np.nanmax(scores,axis=1)})
    for m in range(scores.shape[1]):proj['module'+str(m+1)+'_score']=scores[:,m]
    save(proj,out/'all_sample_projections.tsv')
    y=typ[ek];env=md.loc[ek,'sample_type'].to_numpy();groups=md.loc[ek,'cluster'].to_numpy()
    folds=list(GroupKFold(n_splits=5).split(y,groups=groups));observed_chi=chi_stat(y,env)
    pred=majority_cv(y,env,folds);observed_acc=ba(y,pred,np.unique(y))
    rng=np.random.default_rng(SEED);association=[];prediction=[]
    for b in range(9999):
        yp=shuffle_within(y,groups,rng);association.append(chi_stat(yp,env))
        if b<999:prediction.append(ba(yp,majority_cv(yp,env,folds),np.unique(yp)))
    assocp=(1+sum(v>=observed_chi-1e-12 for v in association))/10000
    predp=(1+sum(v>=observed_acc-1e-12 for v in prediction))/1000
    save(pd.DataFrame({'draw':np.arange(9999),'chi_square':association}),out/'type_environment_association_null.tsv')
    save(pd.DataFrame({'draw':np.arange(999),'balanced_accuracy':prediction}),out/'environment_prediction_null.tsv')
    save(pd.DataFrame({'sample':md.index[ek],'type_code':y,'predicted_type_code':pred,'environment':env,'cluster':groups}),out/'environment_oof_predictions.tsv')
    study=md.dataset.eq('ours').to_numpy();k4=np.all(coverage[study]>=.30)&np.all(has[study])&np.all(np.nanmax(scores[study],axis=1)>0)
    criteria=[{'criterion':'K1','pass':k1p<=.05 and sum(len(c)>=5 for c in communities)>=2,'observed':obs,'p':k1p,'detail':f'{sum(len(c)>=5 for c in communities)} modules >=5 genera; 999 degree-preserving nulls'},
        {'criterion':'K2','pass':min(s['ARI'] for s in stability)>=.7,'observed':min(s['ARI'] for s in stability),'p':np.nan,'detail':'minimum pairwise ARI; fixed common node universe'},
        {'criterion':'K3','pass':assocp<=.05 and predp<=.05,'observed':observed_acc,'p':max(assocp,predp),'detail':f'project-restricted association p={assocp}; prediction p={predp}; grouped CV; Monte Carlo, not exhaustive exact'},
        {'criterion':'K4','pass':bool(k4),'observed':min(coverage[study]),'p':np.nan,'detail':'all ours15: >=30% observed network genera AND defined positive score; d2024 projections separately tabulated'}]
    save(pd.DataFrame(criteria),out/'C5_criteria_v1.tsv')
    save(pd.DataFrame({'feature':nodes,'eukbank_mean_clr':center}),out/'reference_centering.tsv')
    print('C5 COMPLETE',json.dumps(criteria,default=str),flush=True)

def self_test():
    assert np.allclose(clr([[1,2,3],[2,3,4]]).sum(axis=1),0)
    g=nx.path_graph(6);nx.set_edge_attributes(g,1.,'weight')
    p=partition(g,SEED);assert sorted(sum([list(c) for c in p],[]))==list(g)
    yy=np.array([0,0,1,1]);gg=np.array([0,0,1,1]);assert np.array_equal(shuffle_within(yy,gg,np.random.default_rng(1)),yy)
    assert chi_stat(yy,np.array(['a','a','b','b']))>0
    print('PASS CLR, partition, blocked permutation, association tests')

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('--version',action='version',version=VERSION)
    p.add_argument('--input',type=Path);p.add_argument('--out',type=Path);p.add_argument('--mode',choices=['rf','network','self-test'],required=True)
    p.add_argument('--rf-permutations',type=int,default=99)
    a=p.parse_args()
    if a.mode=='self-test':self_test()
    elif a.mode=='rf':rf_run(a.input,a.out,a.rf_permutations)
    else:network_run(a.input,a.out)
