#!/usr/bin/env python
"""Run A ecological-context RF comparison. Reads only frozen Run A inputs."""
from __future__ import annotations

import hashlib, json, platform, sys
from pathlib import Path
from datetime import datetime, timezone
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from sklearn.ensemble import RandomForestClassifier
from sklearn.dummy import DummyClassifier
from sklearn.model_selection import LeaveOneGroupOut, GroupKFold
from sklearn.metrics import accuracy_score, balanced_accuracy_score, f1_score, confusion_matrix, recall_score
import sklearn

from paths import INPUT_ROOT, OUTPUT_ROOT, require_run
ROOT = OUTPUT_ROOT / "rf_lopo"
INPUT = INPUT_ROOT / "comparative/inputs"
SEED = 20260924
N_BOOT = 1999
CLASSES = ["sediment_sand", "terrestrial_soil", "water_column"]
MAP_CONTEXT = {"land_sediment": "sediment_sand", "marine_sediment": "sediment_sand",
               "land_soil": "terrestrial_soil", "land_freshwater": "water_column",
               "marine_water": "water_column"}

def sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as fh:
        for block in iter(lambda: fh.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()

def write_tsv(frame, name):
    frame.to_csv(ROOT / name, sep="\t", index=False)

def metrics_from_cm(cm):
    """Metrics using only true classes present; avoids single-class sklearn warnings."""
    total = cm.sum()
    row = cm.sum(axis=1)
    col = cm.sum(axis=0)
    present = row > 0
    recall = np.divide(np.diag(cm), row, out=np.zeros(3, dtype=float), where=present)
    f1 = np.divide(2 * np.diag(cm), row + col, out=np.zeros(3, dtype=float), where=(row + col) > 0)
    return {"accuracy": float(np.trace(cm) / total),
            "balanced_accuracy": float(recall[present].mean()),
            "macro_f1": float(f1[present].mean())}

def project_weighted(y, pred, projects):
    rows=[]
    for prj in sorted(pd.unique(projects)):
        ix = projects == prj
        cm = confusion_matrix(y[ix], pred[ix], labels=CLASSES)
        vals = metrics_from_cm(cm)
        vals.update(project=prj, n_samples=int(ix.sum()), n_true_classes=int((cm.sum(axis=1) > 0).sum()),
                    balanced_accuracy_definition="mean recall across true classes present in this project")
        rows.append(vals)
    detail=pd.DataFrame(rows)
    return detail, {f"equal_project_{k}": float(detail[k].mean()) for k in ["accuracy","balanced_accuracy","macro_f1"]}

def bootstrap(y, pred, projects, rng):
    """Project-block resampling of fixed OOF predictions; no model refitting."""
    uniques=np.array(sorted(pd.unique(projects)))
    project_cms=np.stack([confusion_matrix(y[projects == p], pred[projects == p], labels=CLASSES) for p in uniques])
    project_metrics=np.array([[metrics_from_cm(cm)[key] for key in ["accuracy", "balanced_accuracy", "macro_f1"]] for cm in project_cms])
    vals=[]
    for b in range(N_BOOT):
        chosen=rng.integers(0, len(uniques), size=len(uniques))
        met=metrics_from_cm(project_cms[chosen].sum(axis=0))
        eq=project_metrics[chosen].mean(axis=0)
        for k,key in enumerate(["accuracy", "balanced_accuracy", "macro_f1"]):
            met[f"equal_project_{key}"]=float(eq[k])
        met["bootstrap_id"]=b+1
        vals.append(met)
    return pd.DataFrame(vals)

def features_for_fold(counts, train_ix):
    train=counts[:, train_ix]
    rel=train / train.sum(axis=0, keepdims=True)
    prev=(rel > 0).mean(axis=1)
    valid=(prev >= 0.05) & (np.var(rel, axis=1) > 0)
    if valid.sum() == 0:
        raise RuntimeError("No feature passed training-fold prevalence/variance filter")
    return valid

def run_cv(rank, counts, feature_names, meta, cv_name, splitter, groups):
    y=meta["ecological_context"].to_numpy()
    project=meta["project"].to_numpy()
    n=len(meta)
    records=[]; folds=[]; importance_rows=[]
    for fold, (train_ix, test_ix) in enumerate(splitter.split(np.zeros(n), y, groups), start=1):
        valid=features_for_fold(counts, train_ix)
        allrel=counts / counts.sum(axis=0, keepdims=True)
        Xtr=allrel[valid][:, train_ix].T; Xte=allrel[valid][:, test_ix].T
        train_y=y[train_ix]; test_y=y[test_ix]
        if set(CLASSES) - set(train_y):
            raise RuntimeError(f"{rank}/{cv_name} fold {fold} lacks training class(es): {set(CLASSES)-set(train_y)}")
        if len(np.intersect1d(meta.iloc[train_ix]["sample"], meta.iloc[test_ix]["sample"])):
            raise RuntimeError("Train/test sample overlap")
        if cv_name == "project_LOPO" and set(project[train_ix]) & set(project[test_ix]):
            raise RuntimeError("Project leakage in LOPO")
        foldbase=dict(rank=rank, cv=cv_name, fold=fold, n_train=len(train_ix), n_test=len(test_ix),
                      n_features=int(valid.sum()), train_projects=int(pd.Series(project[train_ix]).nunique()),
                      test_projects=int(pd.Series(project[test_ix]).nunique()),
                      train_classes=";".join(sorted(pd.unique(train_y))), test_classes=";".join(sorted(pd.unique(test_y))),
                      train_test_sample_overlap=0,
                      train_test_project_overlap=int(len(set(project[train_ix]) & set(project[test_ix]))))
        folds.append(foldbase)
        models={"RF": RandomForestClassifier(n_estimators=100, class_weight="balanced", max_features="sqrt", random_state=SEED, n_jobs=1),
                "Dummy_most_frequent": DummyClassifier(strategy="most_frequent")}
        for model_name, model in models.items():
            model.fit(Xtr, train_y)
            pred=model.predict(Xte)
            raw=model.predict_proba(Xte)
            proba=np.zeros((len(test_ix), len(CLASSES)))
            for col, klass in enumerate(model.classes_): proba[:, CLASSES.index(klass)]=raw[:, col]
            for j, ix in enumerate(test_ix):
                rec={"rank":rank,"cv":cv_name,"model":model_name,"fold":fold,"sample":meta.iloc[ix]["sample"],
                     "dataset":meta.iloc[ix]["dataset"],"project":project[ix],"cluster":meta.iloc[ix]["cluster"],
                     "true_context":y[ix],"predicted_context":pred[j],"correct":int(pred[j]==y[ix])}
                rec.update({f"prob_{c}":proba[j,k] for k,c in enumerate(CLASSES)})
                records.append(rec)
            if model_name == "RF":
                # Fold-level MDI is retained only as diagnostics; full fit below is the requested output.
                pass
    oof=pd.DataFrame(records)
    expected=len(meta)
    check_rows=[]
    for model in sorted(oof.model.unique()):
        z=oof[oof.model==model]
        ps=z[[f"prob_{c}" for c in CLASSES]].sum(axis=1)
        check_rows.append({"rank":rank,"cv":cv_name,"model":model,"expected_samples":expected,"n_predictions":len(z),
                           "unique_samples":z["sample"].nunique(),"duplicate_predictions":int(z["sample"].duplicated().sum()),
                           "missing_samples":expected-z["sample"].nunique(),"probability_sum_min":float(ps.min()),"probability_sum_max":float(ps.max())})
    return oof, pd.DataFrame(folds), pd.DataFrame(check_rows)

def full_importance(rank, counts, feature_names, meta):
    rel=counts / counts.sum(axis=0, keepdims=True)
    valid=(rel > 0).mean(axis=1) >= .05
    valid &= np.var(rel, axis=1) > 0
    model=RandomForestClassifier(n_estimators=100, class_weight="balanced",max_features="sqrt",random_state=SEED,n_jobs=1)
    model.fit(rel[valid].T, meta.ecological_context)
    out=pd.DataFrame({"rank":rank,"feature":np.array(feature_names)[valid],"mdi_importance":model.feature_importances_,
                      "prevalence_all_samples":(rel[valid]>0).mean(axis=1),"full_training_n_features":int(valid.sum()),
                      "label":"EXPLORATORY_BIASED_FULL_TRAINING_MDI_NOT_VALIDATION_NOT_MARKER"})
    return out.sort_values("mdi_importance",ascending=False)

def make_figure(summary, confusions):
    mpl.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Helvetica","DejaVu Sans"],"font.size":8,
                         "svg.fonttype":"none","pdf.fonttype":42,"axes.titleweight":"bold"})
    fig=plt.figure(figsize=(210/25.4,297/25.4), constrained_layout=True)
    gs=fig.add_gridspec(3,2, height_ratios=[1,1,1.25])
    colors={"RF":"#1F5A94","Dummy_most_frequent":"#8C8C8C"}
    metrics=[("accuracy","Accuracy"),("balanced_accuracy","Balanced accuracy"),("macro_f1","Macro-F1")]
    order=[("Genus","project_LOPO"),("Family","project_LOPO"),("Genus","cluster_GroupKFold5"),("Family","cluster_GroupKFold5")]
    for i,(metric,title) in enumerate(metrics):
        # Panel c deliberately spans the row so no near-empty panel is created.
        ax=fig.add_subplot(gs[i//2,i%2] if i < 2 else gs[1,:])
        sub=summary[summary["metric"]==metric].copy()
        for mi,model in enumerate(["RF","Dummy_most_frequent"]):
            vals=[]; low=[]; high=[]
            for rank,cv in order:
                q=sub[(sub["rank"]==rank)&(sub["cv"]==cv)&(sub["model"]==model)].iloc[0]
                vals.append(q["pooled_value"]); low.append(q["bootstrap_ci_low"]); high.append(q["bootstrap_ci_high"])
            x=np.arange(4)+(mi-.5)*.18
            ax.vlines(x, low, high, color=colors[model], lw=1.1)
            ax.scatter(x, vals, color=colors[model], s=18, zorder=3, label=("RF" if i==0 and mi==0 else "Dummy: most frequent" if i==0 else None))
        ax.set_ylim(0,1.05); ax.set_xticks(np.arange(4),["Genus\nLOPO","Family\nLOPO","Genus\nCluster-5","Family\nCluster-5"],fontsize=8)
        ax.set_ylabel(title); ax.set_title(chr(97+i)+"  "+title,loc="left"); ax.grid(axis="y",color="#DDDDDD",lw=.6); ax.spines[["top","right"]].set_visible(False)
        if i==0: ax.legend(frameon=False,loc="upper center",bbox_to_anchor=(.5,1.18),ncol=2,
                           title="Points: pooled OOF; bars: 95% project-bootstrap CI (1,999)")
    for j,rank in enumerate(["Genus","Family"]):
        ax=fig.add_subplot(gs[2,j])
        cm, den=confusions[rank]
        img=ax.imshow(cm, cmap="Greens", vmin=0, vmax=max(1,cm.max()))
        for r in range(3):
            for c in range(3):
                ax.text(c,r,f"{cm[r,c]} / {den[r]}",ha="center",va="center",fontsize=8,
                        color="white" if cm[r,c] > cm.max()*.55 else "#111111")
        ax.set_xticks(range(3),["sediment/\nsand","terrestrial\nsoil","water\ncolumn"]); ax.set_yticks(range(3),["sediment/sand","terrestrial soil","water column"])
        ax.set_xlabel("Predicted operational context"); ax.set_ylabel("True operational context")
        ax.set_title(chr(100+j)+f"  {rank} RF: project-LOPO confusion",loc="left")
        for edge in np.arange(-.5,3,1): ax.axhline(edge,color="white",lw=1); ax.axvline(edge,color="white",lw=1)
    fig.suptitle("Operational ecological-context classification from relative composition\n(not a particle-attachment measurement)",fontsize=11,fontweight="bold")
    for ext in ["png","svg","pdf"]:
        fig.savefig(ROOT / f"RF_ecological_context_LOPO.{ext}",dpi=600 if ext=="png" else None,transparent=False)
    plt.close(fig)

def main():
    ROOT.mkdir(parents=True,exist_ok=True)
    meta=pd.read_csv(INPUT/"metadata.tsv",sep="\t")
    meta["ecological_context"]=meta.sample_type.map(MAP_CONTEXT)
    if meta.ecological_context.isna().any(): raise RuntimeError("Unmapped sample_type")
    expected={"sediment_sand":71,"terrestrial_soil":13,"water_column":56}
    if meta.ecological_context.value_counts().to_dict()!=expected: raise RuntimeError(f"Unexpected contexts: {meta.ecological_context.value_counts().to_dict()}")
    if len(meta)!=140 or meta.project.nunique()!=34 or meta.cluster.nunique()!=54: raise RuntimeError("Input sample/project/cluster contract failed")
    hashes={p.name:sha256(p) for p in [INPUT/"metadata.tsv",INPUT/"counts_Genus.tsv",INPUT/"counts_Family.tsv"]}
    (ROOT/"input_sha256_v1.json").write_text(json.dumps(hashes,indent=2)+"\n",encoding="utf-8")
    write_tsv(meta,"metadata_ecological_context_v1.tsv")
    all_oof=[]; all_folds=[]; all_checks=[]; all_imp=[]
    for rank in ["Genus","Family"]:
        tab=pd.read_csv(INPUT/f"counts_{rank}.tsv",sep="\t")
        samples=tab.columns[1:].tolist()
        if samples != meta["sample"].tolist(): raise RuntimeError(f"{rank}: count columns are not exactly metadata sample order")
        counts=tab.iloc[:,1:].to_numpy(dtype=float); features=tab.iloc[:,0].astype(str).to_numpy()
        if not np.isfinite(counts).all() or (counts<0).any() or (counts.sum(axis=0)<=0).any(): raise RuntimeError(f"{rank}: invalid counts")
        designs=[("project_LOPO",LeaveOneGroupOut(),meta.project.to_numpy()),("cluster_GroupKFold5",GroupKFold(n_splits=5),meta.cluster.to_numpy())]
        for cv,splitter,groups in designs:
            oof,folds,checks=run_cv(rank,counts,features,meta,cv,splitter,groups)
            all_oof.append(oof); all_folds.append(folds); all_checks.append(checks)
        all_imp.append(full_importance(rank,counts,features,meta))
    oof=pd.concat(all_oof,ignore_index=True); folds=pd.concat(all_folds,ignore_index=True); checks=pd.concat(all_checks,ignore_index=True); imp=pd.concat(all_imp,ignore_index=True)
    if not ((checks["n_predictions"]==140)&(checks["unique_samples"]==140)&(checks["duplicate_predictions"]==0)&(checks["missing_samples"]==0)).all(): raise RuntimeError("OOF completeness failure")
    write_tsv(oof,"oof_predictions_probabilities_v1.tsv"); write_tsv(folds,"fold_train_test_checks_v1.tsv"); write_tsv(checks,"oof_completeness_checks_v1.tsv"); write_tsv(imp,"exploratory_full_training_MDI_importance_v1.tsv")
    summaries=[]; project_rows=[]; boot_rows=[]; confusions={}
    rng=np.random.default_rng(SEED)
    for (rank,cv,model),z in oof.groupby(["rank","cv","model"],sort=True):
        y=z["true_context"].to_numpy(); pred=z["predicted_context"].to_numpy(); projects=z["project"].to_numpy()
        pooled=metrics_from_cm(confusion_matrix(y,pred,labels=CLASSES)); pdet,pweighted=project_weighted(y,pred,projects); pdet.insert(0,"model",model);pdet.insert(0,"cv",cv);pdet.insert(0,"rank",rank);project_rows.append(pdet)
        boots=bootstrap(y,pred,projects,rng); boots.insert(0,"model",model);boots.insert(0,"cv",cv);boots.insert(0,"rank",rank);boot_rows.append(boots)
        for metric,value in {**pooled,**pweighted}.items():
            q=boots[metric].quantile([.025,.975]).to_numpy()
            summaries.append({"rank":rank,"cv":cv,"model":model,"metric":metric,"pooled_value":value,"bootstrap_ci_low":q[0],"bootstrap_ci_high":q[1],"bootstrap_replicates":N_BOOT,"interval_definition":"project-resampling of fixed OOF predictions; conditional on fitted models"})
        cm=confusion_matrix(y,pred,labels=CLASSES); rec=recall_score(y,pred,labels=CLASSES,average=None,zero_division=0)
        for k,c in enumerate(CLASSES):
            summaries.append({"rank":rank,"cv":cv,"model":model,"metric":"class_recall","class":c,"pooled_value":rec[k],"bootstrap_ci_low":np.nan,"bootstrap_ci_high":np.nan,"bootstrap_replicates":N_BOOT,"interval_definition":"class recall (interval not requested)"})
        if cv=="project_LOPO" and model=="RF": confusions[rank]=(cm,cm.sum(axis=1))
        cmout=pd.DataFrame(cm,index=CLASSES,columns=CLASSES); cmout.index.name="true_context"; cmout.to_csv(ROOT/f"confusion_{rank}_{cv}_{model}_v1.tsv",sep="\t")
    summary=pd.DataFrame(summaries); write_tsv(summary,"performance_summary_project_bootstrap_v1.tsv"); write_tsv(pd.concat(project_rows,ignore_index=True),"per_project_performance_v1.tsv"); write_tsv(pd.concat(boot_rows,ignore_index=True),"project_bootstrap_OOF_metrics_v1.tsv")
    # Plot source is a compact copy of the plotted pooled metrics/intervals.
    write_tsv(summary[(summary["metric"].isin(["accuracy","balanced_accuracy","macro_f1"]))],"figure_performance_source_data_v1.tsv")
    make_figure(summary,confusions)
    env={"timestamp_utc":datetime.now(timezone.utc).isoformat(),"python":sys.version,"platform":platform.platform(),"numpy":np.__version__,"pandas":pd.__version__,"scikit_learn":sklearn.__version__,"matplotlib":mpl.__version__,"seed":SEED,"bootstrap_replicates":N_BOOT,"rf":{"n_estimators":100,"class_weight":"balanced","max_features":"sqrt","n_jobs":1}}
    (ROOT/"environment_versions_v1.json").write_text(json.dumps(env,indent=2)+"\n",encoding="utf-8")
    manifest=[]
    for p in sorted(ROOT.iterdir()):
        if p.is_file(): manifest.append({"file":p.name,"sha256":sha256(p),"bytes":p.stat().st_size})
    write_tsv(pd.DataFrame(manifest),"output_manifest_sha256_v1.tsv")
    print("SUCCESS",len(oof),"OOF rows",len(summary),"summary rows")

if __name__ == "__main__":
    require_run()
    main()
