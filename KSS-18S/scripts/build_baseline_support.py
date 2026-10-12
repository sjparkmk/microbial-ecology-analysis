"""Explicit recurrence/core definitions and individual coverage, without new hypothesis tests."""
from pathlib import Path
import argparse,hashlib,json
import numpy as np
import pandas as pd

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument('--out',type=Path,required=True);a=p.parse_args()
 from paths import INPUT_ROOT
 src=INPUT_ROOT/'kss';out=a.out;out.mkdir(parents=True,exist_ok=False)
 paths=[src/'inputs/otutab_rare_20000.tsv',src/'inputs/metadata_analysis.tsv',src/'inputs/taxonomy_pr2v511.tsv',src/'figure_source_data/Figure_3/panel_a_asv_abundance_prevalence.tsv']
 x=pd.read_csv(paths[0],sep='\t',index_col=0);m=pd.read_csv(paths[1],sep='\t');tax=pd.read_csv(paths[2],sep='\t',index_col=0)
 ids=m.loc[m.IncludeInGroupTests.eq('yes'),'SampleID'].tolist();x=x[ids];assert len(ids)==15 and (x.sum()==20000).all()
 d=pd.DataFrame({'asv_id':x.index,'detected_libraries_at_ge1_read':(x>0).sum(axis=1),'pooled_reads_study15':x.sum(axis=1)})
 d=d.reset_index(drop=True)
 d['pooled_relative_abundance_percent']=d.pooled_reads_study15/300000*100;d['strict_ubiquity_15_of_15']=d.detected_libraries_at_ge1_read.eq(15)
 d['source_recurrence_ge8_of15']=d.detected_libraries_at_ge1_read.ge(8);d['display_high_abundance_ge02percent']=d.pooled_reads_study15.ge(600)
 d=d.merge(tax.reset_index().rename(columns={tax.index.name or 'index':'asv_id'}),on='asv_id',validate='one_to_one')
 old=pd.read_csv(paths[3],sep='\t').set_index('asv_id');now=d.set_index('asv_id')
 assert (now.loc[old.index,'detected_libraries_at_ge1_read']==old.prevalence_study15).all()
 assert d.strict_ubiquity_15_of_15.sum()==0 and d.source_recurrence_ge8_of15.sum()==1 and d.display_high_abundance_ge02percent.sum()==66
 recurrent=d[d.source_recurrence_ge8_of15];assert recurrent.asv_id.tolist()==['ASV_000018']
 d.to_csv(out/'all_ASV_recurrence_taxonomy.tsv',sep='\t',index=False);recurrent.to_csv(out/'recurrent_ASV000018_taxonomy.tsv',sep='\t',index=False)
 x.loc[recurrent.asv_id].div(20000).mul(100).to_csv(out/'recurrent_ASV_individual_relative_abundance_percent.tsv',sep='\t')
 rel=x/20000;shannon=-(rel.where(rel>0)*np.log(rel.where(rel>0))).sum();observed=(x>0).sum()
 alpha=pd.DataFrame({'SampleID':ids,'depth':20000,'Observed_ASVs':observed,'Shannon':shannon,'Simpson_1_minus_D':1-(rel**2).sum(),
                     'Pielou':shannon/np.log(observed),'singletons':(x==1).sum(),'Goods_coverage':1-(x==1).sum()/20000})
 alpha=alpha.merge(m[['SampleID','AnalysisGroup','sample_type']],on='SampleID',validate='one_to_one')
 alpha.to_csv(out/'individual15_alpha_Pielou_Goods_coverage.tsv',sep='\t',index=False)
 (out/'README.md').write_text('''# Baseline supporting definitions

- Only the15 study libraries are used, each at the existing20,000-read depth. Detection means at least1 read in a rarefied library.
- Strict ubiquity means detection in15/15: no ASV meets this criterion. This is an operational result, not proof that a universal biological core does not exist.
- Source-defined recurrence means detection in at least8/15. Only ASV_000018 meets it (8/15); this does not mean it was detected in only one sample. Its taxonomy and all15 relative abundances are tabulated separately.
- High abundance means at least0.2% of the300,000 pooled reads (at least600 reads), yielding66 ASVs. This is not the recurrence/core criterion. The historical66-ASV tree (V5.13 Figure S8) is therefore a high-abundance short-amplicon display, not a core tree or organismal phylogeny.
- Recurrence is sensitive to denoising/pooling and detection thresholds; consult the historical denoising sensitivity record (not repeated-subsampling recurrence validation). Absence of ubiquitous ASVs does not itself prove microdiversity, dispersal limitation or habitat specificity.
- Good's coverage is1 − singleton reads / total reads; Pielou is Shannon / ln(observed ASVs). These descriptive values add no new tests. Main15 raw minimum is25,798, so28,000-read rarefaction would not retain all15 in this RunA dataset. The current author-approved20,000-read baseline remains fixed.
''',encoding='utf-8')
 (out/'source_manifest.json').write_text(json.dumps([{'path':str(q),'sha256':hashlib.sha256(q.read_bytes()).hexdigest()} for q in paths],indent=2),encoding='utf-8')
 print(f'PASS strict core0; recurrent1 ({recurrent.asv_id.iloc[0]}); display66; alpha15')
if __name__=='__main__':main()
