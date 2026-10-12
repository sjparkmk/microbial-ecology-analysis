"""Numeric-only extraction of the frozen trophic-composition analysis.
Plot layout and private style helpers intentionally omitted. No new rules.
"""
from pathlib import Path
import re, hashlib
import numpy as np
import pandas as pd
from paths import INPUT_ROOT, OUTPUT_ROOT, require_run
ROOT = OUTPUT_ROOT / 'composition'
INPUT = INPUT_ROOT / 'kss/inputs'
def sha(p): return hashlib.sha256(p.read_bytes()).hexdigest()
SOURCE = 'Adl et al. (2019); doi:10.1111/jeu.12691; Appendix S2'

TYPES = ['Protist-enriched', 'Metazoa', 'Fungi', 'Land plants', 'Other / unresolved']

MODES = ['Heterotroph-associated', 'Phototroph-associated', 'Mixotroph-associated', 'Variable / insufficient resolution', 'Unresolved taxonomy']

def tsv(x, name):
    x.to_csv(ROOT / name, sep='\t', index=True)


def named(g):
    return bool(re.fullmatch('[A-Z][a-z]+', g))


def annotate(row):
    d,s,c,o,g = [str(row[k]) for k in ['Division','Subdivision','Class','Order','Genus']]
    if s == 'Metazoa': organism = 'Metazoa'
    elif s == 'Fungi': organism = 'Fungi'
    elif c == 'Embryophyceae': organism = 'Land plants'
    elif not d or (d == 'Opisthokonta' and s not in ['Choanoflagellata','Ichthyosporea','Nucleariida','Rotosphaerida']):
        organism = 'Other / unresolved'
    else: organism = 'Protist-enriched'

    # Conservative affinity labels; no claim of obligate mode or activity.
    mode, rule = 'Variable / insufficient resolution', 'U01_no_specific_mode_rule'
    if organism == 'Other / unresolved': mode,rule = 'Unresolved taxonomy','U00_ambiguous_taxonomic_context'
    elif s in ['Metazoa','Fungi']: mode,rule = 'Heterotroph-associated','H01_animal_fungal_carbon_affinity_not_feeding_guild'
    elif g in ['Helicosporidium','Polytomella','Polytoma','Prototheca']:
        mode,rule = 'Heterotroph-associated','H02_nonphotosynthetic_green_genus'
    elif d == 'Chlorophyta' and named(g):
        mode,rule = 'Phototroph-associated','P01_named_green_genus_exception_screened'
    elif c in ['Zygnemophyceae','Coleochaetophyceae','Klebsormidiophyceae']:
        mode,rule = 'Phototroph-associated','P02_streptophyte_algal_class'
    elif c == 'Embryophyceae' and named(g):
        # Restrict terrestrial plant calls to named genera present in these data.
        if g in ['Nymphaea','Pinus','Phragmites','Typha','Quercus','Carex','Juncus','Persicaria','Potamogeton','Lemna','Sphagnum','Polytrichum']:
            mode,rule = 'Phototroph-associated','P03_named_photosynthetic_plant_genus'
    elif c in ['Mediophyceae','Coscinodiscophyceae']:
        mode,rule = 'Phototroph-associated','P04_diatom_class_affinity'
    elif c == 'Bacillariophyceae' and named(g) and g not in ['Nitzschia','Hantzschia']:
        mode,rule = 'Phototroph-associated','P05_named_pennate_diatom_conservative_screen'
    elif c in ['Eustigmatophyceae','Glaucocystophyceae']:
        mode,rule = 'Phototroph-associated','P06_phototrophic_class_affinity'
    elif c == 'Dinophyceae':
        if g in ['Palatinus','Asulcocephalium','Borghiella','Unruhdinium','Peridinium','Parvodinium']:
            mode,rule = 'Phototroph-associated','P07_dinoflagellate_named_genus_table2_1'
        elif g in ['Nusuttodinium','Baldinia']:
            mode,rule = 'Mixotroph-associated','M01_dinoflagellate_named_genus_table2_1'
        else: rule = 'U02_dinoflagellate_variable_or_unresolved'
    elif c == 'Chrysophyceae':
        if g in ['Spumella','Pedospumella','Paraphysomonas']:
            mode,rule = 'Heterotroph-associated','H03_colourless_chrysophyte_genus'
        elif g in ['Mallomonas','Synura','Chrysodidymus','Tesselaria']:
            mode,rule = 'Phototroph-associated','P08_synuralean_genus'
        else: rule = 'U03_chrysophyte_variable_or_unresolved'
    elif s == 'Ciliophora':
        # Restrict to selected resolved genera with feeding entries inspected in Table2.2.
        if g in ['Loxodes','Geleia','Cryptopharynx','Spirostomum','Anigsteinia','Blepharisma','Condylostoma','Gruberia','Folliculina','Peritromus','Cyrtohymena','Gastrostyla','Stylonychia']:
            mode,rule = 'Heterotroph-associated','H04_ciliate_named_genus_table2_2'
        else: rule = 'U04_ciliate_mode_not_resolved_conservatively'
    elif s == 'Apicomplexa' or c in ['Endomyxa-Phytomyxea','Perkinsida'] or o == 'Pirsoniales':
        mode,rule = 'Heterotroph-associated','H05_parasite_endobiont_lineage_affinity'
    elif s == 'Choanoflagellata' or c in ['Bicoecea','Colpodellidea','Filosa-Sarcomonadea','Heterolobosea','Kinetoplastea'] or d == 'Centroplasthelida':
        mode,rule = 'Heterotroph-associated','H06_phagotroph_associated_lineage'
    elif c in ['Peronosporomycetes','Hyphochytriomyceta']:
        mode,rule = 'Heterotroph-associated','H07_osmotrophic_lineage_strategy_unresolved'
    parasite = 'Not assessed / not specifically assigned'
    if s == 'Apicomplexa': parasite = 'Apicomplexa (parasite / endobiont)'
    elif c == 'Endomyxa-Phytomyxea': parasite = 'Phytomyxea (parasite-associated)'
    elif c == 'Perkinsida': parasite = 'Perkinsida (parasite-associated)'
    elif o == 'Pirsoniales' or c == 'Pirsoniales': parasite = 'Pirsoniales (parasite-associated)'
    note = 'Taxonomic affinity only; within-lineage variation and activity unmeasured.'
    if mode.startswith('Variable'): note = 'Insufficient resolved taxonomic/trait evidence for one mode; not confirmed mixotrophy.'
    if mode == 'Unresolved taxonomy': note = 'Unresolved taxonomic context; not assumed to be protist or functionally absent.'
    return pd.Series(dict(organism_group=organism,potential_trophic_affinity=mode,mode_rule=rule,
                          selected_parasite_lineage=parasite,evidence=SOURCE,confidence='conservative_annotation_not_function_validation',uncertainty=note))


def aggregate(counts, ledger, field, order=None):
    out=counts.groupby(ledger[field], sort=False).sum()
    if order is not None: out=out.reindex(order,fill_value=0)
    return out


def main():
    ROOT.mkdir(parents=True,exist_ok=True)
    inputs=[INPUT/x for x in ['otutab_rare_20000.tsv','taxonomy_pr2v511.tsv','metadata_analysis.tsv']]
    hashes={str(p):sha(p) for p in inputs}
    meta=pd.read_csv(inputs[2],sep='\t',index_col='SampleID')
    meta=meta.loc[meta.IncludeInGroupTests.eq('yes')].copy()
    meta['_order']=range(len(meta));meta['_exposure']=meta.AnalysisGroup.map({'submerged':0,'exposed':1});meta['_substrate']=meta.sample_type.map({'sediment':0,'sand':1,'soil':2})
    meta=meta.sort_values(['_exposure','_substrate','_order']).drop(columns=['_exposure','_substrate','_order'])
    counts=pd.read_csv(inputs[0],sep='\t',index_col=0).loc[:,meta.index]
    tax=pd.read_csv(inputs[1],sep='\t',index_col=0).fillna('').reindex(counts.index)
    assert len(meta)==15 and 'swpw' not in meta.index
    assert counts.index.is_unique and counts.columns.is_unique
    assert (counts.to_numpy()>=0).all() and (counts.to_numpy()==counts.to_numpy().astype(int)).all()
    assert counts.sum().eq(20000).all() and counts.to_numpy().sum()==300000
    assert counts.sum(axis=1).gt(0).sum()==4586
    assert tax.Domain.eq('Eukaryota').all()
    ledger=pd.concat([tax,tax.apply(annotate,axis=1)],axis=1)
    ledger['study_20k_total_reads']=counts.sum(axis=1);ledger['study_20k_prevalence_libraries']=counts.gt(0).sum(axis=1)
    ledger['retained_in_protist_enriched']=ledger.organism_group.eq('Protist-enriched')
    ledger.index.name='ASV'
    tsv(ledger,'ASV_taxonomy_trophic_annotation_ledger_v1.tsv');tsv(meta,'display_sample_order_metadata_v1.tsv')
    observed=ledger.study_20k_total_reads.gt(0)
    tsv(ledger.loc[observed].groupby(['organism_group','potential_trophic_affinity','mode_rule']).agg(reads=('study_20k_total_reads','sum'),observed_ASVs=('study_20k_total_reads','size')),'annotation_rule_totals_v1.tsv')
    groups=aggregate(counts,ledger,'organism_group',TYPES);group_pct=groups/20000*100
    prot=counts.loc[ledger.retained_in_protist_enriched]
    den=prot.sum();assert den.gt(0).all()
    proto_ledger=ledger.loc[prot.index]
    modes=aggregate(counts,ledger,'potential_trophic_affinity',MODES);mode_pct=modes/20000*100
    pmodes=aggregate(prot,proto_ledger,'potential_trophic_affinity',MODES);pmode_pct=pmodes.div(den,axis=1)*100
    assert groups.sum().eq(20000).all() and modes.sum().eq(20000).all() and pmodes.sum().equals(den)
    for name,data in [('organism_group_counts',groups),('organism_group_percent',group_pct),('all_euk_trophic_counts',modes),('all_euk_trophic_percent',mode_pct),('protist_trophic_counts',pmodes),('protist_trophic_percent',pmode_pct)]:tsv(data,name+'_v1.tsv')
    denominators=pd.DataFrame({'all_euk_reads':counts.sum(),'protist_enriched_reads':den,'protist_enriched_percent':den/20000*100,'observed_all_euk_ASVs':counts.gt(0).sum(),'observed_protist_ASVs':prot.gt(0).sum()})
    tsv(denominators,'library_denominators_v1.tsv')
    classes=proto_ledger.Class.where(proto_ledger.Class.ne(''),'Unresolved class')
    ccounts=prot.groupby(classes).sum();cpct=ccounts.div(den,axis=1)*100
    top=ccounts.sum(axis=1).sort_values(ascending=False).head(12).index
    displayed=cpct.loc[top].copy();displayed.loc['All other classes']=cpct.loc[~cpct.index.isin(top)].sum()
    assert np.allclose(displayed.sum(),100)
    tsv(ccounts,'protist_all_class_counts_v1.tsv');tsv(cpct,'protist_all_class_percent_v1.tsv');tsv(displayed,'figure_protist_class_display_matrix_v1.tsv')
    para=counts.groupby(ledger.selected_parasite_lineage).sum().drop(index='Not assessed / not specifically assigned',errors='ignore')
    tsv(para,'selected_parasite_lineage_counts_v1.tsv');tsv(para/20000*100,'selected_parasite_lineage_all_euk_percent_v1.tsv')

if __name__ == '__main__':
    require_run()
    main()
