from __future__ import annotations
import argparse, csv, hashlib
from pathlib import Path
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap, to_rgba
from matplotlib.patches import Rectangle
from scipy.cluster.hierarchy import dendrogram, leaves_list, linkage
from scipy.spatial.distance import pdist, squareform
import numpy as np

MM=1/25.4
A4=(297*MM,210*MM)
SAMPLES=["1st_BFP","1st_BFW1","1st_BFW2","2nd_BFP","2nd_BFW","2nd_wt","3rd_wt"]
COLORS={"P_biofilm":"#1BAF7A","W_biofilm":"#EB6834","water":"#2A78D6"}
LABELS={"P_biofilm":"AMD precipitate-associated biofilm","W_biofilm":"wall-attached biofilm","water":"water"}
STEM="SFig2_BrayCurtis_Heatmap_A4_v3"

def read(path, delim="\t"):
    with path.open(encoding="utf-8-sig",newline="") as h: return list(csv.DictReader(h,delimiter=delim))
def sha(path):
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1048576),b""):h.update(b)
    return h.hexdigest().upper()
def write(path, rows):
    path.parent.mkdir(parents=True,exist_ok=True)
    with path.open("w",encoding="utf-8",newline="") as h:
        w=csv.writer(h,delimiter="\t",lineterminator="\n"); w.writerows(rows)
def main():
    p=argparse.ArgumentParser();p.add_argument("--source-root",type=Path,default=Path("/input/youre/data/"));p.add_argument("--output-root",type=Path,default=Path("/output/your/data"));a=p.parse_args()
    root=a.output_root.resolve(); src=a.source_root.resolve()
    rare=root/"source_data"/"SFig2_rarefied_ASV_table_v1.tsv"; meta=src/"input_data"/"sample_sheet.tsv"
    figdir=root/"build"/"bray_rarefied_stage"/"figures"; logdir=root/"logs"; sourcedir=root/"source_data"
    paths={"pdf":figdir/f"{STEM}.pdf","svg":figdir/f"{STEM}.svg","png":figdir/f"{STEM}.png","tiff":figdir/f"{STEM}.tiff","matrix":sourcedir/f"{STEM}_matrix.tsv","order":sourcedir/f"{STEM}_display_order.tsv","log":logdir/f"{STEM}_run.log","manifest":logdir/f"{STEM}_manifest.tsv"}
    existing=[str(x) for x in paths.values() if x.exists()]
    if existing: raise FileExistsError("Refusing to overwrite: "+"; ".join(existing))
    rows=read(rare); mrows=read(meta); md={r["analysis_name"]:r for r in mrows}
    if len(rows)!=663 or set(md)!=set(SAMPLES):raise ValueError("Unexpected ASV rows or metadata samples")
    mat=np.array([[int(r[s]) for r in rows] for s in SAMPLES])
    if not np.all(mat.sum(axis=1)==88397): raise ValueError("Rarefied depths differ from 88,397")
    dist=squareform(pdist(mat,metric="braycurtis")); cond=squareform(dist,checks=True); cl=linkage(cond,method="average"); idx=leaves_list(cl).tolist(); ordered=[SAMPLES[i] for i in idx]; show=dist[np.ix_(idx,idx)]
    write(paths["matrix"],[["sample"]+SAMPLES]+[[s]+[f"{x:.12f}" for x in row] for s,row in zip(SAMPLES,dist)])
    write(paths["order"],[["display_order","sample","Type_code"]]+[[i+1,s,md[s]["Type"]] for i,s in enumerate(ordered)])
    plt.rcParams.update({"font.family":"sans-serif","font.sans-serif":["Arial","Liberation Sans","DejaVu Sans"],"font.size":8,"axes.edgecolor":"#000000","axes.linewidth":1.15,"svg.fonttype":"none","pdf.fonttype":42})
    fig=plt.figure(figsize=A4)
    left,bottom,height=0.305,0.145,0.615;width=height*210/297
    td=fig.add_axes([left,0.815,width,0.12]);dendrogram(cl,ax=td,color_threshold=0,above_threshold_color="#4A4946",no_labels=True);td.axis("off")
    ts=fig.add_axes([left,0.775,width,0.030]);rgba=np.array([[to_rgba(COLORS[md[s]["Type"]]) for s in ordered]]);ts.imshow(rgba,aspect="auto",interpolation="nearest",extent=(-.5,6.5,-.5,.5));ts.set(xlim=(-.5,6.5),ylim=(-.5,.5),xticks=[],yticks=[])
    for b in np.arange(.5,6.5,1):ts.axvline(b,color="white",lw=.7)
    for spine in ts.spines.values():spine.set(color="#000000",linewidth=.85)
    ld=fig.add_axes([.135,bottom,.120,height]);dendrogram(cl,ax=ld,orientation="left",color_threshold=0,above_threshold_color="#4A4946",no_labels=True);ld.axis("off")
    ls=fig.add_axes([.270,bottom,.023,height]);rgba=np.array([[to_rgba(COLORS[md[s]["Type"]])] for s in ordered]);ls.imshow(rgba,aspect="auto",interpolation="nearest",extent=(-.5,.5,6.5,-.5));ls.set(xlim=(-.5,.5),ylim=(6.5,-.5),xticks=[],yticks=[])
    for b in np.arange(.5,6.5,1):ls.axhline(b,color="white",lw=.7)
    for spine in ls.spines.values():spine.set(color="#000000",linewidth=.85)
    ax=fig.add_axes([left,bottom,width,height]);cmap=LinearSegmentedColormap.from_list("bray",["#F7F7F5","#B8D0EA","#4A8BDA","#174A7E"]);im=ax.imshow(show,vmin=0,vmax=1,cmap=cmap,aspect="equal",interpolation="nearest");ax.set_xticks(range(7),ordered,rotation=43,ha="right",rotation_mode="anchor");ax.set_yticks(range(7),ordered);ax.yaxis.tick_right();ax.tick_params(axis="both",colors="#000000",width=1.0,labelsize=8.2);ax.tick_params(axis="y",pad=5)
    for sp in ax.spines.values():sp.set(color="#000000",linewidth=1.05)
    for i in range(7):
      for j in range(7):
       v=show[i,j];ax.text(j,i,f"{v:.2f}",ha="center",va="center",fontsize=7.3,color="white" if v>=.57 else "#3F3F3F")
    leg=fig.add_axes([left+width+.060,.58,.21,.32]);leg.axis("off");leg.text(0,.99,"Type",fontsize=10.5,fontweight="bold",va="top")
    for code,y in zip(["P_biofilm","W_biofilm","water"],[.74,.47,.20]):
      leg.add_patch(Rectangle((.01,y-.06),.09,.12,transform=leg.transAxes,color=COLORS[code],clip_on=False));txt=f"{code}\n{LABELS[code]}" if code!="water" else "water";leg.text(.15,y,txt,transform=leg.transAxes,va="center",fontsize=8.4,linespacing=1.05)
    ca=fig.add_axes([.865,.17,.024,.27]);cb=fig.colorbar(im,cax=ca);cb.set_ticks([0,.25,.5,.75,1]);cb.ax.tick_params(labelsize=8.2,width=1,colors="#000000");cb.set_label("Bray-Curtis dissimilarity",fontsize=10,fontweight="bold",labelpad=7);cb.outline.set(color="#000000",linewidth=1)
    figdir.mkdir(parents=True,exist_ok=True)
    fig.savefig(paths["pdf"],facecolor="white");fig.savefig(paths["svg"],facecolor="white");fig.savefig(paths["png"],dpi=320,facecolor="white");fig.savefig(paths["tiff"],dpi=400,facecolor="white",pil_kwargs={"compression":"tiff_lzw"});plt.close(fig)
    lines=["PASS: Descriptive rarefied Bray-Curtis heatmap rendered.","Input: R-rarefied ASV table; vegan::rrarefy seed=100; 88,397 reads/sample.","Distance: Bray-Curtis on rarefied ASV counts; no environmental or Type hypothesis test.","Display order: "+", ".join(ordered),"Mean upper-triangle distance: %.6f"%(dist[np.triu_indices(7,1)].mean())]
    paths["log"].write_text("\n".join(lines)+"\n",encoding="utf-8")
    manifest=[("plot_script",Path(__file__).resolve()),("input_rarefied_asv",rare),("input_metadata",meta)]+[("output_"+k,v) for k,v in paths.items() if k!="manifest"]
    write(paths["manifest"],[["role","path","bytes","sha256"]]+[[r,str(p),p.stat().st_size,sha(p)] for r,p in manifest])
    print("\n".join(lines))
if __name__=="__main__":main()
