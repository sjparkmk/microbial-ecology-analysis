\
\
\
\
\
\
\
\
\
\
\
\
\
\

import time
from Bio import AlignIO, Phylo
from Bio.Phylo.TreeConstruction import DistanceCalculator, DistanceTreeConstructor

ALIGN_IN = "/input/youre/data/all_663_asvs.aligned.fasta"
TREE_OUT = "/output/your/data/asv_tree_all663.newick"
SUMMARY_OUT = "/output/your/data/asv_tree_summary.tsv"

t0 = time.time()
aln = AlignIO.read(ALIGN_IN, "fasta")
print(f"Loaded alignment: {len(aln)} seqs x {aln.get_alignment_length()} cols ({time.time()-t0:.1f}s)")

t0 = time.time()
calculator = DistanceCalculator("identity")
dm = calculator.get_distance(aln)
print(f"Computed distance matrix ({time.time()-t0:.1f}s)")

t0 = time.time()
constructor = DistanceTreeConstructor()
tree = constructor.nj(dm)
print(f"Built NJ tree ({time.time()-t0:.1f}s)")


tree.root_with_outgroup({"name": "ASV_1"})


n_clipped = 0
for clade in tree.find_clades():
    if clade.branch_length is not None and clade.branch_length < 0:
        clade.branch_length = 0.0
        n_clipped += 1
print(f"Clipped {n_clipped} negative branch lengths to 0")

Phylo.write(tree, TREE_OUT, "newick")

n_tips = tree.count_terminals()
with open(SUMMARY_OUT, "w") as f:
    f.write("metric\tvalue\n")
    f.write(f"input\t{ALIGN_IN}\n")
    f.write(f"output\t{TREE_OUT}\n")
    f.write(f"sequences\t{len(aln)}\n")
    f.write(f"alignment_columns\t{aln.get_alignment_length()}\n")
    f.write(f"tree_tips\t{n_tips}\n")
    f.write("distance_model\tidentity (p-distance / Hamming on full MSA)\n")
    f.write("rooting\toutgroup = ASV_1 (most abundant ASV, matching original script convention)\n")
    f.write(f"negative_branch_lengths_clipped\t{n_clipped}\n")

print(f"Wrote tree: {TREE_OUT}")
print(f"Wrote summary: {SUMMARY_OUT}")
print("n_tips", n_tips)
