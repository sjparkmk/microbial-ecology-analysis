# Public release manifest

Release scope: `scripts + config + provenance + validation + documentation`

## Included

- In-house Python and R scripts under `scripts/`.
- Analysis contracts and path templates under `config/`.
- Software/session/exclusion records under `provenance/`.
- The read-only public-release checker under `validation/`.
- README, license-boundary, and third-party notices.

## Excluded

- Raw FASTQ/FASTA files and sequencing archives.
- Processed ASV/count tables, sample metadata, and figure-source tables.
- Generated figures, tables, logs, caches, and output packages.
- PR2, SILVA, RDP, FAPROTAX, or any other reference database.
- USEARCH or other executable binaries.
- EasyAmplicon2 or other copied third-party implementations.
- Failed-QC artifacts and superseded RF/network figure renderers that do not
  represent the final manuscript contract.

## Provenance boundary

The final analysis starts from externally supplied common-75 inputs and uses
fixed depths of 20,000 (16S) and 50,000 (18S) reads per sample. The public
package documents these requirements but does not publish the underlying data.
NCM/tNST scripts are retained as result-bearing diagnostics and are gated
against positive assembly interpretation.

## Archival links

- GitHub: <https://github.com/sjparkmk/microbial-ecology-analysis>
- DOI/Zenodo/Figshare: to be added after archival deposit.
