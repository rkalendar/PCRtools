# DigitalGens

**A free, fully client-side web server for nucleic-acid assay design: genotyping, PCR, isothermal
amplification, CRISPR, assembly and sequence analysis.**

[![License](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE.txt)
[![Language](https://img.shields.io/badge/Language-JavaScript-yellow.svg)]()

**Live server:** [https://digitalgens.org/](https://digitalgens.org/)

**Archived release:** [10.5281/zenodo.23498156](https://doi.org/10.5281/zenodo.23498156)

**Author:** Ruslan Kalendar, Helsinki Institute of Life Science (HiLIFE), University of Helsinki

**Version:** v1.0, 15 October 2026  

Everything runs in your browser. A pasted or uploaded sequence never leaves your machine, no account
is required, nothing is installed, and the only outbound traffic is an accession number or rsID you
type into an optional retrieval box.

---

## Contents

- [Getting started](#getting-started)
- [Tools](#tools)
  - [Genotyping](#genotyping)
  - [Amplification and panel design](#amplification-and-panel-design)
  - [CRISPR](#crispr)
  - [Assembly and synthesis](#assembly-and-synthesis)
  - [Oligonucleotide analysis](#oligonucleotide-analysis)
  - [Sequences, repeats and phylogeny](#sequences-repeats-and-phylogeny)
  - [Calculators](#calculators)
- [How the shared engine works](#how-the-shared-engine-works)
- [Citation](#citation)
- [Licence and third-party components](#licence-and-third-party-components)
- [Contact](#contact)

---

## Getting started

### Online

Open [https://digitalgens.org/](https://digitalgens.org/). Nothing to install.

### Offline

```bash
git clone https://github.com/rkalendar/PCRtools.git
cd PCRtools
```

Open `sites/index.html` in any modern browser. There is no build step and no dependency to fetch: all
computation is plain JavaScript and the two third-party packages are vendored in `sites/assets/vendor/`.
Only the optional NCBI and Ensembl retrieval boxes need a network connection.

### A worked test case

Open `sites/kasp.html?example=1` and press **Generate Primers**. Every tool page accepts `?example=1` and
loads with its own example already in place, so you can confirm the installation works without
supplying a sequence. Each page also has a one-click example button.

The output is a tab-delimited table beginning:

```
Location(ID)	Sequence(5'-3')	Length(nt)	Tm(°C)	CG(%)	Linguistic_Complexity(%)
```

Full documentation for every tool, parameter and output column is in `help.html`, which also carries
an end-to-end worked tutorial and a troubleshooting section.

---

## Tools

Twenty-nine tools. Those marked **[G]** are the genotyping modules described in the DigitalGens
paper (in preparation).

### Genotyping

| Tool | Page | What it does |
|---|---|---|
| **KASP / AS-PCR** **[G]** | `kasp.html` | Allele-specific primers for KASP, PACE, ASQ and generic AS-PCR, for SNPs, InDels and haplotypes, with automatic universal-tail selection and Tm balancing. An allele-sizing mode replaces the fluorophore tails with a length ladder so the design can be read on a gel. |
| **Dual-labelled qPCR probes** **[G]** | `qpcr.html` | TaqMan, minor-groove-binder, molecular-beacon and UniQ probes for allelic discrimination. Probe Tm is set relative to the primer window by a format-specific offset (TaqMan +8 °C, MGB −2 °C, UniQ +6 °C). Beacon stem arms are stripped before dimer comparison, so a designed stem is not scored as a defect. |
| **castPCR** **[G]** | `castpcr.html` | Competitive allele-specific TaqMan: an allele-specific primer, a blocker that suppresses the other allele, an allele-specific probe and a common counter-primer, designed in both orientations. For rare somatic variants in cfDNA/ctDNA and FFPE material. |
| **SNaPshot (SBE)** **[G]** | `snapshot.html` | Single-base-extension primers for multiplex SNP and InDel genotyping on capillary electrophoresis, with automatic 5′ length staggering (a poly(T) tail by default, or any non-complementary tail you supply) so several loci resolve in one trace. Also supports BAC fingerprinting and CpG methylation analysis. |
| **MLPA / digitalMLPA** **[G]** | `mlpa.html` | Ligation-dependent half-probe pairs (LPO + RPO) with the ligation junction placed on the variant, one left probe per allele up to four, a shared right probe, and length-tuning stuffers. Exchanging the universal tags for sequencing adapters converts the same arms to the digitalMLPA format. |
| **STR / SSR** **[G]** | `str.html` | Flanking primers for microsatellite fragment analysis, with the M13(−21) universal tail for economical fluorescent labelling and a reverse PIG-tail for complete 3′ adenylation. Panels are packed over four dyes subject to a minimum size gap, and loci that cannot be separated are listed rather than dropped silently. |
| **TSDR probes** **[G]** | `tsdr.html` | Toehold-exchange probes for enzyme-free allele discrimination. Toehold and clamp lengths are searched until the standard reaction free energy approaches zero, using measured internal single-mismatch parameters alongside the Watson–Crick set, and the predicted equilibrium yield is reported per allele. |
| **Insertion & knock-in genotyping** **[G]** | `genotype.html` | A paired knock-in assay (one primer inside the cassette, one in the flanking genome) and empty-site assay (both primers in the flanks) designed from one context sequence, so that together they report zygosity. Both assays are checked in silico against the two alleles supplied. |
| **loxP / FRT / rox assays** **[G]** | `loxp.html` | Locates Cre, Flp and Dre sites (loxP, lox2272, lox511, loxN, FRT, F3, rox), resolves each into arms and spacer, predicts excision or inversion from site orientation, writes out the resulting alleles, and places primers with respect to the site hairpin. |
| **Universal PCR from an alignment** **[G]** | `upcr.html` | Consensus primers for targets that exist as divergent variants. Sequences are aligned in the browser, a per-column conservation profile is computed, and primer windows are accepted only where the 3′ end is conserved in every input; 5′ columns may vary and are encoded as IUPAC degenerate bases. |

### Amplification and panel design

| Tool | Page | What it does |
|---|---|---|
| **PCR and multiplex PCR** | `pcr.html` | Standard, inverse, multiplex and bisulphite primer design, with a preset for recombinase polymerase amplification. Fluorescent STR panels for dosage analysis (QF-PCR) are designed on the STR / SSR page. |
| **Multiplex tiling panel** | `panel.html` | Custom overlapping amplicon panels for targeted sequencing, built as two interleaved pools with pooled cross-dimer checking and a coverage report. |
| **LAMP** | `lamp.html` | Four to six primers recognising six to eight target regions, with a thermodynamic dimer screen over the whole assembled set and bisulphite support. |
| **In silico PCR (ePCR)** **[G]** | `epcr.html` | Mismatch-aware virtual PCR for primers, probes, miRNA or gRNA: both strands and both bisulphite-converted strands are searched, sites are paired into amplicons and the annealing temperature is reported. |

### CRISPR

| Tool | Page | What it does |
|---|---|---|
| **Guide RNA design** | `crispr.html` | Eleven nuclease variants with their PAMs, a transparent composite activity heuristic, MIT/Hsu position-weighted off-target scoring against sequences you supply, HDR/ssODN donors with silent PAM protection, prime-editing pegRNAs with PE3 nicking, and base-editing windows with bystander prediction. |
| **CRISPR diagnostics** | `crisprdx.html` | SHERLOCK and DETECTR assay design: crRNAs for Cas12a, Cas12b and Cas13a, RPA pre-amplification primers, optional T7 promoter, collateral reporters, and allele-specific synthetic-mismatch crRNAs for genotyping. |

### Assembly and synthesis

| Tool | Page | What it does |
|---|---|---|
| **Gibson assembly** | `assembly.html` | Junction primers for seamless multi-fragment isothermal assembly, Tm-matched across overlaps and dimer-screened. |
| **PCA gene assembly** | `pca.html` | Polymerase cycling assembly: the overlapping oligonucleotide set that builds a gene from scratch, with junctions placed for an even overlap Tm, a misannealing check, end primers, a cycling programme, and block splitting for long or repetitive constructs. |
| **Gator GeneSwift probes** | `gator.html` | Fluorescein and biotin probe pairs for AAV titre determination by biolayer interferometry. |

### Oligonucleotide analysis

| Tool | Page | What it does |
|---|---|---|
| **PrimerAnalyser** | `primer.html` | Single-oligonucleotide analysis: Tm by two models with the parameter set named, GC content, molecular weight, extinction coefficient, optical density, linguistic complexity, and self-dimers drawn as alignments with Tm and ΔG. |
| **PrimersList** | `primers.html` | Batch analysis with all-against-all cross-dimer screening. |
| **Secondary structure** | `fold.html` | Minimum-free-energy folding of one sequence or hybridisation of two, with suboptimal structures, a loop-by-loop energy breakdown, and radial, arc and linear drawings exportable as SVG. Full Turner and Mathews parameter sets at any temperature. |

### Sequences, repeats and phylogeny

| Tool | Page | What it does |
|---|---|---|
| **TotalRepeats** | `repeats.html` | De novo identification, masking and clustering of interspersed and tandem repeats at genomic scale, with a zoomable repeat map and TSV, GFF3, BED, FASTA, PNG and SVG export. |
| **MUSCLE-JS** | `msa.html` | Needleman–Wunsch, Smith–Waterman and progressive multiple alignment, with automatic orientation correction and CLUSTAL or FASTA output. |
| **Phylogenetic tree** | `phylo.html` | Neighbour-joining or UPGMA from an alignment or a Newick string, with corrected distances, bootstrap support and an editable on-screen tree exportable as SVG, PNG, Newick or NEXUS. |
| **CAP-JS assembly** | `cap.html` | CAP3-style overlap–layout–consensus assembly of short, repeated or imperfect fragments, calling polymorphic columns as IUPAC codes, parallelised across Web Workers. |
| **Sequence utilities** | `seqtools.html` | Reverse complement, case and bisulphite transforms, clustering, conversion between FASTA, FASTQ, GenBank, EMBL, Clustal, PHYLIP, NEXUS and GFF3, feature and CDS extraction with translation, and FASTQ quality control with adapter screening. |

### Calculators

| Tool | Page | What it does |
|---|---|---|
| **Reaction setup** | `mix.html` | Reagent volumes for PCR, qPCR and LAMP, with presets for more than ten polymerases. |
| **Dilution and mixing** | `uni.html` | Two-solution mixing in molar, percentage, pH or any unit: supply four of six variables and the rest are computed. |

---

## How the shared engine works

Every design page loads the same core, so oligonucleotides designed for different chemistries are
scored under identical criteria and can be combined in a single panel.

- **Parsing** (`js/readseq.js`) accepts FASTA, tab- or space-delimited records, or unformatted text,
  and never throws on malformed input.
- **Validation** (`js/seq-validate.js`) runs before the engine and reports problems in place: empty
  or too-short input, protein pasted in place of nucleotide sequence, unbalanced variant brackets, a
  missing variant marker, and invalid bases listed with their counts. Errors block the run and say
  what to change; warnings allow it to proceed.
- **Repeat masking** (`js/dna.js`) uses a two-bit rolling 15-mer hash counted on both strands with a
  difference array, so cost is linear in sequence length.
- **Complexity** (`js/dna.js`) scores linguistic complexity with bit-set *n*-gram vocabulary counters
  over the four-letter and purine/pyrimidine alphabets; the default acceptance threshold is 70% on
  the design pages and 65% on the STR page. A second detector
  (`js/LowComplexitySequence.js`) merges simple-sequence and telomere-like windows into blocks that
  are excluded as priming sites.
- **Melting temperature** (`js/DNAtm.js`) uses one nearest-neighbour model throughout: Allawi and
  SantaLucia for DNA, Xia and colleagues for RNA, McTigue increments for locked residues, Watkins and
  SantaLucia for inosine, with Owczarzy monovalent and Mg²⁺ corrections. Site-wide conditions are
  55 mM monovalent cation, 1 mM Mg²⁺ and 0.2 µM oligonucleotide. On the design pages a template
  may be written as DNA or RNA: uracil is read as thymine, and inosine and locked residues as their
  parent bases. The RNA, locked-nucleic-acid and inosine parameters apply when an oligonucleotide
  carrying those residues is analysed directly in PrimerAnalyser.
- **Dimers** (`js/dimers.js`, `js/melting2.js`) are found by seeded, mismatch-tolerant
  reverse-complement alignment and re-scored as duplexes, so the output gives the structure, its Tm
  and ΔG, and whether the 3′ end is involved, rather than an unscaled flag.
- **Panel growth**: each accepted oligonucleotide is checked against the pool already accepted, so
  panels of arbitrary size are assembled incrementally.

---

## Citation

If DigitalGens contributed to your work, please cite the platform paper. The other publications
describe individual methods and applications implemented here.

### Software

- Kalendar R. (2026) DigitalGens, version 1.0. https://github.com/rkalendar/PCRtools.
  See `CITATION.cff`; the archived release DOI is listed on the repository page.

### Platform

- Kalendar R. (2025) Comprehensive web-based platform for advanced PCR design, genotyping, synthetic
  biology, molecular diagnostics, and sequence analysis. *Molecular Therapy Nucleic Acids*
  **36**(4):102716.
  [doi:10.1016/j.omtn.2025.102716](https://doi.org/10.1016/j.omtn.2025.102716) · PMID 41069987
- Kalendar R., Khassenov B., Ramankulov Y., Samuilova O., Ivanov K.I. (2017) FastPCR: an in silico
  tool for fast primer and probe design and advanced sequence analysis. *Genomics* **109**:312–319.
  [doi:10.1016/j.ygeno.2017.05.005](https://doi.org/10.1016/j.ygeno.2017.05.005) · PMID 28502701
- Kalendar R., Lee D., Schulman A.H. (2011) Java web tools for PCR, in silico PCR, and
  oligonucleotide assembly and analysis. *Genomics* **98**:137–144.
  [doi:10.1016/j.ygeno.2011.04.009](https://doi.org/10.1016/j.ygeno.2011.04.009) · PMID 21569836
- Kalendar R. (2022) A guide to using FastPCR software for PCR, in silico PCR, and oligonucleotide
  analysis. *Methods in Molecular Biology* **2392**:223–243.
  [doi:10.1007/978-1-0716-1799-1_16](https://doi.org/10.1007/978-1-0716-1799-1_16) · PMID 34773626

### Methods implemented in specific modules

- **In silico PCR** — Kalendar R., Shevtsov A., Otarbay Z., Ismailova A. (2024) In silico PCR
  analysis: a comprehensive bioinformatics tool for enhancing nucleic acid amplification assays.
  *Frontiers in Bioinformatics* **4**:1464197.
  [doi:10.3389/fbinf.2024.1464197](https://doi.org/10.3389/fbinf.2024.1464197) · PMID 39435190
- **TotalRepeats** — Kalendar R., Kairov U. (2024) Genome-wide tool for sensitive de novo
  identification and visualisation of interspersed and tandem repeats. *Bioinformatics and Biology
  Insights* **18**:11779322241306391.
  [doi:10.1177/11779322241306391](https://doi.org/10.1177/11779322241306391) · PMID 39703748
- **Allele-specific PCR** — Kalendar R., Shustov A.V., Akhmetollayev I., Kairov U. (2022) Designing
  allele-specific competitive-extension PCR-based assays for high-throughput genotyping and gene
  characterization. *Frontiers in Molecular Biosciences* **9**:773956.
  [doi:10.3389/fmolb.2022.773956](https://doi.org/10.3389/fmolb.2022.773956) · PMID 35300118
- **Allele-specific qPCR** — Kalendar R., Baidyussen A., Serikbay D., *et al.* (2021) Modified
  "allele-specific qPCR" method for SNP genotyping based on FRET. *Frontiers in Plant Science*
  **12**:747886.
  [doi:10.3389/fpls.2021.747886](https://doi.org/10.3389/fpls.2021.747886) · PMID 35082803
- **KASP allele sizing** — Kalendar R., Shevtsov V., Ismailova A. (2025) An assay for variable
  fragment length allele-specific genotyping. *Methods in Molecular Biology* **2943**:1–17.
  [doi:10.1007/978-1-0716-4642-7_1](https://doi.org/10.1007/978-1-0716-4642-7_1) · PMID 40580281
- **Repeat-based fingerprinting** — Kalendar R., Schulman A.H. (2006) IRAP and REMAP for
  retrotransposon-based genotyping and fingerprinting. *Nature Protocols* **1**:2478–2484.
  [doi:10.1038/nprot.2006.377](https://doi.org/10.1038/nprot.2006.377) · PMID 17406494

### Applications

- **Tiling-PCR panels** — Kalendar R., Kairov U., Karabayev D., *et al.* (2023) Universal
  whole-genome Oxford nanopore sequencing of SARS-CoV-2 using tiled amplicons. *Scientific Reports*
  **13**:10334.
  [doi:10.1038/s41598-023-37588-x](https://doi.org/10.1038/s41598-023-37588-x) · PMID 37365249
- **Ultra-multiplex panels** — Kalendar R., Daniyarov A., Serikzhan A., Kairov U., Basharova D.,
  Toleuzhanova A., Bekbayeva A., Romanova A., Zholdybayeva E. (2026) PCRpanel: automated design of
  ultra-multiplex PCR panels for targeted amplicon sequencing. *BMC Bioinformatics* **27**.
  [doi:10.1186/s12859-026-06648-0](https://doi.org/10.1186/s12859-026-06648-0) · published 10
  September 2026; not yet indexed in PubMed at the time of writing

---

## Licence and third-party components

DigitalGens is free and open-source software released under the **GNU General Public License
version 3** (see [`LICENSE.txt`](LICENSE.txt)), which permits use, study, modification and
redistribution by anyone, including commercially.

The website is free and open to all users and there is no login requirement.

### Vendored third-party packages

| Package | Version | Licence | Location |
|---|---|---|---|
| Bootstrap | 5.3.3 | MIT | `assets/vendor/bootstrap.min.css`, `assets/vendor/bootstrap.bundle.min.js` |
| Phosphor Icons | 2.1.1 | MIT — [`assets/vendor/phosphor/LICENSE`](assets/vendor/phosphor/LICENSE) | `assets/vendor/phosphor/` |

Both are self-hosted so that no visitor request reaches a third-party CDN. Versions, licences and
the one local modification are documented in [`assets/vendor/README.md`](assets/vendor/README.md).

### Published parameter sets embedded in the source

These are third-party scientific data re-encoded in otherwise original code, used under their
original citations:

- Turner 2004 (RNA) and Mathews/SantaLucia 2004 (DNA) folding parameters — `js/fold-params.js`
- DNA_NN3 and DNA_IMM1 nearest-neighbour tables as compiled by Biopython — `js/tsdr-nn.js`
- BLOSUM62 — `js/muscle.js`
- Hsu *et al.* 2013 off-target position weights — `js/crispr.js`
- Allawi and SantaLucia, McTigue, Watkins, Xia and Owczarzy thermodynamic parameters — `js/DNAtm.js`

Algorithms re-implemented independently, with their originals acknowledged in the source headers:
MUSCLE (`js/muscle.js`), CAP3 (`js/cap.js`), ViennaRNA recursions and dangling-end treatment
(`js/fold.js`), and neighbour-joining and UPGMA (`js/phylo.js`).

If you modify or extend this code, please include:

> Based on the work by Ruslan Kalendar — https://github.com/rkalendar/

---

## Contact

**Ruslan Kalendar**
Helsinki Institute of Life Science (HiLIFE), University of Helsinki, Finland
Email: ruslan.kalendar@helsinki.fi · [ORCID 0000-0003-3986-2460](https://orcid.org/0000-0003-3986-2460)

Bug reports are handled fastest through [GitHub Issues](https://github.com/rkalendar/PCRtools/issues).

---

## Browser support

Tested in Chrome, Firefox, Edge and Safari. The code uses no engine-specific API: no WebAssembly, no
`OffscreenCanvas`, no regular-expression lookbehind. Three tools (`fold`, `repeats`, `cap`) use Web
Workers, with a same-thread fallback. Any browser from roughly 2020 onward will run the whole site.
