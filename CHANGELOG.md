# Changelog

## v1.0 — 15 October 2026

First versioned release of DigitalGens at https://digitalgens.org/, superseding the toolset
previously hosted at primerdigital.com/tools/ and described in Kalendar (2025),
*Mol Ther Nucleic Acids* 36:102716. This release replaces the earlier PCRtools files in this
repository with the complete DigitalGens source: 34 pages, 61 JavaScript modules and all assets.

### Tools

Twenty-nine tools, including the genotyping modules KASP / allele-specific PCR, dual-labelled qPCR
probes (TaqMan, MGB, molecular beacon, UniQ), castPCR, SNaPshot single-base extension,
MLPA / digitalMLPA, STR / SSR fragment analysis, toehold-exchange (TSDR) probes, insertion and
knock-in genotyping, loxP / FRT / rox recombinase-site assays, and universal consensus PCR from an
alignment. See `README.md` for the full list.

### Corrections made before this release

- dbSNP retrieval (`js/ncbi.js`) now keeps every alternate allele of a multi-allelic rsID. Earlier
  builds kept only the first listed alternate.
- Bisulphite mode (`js/dna.js`, `js/kasp.js`, `js/qpcr.js`, `js/castpcr.js`, `js/snapshot.js`) now
  derives both primer strands from one converted strand, so a reported pair can amplify the
  converted template.
- Uracil, inosine and locked-nucleic-acid letters in a design template are read as their parent
  bases instead of being dropped.
- Removed a handler copied from the PCR page that raised an error on every load of the SNaPshot page.

### Site

- Input validation before every design run, with blocking errors that say what to change.
- Documentation for every tool in `help.html`, with an end-to-end worked tutorial and a
  troubleshooting section linked from every tool page.
- Version and date shown in the footer of every page.
- Terms page states the GPL-3.0 licence, free access with no login, and a five-year maintenance
  commitment.
