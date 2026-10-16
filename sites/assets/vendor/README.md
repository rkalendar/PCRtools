# Third-party components

Both packages are vendored rather than loaded from a CDN, so that no visitor request reaches a
third-party host. Both are MIT-licensed, which is compatible with the GPL-3.0 under which DigitalGens
is released.

| Package | Version | Licence | Files |
|---|---|---|---|
| [Bootstrap](https://getbootstrap.com/) | 5.3.3 | MIT | `bootstrap.min.css`, `bootstrap.bundle.min.js` |
| [Phosphor Icons](https://phosphoricons.com/) | 2.1.1 (regular weight) | MIT — see [`phosphor/LICENSE`](phosphor/LICENSE) | `phosphor/style.css`, `phosphor/Phosphor.woff2` |

Bootstrap carries its copyright and licence banner inside the minified files.

Phosphor Icons is vendored from `@phosphor-icons/web`. The only change from upstream is in
`phosphor/style.css`: the `@font-face` rule lists woff2 alone, because the woff, ttf and svg
fallbacks target browsers older than 2016, and it adds `font-display: block`.

## Published parameter sets embedded elsewhere in the source

These are third-party scientific data re-encoded in otherwise original code, used under their
original citations rather than under a software licence:

| Data | Source | Where |
|---|---|---|
| RNA and DNA folding parameters | Turner 2004; Mathews / SantaLucia 2004 | `js/fold-params.js` |
| Nearest-neighbour and internal-mismatch tables | DNA_NN3 and DNA_IMM1 as compiled by Biopython | `js/tsdr-nn.js` |
| Protein substitution matrix | BLOSUM62 | `js/muscle.js` |
| CRISPR off-target position weights | Hsu *et al.* 2013, *Nat Biotechnol* 31:827–832 | `js/crispr.js` |
| Oligonucleotide thermodynamics | Allawi & SantaLucia 1997; SantaLucia 1998; Xia *et al.* 1998; McTigue *et al.* 2004; Watkins & SantaLucia 2005; Owczarzy *et al.* 2008 | `js/DNAtm.js` |

Algorithms re-implemented independently, with their originals acknowledged in the file headers:
MUSCLE (`js/muscle.js`), CAP3 (`js/cap.js`), the ViennaRNA recursions and d2 dangling-end treatment
(`js/fold.js`), and neighbour-joining and UPGMA (`js/phylo.js`).
