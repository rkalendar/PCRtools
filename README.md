# PCRtools

**Comprehensive Web-Based Platform for Advanced PCR Design, Molecular Diagnostics, and Sequence Analysis**

[![DOI](https://img.shields.io/badge/DOI-10.1016%2Fj.omtn.2025.102716-blue)](https://doi.org/10.1016/j.omtn.2025.102716)
[![License](https://img.shields.io/badge/License-GPLv3-blue.svg)](LICENSE.txt)
[![Language](https://img.shields.io/badge/Language-JavaScript-yellow.svg)]()


---

## 📋 Table of Contents

- [Overview](#overview)
- [Features](#features)
- [Getting Started](#getting-started)
- [Tools & Applications](#tools--applications)
  - [PCR Design](#1-pcr-multiplex-and-qf-pcr-primer-design)
  - [KASP Genotyping](#2-kasp-primers-genotyping-assay-design)
  - [LAMP Design](#3-loop-mediated-isothermal-amplification-lamp)
  - [Multiplex Tiling](#4-custom-multiplex-tiling-pcr-panel-design)
  - [Virtual PCR](#5-in-silico-pcr)
  - [Repeat Analysis](#6-totalrepeats)
  - [Gibson Assembly](#7-gibson-assembly-primer-design)
  - [Sequence Analysis](#8-primer-analysis-tools)
  - [Calculators](#9-calculators)
- [Citation](#citation)
- [Contact](#contact)
- [System Requirements](#system-requirements)

---

## 🔬 Overview

PCRtools is a comprehensive web-based platform designed for advanced molecular biology applications. It provides state-of-the-art tools for PCR primer design, genotyping, synthetic biology, molecular diagnostics, and sequence analysis. The platform supports researchers in designing robust and efficient molecular assays across various applications.

**Developed by:** Ruslan Kalendar  
**Online Platform:** [https://primerdigital.com/tools/](https://primerdigital.com/tools/)

---

## ✨ Features

- **User-Friendly Interface**: Modern web-based platform accessible from any browser
- **Comprehensive Toolkit**: 10+ specialized tools for molecular biology applications
- **Advanced Algorithms**: Optimized for specificity, sensitivity, and efficiency
- **Multiplex Support**: Design and validate complex multiplex assays
- **Real-Time Validation**: Instant feedback on primer quality and compatibility
- **Cross-Platform**: Works on Windows, macOS, Linux, and mobile devices
- **No Installation Required**: Run directly in your web browser

---

## 🚀 Getting Started

### Basic Usage

1. Clone or download this repository
2. Navigate to the `/sites/` directory
3. Open `index.html` in a modern web browser (Chrome, Firefox, Edge, or Safari)
4. Select the desired tool from the main menu
5. Follow the on-screen instructions for your specific application

### Online Access

Visit **[https://primerdigital.com/tools/](https://primerdigital.com/tools/)** to use PCRtools online without installation.

---

## 🛠️ Tools & Applications

### 1. PCR, Multiplex, and QF-PCR Primer Design

Design primers for a comprehensive range of PCR applications:

- **Standard PCR**: Conventional primer design with optimized parameters
- **Inverse PCR**: Primer design for unknown flanking sequences
- **Multiplex PCR**: Simultaneous amplification of multiple targets
- **Quantitative Fluorescence PCR (QF-PCR)**: TaqMan and MGB probe assay design
- **Bisulfite PCR**: Methylation analysis primer design
- **Real-Time qPCR**: Fluorescence probe-based multiplex assays

**Key Features:**
- SNP and InDel genotyping support
- High-throughput multiplex validation
- Primer quality scoring and optimization
- Secondary structure analysis
- Compatibility checking for multiplex reactions

---

### 2. KASP Primers Genotyping Assay Design

Design allele-specific PCR assays for genotyping applications:

- **Kompetitive Allele Specific PCR (KASP™)**
- **PCR Allele Competitive Extension (PACE™)**
- **Allele-Specific Quantitative PCR (ASQ)**

**Applications:**
- Multiallelic discrimination of SNPs
- Insertion/deletion (InDel) detection
- High-throughput genotyping
- Marker-assisted selection

**Output:**
- Optimized allele-specific forward primers
- Common reverse primer
- Complete KASP Assay Mix specifications

---

### 3. Loop-Mediated Isothermal Amplification (LAMP)

Design LAMP primer sets for rapid, isothermal DNA amplification:

**LAMP Primer Components:**
- Forward Inner Primer (FIP)
- Backward Inner Primer (BIP)
- Forward Outer Primer (F3)
- Backward Outer Primer (B3)
- Loop Primers (LF/LB) - optional

**Advantages:**
- No thermal cycler required
- Isothermal reaction (60-65°C)
- Highly specific (6-8 recognition sites)
- Rapid amplification (30-60 minutes)
- Bisulfite LAMP support for methylation analysis

---

### 4. Custom Multiplex Tiling PCR Panel Design

Design comprehensive amplicon panels for targeted sequencing:

**Target Technologies:**
- Next-Generation Sequencing (NGS)
- Third-Generation Sequencing (TGS)

**Applications:**
- Whole exome sequencing
- Targeted gene panels
- Viral genome sequencing
- Microbial diversity studies
- Plant and animal genomics

**Features:**
- Automated primer tiling across target regions
- Optimized amplicon size distribution
- Multiplex compatibility validation
- Coverage optimization
- Primer pooling strategies

---

### 5. In Silico PCR

Virtual PCR tool for computational primer validation:

**Search Capabilities:**
- Genome-wide primer binding prediction
- Off-target detection and analysis
- Multiple simultaneous target search
- Mismatch tolerance configuration
- Melting temperature prediction

**Applications:**
- PCR primer validation
- Probe specificity checking
- microRNA (miRNA) target prediction
- CRISPR guide RNA (crRNA) off-target analysis
- Multiplex primer compatibility testing

**Output:**
- Predicted PCR products
- Binding site locations
- Mismatch positions and severity
- Expected amplicon sizes

---

### 6. TotalRepeats

Comprehensive repeat sequence analysis tool:

**Detection Capabilities:**
- Interspersed repeats
- Simple sequence repeats (SSRs/microsatellites)
- Tandem repeats
- Telomeric sequences
- Satellite DNA

**Features:**
- De novo repeat identification
- Genomic-scale analysis
- Repeat masking
- Clustering and classification
- Low-complexity region detection

**Applications:**
- Genome annotation
- Evolutionary studies
- Marker development
- Population genetics

---

### 7. Gibson Assembly Primer Design

Design primers for seamless DNA assembly:

**Gibson Assembly Method:**
- Isothermal single-reaction assembly
- Multiple DNA fragment joining
- Scarless cloning
- Plasmid construction

**Applications:**
- Synthetic gene construction
- Genetic pathway assembly
- Genome engineering
- Plasmid construction
- Multi-fragment cloning

**Features:**
- Automated overlap design
- Optimal primer design for assembly junctions
- Fragment order optimization
- Assembly efficiency prediction

---

### 8. Primer Analysis Tools

#### PrimerAnalyser
Comprehensive single-sequence analysis:

**Calculated Properties:**
- Sequence length and composition
- GC content percentage
- Melting temperature (Tm)
  - Standard oligonucleotides
  - Mixed bases
  - LNA modifications
- Molecular weight
- Extinction coefficient
- Optical density (OD)
- Linguistic complexity
- Self-dimer formation potential

**Additional Features:**
- Stock solution calculator
- Dilution calculator
- Resuspension calculator

#### PrimersList
Batch analysis of multiple primers:

**Analysis Features:**
- Simultaneous analysis of primer sets
- Cross-dimer detection
- Multiplex compatibility check
- Side-by-side comparison
- Quality scoring

---

### 9. Calculators

#### PCR Reaction Setup Calculator
- PCR reaction component calculator
- LAMP reaction setup
- Custom reaction mix calculator
- Master mix preparation
- Multiple reaction scaling

#### Universal Dilution Calculator
- Molar concentration conversions
- Percentage solution mixing
- pH adjustment calculations
- Stock solution dilution
- Two-solution mixing calculator

---

## 📚 Citation

If you use PCRtools in your research, please cite:

```
Kalendar R. 2025. Comprehensive web-based platform for advanced PCR design, 
genotyping, synthetic biology, molecular diagnostics, and sequence analysis. 
Molecular Therapy Nucleic Acids, 36(4): 102716.
DOI: 10.1016/j.omtn.2025.102716
```

**Full Text Available:**
- [Cell Press](https://www.cell.com/molecular-therapy-family/nucleic-acids/fulltext/S2162-2531(25)00270-7)
- [ScienceDirect](https://www.sciencedirect.com/science/article/pii/S2162253125002707)

**If you modify or extend this code, please include:**
> Based on the work by Ruslan Kalendar — [https://github.com/rkalendar/](https://github.com/rkalendar/)

Your acknowledgment is greatly appreciated!

---

## 📧 Contact

**Author:** Ruslan Kalendar  
**Email:** ruslan.kalendar@helsinki.fi  
**GitHub:** [https://github.com/rkalendar/](https://github.com/rkalendar/)  
**Website:** [https://primerdigital.com/tools/](https://primerdigital.com/tools/)

---

## 💻 System Requirements

### Software Requirements
- **Programming Language:** JavaScript
- **Web Browser:** Any modern browser
  - Google Chrome (recommended)
  - Mozilla Firefox
  - Microsoft Edge
  - Safari

### Hardware Requirements
- **Minimum:** Any device capable of running a modern web browser
- **Recommended:** 
  - 4GB RAM or higher
  - Modern multi-core processor
  - Stable internet connection (for online version)

### Browser Compatibility
- Chrome 90+
- Firefox 88+
- Edge 90+
- Safari 14+

---

## 📄 License

This project is open source. Please see the repository for license details.

---

## 🙏 Acknowledgments

Special thanks to the molecular biology and bioinformatics communities for their valuable feedback and contributions to the development of PCRtools.

---

<div align="center">

**PCRtools** — Advancing Molecular Biology Research Through Computational Innovation

[Website](https://primerdigital.com/tools/) • [Documentation](https://primerdigital.com/tools/) • [GitHub](https://github.com/rkalendar/)

</div>
