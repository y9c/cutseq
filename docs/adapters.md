# 📚 Adapter Schemes

CutSeq supports various built-in adapter schemes for different NGS library preparation methods. Each scheme follows a general pattern:

## Components

- **p5**: 5' Illumina sequencing adapter (shown in <span style="background-color: #A8E6CF; padding: 5px;">light green</span>)
- **p7**: 3' Illumina sequencing adapter (shown in <span style="background-color: #D1E8D1; padding: 5px;">pale green</span>)
- **inline5**: Fixed DNA barcode sequences in brackets () (shown in <span style="background-color: #B2EBF2; padding: 5px;">light cyan</span>)
- **inline3**: Fixed DNA barcode sequences in brackets () (shown in <span style="background-color: #81D4FA; padding: 5px;">light sky blue</span>)
- **umi5**: 5' Random UMI sequences marked as N (shown in <span style="background-color: #1E90FF; padding: 5px;">dodger blue</span>)
- **umi3**: 3' Random UMI sequences marked as N (shown in <span style="background-color: #4682B4; padding: 5px;">steel blue</span>)
- **mask5**: Sequences to be masked marked as X (shown in <span style="background-color: #DCDCDC; padding: 5px;">gainsboro</span>)
- **mask3**: Sequences to be masked marked as X (shown in <span style="background-color: #D3D3D3; padding: 5px;">light gray</span>)
- **strand**: Direction indicator (shown with <span style="background-color: #FF6F61; padding: 5px; width: 9ch; display: inline-block;">&gt;</span>, <span style="background-color: #FF6F61; padding: 5px; width: 9ch; display: inline-block;">&lt;</span>, or <span style="background-color: #FF6F61; padding: 5px; width: 9ch; display: inline-block;">-</span>)

## Built-in Schemes

### DSLIGATION (dsDNA Ligation)

<div class="adapter-scheme">
<div style="display: flex; align-items: center; font-family: monospace; margin: 20px 0;">
  <span style="background-color: #A8E6CF; padding: 5px;">AGTTCTACAGTCCGACGATC</span>
  <div style="position: relative; width: 30px; height: 30px;">
    <div style="background-color: #FF6F61; width: 100%; height: 100%; position: absolute; clip-path: polygon(0% 0%, 100% 0%, 100% 100%, 0% 100%);"></div>
    <span style="position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%);">></span>
  </div>
  <span style="background-color: #D1E8D1; padding: 5px;">AGATCGGAAGAGCACACGTC</span>
</div>
</div>

- Basic dsDNA ligation with A-tailing
- Forward orientation
- No UMIs or special trimming needed

### SMALLRNA (Small RNA Libraries)

<div class="adapter-scheme">
<div style="display: flex; align-items: center; font-family: monospace; margin: 20px 0;">
  <span style="background-color: #A8E6CF; padding: 5px;">AGTTCTACAGTCCGACGATC</span>
  <div style="position: relative; width: 30px; height: 30px;">
    <div style="background-color: #FF6F61; width: 100%; height: 100%; position: absolute; clip-path: polygon(0% 0%, 100% 0%, 100% 100%, 0% 100%);"></div>
    <span style="position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%);">></span>
  </div>
  <span style="background-color: #D1E8D1; padding: 5px;">AGATCGGAAGAGCACACGTC</span>
</div>
</div>

- Used for small RNA sequencing
- Double ligation method
- Forward orientation
- Optional 2nt trimming on both ends for quality

### INLINE (Custom Barcoded Libraries)

<div class="adapter-scheme">
<div style="display: flex; align-items: center; font-family: monospace; margin: 20px 0;">
  <span style="background-color: #A8E6CF; padding: 5px;">AGTTCTACAGTCCGACGATC</span>
  <span style="background-color: #B2EBF2; padding: 5px;">NNNNN</span>
  <div style="position: relative; width: 30px; height: 30px;">
    <div style="background-color: #FF6F61; width: 100%; height: 100%; position: absolute; clip-path: polygon(0% 0%, 100% 0%, 100% 100%, 0% 100%);"></div>
    <span style="position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%);">></span>
  </div>
  <span style="background-color: #B2EBF2; padding: 5px;">NNNNN</span>
  <span style="background-color: #FFD700; padding: 5px;">(ATCACG)</span>
  <span style="background-color: #D1E8D1; padding: 5px;">AGATCGGAAGAGCACACGTC</span>
</div>
</div>

- Used for libraries with inline barcodes
- Dual UMI design (5nt each)
- Forward orientation
- Contains fixed barcode sequence

### TAKARAV2 (SMARTer® Stranded Protocol V2)

<div class="adapter-scheme">
<div style="display: flex; align-items: center; font-family: monospace; margin: 20px 0;">
  <span style="background-color: #A8E6CF; padding: 5px;">ACACGACGCTCTTCCGATCT</span>
  <span style="background-color: #DCDCDC; padding: 5px;">XXXXXX</span>
  <div style="position: relative; width: 30px; height: 30px;">
    <div style="background-color: #FF6F61; width: 100%; height: 100%; position: absolute; clip-path: polygon(0% 0%, 100% 0%, 100% 100%, 0% 100%);"></div>
    <span style="position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%);">&lt;</span>
  </div>
  <span style="background-color: #DCDCDC; padding: 5px;">XXXXXX</span>
  <span style="background-color: #D1E8D1; padding: 5px;">AGATCGGAAGAGCACACGTC</span>
</div>
</div>

- Handles polyC/G artifacts from random RT priming
- Extended masking for adaptase tail (up to 15bp)
- Reverse orientation
- Uses random polyC tail as pseudo-UMI

### PBAT (Post-Bisulfite Adapter Tagging)

<div class="adapter-scheme">
<div style="display: flex; align-items: center; font-family: monospace; margin: 20px 0;">
  <span style="background-color: #A8E6CF; padding: 5px;">ACACGACGCTCTTCCGATCT</span>
  <span style="background-color: #DCDCDC; padding: 5px;">XXXXXX</span>
  <div style="position: relative; width: 30px; height: 30px;">
    <div style="background-color: #FF6F61; width: 100%; height: 100%; position: absolute; clip-path: polygon(0% 0%, 100% 0%, 100% 100%, 0% 100%);"></div>
    <span style="position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%);">&lt;</span>
  </div>
  <span style="background-color: #DCDCDC; padding: 5px;">XXXXXX</span>
  <span style="background-color: #D1E8D1; padding: 5px;">AGATCGGAAGAGCACACGTC</span>
</div>
</div>

- Used for post-bisulfite DNA sequencing
- Random primer-based adapter addition
- Reverse orientation
- Symmetric masking for random tails

### NEXTERA (ATAC-seq)

<div class="adapter-scheme">
<div style="display: flex; align-items: center; font-family: monospace; margin: 20px 0;">
  <span style="background-color: #A8E6CF; padding: 5px;">AGATGTGTATAAGAGACAG</span>
  <div style="position: relative; width: 30px; height: 30px;">
    <div style="background-color: #FF6F61; width: 100%; height: 100%; position: absolute; clip-path: polygon(0% 0%, 100% 0%, 100% 100%, 0% 100%);"></div>
    <span style="position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%);">></span>
  </div>
  <span style="background-color: #D1E8D1; padding: 5px;">CTGTCTCTTATACACATCT</span>
</div>
</div>

- Used for ATAC-seq libraries
- Simple design without UMIs or barcodes
- Forward orientation
- Standard Nextera adapters

### ILLUMINARNA (Illumina Stranded RNA-Seq)

<div class="adapter-scheme">
<div style="display: flex; align-items: center; font-family: monospace; margin: 20px 0;">
  <span style="background-color: #A8E6CF; padding: 5px;">AGATGTGTATAAGAGACAG</span>
  <div style="position: relative; width: 30px; height: 30px;">
    <div style="background-color: #FF6F61; width: 100%; height: 100%; position: absolute; clip-path: polygon(0% 0%, 100% 0%, 100% 100%, 0% 100%);"></div>
    <span style="position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%);">&lt;</span>
  </div>
  <span style="background-color: #D1E8D1; padding: 5px;">CTGTCTCTTATACACATCT</span>
</div>
</div>

- Standard Illumina stranded RNA-seq protocol
- Reverse orientation
- Simple design without UMIs or masking
- Direct adapter ligation method

### STRANDED (Stranded RNA Library)

<div class="adapter-scheme">
<div style="display: flex; align-items: center; font-family: monospace; margin: 20px 0;">
  <span style="background-color: #A8E6CF; padding: 5px;">ACACGACGCTCTTCCGATCT</span>
  <span style="background-color: #DCDCDC; padding: 5px;">X</span>
  <div style="position: relative; width: 30px; height: 30px;">
    <div style="background-color: #FF6F61; width: 100%; height: 100%; position: absolute;"></div>
    <span style="position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%);">&lt;</span>
  </div>
  <span style="background-color: #DCDCDC; padding: 5px;">XXX</span>
  <span style="background-color: #D1E8D1; padding: 5px;">AGATCGGAAGAGCACACGTC</span>
</div>
</div>

- Reverse orientation
- For SMARTer-type stranded RNA-seq
- Minimal masking on 5’ and 3’ ends

### TAKARAV3 (SMARTer® Pico v3)

<div class="adapter-scheme">
<div style="display: flex; align-items: center; font-family: monospace; margin: 20px 0;">
  <span style="background-color: #A8E6CF; padding: 5px;">ACACGACGCTCTTCCGATCT</span>
  <span style="background-color: #DCDCDC; padding: 5px;">XXX</span>
  <div style="position: relative; width: 30px; height: 30px;">
    <div style="background-color: #FF6F61; width: 100%; height: 100%; position: absolute;"></div>
    <span style="position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%);">&lt;</span>
  </div>
  <span style="background-color: #DCDCDC; padding: 5px;">XXXXXX</span>
  <span style="background-color: #1E90FF; padding: 5px;">NNNNNNNN</span>
  <span style="background-color: #D1E8D1; padding: 5px;">AGATCGGAAGAGCACACGTC</span>
</div>
</div>

- Reverse orientation
- SMARTer v3 protocol with complex UMI linker
- 14-nt UMI used for molecule tracking

### ECLIP6 (eCLIP or SAC-seq with 6nt UMI)

<div class="adapter-scheme">
<div style="display: flex; align-items: center; font-family: monospace; margin: 20px 0;">
  <span style="background-color: #A8E6CF; padding: 5px;">ACACGACGCTCTTCCGATCT</span>
  <span style="background-color: #DCDCDC; padding: 5px;">XX</span>
  <div style="position: relative; width: 30px; height: 30px;">
    <div style="background-color: #FF6F61; width: 100%; height: 100%; position: absolute;"></div>
    <span style="position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%);">&lt;</span>
  </div>
  <span style="background-color: #DCDCDC; padding: 5px;">X</span>
  <span style="background-color: #1E90FF; padding: 5px;">NNNNNN</span>
  <span style="background-color: #D1E8D1; padding: 5px;">AGATCGGAAGAGCACACGTC</span>
</div>
</div>

- Reverse orientation
- 6-nt UMI for identifying PCR duplicates
- Common in SAC-seq or eCLIP protocols

### ECLIP10 (eCLIP with 10nt UMI)

<div class="adapter-scheme">
<div style="display: flex; align-items: center; font-family: monospace; margin: 20px 0;">
  <span style="background-color: #A8E6CF; padding: 5px;">ACACGACGCTCTTCCGATCT</span>
  <span style="background-color: #DCDCDC; padding: 5px;">XX</span>
  <div style="position: relative; width: 30px; height: 30px;">
    <div style="background-color: #FF6F61; width: 100%; height: 100%; position: absolute;"></div>
    <span style="position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%);">&lt;</span>
  </div>
  <span style="background-color: #DCDCDC; padding: 5px;">X</span>
  <span style="background-color: #1E90FF; padding: 5px;">NNNNNNNNNN</span>
  <span style="background-color: #D1E8D1; padding: 5px;">AGATCGGAAGAGCACACGTC</span>
</div>
</div>

- Reverse orientation
- 10-nt UMI improves molecule resolution
- Used in eCLIP or high-depth cDNA methods

### SACSEQV3 (SAC-seq with dual UMIs)

<div class="adapter-scheme">
<div style="display: flex; align-items: center; font-family: monospace; margin: 20px 0;">
  <span style="background-color: #A8E6CF; padding: 5px;">ACACGACGCTCTTCCGATCT</span>
  <span style="background-color: #1E90FF; padding: 5px;">NNNNNNNN</span>
  <span style="background-color: #DCDCDC; padding: 5px;">X</span>
  <div style="position: relative; width: 30px; height: 30px;">
    <div style="background-color: #FF6F61; width: 100%; height: 100%; position: absolute;"></div>
    <span style="position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%);">&gt;</span>
  </div>
  <span style="background-color: #DCDCDC; padding: 5px;">XX</span>
  <span style="background-color: #4682B4; padding: 5px;">NNNNNNNN</span>
  <span style="background-color: #D1E8D1; padding: 5px;">AGATCGGAAGAGCACACGTC</span>
</div>
</div>

- Dual-UMI design (8 nt each side)
- cDNA ligation-based protocols
- Used in SAC-seq v3 and similar workflows (designed by YC)

### XGENRNA (xGen RNA with polyG/C tails)

<div class="adapter-scheme">
<div style="display: flex; align-items: center; font-family: monospace; margin: 20px 0;">
  <span style="background-color: #A8E6CF; padding: 5px;">ACACGACGCTCTTCCGATCT</span>
  <span style="background-color: #DCDCDC; padding: 5px;">XXXXXX</span>
  <div style="position: relative; width: 30px; height: 30px;">
    <div style="background-color: #FF6F61; width: 100%; height: 100%; position: absolute;"></div>
    <span style="position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%);">&lt;</span>
  </div>
  <span style="background-color: #DCDCDC; padding: 5px;">XXXXXXXXXXXXXXX</span>
  <span style="background-color: #D1E8D1; padding: 5px;">AGATCGGAAGAGCACACGTC</span>
</div>
</div>

- Reverse orientation
- Random RT primer artifacts masked
- Adaptase tail trimmed (up to 15bp)

### XGENMETHY (xGen Methyl-seq)

<div class="adapter-scheme">
<div style="display: flex; align-items: center; font-family: monospace; margin: 20px 0;">
  <span style="background-color: #A8E6CF; padding: 5px;">ACACGACGCTCTTCCGATCT</span>
  <span style="background-color: #DCDCDC; padding: 5px;">XX</span>
  <div style="position: relative; width: 30px; height: 30px;">
    <div style="background-color: #FF6F61; width: 100%; height: 100%; position: absolute;"></div>
    <span style="position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%);">&gt;</span>
  </div>
  <span style="background-color: #DCDCDC; padding: 5px;">XXXXXXXXXX</span>
  <span style="background-color: #D1E8D1; padding: 5px;">AGATCGGAAGAGCACACGTC</span>
</div>
</div>

- Forward orientation
- Trim 2nt from R1 start and 10nt from R2 start
- Designed for bisulfite-converted methylation libraries

### XGENSNMC (xGen snmC-seq)

<div class="adapter-scheme">
<div style="display: flex; align-items: center; font-family: monospace; margin: 20px 0;">
  <span style="background-color: #A8E6CF; padding: 5px;">ACACGACGCTCTTCCGATCT</span>
  <span style="background-color: #DCDCDC; padding: 5px;">XXXXXX</span>
  <div style="position: relative; width: 30px; height: 30px;">
    <div style="background-color: #FF6F61; width: 100%; height: 100%; position: absolute;"></div>
    <span style="position: absolute; top: 50%; left: 50%; transform: translate(-50%, -50%);">&gt;</span>
  </div>
  <span style="background-color: #DCDCDC; padding: 5px;">XXXXXXXXXXXXXXX</span>
  <span style="background-color: #D1E8D1; padding: 5px;">AGATCGGAAGAGCACACGTC</span>
</div>
</div>

- Forward orientation
- For single-cell bisulfite sequencing
- 15bp trimming due to RT or ligation artifact

<style>
.adapter-scheme {
  background: #f8f9fa;
  border-radius: 8px;
  padding: 15px;
  margin: 15px 0;
  box-shadow: 0 2px 5px rgba(0, 0, 0, 0.1);
  border: 1px solid #dee2e6;
}
</style>
