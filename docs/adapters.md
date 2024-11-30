# Adapter Schemes

CutSeq supports various built-in adapter schemes for different NGS library preparation methods. Each scheme follows a general pattern:

## Components

- **p5**: 5' Illumina sequencing adapter (shown in <span style="background-color: #90EE90; padding: 5px;">light green</span>)
- **p7**: 3' Illumina sequencing adapter (shown in <span style="background-color: #98FB98; padding: 5px;">pale green</span>)
- **inline5**: Fixed DNA barcode sequences in brackets () (shown in <span style="background-color: #ADD8E6; padding: 5px;">light blue</span>)
- **inline3**: Fixed DNA barcode sequences in brackets () (shown in <span style="background-color: #ADD8E6; padding: 5px;">light blue</span>)
- **umi5**: 5' Random UMI sequences marked as N (shown in <span style="background-color: #4169E1; padding: 5px;">royal blue</span>)
- **umi3**: 3' Random UMI sequences marked as N (shown in <span style="background-color: #6495ED; padding: 5px;">cornflower blue</span>)
- **mask5**: Sequences to be masked marked as X (shown in <span style="background-color: #D3D3D3; padding: 5px;">light gray</span>)
- **mask3**: Sequences to be masked marked as X (shown in <span style="background-color: #D3D3D3; padding: 5px;">light gray</span>)
- **strand**: Direction indicator (shown with <span style="background-color: #FF0000; padding: 5px; width: 10ch; display: inline-block;">&gt;</span>, <span style="background-color: #FF0000; padding: 5px; width: 10ch; display: inline-block;">&lt;</span>, or <span style="background-color: #FF0000; padding: 5px; width: 10ch; display: inline-block;">-</span>)

## Built-in Schemes

### DSLIGATION (dsDNA Ligation)
<div class="adapter-scheme">
<div style="display: flex; align-items: center; font-family: monospace; margin: 20px 0;">
  <span style="background-color: #90EE90; padding: 5px;">AGTTCTACAGTCCGACGATCT</span>
  <span style="background-color: #FF0000; padding: 5px; width: 10ch; display: inline-block;">&gt;</span>
  <span style="background-color: #90EE90; padding: 5px;">AGATCGGAAGAGCACACGTC</span>
</div>
</div>

- Basic dsDNA ligation with A-tailing
- Forward orientation
- No UMIs or special trimming needed

### SMALLRNA (Small RNA Libraries)
<div class="adapter-scheme">
<div style="display: flex; align-items: center; font-family: monospace; margin: 20px 0;">
  <span style="background-color: #90EE90; padding: 5px;">CACGACGCTCTTCCGATCT</span>
  <span style="background-color: #FF0000; padding: 5px; width: 10ch; display: inline-block;">&gt;</span>
  <span style="background-color: #90EE90; padding: 5px;">AGATCGGAAGAGCACACGTC</span>
</div>
</div>

- Used for small RNA sequencing
- Double ligation method
- Forward orientation
- Optional 2nt trimming on both ends for quality

### INLINE (Custom Barcoded Libraries)
<div class="adapter-scheme">
<div style="display: flex; align-items: center; font-family: monospace; margin: 20px 0;">
  <span style="background-color: #90EE90; padding: 5px;">AGTTCTACAGTCCGACGATC</span>
  <span style="background-color: #ADD8E6; padding: 5px;">NNNNN</span>
  <span style="background-color: #FF0000; padding: 5px; width: 10ch; display: inline-block;">&gt;</span>
  <span style="background-color: #ADD8E6; padding: 5px;">NNNNN</span>
  <span style="background-color: #FFD700; padding: 5px;">(ATCACG)</span>
  <span style="background-color: #90EE90; padding: 5px;">AGATCGGAAGAGCACACGTC</span>
</div>
</div>

- Used for libraries with inline barcodes
- Dual UMI design (5nt each)
- Forward orientation
- Contains fixed barcode sequence

### TAKARAV2 (SMARTer® Stranded Protocol V2)
<div class="adapter-scheme">
<div style="display: flex; align-items: center; font-family: monospace; margin: 20px 0;">
  <span style="background-color: #90EE90; padding: 5px;">ACACGACGCTCTTCCGATCT</span>
  <span style="background-color: #D3D3D3; padding: 5px;">XXX</span>
  <span style="background-color: #FF0000; padding: 5px; width: 10ch; display: inline-block;">&lt;</span>
  <span style="background-color: #D3D3D3; padding: 5px;">XXX</span>
  <span style="background-color: #90EE90; padding: 5px;">AGATCGGAAGAGCACACGTC</span>
</div>
</div>

- Earlier version of TAKARA stranded protocol
- Includes masking for template switching artifacts
- Reverse orientation to RNA
- No UMI sequences

### STRANDED (Generic Stranded RNA-seq)
<div class="adapter-scheme">
<div style="display: flex; align-items: center; font-family: monospace; margin: 20px 0;">
  <span style="background-color: #90EE90; padding: 5px;">ACACGACGCTCTTCCGATCT</span>
  <span style="background-color: #D3D3D3; padding: 5px;">X</span>
  <span style="background-color: #FF0000; padding: 5px; width: 10ch; display: inline-block;">&lt;</span>
  <span style="background-color: #D3D3D3; padding: 5px;">XXX</span>
  <span style="background-color: #90EE90; padding: 5px;">AGATCGGAAGAGCACACGTC</span>
</div>
</div>

- Basic stranded RNA-seq protocol
- Minimal masking for ligation artifacts
- Reverse orientation
- No UMI sequences

### TAKARAV3 (SMARTer® Stranded Total RNA-Seq Kit v3)
<div class="adapter-scheme">
<div style="display: flex; align-items: center; font-family: monospace; margin: 20px 0;">
  <span style="background-color: #90EE90; padding: 5px;">ACACGACGCTCTTCCGATCT</span>
  <span style="background-color: #D3D3D3; padding: 5px;">XXX</span>
  <span style="background-color: #FF0000; padding: 5px; width: 10ch; display: inline-block;">&lt;</span>
  <span style="background-color: #D3D3D3; padding: 5px;">XXXXXX</span>
  <span style="background-color: #4169E1; padding: 5px;">NNNNNNNN</span>
  <span style="background-color: #90EE90; padding: 5px;">AGATCGGAAGAGCACACGTC</span>
</div>
</div>

- Used for stranded RNA-seq
- Contains 8nt UMI
- Reverse orientation to RNA
- Includes masking for template switching artifacts

### ECLIP6 (eCLIP Protocol)
<div class="adapter-scheme">
<div style="display: flex; align-items: center; font-family: monospace; margin: 20px 0;">
  <span style="background-color: #90EE90; padding: 5px;">ACACGACGCTCTTCCGATCT</span>
  <span style="background-color: #D3D3D3; padding: 5px;">XX</span>
  <span style="background-color: #FF0000; padding: 5px; width: 10ch; display: inline-block;">&lt;</span>
  <span style="background-color: #D3D3D3; padding: 5px;">X</span>
  <span style="background-color: #FF6B6B; padding: 5px;">NNNNNN</span>
  <span style="background-color: #90EE90; padding: 5px;">AGATCGGAAGAGCACACGTC</span>
</div>
</div>

- Used for eCLIP and similar protocols
- Contains 6nt UMI
- Reverse orientation
- Short masking regions

### ECLIP10 (Extended eCLIP Protocol)
<div class="adapter-scheme">
<div style="display: flex; align-items: center; font-family: monospace; margin: 20px 0;">
  <span style="background-color: #90EE90; padding: 5px;">ACACGACGCTCTTCCGATCT</span>
  <span style="background-color: #D3D3D3; padding: 5px;">XX</span>
  <span style="background-color: #FF0000; padding: 5px; width: 10ch; display: inline-block;">&lt;</span>
  <span style="background-color: #D3D3D3; padding: 5px;">X</span>
  <span style="background-color: #FF6B6B; padding: 5px;">NNNNNNNNNN</span>
  <span style="background-color: #90EE90; padding: 5px;">AGATCGGAAGAGCACACGTC</span>
</div>
</div>

- Extended version of eCLIP protocol
- Contains 10nt UMI for higher complexity
- Reverse orientation
- Short masking regions

### SACSEQV3 (SAC-seq Protocol V3)
<div class="adapter-scheme">
<div style="display: flex; align-items: center; font-family: monospace; margin: 20px 0;">
  <span style="background-color: #90EE90; padding: 5px;">AGTTCTACAGTCCGACGATCT</span>
  <span style="background-color: #FF6B6B; padding: 5px;">NNNNNNNN</span>
  <span style="background-color: #D3D3D3; padding: 5px;">X</span>
  <span style="background-color: #FF0000; padding: 5px; width: 10ch; display: inline-block;">&gt;</span>
  <span style="background-color: #D3D3D3; padding: 5px;">XX</span>
  <span style="background-color: #FF6B6B; padding: 5px;">NNNNNNNN</span>
  <span style="background-color: #90EE90; padding: 5px;">AGATCGGAAGAGCACACGTC</span>
</div>
</div>

- Dual UMI design (8nt each)
- Forward orientation
- Balanced masking on both sides
- Used for high-complexity libraries

### XGENRNA (xGen RNA Library Prep)
<div class="adapter-scheme">
<div style="display: flex; align-items: center; font-family: monospace; margin: 20px 0;">
  <span style="background-color: #90EE90; padding: 5px;">ACACGACGCTCTTCCGATCT</span>
  <span style="background-color: #D3D3D3; padding: 5px;">XXXXXX</span>
  <span style="background-color: #FF0000; padding: 5px; width: 10ch; display: inline-block;">&lt;</span>
  <span style="background-color: #D3D3D3; padding: 5px;">XXXXXXXXXXXXXXX</span>
  <span style="background-color: #90EE90; padding: 5px;">AGATCGGAAGAGCACACGTC</span>
</div>
</div>

- Handles polyC/G artifacts from random RT priming
- Extended masking for adaptase tail (up to 15bp)
- Reverse orientation
- Uses random polyC tail as pseudo-UMI

### XGENMETHY (xGen Methyl-Seq)
<div class="adapter-scheme">
<div style="display: flex; align-items: center; font-family: monospace; margin: 20px 0;">
  <span style="background-color: #90EE90; padding: 5px;">ACACGACGCTCTTCCGATCT</span>
  <span style="background-color: #D3D3D3; padding: 5px;">XX</span>
  <span style="background-color: #FF0000; padding: 5px; width: 10ch; display: inline-block;">&gt;</span>
  <span style="background-color: #D3D3D3; padding: 5px;">XXXXXXXXXX</span>
  <span style="background-color: #90EE90; padding: 5px;">AGATCGGAAGAGCACACGTC</span>
</div>
</div>

- Designed for methylation sequencing
- Trims 10 bases from read ends
- Forward orientation
- Includes random primer artifact removal

### XGENSNMC (snmC-seq Protocol)
<div class="adapter-scheme">
<div style="display: flex; align-items: center; font-family: monospace; margin: 20px 0;">
  <span style="background-color: #90EE90; padding: 5px;">ACACGACGCTCTTCCGATCT</span>
  <span style="background-color: #D3D3D3; padding: 5px;">XXXXXX</span>
  <span style="background-color: #FF0000; padding: 5px; width: 10ch; display: inline-block;">&gt;</span>
  <span style="background-color: #D3D3D3; padding: 5px;">XXXXXXXXXXXXXXX</span>
  <span style="background-color: #90EE90; padding: 5px;">AGATCGGAAGAGCACACGTC</span>
</div>
</div>

- Specialized for single-nucleus methylome sequencing
- Extended 15-base trimming
- Forward orientation
- Heavy masking for protocol artifacts

### PBAT (Post-Bisulfite Adapter Tagging)
<div class="adapter-scheme">
<div style="display: flex; align-items: center; font-family: monospace; margin: 20px 0;">
  <span style="background-color: #90EE90; padding: 5px;">ACACGACGCTCTTCCGATCT</span>
  <span style="background-color: #D3D3D3; padding: 5px;">XXXXXX</span>
  <span style="background-color: #FF0000; padding: 5px; width: 10ch; display: inline-block;">&lt;</span>
  <span style="background-color: #D3D3D3; padding: 5px;">XXXXXX</span>
  <span style="background-color: #90EE90; padding: 5px;">AGATCGGAAGAGCACACGTC</span>
</div>
</div>

- Used for post-bisulfite DNA sequencing
- Random primer-based adapter addition
- Reverse orientation
- Symmetric masking for random tails

### NEXTERA (ATAC-seq)
<div class="adapter-scheme">
<div style="display: flex; align-items: center; font-family: monospace; margin: 20px 0;">
  <span style="background-color: #90EE90; padding: 5px;">AGATGTGTATAAGAGACAG</span>
  <span style="background-color: #FF0000; padding: 5px; width: 10ch; display: inline-block;">&gt;</span>
  <span style="background-color: #90EE90; padding: 5px;">CTGTCTCTTATACACATCT</span>
</div>
</div>

- Used for ATAC-seq libraries
- Simple design without UMIs or barcodes
- Forward orientation
- Standard Nextera adapters

### ILLUMINARNA (Illumina Stranded RNA-Seq)
<div class="adapter-scheme">
<div style="display: flex; align-items: center; font-family: monospace; margin: 20px 0;">
  <span style="background-color: #90EE90; padding: 5px;">AGATGTGTATAAGAGACAG</span>
  <span style="background-color: #FF0000; padding: 5px; width: 10ch; display: inline-block;">&lt;</span>
  <span style="background-color: #90EE90; padding: 5px;">CTGTCTCTTATACACATCT</span>
</div>
</div>

- Standard Illumina stranded RNA-seq protocol
- Reverse orientation
- Simple design without UMIs or masking
- Direct adapter ligation method

<style>
.adapter-scheme {
  background: #f8f9fa;
  border-radius: 4px;
  padding: 10px;
  margin: 10px 0;
}
</style> 