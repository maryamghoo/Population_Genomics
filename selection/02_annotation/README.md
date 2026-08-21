# ANNOVAR Annotation Pipeline

## Overview
This workflow performs gene-based annotation of FDR-significant outlier SNPs identified by BayeScan using ANNOVAR.

The workflow prepares a custom ANNOVAR `refGene` database from the reference genome and GFF3 annotation and uses it to annotate the genomic location, associated genes, and functional consequences of outlier SNPs.

See:
`annovar_annotation.sh`

## Purpose
This workflow is used to characterize SNPs identified as candidates for selection and determine their potential functional consequences.

The workflow performs:

- Convert BayeScan FDR outlier SNPs to ANNOVAR input format
- - Match GFF3 chromosome/scaffold names to the reference genome when sequence identifiers differ
- Prepare the GFF3 annotation for conversion to GenePred format
- Convert GFF3 annotation to ANNOVAR-compatible `refGene` format
- Prepare reference genome FASTA headers
- Generate transcript sequences for the custom ANNOVAR database
- Perform gene-based annotation of FDR outlier SNPs
- Identify genomic locations, associated genes, and functional consequences of candidate SNPs


### Directory Structure

```
annotation/
│
├── bash_data/
│   ├── fdr05_outliers.txt
│   └── gff3_ref_map.txt
│
├── bash_results/
│   ├── fdr_outliers.avinput
│   ├── fsj_annovar_ready.gff3
│   ├── fsj_reference_short_headers.fna
│   └── fdr_outliers_annotated.fsj_multianno.txt
│
├── annovar_annotation.sh
└── README.md
```

The custom ANNOVAR database files generated during the workflow are stored in the ANNOVAR database directory:

```
/path/to/annovar/humandb/fsj_refGene.txt
/path/to/annovar/humandb/fsj_refGeneMrna.fa
```

### Expected Inputs
- BayeScan FDR outlier SNPs
`bash_data/fdr05_outliers.txt`. 
- chromosome/scaffold names map
`bash_data/gff3_ref_map.txt`
- Reference genome annotation file
`/path/to/genome_annotation.gff3`
- Reference genome
`/path/to/reference_genome.fna`

### Output
- Intermediate ANNOVAR input:
`bash_results/fdr_outliers.avinput`
- Prepared GFF3 annotation:
`bash_results/fsj_annovar_ready.gff3`
- Reference genome with shortened FASTA headers:
`bash_results/fsj_reference_short_headers.fna`
- Custom ANNOVAR database:
`/path/to/annovar/humandb/fsj_refGene.txt`
`/path/to/annovar/humandb/fsj_refGeneMrna.fa`
- Main annotated SNP table:
`bash_results/fdr_outliers_annotated.fsj_multianno.txt`


#### Prepare ANNOVAR input format

This step converts `fdr05_outliers.txt` file, containing FDR-significant SNPs generated from the BayeScan analysis, into ANNOVAR input format (`fdr_outliers.avinput`)

The expected column order of `fdr05_outliers.txt` is:
```
ID    CHR    Chr_Label    POS    REF    ALT
```

For example:
```
1234    CM083937.1    Chr1A    672932    C    G
1235    CM083937.1    Chr1A    701620    C    T
1236    CM083937.1    Chr1A    997649    C    A
```

The script converts columns 2, 4, 5, and 6 into ANNOVAR input format:

```
Chr    Start    End    Ref    Alt
```

For example:

```text
CM083937.1    672932    672932    C    G
CM083937.1    701620    701620    C    T
CM083937.1    997649    997649    C    A
```

#### GFF3-to-reference sequence name map

`gff3_ref_map.txt` file maps chromosome/scaffold names used in the GFF3 annotation (`/path/to/genome_annotation.gff3`) to the corresponding sequence names in the reference genome (`/path/to/reference_genome.fna`).

The `gff3_ref_map.txt` expected format is:
```
GFF3_name    Reference_name
```

For example:

```text
Chr1A    CM083937.1
Chr1B    CM083938.1
Chr2A    CM083939.1
```

The file should not contain a header.

If the sequence identifiers in the GFF3 and reference FASTA already match, this renaming step is not required and the script should be adjusted accordingly.

During this step, comma-separated values in GFF3 `Name=` attributes are also simplified by retaining the first name. This was required for the annotation used in this study and may not be necessary for all GFF3 files.

The created file at this step is called `fsj_annovar_ready.gff3` which is an interemediate file and will be used to create `fsj_refGene.txt`

#### Create refGene annotation

The prepared GFF3 file (`fsj_annovar_ready.gff3`) is converted to GenePred/refGene format using `gff3ToGenePred`.

#### Reference genome headers

The script creates a derived FASTA from the reference genome (`/path/to/reference_genome.fna`), in which text following the first space in each sequence header is removed so that sequence identifiers match those used by the annotation.

The resulting file is `fsj_reference_short_headers.fna` which is an intermediate file and is used by `/path/to/annovar/retrieve_seq_from_fasta.pl` to generate transcript sequences (`fsj_refGeneMrna.fa`) in `/path/to/annovar/humandb`:

#### Gene-based annotation

The custom ANNOVAR database consists of:

`fsj_refGene.txt` contains the GenePred/refGene annotation generated from the prepared GFF3 file.

`fsj_refGeneMrna.fa` contains transcript sequences extracted from the reference genome using the generated refGene annotation.

The FDR outlier SNPs in `fdr_outliers.avinput` are then annotated using `/path/to/annovar/table_annovar.pl` with the custom `refGene` database.

The main ANNOVAR output is:

`fdr_outliers_annotated.fsj_multianno.txt`

and contains fields including:
```
Chr
Start
End
Ref
Alt
Func.refGene
Gene.refGene
GeneDetail.refGene
ExonicFunc.refGene
AAChange.refGene
```
These fields can be used to distinguish intergenic, intronic, exonic, synonymous, nonsynonymous, and other functional classes of outlier SNPs.

### Usage
- This script was run using `ANNOVAR` (June 7, 2020 version). Ensure it is installed and available in your environment.
- `ANNOVAR` must include:
```
table_annovar.pl
retrieve_seq_from_fasta.pl
```
- The workflow also requires the UCSC `gff3ToGenePred` utility to convert the GFF3 genome annotation into GenePred/refGene format. Ensure `gff3ToGenePred` is installed and available in your environment.
- Adjust SLURM job settings (partition, CPUs, memory, runtime, output paths, etc.) based on your cluster configuration and dataset size.
- Update the following in `annovar_annotation.sh` as needed:
```
/annovar_annotation/logs/annovar_annotation.out
/annovar_annotation/logs/annovar_annotation.err

DATA_DIR="bash_data"
OUTDIR="bash_results"

GFF3="/path/to/genome_annotation.gff3"
REFERENCE_FASTA="/path/to/reference_genome.fna"

ANNOVAR="/path/to/annovar"
GFF3_TO_GENEPRED="/path/to/gff3ToGenePred"

BUILD="fsj"
```
- Submit the job:
`sbatch annovar_annotation.sh`


## Important Notes

- The reference genome and GFF3 annotation must correspond to the same genome assembly.
- Sequence names used by the GFF3 annotation must match those in the reference FASTA before ANNOVAR annotation.
- `gff3_ref_map.txt` is used to reconcile sequence-name differences between these files.
- The BayeScan outlier input used by this script is expected to have no header.
- The script does not modify the original reference genome or original GFF3 annotation.
- Prepared reference and annotation files are generated as new files.
- The custom build name (`fsj`) can be changed when adapting the workflow to another genome.
- Update all `/path/to/...` locations before submitting the job.
- Check the SLURM output and error logs after the job completes.
