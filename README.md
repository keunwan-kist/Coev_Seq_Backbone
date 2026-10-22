# Coev_Seq_Backbone

**Evolution-guided modular protein design: coevolved backbone (Coev-BB) sequence
re-design combined with non-coevolving fragment (NCF) backbone remodeling.**

This repository contains the scripts that implement the hybrid design protocol
outlined below. Starting from a target structure and its multiple sequence
alignment (MSA), the protocol (i) segments the protein into a *coevolved
backbone* that carries the functional/structural signal and *non-coevolving
fragments* that do not, (ii) re-designs the sequence of the coevolved backbone
under the evolutionary coupling model while keeping binding hot-spot residues
fixed, and (iii) remodels the non-coevolving fragments as flexible backbone
segments. The result is a set of idealized backbone variants of the target
enzyme.

```
  Target structure (PDB)                MSA (PSI-BLAST / HHblits)
            |                                      |
            +------------------+-------------------+
                               |
                 [0] edge_list_generate_simple.pl      structural edge list
                               |
                 [1] pmrf build / stat                 MRF model (.mrf, .pair, .pos)
                               |
       +-----------------------+------------------------+
       |                                                |
  FUNCTIONAL SEGMENTATION                               |
  [2] pair_SMRF_to_readable_simple.pl  -> coev network (.sif)
      pick_2nd_shell_residues.pl       -> hot-spot-connected positions
      coev_seq_backbone.pl             -> Coev-BB / NCF assignment
       |                                                |
       +-------------------------+----------------------+
                                 |
         +-----------------------+------------------------+
         |                                                |
   Coev-BB: SEQUENCE RE-DESIGN                   NCFs: BACKBONE REMODELING
   route A  GA_coev_Opt.pl (global MRF opt.)     coev_seq_backbone.pl
   route B  coev_seq_backbone.pl -> resfile          -> contigmap.contigs
            + Rosetta FastDesign                     -> RFdiffusion
         |                                                |
         +-----------------------+------------------------+
                                 |
                  [4] scoring, filtering, candidate selection
                  select_by_filters.pl / get_pmrf_and_energy.pl
                  candidate_selection_from_summary.pl
                                 |
                   Idealized enzyme backbone variants
```

---

## 1. Requirements

| Component | Used for | Note |
|---|---|---|
| **Perl 5** (tested on v5.16.3, x86_64-linux-thread-multi) | all scripts in `scripts/` | core modules only (`Getopt::Long`, `File::Basename`, `Storable`, `Data::Dumper`); no CPAN dependency |
| **[PMRF](https://github.com/jeongchans/pmrf)** | MRF (coevolution) model build, inference | the `pmrf` binary must be in `$PATH` |
| **PSI-BLAST** (e.g. [MPI Bioinformatics Toolkit](https://toolkit.tuebingen.mpg.de/)) or HHblits | MSA generation | the **first sequence of the alignment must be the query** |
| **Rosetta** (`rosetta_scripts`, `ref2015`) | sequence design / relaxation / filtering | route B and stage 4 |
| **RFdiffusion** | NCF backbone remodeling | optional; consumes the `contigmap.contigs` string |
| **Cytoscape** | visualization of the coevolution network | optional; reads the `.sif` / `.attrs` files |
| **PyMOL** | visualization of coevolving pairs | optional |

## 2. Repository layout

```
Coev_Seq_Backbone/
├── EA/
│   └── GA_COEV_ar.pm                     Genetic-algorithm engine (Perl module)
├── scripts/
│   ├── edge_list_generate_simple.pl      [0] PDB -> residue-pair edge list for PMRF
│   ├── pair_SMRF_to_readable_simple.pl   [2] PMRF .pair -> Cytoscape network (.sif/.attrs)
│   ├── pick_2nd_shell_residues.pl        [2] hot-spot -> N-th shell positions on the network
│   ├── high_zscore_pair_simple.pl        [2] .sif -> coevolving pair list for PyMOL
│   ├── coev_seq_backbone.pl              [2/3] Coev-BB + NCF assignment, resfile, RFdiffusion contigs
│   ├── GA_coev_Opt.pl                    [3A] GA sequence optimization under the MRF score
│   ├── diff_pos2ddgFomat.pl              [3A] GA output sequences -> per-sequence Rosetta resfiles
│   ├── mk_model_refine_job.pl            [4] resfile list -> rosetta_scripts command lines
│   ├── select_by_filters.pl              [4] Rosetta output PDB -> filter values
│   ├── get_pmrf_and_energy.pl            [4] join MRF score and Rosetta energy
│   └── candidate_selection_from_summary.pl [4] top-ranked candidate selection
├── protocols/
│   └── fast_Design_New.xml               RosettaScripts design/minimization protocol
└── examples/
    └── AM1_1557/                         worked example (AM1_1557g2)
        ├── inputs/                       aln, fasta, pdb, allowed_pos, binding_site
        ├── edge_simple.list              output of step 0
        ├── AM1_1557g2.pmrf               output of step 1
        ├── outputs/                      reference output of step 3A
        └── run_example.sh                the commands below, as a runnable script
```

`GA_coev_Opt.pl` loads `EA::GA_COEV_ar`, so it must be run with the repository
root on the include path (`perl -I<repo root> scripts/GA_coev_Opt.pl ...`).

---

## 3. Protocol

All commands below are written as if run from `examples/AM1_1557/`, with
`$REPO` pointing at the repository root.

### Step 0 — Multiple sequence alignment

Run PSI-BLAST against nr90 (3 iterations, E-value 1e-5, inclusion 1e-10,
BLOSUM62, 500 target sequences) and download the raw alignment output.
The query must be the first sequence in the file.

Example alignment: `examples/AM1_1557/inputs/AM1_1557g2.aln`.

### Step 1 — Coevolution (MRF) model

The structural edge list restricts the MRF to residue pairs that are actually in
contact in the target structure, which sharply reduces the number of spurious
couplings.

```bash
perl $REPO/scripts/edge_list_generate_simple.pl inputs/AM1_1557g2.pdb > edge_simple.list

pmrf build inputs/AM1_1557g2.aln --edge edge_simple.list -o AM1_1557g2.pmrf
pmrf stat  AM1_1557g2.pmrf --mode pos  > AM1_1557g2.pos
pmrf stat  AM1_1557g2.pmrf --mode pair > AM1_1557g2.pair
```

`edge_list_generate_simple.pl` connects two residues when any pair of their
heavy atoms other than the backbone N, C and O atoms (i.e. CA plus side-chain
atoms) is within **7.0 Å**.

### Step 2 — Functional segmentation

Convert the coupling table into a network, then read the coevolved backbone off
that network.

```bash
# PMRF .pair -> Cytoscape network, keeping pairs with Z-score >= 1.0
perl $REPO/scripts/pair_SMRF_to_readable_simple.pl \
     inputs/AM1_1557g2.pdb AM1_1557g2.pair 1.0 inputs/binding_site.txt np_structure

# positions connected to the binding hot spot within <depth> steps
perl $REPO/scripts/pick_2nd_shell_residues.pl \
     AM1_1557g2.1.0.sif inputs/binding_site.txt 100 > design.1.0.pos_all_connected
```

Outputs of the first command are `<pdb>.<zcut>.sif` (network),
`<pdb>.<zcut>.edges.attrs` (Z-scores), `<pdb>.<zcut>.edges2.attrs` (long-range
flag, |i-j| > 7) and `<pdb>.nodes.attrs` (binding-site annotation) — all
directly loadable in Cytoscape.

The fifth argument of `pair_SMRF_to_readable_simple.pl` selects whether couplings
are additionally required to be in structural contact (within 5.0 Å):
`structure` applies the filter, any other value (e.g. `np_structure`) does not.

`coev_seq_backbone.pl` performs the same segmentation internally and in one
pass, directly from an a3m alignment and the `.pair` file — see step 3B.

For visual inspection of the surviving couplings:

```bash
perl $REPO/scripts/high_zscore_pair_simple.pl -pair AM1_1557g2.1.0.sif \
     -pymol_script /path/to/coevpair_line.py
```

### Step 3A — Coev-BB sequence re-design, route A (genetic algorithm)

A genetic algorithm maximizes the MRF (coevolution) score of the sequence while
the functional constraints stay fixed. `-fixed_pos` (binding site / active site)
is never mutated and takes priority; `-allowed_pos` (the coevolved backbone
positions from step 2) is the only set of positions that may be mutated.

```bash
perl -I$REPO $REPO/scripts/GA_coev_Opt.pl \
     -pmrf AM1_1557g2.pmrf \
     -ref_seq inputs/AM1_1557g2.fasta \
     -allowed_pos inputs/allowed_pos.txt \
     -fixed_pos   inputs/binding_site.txt \
     -survive 0.2 -mut_rate 0.2 -N 1000 -iter 300 \
     > outputs/test_out.txt 2> outputs/test_out.trace
```

| Option | Meaning | Default |
|---|---|---|
| `-pmrf` | PMRF model used to score sequences | required |
| `-ref_seq` | query sequence (FASTA or plain) | required |
| `-allowed_pos` | mutatable positions, one 1-based index per line | required |
| `-fixed_pos` | positions held fixed (functional constraints) | required |
| `-N` | population size | 100 |
| `-mut_rate` | mutation probability | 0.05 |
| `-survive` | elitism ratio | 0.2 |
| `-iter` | number of generations | 100 |

Outputs: `test_out.txt` — the final population with its MRF scores;
`test_out.trace` — the best sequence and score per generation.
Reference outputs for this example are in `examples/AM1_1557/outputs/`.

Fitness evaluation shells out to `pmrf infer` once per generation, writing a
temporary `pop<pid>.fas` file in the working directory.

Each optimized sequence is then turned into a Rosetta resfile:

```bash
perl $REPO/scripts/diff_pos2ddgFomat.pl outputs/test_out.txt inputs/AM1_1557g2.fasta
# -> CoevOpt_str/test_out/test_out.<n>.resfile  (one per unique sequence)
```

### Step 3B — Coev-BB resfile and NCF remodeling, route B

`coev_seq_backbone.pl` does the segmentation and both design definitions in a
single pass. It reads the alignment in **a3m** format, so a PSI-BLAST `.aln`
has to be converted first (e.g. `reformat.pl` from HH-suite); the bundled
example ships only the `.aln`.

```bash
perl $REPO/scripts/coev_seq_backbone.pl \
     -a3m   inputs/AM1_1557g2.a3m \
     -pdb   inputs/AM1_1557g2.pdb \
     -pair  AM1_1557g2.pair \
     -bsite inputs/binding_site.txt \
     -zscore 0.5 -resf_by_hotspot 1 \
     -resf_outfn AM1_1557g2.cb.resf
```

| Option | Meaning | Default |
|---|---|---|
| `-a3m` | alignment in a3m format, query first | required |
| `-pair` | PMRF `.pair` file | required |
| `-pdb` / `-chain` | target structure / chain | — / `A` |
| `-bsite` | binding hot-spot positions | — |
| `-zscore` | Z-score cutoff for a coevolving pair | 0.5 |
| `-dist` | residue–residue contact cutoff (Å) | 7 |
| `-depth` | network depth from the hot spot defining Coev-BB | 10 |
| `-resf_by_hotspot` | 1 = hot-spot-connected backbone only, 0 = whole Coev-BB | 0 |
| `-exclude_bsite` | 1 = binding-site residues are not designed | 1 |
| `-aa_ratio` | minimum frequency in the evolutionary profile for an amino acid to be allowed | 0.1 |
| `-ncf_len` | minimum length of a non-coevolving fragment | 5 |
| `-ncf_var_std` | allowed NCF length variation factor (L/f … L×f) | 1.5 |
| `-coev_pair_seq_dist` | minimum sequence separation for a pair to be printed for PyMOL | 5 |
| `-pymol_script` | path to `coevpair_line.py`; if given, the matching `run` line is printed | — |
| `-resf_outfn` | resfile output name | `<a3m>.<pdb>.cb.resf` |

It produces:

1. **a Rosetta resfile** (`NATRO` background; `PIKAA` per designed position,
   restricted to the amino acids observed above `-aa_ratio` in the evolutionary
   profile at that position, plus the wild-type residue);
2. **the NCF segment list**, plus a PyMOL `select NCF_regions, ...` string;
3. **two RFdiffusion `contigmap.contigs` strings** — one keeping each NCF at its
   original length, one allowing the `-ncf_var_std` length range.

The resfile is applied with the RosettaScripts protocol in `protocols/`:

```bash
rosetta_scripts.linuxgccrelease \
  -parser:protocol $REPO/protocols/fast_Design_New.xml \
  -s inputs/AM1_1557g2.pdb -resfile AM1_1557g2.cb.resf \
  -ex1 -ex2 -use_input_sc -linmem_ig 10 -nstruct 1
```

`fast_Design_New.xml` runs Cartesian minimization → soft-potential repacking →
rigid-body/side-chain minimization → hard-potential repacking → minimization →
Cartesian minimization, and reports the `score_filter`, `rmsd`, `deg`, `geom`
and `cav` filter values into the output PDB.

The contig strings are passed to RFdiffusion to rebuild the NCFs around the
fixed coevolved backbone.

### Step 4 — Scoring, filtering and candidate selection

```bash
# one rosetta_scripts command line per resfile
find . -name "*.resfile" > resfile_list.txt
perl $REPO/scripts/mk_model_refine_job.pl resfile_list.txt \
     $REPO/protocols/fast_Design_New.xml inputs/AM1_1557g2.pdb > job_refine.txt

# collect filter values from the designed models
for p in CoevOpt_str/test_out/*.pdb; do
    perl $REPO/scripts/select_by_filters.pl "$p"
done > CoevOpt_str/test_out/filter_stat.txt

# join MRF score and Rosetta energy, then pick the top candidates
perl $REPO/scripts/get_pmrf_and_energy.pl \
     outputs/test_out.txt CoevOpt_str/test_out/filter_stat.txt > summary.test_out
perl $REPO/scripts/candidate_selection_from_summary.pl summary.test_out > scores_candidates.txt
```

`select_by_filters.pl` reads `cav`, `deg`, `geom` and `score_filter` from the
`#END_POSE_ENERGIES_TABLE` block of a Rosetta output PDB.
`get_pmrf_and_energy.pl` discards models with `geom == 0` (bad geometry).
`candidate_selection_from_summary.pl` reports the top **1.5 %** by MRF score.

---

## 4. Worked example

`examples/AM1_1557/` is a complete, small example for the AM1_1557g2 protein
(154 residues, 24 binding-site positions, 58 mutatable positions) and includes
the step-0 and step-1 products so that step 3A can be run directly.
`examples/AM1_1557/run_example.sh` lists the commands end to end.

---

## 5. Known limitations

These are properties of the current code, listed so that results can be
interpreted correctly.

* **Hard-coded constants.** The contact cutoff in
  `edge_list_generate_simple.pl` (7.0 Å) and the contact/pocket cutoff in
  `pair_SMRF_to_readable_simple.pl` (5.0 Å, with the long-range threshold
  |i-j| > 7) are set in the source, not exposed as options.
* **`pick_2nd_shell_residues.pl` depth.** The shell set is extended while it is
  being iterated over, so the effective expansion can exceed the requested
  `depth`. The published runs used `depth = 100`, i.e. the full connected
  component, for which this does not matter.
* **`diff_pos2ddgFomat.pl` output location** is fixed to
  `./CoevOpt_str/<prefix>/` relative to the current directory.
* **`candidate_selection_from_summary.pl`** assumes model paths of the form
  `CoevOpt_str/<prefix>/<pdb>` and a fixed 1.5 % selection ratio.
* **Route A and route B take different inputs** — route A reads the PSI-BLAST
  `.aln` and the `.mrf` model, route B reads an `a3m` alignment and the `.pair`
  table. They are not yet driven from a single configuration.
* **RFdiffusion output is not yet fed back** into the resfile-based design step;
  the two halves of the protocol are currently joined by hand.

## 6. Citation

If you use this code, please cite:

> Park, K. et al. *Optimizing protein properties with functional constraints.* (2019)

## 7. License

**To be decided.**

Note for whoever settles this: the genetic-algorithm engine in
`EA/GA_COEV_ar.pm` derives from code by Stefan Kamphausen
(<https://www.perlmonks.org/?node_id=55453>, Copyright 2001) and carries an
explicit GPL notice in its header, so anything distributed together with that
file inherits the GPL terms.

## 8. Acknowledgements

* Genetic algorithm skeleton: Stefan Kamphausen.
* [PMRF](https://github.com/jeongchans/pmrf): Jeongchan Shin et al.
