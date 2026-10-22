#!/bin/bash
# ---------------------------------------------------------------------------
# Coev_Seq_Backbone - worked example (AM1_1557g2, 154 aa)
#
# Run from this directory:   bash run_example.sh
# Requires: perl 5, pmrf in $PATH. Steps 4+ additionally require Rosetta.
#
# Steps 0 and 1 are reproduced here for completeness; their products
# (edge_simple.list, AM1_1557g2.pmrf) are already shipped with the example,
# so step 3A can be run on its own.
# ---------------------------------------------------------------------------
set -e

REPO="$(cd "$(dirname "$0")/../.." && pwd)"
ZCUT=1.0          # Z-score cutoff for a coevolving pair
DEPTH=100         # network depth from the binding site (100 = full component)

# --- [0] structural edge list ----------------------------------------------
perl "$REPO/scripts/edge_list_generate_simple.pl" inputs/AM1_1557g2.pdb > edge_simple.list

# --- [1] coevolution (MRF) model -------------------------------------------
pmrf build inputs/AM1_1557g2.aln --edge edge_simple.list -o AM1_1557g2.pmrf
pmrf stat  AM1_1557g2.pmrf --mode pos  > AM1_1557g2.pos
pmrf stat  AM1_1557g2.pmrf --mode pair > AM1_1557g2.pair

# --- [2] functional segmentation -------------------------------------------
perl "$REPO/scripts/pair_SMRF_to_readable_simple.pl" \
     inputs/AM1_1557g2.pdb AM1_1557g2.pair "$ZCUT" inputs/binding_site.txt np_structure

perl "$REPO/scripts/pick_2nd_shell_residues.pl" \
     "AM1_1557g2.$ZCUT.sif" inputs/binding_site.txt "$DEPTH" \
     > "design.$ZCUT.pos_all_connected"

# --- [3A] Coev-BB sequence re-design, GA route ------------------------------
# allowed_pos.txt is the pre-computed position set shipped with the example;
# to use the freshly derived one, swap in design.$ZCUT.pos_all_connected.
mkdir -p outputs
perl -I"$REPO" "$REPO/scripts/GA_coev_Opt.pl" \
     -pmrf AM1_1557g2.pmrf \
     -ref_seq inputs/AM1_1557g2.fasta \
     -allowed_pos inputs/allowed_pos.txt \
     -fixed_pos   inputs/binding_site.txt \
     -survive 0.2 -mut_rate 0.2 -N 1000 -iter 300 \
     > outputs/test_out.txt 2> outputs/test_out.trace

# optimized sequences -> one Rosetta resfile each, under ./CoevOpt_str/
perl "$REPO/scripts/diff_pos2ddgFomat.pl" outputs/test_out.txt inputs/AM1_1557g2.fasta

# --- [3B] Coev-BB resfile + NCF contigs, single-pass route ------------------
# Needs an a3m alignment; convert inputs/AM1_1557g2.aln first, e.g. with
# HH-suite's reformat.pl, then:
#
# perl "$REPO/scripts/coev_seq_backbone.pl" \
#      -a3m inputs/AM1_1557g2.a3m -pdb inputs/AM1_1557g2.pdb \
#      -pair AM1_1557g2.pair -bsite inputs/binding_site.txt \
#      -zscore 0.5 -resf_by_hotspot 1 -resf_outfn AM1_1557g2.cb.resf

# --- [4] Rosetta design jobs, filtering, candidate selection ----------------
# find . -name "*.resfile" > resfile_list.txt
# perl "$REPO/scripts/mk_model_refine_job.pl" resfile_list.txt \
#      "$REPO/protocols/fast_Design_New.xml" "$PWD/inputs/AM1_1557g2.pdb" > job_refine.txt
# # ... run the command lines in job_refine.txt on your cluster ...
#
# for p in CoevOpt_str/test_out/*.pdb; do
#     perl "$REPO/scripts/select_by_filters.pl" "$p"
# done > CoevOpt_str/test_out/filter_stat.txt
#
# perl "$REPO/scripts/get_pmrf_and_energy.pl" \
#      outputs/test_out.txt CoevOpt_str/test_out/filter_stat.txt > summary.test_out
# perl "$REPO/scripts/candidate_selection_from_summary.pl" summary.test_out > scores_candidates.txt

echo "done. optimized sequences -> outputs/test_out.txt"
