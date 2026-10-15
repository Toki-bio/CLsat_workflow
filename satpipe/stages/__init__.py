"""Stage registry: the contract each stage must keep (HANDOFF.md, "The pipeline as it is").

`scripts` are the published scripts under satpipe/stages/ (the units of work; they are what the
reproduction report tested). `wired` becomes True when the driver (work-plan step 5+) can call a
stage from the config. Nothing is wired yet: `satpipe run` refuses to start until it is.
"""
from __future__ import annotations

import os
from dataclasses import dataclass


@dataclass(frozen=True)
class Stage:
    id: str
    name: str
    scripts: tuple
    inputs: str
    outputs: str
    hand: str = ""          # curated input that cannot be derived automatically
    wired: bool = False
    optional: bool = False
    needs: tuple = ()       # external tools


STAGES = [
    Stage("02", "TRF scan", (), "genome FASTA per species", "<code>.dat", needs=("trf",)),
    Stage("03", "motif extraction", ("01_discovery_clustering/Step3.sh",),
          ".dat files (basename = species code)", "motifs_output/<code>.fa", needs=("gawk",)),
    Stage("04", "simple-repeat filter", ("01_discovery_clustering/Step4.sh",), "motif FASTA",
          "<code>.fa.clean.fas / .repeats.fas", needs=("trf", "gawk")),
    Stage("05", "canonical rotation", ("01_discovery_clustering/Step5.sh",), ".clean.fas", "clean/<file>.rot"),
    Stage("06", "per-species clustering", ("01_discovery_clustering/Step6.sh",), ".rot",
          "clusters_by_file/<code>/<code>_cluster_N.fa", needs=("vsearch",)),
    Stage("07", "MAFFT rotation", ("01_discovery_clustering/Step7.sh",), "cluster FASTA",
          "*_rotated/*.rotated.fa", needs=("mafft",)),
    Stage("08", "cross-species clustering", ("02_seed_expansion_scan/Step8_cross_cluster.sh",),
          "all rotated + rotated reference consensuses", "cross_clusters/clust_N, cross_species_centroids.fa",
          hand="references.rotated_consensus", needs=("vsearch",)),
    Stage("09", "family assignment", ("02_seed_expansion_scan/Step9_family_fasta.sh",),
          "clusters + family_assignment TSV", "<family>.members.fa",
          hand="family_assignment (cluster -> family table from UPGMA tree inspection)"),
    Stage("10", "seed BEDs", ("02_seed_expansion_scan/make_bed_from_clsats_headers.sh",), "members FASTA",
          "<family>.bed (BED6, 0-based)", needs=("gawk",)),
    Stage("11", "expansion + scan", ("02_seed_expansion_scan/clisat_expand_then_scan_merged.sh",),
          "seeds, genomes.tsv, consensus, MONO", "out_<family>.expanded_merged.bed, .bits.tsv, .dropped.tsv",
          needs=("ssearch36", "samtools", "bedtools", "parallel")),
    Stage("12", "MASTER", ("03_master_arbitration/build_master.sh",), "expansion BEDs of every family",
          "MASTER.merged.bed (+ provenance)", needs=("bedtools",)),
    Stage("13", "symmetric rescan", ("02_seed_expansion_scan/clisat_expand_then_scan_merged.sh",),
          "MASTER BED x every consensus (--scan-only)", "MASTER_<family>.bits.tsv",
          needs=("ssearch36", "samtools", "bedtools", "parallel")),
    Stage("14", "segmentation + arbitration", ("03_master_arbitration/arbitration_patched.sh",),
          "MASTER, bits, seeds", "MASTER.loci.segmented.bed, MASTER.final_assignment.tsv",
          needs=("gawk", "bedtools", "python3")),
    Stage("orthology", "flank orthology", ("04_flank_orthology/rebuild_ortho_segmented.py",),
          "segmented loci, genomes", "flanks.json, flank_orthologs.json, matrix",
          needs=("ssearch36", "samtools", "python3")),
    Stage("subfamily", "sub-subfamily, terminal monomers, phase",
          ("06_catalogue_finalize/subfamily_rerun/run_extract.sh",
           "05_subfamily_terminal_boundary/assign_subfamily.sh",
           "05_subfamily_terminal_boundary/terminal_monomers.sh",
           "05_subfamily_terminal_boundary/boundary_phase.sh",
           "05_subfamily_terminal_boundary/build_subfam_term_json.py"),
          "loci, subfamily consensuses, genomes", "subfamily.json, terminal_monomers.json",
          hand="subfamily.consensuses (subx/)", optional=True, needs=("ssearch36", "samtools")),
    Stage("finalisation", "large-scaffold recovery, edge refinement, re-arbitration",
          ("06_catalogue_finalize/85mb_recovery/run_chunked.sh",
           "06_catalogue_finalize/edge_refine/run_refine_edges_parallel.sh",
           "06_catalogue_finalize/rearbitration/arbitration_540.sh"),
          "catalogue", "corrected catalogue", hand="manual_loci", optional=True,
          needs=("trf", "ssearch36", "parallel")),
    Stage("post", "post-catalogue checks",
          ("07_post_catalogue/singleton_screen.sh", "07_post_catalogue/candidate_scan.py",
           "07_post_catalogue/make_fig2_from_data.py"),
          "catalogue, singletons", "C1-C7, figure 2", optional=True, needs=("ssearch36",)),
]


def scripts_dir() -> str:
    return os.path.dirname(os.path.abspath(__file__)).replace("\\", "/")


def script_path(rel: str) -> str:
    return f"{scripts_dir()}/{rel}"


def missing_scripts() -> list:
    """Scripts named in the registry that are not present in the package (should be empty)."""
    out = []
    for st in STAGES:
        for s in st.scripts:
            if not os.path.isfile(script_path(s)):
                out.append(s)
    return out
