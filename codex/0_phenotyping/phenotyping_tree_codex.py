#!/usr/bin/env python3
# -*- coding: UTF-8 -*-
# '''=================================================
# @Project -> File :cd8_mel_codex_phenotyping -> phenotyping_tree_asmaike_draft.py
# @Author : Gemma van der Voort
# @Time   : 16.05.24 14:45
# @Desc   : phenotyping tree for 0_initial_phenotyping.py of codex dataset
# draft phenotyping according to development version of the SPARQ-MI pipeline. Repository version used is
# available as a zip in the zenodo repository. branch: dev, commit: 6b1707ae. Phenotyping result heavily postprocessed,
# see 1_phenotyping_finetuning.py
# '''=================================================

from codex_phenotyping.CompartmentV2 import CompartmentDefinition
import matplotlib.pyplot as plt

PROJECT_ROOT = Path(__file__).resolve()
while PROJECT_ROOT.name != "ACT_melanoma_neutrophil" and PROJECT_ROOT.parent != PROJECT_ROOT:
    PROJECT_ROOT = PROJECT_ROOT.parent
if PROJECT_ROOT.name != "ACT_melanoma_neutrophil":
    raise RuntimeError("ACT_melanoma_neutrophil root not found")
sys.path.append(str(PROJECT_ROOT / "helper_files"))
from paths_parameters import fig_dir

modifier = 'phenotyping_tree_v6'

tree = CompartmentDefinition("Unclassified", [])
non_immune = tree.add_sub_compartment("Non Immune", [], negative=["CD45"])
myeloid = tree.add_sub_compartment("Myeloid", ["CD45", "CD11b"],
                                   negative=["CD31", "CD3", "B220", "NKp46"])
lymphoid = tree.add_sub_compartment("Lymphoid", ["CD45"], negative=["CD31", "CD11b", "CD11c"])

endo = non_immune.add_sub_compartment("Endothelial", ["CD31"])
frc = non_immune.add_sub_compartment("FRC", ["ERTR7", "aSMA"], negative=["CD21-35", "CD31"])
fdc = non_immune.add_sub_compartment("FDC", ["CD21-35"], negative=["aSMA", "CD31", "B220"])

mono = myeloid.add_sub_compartment("Monocyte", ["Ly6C"], negative=["F4-80", "Ly6G"])
cdc = myeloid.add_sub_compartment("cDC2", ["CD11c", "MHCII"])
rpm = myeloid.add_sub_compartment("RPM", ["F4-80"])  # not found.
cd169m = myeloid.add_sub_compartment("CD169+ Macrophage", ["CD169"], negative=["F4-80"])  # also ssm/mzm (ln/spleen)
msm = myeloid.add_sub_compartment("MSM", ["CD169", "F4-80"])  # generic slo macrophage, not found
neutro = myeloid.add_sub_compartment("Neutrophil", ["Ly6G", "Ly6C"])

nk = lymphoid.add_sub_compartment("NK", ["NKp46"], negative=["CD3", "B220"])  # check compartment!

b_cell = lymphoid.add_sub_compartment("B-cell", ["B220"], negative=["CD3", "NKp46"])
t_cell = lymphoid.add_sub_compartment("T-cell", ["CD3"], negative=["B220", "NKp46"])

# Level 3
cdc.add_sub_compartment("cDC1", ["CD103"])
# cdc.add_sub_compartment("cDC2", [])
nk.add_sub_compartment("Mature NK", ["CD11b"])
b_cell.add_sub_compartment("Mature B-cell", ["CD21-35"])

cd4 = t_cell.add_sub_compartment("CD4", ["CD4"], negative=["CD8"])
cd8_endo = t_cell.add_sub_compartment("endo CD8", ["CD8"], negative=["CD4", "CD90.1"])
cd8_pmel = t_cell.add_sub_compartment("trans CD8", ["CD8", "CD90.1"], negative=["CD4"])
apc_neutro = neutro.add_sub_compartment("APC Neutrophil", ["MHCII"])
# Level 4
cd4.add_sub_compartment("Treg", ["FoxP3"])

cd8_endo.add_sub_compartment("endo T cycling", ["Ki67"])  # adapted: proliferating
cd8_endo.add_sub_compartment("endo T exhausted", ["Tox"], negative=["Ki67"])

cd8_pmel.add_sub_compartment("trans T cycling", ["Ki67"])  # adapted
cd8_pmel.add_sub_compartment("trans T exhausted", ["Tox"], negative=["Ki67"])

if __name__ == "__main__":
    height = 10
    aspect = 5
    width = height * aspect
    plt.figure(figsize=(width, height))
    tree.draw_tree(edge_kwargs={"node_size": 5000}, label_kwargs={"clip_on": False})
    plt.tight_layout()
    plt.savefig(f"{fig_dir}/tree_{modifier}_wide.png")
    plt.show()
