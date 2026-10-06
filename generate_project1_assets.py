"""
Script to execute the MolVis-LeadOpt workflow:
Generates analog data, computes properties, and creates publication-grade plots.
"""

import os
import sys
import json
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns

# Ensure local module import
sys.path.append(os.path.abspath('01_molecular_visualization_editor'))
from molvis_toolkit import Molecule, LigandEditor, PluginRegistry

os.makedirs('01_molecular_visualization_editor/figures', exist_ok=True)
os.makedirs('01_molecular_visualization_editor/data', exist_ok=True)
os.makedirs('outputs/plots', exist_ok=True)

# 1. Define parent scaffold and modifications
# Erlotinib scaffold: 4-(3-ethynylanilino)quinazoline with C6 and C7 bis(2-methoxyethoxy) chains
scaffold_template = "C#Cc1cccc(Nc2ncnc3cc([R1])c([R2])cc23)c1"

editor = LigandEditor("Erlotinib_Quinazoline_Core", scaffold_template)
editor.register_site("R1", "C6 substitution vector pointing toward solvent channel", "OCCOC", "Solvent-exposed")
editor.register_site("R2", "C7 substitution vector adjacent to ribose pocket", "OCCOC", "Pocket-lining")

# Define systematic lead optimization analog series:
# Parent Erlotinib and 5 deliberate medicinal chemistry modifications
analogs_def = [
    ("Erlotinib (Parent Lead)", "OCCOC", {
        "mw": 393.44, "logp": 3.41, "tpsa": 74.73, "hbd": 1, "hba": 7, "rotb": 10
    }),
    ("Analog-1 (C6-Fluoro, des-PEG)", "F", {
        "mw": 337.35, "logp": 3.65, "tpsa": 56.27, "hbd": 1, "hba": 5, "rotb": 5
    }),
    ("Analog-2 (C6-Morpholino Solubilized)", "N1CCOCC1", {
        "mw": 404.47, "logp": 2.85, "tpsa": 71.41, "hbd": 1, "hba": 7, "rotb": 6
    }),
    ("Analog-3 (C6-N-Methylpiperazine)", "N1CCN(C)CC1", {
        "mw": 417.51, "logp": 2.45, "tpsa": 66.83, "hbd": 1, "hba": 7, "rotb": 6
    }),
    ("Analog-4 (C6-Methoxy Bioisostere / Gefitinib-like)", "OC", {
        "mw": 349.39, "logp": 3.20, "tpsa": 62.43, "hbd": 1, "hba": 6, "rotb": 6
    }),
    ("Analog-5 (C6-Tetrazole Acid Isostere)", "c1nnn[nH]1", {
        "mw": 387.39, "logp": 1.95, "tpsa": 105.14, "hbd": 2, "hba": 8, "rotb": 5
    }),
    ("Lead-Opt-06 (Optimized MPO Candidate)", "OCC1CCN(C)CC1", {
        "mw": 431.54, "logp": 2.70, "tpsa": 69.95, "hbd": 1, "hba": 7, "rotb": 7
    })
]

analogs = editor.generate_bioisosteres("R1", analogs_def)

# Score all analogs with plugin
scorer = PluginRegistry.get("KinaseHingeBinder_v1.0")
summary_data = []

for mol in analogs:
    p = mol.properties
    scores = scorer.score(mol)
    rec = {
        'Molecule': mol.name,
        'MW (Da)': p.molecular_weight,
        'cLogP': p.logp,
        'TPSA (Å²)': p.tpsa,
        'HBD': p.hbd,
        'HBA': p.hba,
        'RotB': p.rotatable_bonds,
        'QED': p.qed_score,
        'SAScore': p.sascore,
        'Ro5_Violations': p.lipinski_violations,
        'Pfizer_3_75': 'Pass' if not p.pfizer_3_75_alert else 'Alert',
        'MPO_Score': p.mpo_score,
        'Est_DeltaG_kcal_mol': scores['estimated_delta_g_kcal_mol'],
        'Pred_Ki_nM': scores['predicted_ki_nm'],
        'Ligand_Efficiency': scores['ligand_efficiency']
    }
    summary_data.append(rec)

df = pd.DataFrame(summary_data)
df.to_csv('01_molecular_visualization_editor/data/lead_optimization_analogs.csv', index=False)
with open('01_molecular_visualization_editor/data/lead_optimization_analogs.json', 'w') as f:
    json.dump(summary_data, f, indent=2)

print("Saved analog dataset (7 compounds).")

# ============================================================================
# Plot 1: Radar / Spider Plot for Multi-Parameter Lead Optimization Trajectory
# ============================================================================
categories = ['MPO Score', 'QED Drug-likeness', 'Affinity (norm)', 'Ligand Eff.', 'Permeability (LogP)', 'Solubility (TPSA)']
num_vars = len(categories)

# Select 3 representative molecules: Parent, C6-N-Methylpiperazine, and Lead-Opt-06
rep_mols = [df.iloc[0], df.iloc[3], df.iloc[6]]
colors = ['#ef4444', '#3b82f6', '#10b981']

angles = np.linspace(0, 2 * np.pi, num_vars, endpoint=False).tolist()
angles += angles[:1]

fig, ax = plt.subplots(figsize=(7, 7), subplot_kw=dict(polar=True), dpi=300)

for idx, (row, col) in enumerate(zip(rep_mols, colors)):
    # Normalize metrics to 0-1 scale for radar chart
    norm_aff = min(max((-row['Est_DeltaG_kcal_mol'] - 4.0) / 4.0, 0.1), 1.0)
    norm_le = min(max(row['Ligand_Efficiency'] / 0.35, 0.1), 1.0)
    norm_logp = 1.0 - abs(row['cLogP'] - 2.5) / 2.5
    norm_tpsa = 1.0 - abs(row['TPSA (Å²)'] - 70.0) / 50.0
    
    values = [row['MPO_Score'], row['QED'], norm_aff, norm_le, norm_logp, norm_tpsa]
    values += values[:1]
    
    ax.plot(angles, values, color=col, linewidth=2.2, label=row['Molecule'])
    ax.fill(angles, values, color=col, alpha=0.15)

ax.set_theta_offset(np.pi / 2)
ax.set_theta_direction(-1)
ax.set_xticks(angles[:-1])
ax.set_xticklabels(categories, fontsize=10, fontweight='semibold')
ax.set_ylim(0, 1.0)
ax.set_yticks([0.2, 0.4, 0.6, 0.8, 1.0])
ax.set_yticklabels(["0.2", "0.4", "0.6", "0.8", "1.0"], color="grey", size=8)
ax.grid(color='#cbd5e1', linestyle='--', linewidth=0.7)
ax.set_title("Multi-Parameter Lead Optimization Radar Profile\nScaffold: 4-Anilinoquinazoline", fontsize=13, fontweight='bold', pad=25)
ax.legend(loc='upper right', bbox_to_anchor=(1.35, 1.1), frameon=True, facecolor='#f8fafc', edgecolor='#94a3b8')

plt.tight_layout()
fig.savefig('01_molecular_visualization_editor/figures/lead_opt_radar_chart.png', bbox_inches='tight')
fig.savefig('outputs/plots/lead_opt_radar_chart.png', bbox_inches='tight')
plt.close(fig)

# ============================================================================
# Plot 2: Chemical Space & Lipinski Rule of 5 / Pfizer 3/75 Landscape
# ============================================================================
fig, ax = plt.subplots(figsize=(8, 6), dpi=300)

# Shaded desirable drug-like space (Lipinski & Pfizer 3/75 safety)
ax.axvspan(200, 500, color='#f0fdf4', alpha=0.6, label='Lipinski MW Space (≤500 Da)')
ax.axhspan(1.0, 3.0, color='#ecfdf5', alpha=0.6, label='Pfizer 3/75 Safety Zone (cLogP ≤ 3.0)')
ax.axvline(500, color='#dc2626', linestyle='--', linewidth=1.2, label='Ro5 MW Cutoff (500 Da)')
ax.axhline(3.0, color='#ea580c', linestyle='--', linewidth=1.2, label='Pfizer 3/75 Cutoff (cLogP = 3.0)')

# Scatter plot of analogs
scatter = ax.scatter(
    df['MW (Da)'], df['cLogP'],
    s=df['TPSA (Å²)'] * 3.5,
    c=df['QED'], cmap='viridis',
    edgecolors='black', linewidth=1.5, alpha=0.9, zorder=5
)

cbar = plt.colorbar(scatter, ax=ax)
cbar.set_label('QED Drug-Likeness Score', fontsize=11, fontweight='semibold')

for _, row in df.iterrows():
    ax.annotate(
        row['Molecule'].split(' ')[0],
        (row['MW (Da)'], row['cLogP']),
        xytext=(6, 5), textcoords='offset points',
        fontsize=9, fontweight='medium',
        bbox=dict(boxstyle='round,pad=0.2', facecolor='white', alpha=0.8, edgecolor='#cbd5e1')
    )

ax.set_xlabel('Molecular Weight (Da)', fontsize=12, fontweight='bold')
ax.set_ylabel('Calculated LogP (Wildman-Crippen)', fontsize=12, fontweight='bold')
ax.set_title("Medicinal Chemistry Chemical Space & Safety Landscape\nSize ~ TPSA | Color ~ QED Score", fontsize=13, fontweight='bold', pad=12)
ax.set_xlim(300, 480)
ax.set_ylim(1.2, 4.2)
ax.grid(True, linestyle=':', alpha=0.6)
ax.legend(loc='lower left', fontsize=9, frameon=True, facecolor='white', framealpha=0.9)

plt.tight_layout()
fig.savefig('01_molecular_visualization_editor/figures/chemical_space_distribution.png', bbox_inches='tight')
fig.savefig('outputs/plots/chemical_space_distribution.png', bbox_inches='tight')
plt.close(fig)

# ============================================================================
# Plot 3: Structure-Property Normalized Heatmap
# ============================================================================
fig, ax = plt.subplots(figsize=(9, 5), dpi=300)

heatmap_features = ['MW (Da)', 'cLogP', 'TPSA (Å²)', 'RotB', 'QED', 'SAScore', 'MPO_Score', 'Est_DeltaG_kcal_mol']
hm_df = df.set_index('Molecule')[heatmap_features]

# Z-score normalization for balanced comparison
hm_norm = (hm_df - hm_df.mean()) / hm_df.std()

sns.heatmap(hm_norm.T, cmap='coolwarm', annot=hm_df.T.round(2), fmt='g',
            cbar_kws={'label': 'Z-Score Standardized Feature'},
            linewidths=0.8, linecolor='white', ax=ax)

ax.set_title("Physicochemical & In Silico Affinity Profile Matrix", fontsize=13, fontweight='bold', pad=14)
ax.set_xlabel("")
ax.set_ylabel("Properties & Descriptors", fontsize=11, fontweight='bold')
plt.xticks(rotation=25, ha='right', fontsize=9)
plt.yticks(fontsize=10)

plt.tight_layout()
fig.savefig('01_molecular_visualization_editor/figures/structure_property_matrix.png', bbox_inches='tight')
fig.savefig('outputs/plots/structure_property_matrix.png', bbox_inches='tight')
plt.close(fig)

print("All Project 1 plots successfully generated.")
