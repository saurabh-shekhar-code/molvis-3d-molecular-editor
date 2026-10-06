# Project 1: MolVis-LeadOpt — Interactive Molecular Visualization & Ligand Editing Workbench

[![Python 3.10+](https://img.shields.io/badge/python-3.10+-blue.svg)](https://www.python.org/)
[![License: MIT](https://img.shields.io/badge/License-MIT-green.svg)](https://opensource.org/licenses/MIT)
[![3Dmol.js](https://img.shields.io/badge/WebGL-3Dmol.js-orange.svg)](https://3dmol.csb.pitt.edu/)
[![Target: EGFR (1M17)](https://img.shields.io/badge/PDB-1M17-purple.svg)](https://www.rcsb.org/structure/1M17)

A modular, extensible in silico platform for **interactive 3D protein-ligand structural visualization**, **scaffold-based R-group editing**, **bioisosteric replacement**, and **medicinal chemistry multi-parameter optimization (MPO)**.

---

## 🔬 Scientific Context
Lead optimization is the iterative refinement of hit molecules to balance target potency with ADMET properties. This workbench uses the crystallographic complex of the **Epidermal Growth Factor Receptor (EGFR)** tyrosine kinase domain with **Erlotinib** (PDB: **1M17**, 2.60 Å resolution) as an authentic clinical benchmark.

Key structural interactions analyzed:
- **Met793**: Canonical hinge region hydrogen bond to the quinazoline core.
- **Thr790**: Gatekeeper residue governing steric access to the hydrophobic pocket.
- **Lys745 / Glu762**: Catalytic salt bridge maintaining kinase active conformation.

---

## 🌟 Key Features
- **WebGL 3D Viewer (`MolecularVisualizer3D`)**: Renders protein-ligand complexes interactively inside Jupyter notebooks and standalone HTML with cartoons, sticks, atom labels, and solvent-accessible surfaces.
- **In Silico Ligand Editing (`LigandEditor`)**: Rapid R-group enumeration and bioisosteric substitution (e.g. PEG to morpholine/piperazine, carboxylic acid to tetrazole).
- **Comprehensive Property Engine (`Molecule.calculate_properties`)**:
  - Molecular Weight (MW), Wildman-Crippen LogP (cLogP), Topological Polar Surface Area (TPSA).
  - Hydrogen Bond Donors/Acceptors (HBD/HBA), Rotatable Bonds (RotB).
  - Quantitative Estimate of Drug-likeness (**QED**).
  - Synthetic Accessibility Score (**SAScore**).
  - **Lipinski Rule of 5**, **Veber criteria**, and **Pfizer 3/75 safety rule** alerts.
- **Multi-Parameter Optimization (MPO)**: Composite desirability trajectories balancing potency and ADMET profiles.
- **Extensible Plugin Interface (`BaseScoringPlugin`)**: Open-source hook allowing users to inject custom machine learning or docking scoring functions (e.g. AutoDock Vina, GNINA, RF-Score).

---

## 📁 Repository Structure
```text
01_molecular_visualization_editor/
├── Molecular_Visualization_and_Ligand_Editor.ipynb   # Executed Jupyter Notebook
├── Molecular_Visualization_and_Ligand_Editor.html    # Standalone Interactive HTML Report
├── molvis_toolkit.py                                 # Core modular Python library
├── build_notebook_and_html.py                        # Notebook & web report compiler
├── generate_project1_assets.py                       # Plotting and scoring pipeline
├── data/
│   ├── egfr_pocket_1m17.pdb                          # Extracted 1M17 kinase binding cleft
│   ├── erlotinib_aq4.sdf                             # Cleaned ligand with perceived bond orders
│   ├── lead_optimization_analogs.csv                 # Analog dataset with ADMET descriptors
│   └── lead_optimization_analogs.json                # Structured JSON analog records
└── figures/
    ├── lead_opt_radar_chart.png                      # MPO radar plot
    ├── chemical_space_distribution.png               # MW vs cLogP scatter plot
    └── structure_property_matrix.png                 # Z-score property heatmap
```

---

## 🚀 Quick Start

### 1. Launch the Standalone Web Visualizer & Ligand Editor App
Open `local_visualizer_app.html` directly in any web browser (Chrome, Edge, Firefox, Safari).
- **Interactive 3D Viewport:** Rotate, pan, zoom, toggle pocket surface, and change representations.
- **1-Click Presets & File Upload:** Load built-in EGFR-Erlotinib (1M17) or Dopamine D3 (3PBL) complexes, or drag & drop any `.pdb`, `.sdf`, or `.mol2` file.
- **Live Ligand Editing:** Substitute R-groups (C6 solvent channel, C7 ribose pocket) and observe real-time recalculation of MW, cLogP, TPSA, QED, Lipinski Ro5, and Pfizer 3/75 alerts.
- **Export Capabilities:** Download edited 3D ligand coordinate files (.pdb) and high-res PNG renders.

### 2. View the Static/WebGL Publication Report
Open `Molecular_Visualization_and_Ligand_Editor.html` in any web browser to view the executed study and interactive 3D WebGL viewer.

### 3. Run the Jupyter Notebook
```bash
jupyter notebook Molecular_Visualization_and_Ligand_Editor.ipynb
```

### 3. Extend with Custom Plugins
```python
from molvis_toolkit import BaseScoringPlugin, PluginRegistry, Molecule

class MyCustomDockingPlugin(BaseScoringPlugin):
    def name(self) -> str:
        return "CustomVinaEngine"

    def score(self, molecule: Molecule, target_pdb: str = None):
        p = molecule.properties
        # Implement custom docking or ML model scoring here:
        return {"docking_score_kcal_mol": -7.8, "confidence": 0.92}

PluginRegistry.register(MyCustomDockingPlugin())
```

---

## 📊 Sample Results

| Compound | MW (Da) | cLogP | TPSA (Å²) | QED | Pfizer 3/75 | MPO Score | Pred Ki (nM) |
| :--- | :--- | :--- | :--- | :--- | :--- | :--- | :--- |
| **Erlotinib (Parent)** | 393.44 | 3.41 | 74.73 | 0.682 | Alert | 0.685 | 2.4 |
| **Analog-2 (Morpholino)** | 404.47 | 2.85 | 71.41 | 0.768 | **Pass** | 0.782 | 1.8 |
| **Analog-3 (N-Methylpiperazine)** | 417.51 | 2.45 | 66.83 | 0.812 | **Pass** | 0.825 | 3.0 |
| **Lead-Opt-06 (Optimized)** | 431.54 | 2.70 | 69.95 | **0.842** | **Pass** | **0.865** | 2.8 |

---

## 🤝 Contributing
Contributions are welcomed! Feel free to:
1. Fork the repository.
2. Implement a new scoring plugin in `molvis_toolkit.py` (e.g., QM-ESP or machine learning affinity).
3. Submit a Pull Request.
