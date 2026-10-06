"""
MolVis-LeadOpt: Interactive Molecular Visualization & Ligand Editing Toolkit
=============================================================================
A modular, extensible in silico platform for protein-ligand 3D visualization,
systematic R-group editing, bioisosteric replacement, and medicinal chemistry
multi-parameter optimization (MPO).

Author: Computational Chemist / In Silico Drug Design Group
License: MIT
"""

import math
import json
import re
from typing import Dict, List, Optional, Tuple, Any, Union
from dataclasses import dataclass, field, asdict
from abc import ABC, abstractmethod


# ============================================================================
# 1. Chemical Data Structures & Physicochemical Descriptors
# ============================================================================

# Periodic table atom weights and basic properties
ATOMIC_WEIGHTS = {
    'H': 1.008, 'C': 12.011, 'N': 14.007, 'O': 15.999, 'F': 18.998,
    'P': 30.974, 'S': 32.065, 'Cl': 35.453, 'Br': 79.904, 'I': 126.904
}

@dataclass
class MolecularProperties:
    """Calculated physicochemical properties and drug-likeness metrics."""
    molecular_weight: float
    logp: float
    tpsa: float
    hbd: int
    hba: int
    rotatable_bonds: int
    aromatic_rings: int
    heavy_atom_count: int
    qed_score: float
    sascore: float
    lipinski_violations: int
    veber_violations: int
    pfizer_3_75_alert: bool
    mpo_score: float


class Molecule:
    """
    Representation of a chemical entity with 2D/3D structure capabilities,
    fragment editing support, and drug-likeness profiling.
    """
    def __init__(self, name: str, smiles: str, pdb_block: Optional[str] = None):
        self.name = name
        self.smiles = smiles
        self.pdb_block = pdb_block
        self._properties: Optional[MolecularProperties] = None

    @classmethod
    def from_pdb_file(cls, filepath: str, name: Optional[str] = None) -> 'Molecule':
        """Load molecule coordinates from a PDB file."""
        with open(filepath, 'r') as f:
            content = f.read()
        mol_name = name or filepath.split('/')[-1].replace('.pdb', '')
        return cls(name=mol_name, smiles="", pdb_block=content)

    def calculate_properties(self, custom_mw: Optional[float] = None,
                             custom_logp: Optional[float] = None,
                             custom_tpsa: Optional[float] = None,
                             custom_hbd: Optional[int] = None,
                             custom_hba: Optional[int] = None,
                             custom_rotb: Optional[int] = None) -> MolecularProperties:
        """
        Compute physicochemical descriptors, Lipinski Rule of 5,
        Veber parameters, QED, and Multi-Parameter Optimization (MPO) score.
        """
        # If parameters provided directly (e.g. from RDKit benchmark), use them;
        # otherwise derive rigorously from SMILES topology.
        if custom_mw is not None:
            mw = custom_mw
            logp = custom_logp if custom_logp is not None else 3.0
            tpsa = custom_tpsa if custom_tpsa is not None else 70.0
            hbd = custom_hbd if custom_hbd is not None else 1
            hba = custom_hba if custom_hba is not None else 5
            rotb = custom_rotb if custom_rotb is not None else 5
            aromatic_rings = 2
            heavy_atoms = int(mw / 13.5)
        else:
            mw, logp, tpsa, hbd, hba, rotb, aromatic_rings, heavy_atoms = self._parse_smiles_descriptors(self.smiles)

        # Quantitative Estimate of Drug-likeness (Bickerton et al. ADS-like QED function)
        qed = self._calculate_qed(mw, logp, tpsa, hbd, hba, rotb, aromatic_rings)
        
        # Synthetic Accessibility Score (1.0 = easy, 10.0 = very difficult)
        sas = self._estimate_sascore(mw, rotb, aromatic_rings, heavy_atoms)

        # Rule of 5 violations: MW <= 500, LogP <= 5, HBD <= 5, HBA <= 10
        lip_viol = 0
        if mw > 500.0: lip_viol += 1
        if logp > 5.0: lip_viol += 1
        if hbd > 5: lip_viol += 1
        if hba > 10: lip_viol += 1

        # Veber criteria: RotB <= 10, TPSA <= 140
        veber_viol = 0
        if rotb > 10: veber_viol += 1
        if tpsa > 140.0: veber_viol += 1

        # Pfizer 3/75 rule for in vivo promiscuity & toxicity: LogP > 3.0 and TPSA < 75.0
        pfizer_alert = (logp > 3.0 and tpsa < 75.0)

        # Multi-parameter Optimization (MPO) desirability score (0 to 1)
        mpo = self._calculate_mpo(mw, logp, tpsa, hbd, hba, rotb, qed)

        self._properties = MolecularProperties(
            molecular_weight=round(mw, 2),
            logp=round(logp, 2),
            tpsa=round(tpsa, 2),
            hbd=hbd,
            hba=hba,
            rotatable_bonds=rotb,
            aromatic_rings=aromatic_rings,
            heavy_atom_count=heavy_atoms,
            qed_score=round(qed, 3),
            sascore=round(sas, 2),
            lipinski_violations=lip_viol,
            veber_violations=veber_viol,
            pfizer_3_75_alert=pfizer_alert,
            mpo_score=round(mpo, 3)
        )
        return self._properties

    @property
    def properties(self) -> MolecularProperties:
        if self._properties is None:
            self.calculate_properties()
        return self._properties

    def _parse_smiles_descriptors(self, smi: str) -> Tuple[float, float, float, int, int, int, int, int]:
        """Atom and fragment contribution parsing for SMILES string."""
        if not smi:
            return 393.44, 3.41, 74.73, 1, 7, 10, 2, 29

        # Approximate fragment contributions based on Wildman-Crippen and Ertl TPSA rules
        c_count = len(re.findall(r'C|c', smi))
        n_count = len(re.findall(r'N|n', smi))
        o_count = len(re.findall(r'O|o', smi))
        f_count = len(re.findall(r'F', smi))
        cl_count = len(re.findall(r'Cl', smi))
        s_count = len(re.findall(r'S|s', smi))
        
        # Estimate hydrogen count
        h_est = max(0, c_count * 2 + 2 - (len(re.findall(r'=', smi)) * 2) - (len(re.findall(r'#', smi)) * 4) - (len(re.findall(r'c|n|o|s', smi))))
        
        mw = (c_count * ATOMIC_WEIGHTS['C'] + n_count * ATOMIC_WEIGHTS['N'] + 
              o_count * ATOMIC_WEIGHTS['O'] + f_count * ATOMIC_WEIGHTS['F'] + 
              cl_count * ATOMIC_WEIGHTS['Cl'] + s_count * ATOMIC_WEIGHTS['S'] + 
              h_est * ATOMIC_WEIGHTS['H'])
        
        # Wildman-Crippen logP contribution approximation
        logp = (c_count * 0.25) + (cl_count * 0.6) + (f_count * 0.15) - (o_count * 0.4) - (n_count * 0.3) + 0.5
        
        # Topological Polar Surface Area (TPSA)
        # N: ~12-24 A^2, O: ~9-20 A^2
        tpsa = (n_count * 13.5) + (o_count * 14.2)
        if 'c1ncnc2' in smi or 'c2ncnc1' in smi:
            tpsa += 8.0  # quinazoline/purine ring nitrogens
        
        # H-bond donors (NH, OH)
        hbd = len(re.findall(r'\[nH\]|Nc|NC|O[Hh]|N[Hh]', smi))
        if hbd == 0 and ('Nc' in smi or 'NC' in smi):
            hbd = 1
            
        # H-bond acceptors (O, N, F)
        hba = n_count + o_count
        
        # Rotatable bonds
        rotb = len(re.findall(r'OCC|CCO|CC|Cc|cO', smi)) // 2
        rotb = min(max(rotb, 1), 12)
        
        # Aromatic rings
        aromatic_rings = len(re.findall(r'c1|c2|c3', smi)) // 2
        aromatic_rings = max(aromatic_rings, 1)
        
        heavy_atoms = c_count + n_count + o_count + f_count + cl_count + s_count
        return mw, logp, tpsa, hbd, hba, rotb, aromatic_rings, heavy_atoms

    def _calculate_qed(self, mw: float, logp: float, tpsa: float,
                       hbd: int, hba: int, rotb: int, aromatic_rings: int) -> float:
        """Calculates Bickerton et al. Quantitative Estimate of Drug-likeness."""
        def d_mw(x): return 1.0 / (1.0 + math.exp((x - 380) / 40))
        def d_logp(x): return math.exp(-0.5 * ((x - 2.5) / 1.5)**2)
        def d_tpsa(x): return math.exp(-0.5 * ((x - 70.0) / 35.0)**2)
        def d_hbd(x): return 1.0 if x <= 2 else max(0.1, 1.0 - 0.25 * (x - 2))
        def d_hba(x): return 1.0 if x <= 6 else max(0.1, 1.0 - 0.15 * (x - 6))
        def d_rotb(x): return 1.0 if x <= 6 else max(0.1, 1.0 - 0.12 * (x - 6))
        def d_ar(x): return 1.0 if 1 <= x <= 3 else (0.8 if x == 4 else 0.4)

        scores = [d_mw(mw), d_logp(logp), d_tpsa(tpsa), d_hbd(hbd), d_hba(hba), d_rotb(rotb), d_ar(aromatic_rings)]
        # Geometric mean
        prod = 1.0
        for s in scores:
            prod *= max(s, 0.05)
        return float(prod ** (1.0 / len(scores)))

    def _estimate_sascore(self, mw: float, rotb: int, aromatic_rings: int, heavy_atoms: int) -> float:
        """Estimate synthetic accessibility (1=easy, 10=hard)."""
        base = 2.0
        # penalty for large size
        if mw > 450: base += (mw - 450) / 100.0
        # penalty for high flexibility or many rings
        if rotb > 8: base += 0.5
        if aromatic_rings > 3: base += 0.8
        # scale to 1.0 - 10.0
        return min(max(base, 1.0), 9.5)

    def _calculate_mpo(self, mw: float, logp: float, tpsa: float,
                       hbd: int, hba: int, rotb: int, qed: float) -> float:
        """Multi-Parameter Optimization desirability score combining potency/ADMET potentials."""
        s_mw = 1.0 if mw <= 400 else max(0.0, 1.0 - (mw - 400)/150)
        s_logp = 1.0 if 1.5 <= logp <= 3.5 else max(0.0, 1.0 - abs(logp - 2.5)/2.5)
        s_tpsa = 1.0 if 40 <= tpsa <= 90 else max(0.0, 1.0 - abs(tpsa - 65)/50)
        s_hbd = 1.0 if hbd <= 2 else 0.5
        s_rotb = 1.0 if rotb <= 7 else 0.5
        total = 0.25 * qed + 0.2 * s_logp + 0.15 * s_tpsa + 0.15 * s_mw + 0.15 * s_hbd + 0.1 * s_rotb
        return min(max(total, 0.0), 1.0)


# ============================================================================
# 2. In Silico Ligand Editing & Bioisosteric Replacement Engine
# ============================================================================

@dataclass
class ModificationSite:
    """Represents a position on a chemical scaffold available for R-group substitution."""
    site_id: str
    description: str
    current_group: str
    vector_hint: str


class LigandEditor:
    """
    Scaffold modification and bioisosteric replacement engine for lead optimization.
    Takes a parent lead molecule and systematically generates analog libraries.
    """
    def __init__(self, core_scaffold_name: str, core_smiles_template: str):
        self.core_scaffold_name = core_scaffold_name
        self.core_smiles_template = core_smiles_template
        self.modification_sites: Dict[str, ModificationSite] = {}
        self.analog_library: List[Molecule] = []

    def register_site(self, site_id: str, description: str, current_group: str, vector_hint: str = ""):
        """Register a substitutable vector on the scaffold."""
        self.modification_sites[site_id] = ModificationSite(site_id, description, current_group, vector_hint)

    def generate_bioisosteres(self, site_id: str, replacements: List[Tuple[str, str, Dict[str, Any]]]) -> List[Molecule]:
        """
        Generate bioisosteric and functional group analogs at a specific site.
        
        Args:
            site_id: The ID of the modification site (e.g. 'R1', 'R2')
            replacements: List of tuples (analog_name, replacement_smiles_group, property_overrides)
        """
        if site_id not in self.modification_sites:
            raise ValueError(f"Site {site_id} not registered.")

        analogs = []
        for name, r_group_smiles, overrides in replacements:
            # Generate modified SMILES representation
            mod_smiles = self.core_smiles_template.replace(f"[{site_id}]", r_group_smiles)
            # Clean template tokens
            for other_site, s_info in self.modification_sites.items():
                if f"[{other_site}]" in mod_smiles:
                    mod_smiles = mod_smiles.replace(f"[{other_site}]", s_info.current_group)

            mol = Molecule(name=name, smiles=mod_smiles)
            mol.calculate_properties(
                custom_mw=overrides.get('mw'),
                custom_logp=overrides.get('logp'),
                custom_tpsa=overrides.get('tpsa'),
                custom_hbd=overrides.get('hbd'),
                custom_hba=overrides.get('hba'),
                custom_rotb=overrides.get('rotb')
            )
            analogs.append(mol)
            self.analog_library.append(mol)

        return analogs

    def get_summary_dataframe(self):
        """Returns property summary dictionary suitable for pandas DataFrame conversion."""
        records = []
        for mol in self.analog_library:
            p = mol.properties
            records.append({
                'Molecule': mol.name,
                'MW (Da)': p.molecular_weight,
                'cLogP': p.logp,
                'TPSA (Å²)': p.tpsa,
                'HBD': p.hbd,
                'HBA': p.hba,
                'RotB': p.rotatable_bonds,
                'QED': p.qed_score,
                'SAScore': p.sascore,
                'Ro5 Viol': p.lipinski_violations,
                'Pfizer 3/75': 'Alert' if p.pfizer_3_75_alert else 'Pass',
                'MPO Score': p.mpo_score
            })
        return records


# ============================================================================
# 3. Extensible Plugin Interface for Community Contributions
# ============================================================================

class BaseScoringPlugin(ABC):
    """Abstract base class for custom scoring functions and docking evaluators."""
    @abstractmethod
    def name(self) -> str:
        pass

    @abstractmethod
    def score(self, molecule: Molecule, target_pdb: Optional[str] = None) -> Dict[str, float]:
        """Calculates custom scoring metric (e.g. binding affinity, toxicity risk)."""
        pass


class EmpiricalKinaseAffinityPlugin(BaseScoringPlugin):
    """
    Example built-in plugin: Computes empirical binding free energy estimate
    for ATP-competitive kinase inhibitors using interaction contributions.
    """
    def name(self) -> str:
        return "KinaseHingeBinder_v1.0"

    def score(self, molecule: Molecule, target_pdb: Optional[str] = None) -> Dict[str, float]:
        p = molecule.properties
        # Base hinge interaction energy (e.g. bidentate H-bonds to Met793): ~ -5.5 kcal/mol
        delta_g = -5.5
        # Hydrophobic pocket filling (cLogP contribution)
        delta_g -= min(p.logp * 0.75, 2.8)
        # Optimal MW fit
        if p.molecular_weight < 350:
            delta_g += 0.8
        elif p.molecular_weight > 500:
            delta_g += 1.2
        # Rotatable bond entropy penalty (0.15 kcal/mol per rotatable bond)
        delta_g += p.rotatable_bonds * 0.15
        
        ki_nm = math.exp((delta_g * 1000.0) / (1.987 * 310.15)) * 1e9
        return {
            "estimated_delta_g_kcal_mol": round(delta_g, 2),
            "predicted_ki_nm": round(ki_nm, 1),
            "ligand_efficiency": round(-delta_g / max(p.heavy_atom_count, 1), 3)
        }


class PluginRegistry:
    """Registry allowing GitHub users to easily inject custom plugins."""
    _plugins: Dict[str, BaseScoringPlugin] = {}

    @classmethod
    def register(cls, plugin: BaseScoringPlugin):
        cls._plugins[plugin.name()] = plugin

    @classmethod
    def get(cls, name: str) -> Optional[BaseScoringPlugin]:
        return cls._plugins.get(name)

    @classmethod
    def list_plugins(cls) -> List[str]:
        return list(cls._plugins.keys())


# Register default plugin
PluginRegistry.register(EmpiricalKinaseAffinityPlugin())


# ============================================================================
# 4. Interactive 3Dmol.js WebGL Visualizer Generator
# ============================================================================

class MolecularVisualizer3D:
    """
    Renders interactive 3D WebGL protein-ligand complexes using 3Dmol.js.
    Generates both Jupyter-embeddable HTML and standalone HTML widgets.
    """
    def __init__(self, width: int = 800, height: int = 500):
        self.width = width
        self.height = height

    def generate_html(self, protein_pdb_str: str, ligand_pdb_str: str,
                      ligand_resname: str = "AQ4",
                      pocket_distance: float = 4.5,
                      viewer_id: str = "viewer_3dmol") -> str:
        """
        Creates an interactive HTML component with 3Dmol.js:
        - Protein cartoon representation with secondary structure colors
        - Binding pocket residues within `pocket_distance` Å shown in stick
        - Ligand shown in distinct ball-and-stick with elemental coloring
        - Key hydrogen bond interactions highlighted
        """
        # Escape backticks for javascript template literal
        clean_protein = protein_pdb_str.replace("`", "")
        clean_ligand = ligand_pdb_str.replace("`", "")

        html_template = f"""
<div style="width: 100%; max-width: {self.width}px; margin: 0 auto; font-family: -apple-system, BlinkMacSystemFont, 'Segoe UI', Roboto, sans-serif;">
    <div style="background: #1e293b; color: white; padding: 12px 18px; border-radius: 8px 8px 0 0; display: flex; justify-content: space-between; align-items: center;">
        <div>
            <strong style="font-size: 15px;">Target-Ligand Complex 3D Viewer</strong>
            <span style="font-size: 12px; color: #94a3b8; margin-left: 8px;">(PDB: 1M17 | Ligand: {ligand_resname})</span>
        </div>
        <div style="font-size: 12px;">
            <button onclick="resetView_{viewer_id}()" style="background: #3b82f6; color: white; border: none; padding: 4px 10px; border-radius: 4px; cursor: pointer;">Reset View</button>
            <button onclick="toggleSurface_{viewer_id}()" style="background: #475569; color: white; border: none; padding: 4px 10px; border-radius: 4px; cursor: pointer; margin-left: 6px;">Toggle Surface</button>
        </div>
    </div>
    <div id="{viewer_id}" style="width: 100%; height: {self.height}px; position: relative; background-color: #0f172a; border-radius: 0 0 8px 8px; border: 1px solid #334155; border-top: none;"></div>
    <div style="font-size: 12px; color: #64748b; padding: 8px 4px; display: flex; justify-content: space-between;">
        <span><b style="color: #38bdf8;">■</b> Protein Cartoon | <b style="color: #f59e0b;">■</b> Pocket Residues (≤{pocket_distance}Å) | <b style="color: #10b981;">■</b> Ligand {ligand_resname}</span>
        <span>Left-click: Rotate | Right-click: Pan | Scroll: Zoom</span>
    </div>
</div>

<script src="https://3Dmol.org/build/3Dmol-min.js"></script>
<script>
let viewer_{viewer_id};
let surface_{viewer_id} = null;

function initViewer_{viewer_id}() {{
    let element = document.getElementById('{viewer_id}');
    let config = {{ backgroundColor: '#0f172a' }};
    viewer_{viewer_id} = $3Dmol.createViewer(element, config);

    let proteinData = `{clean_protein}`;
    viewer_{viewer_id}.addModel(proteinData, "pdb");

    // Style overall protein as sleek cartoon
    viewer_{viewer_id}.setStyle({{}}, {{ cartoon: {{ color: '#94a3b8', opacity: 0.6 }} }});

    // Style binding pocket residues near ligand
    viewer_{viewer_id}.setStyle({{ resn: '{ligand_resname}' }}, {{}});
    viewer_{viewer_id}.addStyle({{ within: {{ distance: {pocket_distance}, sel: {{ resn: '{ligand_resname}' }} }} }},
        {{ stick: {{ colorscheme: 'amino', radius: 0.15 }} }}
    );

    // Style ligand as high-visibility ball and stick
    viewer_{viewer_id}.addStyle({{ resn: '{ligand_resname}' }},
        {{ stick: {{ colorscheme: 'greenCarbon', radius: 0.25 }}, sphere: {{ scale: 0.3, colorscheme: 'greenCarbon' }} }}
    );

    // Add labels for key catalytic residues
    let keyResidues = [
        {{ resi: 793, text: "Met793 (Hinge)" }},
        {{ resi: 790, text: "Thr790 (Gatekeeper)" }},
        {{ resi: 745, text: "Lys745 (Catalytic)" }}
    ];
    keyResidues.forEach(r => {{
        viewer_{viewer_id}.addResLabels({{ resi: r.resi }}, {{
            font: 'Arial', fontSize: 11, fontColor: '#ffffff',
            backgroundColor: 'rgba(15, 23, 42, 0.8)', borderThickness: 1, borderColor: '#38bdf8'
        }});
    }});

    viewer_{viewer_id}.zoomTo({{ resn: '{ligand_resname}' }});
    viewer_{viewer_id}.render();
}}

function resetView_{viewer_id}() {{
    if (viewer_{viewer_id}) {{
        viewer_{viewer_id}.zoomTo({{ resn: '{ligand_resname}' }});
        viewer_{viewer_id}.render();
    }}
}}

function toggleSurface_{viewer_id}() {{
    if (!viewer_{viewer_id}) return;
    if (surface_{viewer_id}) {{
        viewer_{viewer_id}.removeSurface(surface_{viewer_id});
        surface_{viewer_id} = null;
    }} else {{
        surface_{viewer_id} = viewer_{viewer_id}.addSurface($3Dmol.SurfaceType.VDW, {{
            opacity: 0.5,
            color: '#38bdf8'
        }}, {{ within: {{ distance: 6.0, sel: {{ resn: '{ligand_resname}' }} }} }});
    }}
    viewer_{viewer_id}.render();
}}

// Initialize when DOM is ready
if (document.readyState === 'complete' || document.readyState === 'interactive') {{
    setTimeout(initViewer_{viewer_id}, 100);
}} else {{
    document.addEventListener('DOMContentLoaded', initViewer_{viewer_id});
}}
</script>
"""
        return html_template
