# PLASTRO 🧬

[![PyPI version](https://badge.fury.io/py/plastro.svg)](https://badge.fury.io/py/plastro)
[![Conda version](https://img.shields.io/conda/vn/conda-forge/plastro.svg)](https://anaconda.org/conda-forge/plastro)
[![Documentation Status](https://readthedocs.org/projects/plastro/badge/?version=latest)](https://plastro.readthedocs.io/en/latest/?badge=latest)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)

**PLASTRO** is a Python package for simulating and analyzing cellular plasticity in single-cell data. It provides comprehensive tools for studying how cells transition between different phenotypic states and how these transitions relate to lineage relationships.

## 🚀 Key Features

- **🔬 Plasticity Simulation**: Multiple methods for simulating cellular plasticity including random walk plasticity and discrete cluster switches
- **🧬 Lineage Tracing Integration**: Full integration with CRISPR-based lineage tracing data and Cassiopeia trees  
- **📊 PLASTRO Score**: Novel overlap-based metric for quantifying cellular plasticity from combined lineage and phenotypic data
- **🌳 Phylogenetic Analysis**: Neighbor-joining tree construction and phylogenetic distance calculations
- **🧪 Data Simulation**: Generate realistic synthetic single-cell datasets with branching differentiation trajectories
- **📈 Comprehensive Analysis**: Distance metrics, archetype analysis, and visualization tools for plasticity research

## 📦 Installation

### Quick Install

```bash
pip install plastro
```

### Conda Install

```bash
conda install -c conda-forge plastro
```

### Development Install

```bash
git clone https://github.com/username/plastro.git
cd plastro
pip install -e \".[dev]\"
```

### Install with Phylogenetic Support

```bash
pip install plastro[phylo]
# or for editable install
pip install -e \".[phylo]\"
```

## 🎯 Quick Start

### Generate Synthetic Data

```python
import plastro

# Generate realistic synthetic single-cell data
adata = plastro.simulate_realistic_dataset(
    n_cell_types=6,
    cells_per_type=200, 
    n_genes=25,
    seed=42
)

# Simulate random walk plasticity
walk_params = {
    100: 0.1,   # 10% of cells perform short walks (low plasticity)
    500: 0.05,  # 5% of cells perform medium walks (medium plasticity)
    1000: 0.02  # 2% of cells perform long walks (high plasticity)
}

plastic_adata = plastro.random_walk_plasticity(adata, walk_params)
print(f\"Mean phenotypic change: {plastic_adata.obs['change_in_phenotype'].mean():.3f}\")
```

### Cluster Switch Plasticity

```python
# Simulate discrete cluster switches
degree_props = {
    1: 0.3,  # 30% switch to closest cluster
    2: 0.2,  # 20% switch to 2nd closest cluster  
    3: 0.1   # 10% switch to 3rd closest cluster
}

plastic_adata, plasticity_info = plastro.leiden_switch_plasticity(adata, degree_props)
```

### PLASTRO Score Calculation

```python
# Compute PLASTRO scores (requires lineage tracing data)
plastro_scores = plastro.PLASTRO_score(
    character_matrix=character_matrix,
    ad=adata,
    threshold=0.95,
    latent_space_key='X_pca'
)
```

### Lineage Tree Construction

```python
# Construct phylogenetic tree from single-cell data
cass_tree = plastro.simulate_lineage_tracing(
    sim_ad=full_data,
    terminal_ad=observed_cells,
    latent_space_key='X_dc'
)

# Extract character matrix
character_matrix = cass_tree.character_matrix
```

## 📚 Documentation

Comprehensive documentation is available at [plastro.readthedocs.io](https://plastro.readthedocs.io/).

### Key Documentation Sections

- **[Installation Guide](https://plastro.readthedocs.io/en/latest/installation.html)**: Detailed installation instructions
- **[API Reference](https://plastro.readthedocs.io/en/latest/api/)**: Complete function documentation
- **[Tutorials](https://plastro.readthedocs.io/en/latest/tutorials/)**: Step-by-step guides
- **[Example Notebooks](https://plastro.readthedocs.io/en/latest/examples/)**: Jupyter notebook examples

## 🔧 Core Modules

### `plastro.plasticity`
Core plasticity simulation functions:
- `random_walk_plasticity()`: Simulate plasticity via random walks
- `leiden_switch_plasticity()`: Discrete cluster-based plasticity
- `subclone_switch_plasticity()`: Phylogenetic subclone plasticity
- `visualize_walk()`: Visualize random walk paths

### `plastro.overlap`
PLASTRO score computation:
- `PLASTRO_score()`: Main function for computing plasticity scores
- `PLASTRO_overlaps()`: Compute overlap matrices
- `overlaps_to_score()`: Convert overlaps to final scores

### `plastro.lineage_simulation`
CRISPR-based lineage tracing simulation:
- `simulate_lineage_tracing()`: End-to-end lineage simulation
- `construct_tree()`: Build trees from phenotypic data
- `introduce_crispr_mutations()`: Add CRISPR mutations to trees

### `plastro.distances`
Phenotypic distance calculations:
- `archetype_distance()`: Archetype-based distances
- `euclidean_distance()`: Standard Euclidean distances
- `correlation_distance()`: Gene expression correlations

### `plastro.phylo`
Phylogenetic analysis tools:
- `neighbor_joining()`: Construct phylogenetic trees
- `calculate_tree_distances()`: Compute tree-based distances
- `bootstrap_tree()`: Bootstrap consensus trees

### `plastro.phenotype_simulation`
Synthetic phenotypic data generation:
- `simulate_realistic_dataset()`: Generate complete synthetic datasets
- `generate_ad()`: Create data from tree structures
- `create_random_binary_tree()`: Build differentiation hierarchies

## 🧪 Example Workflow

Here's a complete example workflow:

```python
import plastro
import scanpy as sc
import pandas as pd

# 1. Load and preprocess data
adata = sc.read_h5ad('single_cell_data.h5ad')
sc.pp.neighbors(adata)
sc.tl.umap(adata)
sc.tl.leiden(adata)

# 2. Simulate plasticity
walk_params = {200: 0.15, 800: 0.08}
plastic_adata = plastro.random_walk_plasticity(
    adata, 
    walk_params, 
    save_to='results/'
)

# 3. Construct lineage tree
tree = plastro.construct_tree(
    adata, 
    plastic_adata, 
    latent_space_key='X_pca'
)

# 4. Add CRISPR mutations
cass_tree = plastro.introduce_crispr_mutations(tree)

# 5. Compute PLASTRO scores
scores = plastro.PLASTRO_score(
    cass_tree.character_matrix,
    plastic_adata,
    threshold=0.95
)

# 6. Analyze results
print(f\"Mean PLASTRO score: {scores.mean().values[0]:.3f}\")
print(f\"Plasticity detected in {(scores < 0.5).sum().values[0]} cells\")
```

## 🛠️ Dependencies

**Core requirements:**
- Python ≥ 3.8
- NumPy ≥ 1.20.0
- Pandas ≥ 1.3.0
- SciPy ≥ 1.7.0
- Scanpy ≥ 1.8.0
- NetworkX ≥ 2.6.0

**Optional dependencies:**
- `cassiopeia-lineage` (for lineage tracing simulation)
- `py-pcha` (for archetype analysis)
- `ete3` (for phylogenetic trees)
- `scikit-bio` (for robust neighbor-joining algorithm - install with `pip install plastro[phylo]`)

## 📊 Data Requirements

PLASTRO works with:
- **Single-cell RNA-seq data** in AnnData format
- **Character matrices** from CRISPR lineage tracing
- **Phylogenetic trees** in Newick format or ETE3 objects
- **Dimensionality reduction** coordinates (PCA, UMAP, diffusion maps)

## 🤝 Contributing

We welcome contributions! Please see our [Contributing Guide](CONTRIBUTING.md) for details.

### Development Setup

```bash
git clone https://github.com/username/plastro.git
cd plastro
pip install -e \".[dev]\"
pre-commit install
```

### Running Tests

```bash
pytest tests/
```

## 📄 Citation

If you use PLASTRO in your research, please cite:

```bibtex
@software{plastro2024,
  author = {Persad, Sitara},
  title = {PLASTRO: A Python package for simulating cellular plasticity},
  url = {https://github.com/username/plastro},
  version = {0.1.0},
  year = {2024}
}
```

## 📝 License

This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.

## 🆘 Support

- **Documentation**: [plastro.readthedocs.io](https://plastro.readthedocs.io/)
- **Issues**: [GitHub Issues](https://github.com/username/plastro/issues)
- **Discussions**: [GitHub Discussions](https://github.com/username/plastro/discussions)

## 🔗 Related Projects

- **[Scanpy](https://scanpy.readthedocs.io/)**: Single-cell analysis in Python
- **[Cassiopeia](https://cassiopeia-lineage.readthedocs.io/)**: Lineage tracing analysis
- **[AnnData](https://anndata.readthedocs.io/)**: Annotated data structures
- **[NetworkX](https://networkx.org/)**: Network analysis in Python

---

**PLASTRO** - Enabling comprehensive analysis of cellular plasticity in single-cell data 🧬✨