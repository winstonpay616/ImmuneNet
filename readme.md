# ImmuneNet

ImmuneNet builds and analyzes immune interaction networks from an InnateDB MITAB dataset. The pipeline converts raw MITAB interactions into a cleaned edge list, derives cytokine- and receptor-centered subnetworks, computes network metrics/communities, and generates summary figures and Gephi exports files.

The core behavior is the following: start from curated human–human interactions, infer coarse node types (cytokine, receptor, protein), and analyze immune-specific subnetworks.

## Pipeline

1. `convert_to_edgelist.py`
   - Reads `data/innatedb_all.mitab.txt`
   - Filters to human–human interactions
   - Resolves gene/protein symbols
   - Deduplicates undirected pairs
   - Parses confidence into edge weights
   - Writes `data/edges.csv` and `data/nodes.csv`

2. `build_subnetworks.py`
   - Builds cytokine- and receptor-centered subnetworks by keeping edges that touch seed node types
   - Writes `data/edges_cytokine.csv`, `data/nodes_cytokine.csv`, `data/edges_receptor.csv`, `data/nodes_receptor.csv`

3. `analyze_network.py`
   - Loads a chosen subnetwork (cytokine or receptor)
   - Computes metrics (degree, betweenness, closeness, PageRank, k-core)
   - Runs community detection (greedy modularity)
   - Exports CSVs and a Gephi `.gexf`

4. `make_figures.py`
   - Generates:
     - Degree distribution (log–log)
     - Top-10 hubs barplot
     - Community × node-type heatmap

## Quickstart

Install Python dependencies (if not already available):

```bash
pip install pandas networkx numpy matplotlib
```

Run the full pipeline:

```bash
python3 convert_to_edgelist.py
python3 build_subnetworks.py

# Edit analyze_network.py to pick cytokine vs receptor (see below)
python3 analyze_network.py

# Edit make_figures.py to match the analyzed output dir
python3 make_figures.py
```

## Configuration Tips

- **Input MITAB file**: update `MITAB_PATH` in `convert_to_edgelist.py` if you replace the dataset.
- **Choose subnetwork to analyze** (in `analyze_network.py`):
  - Set `OUTPUT_DIR` to `output/cytokine` or `output/receptor`
  - Set `EDGES_FILE` / `NODES_FILE` to the matching `data/edges_*.csv` and `data/nodes_*.csv`
- **Match figure generation** (in `make_figures.py`):
  - Set `BASE_DIR` to the same output dir used in `analyze_network.py`.

## Outputs

- **Global network**
  - `data/edges.csv`
  - `data/nodes.csv`
- **Subnetworks**
  - `data/edges_cytokine.csv`, `data/nodes_cytokine.csv`
  - `data/edges_receptor.csv`, `data/nodes_receptor.csv`
- **Analysis (per subnetwork)**
  - `output/<subnet>/nodes_metrics.csv`
  - `output/<subnet>/edges_clean.csv`
  - `output/<subnet>/report_summary.csv`
  - `output/<subnet>/network.gexf`
- **Figures (per subnetwork)**
  - `output/<subnet>/figures/*_degree_distribution.png`
  - `output/<subnet>/figures/*_top_hubs.png`
  - `output/<subnet>/figures/*_community_type_heatmap.png`

## Repository Layout

- `convert_to_edgelist.py` – MITAB → edges/nodes CSVs
- `build_subnetworks.py` – cytokine/receptor subnetworks
- `analyze_network.py` – network metrics + GEXF export
- `make_figures.py` – summary plots
- `data/` – input MITAB and generated CSVs
- `output/` – analysis outputs and figures

## Notes

- The graph is treated as **undirected** by default (`IS_DIRECTED = False` in `analyze_network.py`).
- Node types are inferred by regex rules (cytokine vs receptor), otherwise defaulting to `protein`.
- A few higher-level nodes (e.g., `Macrophage`, `T_cell`) are added to `nodes.csv` for convenience but are not connected by edges unless added manually.
