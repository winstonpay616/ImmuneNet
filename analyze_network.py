#!/usr/bin/env python3
# analyze_network.py

import os
import math
import pandas as pd
import networkx as nx

# ---------------------
# Configuration
# ---------------------
# Path can be changed depending on the desired subgraph (cytokine/receptor) to analyse.

DATA_DIR = "data"
OUTPUT_DIR = "output/cytokine"
#OUTPUT_DIR = "output/receptor"

GLOBAL_EDGES_FILE = os.path.join(DATA_DIR, "edges.csv")
GLOBAL_NODES_FILE = os.path.join(DATA_DIR, "nodes.csv")

CYTOKINE_EDGES_FILE = os.path.join(DATA_DIR, "edges_cytokine.csv")
CYTOKINE_NODES_FILE = os.path.join(DATA_DIR, "nodes_cytokine.csv")

RECEPTOR_EDGES_FILE = os.path.join(DATA_DIR, "edges_receptor.csv")
RECEPTOR_NODES_FILE = os.path.join(DATA_DIR, "nodes_receptor.csv")

if OUTPUT_DIR.endswith("cytokine"):
   EDGES_FILE = CYTOKINE_EDGES_FILE
   NODES_FILE = CYTOKINE_NODES_FILE
if OUTPUT_DIR.endswith("receptor"):
   EDGES_FILE = RECEPTOR_EDGES_FILE
   NODES_FILE = RECEPTOR_NODES_FILE

print(f"\n[INFO] Using edges: {EDGES_FILE}\n")

# Se la rete è PPI/metabolica ecc. lasciala come non diretta
IS_DIRECTED = False   

# Column names expected in edges.csv
EDGE_COL_SOURCE = "source"
EDGE_COL_TARGET = "target"
EDGE_COL_WEIGHT = "weight"        # optional
EDGE_COL_INTERACTION = "interaction"  # optional
EDGE_COL_EVIDENCE = "evidence"        # optional

# Column names expected in nodes.csv
NODE_COL_ID = "id"
NODE_COL_TYPE = "type"  # optional (e.g., "cell", "cytokine", "protein", "gene")

# Visual mappings
COLOR_BY_TYPE = {
    "cell": "#4CAF50",
    "cytokine": "#2196F3",
    "protein": "#9C27B0",
    "gene": "#9C27B0",
    "receptor": "#FF9800",
    "pathogen": "#F44336",
    "other": "#9E9E9E"
}
DEFAULT_NODE_COLOR = "#9E9E9E"

NODE_SIZE_MIN = 10
NODE_SIZE_MAX = 40

# ---------------------
# Helpers
# ---------------------
def ensure_outdir():
    """Create the output directory if it does not exist."""
    os.makedirs(OUTPUT_DIR, exist_ok=True)

def load_edges():
    """Load the edges CSV and validate required columns."""
    if not os.path.exists(EDGES_FILE):
        raise FileNotFoundError(f"Missing {EDGES_FILE}")
    df = pd.read_csv(EDGES_FILE)
    for col in [EDGE_COL_SOURCE, EDGE_COL_TARGET]:
        if col not in df.columns:
            raise ValueError(f"edges.csv must contain column '{col}'")
    return df

def load_nodes():
    """Load the nodes CSV if present and validate required columns."""
    if os.path.exists(NODES_FILE):
        df = pd.read_csv(NODES_FILE)
        if NODE_COL_ID not in df.columns:
            raise ValueError(f"nodes.csv must contain column '{NODE_COL_ID}'")
        return df
    return None

def build_graph(edges_df, nodes_df=None, is_directed=True):
    """Build a NetworkX graph from edge and optional node tables."""
    G = nx.DiGraph() if is_directed else nx.Graph()

    # Add nodes (from nodes_df, if provided)
    if nodes_df is not None:
        attrs = nodes_df.set_index(NODE_COL_ID).to_dict(orient="index")
        for nid, data in attrs.items():
            G.add_node(nid, **data)

    # Add edges
    has_weight = EDGE_COL_WEIGHT in edges_df.columns
    for _, row in edges_df.iterrows():
        u = str(row[EDGE_COL_SOURCE])
        v = str(row[EDGE_COL_TARGET])
        data = {}
        if has_weight and pd.notnull(row[EDGE_COL_WEIGHT]):
            data["weight"] = float(row[EDGE_COL_WEIGHT])
        if EDGE_COL_INTERACTION in edges_df.columns and pd.notnull(row.get(EDGE_COL_INTERACTION, None)):
            data["interaction"] = str(row[EDGE_COL_INTERACTION])
        if EDGE_COL_EVIDENCE in edges_df.columns and pd.notnull(row.get(EDGE_COL_EVIDENCE, None)):
            data["evidence"] = str(row[EDGE_COL_EVIDENCE])

        G.add_edge(u, v, **data)

        # If nodes weren’t predeclared, ensure existence
        if u not in G.nodes:
            G.add_node(u)
        if v not in G.nodes:
            G.add_node(v)

    return G

def safe_normalize(values, new_min=NODE_SIZE_MIN, new_max=NODE_SIZE_MAX):
    """Normalize numeric values to a target range, guarding against NaN/inf."""
    vals = list(values)
    finite_vals = [v for v in vals if pd.notnull(v) and math.isfinite(v)]
    if not finite_vals:
        return [new_min for _ in vals]
    vmin, vmax = min(finite_vals), max(finite_vals)
    if vmax == vmin:
        return [int((new_min + new_max) / 2)] * len(vals)
    out = []
    for v in vals:
        if pd.isnull(v) or not math.isfinite(v):
            out.append(new_min)
        else:
            x = (v - vmin) / (vmax - vmin)
            out.append(int(new_min + x * (new_max - new_min)))
    return out

def node_color(node_attrs):
    """Map node attributes to a display color using the type field."""
    t = (node_attrs or {}).get(NODE_COL_TYPE, None)
    if pd.isnull(t) or t is None:
        return DEFAULT_NODE_COLOR
    t = str(t).lower()
    return COLOR_BY_TYPE.get(t, DEFAULT_NODE_COLOR)

def undirected_view(G):
    """Return an undirected copy for algorithms that require it (e.g., path length, clustering, communities)."""
    return G.to_undirected(as_view=False)

# ---------------------
# Main analysis steps
# ---------------------
def analyze(G):
    """Compute network metrics, per-node stats, and community assignments."""
    from networkx.algorithms.community import greedy_modularity_communities, modularity

    report = {}

    # Basic stats
    report["n_nodes"] = G.number_of_nodes()
    report["n_edges"] = G.number_of_edges()
    report["is_directed"] = G.is_directed()
    report["density"] = nx.density(G)

    # Undirected projection for many metrics
    U = undirected_view(G)

    # Connected components (sempre sulla versione non diretta)
    components = list(nx.connected_components(U))
    n_components = len(components)
    report["n_components_undirected"] = n_components

    if n_components > 0:
        giant = max(components, key=len)
        giant_size = len(giant)
        report["giant_component_size"] = giant_size
        report["giant_component_frac"] = giant_size / G.number_of_nodes()
    else:
        giant = set()
        report["giant_component_size"] = 0
        report["giant_component_frac"] = 0.0

    # Per grafi diretti, tieni anche weak/strong component count
    if G.is_directed():
        report["n_weak_components"] = nx.number_weakly_connected_components(G)
        report["n_strong_components"] = nx.number_strongly_connected_components(G)

    # Degree stats
    deg = dict(G.degree())
    indeg = dict(G.in_degree()) if G.is_directed() else None
    outdeg = dict(G.out_degree()) if G.is_directed() else None

    # Centralities (use weighted if available)
    weight_kw = "weight" if any("weight" in e for _, _, e in G.edges(data=True)) else None
    bet = nx.betweenness_centrality(G, weight=weight_kw, normalized=True)

    try:
        clo = nx.closeness_centrality(G, distance=(weight_kw if weight_kw else None))
    except Exception:
        clo = {n: float("nan") for n in G.nodes}

    try:
        if G.is_directed():
            pg = nx.pagerank(G, weight=weight_kw)
        else:
            pg = nx.pagerank(U, weight=weight_kw)
    except Exception:
        pg = {n: float("nan") for n in G.nodes}

    # k-core (su U)
    try:
        core = nx.core_number(U)
    except Exception:
        core = {n: 0 for n in G.nodes}

    # Assortativity (by degree) su U
    try:
        assort = nx.degree_assortativity_coefficient(U)
    except Exception:
        assort = float("nan")
    report["degree_assortativity"] = assort

    # Average clustering coefficient su U
    try:
        if weight_kw:
            avg_clust = nx.average_clustering(U, weight=weight_kw)
        else:
            avg_clust = nx.average_clustering(U)
    except Exception:
        avg_clust = float("nan")
    report["avg_clustering"] = avg_clust

    # Path length e diameter sulla giant component
    if len(giant) > 1:
        S = U.subgraph(giant).copy()
        try:
            avg_path = nx.average_shortest_path_length(S, weight=weight_kw)
        except Exception:
            avg_path = float("nan")
        try:
            diam = nx.diameter(S)
        except Exception:
            diam = float("nan")
    else:
        avg_path = float("nan")
        diam = float("nan")
    report["avg_path_length_giant"] = avg_path
    report["diameter_giant"] = diam

    # Community detection (greedy modularity) su U
    try:
        communities = list(greedy_modularity_communities(U, weight=weight_kw))
    except Exception:
        communities = [set(U.nodes())]

    report["n_communities_greedy"] = len(communities)
    try:
        mod_value = modularity(U, communities, weight=weight_kw)
    except Exception:
        mod_value = float("nan")
    report["modularity_greedy"] = mod_value

    # Map node -> community id
    node_to_comm = {}
    for cid, com in enumerate(communities):
        for n in com:
            node_to_comm[n] = cid

    # Build node table
    rows = []
    for n, attrs in G.nodes(data=True):
        rows.append({
            "node": n,
            "type": attrs.get(NODE_COL_TYPE, None),
            "degree": deg.get(n, 0),
            "in_degree": indeg.get(n, None) if indeg else None,
            "out_degree": outdeg.get(n, None) if outdeg else None,
            "betweenness": bet.get(n, float("nan")),
            "closeness": clo.get(n, float("nan")),
            "pagerank": pg.get(n, float("nan")),
            "kcore": core.get(n, 0),
            "community": node_to_comm.get(n, -1),
        })
    nodes_df = pd.DataFrame(rows)

    # Edge table (preserve attributes)
    e_rows = []
    for u, v, d in G.edges(data=True):
        e_rows.append({
            "source": u, "target": v,
            **({"weight": d.get("weight")} if "weight" in d else {}),
            **({"interaction": d.get("interaction")} if "interaction" in d else {}),
            **({"evidence": d.get("evidence")} if "evidence" in d else {})
        })
    edges_df = pd.DataFrame(e_rows)

    # Attach report
    report_df = pd.DataFrame([report])

    return nodes_df, edges_df, report_df, node_to_comm

def export_results(nodes_df, edges_df, report_df):
    """Write analysis tables to CSV files under the output directory."""
    ensure_outdir()
    nodes_path = os.path.join(OUTPUT_DIR, "nodes_metrics.csv")
    edges_path = os.path.join(OUTPUT_DIR, "edges_clean.csv")
    report_path = os.path.join(OUTPUT_DIR, "report_summary.csv")
    nodes_df.to_csv(nodes_path, index=False)
    edges_df.to_csv(edges_path, index=False)
    report_df.to_csv(report_path, index=False)
    print(f"[OK] Saved: {nodes_path}\n[OK] Saved: {edges_path}\n[OK] Saved: {report_path}")

def export_gephi(G, nodes_df, out_path):
    """Export the graph to GEXF with node metrics attached for Gephi."""
    # Attach node attributes from nodes_df to the graph
    metrics = nodes_df.set_index("node").to_dict(orient="index")

    # Make a copy to avoid mutating original G
    H = G.copy()

    for n in H.nodes():
        attrs = metrics.get(n, {})
        clean_attrs = {}
        for k, v in attrs.items():
            # Cast numpy types to native Python or strings (GEXF is picky)
            if pd.isna(v):
                continue
            if isinstance(v, (int, float, str, bool)):
                clean_attrs[k] = v
            else:
                clean_attrs[k] = str(v)
        nx.set_node_attributes(H, {n: clean_attrs})

    nx.write_gexf(H, out_path)
    print(f"[OK] Saved Gephi GEXF file: {out_path}")

def main():
    """Load inputs, run analysis, and export CSV/GEXF outputs."""
    edges_df = load_edges()
    nodes_df = load_nodes()

    # If all edges are undirected (PPI), set IS_DIRECTED = False at the top.
    G = build_graph(edges_df, nodes_df, is_directed=IS_DIRECTED)

    nodes_metrics, edges_clean, report_df, _ = analyze(G)
    export_results(nodes_metrics, edges_clean, report_df)
    
    gephi_path = os.path.join(OUTPUT_DIR, "network.gexf")
    export_gephi(G, nodes_metrics, gephi_path)

if __name__ == "__main__":
    main()
