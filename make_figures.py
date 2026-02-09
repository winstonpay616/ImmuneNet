#!/usr/bin/env python3
"""
make_figures.py

Generate:
 - degree distribution (log-log)
 - top-10 hubs barplot
 - community x node-type heatmap

for a given subnetwork (cytokine-centered or receptor-centered),
starting from nodes_metrics.csv produced by analyze_network.py.
"""

import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

# -------------------------------------------------------------------
# CONFIG
# -------------------------------------------------------------------
# Path can be changed depending on the desired subgraph (cytokine/receptor) to plot.

BASE_DIR = "output/cytokine"   
#BASE_DIR = "output/receptor"   

NODES_METRICS_FILE = os.path.join(BASE_DIR, "nodes_metrics.csv")
FIG_DIR = os.path.join(BASE_DIR, "figures")

os.makedirs(FIG_DIR, exist_ok=True)


# -------------------------------------------------------------------
# 1) Degree distribution (log–log)
# -------------------------------------------------------------------
def plot_degree_distribution(nodes_df, out_path, title="Degree distribution"):
    """
    Plot the degree distribution on log-log scale.

    nodes_df: DataFrame with a 'degree' column
    out_path: path to save the PNG
    """
    deg = nodes_df["degree"].astype(int)

    # Compute P(k)
    unique_deg, counts = np.unique(deg, return_counts=True)
    prob = counts / counts.sum()

    # Remove zeros if present (non ha senso log(0))
    mask = unique_deg > 0
    unique_deg = unique_deg[mask]
    prob = prob[mask]

    plt.figure()
    plt.scatter(unique_deg, prob)
    plt.xscale("log")
    plt.yscale("log")
    plt.xlabel("Degree k")
    plt.ylabel("P(k)")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"[OK] Saved degree distribution plot -> {out_path}")


# -------------------------------------------------------------------
# 2) Top-10 hubs barplot (by degree or betweenness)
# -------------------------------------------------------------------
import matplotlib.pyplot as plt
import pandas as pd
import os

# Same color map as in Gephi
COLOR_BY_TYPE = {
    "cell": "#4CAF50",
    "cytokine": "#FF9800",
    "protein": "#9C27B0",
    "receptor": "#F44336",
    "other": "#9E9E9E"
}

def get_node_color(node_type: str) -> str:
    """Map a node type string to a display color."""
    if pd.isna(node_type):
        return COLOR_BY_TYPE["other"]
    t = str(node_type).strip().lower()
    return COLOR_BY_TYPE.get(t, COLOR_BY_TYPE["other"])


def plot_top_hubs(nodes_df,
                  out_path,
                  metric="degree",
                  n=10,
                  title="Top hubs by degree"):
    """
    Plot a bar chart of top-n nodes according to a centrality metric,
    coloring bars and labels according to node type (cytokine, receptor, protein, ...).

    nodes_df must contain:
      - 'node'  (node id)
      - 'type'  (biological type)
      - metric  (e.g. 'degree', 'betweenness', 'pagerank')
    """
    if metric not in nodes_df.columns:
        raise ValueError(f"Metric '{metric}' not found in nodes_df columns")

    # Sort and take top-n
    df_sorted = nodes_df.sort_values(by=metric, ascending=False).head(n)

    # Reverse so the highest is at the top of the horizontal bar chart
    df_sorted = df_sorted.iloc[::-1]

    # Build color list based on 'type'
    types = df_sorted["type"].fillna("other")
    colors = [get_node_color(t) for t in types]

    plt.figure(figsize=(6, 4))
    bars = plt.barh(df_sorted["node"], df_sorted[metric])

    # Color bars
    for bar, c in zip(bars, colors):
        bar.set_color(c)

    # Axis labels / title
    plt.xlabel(metric.capitalize())
    plt.ylabel("Node")
    plt.title(title)

    # Color y-tick labels according to node type
    ax = plt.gca()
    ax.set_yticks(range(len(df_sorted)))
    ax.set_yticklabels(df_sorted["node"])

    for label, c in zip(ax.get_yticklabels(), colors):
        label.set_color(c)

    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"[OK] Saved top hubs plot -> {out_path}")


# -------------------------------------------------------------------
# 3) Heatmap community × node type
# -------------------------------------------------------------------
def plot_community_type_heatmap(nodes_df, out_path,
                                max_communities=10,
                                title="Community vs node type"):
    """
    Heatmap delle counts (node type x community).
    Filtra alle comunità più grandi (max_communities) per leggibilità.
    Assumes:
      - 'community' column with integer ids
      - 'type' column with biological type (cytokine, receptor, protein, ...)
    """

    if "community" not in nodes_df.columns:
        raise ValueError("nodes_df must contain a 'community' column")

    if "type" not in nodes_df.columns:
        raise ValueError("nodes_df must contain a 'type' column for the heatmap")

    # Prendi comunità più grandi
    comm_counts = nodes_df["community"].value_counts().sort_values(ascending=False)
    top_communities = comm_counts.head(max_communities).index.tolist()

    df_sub = nodes_df[nodes_df["community"].isin(top_communities)].copy()

    # Sostituisci NaN o vuoti in 'type' con "other"
    df_sub["type"] = df_sub["type"].fillna("other")

    # Crea tabella type x community
    table = pd.crosstab(df_sub["type"], df_sub["community"])

    plt.figure(figsize=(8, 5))
    im = plt.imshow(table, aspect="auto")
    plt.colorbar(im, label="Node count")

    plt.yticks(range(len(table.index)), table.index)
    plt.xticks(range(len(table.columns)), table.columns)
    plt.xlabel("Community ID")
    plt.ylabel("Node type")
    plt.title(title)
    plt.tight_layout()
    plt.savefig(out_path, dpi=300)
    plt.close()
    print(f"[OK] Saved community-type heatmap -> {out_path}")


# -------------------------------------------------------------------
# MAIN
# -------------------------------------------------------------------
def main():
    """Load node metrics and render all summary figures."""
    if not os.path.exists(NODES_METRICS_FILE):
        raise FileNotFoundError(f"nodes_metrics.csv not found at {NODES_METRICS_FILE}")

    nodes_df = pd.read_csv(NODES_METRICS_FILE)

    # --- Decide i nomi in base al tipo di rete (cytokine / receptor) ---
    if "cytokine" in BASE_DIR.lower():
        prefix = "cytokine"
        degree_title = "Degree distribution of the cytokine-centered subnetwork"
        hubs_title = "Top-10 hubs in the cytokine-centered subnetwork"
        heatmap_title = "Community vs node type (cytokine-centered)"
    elif "receptor" in BASE_DIR.lower():
        prefix = "receptor"
        degree_title = "Degree distribution of the receptor-centered subnetwork"
        hubs_title = "Top-10 hubs in the receptor-centered subnetwork"
        heatmap_title = "Community vs node type (receptor-centered)"
    else:
        prefix = "network"
        degree_title = "Degree distribution"
        hubs_title = "Top-10 hubs"
        heatmap_title = "Community vs node type"

    # 1) Degree distribution
    out_deg = os.path.join(FIG_DIR, f"{prefix}_degree_distribution.png")
    plot_degree_distribution(nodes_df, out_deg, title=degree_title)

    # 2) Top-10 hubs by degree
    out_hubs = os.path.join(FIG_DIR, f"{prefix}_top_hubs.png")
    plot_top_hubs(nodes_df, out_hubs, metric="degree", n=10, title=hubs_title)

    # 3) Heatmap community x type
    out_heatmap = os.path.join(FIG_DIR, f"{prefix}_community_type_heatmap.png")
    plot_community_type_heatmap(nodes_df, out_heatmap, title=heatmap_title)


if __name__ == "__main__":
    main()
