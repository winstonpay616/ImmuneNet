#!/usr/bin/env python3
import os
import pandas as pd

DATA_DIR = "data"

FULL_EDGES = os.path.join(DATA_DIR, "edges.csv")
FULL_NODES = os.path.join(DATA_DIR, "nodes.csv")


def build_subgraph(label, seed_types):
    """
    Build and write a subgraph centered on seed node types.

    label: 'cytokine' or 'receptor'
    seed_types: list of node types to use as seeds, e.g. ['cytokine'] or ['receptor']
    """
    # Load full graph
    if not os.path.exists(FULL_EDGES) or not os.path.exists(FULL_NODES):
        raise FileNotFoundError("Run convert_to_edgelist.py first to create data/edges.csv and data/nodes.csv")

    edges = pd.read_csv(FULL_EDGES)
    nodes = pd.read_csv(FULL_NODES)

    if "id" not in nodes.columns or "type" not in nodes.columns:
        raise ValueError("nodes.csv must contain columns 'id' and 'type'")

    # Map id -> type
    nodes_idx = nodes.set_index("id")
    types = nodes_idx["type"].to_dict()

    # Seed nodes: nodes whose type is in seed_types
    seed_nodes = {nid for nid, t in types.items() if t in seed_types}
    print(f"[{label}] seed nodes: {len(seed_nodes)}")

    if not seed_nodes:
        print(f"[{label}] WARNING: no seed nodes found with types: {seed_types}")
    
    # Keep edges that touch at least one seed node
    mask = edges["source"].isin(seed_nodes) | edges["target"].isin(seed_nodes)
    edges_sub = edges[mask].copy()
    print(f"[{label}] edges after seed filter: {len(edges_sub)}")

    if edges_sub.empty:
        print(f"[{label}] WARNING: no edges in subgraph, nothing to write.")
        return

    # Nodes in this subgraph = all endpoints of remaining edges
    node_ids = sorted(set(edges_sub["source"]) | set(edges_sub["target"]))
    nodes_sub = nodes[nodes["id"].isin(node_ids)].copy()

    # Keep only the columns we need in correct order
    nodes_sub = nodes_sub[["id", "type"]]

    # Write out
    edges_out = os.path.join(DATA_DIR, f"edges_{label}.csv")
    nodes_out = os.path.join(DATA_DIR, f"nodes_{label}.csv")
    edges_sub.to_csv(edges_out, index=False)
    nodes_sub.to_csv(nodes_out, index=False)

    print(f"[{label}] wrote {edges_out} ({len(edges_sub)} edges) and {nodes_out} ({len(nodes_sub)} nodes)\n")


def main():
    """Generate cytokine- and receptor-centered subnetworks."""
    # Cytokine-centered subgraph
    build_subgraph(label="cytokine", seed_types=["cytokine"])

    # Receptor-centered subgraph
    build_subgraph(label="receptor", seed_types=["receptor"])


if __name__ == "__main__":
    main()
