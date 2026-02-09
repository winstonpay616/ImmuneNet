
#!/usr/bin/env python3
import os
import re
import pandas as pd

MITAB_PATH = "data/innatedb_all.mitab.txt"   # update if needed
OUT_DIR = "data"
EDGES_OUT = os.path.join(OUT_DIR, "edges.csv")
NODES_OUT = os.path.join(OUT_DIR, "nodes.csv")

# --- Helpers -----------------------------------------------------------------

HGNC_RE = re.compile(r"\bhgnc:([A-Za-z0-9\-]+)")
UNIPROT_HUMAN_RE = re.compile(r"\buniprotkb:([A-Za-z0-9]+_HUMAN)\b")
ENSEMBL_RE = re.compile(r"\bensembl:(ENSG[0-9]+)")
SYMBOL_FROM_UNIPROT_RE = re.compile(r"^([A-Za-z0-9\-]+)_HUMAN$")

def first_match(pattern, text):
    """Return the first regex group match from text, or None if no match."""
    if not isinstance(text, str):
        return None
    m = pattern.search(text)
    return m.group(1) if m else None

def choose_symbol(row_alias, row_altid):
    """
    Pick a stable gene/protein symbol from MITAB alias/alt id fields.

    Priority:
      1) hgnc:SYMBOL from alias
      2) UniProt NAME_HUMAN from alias → SYMBOL
      3) Ensembl ENSG... from alt_identifier
      4) any 'hgnc:SYMBOL' appearing in alt_identifier
    """
    # 1) HGNC in alias
    s = first_match(HGNC_RE, row_alias)
    if s:
        return s

    # 2) UniProt *_HUMAN in alias
    u = first_match(UNIPROT_HUMAN_RE, row_alias)
    if u:
        m = SYMBOL_FROM_UNIPROT_RE.match(u)
        if m:
            return m.group(1)

    # 3) Ensembl in alt_identifier
    e = first_match(ENSEMBL_RE, row_altid)
    if e:
        return e

    # 4) HGNC in alt_identifier as last resort
    s2 = first_match(HGNC_RE, row_altid)
    if s2:
        return s2

    # If nothing found, return None; caller will drop the row
    return None

def parse_confidence(conf_value):
    """
    Parse MITAB confidence strings into a single numeric weight.

    MITAB 'confidence_score' often like: 'lpr:5|hpr:5|np:1|'
    We’ll take mean of lpr and hpr if present; else fallback to any numeric found.
    """
    if not isinstance(conf_value, str) or not conf_value:
        return None
    nums = {}
    for part in conf_value.strip("|").split("|"):
        if ":" in part:
            key, val = part.split(":", 1)
            try:
                nums[key.strip()] = float(val.strip())
            except:
                pass
    if "lpr" in nums and "hpr" in nums:
        return 0.5 * (nums["lpr"] + nums["hpr"])
    # fallback: first numeric
    for v in nums.values():
        return v
    return None

CYTOKINE_RE = re.compile(r"^(IL[0-9]+|IL-\d+|TNF|IFN[\-\w]*|CSF\d*|CXCL\d+|CCL\d+)$", re.IGNORECASE)
RECEPTOR_RE = re.compile(r"(?:^TLR\d+$|^CD\d+$|^TNFR|^IFNAR|^IFNGR|^IL\d?R|^CCR\d+|^CXCR\d+)", re.IGNORECASE)

def infer_type(symbol):
    """Infer a coarse node type (cytokine/receptor/protein) from a symbol."""
    if symbol is None:
        return "protein"
    sym = symbol.upper()
    if CYTOKINE_RE.match(sym):
        return "cytokine"
    if RECEPTOR_RE.search(sym):
        return "receptor"
    return "protein"

# --- Load and convert ---------------------------------------------------------

def main():
    """Convert MITAB interactions into edges/nodes CSVs for analysis."""
    os.makedirs(OUT_DIR, exist_ok=True)

    usecols = [
        "alias_A","alias_B",
        "alt_identifier_A","alt_identifier_B",
        "interaction_type","confidence_score",
        "ncbi_taxid_A","ncbi_taxid_B"
    ]
    df = pd.read_csv(MITAB_PATH, sep="\t", dtype=str, usecols=lambda c: True)
    # Some dumps have different exact headers; ensure aliases exist
    for col in usecols:
        if col not in df.columns:
            df[col] = None

    # Keep human-human only (taxid:9606)
    def is_human(x):
        """Return True if the taxid field indicates a human interactor."""
        return isinstance(x, str) and ("taxid:9606" in x or "9606(Human)" in x)
    df = df[ df["ncbi_taxid_A"].apply(is_human) & df["ncbi_taxid_B"].apply(is_human) ].copy()

    # Extract clean symbols
    df["sym_A"] = df.apply(lambda r: choose_symbol(r.get("alias_A"), r.get("alt_identifier_A")), axis=1)
    df["sym_B"] = df.apply(lambda r: choose_symbol(r.get("alias_B"), r.get("alt_identifier_B")), axis=1)

    # Drop rows we cannot name
    df = df[df["sym_A"].notna() & df["sym_B"].notna()].copy()

    # Remove self-loops
    df = df[df["sym_A"] != df["sym_B"]].copy()

    # Weight from confidence
    df["weight"] = df["confidence_score"].apply(parse_confidence)

    # Interaction label (short)
    def short_inter_type(s):
        """Extract a short interaction label from the MITAB interaction type."""
        if not isinstance(s, str):
            return None
        # e.g. psi-mi:"MI:0915"(physical association) → physical association
        m = re.search(r"\(([^)]+)\)\s*$", s)
        return m.group(1) if m else s
    df["interaction"] = df["interaction_type"].apply(short_inter_type)

    # Deduplicate (undirected – sort endpoints)
    endpoints = df[["sym_A","sym_B"]].values
    undirected_keys = [tuple(sorted(pair)) for pair in endpoints]
    df["pair"] = undirected_keys
    df = df.sort_values(["pair","weight"], ascending=[True, False]).drop_duplicates("pair", keep="first")

    # Build edges.csv
    edges = pd.DataFrame({
        "source": df["sym_A"].values,
        "target": df["sym_B"].values,
        "weight": df["weight"].values,
        "interaction": df["interaction"].values
    })
    edges.to_csv(EDGES_OUT, index=False)

    # Build nodes.csv with simple type inference
    all_syms = pd.unique(edges[["source","target"]].values.ravel("K"))
    nodes = pd.DataFrame({"id": all_syms})
    nodes["type"] = nodes["id"].apply(infer_type)
    nodes.to_csv(NODES_OUT, index=False)
    
    # --- define new higher-level nodes ---------------------------------
    new_nodes = [
        {"id": "Macrophage",     "type": "cell"},
        {"id": "Dendritic_cell", "type": "cell"},
        {"id": "T_cell",         "type": "cell"},
        {"id": "B_cell",         "type": "cell"},
        {"id": "Pathogen",       "type": "pathogen"},
    ]
    
    # Only keep those that aren't already in the graph
    existing_ids = set(nodes["id"].astype(str))
    new_nodes = [n for n in new_nodes if n["id"] not in existing_ids]
    
    nodes = pd.concat([nodes, pd.DataFrame(new_nodes)], ignore_index=True)
    
    edges.to_csv("data/edges.csv", index=False)
    nodes.to_csv("data/nodes.csv", index=False)
    
    print("Wrote data/edges.csv and data/nodes.csv")

if __name__ == "__main__":
    main()
