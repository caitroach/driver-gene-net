'''
builds the protein-protein interaction (PPI) graph that everything else runs on
gets a simple undirected unweighted graph G=(V,E) where 
- V: human gene
- E = {u,v} if proteins made by u and v interact
also translates STRING ENSEMBL protein IDs to gene symbols
because our mutation data is keyed by gene symbol 
also only keeping high-confidence edges from STRING
+ only using the physical subnetwork to reduce noise
'''

import gzip 
import os 
import urllib.request 

import networkx as nx 
import pandas as pd 

DATA_DIR = "data"
SCORE_CUTOFF = 700 # high confidence from STRING, may change l8r

# grab data
STRING_BASE = "https://stringdb-downloads.org/download"
LINKS_URL = f"{STRING_BASE}/protein.physical.links.v12.0/9606.protein.physical.links.v12.0.txt.gz"
INFO_URL = f"{STRING_BASE}/protein.info.v12.0/9606.protein.info.v12.0.txt.gz"

LINKS_FILE = os.path.join(DATA_DIR, "9606.protein.physical.links.v12.0.txt.gz")
INFO_FILE = os.path.join(DATA_DIR, "9606.protein.info.v12.0.txt.gz")
EDGE_OUT = os.path.join(DATA_DIR, "ppi_edges.tsv")

# fetches file once and caches it, never redownloads STRING file bcoz static
def download_ifmissing(url, path):
    if os.path.exists(path):
        print(f"[cached] {os.path.basename(path)}")
        return
    print(f"[downloading] {os.path.basename(path)}")
    urllib.request.urlretrieve(url,path)
    return

# STRING protein ID -> gene symbpl
def load_id_map():
    with gzip.open(INFO_FILE, "rt") as fh:
        info = pd.read_csv(fh, sep="\t")
    info.columns = ["string_id", "symbol", "size", "annotation"]
    return dict(zip(info["string_id"], info["symbol"]))

def build_graph(): 
    id2sym = load_id_map()
    print(f"mapped {len(id2sym):,} STRING protein IDs to gene symbols")

    # links file is space separated
    with gzip.open(LINKS_FILE, "rt") as fh:
        links = pd.read_csv(fh, sep=" ")
    print(f" {len(links):,} raw physical-interaction edges")

    # filters by confidence 
    links = links[links["combined_score"] >= SCORE_CUTOFF]
    print(f"{len(links):,} edges at score >= {SCORE_CUTOFF}")

    # IDs -> symbols 
    links["g1"] = links["protein1"].map(id2sym)
    links["g2"] = links["protein2"].map(id2sym)
    links = links.dropna(subset=["g1", "g2"])

    # networkx.Graph collapses duplicates but im dropping self loops explicitly
    # so the degrees dont get messed up 
    links = links[links["g1"] != links["g2"]]
    G = nx.Graph()
    G.add_edges_from(zip(links["g1"], links["g2"]))
    print(f"graph: {G.number_of_nodes():,} nodes, {G.number_of_edges():,} edges")

    # keeping giant component
    components = sorted(nx.connected_components(G), key=len, reverse=True)
    print(f"{len(components)} connected components," 
          f"largest holds {len(components[0]):,} nodes"
          f"({len(components[0]) / G.number_of_nodes():.1%} of all nodes)")
    G = G.subgraph(components[0]).copy() # giant

    return G

# sanity check, if these look wrong then stop
def describe(G): 
    degrees = [d for _, d in G.degree()]
    print("\nNetwork summary:")
    print(f"nodes:\t\t {G.number_of_nodes():,}")
    print(f"edges:\t\t {G.number_of_edges():,}")
    print(f"mean degree\t\t: {sum(degrees) / len(degrees):.2f}")
    print(f"median degree\t\t: {sorted(degrees)[len(degrees) // 2]}")
    print(f"max degree\t\t: {max(degrees)}")
    print(f"density\t\t: {nx.density(G):.6f}")

    hubs = sorted(G.degree(), key=lambda kv: kv[1], reverse=True)[:10]
    print(f"top hubs\t\t: {', '.join(f'{g}({d})' for g, d in hubs)}")
    return

def main():
    os.makedirs(DATA_DIR, exist_ok=True)
 
    print("Fetching STRING v12.0 ...")
    download_ifmissing(LINKS_URL, LINKS_FILE)
    download_ifmissing(INFO_URL, INFO_FILE)
 
    print("\nBuilding graph ...")
    G = build_graph()
    describe(G)
 
    with open(EDGE_OUT, "w") as fh:
        fh.write("gene1\tgene2\n")
        for u, v in G.edges():
            fh.write(f"{u}\t{v}\n")
    print(f"\nWrote {EDGE_OUT}")
 
 
if __name__ == "__main__":
    main()