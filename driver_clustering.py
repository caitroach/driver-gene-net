'''
do known driver genes cluster in the PPI network? 
take the driver genes D (subset of V) and look at the induced subgraph G[D]
which is the graph you get if u keep only driver nodes and only edges where both endpoints are drivers. 

the number of edges n = |E(G[D])| is how many driver-driver interactions exist. 
if drivers were scattered at random across the PPI network, n will be very small

this script basically 
- samples a null set with the same degree profile as D
- bins every node in V by degree 
- for each driver draw a random replacement from its own bin 
- result is fake driver set with same hubs as real one 
- SO if the real drivers still beat this null then clustering is more than a hub artifact
'''

import os 
import random 
import networkx as nx 
import numpy as np 
import pandas as pd 

DATA_DIR = "data"
EDGE_FILE = os.path.join(DATA_DIR, "ppi_edges.tsv") # u get this from running build_network.py, do that first

N_PERMUTATIONS = 1000 # somewhat arbitrary see docstring 
N_DEGREE_BINS = 100 
SEED = 0 


FALLBACK_DRIVERS = """
TP53 KRAS PIK3CA PTEN APC EGFR BRAF NRAS RB1 CDKN2A NF1 ARID1A KMT2D KMT2C
SMAD4 ATM BRCA1 BRCA2 MYC CTNNB1 FBXW7 IDH1 IDH2 VHL MET ERBB2 ALK NOTCH1
FAT1 CREBBP EP300 SETD2 BAP1 PBRM1 STK11 KEAP1 NFE2L2 CDH1 GATA3 MAP3K1
RUNX1 TET2 DNMT3A ASXL1 JAK2 FLT3 NPM1 CEBPA WT1 SF3B1 U2AF1 SRSF2 CBL EZH2
RAD21 STAG2 SMARCA4 SMARCB1 TSC1 TSC2 MTOR AKT1 CCND1 CDK4 MDM2 ATRX MEN1
RET KIT PDGFRA ABL1 POLE MSH2 MLH1 MSH6 PMS2 PPP2R1A ERBB3 FGFR2 FGFR3 HRAS
RHOA TGFBR2 CASP8 B2M NFKBIA TRAF3 CYLD SOX2 CCNE1 CHEK2 PALB2 RAD51C BARD1
""".split() 
# in case the driver file fails, these are the most famous cancer drivers IN GENERAL...
# my family carries one of the mutations on this list !!! so weird to see those letters here in a different context haha 

def load_graph(): 
    G = nx.Graph() 
    with open(EDGE_FILE) as fh: 
        next(fh) # header 
        for line in fh: 
            u,v = line.rstrip("\n").split("\t")
            G.add_edge(u,v)
    return G 


# # here we grab the data and just get a list of drivers for ovarian epithelial tumors 
def load_drivers(G):
    try:
        df = pd.read_csv("data/IntOGen-DriverGenes_OVT.tsv", sep="\t")
        drivers = df.iloc[:,0].tolist()
        print(f"drivers loaded.\n{len(drivers)} mutational cancer genes for ovarian epithelial tumors (OVT) via IntOGen:")
        print("---------------------------------------------------------------------------------------------")
        for driver in drivers: 
            print(driver,end=" ")
        print("\n---------------------------------------------------------------------------------------------")
    except FileNotFoundError: 
        print("file not found. using fallback drivers.")
        drivers = FALLBACK_DRIVERS
        print(f"fallback drivers loaded. {len(drivers)} well-documented general mutational cancer genes:")
        print("---------------------------------------------------------------------------------------------")
        for driver in drivers: 
            print(driver,end=" ")
        print("\n---------------------------------------------------------------------------------------------")

    present = [g for g in drivers if g in G]
    print(f"{len(present)}/{len(drivers)} present in the network")
    return set(present)
# yes this is inefficient... don't @ me i don't care right now it's not a priority


# returns edge count, largest component size of induced subgraph 
def clustering_stats(G, node_set):
    H = G.subgraph(node_set)
    n_edges = H.number_of_edges()
    if n_edges == 0: 
        return 0,1 if len(node_set) else 0 
    lcc = max(len(c) for c in nx.connected_components(H))
    return n_edges, lcc 

def uniform_null(G, k, rng):
    return rng.sample(list(G.nodes()), k)

 
# partition nodes into equal size bins ordered by degree. size is better than width because that''d put most genes in bin 0 and leave hub bins mostly empty 
def build_degree_bins(G, n_bins=N_DEGREE_BINS): 
    nodes_by_degree = sorted(G.degree(), key=lambda kv: kv[1])
    per_bin = max(1, len(nodes_by_degree)//n_bins)
    bin_of_node, bins = {}, [] 
    for i in range(0, len(nodes_by_degree), per_bin):
        chunk = [g for g,_ in nodes_by_degree[i:i+per_bin]]
        idx = len(bins)
        bins.append(chunk)
        for g in chunk: 
            bin_of_node[g] = idx 
    return bin_of_node, bins 

def degree_matched_null(drivers, bin_of_node, bins, rng): 
    needed = {} 
    for g in drivers: 
        needed[bin_of_node[g]] = needed.get(bin_of_node[g], 0)+1 

    sampled = [] 
    for bin_idx, count in needed.items(): 
        pool = bins[bin_idx]
        sampled.extend(rng.sample(pool, min(count, len(pool))))
    return sampled 

# one sided p value with +1 so p is never exactly 0 
def summarize(label, observed, null_values):
    null_values = np.array(null_values, dtype=float)
    n_ge = int((null_values >= observed).sum())
    p = (n_ge + 1) / (len(null_values) + 1)
    mu, sd = null_values.mean(), null_values.std() 
    z = (observed - mu)/sd if sd > 0 else float("nan")

    print(f"\t\t{label:<22} null mean {mu:8.1f}+/-{sd:6.1f}"
          f"fold {observed/mu if mu else float('inf'):5.2f}x"
          f"z = {z:6.2f} p = {p:.4f}")

def main(): 
    rng = random.Random(SEED)
    print("Loading network...")
    G = load_graph() 
    print(f"{G.number_of_nodes():,} nodes, {G.number_of_edges():,} edges")
    print("\nLoading drivers...")
    drivers = load_drivers(G) 
    k = len(drivers)

    obs_edges, obs_lcc = clustering_stats(G, drivers)
    driver_degrees = [G.degree(g) for g in drivers]
    all_degrees = [d for _, d in G.degree()]
    print(f"\nObserved:")
    print(f"largest connected blob:\t\t {obs_lcc}/{k} drivers")
    print(f"driver-driver edges:\t\t {obs_edges}")
    print(f"mean driver degree:\t\t{np.mean(driver_degrees):.1f}")
    print(f"network-wide mean:\t\t {np.mean(all_degrees):.1f}")

    print(f"\nRunning {N_PERMUTATIONS} permutations per null...")

    bin_of_node, bins = build_degree_bins(G)
    null_uniform_edges, null_uniform_lcc = [], []
    null_matched_edges, null_matched_lcc = [],[] 

    for _ in range(N_PERMUTATIONS):
        e, l = clustering_stats(G, uniform_null(G,k,rng))
        null_uniform_edges.append(e)
        null_uniform_lcc.append(l)

        e, l = clustering_stats(G, degree_matched_null(drivers, bin_of_node, bins, rng))
        null_matched_edges.append(e)
        null_matched_lcc.append(l)

    print(f"\nedges in induced subgraph (observed={obs_edges}):")
    summarize("vs uniform", obs_lcc, null_uniform_lcc)
    summarize("vs degree-matched", obs_lcc, null_matched_lcc)

    print("""
    interpreting this: 
    - big effect vs uniform, big effect vs degree-matched: drivers genuinely concentrate in shared network neighbourhoods
    - big effect vs uniform, effect gone vs degree-matched: driver clusteting is mostly the hub effect
    """)

if __name__ == "__main__":
    main()