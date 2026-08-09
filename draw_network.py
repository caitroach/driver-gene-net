"""
03_plot_driver_subgraph.py
==========================

Draw the induced driver subgraph G[D] -- the figure that goes at the top of
your README.

THE LAYOUT PROBLEM
------------------
You have 49 driver nodes: 33 form one connected blob, a handful form small
pairs/triples, and the rest are isolated (no driver-driver edges at all).

If you hand all of that to nx.spring_layout(), you get a mess. Spring layout is
a physics simulation -- edges act as springs pulling nodes together, and all
nodes repel each other. Isolated nodes have no springs attached, so nothing
pulls them anywhere; they just get shoved to the outer rim by repulsion and
sprayed around the border in a meaningless ring. Worse, they inflate the
bounding box, squashing the part you actually care about into the middle.

So we lay out ONLY the connected part with springs, and place the isolated
genes deliberately in a labelled strip underneath. This is honest -- the
isolates are a real part of the result (16 of your drivers have no driver
neighbour) and hiding them would overstate the clustering.

VISUAL ENCODING
---------------
    node size   = degree in the FULL network (how hubby the gene is overall)
    node colour = degree within the driver subgraph (how many driver partners)
    edge width  = fixed, thin

Size and colour deliberately encode different things. A big pale node is a
famous hub that happens not to touch other drivers; a small dark node is an
unremarkable gene sitting in a driver-dense pocket. That contrast is the
degree-vs-clustering story your permutation test is about, made visible.

Run:  python 03_plot_driver_subgraph.py
Out:  figures/driver_subgraph.png  (and .svg for editing)
"""

import os

import matplotlib.pyplot as plt
import networkx as nx
from matplotlib.lines import Line2D

DATA_DIR = "data"
FIG_DIR = "figures"
EDGE_FILE = os.path.join(DATA_DIR, "ppi_edges.tsv")
DRIVER_FILE = os.path.join(DATA_DIR, "ovtdrivers.txt")

CANCER_LABEL = "ovarian epithelial tumours (IntOGen OVT)"
SEED = 5  # change this if the layout looks cramped -- it's just the RNG start


def load_graph():
    G = nx.Graph()
    with open(EDGE_FILE) as fh:
        next(fh)
        for line in fh:
            u, v = line.rstrip("\n").split("\t")
            G.add_edge(u, v)
    return G


def load_drivers(G):
    with open(DRIVER_FILE) as fh:
        genes = {ln.strip() for ln in fh if ln.strip() and not ln.startswith("#")}
    return {g for g in genes if g in G}


def build_layout(H, seed=SEED):
    """
    Position nodes: springs for the connected components, a tidy grid for isolates.

    Returns (pos dict, list of isolate names).
    """
    isolates = sorted(n for n in H.nodes() if H.degree(n) == 0)
    connected = [n for n in H.nodes() if H.degree(n) > 0]

    sub = H.subgraph(connected)

    # k controls the target distance between nodes. The default is
    # 1/sqrt(n), which packs labelled nodes far too tightly to read --
    # roughly tripling it gives labels room to breathe.
    pos = nx.spring_layout(sub, k=2.5 / (len(connected) ** 0.5),
                           iterations=300, seed=seed)

    # Normalise the connected part into y in [0, 1] so the isolate strip
    # below it always sits clear of the blob, whatever the layout produced.
    if pos:
        ys = [p[1] for p in pos.values()]
        xs = [p[0] for p in pos.values()]
        y_min, y_span = min(ys), (max(ys) - min(ys)) or 1
        x_min, x_span = min(xs), (max(xs) - min(xs)) or 1
        pos = {n: ((x - x_min) / x_span, (y - y_min) / y_span * 0.8 + 0.2)
               for n, (x, y) in pos.items()}

    # Isolates: even row(s) across the bottom, on a FIXED column grid so that
    # a short final row lines up under the row above instead of stretching to
    # fill the width (which reads as meaningful spacing when it isn't).
    per_row = 9
    for i, gene in enumerate(isolates):
        row, col = divmod(i, per_row)
        pos[gene] = ((col + 0.5) / per_row, 0.08 - row * 0.09)

    return pos, isolates


def draw(G, H, pos, isolates, out_stem):
    fig, ax = plt.subplots(figsize=(14, 12))

    # --- encodings ---------------------------------------------------------
    full_degree = {n: G.degree(n) for n in H.nodes()}
    sub_degree = {n: H.degree(n) for n in H.nodes()}

    nodes = list(H.nodes())
    sizes = [120 + 12 * full_degree[n] for n in nodes]
    colours = [sub_degree[n] for n in nodes]

    nx.draw_networkx_edges(H, pos, ax=ax, edge_color="#9aa5b1", width=1.2, alpha=0.8)

    drawn = nx.draw_networkx_nodes(
        H, pos, ax=ax, nodelist=nodes,
        node_size=sizes, node_color=colours,
        cmap=plt.cm.viridis, vmin=0, edgecolors="white", linewidths=1.2,
    )

    # Labels sit just BELOW each node rather than on top of it. Centred labels
    # are unreadable against the dark end of viridis, and long gene symbols
    # spill over the node edge. The offset scales with node size so big hubs
    # push their label further out.
    for n in nodes:
        x, y = pos[n]
        offset = 0.018 + 0.00016 * full_degree[n]
        ax.text(x, y - offset, n,
                ha="center", va="top", fontsize=8.5, fontweight="bold",
                color="#1f2933", zorder=5,
                bbox=dict(boxstyle="round,pad=0.15", facecolor="white",
                          edgecolor="none", alpha=0.75))

    # --- annotation --------------------------------------------------------
    if isolates:
        ax.axhline(0.16, color="#c0c8d0", linestyle="--", linewidth=1)
        ax.text(0.5, 0.175, f"no driver-driver edge  (n = {len(isolates)})",
                ha="center", va="bottom", fontsize=10,
                color="#5a6470", style="italic", transform=ax.transData)

    cbar = fig.colorbar(drawn, ax=ax, fraction=0.03, pad=0.02)
    cbar.set_label("driver neighbours (degree within subgraph)", fontsize=10)

    # Legend for node size -- a colourbar handles colour, but size needs its own
    # key or readers will guess wrong about what big means.
    size_key = [
        Line2D([], [], marker="o", linestyle="", markerfacecolor="#b9c2cc",
               markeredgecolor="white", markersize=(120 + 12 * d) ** 0.5 / 2.4,
               label=f"{d} PPI partners")
        for d in (10, 50, 150)
    ]
    ax.legend(handles=size_key, loc="upper left", frameon=False,
              fontsize=9, labelspacing=1.4, title="node size",
              title_fontsize=9, borderpad=1)

    n_edges = H.number_of_edges()
    lcc = max((len(c) for c in nx.connected_components(H)), default=0)
    ax.set_title(
        f"Driver genes in the human protein interaction network\n"
        f"{CANCER_LABEL}  |  {H.number_of_nodes()} drivers, {n_edges} "
        f"driver-driver edges, largest component {lcc}",
        fontsize=13, pad=18,
    )

    ax.set_axis_off()
    ax.margins(0.08)
    fig.tight_layout()

    os.makedirs(FIG_DIR, exist_ok=True)
    for ext in ("png", "svg"):
        path = f"{out_stem}.{ext}"
        fig.savefig(path, dpi=200, bbox_inches="tight",
                    facecolor="white")
        print(f"  wrote {path}")
    plt.close(fig)


def main():
    G = load_graph()
    drivers = load_drivers(G)
    H = G.subgraph(drivers).copy()

    print(f"drivers in network : {H.number_of_nodes()}")
    print(f"driver-driver edges: {H.number_of_edges()}")

    pos, isolates = build_layout(H)
    print(f"isolated drivers   : {len(isolates)}  ({', '.join(isolates)})")

    draw(G, H, pos, isolates, os.path.join(FIG_DIR, "driver_subgraph"))


if __name__ == "__main__":
    main()