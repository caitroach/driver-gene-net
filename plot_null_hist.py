import random
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

import driver_clustering as dc

# colours: blue / orange / aqua = first three slots of the validated categorical palette
BLUE, ORANGE, AQUA = "#2a78d6", "#eb6834", "#1baf7a"
INK, INK2, GRID = "#0b0b0b", "#52514e", "#e6e5e1"

rng = random.Random(dc.SEED)
G = dc.load_graph()
drivers = dc.load_drivers(G)
k = len(drivers)
obs_edges, obs_lcc = dc.clustering_stats(G, drivers)
bin_of_node, bins = dc.build_degree_bins(G)

uni_edges, mat_edges, mat_mean_deg = [], [], []
for _ in range(dc.N_PERMUTATIONS):
    s = dc.uniform_null(G, k, rng)
    uni_edges.append(dc.clustering_stats(G, s)[0])
    s = dc.degree_matched_null(drivers, bin_of_node, bins, rng)
    mat_edges.append(dc.clustering_stats(G, s)[0])
    mat_mean_deg.append(np.mean([G.degree(g) for g in s]))

uni_edges, mat_edges = np.array(uni_edges), np.array(mat_edges)
mu, sd = mat_edges.mean(), mat_edges.std()
p = (int((mat_edges >= obs_edges).sum()) + 1) / (len(mat_edges) + 1)

real_deg = np.mean([G.degree(g) for g in drivers])
print(f"\nobserved edges {obs_edges} | degree-matched null {mu:.1f} +/- {sd:.1f} | p = {p:.4f}")
print(f"[check 1] mean degree: real drivers {real_deg:.1f} vs degree-matched null "
      f"{np.mean(mat_mean_deg):.1f}  (should be close; a big gap means the null is under-matched)")
print(f"[check 2] real driver degrees, top 5: "
      f"{sorted([G.degree(g) for g in drivers], reverse=True)[:5]}")
print(f"[check 3] widest degree bin spans degrees "
      f"{min(G.degree(g) for g in bins[-1])}-{max(G.degree(g) for g in bins[-1])} "
      f"(top bin is wide, so matching is loose for the biggest hubs)")

plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 15})
fig, ax = plt.subplots(figsize=(10, 5.6), dpi=200)
xmax = max(obs_edges, mat_edges.max(), uni_edges.max()) + 6
bin_edges = np.arange(-0.5, xmax + 1.5, 1)

ax.hist(uni_edges, bins=bin_edges, color=AQUA, alpha=0.55, edgecolor="white", linewidth=1,
        label="random gene sets")
ax.hist(mat_edges, bins=bin_edges, color=BLUE, alpha=0.9, edgecolor="white", linewidth=1,
        label="degree-matched random sets")
ax.axvline(obs_edges, color=ORANGE, linewidth=3)

ptxt = "p ≤ 0.001" if p <= 0.0011 else f"p = {p:.3f}"
ymax = ax.get_ylim()[1]
ax.text(obs_edges - 0.8, ymax * 0.93, f"real drivers\n{obs_edges} edges", color=INK,
        ha="right", va="top", fontsize=16, fontweight="bold")
ax.text(mu + 7.5, ymax * 0.62, f"degree-matched null\n{mu:.1f} ± {sd:.1f}\n"
        f"real set ≈ {obs_edges/mu:.1f}× more\n{ptxt}", color=INK2, ha="left", va="top")

ax.set_xlabel(f"driver–driver edges in a set of {k} genes")
ax.set_ylabel(f"random sets (of {dc.N_PERMUTATIONS:,})")
ax.set_title(f"Real ovarian drivers vs {dc.N_PERMUTATIONS:,} random gene sets", loc="left",
             fontsize=18, fontweight="bold", color=INK)
for s in ("top", "right"):
    ax.spines[s].set_visible(False)
ax.spines["left"].set_color(GRID); ax.spines["bottom"].set_color(INK2)
ax.tick_params(colors=INK2)
ax.yaxis.grid(True, color=GRID, linewidth=1); ax.set_axisbelow(True)
ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.28, 0.98), fontsize=13)
fig.tight_layout()
fig.savefig("figures/null_histogram.png", facecolor="white")
print("saved figures/null_histogram.png")
