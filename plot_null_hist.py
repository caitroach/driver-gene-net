import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

BLUE, ORANGE = "#2a78d6", "#eb6834"
INK, INK2, GRID, FAINT = "#0b0b0b", "#52514e", "#e6e5e1", "#bdbcb7"


def plot_null(uni_edges, mat_edges, obs_edges, k, n_perm, p, out="figures/null_histogram.png"):
    uni_edges, mat_edges = np.asarray(uni_edges), np.asarray(mat_edges)
    mu, sd = mat_edges.mean(), mat_edges.std()

    plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 15})
    fig, ax = plt.subplots(figsize=(10, 5.6), dpi=200)

    xmax = max(obs_edges, mat_edges.max(), uni_edges.max()) + 6
    bins = np.arange(-0.5, xmax + 1.5, 1)

    # fair null = hero; naive null = faint gray context
    ax.hist(uni_edges, bins=bins, color=FAINT, edgecolor="white", linewidth=1,
            label="random gene sets (ignores degree)")
    ax.hist(mat_edges, bins=bins, color=BLUE, edgecolor="white", linewidth=1,
            label="degree-matched random sets")
    ax.axvline(obs_edges, color=ORANGE, linewidth=3)

    ptxt = "p < 0.001" if p <= 0.0011 else f"p = {p:.3f}"
    ymax = ax.get_ylim()[1]

    # label sits right next to the thing it describes
    ax.text(obs_edges + 0.8, ymax * 0.93,
            f"real drivers\n{obs_edges} edges\n{obs_edges/mu:.1f}× the null\n{ptxt}",
            color=INK, ha="left", va="top", fontsize=15, fontweight="bold")
    ax.text(mu + sd + 5, ymax * 0.42, f"degree-matched null\n{mu:.1f} ± {sd:.1f}",
            color=BLUE, fontsize=14, fontweight="bold", va="center")

    ax.set_xlim(-0.5, obs_edges + 14)
    ax.set_xlabel(f"driver–driver edges in a set of {k} genes")
    ax.set_ylabel(f"random sets (of {n_perm:,})")
    ax.set_title(f"Real ovarian drivers vs {n_perm:,} random gene sets",
                 loc="left", fontsize=18, fontweight="bold", color=INK)
    for s in ("top", "right"):
        ax.spines[s].set_visible(False)
    ax.spines["left"].set_color(GRID); ax.spines["bottom"].set_color(INK2)
    ax.tick_params(colors=INK2)
    ax.yaxis.grid(True, color=GRID, linewidth=1); ax.set_axisbelow(True)
    ax.legend(frameon=False, loc="upper left", bbox_to_anchor=(0.2, 1.0), fontsize=13)
    fig.tight_layout()
    fig.savefig(out, facecolor="white")
    plt.close(fig)
    print("saved", out)


if __name__ == "__main__":
    # layout test with fake data shaped like the 49-edge run
    rng = np.random.default_rng(0)
    plot_null(rng.normal(2.5, 1.6, 1000).round().clip(0), rng.normal(12.6, 3.4, 1000).round().clip(0),
              49, 49, 1000, 1 / 1001, "test.png")