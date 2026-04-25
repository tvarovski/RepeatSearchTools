from __future__ import annotations
from typing import Sequence
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.path import Path
from .repeat_parsing import RepeatRecord

def visualizeRepeats(
    sequence: str,
    repeats: Sequence[RepeatRecord],
    title: str = "Repeat Position Map",
    figsize: tuple[int, int] = (20, 8),
    save_path: str = "repeat_visualization.png",
) -> None:
    """Visualize repeat-pair locations and per-position repeat density.

    Args:
        sequence: Input DNA sequence.
        repeats: Position-aware repeat records.
        title: Figure title.
        figsize: Matplotlib figure size in inches.
        save_path: Output image path.

    Returns:
        None.
    """
    if not repeats:
        print("No repeats to visualize.")
        return

    seq_len = len(sequence)

    color_map: dict[tuple[bool, bool], tuple[str, str]] = {
        (True, True): ("#1565C0", "-"),
        (False, True): ("#E65100", "--"),
        (True, False): ("#2E7D32", "-"),
        (False, False): ("#6A1B9A", "--"),
    }

    lengths = [record["length"] for record in repeats]
    min_len = min(lengths)
    max_len = max(lengths)
    len_range = max(max_len - min_len, 1)

    def get_height(length: int) -> float:
        return 0.4 + 0.6 * (length - min_len) / len_range

    density = np.zeros(seq_len)
    for record in repeats:
        density[record["query_start"] : record["query_end"]] += 1
        density[record["match_start"] : record["match_end"]] += 1

    fig = plt.figure(figsize=figsize)
    gs = gridspec.GridSpec(2, 1, height_ratios=[5, 1], hspace=0.05)
    ax_main = fig.add_subplot(gs[0])
    ax_heat = fig.add_subplot(gs[1], sharex=ax_main)
    fig.suptitle(title, fontsize=15, fontweight="bold", y=0.98)

    ax_main.add_patch(
        mpatches.FancyArrow(
            0,
            0,
            seq_len * 0.99,
            0,
            width=0.03,
            length_includes_head=True,
            head_width=0.07,
            head_length=seq_len * 0.01,
            color="#455A64",
            zorder=2,
        )
    )

    for record in repeats:
        color, line_style = color_map[(record["is_perfect"], record["inverted"])]
        height = get_height(record["length"])
        alpha = 0.75 if record["is_perfect"] else 0.55
        linewidth = 0.8 + 1.2 * (record["length"] - min_len) / len_range

        for segment_start in (record["query_start"], record["match_start"]):
            ax_main.add_patch(
                mpatches.Rectangle(
                    (segment_start, -0.06),
                    record["length"],
                    0.12,
                    color=color,
                    alpha=min(alpha + 0.1, 1.0),
                    zorder=3,
                    linewidth=0,
                )
            )

        x1 = record["query_start"] + record["length"] / 2
        x2 = record["match_start"] + record["length"] / 2
        path = Path(
            [(x1, 0.0), (x1, height), (x2, height), (x2, 0.0)],
            [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4],
        )
        ax_main.add_patch(
            mpatches.PathPatch(
                path,
                facecolor="none",
                edgecolor=color,
                linestyle=line_style,
                alpha=alpha,
                lw=linewidth,
                zorder=4,
            )
        )

    max_height = get_height(max_len)
    ax_main.set_xlim(-seq_len * 0.01, seq_len * 1.02)
    ax_main.set_ylim(-0.25, max_height + 0.3)
    ax_main.set_ylabel("Relative Repeat Length", fontsize=11)
    ax_main.set_yticks([])
    for spine in ["top", "right", "left"]:
        ax_main.spines[spine].set_visible(False)
    ax_main.tick_params(axis="x", labelbottom=False)

    legend_handles = [
        mpatches.Patch(color="#1565C0", label="Perfect Inverted"),
        mpatches.Patch(color="#E65100", label="Imperfect Inverted"),
        mpatches.Patch(color="#2E7D32", label="Perfect Direct"),
        mpatches.Patch(color="#6A1B9A", label="Imperfect Direct"),
    ]
    ax_main.legend(
        handles=legend_handles,
        loc="upper right",
        fontsize=9,
        framealpha=0.85,
        ncol=2,
    )

    heatmap = ax_heat.imshow(
        density.reshape(1, -1),
        aspect="auto",
        cmap="YlOrRd",
        extent=[0, seq_len, 0, 1],
        vmin=0,
        vmax=max(density.max(), 1),
    )
    ax_heat.set_yticks([0.5])
    ax_heat.set_yticklabels(["Density"], fontsize=9)
    ax_heat.set_xlabel("Position in Sequence (bp)", fontsize=11)

    color_bar = plt.colorbar(
        heatmap,
        ax=ax_heat,
        orientation="vertical",
        pad=0.01,
        fraction=0.02,
        shrink=0.8,
    )
    color_bar.set_label("Repeat\nCoverage", fontsize=8)

    plt.savefig(save_path, dpi=150, bbox_inches="tight")
    plt.show()
    print(f"Visualization saved as '{save_path}'")
