from __future__ import annotations
from typing import Sequence
import matplotlib.gridspec as gridspec
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D
from matplotlib.path import Path
from .repeat_parsing import RepeatRecord

def visualizeRepeats(
    sequence: str,
    repeats: Sequence[RepeatRecord],
    title: str = "Repeat Position Map",
    figsize: tuple[int, int] = (20, 8),
    save_path: str = "repeat_visualization.png",
    density_cmap: str = "magma",
    arc_alpha: float = 0.28,
    min_interval_length: int | None = None,
) -> None:
    """Visualize repeat-pair locations and per-position repeat density.

    Args:
        sequence: Input DNA sequence.
        repeats: Position-aware repeat records.
        title: Figure title.
        figsize: Matplotlib figure size in inches.
        save_path: Output image path.
        density_cmap: Matplotlib colormap used for repeat-density heatmap.
        arc_alpha: Base transparency for repeat arcs.
        min_interval_length: Minimum query length for drawing baseline interval
            rectangles. Defaults to the midpoint of the searched repeat lengths.

    Returns:
        None.
    """
    if not repeats:
        print("No repeats to visualize.")
        return

    seq_len = len(sequence)

    color_map: dict[tuple[bool, bool], str] = {
        (True, True): "#0072B2",
        (False, True): "#D55E00",
        (True, False): "#009E73",
        (False, False): "#CC79A7",
    }
    baseline_y = 0.0
    interval_y = -0.045
    interval_height = 0.09
    arc_base_y = interval_y + interval_height

    lengths = [record["length"] for record in repeats]
    min_len = min(lengths)
    max_len = max(lengths)
    len_range = max(max_len - min_len, 1)
    if min_interval_length is None:
        min_interval_length = round((min_len + max_len) / 2)

    def get_height(length: int) -> float:
        return 0.25 + 0.9 * (length - min_len) / len_range

    def get_center(start: int, end: int) -> float:
        return (start + end) / 2

    def add_center_weighted_density(start: int, end: int) -> None:
        clipped_start = min(max(start, 0), seq_len)
        clipped_end = min(max(end, 0), seq_len)
        span = clipped_end - clipped_start
        if span <= 0:
            return

        positions = np.arange(clipped_start, clipped_end) + 0.5
        center = get_center(start, end)
        half_width = max((end - start) / 2, 0.5)
        weights = 1 - np.abs(positions - center) / half_width
        weights = np.clip(weights, 0, None)
        if weights.sum() > 0:
            weights *= span / weights.sum()
        density[clipped_start:clipped_end] += weights

    density = np.zeros(seq_len)
    for record in repeats:
        add_center_weighted_density(record["query_start"], record["query_end"])
        add_center_weighted_density(record["match_start"], record["match_end"])

    fig = plt.figure(figsize=figsize, facecolor="white", constrained_layout=True)
    gs = gridspec.GridSpec(2, 1, height_ratios=[5, 1.1], hspace=0.08, figure=fig)
    ax_main = fig.add_subplot(gs[0])
    ax_heat = fig.add_subplot(gs[1], sharex=ax_main)
    fig.suptitle(title, fontsize=16, fontweight="bold")
    ax_main.set_facecolor("#FAFAFA")
    ax_heat.set_facecolor("#FAFAFA")

    ax_main.add_patch(
        mpatches.FancyArrow(
            0,
            baseline_y,
            seq_len,
            0,
            width=0.02,
            length_includes_head=True,
            head_width=0.06,
            head_length=max(seq_len * 0.012, 1),
            color="#263238",
            alpha=0.85,
            zorder=8,
        )
    )

    def get_record_style(record: RepeatRecord) -> tuple[str, float, float, float]:
        color = color_map[(record["is_perfect"], record["inverted"])]
        height = get_height(record["length"])
        alpha = arc_alpha if record["is_perfect"] else arc_alpha * 0.75
        linewidth = 0.5 + 1.1 * (record["length"] - min_len) / len_range
        return color, height, alpha, linewidth

    interval_segments = []
    for record in repeats:
        if record["length"] < min_interval_length:
            continue
        color, _, alpha, _ = get_record_style(record)
        interval_segments.extend(
            [
                (record["query_start"], record["query_end"], color, alpha),
                (record["match_start"], record["match_end"], color, alpha),
            ]
        )

    for segment_start, segment_end, color, alpha in sorted(
        interval_segments,
        key=lambda segment: segment[1] - segment[0],
    ):
        ax_main.add_patch(
            mpatches.Rectangle(
                (segment_start, interval_y),
                segment_end - segment_start,
                interval_height,
                color=color,
                alpha=min(alpha + 0.2, 0.7),
                zorder=5,
                linewidth=0,
            )
        )

    for record in repeats:
        color, height, alpha, linewidth = get_record_style(record)

        x1 = get_center(record["query_start"], record["query_end"])
        x2 = get_center(record["match_start"], record["match_end"])
        path = Path(
            [(x1, arc_base_y), (x1, height), (x2, height), (x2, arc_base_y)],
            [Path.MOVETO, Path.CURVE4, Path.CURVE4, Path.CURVE4],
        )
        ax_main.add_patch(
            mpatches.PathPatch(
                path,
                facecolor="none",
                edgecolor=color,
                linestyle="-",
                alpha=alpha,
                lw=linewidth,
                capstyle="round",
                joinstyle="round",
                zorder=6,
            )
        )

    max_height = get_height(max_len)
    ax_main.set_xlim(-seq_len * 0.01, seq_len * 1.03)
    ax_main.set_ylim(-0.25, max_height + 0.3)
    ax_main.set_ylabel("Relative Repeat Length", fontsize=11)
    ax_main.set_yticks([])
    ax_main.grid(axis="x", color="#E0E0E0", linewidth=0.6, alpha=0.7)
    for spine in ["top", "right", "left", "bottom"]:
        ax_main.spines[spine].set_visible(False)
    ax_main.tick_params(axis="x", labelbottom=False)

    legend_handles = [
        Line2D([0], [0], color=color_map[(True, True)], lw=2, label="Perfect Inverted"),
        Line2D([0], [0], color=color_map[(False, True)], lw=2, label="Imperfect Inverted"),
        Line2D([0], [0], color=color_map[(True, False)], lw=2, label="Perfect Direct"),
        Line2D([0], [0], color=color_map[(False, False)], lw=2, label="Imperfect Direct"),
    ]
    ax_main.legend(
        handles=legend_handles,
        loc="upper right",
        fontsize=9,
        frameon=False,
        ncol=2,
    )

    heatmap = ax_heat.imshow(
        density.reshape(1, -1),
        aspect="auto",
        cmap=density_cmap,
        extent=[0, seq_len, 0, 1],
        vmin=0,
        vmax=max(density.max(), 1),
    )
    ax_heat.set_yticks([0.5])
    ax_heat.set_yticklabels(["Density"], fontsize=9)
    ax_heat.set_xlabel("Position in Sequence (bp)", fontsize=11)
    ax_heat.tick_params(axis="both", length=0)
    for spine in ["top", "right", "left", "bottom"]:
        ax_heat.spines[spine].set_visible(False)

    color_bar = plt.colorbar(heatmap, ax=ax_heat, orientation="horizontal", pad=0.35)
    color_bar.set_label("Repeat Density", fontsize=8)
    color_bar.outline.set_visible(False)
    color_bar.ax.tick_params(axis="x", labelsize=8, length=0)
    heat_pos = ax_heat.get_position()
    cbar_width = heat_pos.width * 0.35
    cbar_height = heat_pos.height * 0.22
    cbar_left = heat_pos.x0 + (heat_pos.width - cbar_width) / 2
    cbar_bottom = max(heat_pos.y0 - cbar_height * 1.8, 0.02)
    color_bar.ax.set_position([cbar_left, cbar_bottom, cbar_width, cbar_height])

    plt.savefig(save_path, dpi=600, bbox_inches="tight")
    plt.show()
    print(f"Visualization saved as '{save_path}'")
