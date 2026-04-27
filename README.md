# RepeatSearchTools

Tools for finding imperfect or perfect direct and inverted repeats in DNA sequences.

The notebook is intentionally kept short. Use it as the runnable workflow, and use this README for setup, inputs, parameter notes, and method context.

## Repository Structure

- `RepeatSearch.ipynb`: minimal repeat-search workflow
- `utils/repeat_parsing.py`: FASTA parsing, repeat search, coordinate tracking, and redundancy filtering
- `utils/repeat_visualization.py`: repeat-map plotting helpers
- `method_overview.png`: method overview figure

## Requirements

Install third-party dependencies with:

```bash
pip install regex matplotlib numpy
```

`json` is part of the Python standard library and does not need to be installed.

## Basic Workflow

1. Open `RepeatSearch.ipynb`.
2. Set `sequence` directly, or load a FASTA file with `extractSeqFromFileToList`.
3. Adjust `search_params`.
4. Run the combined search section once to create `repeats`, save `repeat_records.json`, and plot `repeat_visualization.png`.

## Search Parameters

- `repeat_type`: `"direct"`, `"inverted"`, or `"both"`
- `min_query_length`: shortest repeat length to test
- `max_query_length`: longest repeat length to test
- `min_spacer`: minimum distance between the query and matched repeat
- `window_size`: local sequence window used during scanning
- `imperfect_homology`: allow approximate matches when `True`
- `min_homology`: minimum fractional homology for approximate matching
- `fixed_errors`: explicit edit-distance limit; set to `False` to use `min_homology`

For imperfect matches, terminal bases must match exactly. Allowed errors are only accepted inside the repeat, which avoids greedy edge matches. Candidate match spans may be shorter or longer than the query when insertions or deletions are used, and `match_end` records the actual matched span.

Search time is very sensitive to `window_size`. Smaller values run faster because distant repeat pairs are not evaluated.

## Visualization

The notebook runs `searchSequenceForRepeats` once. It returns position-aware repeat records, and those same records are used for both JSON output and visualization.

Each JSON record contains:

- `query_seq` and `match_seq`
- `query_start` and `query_end`
- `match_start` and `match_end`
- `is_perfect`
- `length`
- `repeat_type`
- `inverted`

By default, the plot is saved as `repeat_visualization.png`.

`visualizeRepeats` uses translucent solid arcs and distinct colors for perfect/imperfect direct/inverted repeats. You can tune the figure with optional arguments such as `arc_alpha=0.2` for lighter arcs, `density_cmap="viridis"` for a different density colormap, or `min_interval_length=18` to show baseline interval rectangles only for longer repeats.

Arcs are anchored at each repeat segment midpoint. Heatmap density is center-weighted across each repeat span, so high-activity regions remain visible without treating start and end positions as separate signals.

Example visualization:
![](https://github.com/tvarovski/RepeatSearchTools/blob/main/repeat_visualization.png)

## Method Overview

![](https://github.com/tvarovski/RepeatSearchTools/blob/main/method_overview.png)
