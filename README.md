# RepeatSearchTools
This repo contains tools allowing for search of imperfect repeats (Inverted and Direct) in a DNA sequence.

## Refactored Structure

Core logic is now split into dedicated utility modules:

- `utils/repeat_parsing.py`: sequence parsing, repeat-search logic, and positional filtering
- `utils/repeat_visualization.py`: plotting and repeat-map visualization helpers
- `RepeatSearch.ipynb`: workflow notebook that imports and uses the utility modules

All utility functions include Python type hints and Google-style docstrings.

## Requirements

Install dependencies with:

```bash
pip install regex matplotlib numpy
```

`json` is part of the Python standard library and does not need installation.

## Method Overview

![](https://github.com/tvarovski/RepeatSearchTools/blob/main/method_overview.png)
