# 🧬 EI Genetic Algorithm

## Project Overview

This code was developed as part of a **one-week project** in collaboration with **Université d’Évry Val d’Essonne** and **Université Paris-Saclay**, focused on **optimization in a biological context**.

The main research question was:
> *How can we optimize a model originally designed for short segments of naked DNA so that it also accounts for the properties of long DNA chains, particularly their circular structure?*

The project leverages **optimization algorithms** (including simulated annealing and genetic algorithms) to adjust a DNA conformation model.

---

## Installation

A Python environment is required. This repository was tested with Python 3.12, but it should work with most Python 3 versions.

- `pip install -r requirements.txt` installs **all required packages**.

## Usage

Everything is launched from a terminal command: `python -m dna`  
This command will build and display the DNA strand extracted from _data/plasmid_8k.fasta_ using the 1993 conformation model: _dna/table.json_.

Several parameters allow you to change this behavior. They can be listed using `python -m dna --help`, and are detailed here:

- `-m [traditional | recuit | genetic] (default: traditional)` allows you to choose the operating mode. See the "Modes" section for more details.
- `-d [path_to_file] (default: data/plasmid_8k.fasta)` lets you choose a DNA sequence. Not needed for the `recuit` mode, which automatically trains on both sequences.
- `-j [path_to_file] (default: dna/table.json)` lets you choose a conformation model. In `traditional` mode, it's used for plotting; in other modes, it's the starting model for optimization.
- `-i [positive integer] (default: 100)` defines the maximum number of iterations. Useful for `recuit` mode.
- `-p [positive integer] (default: 10)` defines the population size. Useful for `genetic` mode.
- `-s` enables plotting of comparative curves across different selection modes. Requires `genetic` mode to work.

## Tests

Code coverage is performed using `pytest`.  
The latest version of the coverage report is accessible via the _Coverage.lnk_ shortcut.  
To regenerate the report, run the following command:

```bash
pytest --cov=dna --cov-report html test_*.py
```


## Modes

- **traditional** : Calculates and displays the spatial trajectory of a DNA sequence based on the provided conformation model.

- **recuit**: Optimizes the input conformation model using a simulated annealing algorithm. The goal is to “close” the DNA sequence by minimizing the distance between its start and end points.

- **genetic** : Uses a genetic algorithm to improve the input conformation model for the same purpose: promoting the circularization of the DNA chain.
