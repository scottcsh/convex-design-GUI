<p align="center">
<strong>Convex Design GUI</strong>
</p>

# Convex Design GUI

A lightweight PySide6 desktop GUI for scaffold-guided convex/concave protein binder design using RFdiffusion and the 5HCS scaffold library.

The GUI is designed for Rocky Linux workstation/server environments and focuses on generating and running RFdiffusion jobs from a single target structure, then ranking output binders by paper-style sphere convexity.

---

## Table of Contents
- [Installation](#installation)
- [System Requirements](#system-requirements)
- [First-run Configuration](#first-run-configuration)
- [Running the App](#running-the-app)
- [Input and Output](#input-and-output)
- [Design Workflow](#design-workflow)
- [Convexity Calculation](#convexity-calculation)
- [Acknowledgements](#acknowledgements)

---

## Installation

Create a virtual environment and install Python dependencies:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install --upgrade pip
pip install -r requirements.txt
```

If your system does not already provide virtualenv support:

```bash
sudo dnf install -y python3 python3-pip python3-virtualenv
```

---

## System Requirements

Required external components:

- RFdiffusion repository and runtime environment
- 5HCS scaffold library
- A pre-extracted 5HCS structure directory
- CUDA-capable GPU for RFdiffusion inference

Python packages are listed in:

```text
requirements.txt
```

Current Python dependencies:

- PySide6
- Biopython

---

## First-run Configuration

Open the GUI and click **Settings**.

Configure the following paths:

- `HCS_DIR`
  - Root directory containing the downloaded 5HCS archive files.
- `HCS_EXTRACTED_DIR`
  - Directory containing extracted 5HCS PDB/CIF structures.

Settings are stored locally in:

```text
resources/settings.json
```

This file is intentionally ignored by Git because it contains machine-specific paths.

---

## Running the App

From this directory:

```bash
source .venv/bin/activate
python main.py
```

If Qt reports DBus warnings in SSH/X11 or non-GNOME sessions, the app may still work normally. If the GUI does not open, try:

```bash
dbus-run-session -- python main.py
```

---

## Input and Output

Input structure:

- `.pdb`
- `.cif`
- `.mmcif`

The GUI requires:

- Input PDB/CIF file
- Target chain
- Library search length
- Hotspot residues
- Number of designs
- Top N designs to select
- GPU ID

If the input structure contains multiple chains, the selected target chain is extracted to a target-only PDB before RFdiffusion is run. The original input structure is not modified.

Runtime files are written under:

```text
<input_directory>/_convex_design_gui/
```

---

## Design Workflow

For each run:

1. Validate the input target structure and selected target chain.
2. Calculate paper-style sphere convexity for the selected target-chain hotspot patch.
3. Extract the selected target chain into a target-only PDB.
4. Generate an RFdiffusion shell script.
5. Filter extracted 5HCS scaffolds by residue length.
6. Generate scaffold-guided secondary-structure and adjacency files.
7. Run RFdiffusion scaffold-guided inference.
8. Monitor generated PDB files while the job is running.
9. Calculate output binder convexity, convexity difference, shape match, and binder length.
10. Copy the top N designs and write `selected_top_designs.csv`.

---

## Convexity Calculation

The GUI implements a paper-style sphere convexity approximation:

1. Collect interface heavy atoms.
2. Fit a sphere to interface heavy atoms using RANSAC.
3. Calculate convexity as the reciprocal of the fitted sphere radius.
4. Assign sign by comparing the fitted sphere center direction to the protein center direction.

Current sign convention:

- Negative value: convex surface
- Positive value: concave surface

---

## Acknowledgements

This GUI is a thin interface layer around external research software and datasets including:

- RFdiffusion
- 5HCS scaffold library
- Biopython
- PySide6

Please follow the licenses, citation policies, and usage requirements of each upstream project.

</br>
</br>

[Return to top](#table-of-contents)
