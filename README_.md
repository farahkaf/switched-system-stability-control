# Switched-System Stability and Control for Power Systems

This repository contains the MATLAB implementations associated with the numerical examples and control methods presented in my Master 2 internship report on **switched-system stability and control for power-grid applications**.

The repository is intended to make the main numerical results reproducible without reproducing every intermediate or experimental script developed during the internship.

## Repository content

The codes are organized around two groups of examples:

- **Chapter 3 — Stability and control-design tools:** small two-mode examples used to illustrate common quadratic, polyquadratic and polyhedral methods.
- **Chapter 4 — Five-node power-grid example:** switching analysis, common-certificate groups, additional controller design and path-dependent Lyapunov verification.

A typical organization is:

```text
switched-system-stability-control/
│
├── README.md
├── data/
│   └── ss_matrices_A_conns.mat
│
├── examples/
│   ├── common_quadratic_example.m
│   ├── polyquadratic_example.m
│   └── polyhedral_control_example.m
│
├── power_grid/
│   ├── five_node_certificate_groups.m
│   └── path_dependent_bridge.m
│
└── utils/
    ├── bitsoris_test.m
    ├── isSubsetScaled.m
    ├── make_set.m
    ├── image_set.m
    └── plot_invariance.m
```

The exact folder names can be changed, but the same dependencies between the scripts should be preserved.

## Main scripts

### `common_quadratic_example.m`

Two-mode numerical example for the **common quadratic Lyapunov approach**.

The script:

- designs mode-dependent feedback gains using LMIs;
- searches for a common quadratic Lyapunov function;
- verifies the Lyapunov decrease inequalities directly;
- plots the common quadratic level set and its images;
- constructs a common invariant polyhedral set;
- performs a Bitsoris-type invariance/contraction verification.

This script corresponds to the common quadratic and invariant-set illustrations used in Chapter 3.

### `polyquadratic_example.m`

Two-mode numerical example for **polyquadratic stability**.

The script:

- defines two stable closed-loop modes;
- solves the polyquadratic LMIs;
- associates a different Lyapunov matrix with each mode;
- checks the transition inequalities;
- constructs and plots a common invariant polyhedral set used for numerical illustration.

This script is used to illustrate why a mode-dependent Lyapunov description can be less conservative than a single common quadratic certificate.

### `polyhedral_control_example.m`

Implementation of the constructive **polyhedral control design** used in the report.

The script:

- starts from an initial polyhedral set;
- computes mode-dependent predecessor sets;
- projects them onto the state space;
- intersects the resulting sets;
- iterates the construction until the stopping condition is satisfied;
- computes the controller matrices;
- simulates the resulting closed-loop switched system.

### `five_node_certificate_groups.m`

Main preliminary analysis for the **five-node, nine-mode power-grid model**.

The script:

- loads the system matrices;
- builds the initial mode-dependent controllers;
- constructs the closed-loop matrices;
- searches for subsets of modes sharing common quadratic certificates;
- identifies the main certificate groups used in the report;
- designs the additional controller associated with topology 2;
- compares direct and modified switching sequences.

The main groups used in the numerical study are

\[
\mathcal{I}_1 = \{1,4,7\},
\qquad
\mathcal{I}_2 = \{2,3,5,6,8,9\}.
\]

### `path_dependent_bridge.m`

Final path-dependent implementation used for the Chapter 4 numerical study.

The considered periodic sequence is

\[
1 \rightarrow 4 \rightarrow 7 \rightarrow \widetilde{2}
\rightarrow 2 \rightarrow 3 \rightarrow 2
\rightarrow \widetilde{2} \rightarrow 1.
\]

Here, \(\widetilde{2}\) represents the **same physical topology as mode 2**, but with a different feedback controller.

The script:

- computes the controller associated with \(\widetilde{2}\);
- forms the prescribed periodic switching path;
- checks the stability of the complete cycle;
- solves path-dependent Lyapunov LMIs along the selected transitions;
- recovers the Lyapunov matrices;
- verifies the original decrease inequalities directly;
- simulates the repeated switching path;
- evaluates the active path-dependent Lyapunov function.

## Utility functions

### `bitsoris_test.m`

Checks a polyhedral invariance or contraction condition of the form

\[
FH = H(A+BK),
\qquad
Fw \leq \lambda w,
\qquad
F \geq 0
\]

with elementwise non-negativity imposed on \(F\).

### `isSubsetScaled.m`

Tests whether one polyhedral set is contained in a scaled version of another set.

### `make_set.m`

Creates a polyhedron from its half-space representation and computes a minimal representation.

### `image_set.m`

Computes the image of a polyhedral set under a linear map.

### `plot_invariance.m`

Plots a polyhedral set together with its image in order to visualize invariance.

## Requirements

The codes were developed in MATLAB and require:

- MATLAB
- [YALMIP](https://yalmip.github.io/)
- [MOSEK](https://www.mosek.com/)
- [MPT3](https://www.mpt3.org/)
- MATLAB Control System Toolbox
- MATLAB Optimization Toolbox

A valid MOSEK installation and license are required for the scripts that solve semidefinite programs.

## How to use the repository

### 1. Download the repository

Clone it with Git:

```bash
git clone https://github.com/farahkaf/switched-system-stability-control.git
```

or use **Code → Download ZIP** from the GitHub page.

### 2. Open MATLAB in the repository folder

Set the repository as the current MATLAB folder and add all subfolders to the MATLAB path:

```matlab
addpath(genpath(pwd))
```

### 3. Check the required toolboxes

Make sure YALMIP, MOSEK and MPT3 are installed correctly before running the examples.

For example:

```matlab
yalmiptest
```

can be used to check the YALMIP installation.

### 4. Check the data file

The five-node examples require:

```text
ss_matrices_A_conns.mat
```

This file contains the plant matrices `A_matrices` and the input matrix `B`.

If the data file is stored in the `data/` folder, use a relative path such as:

```matlab
dataFile = fullfile('data','ss_matrices_A_conns.mat');
```

Avoid absolute paths such as `C:\Users\...\Downloads\...`, since these will not work on another computer.

### 5. Run the examples

A suggested order is:

```text
1. common_quadratic_example.m
2. polyquadratic_example.m
3. polyhedral_control_example.m
4. five_node_certificate_groups.m
5. path_dependent_bridge.m
```

The first three scripts reproduce the small numerical examples used to illustrate the theoretical tools. The last two scripts correspond to the five-node power-grid study.

## Relation with the internship report

The repository is intended as a companion to the internship report:

- the **two-mode scripts** reproduce the numerical illustrations of the stability and control-design tools;
- the **five-node scripts** reproduce the switching analysis and controller-design results for the power-grid benchmark;
- the **utility functions** contain the numerical routines used by the main scripts.

The repository contains the implementations needed to reproduce the main examples rather than every intermediate code written during the development of the project.

## Notes

Numerical results can depend slightly on solver versions, tolerances and installed MATLAB/YALMIP versions.

For the semidefinite programs, solver feasibility is complemented by direct numerical checks of the recovered Lyapunov inequalities whenever possible.

## Author

**Farah Kafnemer**

Master 2 internship — switched-system stability and control for power-grid applications.

## Repository

https://github.com/farahkaf/switched-system-stability-control
