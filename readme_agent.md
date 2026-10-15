# Agent runbook: diffmol.py, quantmol.py, solvmol.py

## 0. Audience, authority and scope

This document is an operational contract for an AI agent, not a promotional
feature list. Use it to choose tools, build noninteractive commands, validate
outputs and restrict scientific claims. It covers ONLY these three scripts.

The inspected implementations are diffmol.py (the supplied diffmol-3.py),
quantmol.py 0.3.0 with single-conformer mode, and solvmol.py (the supplied
solvmol-4.py). Upload suffixes are not required runtime names. Resolve the actual
files in the workspace; do not rename or replace scripts merely to match examples.

The installed source and `--help` are the final authority for execution. If they
contradict this runbook, stop relying on that part of the runbook, inspect the
relevant function, and record the version mismatch. Do not invent CLI switches,
auto-detect undocumented modes or change source merely to make a command succeed.

No host-specific directory is part of this contract. Resolve script paths,
interpreter, xTB executable, density library and output locations from the current
project configuration or user instructions. Pass resource paths explicitly.
Do not rely on a developer's shebang, executable default, or density-directory
default being valid on a different machine.

## 1. Task routing

| User question | Tool/mode | What it can establish |
|---|---|---|
| Where/how do two geometries differ? | diffmol, paired inputs | Connectivity heuristics, distance patterns, IC changes, coordination/rings/shape |
| How do corresponding fragments differ geometrically? | diffmol, paired selections | Selected-fragment geometry, with full-REF branch/topology context |
| Is a single structure's AO basis near-redundant? | quantmol, one input | Numerical overlap health in the specified basis |
| Which contacts and model dispersion contributions exist in one structure? | quantmol, one input | Absolute D4 contribution tables and reference-density pair-IGM contacts |
| Where do model contact/dispersion contributions change? | quantmol, paired inputs | D4 and IGM changes, overlap signatures and AO-metric distortion |
| Does solvent alter the preference between whole conformers? | solvmol, paired complete molecules | Fixed-geometry GFN-xTB/ALPB gas/solution comparisons and SASA/dipole descriptors |
| Need total conformational Gibbs energies, MR diagnosis, spin ranking, reaction barriers? | None of these alone | Explain the missing evidence; do not manufacture it |

Do not automatically run all three. For a broad two-conformer study, use diffmol
first, then quantmol for model contributions, then solvmol only when solvent and
electronic-state inputs are established. For one structure, use quantmol directly;
do not fabricate a comparison by repeating the same XYZ.

### Computational contracts

- diffmol is coordinate/topology analysis; no molecular electronic calculation.
- quantmol analyzes frozen input geometries without molecular SCF. Its explicit
  reference-density preparation DOES run separate atomic HF calculations.
- solvmol DOES launch molecular GFN1/2-xTB SCC single points: REF gas, REF solution,
  PROBE gas, PROBE solution. Do not describe it as no-SCF or merely geometric.
- None optimizes these input geometries or supplies vibrational thermochemistry.

## 2. Mandatory preflight

Before calculations:

1. Resolve an explicit Python interpreter and actual script files. Prefer
   `"$PYTHON" "$SCRIPT" ...` to assuming executable permissions/shebangs.
2. Read the relevant `--help`. quantmol also supports `--version`; the other two
   do not implement a general version option. Record source hashes when needed.
3. Require one complete finite-coordinate XYZ frame per input, coordinates in
   Angstrom, correct explicit hydrogens, and no coincident atoms. Prevalidate frames:
   quantmol's reader does not reliably reject trailing trajectory frames, whereas
   diffmol/solvmol reject extra nonblank content. Convert trajectories upstream.
4. Preserve original input files. Do not reorder, add hydrogens, cap, optimize or
   change chemical identity without explicit user authorization. Record any
   authorized transformation and its atom mapping.
5. Establish the task: single vs comparison, full molecule vs geometric fragment,
   intended atom correspondence, required modules and resource budget.
6. For energy comparisons, establish chemical identity, composition, protonation,
   charge and intended electronic state. Same element counts are insufficient.
7. Confirm required dependencies and data coverage before launching a large job.
8. Create a NEW, task-specific working directory. Resolve input/script paths to
   absolute paths before changing directory. Use unique prefixes. This matters
   particularly because diffmol always writes a fixed-name response CSV by default.
9. Specify inputs that would otherwise prompt. Do not let unattended jobs hang.
10. Save exact command, interpreter/backend versions, stdout, stderr, exit status,
    mapping, parameters and reasons for intentionally skipped modules.

If a scientific input such as charge or xTB UHF is unknown, ask the user; do not
infer neutrality/closed shell merely from a successful parser default.

### Noninteractive interface differences

| Script | No-prompt full geometry selection | Fragment selection | Other potential prompts |
|---|---|---|---|
| diffmol | `--all` | BOTH `--ref-atoms` and `--probe-atoms` | Selections if omitted |
| quantmol, one input | `--all` | `--ref-atoms` ONLY | One selection if omitted |
| quantmol, two inputs | `--all` | BOTH paired selection flags | Two selections if omitted |
| solvmol | Entire input molecules automatically; NO `--all` | Not supported | Charge, UHF, missing probe radius |

C/c means all at applicable selection prompts. Indices/ranges are 1-based,
inclusive, ascending within a range, unique, and preserve supplied order.
Paired selection order defines correspondence; there is no atom-mapping search.
`--all` means identical file-order mapping for diffmol/quantmol, NOT an automatic
matching operation. Do not combine it with selection flags.

Solvmol checks elemental composition but does not require identical atom ordering;
therefore do not compare per-atom SASA records across reordered files without a
separate mapping. It does not prove the inputs are conformers rather than isomers.

### Dependencies

- diffmol: NumPy and QCElemental; RDKit for CORE/shape; NetworkX for whole-molecule
  topology. For a reduced geometric run, BOTH `--skip-advanced --skip-whole` avoid
  those advanced modules, not NumPy/QCElemental.
- quantmol: NumPy/SciPy; PySCF for S/T and explicit atomic-density preparation;
  dftd4 Python API for dispersion; radial density JSON library for IGM; complete
  compatible Thole parameters for polarization. No RDKit is required by this script.
- solvmol: NumPy, RDKit, FreeSASA Python module (`import freesasa`), and an xTB
  executable supporting the requested GFN/ALPB/JSON options. Source documentation
  specifies xTB >=6.3.3; verify actual solvent/features in the installed version.

Do not install packages globally, edit environments or download reference libraries
silently. A successful import is not a scientific parameter-coverage check.

## 3. Output isolation and reusable command templates

The following are templates. Bind variables from verified project configuration;
placeholder names are not literal paths and have no automatic meaning to scripts.
`CHARGE` and `UHF` below must be supplied scientifically, not guessed.

Create a new working directory before redirection so log files do not collide.
Capture exit status before launching any subsequent command. If using `tee`, use
`set -o pipefail` or inspect the Python process status; tee's success is not enough.

### Geometry comparison

```bash
"$PYTHON" "$DIFFMOL" "$REF_XYZ" "$PROBE_XYZ" \
  --all --csv-prefix geometry > geometry.stdout.log 2> geometry.stderr.log
```

### Single-structure quantum/model diagnostics

For a metal system without validated Thole parameters, explicitly skip polar:

```bash
"$PYTHON" "$QUANTMOL" "$REF_XYZ" --all \
  --basis def2-svp --charge "$CHARGE" \
  --density-library "$DENSITY_LIBRARY" --skip-polar \
  --strict --csv-prefix single > single.stdout.log 2> single.stderr.log
```

### Two-conformer model diagnostics

```bash
"$PYTHON" "$QUANTMOL" "$REF_XYZ" "$PROBE_XYZ" --all \
  --basis def2-svp --charge "$CHARGE" \
  --density-library "$DENSITY_LIBRARY" --skip-polar \
  --strict --csv-prefix model > model.stdout.log 2> model.stderr.log
```

If IGM is not required or no validated density library exists, explicitly use
`--skip-igm`, explaining the omission. Do not present the run as an IGM analysis.
If polar response is required, skipping it is not a substitute for obtaining valid
parameters. `--strict` tests module availability/coverage, not every numerical
contact's convergence.

### Solvent comparison

```bash
"$PYTHON" "$SOLVMOL" "$REF_XYZ" "$PROBE_XYZ" \
  --solvent water --charge "$CHARGE" --uhf "$UHF" \
  --xtb "$XTB_EXECUTABLE" --threads 1 \
  --output-dir xtb_records > solvent.stdout.log 2> solvent.stderr.log
```

`xtb_records` must not already exist. If the solvent lacks a radius preset, add
`--probe-radius "$PROBE_RADIUS_A"` using a documented/user-approved value.
Water in this example has a built-in preset; it is not a default solvent choice.

## 4. diffmol.py execution and interpretation contract

### 4.1 Scope and selection semantics

Exactly two positional XYZ files are required. Fragment comparisons are geometric
subsets; the script uses full-molecule context for some interaction/drive-response
operations. This is NOT the same as quantmol's isolated selected-subsystem model.

Full-molecule report is triggered by `--all` or by entering C/c for BOTH selections
(including supplying C/c in both flags). An explicit numeric list covering all atoms
does not necessarily trigger WHOLE_MOLECULE. Use --all if the global report is wanted.

Example of an intentional paired geometric fragment mapping:

```bash
"$PYTHON" "$DIFFMOL" "$REF_XYZ" "$PROBE_XYZ" \
  --ref-atoms 1-6,12 --probe-atoms 21-26,30 --csv-prefix fragment
```

Validate element identity and chemistry at each paired position beforehand.

### 4.2 Sections and admissible conclusions

- CONNECTION: compares distance-inferred covalent connectivity and candidate
  hydrogen/halogen/noncovalent contacts. Labels are geometric heuristics, not
  electronic bond orders, reaction events or interaction energies.
- DISPLACEMENT: DDM = D(PROBE)-D(REF), no alignment needed. dRMSD describes pair
  distance changes, not Cartesian atom displacement. DDM spectral fractions are
  squared-matrix-change fractions, not PCA statistical variance, energies or
  fractions of independent physical motion. Mirror inversions can preserve all
  pair distances. Do not export/describe DDM patterns as normal modes or pathways.
- CORE: inferred stereocenters, metal coordination shape and ring puckering.
  RDKit stereochemistry depends on inferred chemistry; metal shape changes do not
  determine oxidation state, spin, coordination strength or reaction mechanism.
- INTERNAL_COORD: selected covalently connected REF-topology bond/angle/dihedral
  changes. REF topology defines enumerated coordinates; connectivity warnings
  reduce interpretation confidence. Periodic angle differences need their own
  wrapped conventions; consume exported delta rather than subtracting naively.
- IC DRIVE-RESPONSE: rigid-branch constructions drive one coordinate toward the
  PROBE value and observe others. This is not an energy-relaxed response, statistical
  correlation, force constant or causality. IC_impact is not an additive percentage.
- WHOLE_MOLECULE: graph/shape/mass-spread descriptors and mapped proper alignment.
  Heuristic shape bands are not calibrated chemical significance thresholds.
  Near-degenerate principal axes have ambiguous orientation; inspect stability flags.

Response ratio = induced response / observed REF-to-PROBE change of the RESPONSE
coordinate. MATCH is 0.8-1.2, PARTIAL positive <0.8, OVERSHOOT >1.2, OPPOSITE negative,
after visibility filters. Default candidate selection is up to 10 changes each
for bonds, angles and dihedrals, at most 870 ordered nonself pairs; only 20 prominent
relations are printed by default. B1/A1/D1-style IDs must be mapped through candidate
records; they are not atom numbers. WEAK/N/A is not absence of physical coupling.

IMPORTANT: the current console suggests not constraining MATCH pairs simultaneously.
Treat that as a heuristic redundancy warning, not an instruction to edit a real
optimization. Do not choose/freeze/release constraints based only on this branch
model; verify the actual optimization strategy and obtain authorization.

### 4.3 Output and overwrite behavior

With --csv-prefix P, expected conditional files include:

- P_all.csv / P_top.csv / P_excluded.csv: internal-coordinate records;
- P_covalent_changes.csv / P_contacts.csv;
- P_chirality.csv / P_metal_coordination.csv;
- P_ring_puckering.csv / P_ring_modes.csv;
- P_ddm_spectrum.csv / P_ddm_atoms.csv / P_ddm_pairs.csv / P_ddm_loadings.csv;
- P_ddm_matrix.npz only with --ddm-save-matrix (requires a prefix);
- P_ic_response_candidates.csv / P_ic_response_all.csv / P_ic_response_top.csv;
- P_whole_shape.csv / P_whole_axes.csv / P_whole_summary.csv when applicable.

CSV writes overwrite matching paths without a universal overwrite guard. Output
parents are not comprehensively created by diffmol: create directories beforehand.
There is no native results JSON or --output-dir flag.

CRITICAL SIDE EFFECT: when IC response is enabled, `ic_response_prominent.csv` is
ALWAYS written in the working directory, replacing any existing file, even without
--csv-prefix. It contains ALL prominent responses, not just displayed rows. A unique
working directory is mandatory for parallel/batch jobs. --skip-ic-response disables
that analysis/export but does not necessarily eliminate IC_impact computation.

### 4.4 Resource/failure handling

- --ddm-max-atoms default 500 can skip DDM while the script otherwise succeeds.
- Ring enumeration defaults to all simple cycles with a 10,000-ring safeguard;
  inspect truncation. --ring-mode symmsssr changes the ring definition and should
  be recorded, not used as if identical. --ring-max-size restricts scope.
- --skip-advanced does not disable WHOLE_MOLECULE dependencies; use --skip-whole too
  when intentionally doing reduced analysis.
- --top is a reporting limit, not a universal compute limit.
- Handled input/dependency errors exit 2. A zero exit code does not guarantee every
  optional section ran: parse skipped/truncated/undefined and warning messages.

## 5. quantmol.py execution and interpretation contract

### 5.1 Modes, mapping and coordinates

One positional XYZ -> single mode. Two -> comparison mode. No --single flag exists.
Single mode accepts --ref-atoms, rejects --probe-atoms and --align-atoms, and keeps
the input coordinate frame. Its metric module is intentionally NOT_APPLICABLE.

Comparison mode uses file-order/explicit paired mapping and aligns PROBE to REF
BEFORE integral construction. --align-atoms specifies original selected REF indices
for alignment; default uses all selected atoms. Use a common rigid scaffold when
appropriate. Collinear/single-atom alignment makes axis-dependent results ambiguous.
Changing alignment rules changes interpretation of metric/tensor comparisons.

Selected atoms form the COMPLETE modeled subsystem. There are no caps and omitted
atoms contribute no field/density/dispersion. Supply the selected subsystem's charge.
For chemistry use --all plus named reporting fragments rather than casually cutting
an organic ligand from a metal complex. PySCF neutral/parity Mole initialization
is not an electronic-state solution and does not infer oxidation state or spin.

### 5.2 Fragment and topology input schemas

`--fragments FILE` is a JSON object of disjoint, nonempty named lists of ORIGINAL
REF atom indices. All indices must be selected. Unlisted atoms form unassigned.
The name unassigned is reserved if a remainder needs that name. Example structure:

```json
{"metal": [1], "ligand_A": [2,3,4,5], "substrate": [6,7,8]}
```

Example only: choose actual chemically meaningful atom lists, never assume these
numbers fit the user's molecule. In single mode the same original-input convention
applies. Without fragments each selected atom is its own reporting fragment.

`--igm-bonds FILE` is a COMPLETE replacement nonmetal covalent graph (not an
additive patch), e.g. `[[2,3],[3,4]]`, using original selected REF indices. Metal
edges are deliberately omitted from covalent exclusion paths. Default inference
uses the union of REF/PROBE distance-based nonmetal graphs. Inspect changed-edge
warnings; unsupported radius coverage needs an explicit graph. These graph rules
are not an electronic bond-order analysis.

Singleton IGM excludes 1-2/1-3 by default; 1-4 remains. Metal singleton pairs are
excluded unless --igm-include-coordination. Metal-centered paths never suppress
ligand-ligand pairs. Named fragments crossing bonds retain a COVALENT_BOUNDARY label;
nonlinear fragment IGM is NOT decomposed and subtracted to remove such contacts.
Named fragment metal pairs are COORDINATION_OR_METAL, not automatically actual bonds.

### 5.3 Numerical modules

BASIS_HEALTH: normalized AO overlap eigenvalues/condition, low-subspace shell
weights and retained-space kinetic spectra. Default thresholds 1e-6,1e-8,1e-10.
Small modes indicate numerical near-redundancy, not high correlation. Shell weights
locate participation, not unique causes. Tight functions legitimately have large
kinetic energies. Never delete basis functions automatically. Check sensitivity
and intended calculation settings before proposing a change.

FINGERPRINT: uses --basis or optional --fingerprint-basis, not a true density.
Comparison global/local spectral distances localize overlap-environment changes.
Local neighborhoods use the union in either geometry under a hard radius cutoff;
identical neighbor sets can yield duplicate local entries. Do not count duplicates
as independent chemical evidence. Default centers: all selected atoms; local radius
4 A. Single mode saves a global spectrum only, no local difference ranking.
Across molecules/bases/elements, magnitudes are not directly calibrated.

METRIC: generalized eigenvalues SB v = gamma SA v; log distance = norm(log gamma),
log RMS divides by sqrt(dimension). Gamma extrema bound squared-norm changes under
unchanged AO coefficients. They do NOT measure actual occupied-state change,
charge transfer or guaranteed DIIS failure. Below --metric-threshold 1e-8, full-space
results are withheld. --metric-project explicitly uses a common retained space,
rechecks positivity and labels its dimension. Do not silently jitter or separately
prune the two spaces to force a number.

Single mode builds S/T once per basis; a different fingerprint basis adds another
build. --max-ao is a safeguard, not an accuracy knob or speedup. Raising it to 5000
permits larger dense matrices with substantial CPU/memory cost. Do not promise
universal subsecond analysis.

### 5.4 D4 module

--d4-method (default pbe0) selects damping parameters, NOT a molecular DFT calculation.
Input charge affects the model; verify it explicitly. Coordinates are converted
internally to Bohr. Pair-matrix sums are checked against backend total energy.

Comparison signs: delta = PROBE - REF. Negative total delta favors PROBE IN THE D4
COMPONENT ONLY. Pair/fragment changes refer to the additive two-body contribution;
the three-body difference is reported separately. Do not interpret its pairwise
representation as a unique physical pair decomposition.

Single: total/two-body/three-body, absolute atom pairs and fragment sums. Pair
categories are COVALENT_12, COVALENT_13, NONBONDED, METAL_INVOLVING. All remain in the
total; only reporting is partitioned. The single-mode category graph uses
--igm-bonds/inference even when IGM itself is skipped. Unsupported topology can
therefore affect single D4 categorization.

Do not call a total D4 correction an unfolding energy or binding energy. No separated
fragment/reference calculation is defined. Do not add D4 energy onto a solvmol total
or an externally computed dispersion-inclusive energy without an explicitly justified
method; these outputs are not additive modules of a universal total energy.

The current dual console may state a preference for a tiny nonzero difference.
Report its numerical value but call near-cancelling values practically similar
within this model rather than treating many printed decimals as accuracy. A chosen
reporting band must be disclosed and is not a model uncertainty estimate.

### 5.5 Reference density library

Always pass `--density-library "$DENSITY_LIBRARY"` in reproducible commands; do
not encode a specific computer's directory. The library is reusable across runs.
Each required element needs Element.json with element, units='atomic', source,
r_bohr, rho, drho. Radii are Bohr, density and radial derivative atomic units;
values must be finite, radial grid strictly increasing from zero, density nonnegative,
and tail negligible. The current reader validates structure/tail but does NOT fully
certify imported normalization or physical provenance. Validate external data separately.

Do not invent Gaussian densities to bypass missing physical data. Synthetic tables
are appropriate only for tests. Reference densities are frozen neutral atoms, not
actual molecular density/spin density; no chemical polarization/charge transfer.

Explicit preparation (not routine target-molecule analysis):

```bash
"$PYTHON" "$QUANTMOL" --prepare-density "$NEW_DENSITY_LIBRARY" \
  --density-elements H,C,N,O --density-basis def2-tzvp
```

This invokes atomic HF, convergence checks, total alpha+beta density and radial
normalization. It can take time and fail for individual atoms. It never occurs
implicitly in an IGM run. Existing element files require --overwrite; do not
overwrite a shared library without authorization. Prefer isolated per-element
attempts/logs and retain successful references.

The builder assumes ALL-ELECTRON atoms. Standard def2 heavy-element bases from Rb
are ECP-oriented: do not use them as full-electron references or cure occupancy
index errors by truncating occupations. It does not implement an ECP-core density
reconstruction. Select and validate a suitable reference method separately. Do not
mix physical reference conventions without recording per-element provenance.

### 5.6 IGM coverage, ranking and numerical status

Pair descriptor: |grad rhoF| + |grad rhoG| - |grad(rhoF+rhoG)|. Its integral has unit
bohr^-1, NOT energy. Threshold volume is A^3; threshold value is in atomic gradient
units. Pair fields overlap and are nonadditive. No NCI density-Hessian classification,
connected-component tracking, surface area or 3D rendering is implemented.

Defaults: 4 A candidate cutoff in either geometry, .25 A coarse step, 2 A padding,
thresholds 1e-4/1e-3/1e-2, all candidates coarse-evaluated, fast top-12 refinement.
'Complete coverage' means candidate coarse coverage after filters, not all physical
interactions or a global fine-grid proof. --igm-max-pairs 0 means no cap; a positive
cap is an incomplete deterministic lexical subset, NOT an importance ranking.

Two modes:
- comparison fast ranks coarse |delta J|; checks 12 largest;
- single fast ranks coarse absolute J; checks 12 largest.
--igm-mode complete or --igm-refine checks all evaluated candidates. Costs can grow
substantially. --igm-max-points limits each grid, default 1,000,000.

Checks use nested half-step cells in the same physical box, then a larger fine-grid
box. This separates spacing and domain effects. It does not refine the source radial
density tables or prove quadrature convergence with only one refinement step.

Comparison labels:
- RESOLVED_INCREASE / RESOLVED_DECREASE: sign-consistent delta passes baseline,
  boundary and delta-sensitivity tests under configured heuristics.
- NUMERICALLY_UNRESOLVED: no reliable change sign/magnitude established, NOT no contact.
- NOT_CHECKED: coarse only.
- RESOURCE_LIMIT: retain coarse result; no validation claim.

Single labels:
- CHECKED_CONTACT: absolute J passes integral and box sensitivity tests;
- GRID_SENSITIVE / BOX_SENSITIVE;
- NOT_CHECKED / RESOURCE_LIMIT.
No increase/decrease can be inferred in single mode. Self-comparison cancellation
cannot validate an absolute integral.

Default baseline tolerances: rtol=.02, atol=1e-5 bohr^-1. Comparison additionally
uses delta_atol=1e-4 bohr^-1 and safety=3 times summed grid/padding delta sensitivities;
extra padding 1 A. These are numerical heuristics, not statistical confidence or
rigorous error bounds. Single mode does not use the delta floor/safety settings.
Threshold volumes have sensitivity data but no automatic formed/lost-contact test.

Do not sum pair J as a global interaction energy. A larger J does not uniquely mean
stronger attraction, larger Pauli repulsion or a stronger hydrogen bond. To name a
contact mechanism, inspect element identities, topology, geometry and independent
energetic/electronic evidence. No force or barrier is computed.

### 5.7 Polarization parameter contract

Current model is the empirical LINEAR damped Thole set, a=1.7278, with effective
isotropic site values for H,C,N,O,F,S,Cl,Br,I. Fe,Zn,P and other absent species do
NOT receive fallback parameters. A supported atomic symbol does not guarantee
accuracy in every molecular environment.

--polar-params FILE REPLACES the default parameter source; it is not an incremental
one-element patch. Schema:

```json
{
  "model": "thole_linear",
  "units": "bohr^3",
  "source": "Documented compatible parameterization",
  "a": 1.7278,
  "elements": {"H": 3.5020, "O": 6.3940},
  "atoms": {}
}
```

This is a schema example for covered H/O, NOT a complete metal parameterization.
Cover every selected atom. Optional atoms keys are original REF atom numbers and
override element values. Require positive finite values and compatible damping.
Do not transplant an AMOEBA/Drude value into this different formula without validation.
Do not use free-neutral-metal polarizability, zero, or a tiny epsilon to bypass coverage.
Do not delete the metal silently. To run without response, explicitly --skip-polar
and disclose it; if polarization is the requested result, parameters are required.

Output is a classical frequency-independent response tensor, mean, anisotropy and
optional axis projection. It is not a spin-resolved quantum polarizability or Raman
intensity (which requires derivatives). Response positivity/conditioning is checked;
near-instability is not a real 'polarization hotspot'. Fragment sums are coupled
contributions under a common field, not isolated-fragment polarizabilities.

--response-axis i:j uses original selected REF indices. Dual mode uses the SAME REF
axis for both aligned tensors, not two independently rotating bond axes. Single
mode uses the input geometry. Principal directions need degeneracy caution.

### 5.8 Exports, statuses and error policy

With --csv-prefix P:
- P_quantmol.json: settings, selections/mapping, module statuses, conclusions and data;
- P_health_shells.csv / P_kinetic.csv;
- P_om_spectrum.csv; dual also P_om_local.csv;
- P_metric_spectrum.csv only for applicable valid dual metric;
- P_d4_pairs.csv / P_d4_fragments.csv; single also P_d4_categories.csv;
- P_polar_tensor.csv / P_polar_fragments.csv;
- P_igm_pairs.csv / P_igm_volumes.csv / P_igm_filtered.csv / P_igm_excluded.csv.
Empty datasets need not produce a CSV. Missing file is not automatically zero.
No prefix -> no module CSV/JSON export. Parents are created at export. Existing
matching outputs require --overwrite; check is late, so unique prefixes avoid wasted
computation. No NPZ/cube/surface output is implemented by this quantmol version.

Single JSON analysis_mode='single', mapping.atoms, health.value, fingerprint.spectrum,
dispersion.total_kcal, polarization.tensor_A3, igm.pairs[].integral. No fake deltas.
Dual mapping.ref/probe, health.ref/probe, fingerprint.spectrum_ref/probe,
dispersion.ref_kcal/probe_kcal/delta_kcal, igm.pairs[].delta_integral.
Single om_spectrum CSV is mode/eigenvalue; dual is mode/ref/probe/delta. Do not parse
both modes with the same assumed schema. Internal topology indices are selected-array
zero-based; translate via stored mapping before reporting original atom numbers.

Module statuses include OK, PARTIAL, UNAVAILABLE, ILL_CONDITIONED, and intentional
single METRIC NOT_APPLICABLE. A non-OK requested module makes overall run PARTIAL;
NOT_APPLICABLE does not. Default PARTIAL can exit 0; --strict exits 3. Handled input
errors exit 2; unexpected errors may produce other nonzero failures/tracebacks.

Strict does NOT ensure every contact passed numerical checks: RESOURCE_LIMIT or
NOT_CHECKED pairs can coexist with module OK when coarse coverage is complete.
Read pair-level labels, coverage, exclusions and warnings even after exit 0.
Do not confuse intentionally skipped (absent) modules with unavailable modules.

## 6. solvmol.py execution and interpretation contract

### 6.1 Inputs and calculation model

Exactly two COMPLETE same-species molecules, one solvent, one COMMON integer charge
and COMMON nonnegative unpaired-electron count. No selections, --all, per-conformer
states, or single-input mode. To compare different spin states, this CLI does not
express that task; do not misuse it.

--uhf is the xTB unpaired-electron count, not charge or multiplicity. Do not pass a
multiplicity value as --uhf. Source checks electron parity and maximum count, but
passing parity does not establish the correct physical state.

Coordinates are centered to COM, not optimized. Charged-molecule dipoles are origin
dependent; the COM convention is part of the result. Compare model total dipole
magnitudes with this caveat, not as invariant observables for ions.

Four isolated xTB single points are launched, each with --sp --json, matching
--gfn/charge/UHF/accuracy and gas or --alpb solvent. --state controls the ALPB state
convention; keep identical across conformers and compatible across any solvent study.
No explicit solvent molecules, coordination bridge, ZPE, entropy or ensemble.

### 6.2 Solvent and probe resolution

Use `--list-solvents` and `--list-probe-radii` first. --solvent NAME and --NAME are
mutually exclusive alternatives. Aliases are normalized for --solvent, but only
implemented ALPB names are accepted; a SASA probe preset does not add ALPB support.
Installed xTB must also support the name.

Radius priority: explicit --probe-radius > built-in source-tagged preset > prompt.
If no preset and stdin is not a TTY, the script fails rather than using water's
radius. Always supply a documented radius explicitly for such a noninteractive run.
Do not choose an arbitrary number just to avoid prompting.

SASA uses FreeSASA Shrake-Rupley with RDKit elemental VDW radii, explicit atoms and
by default 960 points/atom. Probe presets are imported PBF conventions, not actual
ALPB cavity radii. Changing --probe-radius affects ONLY FreeSASA, not xTB ALPB.
Hold radius fixed for comparisons; surface values at different radii are different
geometric definitions. C/H vs other surface is an element split, not a rigorous
polar/nonpolar chemical partition. --threads controls xTB; FreeSASA code fixes one
thread independently.

### 6.3 Energy algebra: always preserve these signs

Let R=REF, P=PROBE, g=gas and s=solution. All returned relative energies are P-R:

- A_R = E_R,s - E_R,g; A_P = E_P,s - E_P,g.
- delta_g = E_P,g - E_R,g.
- delta_s = E_P,s - E_R,s.
- delta_A = A_P - A_R = delta_s - delta_g.

JSON relative keys:
ref_stabilization_kcal, probe_stabilization_kcal,
relative_gas_kcal, relative_solution_kcal, relative_stabilization_kcal.

Negative delta_A: solvent preferentially stabilizes PROBE relative to REF in this
fixed-geometry model. It does NOT imply PROBE has lower solution energy; inspect
delta_s. Positive delta_s means REF remains lower in solution. Check algebraic
closure numerically and quote both quantities when discussing solvent effects.

Parsed model terms: Gsolv, Gelec, Gsasa, Ghb, Gshift (stored in Eh, printed kcal/mol).
Missing decomposition terms stay unavailable, not zero. Gsolv is already included
in the xTB solution total: NEVER add it again. FreeSASA area is not an extra energy
correction. Printed Gsolv need not equal E_solution-E_gas, which also includes the
model electronic relaxation between separate single points.

Default reporting band is --energy-tolerance 1 kcal/mol: a heuristic communication
band, not validated uncertainty or a rigorous significance test. SASA's 2% exposure
label is likewise heuristic. Do not infer solvent stabilization from dipole or
SASA alone when calculated model energies are available.

### 6.4 Records and failure handling

--output-dir must be NEW; existing directories are never overwritten. Without it,
a new solvmol_* directory is created in the working directory. No --csv-prefix,
--overwrite, --strict, or resume switch exists.

Success directory contains:
- settings.json (including resolved radius/provenance/executable);
- xtb_version.txt;
- ref_gas/, ref_solution/, probe_gas/, probe_solution/;
- each subdirectory: input.xyz, command.json, xtb.out, xtb.err, backend files;
- results.json with settings, conformers and relative;
- summary.csv (one row of settings and relative energy values).

Failure may leave settings/subdirectories without results.json. Preserve them for
diagnosis. Directory existence is not success. Check all four normal terminations,
SCC convergence, JSON energy/dipole fields, finite results, required Gsolv and
relative-energy closure. Script handles these checks but agents should audit logs
for relevant warnings. Do not synthesize missing final results from partial runs.

An xTB subprocess timeout is converted into an error; --timeout is PER single
point, not the entire four-job workflow. --accuracy is a numerical setting, not
chemical accuracy. Record --gfn, solvent, --state, probe radius and electron state.
Changing any for recovery requires recording that the comparison method changed.
Do not compare a GFN1 recovered conformer against a GFN2 original conformer.

Handled errors exit 2. On failure inspect the specific xtb.out/xtb.err and command,
not just the final one-line message. Use a fresh directory for a rerun. Do not
silently optimize geometries, change charge/spin, or remove atoms to force SCC.

## 7. Cross-tool evidence and reporting rules

1. Keep REF/PROBE meanings fixed throughout. Record swapped inputs explicitly.
2. Use original atom identifiers and element/fragment identities, not only IDs.
   Consecutive atom numbers do not establish chemical adjacency.
3. Classify every claim: direct numerical result; parameterized model result;
   interpretation hypothesis; or not established. Tie it to a file/row/module.
4. A named mechanism (H bond, pi stack, ligand exchange) needs identity/topology and
   relevant geometry. IGM size alone does not name it.
5. Agreement between reduced contact J and less favorable pair D4 strengthens a
   local contact/dispersion interpretation, not proof of total attraction loss.
6. Large local changes with nearly zero total D4 indicate possible redistribution/
   cancellation; inspect positive AND negative contributions. Do not use pair counts
   as net energy or sum nonadditive IGM values as energy.
7. A larger metric can coexist with healthy S spectra. It describes a different
   mathematical question, not a contradiction or automatic electronic transition.
8. No suite member diagnoses MR, spin ordering, oxidation state or barrier reliably
   from these descriptors alone. Acknowledge the missing evidence.
9. Do not add quantmol D4 to solvmol energy or append a FreeSASA 'energy' to ALPB.
10. Do not claim precise stability from tiny model differences; preserve numbers and
    explain model scope. Absolute D4 does not define a binding/unfolding reference.
11. Numerical checks are not model validation. Self-comparison zero is a consistency
    check, not proof that physical predictions are accurate.
12. If a partial result answers a subquestion, report it with scope; do not present
    unavailable response/IGM as computed. A neutral input default is not a charge
    determination. Metal parameter gaps are never cured by inventing constants.

### Agent final-response schema

Use concise prose plus a small evidence table when useful:

- Task and scope: modes, full/subset, REF/PROBE, charge/state, solvent/model.
- Execution: commands/settings, completed/skipped/failed modules, coverage/checks.
- Findings: 3-5 strongest supported results with units, signs, original atoms/fragments.
- Interpretation: connect multiple results without turning descriptors into energies.
- Limitations: unrefined contacts, topology assumptions, parameter gaps, model scope.
- Artifacts: identify recorded outputs and their purpose.
- Next computation ONLY if it materially resolves a stated uncertainty; no gratuitous
  full scans, library rebuilds or parameter changes.

## 8. Recovery decision table

| Symptom | Inspect | Allowed response |
|---|---|---|
| Agent is waiting for stdin | Missing selection/charge/UHF/radius | Stop safely, supply verified inputs explicitly; no guessed defaults |
| Import/executable error | Active interpreter and configured executable | Diagnose environment; do not edit chemistry |
| quantmol IGM library missing | Explicit directory + all element files | Select existing validated library or request explicit preparation |
| Atomic occupancy IndexError | Element, basis core/ECP conventions | Do not truncate occupation; use validated all-electron setup |
| Fe/Zn/P Thole unavailable | Compatible parameter coverage | Report/explicitly skip if optional; no neutral-metal fallback |
| IGM RESOURCE_LIMIT | Grid count, box/step, pair size | Preserve coarse label; raise budget only deliberately; do not lower accuracy silently |
| IGM unrefined/unresolved | Baseline/delta/box checks and coverage | Report limited claim; targeted refinement if needed |
| Whole/CORE/DDM missing in diffmol | Selection trigger, dependencies, safeguards | Record skips; choose supported reduced scope only if appropriate |
| diffmol output overwritten | CWD fixed-name response CSV/prefix | Use a new workdir; do not claim recovery of lost files |
| quantmol output collision | Prefix and --overwrite policy | New prefix by default; overwriting requires authorization |
| solvmol existing output dir | Directory existence | New directory; no built-in overwrite/resume |
| xTB SCC/normal-termination/JSON failure | Individual job stdout/stderr/command | Preserve diagnostics; no completed scientific result |
| Solvent has no radius preset | Both ALPB list and SASA radius list | Ask/use documented approved radius; no implicit water fallback |
| Tiny self-comparison differences | Magnitude relative to precision/checks | Treat as numerical residuals, do not rank hotspots |

## 9. Capability checks and testing

Use installed --help to confirm CLI before execution. If supplied, run
`test_quantmol.py`, `test_quantmol_igm.py`, and `test_quantmol_single.py` in an
isolated directory before depending on a modified installation. Missing test files
are not failures of physical calculations. Report passed/skipped/backend-mocked
counts separately. A test name containing 'real change' can refer to a nonzero
synthetic Gaussian change, NOT a real molecule/backend test.

For optional installed-backend smoke tests, check whether PySCF and dftd4 actually
ran; absence-based skips do not validate them. Do not regenerate synthetic test
references as production element density data.

Useful task checks:
- comparison self-test: difference quantities near numerical zero;
- single-vs-self comparison: absolute D4/polar values agree, absolute IGM requires
  its own baseline checks;
- reordered/mapped structures: consistent original atom identifiers;
- coordinate rotation: scalar quantities stable within numerical/grid tolerance,
  tensors transform; a finite Cartesian grid is only approximately invariant;
- solvmol: delta_s = delta_g + delta_A, no Gsolv double-counting.

## 10. Source-checked CLI inventory

This inventory records implemented parser options, not proposed features. The
script help/source wins after upgrades. Site-specific resource defaults are
intentionally omitted: resolve them at runtime. Positional arguments and dynamically
generated flags are described separately. '--' option defaults below are parser
values, not recommended scientific assumptions.

