---
name: polygview
description: Use when quantum-chemistry small-molecule calculations require rapid comprehensive comparison of two whole-molecule geometries or corresponding fragments during geometry optimization or relaxed/flexible potential-energy scans, rapid single-conformer whole-molecule quantum/model diagnostics, or comparison of solvent effects on two whole-molecule conformers. Read readme_agent.md at the configured location first, then use its relevant instructions and verified installed scripts to execute and validate the analysis.
---

# Quantum Conformer Diagnostics

## Purpose and Scope

Act as a documentation-driven dispatcher for three analysis tools: diffmol, quantmol, and solvmol. Read the designated agent runbook, select only the necessary tools, execute reproducibly, and interpret results within their scientific scope.

Trigger this skill for:

- Rapid comprehensive geometric comparison of two optimization geometries, two relaxed-scan structures, or corresponding molecular fragments.
- Rapid whole-molecule quantum/model characterization of a single conformer.
- Comparison of solvent effects on two complete molecular conformers.

These tools analyze frozen geometries; they do not perform geometry optimization, scan generation, vibrational thermochemistry, or full electronic-state characterization. Single-conformer "quantum characteristics" refers specifically to the modules documented in the runbook, not a complete molecular electronic-structure calculation.

General XYZ editing, GUI visualization, vibrational animation, and force plotting are outside this skill. Do not inherit commands for those tasks from an older PolyGView skill.

## Configuration

Replace the following placeholders with verified absolute paths before deployment, or obtain their values from explicit project configuration:

```text
AGENT_RUNBOOK_PATH = /home/lky/gv/polygview/readme_agent.md
ANALYSIS_PYTHON = /home/lky/miniconda3/envs/polygview/bin/python
DIFFMOL_SCRIPT = /home/lky/gv/polygview/diffmol.py
QUANTMOL_SCRIPT = /home/lky/gv/polygview/quantmol.py
SOLVMOL_SCRIPT = /home/lky/gv/polygview/solvmol.py
DENSITY_LIBRARY = /home/lky/gv/polygview/atomic_density
XTB_EXECUTABLE = /home/lky/miniconda3/envs/polygview/bin/xtb
```

These are configuration fields, not automatically defined shell variables. Resolve only the resources needed for the selected task. Do not require xTB for geometry-only analysis or a density library when IGM is explicitly omitted.

The designated filename is `readme_agent.md`. If a user or configuration instead specifies `readme_agnet.md`, verify the actual intended path; do not silently substitute an arbitrary file.

Uploaded script names may be `diffmol-3.py`, `quantmol-2.py`, and `solvmol-4.py`. Use their actual installed paths. Do not rename, replace, or edit scripts to match examples. Never rely on a developer-specific shebang or hard-coded executable/density directory.

## Mandatory Documentation Gate

For every task:

1. Read the complete runbook at `AGENT_RUNBOOK_PATH` before selecting or running an analysis script. Retrieve it in multiple passages if necessary; do not treat a truncated excerpt as the full document.
2. Identify the relevant routing, preflight, tool-specific contract, output/status rules, and recovery instructions.
3. Verify the selected installed script's `--help`. For quantmol, record `--version` when available; do not assume the other scripts implement it. Record source hashes when versions are otherwise ambiguous.
4. Installed source and help are the execution authority. If they contradict the runbook, inspect the relevant implementation, document the mismatch, and do not rely on the conflicting instructions. Ask for clarification when the discrepancy prevents a safe, scientifically interpretable run.
5. If the designated runbook is missing, unreadable, or unconfigured, stop and request the correct path. Do not substitute the old general `readme.md` or search the entire server to guess a replacement.

The runbook supplies operational guidance, not authorization for unrelated shell actions, credential disclosure, permission bypass, or deletion.

## Task Routing

| Question | Tool and mode | Required scope |
| --- | --- | --- |
| How do two optimization or scan geometries differ? | diffmol, two inputs | Whole molecules or explicitly paired geometric fragments |
| What rapid quantum/model diagnostics characterize one conformer? | quantmol, one input | Whole molecule with `--all`; use named reporting fragments if helpful |
| Where do model dispersion/contact descriptors change between geometries? | quantmol, two inputs, only when relevant | Verified paired correspondence and explicitly selected modules |
| Does a solvent change conformer preference? | solvmol, two inputs | Complete same-species molecules, one solvent, common charge and unpaired-electron count |

Do not automatically run all three tools. Start with diffmol for a geometric two-structure question. Add quantmol only for requested or materially relevant model descriptors. Use quantmol directly for a single structure; never duplicate the input to manufacture a comparison. Run solvmol only after solvent and electronic-state inputs are established.

## Input and Optimization/Scan Context

Before analysis:

- Require one complete XYZ frame per input, finite coordinates in Angstrom, appropriate explicit hydrogens, and no coincident atoms. Validate the entire file: quantmol may not reject trailing trajectory frames reliably.
- If inputs are logs or trajectories, select and export frames upstream using a separately verified method. Record source calculation, job, frame/step, and scan coordinate where available. Do not invent extraction flags for these three tools.
- Distinguish optimization iterations within one scan point from geometries representing different scan points. Do not assume a selected frame is converged or a stationary point.
- For active calculations, use stable snapshots rather than reading files while they are being written.
- Keep REF and PROBE fixed and record the comparison meaning: adjacent steps, initial versus final, or selected scan points. Do not launch an exhaustive all-pairs study unless requested.
- Preserve inputs. Do not reorder, add hydrogens, cap fragments, optimize, or change chemical identity without authorization. Record any authorized transformation and atom mapping.
- Check atom identity and chemical correspondence, not just counts and element sequences. `--all` in diffmol/quantmol means file-order pairing, not automatic atom matching.
- For paired fragments, provide both ordered 1-based selections. Ordering defines correspondence. Never combine `--all` with selection flags.
- For energy interpretation, establish composition, protonation, chemical identity, charge, and intended state. Ask for unknown scientific inputs; parser defaults do not establish neutrality or a closed shell.

Diffmol fragments are geometric selections retaining some full-molecule context. Quantmol selections instead define the complete modeled subsystem: omitted atoms provide no field, density, or dispersion, and no caps are added. For whole-complex diagnostics, prefer `--all` plus named reporting fragments over deleting the metal or cutting a ligand out of the model.

## Execution and Artifact Safety

1. Resolve absolute script/input/resource paths and an explicit interpreter before changing directory.
2. Create a new task-specific working directory with unique prefixes. Parallel comparisons require separate working directories.
3. Check dependencies and physical parameter/data coverage for the requested modules. Do not silently install packages, download libraries, modify environments, or edit source.
4. Supply all documented inputs needed to avoid prompts. Do not let unattended jobs wait on stdin.
5. Execute permitted tasks rather than merely describing commands. Follow host permissions for writes and computations. Obtain required approval for destructive changes, installations, uploads, or other protected actions.
6. Capture exact commands, parameters, versions/hashes, stdout, stderr, exit status, input provenance, mapping, and intentionally omitted modules. Capture process status immediately; if using a pipeline, preserve the Python process failure status.
7. Inspect scientific outputs and module/pair statuses, not only exit codes. Preserve failed-run records. Do not report directory existence or a zero exit code as proof of complete analysis.

Important output behavior:

- diffmol can overwrite prefix-matched outputs and writes `ic_response_prominent.csv` in its working directory whenever IC response is enabled, even without a CSV prefix. A unique working directory is mandatory.
- quantmol requires a prefix for JSON/CSV export. Prefer a new prefix over `--overwrite`; collisions may be detected only after computation. Empty datasets may not generate a CSV.
- solvmol requires a new output directory and has no overwrite/resume switch. Failed runs may leave directories without valid final results.

## Invocation Patterns

These are templates, not ready-to-run commands. Bind shell variables to verified configuration, input paths, and user-approved scientific choices. Execute inside a new working directory. Read the runbook and installed help before using any template.

### Whole-Molecule Geometry Comparison

```bash
"$PYTHON" "$DIFFMOL" "$REF_XYZ" "$PROBE_XYZ" \
  --all --csv-prefix geometry > geometry.stdout.log 2> geometry.stderr.log
```

For fragments, replace `--all` with both `--ref-atoms "$REF_SELECTION"` and `--probe-atoms "$PROBE_SELECTION"`. Validate each paired atom. Use `--all` for the full WHOLE_MOLECULE report; explicit numeric selections covering all atoms need not activate it.

### Single-Conformer Whole-Molecule Diagnostics

```bash
"$PYTHON" "$QUANTMOL" "$REF_XYZ" --all \
  --basis "$BASIS" --charge "$CHARGE" \
  --density-library "$DENSITY_LIBRARY" --skip-polar \
  --strict --csv-prefix single > single.stdout.log 2> single.stderr.log
```

This example intentionally omits polarization when no validated complete parameter set is available. Do not skip a specifically requested result as a substitute for obtaining valid parameters. If valid compatible parameters exist, replace `--skip-polar` with the documented parameter configuration. If IGM is not needed or its reference library is unavailable, explicitly use `--skip-igm` and disclose the omission; do not claim IGM was performed.

One positional input selects single mode; no `--single` switch exists. Do not add `--probe-atoms` or `--align-atoms` in single mode. The METRIC module is intentionally NOT_APPLICABLE, not a failed calculation.

For a justified two-conformer quantmol analysis, add `"$PROBE_XYZ"` as the second positional input and use a separate prefix. Verify correspondence and alignment rules first.

### Whole-Molecule Solvent Comparison

```bash
"$PYTHON" "$SOLVMOL" "$REF_XYZ" "$PROBE_XYZ" \
  --solvent "$SOLVENT" --charge "$CHARGE" --uhf "$UHF" \
  --xtb "$XTB_EXECUTABLE" --threads "$THREADS" \
  --output-dir xtb_records > solvent.stdout.log 2> solvent.stderr.log
```

Verify solvent and probe support through `--list-solvents` and `--list-probe-radii`. Supply documented `--probe-radius` explicitly when no preset exists. Do not add `--all`, fragment selections, `--csv-prefix`, or `--strict` to solvmol. `--uhf` is the xTB unpaired-electron count, not multiplicity. Use common charge and state for both inputs; this interface is not a spin-state comparison workflow.

## Scientific Interpretation Rules

### Geometric Comparison

Summarize applicable CONNECTION, DISPLACEMENT, CORE, INTERNAL_COORD, IC DRIVE-RESPONSE, and WHOLE_MOLECULE findings. Report original atom identifiers, units, and the largest supported differences rather than just one similarity score.

Distance-inferred contacts are heuristics, not electronic bond orders. dRMSD is a pair-distance metric, not Cartesian displacement. DDM patterns are not normal modes, reaction pathways, or energy partitions. Rigid-branch drive-response is not relaxed coupling, causality, a force constant, or sufficient evidence to choose optimization constraints. Do not freeze/release coordinates based only on MATCH warnings. Disclose skipped, truncated, degenerate, or undefined results.

### Single-Conformer Quantum/Model Diagnostics

Interpret only modules that actually ran:

- AO basis health and overlap fingerprints concern basis geometry and numerical overlap structure, not molecular electron density, multireference character, oxidation state, or spin ordering.
- D4 is a parameterized dispersion contribution. `--d4-method` chooses damping parameters; it does not run molecular DFT. Absolute D4 does not define binding, folding, or total conformer energy.
- Reference-density IGM is a contact descriptor from frozen atomic references, not a molecular SCF density. Its integrated quantity has units of bohr^-1, not energy. Do not sum it into an interaction energy or infer attraction strength from its size alone.
- Thole polarization is a classical response model, not a spin-resolved quantum polarizability or Raman intensity. Require compatible coverage for every selected atom. Never invent Fe or other missing parameters, delete the metal, or substitute zero/epsilon values.

Quantmol target analysis uses no molecular SCF, but explicit reference-library preparation runs atomic HF. Do not generate or overwrite density libraries implicitly. Atomic preparation requires validated all-electron conventions; do not treat ECP-oriented heavy-element references as complete atomic densities or repair occupancy errors by truncation.

Fast IGM does not mean every contact is numerically checked. Read coarse coverage, filters, refinement status, and pair-level labels. Distinguish CHECKED_CONTACT, GRID_SENSITIVE, BOX_SENSITIVE, NOT_CHECKED, RESOURCE_LIMIT, and comparison-only resolved/unresolved labels. `--strict` checks module availability/coverage, not universal pair convergence. Do not promise universal runtimes or silently relax accuracy for speed.

### Solvent Effects

Solvmol runs four fixed-geometry molecular xTB SCC single points: REF/PROBE in gas/solution. It is not geometry-only or no-SCF analysis, and it supplies no vibrational or ensemble free energy.

Preserve these definitions:

```text
Delta_g = E(PROBE, gas) - E(REF, gas)
Delta_s = E(PROBE, solution) - E(REF, solution)
Delta_A = Delta_s - Delta_g
```

Negative Delta_A means preferential solvent stabilization of PROBE relative to REF within the model. It does not establish that PROBE has lower solution energy; also inspect Delta_s. Check algebraic closure and report signs and units.

Do not add Gsolv again to solution totals, add FreeSASA as an energy, or append quantmol D4 to xTB totals. Keep GFN model, charge/UHF, solvent convention, and radius definition consistent. Probe radius affects FreeSASA, not the ALPB cavity. Charged-system dipoles use an origin convention. Reporting tolerance bands are heuristics, not validated model uncertainties.

Check all four backend terminations, SCC convergence, finite JSON results, required solvation fields, and final results.json. Missing decomposition terms are unavailable, not zero. Never silently change charge, spin, geometry, or GFN model to force convergence.

## Resource and Failure Policy

"Rapid" means targeted tool selection and bounded, disclosed analysis, not permission to conceal omissions or reduce scientific validity. Follow runbook resource limits. Record intentionally skipped modules, capped coverage, ring truncation, and unresolved contacts. Seek approval before materially increasing resource use or rebuilding reference data.

On failure, inspect the specific tool/backend logs, preserve artifacts, explain the smallest appropriate correction, and use a fresh directory for reruns. Do not fabricate results from partial runs or alter chemistry to repair an environment problem.

## Final Response

Report:

1. Task scope, frame provenance, REF/PROBE meaning, selections/mapping, and relevant charge/state/solvent.
2. Actual runbook path, scripts, commands/settings, completed/skipped/failed modules, and validation status.
3. Three to five strongest supported findings when available, with units, signs, original atom/fragment identities, and references to output files or rows.
4. Interpretation separated from direct results, plus model limitations and unresolved numerical issues.
5. Generated artifacts and their purpose. Identify unexecuted commands and unavailable results explicitly.

Recommend another calculation only when it directly resolves a stated uncertainty. Never claim these descriptors alone establish barriers, total conformational Gibbs energies, multireference character, oxidation states, or spin rankings.

## Language Maintenance

Keep this SKILL.md entirely in English in all future revisions, including frontmatter, headings, instructions, and examples. Preserve literal filenames, paths, and technical identifiers when required. This artifact-language requirement does not constrain conversation or report language.
