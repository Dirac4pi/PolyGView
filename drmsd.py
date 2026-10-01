#!/home/lky/miniconda3/envs/polygview/bin/python
"""
Compare XYZ fragments using fixed REF topology and paired atom selections.
Author: Dirac4pi
env:polygview
"""

import argparse
import csv
import sys
from itertools import combinations, permutations
import numpy as np

#------------------------------------------------------------------------------
def _parse_indices(text, label, num_atoms):
  """Expand inclusive 1-based ranges while preserving input order."""
  text = text.strip()
  if text.lower() == 'c':
    return list(range(1, num_atoms + 1))
  indices = []
  for item in text.split(','):
    item = item.strip()
    if not item:
      raise ValueError(f'{label}: empty selection item.')
    if '-' in item:
      parts = item.split('-')
      if len(parts) != 2 or not all(p.strip().isdigit() for p in parts):
        raise ValueError(f'{label}: invalid range {item!r}.')
      first, last = (int(p.strip()) for p in parts)
      if first < 1 or last < first:
        raise ValueError(f'{label}: ranges must be positive and ascending.')
      indices.extend(range(first, last + 1))
    else:
      if not item.isdigit() or int(item) < 1:
        raise ValueError(f'{label}: invalid atom number {item!r}.')
      indices.append(int(item))
  if len(indices) != len(set(indices)):
    raise ValueError(f'{label}: duplicate atom numbers are not allowed.')
  return indices

#------------------------------------------------------------------------------
def read_xyz(filepath):
  """Read a single XYZ frame with finite Cartesian coordinates."""
  with open(filepath, encoding='utf-8') as f:
    lines = f.readlines()
  if not lines:
    raise ValueError(f'Empty XYZ: {filepath}')
  n = int(lines[0].strip())
  if n < 1 or len(lines) < n + 2:
    raise ValueError(f'Invalid or incomplete XYZ: {filepath}')
  symbols, coords = [], []
  for number, line in enumerate(lines[2:n+2], 1):
    parts = line.split()
    if len(parts) < 4:
      raise ValueError(f'{filepath}: incomplete atom {number}')
    symbols.append(parts[0].capitalize())
    coords.append([float(v.replace('D', 'E').replace('d', 'e')) for v in parts[1:4]])
  coords = np.asarray(coords)
  if not np.all(np.isfinite(coords)):
    raise ValueError(f'Non-finite coordinates: {filepath}')
  if any(line.strip() for line in lines[n+2:]):
    raise ValueError(f'{filepath}: provide a single-frame XYZ.')
  return symbols, coords

#------------------------------------------------------------------------------
def distance_metrics(a, b):
  """Calculate dRMSD and dMAXD without full pairwise matrices."""
  if a.shape != b.shape:
    raise ValueError('Coordinate shapes differ.')
  if len(a) < 2:
    return 0.0, 0.0
  ss, maximum, count = 0.0, 0.0, 0
  for i in range(len(a)-1):
    delta = np.linalg.norm(b[i+1:]-b[i], axis=1) - np.linalg.norm(a[i+1:]-a[i], axis=1)
    ss += float(delta @ delta)
    maximum = max(maximum, float(np.abs(delta).max()))
    count += len(delta)
  return float(np.sqrt(ss/count)), maximum

#------------------------------------------------------------------------------
def drmsd(coords_a, coords_b):
  return distance_metrics(coords_a, coords_b)[0]

#------------------------------------------------------------------------------
def dmaxd(coords_a, coords_b):
  return distance_metrics(coords_a, coords_b)[1]

#------------------------------------------------------------------------------
def infer_bonds(coords, radii):
  """Apply the original 1.15-times-covalent-radius-sum criterion."""
  bonds = set()
  for i in range(len(coords)-1):
    distances = np.linalg.norm(coords[i+1:] - coords[i], axis=1)
    if np.any(distances <= 1e-12):
      raise ValueError('Coincident selected atoms.')
    for offset in np.flatnonzero(distances <= 1.15*(radii[i] + radii[i+1:])):
      bonds.add((i, i+1+int(offset)))
  return bonds

#------------------------------------------------------------------------------
def dihedral_patterns(atoms, bonds):
  """Return satisfied chain and coordination edge patterns."""
  a, b, c, d = atoms
  edge = lambda i, j: tuple(sorted((i, j)))
  patterns = {
    'chain': [(edge(a, b), edge(b, c), edge(c, d))],
    'coordination': [
      (edge(a, b), edge(b, c), edge(b, d)),
      (edge(a, c), edge(b, c), edge(c, d)),
    ],
  }
  return {label: [p for p in options if all(e in bonds for e in p)]
    for label, options in patterns.items()}

#------------------------------------------------------------------------------
def coordination_key(atoms):
  """Identify the same plane pair, ignoring axis and plane ordering."""
  a, b, c, d = atoms
  return (tuple(sorted((b, c))), tuple(sorted((a, d))))

#------------------------------------------------------------------------------
def enumerate_coordinates(n, bonds):
  """Enumerate internal coordinates and deduplicate coordination plane pairs."""
  neighbors = [set() for _ in range(n)]
  for i, j in bonds:
    neighbors[i].add(j)
    neighbors[j].add(i)
  angles = []
  for j in range(n):
    angles.extend((i, j, k) for i, k in combinations(sorted(neighbors[j]), 2))
  torsions = set()
  for j, k in sorted(bonds):
    for i in neighbors[j] - {k}:
      for l in neighbors[k] - {j}:
        t = (i, j, k, l)
        if len(set(t)) == 4:
          torsions.add(min(t, t[::-1]))
  for center in range(n):
    for a, c, d in permutations(sorted(neighbors[center]), 3):
      t = (a, center, c, d)
      torsions.add(min(t, t[::-1]))
  unique, seen = [], set()
  for t in sorted(torsions):
    if dihedral_patterns(t, bonds)['coordination']:
      key = coordination_key(t)
      if key in seen:
        continue
      seen.add(key)
    unique.append(t)
  return {'bond': sorted(bonds), 'angle': sorted(angles), 'dihedral': unique}

#------------------------------------------------------------------------------
def coordinate_value(coords, atoms, kind, linear_tol):
  """Return lengths in Angstrom and angles in degrees."""
  p = coords[list(atoms)]
  if kind == 'bond':
    return float(np.linalg.norm(p[1]-p[0]))
  if kind == 'angle':
    u,v = p[0]-p[1], p[2]-p[1]
    if np.linalg.norm(u)*np.linalg.norm(v) <= 1e-24:
      raise ValueError('Zero-length angle arm')
    return float(np.degrees(np.arctan2(np.linalg.norm(np.cross(u,v)), u@v)))
  b0,b1,b2 = p[0]-p[1], p[2]-p[1], p[3]-p[2]
  norms = [np.linalg.norm(v) for v in (b0,b1,b2)]
  if min(norms) <= 1e-12:
    raise ValueError('Zero-length dihedral arm')
  axis = b1/norms[1]
  v = b0 - (b0@axis)*axis
  w = b2 - (b2@axis)*axis
  if np.linalg.norm(v)/norms[0] <= linear_tol or np.linalg.norm(w)/norms[2] <= linear_tol:
    raise ValueError('Undefined/near-collinear dihedral')
  return float(np.degrees(np.arctan2(np.cross(axis,v)@w, v@w)))

#------------------------------------------------------------------------------
def compare_internal(a, b, symbols, ref_ids, probe_ids, bonds, probe_bonds, linear_tol=1e-3):
  """Rank changes independently; keep REF bonds that fail the PROBE criterion."""
  rows, excluded = [], []
  candidates = enumerate_coordinates(len(a), bonds)
  for kind, coordinates in candidates.items():
    for atoms in coordinates:
      subtype = ''
      required_patterns = []
      if kind == 'dihedral':
        patterns = dihedral_patterns(atoms, bonds)
        subtype = 'coordination' if patterns['coordination'] else 'chain'
        required_patterns = patterns[subtype]
      base = dict(kind=kind, dihedral_type=subtype, ref_atoms='-'.join(str(ref_ids[i]) for i in atoms),
        probe_atoms='-'.join(str(probe_ids[i]) for i in atoms),
        elements='-'.join(symbols[i] for i in atoms))
      try:
        av = coordinate_value(a, atoms, kind, linear_tol)
        bv = coordinate_value(b, atoms, kind, linear_tol)
      except ValueError as exc:
        excluded.append(dict(base, reason=str(exc)))
        continue
      delta = bv-av
      if kind == 'dihedral':
        delta = (delta+180.0)%360.0-180.0
      if kind == 'dihedral':
        # Validate the same REF topology type in PROBE.
        valid_b = any(all(e in probe_bonds for e in p) for p in required_patterns)
      else:
        valid_b = all(tuple(sorted((i,j))) in probe_bonds for i,j in zip(atoms, atoms[1:]))
      rows.append(dict(base, unit='angstrom' if kind=='bond' else 'degree',
        ref_value=av, probe_value=bv, delta=delta, abs_delta=abs(delta),
        probe_bonds_ok=valid_b))
  ranked = {kind: sorted([r for r in rows if r['kind']==kind],
    key=lambda r:(-r['abs_delta'], r['ref_atoms'])) for kind in candidates}
  return ranked, excluded

#------------------------------------------------------------------------------
def write_csv(filename, rows, fields):
  with open(filename, 'w', newline='', encoding='utf-8') as f:
    writer = csv.DictWriter(f, fieldnames=fields)
    writer.writeheader()
    writer.writerows(rows)


#==============================================================================
def main(argv=None):
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument('ref_xyz')
  parser.add_argument('probe_xyz')
  parser.add_argument('--ref-atoms', help='1-based selection, e.g. 1-9,16-23,25')
  parser.add_argument('--probe-atoms', help='Paired in the exact order of REF selection')
  parser.add_argument('--all', action='store_true', help='Compare all atoms in file order without prompts')
  parser.add_argument('--top', type=int, default=10)
  parser.add_argument('--csv-prefix', help='Write PREFIX_all.csv, PREFIX_top.csv, PREFIX_excluded.csv')
  parser.add_argument('--linear-tol', type=float, default=1e-3, help='Minimum sine of either dihedral adjacent angle')
  args = parser.parse_args(argv)
  if args.top < 1 or not 0 < args.linear_tol < 1:
    parser.error('--top must be positive; --linear-tol must lie between 0 and 1')
  if args.all and (args.ref_atoms is not None or args.probe_atoms is not None):
    parser.error('--all cannot be combined with fragment selections')
  if (args.ref_atoms is None) != (args.probe_atoms is None):
    parser.error('Provide both --ref-atoms and --probe-atoms')
  sa, ca = read_xyz(args.ref_xyz)
  sb, cb = read_xyz(args.probe_xyz)
  if args.all:
    ra, pb = list(range(1,len(sa)+1)), list(range(1,len(sb)+1))
  else:
    print(f'Reference XYZ (first argument): {args.ref_xyz}')
    print(f'Probe XYZ (second argument): {args.probe_xyz}')
    rt = args.ref_atoms if args.ref_atoms is not None else input('Enter REF atoms (1-based; e.g. 1-9,16-23,25): ')
    pt = args.probe_atoms if args.probe_atoms is not None else input('Enter PROBE atoms (1-based; same order as REF): ')
    ra = _parse_indices(rt, 'ref_indices', len(sa))
    pb = _parse_indices(pt, 'probe_indices', len(sb))
  if len(ra) != len(pb) or not ra:
    raise ValueError('Provide equal, nonempty paired selections.')
  for i,j in zip(ra,pb):
    if i > len(sa) or j > len(sb):
      raise ValueError(f'Atom index out of range: REF {i}, PROBE {j}')
    if sa[i-1] != sb[j-1]:
      raise ValueError(f'Element mismatch: REF {i} {sa[i-1]}, PROBE {j} {sb[j-1]}')
  a,b = ca[np.array(ra)-1], cb[np.array(pb)-1]
  symbols = [sa[i-1] for i in ra]
  import qcelemental as qcel
  radii = np.array([float(qcel.covalentradii.get(s, units='angstrom')) for s in symbols])
  if not np.all(np.isfinite(radii)) or np.any(radii <= 0):
    raise ValueError('Invalid covalent radii')
  bonds, probe_bonds = infer_bonds(a,radii), infer_bonds(b,radii)
  print('Atom mapping (REF:PROBE):', ', '.join(f'{i}:{j}' for i,j in zip(ra,pb)))
  print('Topology: fixed REF; d <= 1.15 * (r_cov_i + r_cov_j), Angstrom')
  print(f'REF bonds: {len(bonds)}; missing in PROBE: {len(bonds-probe_bonds)}; PROBE-only: {len(probe_bonds-bonds)}')
  for label, pairs in [('REF bonds failing PROBE criterion', bonds-probe_bonds), ('PROBE-only bonds (not enumerated)', probe_bonds-bonds)]:
    if pairs:
      print(label+':', ', '.join(f'REF {ra[i]}-{ra[j]} / PROBE {pb[i]}-{pb[j]}' for i,j in sorted(pairs)))
  rms, maximum = distance_metrics(a,b)
  print(f'Distance RMSD: {rms:.6f} Angstroms')
  print(f'Distance MAXD: {maximum:.6f} Angstroms')
  ranked, excluded = compare_internal(a,b,symbols,ra,pb,bonds,probe_bonds,args.linear_tol)
  all_rows, top_rows = [], []
  for kind, rows in ranked.items():
    print(f'\n{kind.upper()}: top {min(args.top,len(rows))} / {len(rows)} valid')
    print('Rank  REF atoms         PROBE atoms       Elements       REF value    PROBE value    Delta    |delta|   bonds     Type')
    for rank,row in enumerate(rows,1):
      record = dict(rank=rank, **row)
      all_rows.append(record)
      if rank <= args.top:
        top_rows.append(record)
        print(f"{rank:4d}  {row['ref_atoms']:17s} {row['probe_atoms']:17s} {row['elements']:12s} {row['ref_value']:11.6f} {row['probe_value']:13.6f} {row['delta']:10.6f} {row['abs_delta']:10.6f} {'  OK  ' if row['probe_bonds_ok'] else 'WARN  '} {row['dihedral_type'] or '-'}")
  print(f'\nExcluded coordinates: {len(excluded)}')
  for row in excluded:
    print(f"  {row['kind']} REF {row['ref_atoms']} / PROBE {row['probe_atoms']}: {row['reason']} (type={row['dihedral_type'] or '-'})")
  if args.csv_prefix:
    fields = ['rank','kind','dihedral_type','ref_atoms','probe_atoms','elements','unit','ref_value','probe_value','delta','abs_delta','probe_bonds_ok']
    write_csv(args.csv_prefix+'_all.csv',all_rows,fields)
    write_csv(args.csv_prefix+'_top.csv',top_rows,fields)
    write_csv(args.csv_prefix+'_excluded.csv',excluded,['kind','dihedral_type','ref_atoms','probe_atoms','elements','reason'])
  return ranked, excluded


if __name__ == '__main__':
  try:
    main()
  except (ValueError, OSError, ImportError, EOFError) as exc:
    print(f'{sys.argv[0]}: {exc}', file=sys.stderr)
    sys.exit(2)
