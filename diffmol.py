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


#------------------------------------------------------------------------------
# Bondi-style van der Waals radii in Angstrom for the supported XB elements.
VDW_RADII = {'N': 1.55, 'O': 1.52, 'S': 1.80,
  'Cl': 1.75, 'Br': 1.85, 'I': 1.98}


def interaction_neighbors(n, bonds):
  """Build a full-molecule neighbor graph for chemical-context checks."""
  neighbors = [set() for _ in range(n)]
  for i, j in bonds:
    neighbors[i].add(j)
    neighbors[j].add(i)
  return neighbors


def contact_angle(coords, i, j, k):
  """Return the angle i-j-k in degrees, or None for a zero-length arm."""
  u, v = coords[i] - coords[j], coords[k] - coords[j]
  if min(np.linalg.norm(u), np.linalg.norm(v)) <= 1e-12:
    return None
  return float(np.degrees(np.arctan2(np.linalg.norm(np.cross(u, v)), u @ v)))


def detect_contacts(coords, symbols, selected_ids, full_bonds, options):
  """Screen internal geometric candidates; do not claim electronic bonding."""
  selected = [i - 1 for i in selected_ids]
  local = {original: paired for paired, original in enumerate(selected)}
  neighbors = interaction_neighbors(len(symbols), full_bonds)
  acceptors = []
  for a in selected:
    if symbols[a] not in options.acceptors:
      continue
    # Exclude tetra-coordinate N; XYZ cannot resolve all acceptor chemistry.
    if symbols[a] == 'N' and len(neighbors[a]) >= 4:
      continue
    acceptors.append(a)
  hbonds, xbonds = {}, {}
  for h in selected:
    if symbols[h] != 'H' or len(neighbors[h]) != 1:
      continue
    d = next(iter(neighbors[h]))
    if d not in local or symbols[d] not in options.hb_donors:
      continue
    for a in acceptors:
      if a in (d, h) or a in neighbors[d] or a in neighbors[h]:
        continue
      da = float(np.linalg.norm(coords[d] - coords[a]))
      ha = float(np.linalg.norm(coords[h] - coords[a]))
      angle = contact_angle(coords, d, h, a)
      if angle is None:
        continue
      if da <= options.hb_da and angle >= options.hb_angle:
        hbonds[(local[d], local[h], local[a])] = {
          'primary_distance_A': da, 'secondary_distance_A': ha, 'angle_deg': angle}
  for x in selected:
    if symbols[x] not in ('Cl', 'Br', 'I') or len(neighbors[x]) != 1:
      continue
    r = next(iter(neighbors[x]))
    if r not in local:
      continue
    for a in acceptors:
      if a in (r, x) or a in neighbors[x] or a in neighbors[r]:
        continue
      xa = float(np.linalg.norm(coords[x] - coords[a]))
      rx = float(np.linalg.norm(coords[r] - coords[x]))
      cutoff = options.xb_scale * (VDW_RADII[symbols[x]] + VDW_RADII[symbols[a]])
      angle = contact_angle(coords, r, x, a)
      if angle is None:
        continue
      if xa <= cutoff and angle >= options.xb_angle:
        xbonds[(local[r], local[x], local[a])] = {
          'primary_distance_A': xa, 'secondary_distance_A': rx, 'angle_deg': angle}
  return {'hydrogen': hbonds, 'halogen': xbonds}


def interaction_report(a, b, symbols, ref_ids, probe_ids, bonds, probe_bonds,
  full_ref_coords, full_probe_coords, full_ref_symbols, full_probe_symbols,
                       full_ref_bonds, full_probe_bonds, options):
  """Print all internal HB/XB candidates and presence changes under atom mapping."""
  covalent = []
  print('\n==============================CONNECTION==================================')
  print('COVALENT CONNECTIVITY DIFFERENCES (distance criterion)')
  print('REF atoms / PROBE atoms   Elements   REF(A) PROBE(A) REF bonded PROBE bonded')
  for i, j in sorted(bonds ^ probe_bonds):
    row = dict(ref_atoms=f'{ref_ids[i]}-{ref_ids[j]}',
      probe_atoms=f'{probe_ids[i]}-{probe_ids[j]}',
      elements=f'{symbols[i]}-{symbols[j]}',
      ref_distance_A=float(np.linalg.norm(a[i]-a[j])),
      probe_distance_A=float(np.linalg.norm(b[i]-b[j])),
      ref_bonded=(i, j) in bonds, probe_bonded=(i, j) in probe_bonds,
      status='lost' if (i, j) in bonds else 'gained')
    covalent.append(row)
    print(f"{row['ref_atoms']} / {row['probe_atoms']}  {row['elements']}  "
      f"{row['ref_distance_A']:.6f} {row['probe_distance_A']:.6f} "
      f"{row['ref_bonded']} {row['probe_bonded']} {row['status']}")
  if not covalent:
    print('None.')
  ref_contacts = detect_contacts(full_ref_coords, full_ref_symbols, ref_ids,
    full_ref_bonds, options)
  probe_contacts = detect_contacts(full_probe_coords, full_probe_symbols, probe_ids,
    full_probe_bonds, options)
  print('\nNONCOVALENT GEOMETRIC CANDIDATES (not electronic bond assignments)')
  print(f'HB: donors={options.hb_donors}; acceptors={options.acceptors}; '
    f'D-A <= {options.hb_da:g} A; D-H-A >= {options.hb_angle:g} deg')
  print(f'XB: X=Cl,Br,I; acceptors={options.acceptors}; '
    f'X-A <= {options.xb_scale:g} * VDW sum; R-X-A >= {options.xb_angle:g} deg')
  print('All defining atoms must be selected. Covalent/1-3 contacts are excluded.')
  print('HB distances: D-A, H-A. XB distances: X-A, R-X.')
  contacts = []
  for kind in ('hydrogen', 'halogen'):
    for label, mapping, ids in (('REF', ref_contacts[kind], ref_ids),
      ('PROBE', probe_contacts[kind], probe_ids)):
      print(f'\n{label} {kind.upper()} CANDIDATES: {len(mapping)}')
      for key, value in sorted(mapping.items()):
        names = '-'.join(str(ids[i]) for i in key)
        elements = '-'.join(symbols[i] for i in key)
        print(f"{names:18s} {elements:12s} "
          f"{value['primary_distance_A']:.6f} {value['secondary_distance_A']:.6f} "
          f"{value['angle_deg']:.3f}")
      if not mapping:
        print('None.')
    lost = ref_contacts[kind].keys() - probe_contacts[kind].keys()
    gained = probe_contacts[kind].keys() - ref_contacts[kind].keys()
    print(f'\n{kind.upper()} DIFFERENCES: lost={len(lost)}, gained={len(gained)}')
    for key in sorted(ref_contacts[kind].keys() | probe_contacts[kind].keys()):
      ref_data, probe_data = ref_contacts[kind].get(key), probe_contacts[kind].get(key)
      status = 'lost' if probe_data is None else ('gained' if ref_data is None else 'retained')
      row = dict(kind=kind, status=status,
        ref_atoms='-'.join(str(ref_ids[i]) for i in key),
        probe_atoms='-'.join(str(probe_ids[i]) for i in key),
        elements='-'.join(symbols[i] for i in key),
        ref_present=ref_data is not None, probe_present=probe_data is not None)
      # Evaluate both geometries even when only one passes the screening criteria.
      for label, xyz in (('ref', a), ('probe', b)):
        d, h, acceptor = key
        pair = (d, acceptor) if kind == 'hydrogen' else (h, acceptor)
        secondary = (h, acceptor) if kind == 'hydrogen' else (d, h)
        row[label+'_primary_distance_A'] = float(np.linalg.norm(xyz[pair[0]]-xyz[pair[1]]))
        row[label+'_secondary_distance_A'] = float(np.linalg.norm(xyz[secondary[0]]-xyz[secondary[1]]))
        angle = contact_angle(xyz, d, h, acceptor)
        row[label+'_angle_deg'] = '' if angle is None else angle
      contacts.append(row)
      if status != 'retained':
        print(f"{status.upper():6s} REF {row['ref_atoms']} / PROBE {row['probe_atoms']} "
          f"{row['elements']} primary(A)={row['ref_primary_distance_A']:.6f} -> "
          f"{row['probe_primary_distance_A']:.6f}; angle(deg)="
          f"{row['ref_angle_deg']} -> {row['probe_angle_deg']}")
    if not lost and not gained:
      print('No presence differences.')
  print('==========================================================================\n')
  return covalent, contacts


#------------------------------------------------------------------------------
# Metal recognition covers alkali, alkaline-earth, transition, f-block,
# and common post-transition metals. Metalloids are not included.
METALS = set(('Li Na K Rb Cs Fr Be Mg Ca Sr Ba Ra Al Ga In Tl Sn Pb Bi Po '
  'Sc Ti V Cr Mn Fe Co Ni Cu Zn Y Zr Nb Mo Tc Ru Rh Pd Ag Cd '
  'Hf Ta W Re Os Ir Pt Au Hg La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu '
  'Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv').split())


def rdkit_graph(symbols, coords, bonds):
  """Build a connectivity graph without assuming chemically valid bond orders."""
  from rdkit import Chem
  graph = Chem.RWMol()
  for symbol in symbols:
    atom = Chem.Atom(symbol)
    atom.SetNoImplicit(True)
    graph.AddAtom(atom)
  for i, j in sorted(bonds):
    graph.AddBond(i, j, Chem.BondType.SINGLE)
  mol = graph.GetMol()
  conf = Chem.Conformer(len(symbols))
  conf.Set3D(True)
  for i, point in enumerate(coords):
    conf.SetAtomPosition(i, tuple(float(x) for x in point))
  mol.AddConformer(conf)
  mol.UpdatePropertyCache(strict=False)
  Chem.GetSymmSSSR(mol)
  return mol


def find_chiral_centers(symbols, coords, bonds):
  """Use RDKit to find centers and configurations on the inferred graph.

  XYZ supplies no bond orders: all inferred edges are SINGLE bonds.
  Labels refer to this connectivity model, not a certified chemical topology.
  """
  from rdkit import Chem
  from rdkit.Chem import rdCIPLabeler
  mol = rdkit_graph(symbols, coords, bonds)
  Chem.SetHybridization(mol)
  Chem.AssignStereochemistryFrom3D(mol, confId=0, replaceExistingTags=True)
  centers = Chem.FindMolChiralCenters(mol, force=True,
    includeUnassigned=True, includeCIP=True, useLegacyImplementation=False)
  rdCIPLabeler.AssignCIPLabels(mol)
  return {i: mol.GetAtomWithIdx(i).GetProp('_CIPCode')
    if mol.GetAtomWithIdx(i).HasProp('_CIPCode') else label
    for i, label in centers}


def chirality_report(sa, ca, sb, cb, ra, pb, fa, fb, args):
  """List full-system centers; use only explicitly supplied atom correspondence."""
  ref = find_chiral_centers(sa, ca, fa)
  probe = find_chiral_centers(sb, cb, fb)
  mapping = {i-1: j-1 for i, j in zip(ra, pb)}
  rows, used_probe = [], set()
  print('\n=================================CORE=====================================')
  print('CHIRAL CENTERS (RDKit; complete REF/PROBE structures)')
  print('R/S labels use the inferred SINGLE-bond graph; ? means unassigned.')
  print('REF atom  PROBE atom  Element  REF config  PROBE config')
  for i in sorted(set(ref) | {r for r,p in mapping.items() if p in probe}):
    j = mapping.get(i)
    if j is not None:
      used_probe.add(j)
    rows.append(dict(ref_atom=i+1, probe_atom=j+1 if j is not None else '',
      element=sa[i], ref_configuration=ref.get(i,'-'),
      probe_configuration=probe.get(j,'-') if j is not None else 'unmapped'))
  for j in sorted(set(probe)-used_probe):
    rows.append(dict(ref_atom='', probe_atom=j+1, element=sb[j],
      ref_configuration='unmapped', probe_configuration=probe[j]))
  for row in rows:
    print(f"{str(row['ref_atom'] or '-'):>8s}  {str(row['probe_atom'] or '-'):>10s}  "
      f"{row['element']:7s}  {row['ref_configuration']:10s}  {row['probe_configuration']}")
  if not rows:
    print('No chiral centers identified by RDKit under this connectivity model.')
  return rows


def ideal_coordination():
  """Provide common CN 2-6 unit-direction templates, not a ChemEnv CSM."""
  planar = lambda n: np.column_stack((np.cos(2*np.pi*np.arange(n)/n),
    np.sin(2*np.pi*np.arange(n)/n), np.zeros(n)))
  equator = planar(3)
  axis = np.array([[0.,0,1],[0,0,-1]])
  return {
    2: {'linear': np.array([[1.,0,0],[-1,0,0]])},
    3: {'trigonal_planar': planar(3)},
    4: {'tetrahedral': np.array([[1.,1,1],[1,-1,-1],[-1,1,-1],[-1,-1,1]])/np.sqrt(3),
      'square_planar': planar(4)},
    5: {'trigonal_bipyramidal': np.vstack((equator, axis)),
      'square_pyramidal': np.vstack((planar(4),axis[:1]))},
    6: {'octahedral': np.vstack((np.eye(3),-np.eye(3))),
      'trigonal_prismatic': np.vstack((equator + [0,0,np.sqrt(3)/2],
        equator - [0,0,np.sqrt(3)/2])) / np.sqrt(1.75)},
  }


def fit_unit_geometry(vectors, templates):
  """Exhaustively fit CN<=6 permutations using proper Kabsch rotations.

  Score is RMS Euclidean distance between unit directions, not CSM.
  Angular RMS describes the chosen least-chord-distance mapping.
  """
  unit = vectors / np.linalg.norm(vectors, axis=1)[:, None]
  results = []
  for name, ideal in templates.items():
    best = None
    for perm in permutations(range(len(unit))):
      target = ideal[list(perm)]
      u, _, vt = np.linalg.svd(unit.T @ target)
      correction = np.eye(3)
      correction[-1,-1] = 1 if np.linalg.det(u @ vt) >= 0 else -1
      rotated = unit @ (u @ correction @ vt)
      score = float(np.sqrt(np.mean(np.sum((rotated-target)**2, axis=1))))
      angles = np.degrees(np.arccos(np.clip(np.sum(rotated*target, axis=1), -1, 1)))
      candidate = (score, float(np.sqrt(np.mean(angles**2))))
      if best is None or candidate < best:
        best = candidate
    results.append((name, *best))
  return sorted(results, key=lambda x: (x[1], x[0]))


def shape_interpretation(fits):
  """Use transparent reporting thresholds, not universal chemistry cutoffs."""
  if not fits:
    return 'unsupported coordination number'
  best = fits[0]
  angle = best[2]
  severity = ('near-ideal' if angle <= 5 else 'mildly distorted' if angle <= 10
    else 'moderately distorted' if angle <= 20 else 'strongly distorted')
  ambiguous = len(fits)>1 and (fits[1][1]-best[1] < 0.03 or abs(fits[1][2]-best[2]) < 2.)
  if ambiguous:
    return f'{severity}; mixed/intermediate geometry between {best[0]} and {fits[1][0]}'
  return f'{severity}; clear preference for {best[0]}'


def metal_report(sa, ca, sb, cb, ra, pb, args):
  """Use complete coordination shells, including ligands outside selected fragments."""
  import qcelemental as qcel
  templates, rows = ideal_coordination(), []
  shells = {'REF': {}, 'PROBE': {}}
  print('\nMETAL COORDINATION (selected centers; full-molecule shells)')
  print(f'Neighbor criterion: d <= {args.coord_scale:g} * covalent-radius sum.')
  print('shape_rms: dimensionless unit-direction mismatch, 0 is ideal; angle_rms: directional RMS in degrees.')
  print('Direction scores exclude bond-length stretching; NOT ChemEnv CSM.')
  print('Heuristic angle_rms bands: <=5 near-ideal; <=10 mild; <=20 moderate; >20 strong distortion.')
  print('Close fits (shape gap <0.03 or angle gap <2 deg) are reported as mixed geometry.')
  for label, symbols, coords, ids in (('REF', sa, ca, ra), ('PROBE', sb, cb, pb)):
    radii = np.array([qcel.covalentradii.get(x, units='angstrom') for x in symbols])
    for number in ids:
      m = number-1
      if symbols[m] not in METALS:
        continue
      distances = np.linalg.norm(coords-coords[m], axis=1)
      neighbors = [i for i in range(len(coords)) if i != m and distances[i] > 1e-12
        and distances[i] <= args.coord_scale*(radii[m]+radii[i])]
      shells[label][m] = set(neighbors)
      length = distances[neighbors]
      fits = fit_unit_geometry(coords[neighbors]-coords[m], templates[len(neighbors)]) if len(neighbors) in templates else []
      row = dict(structure=label, atom=number, element=symbols[m], CN=len(neighbors),
        neighbors=','.join(str(i+1) for i in neighbors),
        distances_A=','.join(f'{distances[i]:.6f}' for i in neighbors),
        mean_distance_A=float(length.mean()) if len(length) else '',
        distance_std_A=float(length.std()) if len(length) else '',
        best_model=fits[0][0] if fits else 'unsupported_CN',
        shape_rms=fits[0][1] if fits else '', angular_rms_deg=fits[0][2] if fits else '',
        interpretation=shape_interpretation(fits),
        model_scores=';'.join(f'{name}:{score:.8f}:{angular:.6f}' for name,score,angular in fits))
      rows.append(row)
      print(f"{label} {symbols[m]}({number}) CN={len(neighbors)} "
        f"neighbors={row['neighbors']} distances(A)={row['distances_A']}")
      for name, score, angular in fits:
        print(f'  {name}: shape_rms={score:.6f}, angle_rms={angular:.3f} deg')
      print('  Interpretation:', row['interpretation'])
  mapping = {i-1:j-1 for i,j in zip(ra,pb)}
  for i,j in mapping.items():
    if i not in shells['REF'] and j not in shells['PROBE']:
      continue
    ref_shell, probe_shell = shells['REF'].get(i,set()), shells['PROBE'].get(j,set())
    mapped = {mapping[x] for x in ref_shell if x in mapping}
    covered = all(x in mapping for x in ref_shell) and all(x in mapping.values() for x in probe_shell)
    print(f'Metal correspondence REF {i+1} / PROBE {j+1}: CN {len(ref_shell)} -> {len(probe_shell)}')
    ref_row = next((r for r in rows if r['structure']=='REF' and r['atom']==i+1),None)
    probe_row = next((r for r in rows if r['structure']=='PROBE' and r['atom']==j+1),None)
    if ref_row and probe_row and len(ref_shell)==len(probe_shell) and ref_row['model_scores'] and probe_row['model_scores']:
      ref_scores = {x.split(':')[0]:tuple(map(float,x.split(':')[1:])) for x in ref_row['model_scores'].split(';')}
      probe_scores = {x.split(':')[0]:tuple(map(float,x.split(':')[1:])) for x in probe_row['model_scores'].split(';')}
      model = ref_row['best_model']
      ds = probe_scores[model][0]-ref_scores[model][0]
      da = probe_scores[model][1]-ref_scores[model][1]
      trend = 'farther from' if ds>0.01 else 'closer to' if ds<-0.01 else 'similar fit to'
      print(f'  PROBE is {trend} REF-best {model}: delta shape_rms={ds:+.6f}, delta angle_rms={da:+.3f} deg.')
      print('  This comparison uses the same template; it does not establish a reaction pathway.')

    if covered:
      print('  lost PROBE-mapped neighbors:', ','.join(str(x+1) for x in sorted(mapped-probe_shell)) or 'none')
      print('  gained PROBE neighbors:', ','.join(str(x+1) for x in sorted(probe_shell-mapped)) or 'none')
    else:
      print('  Shell correspondence incomplete: some ligand atoms are outside the selected mapping.')
  if not rows:
    print('No selected metal centers.')
  return rows


def canonical_ring(ring):
  """Preserve adjacency while making cyclic ordering deterministic."""
  ring = tuple(ring)
  choices = []
  for seq in (ring, ring[::-1]):
    pivot = seq.index(min(seq))
    choices.append(seq[pivot:] + seq[:pivot])
  return min(choices)


def all_simple_rings(neighbors, max_size=0, max_count=10000, search_limit=1000000):
  """Enumerate undirected simple cycles; never silently claim exhaustive output."""
  rings, count, truncated = set(), 0, False
  for start in range(len(neighbors)):
    stack = [(start, [start])]
    while stack:
      vertex, path = stack.pop()
      count += 1
      if count > search_limit:
        truncated = True
        break
      for nxt in sorted(neighbors[vertex]):
        if nxt == start and len(path) >= 3:
          rings.add(canonical_ring(path))
          if len(rings) >= max_count:
            truncated = True
            break
        elif nxt > start and nxt not in path and (not max_size or len(path) < max_size):
          stack.append((nxt, path+[nxt]))
      if truncated:
        break
    if truncated:
      break
  return sorted(rings), truncated


def puckering(coords):
  """Calculate Cremer-Pople modes and independent best-fit-plane deviations.

  Normal convention: cross(sum(r*sin(t)), sum(r*cos(t))) with ordered atoms.
  Phase is undefined when the corresponding amplitude approaches zero.
  """
  p = np.asarray(coords) - np.mean(coords,axis=0)
  n = len(p)
  _, _, vt = np.linalg.svd(p,full_matrices=False)
  plane_z = p @ vt[-1]
  base = dict(plane_rms_A=float(np.sqrt(np.mean(plane_z**2))),
    plane_max_A=float(np.max(np.abs(plane_z))))
  if n == 3:
    return dict(base, Q_A=0., theta_deg='', modes={}, status='three_member_planar')
  t = 2*np.pi*np.arange(n)/n
  r1 = np.sum(p*np.sin(t)[:,None],axis=0)
  r2 = np.sum(p*np.cos(t)[:,None],axis=0)
  normal = np.cross(r1,r2)
  if np.linalg.norm(normal) <= 1e-10 * max(np.linalg.norm(r1)*np.linalg.norm(r2),1.):
    return dict(base,Q_A='',theta_deg='',modes={},status='degenerate_CP_plane')
  z = p @ (normal/np.linalg.norm(normal))
  modes = {}
  for m in range(2,(n-1)//2+1):
    cosine = float(np.sqrt(2/n)*np.sum(z*np.cos(m*t)))
    sine = float(-np.sqrt(2/n)*np.sum(z*np.sin(m*t)))
    amplitude = float(np.hypot(cosine,sine))
    phase = float(np.degrees(np.arctan2(sine,cosine))%360) if amplitude > 1e-8 else ''
    modes[m] = dict(q_A=amplitude,phi_deg=phase)
  if n%2 == 0:
    modes[n//2] = dict(q_A=float(np.sum(z*(-1.)**np.arange(n))/np.sqrt(n)),phi_deg='')
  Q = float(np.sqrt(np.sum(z*z)))
  theta = float(np.degrees(np.arccos(np.clip(modes[3]['q_A']/Q,-1,1)))) if n==6 and Q>1e-8 else ''
  return dict(base,Q_A=Q,theta_deg=theta,modes=modes,status='ok')


def ring_set(symbols, coords, ids, full_bonds, args):
  """Find internal nonmetal rings; metal-containing cycles are opt-in."""
  selected = [x-1 for x in ids]
  local = {x:i for i,x in enumerate(selected)}
  bonds = {tuple(sorted((local[i],local[j]))) for i,j in full_bonds
    if i in local and j in local and (args.include_metal_rings or
      (symbols[i] not in METALS and symbols[j] not in METALS))}
  neighbors = interaction_neighbors(len(ids),bonds)
  if args.ring_mode == 'all':
    return all_simple_rings(neighbors,args.ring_max_size,args.ring_max_count)
  from rdkit import Chem
  mol = rdkit_graph([symbols[i] for i in selected],coords[selected],bonds)
  rings = {canonical_ring(tuple(int(x) for x in ring)) for ring in Chem.GetSymmSSSR(mol)}
  rings = sorted(r for r in rings if not args.ring_max_size or len(r)<=args.ring_max_size)
  return rings[:args.ring_max_count], len(rings)>args.ring_max_count


def ring_interpretation(data, size, planar_cutoff):
  """Explain dimensional puckering and six-member shape with heuristic bands."""
  rms = data['plane_rms_A']
  level = ('approximately planar' if rms <= planar_cutoff else 'weakly puckered' if rms <= 0.10
    else 'moderately puckered' if rms <= 0.20 else 'strongly puckered')
  if data['status']=='degenerate_CP_plane':
    return level + '; CP reference plane is degenerate'
  if size==3:
    return 'three ring atoms define a plane; puckering is not an independent coordinate'
  q = data['Q_A']
  text = f'{level}; CP RMS height Q/sqrt(N)={q/np.sqrt(size):.4f} A'
  if q/np.sqrt(size) <= planar_cutoff:
    return text + '; amplitude too small for a stable phase/shape label'
  if size==6 and isinstance(data['theta_deg'],float):
    theta = data['theta_deg']
    if theta<=15 or theta>=165:
      text += '; chair-like'
    elif 75<=theta<=105:
      text += '; boat/twist-boat family (theta alone does not distinguish them)'
    else:
      text += '; intermediate six-membered pucker'
  elif size==5:
    text += '; phi locates the pseudorotation phase, not puckering strength'
  return text


def brief_number(value):
  return f'{value:.4f}' if isinstance(value,(float,np.floating)) else str(value or 'undefined')


def ring_report(sa, ca, sb, cb, ra, pb, fa, fb, args):
  """Match cyclic adjacency under the supplied atom map, not only atom sets."""
  ar, at = ring_set(sa,ca,ra,fa,args)
  br, bt = ring_set(sb,cb,pb,fb,args)
  print(f'\nRING PUCKERING: mode={args.ring_mode}; REF={len(ar)}, PROBE={len(br)}')
  print('Plane RMS/max (A): average/largest departure from the best-fit plane.')
  print('Q (A): total CP out-of-plane amplitude; Q/sqrt(N) permits ring-size-normalized comparison.')
  print('q_m (A): mode amplitude; phi_m (deg): phase, not strength; theta (6-rings): chair vs equatorial pucker.')
  print(f'Heuristic plane RMS bands: <={args.ring_planar_rms:g} planar; <=0.10 weak; <=0.20 moderate; >0.20 strong.')
  if at or bt:
    print('WARNING: ring count/search limit reached; ring output and presence differences may be incomplete.')
  if args.ring_max_size:
    print(f'Only rings of size <= {args.ring_max_size} were requested.')
  if not args.include_metal_rings:
    print('Metal-containing cycles excluded (use --include-metal-rings to include).')
  rows, modes = [], []
  a,b = ca[np.array(ra)-1],cb[np.array(pb)-1]
  for ring in sorted(set(ar)|set(br)):
    status = 'retained' if ring in ar and ring in br else ('lost' if ring in ar else 'gained')
    pa,pb_data = puckering(a[list(ring)]),puckering(b[list(ring)])
    # Suppress practically unstable phases when total per-atom puckering is tiny.
    for data in (pa,pb_data):
      if isinstance(data['Q_A'],float) and data['Q_A']/np.sqrt(len(ring)) <= args.ring_planar_rms:
        data['theta_deg'] = ''
        for mode_data in data['modes'].values():
          mode_data['phi_deg'] = ''

    row = dict(ref_atoms='-'.join(str(ra[i]) for i in ring),
      probe_atoms='-'.join(str(pb[i]) for i in ring), size=len(ring), status=status)
    for label,data in (('ref',pa),('probe',pb_data)):
      for key in ('plane_rms_A','plane_max_A','Q_A','theta_deg','status'):
        row[label+'_'+key] = data[key]
    row['delta_Q_A'] = pb_data['Q_A']-pa['Q_A'] if isinstance(pa['Q_A'],float) and isinstance(pb_data['Q_A'],float) else ''
    row['ref_analysis'] = ring_interpretation(pa,len(ring),args.ring_planar_rms)
    row['probe_analysis'] = ring_interpretation(pb_data,len(ring),args.ring_planar_rms)
    rows.append(row)
    print(f"{status:8s} REF {row['ref_atoms']} / PROBE {row['probe_atoms']} "
      f"size={len(ring)} Q(A)={brief_number(pa['Q_A'])} -> {brief_number(pb_data['Q_A'])}; "
      f"plane RMS(A)={pa['plane_rms_A']:.6f} -> {pb_data['plane_rms_A']:.6f}")
    print('  REF:', row['ref_analysis'])
    print('  PROBE:', row['probe_analysis'])
    if status == 'retained':
      delta_rms = pb_data['plane_rms_A']-pa['plane_rms_A']
      trend = 'more puckered' if delta_rms>0.01 else 'flatter' if delta_rms<-0.01 else 'similar puckering strength'
      print(f'  Comparison: PROBE {trend}; delta plane RMS={delta_rms:+.4f} A.')
    if status != 'retained':
      print('  WARNING: opposite-side metrics describe an open/missing cycle, not a confirmed ring.')
    for m in sorted(set(pa['modes'])|set(pb_data['modes'])):
      am,bm = pa['modes'].get(m,{}),pb_data['modes'].get(m,{})
      aq,bq = am.get('q_A',''),bm.get('q_A','')
      ap,bp = am.get('phi_deg',''),bm.get('phi_deg','')
      phase_delta = (bp-ap+180)%360-180 if isinstance(ap,float) and isinstance(bp,float) else ''
      mode = dict(ref_atoms=row['ref_atoms'],probe_atoms=row['probe_atoms'],mode=m,
        ref_q_A=aq,probe_q_A=bq,ref_phi_deg=ap,probe_phi_deg=bp,
        delta_q_A=bq-aq if isinstance(aq,float) and isinstance(bq,float) else '',
        delta_phi_deg=phase_delta)
      modes.append(mode)
      print(f'  mode {m}: q(A)={brief_number(aq)} -> {brief_number(bq)}; phi(deg)={brief_number(ap)} -> {brief_number(bp)}; delta_phi={brief_number(phase_delta)}')
    if len(ring)==6:
      print(f"  theta(deg)={brief_number(pa['theta_deg'])} -> {brief_number(pb_data['theta_deg'])}")
  if not rows:
    print('No internal rings found.')
  print('==========================================================================\n')
  return rows,modes


def advanced_report(sa, ca, sb, cb, ra, pb, fa, fb, args):
  """Run RDKit stereocenter, coordination shape and interpreted puckering reports."""
  from rdkit import Chem
  chiral = chirality_report(sa,ca,sb,cb,ra,pb,fa,fb,args)
  metals = metal_report(sa,ca,sb,cb,ra,pb,args)
  rings,modes = ring_report(sa,ca,sb,cb,ra,pb,fa,fb,args)
  if args.csv_prefix:
    tables = [('_chirality.csv',chiral,['ref_atom','probe_atom','element',
      'ref_configuration','probe_configuration']),
      ('_metal_coordination.csv',metals,['structure','atom','element','CN','neighbors','distances_A',
        'mean_distance_A','distance_std_A','best_model','shape_rms','angular_rms_deg','model_scores','interpretation']),
      ('_ring_puckering.csv',rings,['ref_atoms','probe_atoms','size','status','ref_plane_rms_A',
        'ref_plane_max_A','ref_Q_A','ref_theta_deg','ref_status','probe_plane_rms_A','probe_plane_max_A',
        'probe_Q_A','probe_theta_deg','probe_status','delta_Q_A','ref_analysis','probe_analysis']),
      ('_ring_modes.csv',modes,['ref_atoms','probe_atoms','mode','ref_q_A','probe_q_A','ref_phi_deg',
        'probe_phi_deg','delta_q_A','delta_phi_deg'])]
    for suffix,rows,fields in tables:
      write_csv(args.csv_prefix+suffix,rows,fields)


#------------------------------------------------------------------------------
def analyze_ddm(a, b, tolerance=0.01):
  """Analyze signed distance changes without alignment or Cartesian projection.

  For a symmetric matrix, U=Q, S=abs(lambda), V=Q*sign(lambda) is a valid SVD.
  Squared singular-value fractions describe matrix change, not PCA variance.
  Atom loadings describe distance-pattern participation, not atom displacement.
  """
  a,b = np.asarray(a,dtype=float),np.asarray(b,dtype=float)
  n = len(a)
  da = np.linalg.norm(a[:,None,:]-a[None,:,:],axis=2)
  db = np.linalg.norm(b[:,None,:]-b[None,:,:],axis=2)
  delta = (db-da)
  delta = (delta+delta.T)/2
  np.fill_diagonal(delta,0.)
  values,vectors = np.linalg.eigh(delta)
  order = np.argsort(-np.abs(values),kind='stable')
  values,vectors = values[order],vectors[:,order]
  # Choose a reproducible sign for display; a sign flip does not change a mode.
  for k in range(n):
    pivot = int(np.argmax(np.abs(vectors[:,k])))
    if vectors[pivot,k]<0:
      vectors[:,k] *= -1
  singular = np.abs(values)
  total = float(np.sum(singular**2))
  fraction = singular**2/total if total>1e-20 else np.zeros(n)
  cumulative = np.cumsum(fraction)
  ii,jj = np.triu_indices(n,1)
  differences = delta[ii,jj]
  increased = differences>tolerance
  decreased = differences<-tolerance
  unchanged = ~(increased|decreased)
  rank = int(np.sum(singular>max(singular[0]*1e-10,1e-12))) if n else 0
  effective_rank = float(np.exp(-np.sum(fraction[fraction>0]*np.log(fraction[fraction>0])))) if rank else 0.
  k90 = int(np.searchsorted(cumulative,.90)+1) if rank else 0
  k95 = int(np.searchsorted(cumulative,.95)+1) if rank else 0
  k99 = int(np.searchsorted(cumulative,.99)+1) if rank else 0
  row_rms = np.sqrt(np.sum(delta**2,axis=1)/(n-1)) if n>1 else np.zeros(n)
  row_energy = np.sum(delta**2,axis=1)/total if total>1e-20 else np.zeros(n)
  row_mean = np.sum(delta,axis=1)/(n-1) if n>1 else np.zeros(n)
  return dict(matrix=delta,ref_distances=da,probe_distances=db,eigenvalues=values,
    singular_values=singular,U=vectors,Vt=(vectors*np.sign(values)).T,
    fractions=fraction,cumulative=cumulative,rank=rank,effective_rank=effective_rank,
    k90=k90,k95=k95,k99=k99,pair_i=ii,pair_j=jj,pair_changes=differences,
    increased=increased,decreased=decreased,unchanged=unchanged,
    row_rms=row_rms,row_energy=row_energy,row_mean=row_mean)


def ddm_report(a,b,symbols,ref_ids,probe_ids,args):
  """Print concise matrix-pattern interpretation, never a motion reconstruction."""
  print('\n=============================DISPLACEMENT=================================')
  rms,maximum = distance_metrics(a,b)
  print(f'Distance RMSD: {rms:.6f} Angstroms')
  print(f'Distance MAXD: {maximum:.6f} Angstroms')
  if args.skip_ddm or len(a)<2:
    print('DDM analysis skipped.' if args.skip_ddm else 'DDM requires at least two atoms.')
    print('==========================================================================\n')
    return None
  if len(a)>args.ddm_max_atoms:
    print(f'WARNING: selected atoms exceed --ddm-max-atoms={args.ddm_max_atoms}; dense decomposition skipped.')
    print('==========================================================================\n')
    return None
  data=analyze_ddm(a,b,args.ddm_tolerance)
  n=len(a)
  print('DDM = D(PROBE)-D(REF), Angstrom. SVD obtained from symmetric eigendecomposition.')
  print('Explained = sigma^2/sum(sigma^2): matrix-change fraction, NOT statistical PCA variance or physical energy.')
  print('No alignment, Cartesian mapping, pseudo-vibrations or motion-path interpretation.')
  if data['rank']==0:
    print('No resolvable internal-distance change; rigid motion and mirror inversion are invisible to DDM.')
    print('==========================================================================\n')
    return data
  p=data['pair_changes']
  expanded=int(data['increased'].sum());contracted=int(data['decreased'].sum());stable=int(data['unchanged'].sum())
  print(f'Pair changes (|delta| <= {args.ddm_tolerance:g} A treated as unchanged): '
    f'increased={expanded}, decreased={contracted}, unchanged={stable} / {len(p)}')
  positive_energy=float(np.sum(p[p>0]**2));negative_energy=float(np.sum(p[p<0]**2))
  share=positive_energy/(positive_energy+negative_energy)
  print(f'Squared-change share: increasing distances {100*share:.1f}%; decreasing distances {100*(1-share):.1f}%.')
  trend='mainly distance increases' if share>.65 else 'mainly distance decreases' if share<.35 else 'mixed increases and decreases'
  print(f'Overall: {trend}; this does not imply uniform molecular expansion/contraction.')
  print(f"SVD rank={data['rank']}; effective rank={data['effective_rank']:.2f}; "
    f"modes for 90/95/99%={data['k90']}/{data['k95']}/{data['k99']}")
  concentration='strongly concentrated' if data['k90']<=3 else 'concentrated in a few patterns' if data['k90']<=6 else 'distributed across multiple patterns'
  print('Spectrum:',concentration,'(reporting heuristic; patterns are not independent motions).')
  pair_rows=[]
  for k in np.argsort(-np.abs(p),kind='stable'):
    i,j=int(data['pair_i'][k]),int(data['pair_j'][k])
    pair_rows.append(dict(ref_atoms=f'{ref_ids[i]}-{ref_ids[j]}',probe_atoms=f'{probe_ids[i]}-{probe_ids[j]}',
      elements=f'{symbols[i]}-{symbols[j]}',ref_distance_A=float(data['ref_distances'][i,j]),
      probe_distance_A=float(data['probe_distances'][i,j]),delta_A=float(p[k]),abs_delta_A=float(abs(p[k]))))
  print('Largest actual pair-distance changes (not necessarily bonded pairs):')
  for row in pair_rows[:args.ddm_top_pairs]:
    print(f"  REF {row['ref_atoms']} / PROBE {row['probe_atoms']} {row['elements']}: "
      f"{row['ref_distance_A']:.4f} -> {row['probe_distance_A']:.4f} A ({row['delta_A']:+.4f})")
  print('Atoms with most changed distance environments (not greatest Cartesian displacement):')
  atom_rows=[]
  for i in np.argsort(-data['row_rms'],kind='stable'):
    atom_rows.append(dict(ref_atom=ref_ids[i],probe_atom=probe_ids[i],element=symbols[i],
      distance_environment_rms_A=float(data['row_rms'][i]),signed_mean_change_A=float(data['row_mean'][i]),
      row_energy_fraction=float(data['row_energy'][i])))
  for row in atom_rows[:args.ddm_top_atoms]:
    print(f"  {row['element']} REF({row['ref_atom']})/PROBE({row['probe_atom']}): "
      f"row RMS={row['distance_environment_rms_A']:.4f} A; row share={100*row['row_energy_fraction']:.2f}%")
  print('Mode  Sigma(A)  Signed lambda(A)  Explained(%)  Cumulative(%)  Participation atoms')
  spectrum,loadings=[],[]
  for k in range(n):
    u=data['U'][:,k]
    participation=float(1/np.sum(u**4))
    row=dict(mode=k+1,singular_value_A=float(data['singular_values'][k]),
      signed_eigenvalue_A=float(data['eigenvalues'][k]),explained_fraction=float(data['fractions'][k]),
      cumulative_fraction=float(data['cumulative'][k]),participation_atoms=participation)
    spectrum.append(row)
    if k>=min(args.ddm_modes,data['rank']):
      continue
    print(f"{k+1:4d} {row['singular_value_A']:9.5f} {row['signed_eigenvalue_A']:17.5f} "
      f"{100*row['explained_fraction']:13.3f} {100*row['cumulative_fraction']:14.3f} {participation:19.2f}")
    fraction=participation/n
    locality='localized participation' if fraction<=.25 else 'broad participation' if fraction>=.60 else 'intermediate participation'
    top=np.argsort(-(u*u),kind='stable')[:args.ddm_top_atoms]
    names=', '.join(f'{symbols[i]}({ref_ids[i]})[{100*u[i]**2:.1f}%]' for i in top)
    print(f'  {locality}; leading loadings: {names}')
    # Qualify sign groups as a rank-one distance component, never a net distance verdict.
    cutoff=.25*np.max(np.abs(u))
    pos=sorted(np.flatnonzero(u>=cutoff),key=lambda i:-abs(u[i]))[:3]
    neg=sorted(np.flatnonzero(u<=-cutoff),key=lambda i:-abs(u[i]))[:3]
    fmt=lambda group:','.join(f'{symbols[i]}({ref_ids[i]})' for i in group) or 'none'
    within='increase' if data['eigenvalues'][k]>0 else 'decrease'
    between='decrease' if data['eigenvalues'][k]>0 else 'increase'
    if pos and neg:
      print(f'  Component groups +[{fmt(pos)}], -[{fmt(neg)}]: same-sign pairs {within}; opposite-sign pairs {between}.')
    else:
      print(f'  Leading same-sign group [{fmt(pos or neg)}]: component pair distances {within}.')
    print('  Group signs describe this matrix component only; other components can cancel it.')
    neighbors=[j for j in range(n) if j!=k and abs(data['eigenvalues'][j]-data['eigenvalues'][k])
      <=1e-3*max(abs(data['eigenvalues'][k]),1e-12)]
    if neighbors:
      print('  Near-degenerate eigenvalue: individual atom loadings may rotate within the shared subspace.')
    for i in range(n):
      loadings.append(dict(mode=k+1,ref_atom=ref_ids[i],probe_atom=probe_ids[i],element=symbols[i],
        signed_loading=float(u[i]),participation_fraction=float(u[i]**2)))
  selected=min(args.ddm_modes,data['rank'])
  coverage=float(data['cumulative'][selected-1])
  print(f'Top {selected} patterns cover {100*coverage:.3f}% of squared matrix change; '
    f'relative full-matrix reconstruction error={100*np.sqrt(max(0.,1-coverage)):.3f}%.')
  print('Participation counts/loadings describe distance patterns, not moving atom groups or bond strengths.')
  if args.csv_prefix:
    prefix=args.csv_prefix+'_ddm'
    write_csv(prefix+'_spectrum.csv',spectrum,['mode','singular_value_A','signed_eigenvalue_A',
      'explained_fraction','cumulative_fraction','participation_atoms'])
    write_csv(prefix+'_atoms.csv',atom_rows,['ref_atom','probe_atom','element',
      'distance_environment_rms_A','signed_mean_change_A','row_energy_fraction'])
    write_csv(prefix+'_pairs.csv',pair_rows,['ref_atoms','probe_atoms','elements',
      'ref_distance_A','probe_distance_A','delta_A','abs_delta_A'])
    write_csv(prefix+'_loadings.csv',loadings,['mode','ref_atom','probe_atom','element',
      'signed_loading','participation_fraction'])
    if args.ddm_save_matrix:
      np.savez_compressed(prefix+'_matrix.npz',DDM=data['matrix'],U=data['U'],
        singular_values=data['singular_values'],Vt=data['Vt'],signed_eigenvalues=data['eigenvalues'],
        ref_atoms=np.array(ref_ids),probe_atoms=np.array(probe_ids),elements=np.array(symbols))
    print('DDM tables saved with prefix:',prefix)
  print('==========================================================================\n')
  return data


#------------------------------------------------------------------------------
def molecular_shape_metrics(coords, masses):
  """Calculate mass-weighted gyration and inertia about the center of mass."""
  coords,masses = np.asarray(coords,dtype=float),np.asarray(masses,dtype=float)
  if coords.shape!=(len(masses),3) or not np.all(np.isfinite(masses)) or np.any(masses<=0):
    raise ValueError('Shape analysis requires positive finite atomic masses.')
  total=float(masses.sum())
  com=np.sum(coords*masses[:,None],axis=0)/total
  x=coords-com
  gyration=(x*masses[:,None]).T@x/total
  inertia=total*(np.trace(gyration)*np.eye(3)-gyration)
  moments,axes=np.linalg.eigh(inertia)
  moments=np.maximum(moments,0.)
  geig=np.maximum(np.linalg.eigvalsh(gyration),0.)
  rg=float(np.sqrt(geig.sum()))
  geom=coords-coords.mean(axis=0)
  rg_geom=float(np.sqrt(np.mean(np.sum(geom*geom,axis=1))))
  anisotropy=float(1.5*np.sum(geig**2)/geig.sum()**2-.5) if geig.sum()>1e-20 else 0.
  anisotropy=float(np.clip(anisotropy,0.,1.))
  for k in range(3):
    pivot=int(np.argmax(np.abs(axes[:,k])))
    if axes[pivot,k]<0:
      axes[:,k]*=-1
  return dict(total_mass_u=total,com=com,centered=x,gyration=gyration,
    inertia=inertia,moments=moments,axes=axes,gyration_eigenvalues=geig,
    Rg_A=rg,Rg_geometric_A=rg_geom,kappa2=anisotropy)


def whole_align(a,b,masses):
  """Align corresponding atoms by proper mass-weighted Kabsch rotation."""
  mass=np.asarray(masses)
  ac=a-np.sum(a*mass[:,None],axis=0)/mass.sum()
  bc=b-np.sum(b*mass[:,None],axis=0)/mass.sum()
  u,_,vt=np.linalg.svd((bc*mass[:,None]).T@ac)
  fix=np.eye(3)
  fix[-1,-1]=1. if np.linalg.det(u@vt)>=0 else -1.
  rotation=u@fix@vt
  fitted=bc@rotation
  rms=float(np.sqrt(np.sum(mass*np.sum((fitted-ac)**2,axis=1))/mass.sum()))
  return ac,fitted,rotation,rms


def topology_comparison(symbols,ref_bonds,probe_bonds):
  """Compare element-labeled connectivity, not bond orders or stereochemistry."""
  import networkx as nx
  graphs=[]
  for edges in (ref_bonds,probe_bonds):
    graph=nx.Graph()
    graph.add_nodes_from((i,dict(element=symbol)) for i,symbol in enumerate(symbols))
    graph.add_edges_from(edges)
    graphs.append(graph)
  ga,gb=graphs
  lost,gained=ref_bonds-probe_bonds,probe_bonds-ref_bonds
  union=ref_bonds|probe_bonds
  mapped_jaccard=len(ref_bonds&probe_bonds)/len(union) if union else 1.
  isomorphic=not lost and not gained or nx.is_isomorphic(ga,gb,
    node_match=nx.algorithms.isomorphism.categorical_node_match('element',''))
  ca,cb=nx.number_connected_components(ga),nx.number_connected_components(gb)
  degrees=[i for i in range(len(symbols)) if ga.degree[i]!=gb.degree[i]]
  return dict(ref_edges=len(ref_bonds),probe_edges=len(probe_bonds),
    lost_edges=lost,gained_edges=gained,mapped_edge_jaccard=float(mapped_jaccard),
    element_graph_isomorphic=bool(isomorphic),ref_components=ca,probe_components=cb,
    ref_cycle_rank=len(ref_bonds)-len(symbols)+ca,
    probe_cycle_rank=len(probe_bonds)-len(symbols)+cb,degree_changed_atoms=degrees)


def whole_shape_tanimoto(symbols,ref_coords,probe_coords,ref_bonds,probe_bonds,spacing):
  """Compare RDKit grids after atom-mapped alignment, including explicit H.

  No additional shape-based optimization is performed: this is not the maximum
  possible overlap. The similarity is 1 minus the RDKit shape Tanimoto distance.
  """
  from rdkit.Chem import rdShapeHelpers
  ref=rdkit_graph(symbols,ref_coords,ref_bonds)
  probe=rdkit_graph(symbols,probe_coords,probe_bonds)
  return float(np.clip(1.-rdShapeHelpers.ShapeTanimotoDist(ref,probe,
    gridSpacing=spacing,ignoreHs=False),0.,1.))


#------------------------------------------------------------------------------
def whole_similarity_text(metric, value):
  """Use explicit reporting bands; none are universal chemical cutoffs."""
  if value is None or not np.isfinite(value):
    return 'unavailable'
  if metric == 'mapped_edges':
    if value == 1.:
      return 'identical mapped edge sets'
    if value >= .95:
      return 'high mapped edge overlap; a small fraction of edges differs'
    if value >= .80:
      return 'partial mapped edge differences; most edges are shared'
    return 'substantial mapped edge differences'
  if metric == 'shape':
    if value >= .90:
      return 'very high spatial overlap'
    if value >= .75:
      return 'high spatial overlap with visible contour differences'
    if value >= .50:
      return 'moderate spatial overlap; non-negligible contour differences'
    return 'low spatial overlap; pronounced contour differences at this alignment'
  if metric == 'gyration':
    if value >= .98:
      return 'nearly unchanged normalized second-moment shape'
    if value >= .95:
      return 'closely similar normalized second-moment shape'
    if value >= .90:
      return 'moderate redistribution of normalized second moments'
    return 'substantial redistribution of normalized second moments'
  raise ValueError('Unknown whole-molecule similarity metric.')


def whole_anisotropy_analysis(ref, probe):
  """Interpret dimensionless differences without confusing inertia with length."""
  delta=probe['kappa2']-ref['kappa2']
  magnitude=abs(delta)
  extent='minor' if magnitude<=.01 else 'mild' if magnitude<=.05 else 'clear'
  trend='approximately unchanged directional anisotropy' if magnitude<=.01 else (
    'more directionally uneven mass distribution' if delta>0 else 'more balanced directional mass distribution')
  print(f'  Anisotropy analysis: delta kappa2={delta:+.4f}; {extent} change, {trend}.')
  if min(ref['moments'][-1],probe['moments'][-1])<=1e-20:
    print('  Normalized-PMI interpretation unavailable: zero spatial extent.')
    return
  a=ref['moments'][:2]/ref['moments'][-1]
  b=probe['moments'][:2]/probe['moments'][-1]
  change=b-a
  largest=float(np.max(np.abs(change)))
  extent='minor' if largest<=.02 else 'mild' if largest<=.05 else 'clear'
  print(f'  PMI analysis: delta(I1/I3)={change[0]:+.4f}, delta(I2/I3)={change[1]:+.4f}; '
    f'{extent} normalized-inertia redistribution.')
  if magnitude<=.01 and largest<=.02:
    print('  Overall shape: coarse mass-distribution shape is essentially retained.')
  elif magnitude<=.01 and largest<=.05:
    print('  Overall shape: broad anisotropy retained, with mild redistribution among principal directions; no clear shape-class change.')
  else:
    print('  Overall shape: directional mass distribution changes; no unique rod/disk/sphere assignment follows from these deltas.')
  # Recover ordered gyration fractions from inertia, avoiding simple axis-length claims.
  shares=[]
  for data in (ref,probe):
    values=np.maximum(.5*data['moments'].sum()-data['moments'],0.)
    shares.append(values/values.sum())
  for index,label in ((0,'longest-spread direction'),(2,'weakest-spread direction')):
    diff=float(shares[1][index]-shares[0][index])
    trend='similar' if abs(diff)<=.01 else 'increased' if diff>0 else 'decreased'
    print(f'  Relative squared-spread share along {label}: '
      f'{shares[0][index]:.4f} -> {shares[1][index]:.4f} ({trend}; 0.01 reporting tolerance).')
  print('  These are normalized mass-spread shares, not axis lengths, thickness changes in Angstrom, or principal-axis rotations.')
  print('  Reporting cutoffs: |delta kappa2|=0.01/0.05; max |delta normalized PMI|=0.02/0.05 (heuristic).')


def whole_molecule_report(symbols,a,b,ref_ids,probe_ids,ref_bonds,probe_bonds,args):
  """Summarize full-molecule topology and shape using transparent definitions."""
  from rdkit import Chem
  table=Chem.GetPeriodicTable()
  masses=np.array([table.GetAtomicWeight(symbol) for symbol in symbols],dtype=float)
  ac,bc,rotation,aligned_rms=whole_align(a,b,masses)
  ref,probe=molecular_shape_metrics(ac,masses),molecular_shape_metrics(bc,masses)
  topology=topology_comparison(symbols,ref_bonds,probe_bonds)
  print('\n============================WHOLE_MOLECULE================================')
  print('Full-molecule selection: both C/c or --all; atom correspondence is preserved.')
  print('TOPOLOGY: element-labeled distance-inferred graph; metal contacts included; no bond orders/stereo.')
  print(f"Edges: {topology['ref_edges']} -> {topology['probe_edges']}; "
    f"lost={len(topology['lost_edges'])}, gained={len(topology['gained_edges'])}; "
    f"mapped edge Jaccard={topology['mapped_edge_jaccard']:.4f}")
  print(f"Connected components: {topology['ref_components']} -> {topology['probe_components']}; "
    f"cycle rank E-N+C: {topology['ref_cycle_rank']} -> {topology['probe_cycle_rank']} (not total ring count).")
  print('Element-labeled graph isomorphism (renumbering allowed):',topology['element_graph_isomorphic'])
  print('Mapped edge Jaccard analysis:',whole_similarity_text('mapped_edges',topology['mapped_edge_jaccard']))
  print('  Bands: 1 identical; >=0.95 high overlap; >=0.80 mostly shared; <0.80 substantial differences (heuristic).')
  print('  A high score does not make changed edges chemically unimportant; inspect lost/gained edges and distance cutoffs.')

  for label,edges in (('lost',topology['lost_edges']),('gained',topology['gained_edges'])):
    if edges:
      print('  '+label+': '+', '.join(f"REF {ref_ids[i]}-{ref_ids[j]} / PROBE {probe_ids[i]}-{probe_ids[j]} "
        f"({symbols[i]}-{symbols[j]})" for i,j in sorted(edges)))
  changed=topology['degree_changed_atoms']
  if changed:
    print('Degree-changed corresponding atoms:',','.join(str(ref_ids[i]) for i in changed))
  if not topology['lost_edges'] and not topology['gained_edges']:
    print('Topology conclusion: identical connectivity under the supplied atom mapping.')
  elif topology['element_graph_isomorphic']:
    print('Topology conclusion: mapped edges differ, but the element graph is equivalent under another permutation.')
  else:
    print('Topology conclusion: different inferred element/connectivity graphs; check borderline distances before claiming bond changes.')
  if topology['ref_components']>1 or topology['probe_components']>1:
    print('WARNING: disconnected structures; global shape reflects inter-fragment separation as well as fragment shape.')
  print('\nSHAPE: all explicit atoms; RDKit standard atomic weights; COM-centered, mass-weighted quantities.')
  print('Shape overlap/axis angles use proper mass-weighted atom-mapped alignment, without reflection or scale fitting.')
  print('Metric                         REF                PROBE              Delta')
  def value_line(name,x,y):
    print(f'{name:26s} {x:17.6f} {y:20.6f} {y-x:+17.6f}')
  value_line('Rg mass-weighted (A)',ref['Rg_A'],probe['Rg_A'])
  value_line('Rg geometric (A)',ref['Rg_geometric_A'],probe['Rg_geometric_A'])
  for k in range(3):
    value_line(f'Principal moment I{k+1} (u A^2)',ref['moments'][k],probe['moments'][k])
  value_line('Relative anisotropy kappa2',ref['kappa2'],probe['kappa2'])
  if ref['moments'][-1]>1e-20 and probe['moments'][-1]>1e-20:
    value_line('Normalized PMI I1/I3',ref['moments'][0]/ref['moments'][2],probe['moments'][0]/probe['moments'][2])
    value_line('Normalized PMI I2/I3',ref['moments'][1]/ref['moments'][2],probe['moments'][1]/probe['moments'][2])
  print('Rg: mass spread/compactness. kappa2: 0 isotropic second moments, 1 line-like; neither implies bonding or energetics.')
  whole_anisotropy_analysis(ref,probe)
  print('PRINCIPAL AXES (both in COM-centered REF-aligned coordinates; signs are arbitrary):')
  axes=[]
  for k in range(3):
    stable=[]
    for metrics in (ref,probe):
      gaps=[abs(metrics['moments'][k]-metrics['moments'][j]) for j in range(3) if j!=k]
      stable.append(min(gaps)>args.whole_axis_tol*max(metrics['moments'][-1],1e-15))
    angle=float(np.degrees(np.arccos(np.clip(abs(ref['axes'][:,k]@probe['axes'][:,k]),0.,1.)))) if all(stable) else None
    rvec=','.join(f'{v:+.4f}' for v in ref['axes'][:,k])
    pvec=','.join(f'{v:+.4f}' for v in probe['axes'][:,k])
    text=f'{angle:.3f} deg' if angle is not None else 'undefined/unstable: near-degenerate moments'
    print(f'  Axis {k+1}: REF({rvec}) PROBE({pvec}); unsigned angle={text}')
    axes.append(dict(axis=k+1,ref_x=float(ref['axes'][0,k]),ref_y=float(ref['axes'][1,k]),ref_z=float(ref['axes'][2,k]),
      probe_x=float(probe['axes'][0,k]),probe_y=float(probe['axes'][1,k]),probe_z=float(probe['axes'][2,k]),
      unsigned_angle_deg=angle if angle is not None else '',stable=all(stable)))
  shape_similarity=None
  try:
    shape_similarity=whole_shape_tanimoto(symbols,ac,bc,ref_bonds,probe_bonds,args.whole_shape_grid)
    print(f'Shape Tanimoto similarity: {shape_similarity:.4f} (1 identical overlap; grid={args.whole_shape_grid:g} A; H included).')
    print('This is overlap at the mapped alignment, NOT a globally maximized shape score.')
    print('  Shape Tanimoto analysis:',whole_similarity_text('shape',shape_similarity))
    print('  Bands: >=0.90 very high; >=0.75 high; >=0.50 moderate; <0.50 low overlap (heuristic).')
    print('  Grid spacing, hydrogen inclusion and alignment affect this score; it is not a percentage of unchanged atoms.')

  except (RuntimeError,ValueError) as exc:
    print('WARNING: RDKit shape overlap unavailable:',exc)
  g1,g2=ref['gyration_eigenvalues'],probe['gyration_eigenvalues']
  descriptor_similarity=float(1.-.5*np.sum(np.abs(g1/g1.sum()-g2/g2.sum()))) if min(g1.sum(),g2.sum())>1e-20 else None
  if descriptor_similarity is not None:
    print(f'Normalized gyration-spectrum similarity: {descriptor_similarity:.4f} (size-independent second-moment descriptor, not full shape).')
    print('  Gyration-spectrum analysis:',whole_similarity_text('gyration',descriptor_similarity))
    print('  Bands: >=0.98 nearly unchanged; >=0.95 close; >=0.90 moderate redistribution; <0.90 substantial (heuristic).')
    print('  High similarity constrains coarse directional proportions only, not local contours, chirality or internal atom arrangement.')
    if descriptor_similarity>=.95 and shape_similarity is not None and shape_similarity<.75:
      print('  Combined interpretation: coarse mass-distribution proportions are retained, but spatial contours differ at the mapped alignment.')

  print(f'Atom-mapped mass-weighted aligned RMSD: {aligned_rms:.6f} A (no mirror fit).')
  rms,maximum=distance_metrics(a,b)
  mean_rg=(ref['Rg_A']+probe['Rg_A'])/2.
  geometric_mean=(ref['Rg_geometric_A']+probe['Rg_geometric_A'])/2.
  ratio=rms/mean_rg if mean_rg>1e-12 else None
  geometric_ratio=rms/geometric_mean if geometric_mean>1e-12 else None
  thresholds=args.whole_thresholds
  levels=('negligible','small','moderate','large')
  level=levels[int(np.searchsorted(thresholds,ratio,side='left'))] if ratio is not None else 'undefined'
  rg_relative=(probe['Rg_A']-ref['Rg_A'])/ref['Rg_A'] if ref['Rg_A']>1e-12 else None
  print('\nGLOBAL DIFFERENCE SUMMARY')
  if ratio is not None:
    print(f'dRMSD / mean mass-weighted Rg = {rms:.6f}/{mean_rg:.6f} = {ratio:.5f}; classification: {level}.')
    print(f'dRMSD / mean geometric Rg = {geometric_ratio:.5f} (equal-atom companion normalization).')
  else:
    print('dRMSD/Rg undefined: molecular spatial extent is zero.')
  print('Heuristic ratio bands: <= '+', <= '.join(f'{x:g}' for x in thresholds)+
    ' correspond to negligible/small/moderate; above the last cutoff is large.')
  if rg_relative is not None:
    size='more extended' if rg_relative>.02 else 'more compact' if rg_relative<-.02 else 'similar overall mass spread'
    print(f'PROBE Rg change={100*rg_relative:+.2f}%: {size} (2% reporting threshold).')
  if ratio is not None and ratio<=thresholds[0] and aligned_rms/max(mean_rg,1e-12)>.10:
    print('WARNING: near-equal distances but poor proper alignment; inspect mirror inversion or atom mapping.')
  print('Conclusion refers to global internal-distance differences; small global change can coexist with large local changes.')
  print('No universal chemistry thresholds; topology, scalar shape and spatial overlap answer different questions.')
  if args.csv_prefix:
    records=[]
    for label,metrics in (('REF',ref),('PROBE',probe)):
      records.append(dict(structure=label,total_mass_u=metrics['total_mass_u'],Rg_A=metrics['Rg_A'],
        Rg_geometric_A=metrics['Rg_geometric_A'],I1_uA2=float(metrics['moments'][0]),
        I2_uA2=float(metrics['moments'][1]),I3_uA2=float(metrics['moments'][2]),kappa2=metrics['kappa2']))
    write_csv(args.csv_prefix+'_whole_shape.csv',records,['structure','total_mass_u','Rg_A','Rg_geometric_A',
      'I1_uA2','I2_uA2','I3_uA2','kappa2'])
    write_csv(args.csv_prefix+'_whole_axes.csv',axes,['axis','ref_x','ref_y','ref_z','probe_x','probe_y','probe_z','unsigned_angle_deg','stable'])
    summary=dict(mapped_edge_jaccard=topology['mapped_edge_jaccard'],element_graph_isomorphic=topology['element_graph_isomorphic'],
      ref_edges=topology['ref_edges'],probe_edges=topology['probe_edges'],lost_edges=len(topology['lost_edges']),
      gained_edges=len(topology['gained_edges']),ref_components=topology['ref_components'],probe_components=topology['probe_components'],
      ref_cycle_rank=topology['ref_cycle_rank'],probe_cycle_rank=topology['probe_cycle_rank'],
      shape_tanimoto=shape_similarity if shape_similarity is not None else '',
      gyration_spectrum_similarity=descriptor_similarity if descriptor_similarity is not None else '',
      aligned_mass_RMSD_A=aligned_rms,dRMSD_A=rms,dMAXD_A=maximum,mean_mass_Rg_A=mean_rg,
      normalized_dRMSD=ratio if ratio is not None else '',geometric_normalized_dRMSD=geometric_ratio if geometric_ratio is not None else '',
      difference_class=level)
    write_csv(args.csv_prefix+'_whole_summary.csv',[summary],list(summary))
  print('==========================================================================\n')
  return dict(ref=ref,probe=probe,topology=topology,axes=axes,shape_tanimoto=shape_similarity,
    normalized_dRMSD=ratio,difference_class=level,alignment_rotation=rotation)


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
  parser.add_argument('--hb-da', type=float, default=3.0,
    help='Maximum hydrogen-bond D-A distance in Angstrom (default: 3.0)')
  parser.add_argument('--hb-angle', type=float, default=150.0,
    help='Minimum D-H-A angle in degrees (default: 150)')
  parser.add_argument('--hb-donors', default='N,O,S', help='Comma-separated HB donor elements')
  parser.add_argument('--acceptors', default='N,O,S', help='Candidate acceptor elements: N,O,S')
  parser.add_argument('--xb-scale', type=float, default=1.0,
    help='Halogen-contact cutoff as a multiple of the VDW radius sum')
  parser.add_argument('--xb-angle', type=float, default=150.0,
    help='Minimum R-X-A angle in degrees (default: 150)')
  parser.add_argument('--ring-planar-rms', type=float, default=0.02,
    help='Best-fit-plane RMS threshold in Angstrom for approximately planar rings')
  parser.add_argument('--coord-scale', type=float, default=1.15,
    help='Metal neighbor distance cutoff as a multiple of covalent-radius sum')
  parser.add_argument('--ring-mode', choices=('all','symmsssr'), default='all',
    help='All simple cycles or RDKit SymmSSSR (default: all)')
  parser.add_argument('--ring-max-size', type=int, default=0,
    help='Maximum ring size; 0 means no size restriction')
  parser.add_argument('--ring-max-count', type=int, default=10000,
    help='Ring count safeguard; truncation is explicitly reported')
  parser.add_argument('--include-metal-rings', action='store_true',
    help='Include coordination cycles containing metal atoms in ring analysis')
  parser.add_argument('--skip-advanced', action='store_true',
    help='Skip CORE analysis; whole-molecule shape still requires RDKit unless --skip-whole is used')
  parser.add_argument('--ddm-modes',type=int,default=5,help='Leading DDM patterns to print (default: 5)')
  parser.add_argument('--ddm-top-atoms',type=int,default=5,help='Number of leading atoms to print')
  parser.add_argument('--ddm-top-pairs',type=int,default=5,help='Number of actual changed atom pairs to print')
  parser.add_argument('--ddm-tolerance',type=float,default=0.01,help='Near-zero pair change tolerance in Angstrom')
  parser.add_argument('--ddm-max-atoms',type=int,default=500,help='Dense SVD resource safeguard')
  parser.add_argument('--skip-ddm',action='store_true',help='Print distance metrics without DDM decomposition')
  parser.add_argument('--ddm-save-matrix',action='store_true',help='Save SVD matrices to NPZ when --csv-prefix is supplied')
  parser.add_argument('--skip-whole',action='store_true',help='Disable the automatic whole-molecule report')
  parser.add_argument('--whole-axis-tol',type=float,default=0.02,
    help='Relative principal-moment gap threshold for unstable axes')
  parser.add_argument('--whole-shape-grid',type=float,default=0.25,
    help='RDKit shape-overlap grid spacing in Angstrom')
  parser.add_argument('--whole-thresholds',default='0.02,0.10,0.25',
    help='Three increasing dRMSD/mean-Rg reporting cutoffs')
  args = parser.parse_args(argv)
  try:
    args.whole_thresholds=tuple(float(x.strip()) for x in args.whole_thresholds.split(','))
  except ValueError:
    parser.error('--whole-thresholds requires three numerical cutoffs')
  if len(args.whole_thresholds)!=3 or not all(np.isfinite(x) and x>0 for x in args.whole_thresholds) or not all(
    x<y for x,y in zip(args.whole_thresholds,args.whole_thresholds[1:])):
    parser.error('--whole-thresholds requires three positive increasing cutoffs')
  if not np.isfinite(args.whole_axis_tol) or not 0<args.whole_axis_tol<1:
    parser.error('--whole-axis-tol must lie between 0 and 1')
  if not np.isfinite(args.whole_shape_grid) or args.whole_shape_grid<=0:
    parser.error('--whole-shape-grid must be positive and finite')

  if min(args.ddm_modes,args.ddm_top_atoms,args.ddm_top_pairs)<1 or args.ddm_max_atoms<2:
    parser.error('DDM mode/atom/pair counts must be positive; atom safeguard must be >=2')
  if not np.isfinite(args.ddm_tolerance) or args.ddm_tolerance<0:
    parser.error('--ddm-tolerance must be nonnegative and finite')
  if args.ddm_save_matrix and not args.csv_prefix:
    parser.error('--ddm-save-matrix requires --csv-prefix')

  if not np.isfinite(args.ring_planar_rms) or not 0 < args.ring_planar_rms < 0.10:
    parser.error('--ring-planar-rms must lie between 0 and 0.10 Angstrom')
  if not np.isfinite(args.coord_scale) or args.coord_scale <= 0:
    parser.error('--coord-scale must be positive and finite')
  if args.ring_max_count < 1 or args.ring_max_size < 0 or args.ring_max_size in (1,2):
    parser.error('--ring-max-count must be positive; --ring-max-size must be 0 or >=3')
  if not args.skip_advanced:
    try:
      from rdkit import Chem
    except ImportError as exc:
      raise ImportError('Advanced analysis requires RDKit; install it or use --skip-advanced.') from exc

  args.hb_donors = tuple(x.strip().capitalize() for x in args.hb_donors.split(','))
  args.acceptors = tuple(x.strip().capitalize() for x in args.acceptors.split(','))
  if not args.hb_donors or any(x not in ('N', 'O', 'S') for x in args.hb_donors):
    parser.error('--hb-donors supports N,O,S only')
  if not args.acceptors or any(x not in ('N', 'O', 'S') for x in args.acceptors):
    parser.error('--acceptors supports N,O,S only')
  if not np.isfinite(args.hb_da) or args.hb_da <= 0 or not np.isfinite(args.xb_scale) or args.xb_scale <= 0:
    parser.error('Distance cutoffs must be positive and finite')
  if not 0 <= args.hb_angle <= 180 or not 0 <= args.xb_angle <= 180:
    parser.error('Angle cutoffs must lie between 0 and 180 degrees')

  if args.top < 1 or not 0 < args.linear_tol < 1:
    parser.error('--top must be positive; --linear-tol must lie between 0 and 1')
  if args.all and (args.ref_atoms is not None or args.probe_atoms is not None):
    parser.error('--all cannot be combined with fragment selections')
  if (args.ref_atoms is None) != (args.probe_atoms is None):
    parser.error('Provide both --ref-atoms and --probe-atoms')
  sa, ca = read_xyz(args.ref_xyz)
  sb, cb = read_xyz(args.probe_xyz)
  whole_selected=bool(args.all)
  if args.all:
    ra, pb = list(range(1,len(sa)+1)), list(range(1,len(sb)+1))
  else:
    print(f'Reference XYZ (first argument): {args.ref_xyz}')
    print(f'Probe XYZ (second argument): {args.probe_xyz}')
    rt = args.ref_atoms if args.ref_atoms is not None else input('Enter REF atoms (1-based; e.g. 1-9,16-23,25): ')
    pt = args.probe_atoms if args.probe_atoms is not None else input('Enter PROBE atoms (1-based; same order as REF): ')
    whole_selected=rt.strip().lower()=='c' and pt.strip().lower()=='c'
    ra = _parse_indices(rt, 'ref_indices', len(sa))
    pb = _parse_indices(pt, 'probe_indices', len(sb))
  if whole_selected and not args.skip_whole:
    try:
      import networkx
      from rdkit import Chem
    except ImportError as exc:
      raise ImportError('Whole-molecule analysis requires RDKit and NetworkX; install them or use --skip-whole.') from exc
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
  # Use the full-molecule topology for donor/acceptor context, but report only selected atoms.
  ref_radii = np.array([float(qcel.covalentradii.get(x, units='angstrom')) for x in sa])
  probe_radii = np.array([float(qcel.covalentradii.get(x, units='angstrom')) for x in sb])
  for array in (ref_radii, probe_radii):
    if not np.all(np.isfinite(array)) or np.any(array <= 0):
      raise ValueError('Invalid full-molecule covalent radii')
  covalent_changes, contact_rows = interaction_report(a, b, symbols, ra, pb,
    bonds, probe_bonds, ca, cb, sa, sb, infer_bonds(ca, ref_radii),
    infer_bonds(cb, probe_radii), args)
  ddm_report(a,b,symbols,ra,pb,args)
  if not args.skip_advanced:
    advanced_report(sa, ca, sb, cb, ra, pb, infer_bonds(ca, ref_radii),
      infer_bonds(cb, probe_radii), args)
  ranked, excluded = compare_internal(a,b,symbols,ra,pb,bonds,probe_bonds,args.linear_tol)
  all_rows, top_rows = [], []
  print('\n============================INTERNAL_COORD================================')

  for kind, rows in ranked.items():
    if kind.upper() == 'BOND':
      print(f'{kind.upper()}: top {min(args.top,len(rows))} / {len(rows)} valid')
    else:
      print(f'\n{kind.upper()}: top {min(args.top,len(rows))} / {len(rows)} valid')
    print('Rank  REF atoms         PROBE atoms       Elements       REF value    PROBE value    Delta    |delta|   bonds     Type')
    for rank,row in enumerate(rows,1):
      record = dict(rank=rank, **row)
      all_rows.append(record)
      if rank <= args.top:
        top_rows.append(record)
        print(f"{rank:4d}  {row['ref_atoms']:17s} {row['probe_atoms']:17s} {row['elements']:12s} {row['ref_value']:11.6f} {row['probe_value']:13.6f} {row['delta']:10.6f} {row['abs_delta']:10.6f} {'  OK  ' if row['probe_bonds_ok'] else 'WARN  '} {row['dihedral_type'] or '     -'}")
  print(f'\nExcluded coordinates: {len(excluded)}')
  for row in excluded:
    print(f"  {row['kind']} REF {row['ref_atoms']} / PROBE {row['probe_atoms']}: {row['reason']} (type={row['dihedral_type'] or '-'})")
  print('==========================================================================\n')
  if whole_selected and not args.skip_whole:
    whole_molecule_report(symbols,a,b,ra,pb,bonds,probe_bonds,args)
  if args.csv_prefix:
    fields = ['rank','kind','dihedral_type','ref_atoms','probe_atoms','elements','unit','ref_value','probe_value','delta','abs_delta','probe_bonds_ok']
    write_csv(args.csv_prefix+'_covalent_changes.csv', covalent_changes,
      ['ref_atoms', 'probe_atoms', 'elements', 'ref_distance_A', 'probe_distance_A',
        'ref_bonded', 'probe_bonded', 'status'])
    write_csv(args.csv_prefix+'_contacts.csv', contact_rows,
      ['kind', 'status', 'ref_atoms', 'probe_atoms', 'elements', 'ref_present',
        'probe_present', 'ref_primary_distance_A', 'ref_secondary_distance_A',
        'ref_angle_deg', 'probe_primary_distance_A', 'probe_secondary_distance_A',
        'probe_angle_deg'])
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
