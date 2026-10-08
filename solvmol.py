#!/home/lky/miniconda3/envs/polygview/bin/python
"""Compare fixed-geometry solvation responses of two whole-molecule conformers.
Author: Dirac4pi
Dependencies: numpy, rdkit, freesasa; external executable: xtb >= 6.3.3.
No geometry optimization, frequency calculation, or conformer ensemble is used.
"""
import argparse
from collections import Counter
import csv
import json
import os
from pathlib import Path
import re
import shutil
import subprocess
import sys
import tempfile
import numpy as np

HARTREE_TO_KCAL = 627.5094740631
AU_TO_DEBYE = 2.541746473
SOLVENTS = ('acetone acetonitrile aniline benzaldehyde benzene ch2cl2 chcl3 cs2 '
  'dioxane dmf dmso ethanol ether ethylacetate furane hexadecane hexane methanol '
  'nitromethane octanol woctanol phenol toluene thf water').split()
NUMBER = r'[+-]?(?:\d+\.?\d*|\.\d+)(?:[EeDd][+-]?\d+)?'

#------------------------------------------------------------------------------
# Source: Jaguar 7.6 User Manual (2009), Table 3.3, PBF solvent parameters.
# These imported probe conventions affect FreeSASA geometry only, not xTB ALPB.
PROBE_RADIUS_SOURCE = 'Jaguar 7.6 User Manual (2009), Table 3.3: PBF probe radii'
SOLVENT_PROBE_RADII = {
  'water': 1.40,
  'methanol': 2.00,
  'ethanol': 2.27,
  'ch2cl2': 2.33,
  'dmso': 2.41,
  'dmf': 2.49,
  '1,2-dichloroethane': 2.51,
  'chcl3': 2.52,
  'thf': 2.52,
  'benzene': 2.60,
  'chlorobenzene': 2.72,
  'nitrobenzene': 2.73,
  'cyclohexane': 2.78,
}
SOLVENT_ALIASES = {
  'h2o': 'water',
  'ch3oh': 'methanol',
  'c2h5oh': 'ethanol',
  'dichloromethane': 'ch2cl2',
  'methylene chloride': 'ch2cl2',
  'chloroform': 'chcl3',
  'dimethyl sulfoxide': 'dmso',
  'dimethylsulfoxide': 'dmso',
  'dimethylformamide': 'dmf',
  'n,n-dimethylformamide': 'dmf',
  'tetrahydrofuran': 'thf',
  'dce': '1,2-dichloroethane',
  '1,2-dce': '1,2-dichloroethane',
}

#------------------------------------------------------------------------------
def canonical_solvent(name):
  """Normalize lookup aliases; ALPB choices remain independently validated."""
  key=str(name).strip().lower()
  return SOLVENT_ALIASES.get(key,key)

#------------------------------------------------------------------------------
def positive_probe_radius(value):
  """Require a finite positive probe radius in Angstrom."""
  radius=float(value)
  if not np.isfinite(radius) or radius<=0:
    raise ValueError('SASA probe radius must be positive and finite, in Angstrom.')
  return radius

#------------------------------------------------------------------------------
def resolve_probe_radius(solvent,manual_radius=None):
  """Prefer an explicit override, then a preset, otherwise require user input."""
  name=canonical_solvent(solvent)
  if manual_radius is not None:
    radius=positive_probe_radius(manual_radius)
    source='explicit --probe-radius override'
    method='manual_override'
  elif name in SOLVENT_PROBE_RADII:
    radius=SOLVENT_PROBE_RADII[name]
    source=PROBE_RADIUS_SOURCE
    method='source_documented_preset'
  else:
    print(f'No source-documented SASA probe radius in the local library for {solvent!r}.')
    if not sys.stdin.isatty():
      raise ValueError(f'Noninteractive run: supply --probe-radius for solvent {solvent!r}; no water-radius fallback is used.')
    print('Enter an effective spherical SASA probe radius in Angstrom.')
    while True:
      try:
        radius=positive_probe_radius(input('SASA probe radius (Angstrom): ').strip())
        break
      except ValueError as exc:
        print(exc)
    source='interactive user input'
    method='interactive_input'
  metadata=dict(solvent_key=name,radius_A=radius,selection_method=method,source=source,
    scope='FreeSASA geometry only; xTB ALPB parameters are unchanged')
  print(f'SASA probe radius: {radius:.3f} A; source: {source}')
  return radius,metadata

#------------------------------------------------------------------------------
def print_probe_library():
  """Show provenance and distinguish probe presets from ALPB availability."""
  print('LOCAL SASA PROBE LIBRARY (Angstrom)')
  print('Source:',PROBE_RADIUS_SOURCE)
  print('Solvent                  Probe(A)  In this script ALPB list')
  for name,radius in sorted(SOLVENT_PROBE_RADII.items()):
    print(f'{name:24s} {radius:7.2f}  {"yes" if name in SOLVENTS else "no"}')
  print('These are imported model conventions, not FreeSASA defaults or ALPB cavity parameters.')
  print('A probe preset does not imply support by the installed xTB version.')
  missing=sorted(set(SOLVENTS)-set(SOLVENT_PROBE_RADII))
  print('ALPB names without a local probe preset (manual radius required):',', '.join(missing))


#------------------------------------------------------------------------------
def read_xyz(filename):
  """Read one complete XYZ molecule, preserving all explicitly listed atoms."""
  lines=Path(filename).read_text(encoding='utf-8').splitlines()
  if len(lines)<3:
    raise ValueError(f'Incomplete XYZ: {filename}')
  n=int(lines[0].strip())
  if n<1 or len(lines)<n+2 or any(x.strip() for x in lines[n+2:]):
    raise ValueError(f'Provide a single complete XYZ frame: {filename}')
  symbols,coords=[],[]
  for i,line in enumerate(lines[2:n+2],1):
    parts=line.split()
    if len(parts)<4:
      raise ValueError(f'Incomplete XYZ atom {i}: {filename}')
    symbols.append(parts[0].capitalize())
    coords.append([float(x.replace('D','E').replace('d','e')) for x in parts[1:4]])
  coords=np.array(coords)
  if not np.all(np.isfinite(coords)):
    raise ValueError('Coordinates must be finite.')
  for i in range(n-1):
    if np.any(np.linalg.norm(coords[i+1:]-coords[i],axis=1)<1e-8):
      raise ValueError('Coincident atoms found in XYZ.')
  return symbols,coords

#------------------------------------------------------------------------------
def centered_geometry(symbols,coords):
  """Use a consistent center-of-mass origin, especially for molecular ions."""
  from rdkit import Chem
  table=Chem.GetPeriodicTable()
  masses=np.array([table.GetAtomicWeight(x) for x in symbols])
  return coords-np.sum(coords*masses[:,None],axis=0)/masses.sum()

#------------------------------------------------------------------------------
def write_xyz(filename,symbols,coords):
  with open(filename,'w',encoding='utf-8') as f:
    f.write(str(len(symbols))+'\nFixed geometry; center-of-mass origin; Angstrom.\n')
    for symbol,point in zip(symbols,coords):
      f.write(symbol+' '+' '.join(f'{x:.12f}' for x in point)+'\n')

#------------------------------------------------------------------------------
def sasa_analysis(symbols,coords,probe_radius=1.4,points=960):
  """Use explicit small-molecule atomic radii; no protein atom classifier."""
  import freesasa
  from rdkit import Chem
  table=Chem.GetPeriodicTable()
  radii=[float(table.GetRvdw(x)) for x in symbols]
  if any(not np.isfinite(x) or x<=0 for x in radii):
    raise ValueError('Invalid RDKit van der Waals radius; supply a validated radius scheme.')
  parameters=freesasa.Parameters({'algorithm':freesasa.ShrakeRupley,
    'probe-radius':probe_radius,'n-points':points,'n-threads':1})
  result=freesasa.calcCoord(coords.ravel().tolist(),radii,parameters)
  areas=[float(result.atomArea(i)) for i in range(len(symbols))]
  total=float(result.totalArea())
  if not np.isfinite(total) or total<0 or not np.isclose(sum(areas),total,rtol=1e-6,atol=1e-6):
    raise ValueError('Invalid FreeSASA result.')
  ch=float(sum(area for element,area in zip(symbols,areas) if element in ('C','H')))
  return dict(total_A2=total,CH_area_A2=ch,other_area_A2=total-ch,
    atom_areas_A2=areas,radii_A=radii,probe_radius_A=probe_radius,points=points)

#------------------------------------------------------------------------------
def parse_solvation_terms(text):
  """Read xTB model contributions, preserving missing fields as unavailable."""
  terms={}
  for name in ('Gsolv','Gelec','Gsasa','Ghb','Gshift'):
    matches=re.findall(r'\b'+name+r'\s+('+NUMBER+r')\s+Eh\b',text)
    terms[name]=float(matches[-1].replace('D','E').replace('d','e')) if matches else None
  return terms

#------------------------------------------------------------------------------
def run_xtb(symbols,coords,directory,args,solvent=None):
  """Run an isolated, fixed-geometry single point and retain all input/output."""
  directory.mkdir()
  write_xyz(directory/'input.xyz',symbols,coords)
  command=[args.xtb,'input.xyz','--gfn',str(args.gfn),'--sp','--json',
    '--chrg',str(args.charge),'--uhf',str(args.uhf),'--acc',str(args.accuracy)]
  if solvent is not None:
    command.extend(['--alpb',solvent])
    if args.state!='gsolv':
      command.append(args.state)
  env=dict(os.environ)
  env.update(OMP_NUM_THREADS=str(args.threads),MKL_NUM_THREADS=str(args.threads),
    OMP_STACKSIZE=env.get('OMP_STACKSIZE','1G'))
  (directory/'command.json').write_text(json.dumps(command,indent=2)+'\n')
  with open(directory/'xtb.out','w') as stdout,open(directory/'xtb.err','w') as stderr:
    try:
      completed=subprocess.run(command,cwd=directory,env=env,stdout=stdout,
        stderr=stderr,timeout=args.timeout,check=False)
    except subprocess.TimeoutExpired as exc:
      raise RuntimeError(f'xTB timed out; inspect {directory}') from exc
  text=(directory/'xtb.out').read_text(errors='replace')
  errors=(directory/'xtb.err').read_text(errors='replace')
  combined=text+'\n'+errors
  if completed.returncode!=0 or not re.search(r'(?<!ab)normal termination of xtb',combined,re.I):
    raise RuntimeError(f'xTB did not terminate normally; inspect {directory}')
  if re.search(r'SCC.{0,30}(?:not converged|failed to converge)|convergence failure',combined,re.I):
    raise RuntimeError(f'xTB SCC convergence failure; inspect {directory}')
  path=directory/'xtbout.json'
  if not path.exists():
    raise RuntimeError(f'xTB JSON output missing; inspect {directory}')
  data=json.loads(path.read_text())
  energy=float(data['total energy'])
  dipole=np.asarray(data['dipole / a.u.'],dtype=float)
  if dipole.shape!=(3,) or not np.isfinite(energy) or not np.all(np.isfinite(dipole)):
    raise ValueError(f'Invalid xTB energy/dipole data: {directory}')
  terms=parse_solvation_terms(text)
  if solvent is not None and terms['Gsolv'] is None:
    raise RuntimeError(f'ALPB Gsolv contribution not found; inspect {directory}')
  if solvent is not None and all(v is not None for v in terms.values()):
    discrepancy=terms['Gsolv']-sum(terms[x] for x in ('Gelec','Gsasa','Ghb','Gshift'))
    if abs(discrepancy)>1e-7:
      print(f'WARNING: ALPB printed decomposition does not close: {discrepancy:.3e} Eh.')
  return dict(energy_Eh=energy,dipole_vector_au=dipole.tolist(),
    dipole_D=float(np.linalg.norm(dipole)*AU_TO_DEBYE),model_terms_Eh=terms,
    directory=str(directory))

#------------------------------------------------------------------------------
def relative_energies(ref,probe):
  """Keep conformer preference separate from solvent stabilization preference."""
  rg,rs=ref['gas']['energy_Eh'],ref['solution']['energy_Eh']
  pg,ps=probe['gas']['energy_Eh'],probe['solution']['energy_Eh']
  return dict(ref_stabilization_kcal=(rs-rg)*HARTREE_TO_KCAL,
    probe_stabilization_kcal=(ps-pg)*HARTREE_TO_KCAL,
    relative_gas_kcal=(pg-rg)*HARTREE_TO_KCAL,
    relative_solution_kcal=(ps-rs)*HARTREE_TO_KCAL,
    relative_stabilization_kcal=((ps-pg)-(rs-rg))*HARTREE_TO_KCAL)

#------------------------------------------------------------------------------
def print_report(ref,probe,args):
  """Interpret descriptors cautiously and use calculated energies for preference."""
  values=relative_energies(ref,probe)
  print('\n===============================SOLVMOL====================================')
  print(f'GFN{args.gfn}-xTB / ALPB({args.solvent}); charge={args.charge}, UHF={args.uhf}, state={args.state}')
  print('Fixed input geometries. No optimization, ZPE, vibrational entropy or conformer ensemble.')
  print('\nPOLARITY DESCRIPTOR')
  for phase in ('gas','solution'):
    a,b=ref[phase]['dipole_D'],probe[phase]['dipole_D']
    print(f'Dipole {phase} (Debye): REF={a:.4f}, PROBE={b:.4f}, PROBE-REF={b-a:+.4f}')
  print('A larger dipole is a polarity descriptor, not proof of stronger solvent stabilization.')
  if args.charge:
    print('WARNING: charged-molecule dipoles are origin dependent; both geometries use their COM origin.')
  print('\nACCESSIBLE SURFACE (geometric descriptor, not nonpolar energy)')
  print(f'FreeSASA Shrake-Rupley; RDKit elemental VDW radii; probe={args.probe_radius:g} A; points={args.sasa_points}.')
  for name,label in [('total_A2','Total SASA'),('CH_area_A2','C/H atom surface'),('other_area_A2','Other atom surface')]:
    a,b=ref['sasa'][name],probe['sasa'][name]
    print(f'{label} (A^2): REF={a:.4f}, PROBE={b:.4f}, delta={b-a:+.4f}')
  area_ref,area_probe=ref['sasa']['total_A2'],probe['sasa']['total_A2']
  percent=100*(area_probe-area_ref)/area_ref if area_ref>1e-12 else None
  if percent is not None:
    label='similar exposure' if abs(percent)<=2 else 'greater exposure in PROBE' if percent>0 else 'lower exposure in PROBE'
    print(f'SASA interpretation: {label} (change={percent:+.2f}%; 2% reporting threshold).')
  print('Probe-radius source:',args.probe_radius_metadata['source'])
  print('Local probe presets/manual overrides affect FreeSASA only; xTB ALPB parameters are unchanged.')
  print('SASA comparisons across different probe radii use different surface definitions; compare conformers at the same radius.')
  print('C/H partition is element-based, not a rigorous polar/apolar chemical classification.')
  print('\nALPB MODEL CONTRIBUTIONS (kcal/mol)')
  for name in ('Gsolv','Gelec','Gsasa','Ghb','Gshift'):
    a,b=ref['solution']['model_terms_Eh'][name],probe['solution']['model_terms_Eh'][name]
    if a is None or b is None:
      print(name+': unavailable in one or both xTB outputs')
    else:
      print(f'{name:7s} REF={a*HARTREE_TO_KCAL:+.6f}, PROBE={b*HARTREE_TO_KCAL:+.6f}, delta={(b-a)*HARTREE_TO_KCAL:+.6f}')
  print('ALPB Gsasa is the model surface-area contribution; FreeSASA area is an independent diagnostic, not an extra energy term.')
  print('Gsolv is already included in xTB solution total energy; never add it again.')
  print('\nFIXED-GEOMETRY ENERGIES')
  for label,record in (('REF',ref),('PROBE',probe)):
    print(f"{label}: gas={record['gas']['energy_Eh']:.12f} Eh; solution={record['solution']['energy_Eh']:.12f} Eh")
  print(f"Solvent stabilization E_solution-E_gas (kcal/mol): REF={values['ref_stabilization_kcal']:+.6f}, PROBE={values['probe_stabilization_kcal']:+.6f}")
  print(f"Relative gas energy PROBE-REF: {values['relative_gas_kcal']:+.6f} kcal/mol")
  print(f"Relative solution energy PROBE-REF: {values['relative_solution_kcal']:+.6f} kcal/mol")
  print(f"Relative solvent stabilization PROBE-REF: {values['relative_stabilization_kcal']:+.6f} kcal/mol")
  for key,label in [('relative_stabilization_kcal','Solvent effect'),('relative_solution_kcal','Solution energy')]:
    difference=values[key]
    if abs(difference)<=args.energy_tolerance:
      print(f'{label}: close within the {args.energy_tolerance:g} kcal/mol reporting band; do not infer a robust preference.')
    elif key=='relative_stabilization_kcal':
      print('Solvent effect: preferentially stabilizes '+('PROBE' if difference<0 else 'REF')+' within this fixed-geometry model.')
    else:
      print('Solution energy: '+('PROBE' if difference<0 else 'REF')+' is lower in the model total-energy comparison.')
  print('Solution preference = gas-phase preference + differential solvent stabilization.')
  print('The solvent-minus-gas total-energy difference includes electronic relaxation; it need not equal printed Gsolv.')
  print('These are model fixed-geometry comparisons, not full conformational Gibbs free energies or direct solvent binding energies.')
  print('Transition-metal/spin-state rankings require validation; implicit solvent cannot describe explicit coordination or solvent bridges.')
  print('==========================================================================')
  return values

#==============================================================================
def main(argv=None):
  parser=argparse.ArgumentParser(description=__doc__)
  parser.add_argument('ref_xyz',nargs='?')
  parser.add_argument('probe_xyz',nargs='?')
  selection=parser.add_mutually_exclusive_group()
  selection.add_argument('--solvent',type=canonical_solvent,choices=SOLVENTS)
  for name in SOLVENTS:
    selection.add_argument('--'+name,dest='solvent',action='store_const',const=name)
  parser.add_argument('--list-solvents',action='store_true')
  parser.add_argument('--list-probe-radii',action='store_true',help='List local SASA probe radii and their source')
  parser.add_argument('--gfn',type=int,choices=(1,2),default=2)
  parser.add_argument('--charge',type=int,default=None,help='Common total charge; prompted if omitted')
  parser.add_argument('--uhf',type=int,default=None,help='Common number of unpaired electrons; prompted if omitted')
  parser.add_argument('--state',choices=('gsolv','bar1mol','reference'),default='gsolv')
  parser.add_argument('--probe-radius',type=float,default=None,
    help='Override the solvent-specific SASA probe radius, in Angstrom')
  parser.add_argument('--sasa-points',type=int,default=960)
  parser.add_argument('--threads',type=int,default=1)
  parser.add_argument('--accuracy',type=float,default=0.1)
  parser.add_argument('--timeout',type=float,default=3600.)
  parser.add_argument('--energy-tolerance',type=float,default=1.)
  parser.add_argument(
    '--xtb',
    default='/home/lky/miniconda3/envs/polygview/bin/xtb',
    help='xTB executable name or path'
  )
  parser.add_argument('--output-dir',help='New directory; existing directories are never overwritten')
  args=parser.parse_args(argv)
  if args.list_probe_radii:
    print_probe_library()
    if not args.list_solvents:
      return
  if args.list_solvents:
    print('Documented GFN1/2-ALPB solvent names (availability depends on installed xTB version):')
    print(', '.join(SOLVENTS))
    print('SASA radii use a separate local library; inspect --list-probe-radii or supply --probe-radius.')
    return
  if not args.ref_xyz or not args.probe_xyz or not args.solvent:
    parser.error('Provide ref.xyz probe.xyz and --solvent NAME or --NAME')
  if args.charge is None:
    args.charge=int(input('Enter common whole-molecule charge: ').strip())
  if args.uhf is None:
    args.uhf=int(input('Enter common number of unpaired electrons (xTB --uhf): ').strip())
  if args.uhf<0 or args.threads<1 or args.sasa_points<20:
    parser.error('UHF must be nonnegative, threads positive, SASA points >=20')
  args.probe_radius,args.probe_radius_metadata=resolve_probe_radius(args.solvent,args.probe_radius)
  for name in ('probe_radius','accuracy','timeout'):
    value=getattr(args,name)
    if not np.isfinite(value) or value<=0:
      parser.error(name+' must be positive and finite')
  if not np.isfinite(args.energy_tolerance) or args.energy_tolerance<0:
    parser.error('Energy reporting tolerance must be nonnegative and finite')
  executable=shutil.which(args.xtb)
  if executable is None:
    raise FileNotFoundError('xTB executable not found; install xtb or set --xtb /path/to/xtb')
  args.xtb=str(Path(executable).resolve())
  import freesasa
  from rdkit import Chem
  ra,ca=read_xyz(args.ref_xyz);pb,cb=read_xyz(args.probe_xyz)
  if Counter(ra)!=Counter(pb):
    raise ValueError('Conformer energy comparison requires identical elemental composition.')
  table=Chem.GetPeriodicTable()
  electrons=sum(table.GetAtomicNumber(x) for x in ra)-args.charge
  if electrons<1 or args.uhf>electrons or (electrons-args.uhf)%2:
    raise ValueError('Charge/UHF is inconsistent with total electron count.')
  ca,cb=centered_geometry(ra,ca),centered_geometry(pb,cb)
  # Fail on SASA/radius problems before launching external calculations.
  sasa_ref=sasa_analysis(ra,ca,args.probe_radius,args.sasa_points)
  sasa_probe=sasa_analysis(pb,cb,args.probe_radius,args.sasa_points)
  if args.output_dir:
    root=Path(args.output_dir).resolve()
    root.mkdir(parents=True,exist_ok=False)
  else:
    root=Path(tempfile.mkdtemp(prefix='solvmol_',dir=Path.cwd())).resolve()
  settings=dict(vars(args));settings['output_directory']=str(root)
  settings['note']='Same composition is checked; chemical identity and intended electronic state must be verified by the user.'
  (root/'settings.json').write_text(json.dumps(settings,indent=2)+'\n')
  version=subprocess.run([args.xtb,'--version'],capture_output=True,text=True,timeout=30,check=False)
  (root/'xtb_version.txt').write_text(version.stdout+'\n'+version.stderr)
  print('Calculation records:',root)
  print('Assumption: the two whole molecules are conformers of the same chemical species and intended spin state.')
  records=[]
  for label,symbols,coords,area in [('REF',ra,ca,sasa_ref),('PROBE',pb,cb,sasa_probe)]:
    print(f'Running {label} gas and ALPB({args.solvent}) single points...')
    gas=run_xtb(symbols,coords,root/(label.lower()+'_gas'),args)
    solution=run_xtb(symbols,coords,root/(label.lower()+'_solution'),args,args.solvent)
    records.append(dict(label=label,gas=gas,solution=solution,sasa=area))
  values=print_report(records[0],records[1],args)
  (root/'results.json').write_text(json.dumps(dict(settings=settings,conformers=records,relative=values),indent=2,allow_nan=False)+'\n')
  with open(root/'summary.csv','w',newline='',encoding='utf-8') as f:
    row=dict(solvent=args.solvent,gfn=args.gfn,charge=args.charge,uhf=args.uhf,
      probe_radius_A=args.probe_radius,probe_radius_source=args.probe_radius_metadata['source'],**values)
    writer=csv.DictWriter(f,fieldnames=list(row));writer.writeheader();writer.writerow(row)
  print('Results:',root/'results.json')
  return records,values


if __name__=='__main__':
  try:
    main()
  except (ValueError,RuntimeError,OSError,ImportError,KeyError,EOFError,subprocess.TimeoutExpired) as exc:
    print(f'{sys.argv[0]}: {exc}',file=sys.stderr)
    sys.exit(2)
