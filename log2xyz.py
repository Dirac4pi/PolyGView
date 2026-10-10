#!/home/lky/miniconda3/envs/polygview/bin/python
"""Export a forward frame from Gaussian optimization or relaxed-scan logs.
Optimization: all geometry frames, unchanged. Relaxed scan: completed scan points.
Default filenames are N_INPUT-STEM.xyz; use -o to choose another name.
Author: Dirac4pi
env: polygview
"""
import argparse
import math
import re
import sys
from pathlib import Path

ELEMENTS = ('X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca '
  'Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo '
  'Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu '
  'Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po '
  'At Rn Fr Ra Ac Th Pa U Np Pu Am Cm Bk Cf Es Fm Md No Lr Rf '
  'Db Sg Bh Hs Mt Ds Rg Cn Nh Fl Mc Lv Ts Og').split()

#------------------------------------------------------------------------------
def split_jobs(lines):
  """Split at termination records without merging subsequent jobs."""
  jobs, current = [], []
  for number, line in enumerate(lines, 1):
    current.append((number, line))
    if 'Normal termination of Gaussian' in line or 'Error termination' in line:
      jobs.append(current)
      current = []
  if current and any('orientation:' in line for _, line in current):
    jobs.append(current)
  return jobs

#------------------------------------------------------------------------------
def read_orientation(section, start):
  """Read one complete orientation table; reject unsupported centers."""
  number, title = section[start]
  i, separators = start + 1, 0
  while i < len(section) and separators < 2:
    line = section[i][1].strip()
    if line and set(line) == {'-'}:
      separators += 1
    i += 1
  if separators != 2:
    raise ValueError(f'Incomplete orientation header at line {number}.')
  atoms = []
  while i < len(section):
    line = section[i][1].strip()
    if line and set(line) == {'-'}:
      if not atoms:
        raise ValueError(f'Empty orientation table at line {number}.')
      return {'line': number, 'orientation': title.strip().rstrip(':'),
        'atoms': atoms}
    parts = line.split()
    if len(parts) != 6:
      raise ValueError(f'Invalid orientation row at line {section[i][0]}.')
    center, atomic_number = int(parts[0]), int(parts[1])
    if center != len(atoms) + 1:
      raise ValueError(f'Nonsequential center numbering at line {section[i][0]}.')
    if not 1 <= atomic_number < len(ELEMENTS):
      raise ValueError('Dummy, ghost or translation-vector centers are not supported; '
        'refusing to change atom numbering silently.')
    xyz = tuple(float(v.replace('D', 'E').replace('d', 'e')) for v in parts[3:])
    if not all(math.isfinite(v) for v in xyz):
      raise ValueError(f'Non-finite coordinate at line {section[i][0]}.')
    atoms.append((ELEMENTS[atomic_number], *xyz))
    i += 1
  raise ValueError(f'Incomplete orientation table at line {number}.')

#------------------------------------------------------------------------------
def scan_metadata(section):
  """Detect actual scan records; Opt=ModRedundant alone is not a scan."""
  text = ''.join(line for _,line in section)
  points = list(re.finditer(r'\bon scan point\s+(\d+)\s+out of\s+(\d+)',text,re.I))
  declared = re.search(r'Number of optimizations in scan\s*=\s*(\d+)',text,re.I)
  directive = re.search(
    r'^\s*(?:[BADLXYZ]\s+)?(?:[0-9*]+\s+){1,4}S\s+\d+\s+[-+0-9.DEded]+\s*$',text,re.M|re.I)
  parameter = re.search(r'^\s*!\s*\w+\s+\w+\([^)]*\).*?\bScan\b',text,re.M|re.I)
  is_scan = bool(points or declared or directive or parameter)
  planned = int(declared.group(1)) if declared else int(points[0].group(2)) if points else None
  return is_scan, planned


def collect_orientations(section, scan=False):
  """Preserve orientation priority; tolerate only an unfinished scan tail."""
  groups = {'Standard orientation': [], 'Input orientation': [], 'Z-Matrix orientation': []}
  ends = [number for number,line in section if 'Optimization completed.' in line]
  last_end = max(ends,default=0)
  for i,(number,line) in enumerate(section):
    title = line.strip().rstrip(':')
    if title not in groups:
      continue
    try:
      item = read_orientation(section,i)
    except ValueError as exc:
      if scan and number > last_end and str(exc).startswith('Incomplete orientation'):
        continue
      raise
    groups[title].append(item)
  return groups


def check_elements(frames):
  """Reject changing atom counts/order without renumbering or deduplication."""
  if not frames:
    raise ValueError('No complete orientation tables found.')
  elements = [atom[0] for atom in frames[0]['atoms']]
  if any([atom[0] for atom in item['atoms']] != elements for item in frames):
    raise ValueError('Atom count or element order changes within the selected job.')


def completed_scan_frames(section,groups,planned):
  """Take the last preceding geometry for each successfully optimized scan point."""
  has_completed = any('Optimization completed.' in line for _,line in section)
  completions, point = [], None
  for number,line in section:
    match = re.search(r'\bon scan point\s+(\d+)\s+out of\s+(\d+)',line,re.I)
    if match:
      point = int(match.group(1))
      if planned is None:
        planned = int(match.group(2))
    done = ('Optimization completed.' in line if has_completed
      else '-- Stationary point found.' in line)
    if done:
      completions.append((number,point))
  if not completions:
    raise ValueError('Relaxed scan detected, but no successfully completed scan points were found.')
  if any(point is None for _,point in completions):
    if any(point is not None for _,point in completions):
      raise ValueError('Mixed missing/explicit scan-point numbers; refusing ambiguous scan indexing.')
    completions = [(number,index) for index,(number,_) in enumerate(completions,1)]
  identifiers = [point for _,point in completions]
  if len(set(identifiers)) != len(identifiers) or identifiers != sorted(identifiers):
    raise ValueError('Repeated or nonascending scan-point numbers; split separate scan sequences before exporting.')
  for orientation,raw in groups.items():
    if not raw:
      continue
    selected, cursor, previous_end = [], 0, 0
    for completion,point in completions:
      latest = None
      while cursor < len(raw) and raw[cursor]['line'] < completion:
        if raw[cursor]['line'] > previous_end:
          latest = raw[cursor]
        cursor += 1
      if latest is None:
        break
      item = dict(latest,kind='relaxed_scan',scan_point=point,
        completion_line=completion,planned_scan_points=planned,
        raw_orientation_count=len(raw))
      selected.append(item)
      previous_end = completion
    if len(selected) == len(completions):
      check_elements(selected)
      return selected
  raise ValueError('No consistent orientation convention supplies a geometry for every completed scan point.')


#------------------------------------------------------------------------------
def read_frames(logfile, job=1, mode='auto'):
  """Use raw optimization frames or completed relaxed-scan points automatically."""
  if mode not in ('auto','optimization','scan'):
    raise ValueError('Mode must be auto, optimization or scan.')
  with open(logfile,encoding='utf-8',errors='replace') as f:
    jobs = split_jobs(f.readlines())
  if not 1 <= job <= len(jobs):
    raise ValueError(f'Job {job} is out of range; found {len(jobs)} job sections.')
  section = jobs[job-1]
  text = ''.join(line for _,line in section)
  detected,planned = scan_metadata(section)
  scan = mode == 'scan' or (mode == 'auto' and detected)
  if not scan and not re.search(r'Step number\s+\d+',text) and 'Optimization completed.' not in text:
    raise ValueError('Selected job has no recognized optimization records.')
  groups = collect_orientations(section,scan)
  if scan:
    frames = completed_scan_frames(section,groups,planned)
  else:
    frames = next((groups[key] for key in groups if groups[key]),[])
    check_elements(frames)
  return frames,len(jobs)


#------------------------------------------------------------------------------
def log2xyz(logfile, frame, job=1, output=None, overwrite=False, mode='auto'):
  """Export a forward optimization frame or an explicitly numbered scan point."""
  frames,jobs = read_frames(logfile,job,mode)
  scan = frames[0].get('kind') == 'relaxed_scan'
  if scan:
    selected = next((item for item in frames if item['scan_point'] == frame),None)
    if selected is None:
      available = ','.join(str(item['scan_point']) for item in frames)
      raise ValueError(f'Scan point {frame} is not available as a completed point; available points: {available}. '
        'Unfinished or failed points are not exported and later points are not renumbered.')
  else:
    if not 1 <= frame <= len(frames):
      raise ValueError(f'Frame {frame} is out of range; available frames: 1-{len(frames)}.')
    selected = frames[frame-1]
  name = f'{frame}_{Path(logfile).stem}.xyz'
  destination = Path(output) if output else Path.cwd()/name
  if destination.resolve() == Path(logfile).resolve():
    raise ValueError('Output must not overwrite the input log.')
  with destination.open('w' if overwrite else 'x',encoding='utf-8') as f:
    f.write(f"{len(selected['atoms'])}\n")
    if scan:
      planned = selected['planned_scan_points']
      f.write(f"Source={Path(logfile).name} Job={job} Mode=RelaxedScan ScanPoint={frame} "
        f"CompletedPoints={len(frames)} PlannedPoints={planned if planned is not None else 'unknown'} "
        f"Orientation={selected['orientation']} LogLine={selected['line']} "
        f"CompletionLine={selected['completion_line']} Units=Angstrom\n")
    else:
      f.write(f"Source={Path(logfile).name} Job={job} Frame={frame}/{len(frames)} "
        f"Orientation={selected['orientation']} LogLine={selected['line']} Units=Angstrom\n")
    for symbol,x,y,z in selected['atoms']:
      f.write(f'{symbol:<3s} {x: .8f} {y: .8f} {z: .8f}\n')
  if scan:
    planned = selected['planned_scan_points']
    print(f'Job {job}/{jobs}; mode: relaxed scan; completed points: {len(frames)}; '
      f'planned points: {planned if planned is not None else "unknown"}; '
      f'orientation: {selected["orientation"]}')
    print(f'Exported scan point {frame}, log line {selected["line"]}: {destination}')
    if planned is not None and len(frames) < planned:
      print('Only completed scan points are available; intermediate/unfinished optimization geometries are excluded.')
  else:
    print(f'Job {job}/{jobs}; frames: {len(frames)}; orientation: {selected["orientation"]}')
    print(f'Exported frame {frame}, log line {selected["line"]}: {destination}')
  if jobs > 1:
    print('Other job sections are excluded from frame numbering (use --job to select one).')
  return destination

#==============================================================================
def main(argv=None):
  parser = argparse.ArgumentParser(description=__doc__)
  parser.add_argument('logfile')
  parser.add_argument('frame', help='Forward frame selector, e.g. -12 or 12')
  parser.add_argument('--job', type=int, default=1)
  parser.add_argument('--mode', choices=('auto','optimization','scan'), default='auto',
    help='Auto-detect by default; optimization forces the original raw-frame logic')
  parser.add_argument('-o', '--output')
  parser.add_argument('--overwrite', action='store_true')
  parser.add_argument('--list', action='store_true', help='List optimization frames or completed scan points without exporting')
  args = parser.parse_args(argv)
  if not re.fullmatch(r'-?[1-9][0-9]*', args.frame):
    parser.error('Frame must be a positive index, optionally prefixed with a minus sign.')
  frame = abs(int(args.frame))
  if args.list:
    frames, jobs = read_frames(args.logfile, args.job, args.mode)
    scan = frames[0].get('kind') == 'relaxed_scan'
    print(f'Job {args.job}/{jobs}; mode: {"relaxed scan" if scan else "optimization"}; total frames: {len(frames)}')
    for i,item in enumerate(frames,1):
      if scan:
        print(f"{item['scan_point']:5d} line={item['line']:8d} atoms={len(item['atoms']):5d} "
          f"{item['orientation']} completed_at={item['completion_line']}")
      else:
        print(f"{i:5d} line={item['line']:8d} atoms={len(item['atoms']):5d} {item['orientation']}")
    return
  log2xyz(args.logfile, frame, args.job, args.output, args.overwrite, args.mode)


if __name__ == '__main__':
  try:
    main()
  except (OSError, ValueError) as exc:
    print(f'{sys.argv[0]}: {exc}', file=sys.stderr)
    sys.exit(2)
