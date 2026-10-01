#!/home/lky/miniconda3/envs/polygview/bin/python
"""Export a 1-based geometry frame from a Gaussian optimization log.
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
def read_frames(logfile, job=1):
  """Keep one orientation convention within the selected job, without deduplication."""
  with open(logfile, encoding='utf-8', errors='replace') as f:
    jobs = split_jobs(f.readlines())
  if not 1 <= job <= len(jobs):
    raise ValueError(f'Job {job} is out of range; found {len(jobs)} job sections.')
  section = jobs[job-1]
  text = ''.join(line for _, line in section)
  if not re.search(r'Step number\s+\d+', text) and 'Optimization completed.' not in text:
    raise ValueError('Selected job has no recognized optimization records.')
  groups = {'Standard orientation': [], 'Input orientation': [], 'Z-Matrix orientation': []}
  for i, (_, line) in enumerate(section):
    title = line.strip().rstrip(':')
    if title in groups:
      groups[title].append(read_orientation(section, i))
  frames = next((groups[k] for k in groups if groups[k]), [])
  if not frames:
    raise ValueError('No complete orientation tables found.')
  elements = [atom[0] for atom in frames[0]['atoms']]
  if any([atom[0] for atom in frame['atoms']] != elements for frame in frames):
    raise ValueError('Atom count or element order changes within the selected job.')
  return frames, len(jobs)

#------------------------------------------------------------------------------
def log2xyz(logfile, frame, job=1, output=None, overwrite=False):
  """Export a forward-counted frame while preserving Gaussian atom order."""
  frames, jobs = read_frames(logfile, job)
  if not 1 <= frame <= len(frames):
    raise ValueError(f'Frame {frame} is out of range; available frames: 1-{len(frames)}.')
  selected = frames[frame-1]
  destination = Path(output) if output else Path.cwd() / f'{frame}_{Path(logfile).stem}.xyz'
  if destination.resolve() == Path(logfile).resolve():
    raise ValueError('Output must not overwrite the input log.')
  mode = 'w' if overwrite else 'x'
  with destination.open(mode, encoding='utf-8') as f:
    f.write(f"{len(selected['atoms'])}\n")
    f.write(f"Source={Path(logfile).name} Job={job} Frame={frame}/{len(frames)} "
      f"Orientation={selected['orientation']} LogLine={selected['line']} Units=Angstrom\n")
    for symbol, x, y, z in selected['atoms']:
      f.write(f'{symbol:<3s} {x: .8f} {y: .8f} {z: .8f}\n')
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
  parser.add_argument('-o', '--output')
  parser.add_argument('--overwrite', action='store_true')
  parser.add_argument('--list', action='store_true', help='List frames without exporting')
  args = parser.parse_args(argv)
  if not re.fullmatch(r'-?[1-9][0-9]*', args.frame):
    parser.error('Frame must be a positive index, optionally prefixed with a minus sign.')
  frame = abs(int(args.frame))
  if args.list:
    frames, jobs = read_frames(args.logfile, args.job)
    print(f'Job {args.job}/{jobs}; total frames: {len(frames)}')
    for i, item in enumerate(frames, 1):
      print(f"{i:5d}  line={item['line']:8d}  atoms={len(item['atoms']):5d}  {item['orientation']}")
    return
  log2xyz(args.logfile, frame, args.job, args.output, args.overwrite)


if __name__ == '__main__':
  try:
    main()
  except (OSError, ValueError) as exc:
    print(f'{sys.argv[0]}: {exc}', file=sys.stderr)
    sys.exit(2)
