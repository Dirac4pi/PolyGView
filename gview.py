'''
GaussView Laucher with support for xyz editing
Author: Dirac4pi
env:base
'''

from os import path, remove
from re import search
from subprocess import call

elements = {
  'H' : 1, 'He': 2, 'Li': 3, 'Be': 4, 'B' : 5, 'C' : 6, 'N' : 7, 'O' : 8,
  'F' : 9, 'Ne':10, 'Na':11, 'Mg':12, 'Al':13, 'Si':14, 'P' :15, 'S' :16,
  'Cl':17, 'Ar':18, 'K' :19, 'Ca':20, 'Sc':21, 'Ti':22, 'V' :23, 'Cr':24,
  'Mn':25, 'Fe':26, 'Co':27, 'Ni':28, 'Cu':29, 'Zn':30, 'Ga':31, 'Ge':32,
  'As':33, 'Se':34, 'Br':35, 'Kr':36, 'Rb':37, 'Sr':38, 'Y' :39, 'Zr':40,
  'Nb':41, 'Mo':42, 'Tc':43, 'Ru':44, 'Rh':45, 'Pd':46, 'Ag':47, 'Cd':48,
  'In':49, 'Sn':50, 'Sb':51, 'Te':52, 'I' :53, 'Xe':54, 'Cs':55, 'Ba':56,
  'La':57, 'Ce':58, 'Pr':59, 'Nd':60, 'Pm':61, 'Sm':62, 'Eu':63, 'Gd':64,
  'Tb':65, 'Dy':66, 'Ho':67, 'Er':68, 'Tm':69, 'Yb':70, 'Lu':71, 'Hf':72,
  'Ta':73, 'W' :74, 'Re':75, 'Os':76, 'Ir':77, 'Pt':78, 'Au':79, 'Hg':80,
  'Tl':81, 'Pb':82, 'Bi':83, 'Po':84, 'At':85, 'Rn':86, 'Fr':87, 'Ra':88,
  'Ac':89, 'Th':90, 'Pa':91, 'U' :92, 'Np':93, 'Pu':94, 'Am':95, 'Cm':96,
  'Bk':97, 'Cf':98, 'Es':99
}

#-------------------------------------------------------------------------------
def isint(string: str) -> bool:
  '''
  Determine whether the string is an integer or not
  --
  :string: input string\n
  return: true/false
  '''
  try:
    int(string)
  except ValueError:
    return False
  else:
    return True

#-------------------------------------------------------------------------------
xyzlines = None
Nframe = 0
framestarts = []
frameends = []
def load_xyz(xyz_file:str, frame:int):
  '''
  Load specific frame information from .xyz/.trj.
  --
  :xyz_file: input file name\n
  :frame: index of specific frame\n
  :Natom: number of atoms of specific frame\n
  :title: title of specific frame\n
  :coord: atom coordinates of specific frame\n
  return: Natom, title, coord, fcoords
  '''
  import numpy as np
  global xyzlines, Nframe, framestarts, frameends
  if xyzlines is None:
    try:
      open(xyz_file, 'r')
    except:
      print(f"Error: File {xyz_file} can not be opened.")
      exit(1)
    else:
      with open(xyz_file, 'r') as f:
        xyzlines = f.readlines()
    i = 0
    while i < len(xyzlines):
      if isint(xyzlines[i].strip()):
        framestarts.append(i)
        i += int(xyzlines[i].strip()) + 2
        frameends.append(i)
        Nframe += 1
      else:
        i += 1
  if frame <= 0 or frame > Nframe:
    print("Error: Frame index out of range.")
    exit(1)
  ii = framestarts[frame-1]
  jj = frameends[frame-1]
  Natom = int(xyzlines[ii].strip())
  title = xyzlines[ii+1].strip()
  coord = xyzlines[ii+2:jj]
  fcoords = np.zeros((Natom, 3), dtype=float)
  # Parse coordinates starting from the third line
  for i in range(Natom):
    parts = coord[i].split()
    fcoords[i] = [float(parts[1]), float(parts[2]), float(parts[3])]
  return Natom, title, coord, fcoords

#-------------------------------------------------------------------------------
def single_xyz_edit(xyz_file:str) -> None:
  '''
  Visualise and edit single-frame xyz file via GaussView, support lattice info.
  --
  :xyz_file: input file name\n
  return: None
  '''
  Natom, title, coord, fcoords = load_xyz(xyz_file, 1)
  # Generate .gjf (coordinates only)
  gjf_file = path.splitext(xyz_file)[0] + "_fromxyz.gjf"
  with open(gjf_file, 'w') as f:
    gjfheader = f'''
%chk=fromxyz.chk
%mem=1500MB
%nprocshared=4
#p opt freq upbe1pbe em=gd3bj def2svp

Converted from {xyz_file}

0 1
'''
    f.write(gjfheader)
    for line in coord:
      f.write(line.strip() + '\n')
    if search(r'Lattice="([^"]*)"', title):
      print('this is a periodic system')
      lattice_info = search(r'Lattice="([^"]*)"', title).group(1)
      lattice_data = lattice_info.split()
      f.write('Tv               '+lattice_data[0]+'   '+lattice_data[1]+'   '+\
              lattice_data[2]+'\n')
      f.write('Tv               '+lattice_data[3]+'   '+lattice_data[4]+'   '+\
              lattice_data[5]+'\n')
      f.write('Tv               '+lattice_data[6]+'   '+lattice_data[7]+'   '+\
              lattice_data[8]+'\n')
    f.write('\n')
  # Launch gview
  try:
    call(['/home/lky/gv/gview.sh', gjf_file])
  except Exception as e:
    print(f"GaussView error: {e}")
    return
  # Read edited coords
  with open(gjf_file, 'r') as f:
    edited_lines = f.readlines()
  edited_coords = []
  in_coords = False
  for line in edited_lines:
    stripped = line.strip()
    if stripped and stripped.split()[0].isdigit() and len(stripped.split())==2:
      in_coords = True
      continue
    if in_coords and stripped and stripped[0].isalpha() and \
      stripped.find('Tv')!=0:
      edited_coords.append(stripped)
    if in_coords and not stripped:
      break
  # Write output xyz (with original metadata)
  edited_xyz = path.splitext(xyz_file)[0] + ".xyz"
  with open(edited_xyz, 'w') as f:
    f.write(f"{len(edited_coords)}\n")
    f.write(f"{title}\n")
    for coord in edited_coords:
      f.write(coord + "\n")
  remove(gjf_file)
  print(f"Edited structure saved to {edited_xyz}")

#-------------------------------------------------------------------------------
def log_geom(xyz_coords:str):
  '''
  Docstring for log_geom
  --
  :xyz_coords: atom coordinates in .xyz/.trj file\n
  return: log_coords
  '''
  global elements
  log_coords = ''
  for i in range(len(xyz_coords)):
    index = elements.get(xyz_coords[i].split()[0], 0)
    x = float(xyz_coords[i].split()[1])
    y = float(xyz_coords[i].split()[2])
    z = float(xyz_coords[i].split()[3])

    line = f'     {i+1:3d}         {index:3d}' + \
           f'           0       {x:10.6f}   {y:10.6f}   {z:10.6f}'
    if i != len(xyz_coords)-1:
      log_coords = log_coords + line + '\n'
    else:
      log_coords = log_coords + line
  return log_coords

#-------------------------------------------------------------------------------
def multi_xyz_visual(xyz_file:str) -> None:
  '''
  Visualise multi-frame xyz file via GaussView.
  --
  :xyz_file: input file name\n
  return: None
  '''
  from drmsd import drmsd, dmaxd
  log_file = path.splitext(xyz_file)[0] + "_fromxyz.log"
  with open(log_file, 'w') as f:
    header = \
f''' ! converted from {xyz_file}

 0 basis functions
 0 alpha electrons
 0 beta electrons
GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad'''
    f.write(header)
    f.write("\n")
    for i in range(1, Nframe+1):
      if i != 1:
        fcoords_prev = fcoords
      Natom, title, coords, fcoords = load_xyz(xyz_file, i)
      if i == 1:
        vdrmsd = 0.0
        vdmaxd = 0.0
      else:
        vdrmsd = drmsd(fcoords_prev, fcoords)
        vdmaxd = dmaxd(fcoords_prev, fcoords)
      try:# for ORCA and CP2K, last value of title line is energy.
        energy = float(title.split()[-1])
      except ValueError:
        energy = -675.0
      log_coords = log_geom(coords)
      framer = \
f'''GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad
 
                         Standard orientation:
 ---------------------------------------------------------------------
 Center     Atomic      Atomic             Coordinates (Angstroms)
 Number     Number       Type             X           Y           Z
 ---------------------------------------------------------------------
{log_coords}
 ---------------------------------------------------------------------
 
 SCF Done:      {energy:.7f}     A.U. after   10 cycles
 
 GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad
 Step number   {i}
         Item               Value     Threshold  Converged?
 Maximum Force            1.000000     1.000000     NO
 RMS     Force            1.000000     1.000000     NO
 Maximum Displacement     {vdmaxd:.7f}     1.000000     NO
 RMS     Displacement     {vdrmsd:.7f}     1.000000     NO'''
      f.write(framer)
      f.write("\n")
    footer = \
'''GradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGradGrad
 Normal termination of Gaussian'''
    f.write(footer)
  # Launch gview, but editing is not allowed.
  try:
    call(['/home/lky/gv/gview.sh', log_file])
  except Exception as e:
    print(f"GaussView error: {e}")
    return
  remove(log_file)

#-------------------------------------------------------------------------------
def out_freq_visual(outfile:str) -> None:
  '''
  Visualise vibrational modes from Gaussian/ORCA output file via GaussView.
  --
  :outfile: input output file name\n
  return: None
  '''
  try:
    open(outfile, 'r')
  except:
    print(f"Error: File {outfile} can not be opened.")
    exit(1)
  else:
    found = ''
    with open(outfile, 'r') as f:
      line_count = 0
      for line in f:
        if 'Gaussian' in line:
          print('This is Gaussian output file')
          found = 'Gaussian'
          break
        elif 'O   R   C   A' in line or 'Frank Neese' in line:
          print('This is ORCA output file')
          found = 'ORCA'
          break
        line_count += 1
        if line_count >= 100 and found == '':
          break
    if found == '' and line_count > 0:
      print("Unrecognize output file type.")
      exit(1)
    elif found == 'Gaussian':
      call(['/home/lky/gv/gview.sh', outfile])
    elif found == 'ORCA':
      call(['/home/lky/gv/OfakeG', outfile])  # convert ORCA output by OfakeG
      outfile_gau = outfile.removesuffix('.out') + "_fake.out"
      call(['/home/lky/gv/gview.sh', outfile_gau])
      remove(outfile_gau)

#-------------------------------------------------------------------------------
def MOLDEN_freq_visual(MOLDEN_file:str) -> None:
  '''
  Visualise vibrational modes from CP2K MOLDEN file via GaussView
  --
  :MOLDEN_file: input CP2K MOLDEN file name\n
  return: None
  '''
  try:
    open(MOLDEN_file, 'r')
  except:
    print(f"Error: File {MOLDEN_file} can not be opened.")
    exit(1)
  else:
    call(['/home/lky/gv/MfakeG', MOLDEN_file])
    outfile_gau = MOLDEN_file.removesuffix('.mol') + "_fake.out"
    call(['/home/lky/gv/gview.sh', outfile_gau])
    remove(outfile_gau)

#-------------------------------------------------------------------------------
def cdxml_visual(cdxml_file:str) -> None:
  '''
  Visualise ChemDraw XML file via GaussView
  --
  :cdxml_file: input ChemDraw XML file name\n
  return: None
  '''
  try:
    open(cdxml_file, 'r')
  except:
    print(f"Error: File {cdxml_file} can not be opened.")
    exit(1)
  else:
    gjf_file = cdxml_file.removesuffix('.cdxml') + ".gjf"
    call(['obabel',
            cdxml_file,
            '-O', gjf_file,
            '--gen3d',
            '--minimize',
            '--ff', 'MMFF94',
            '-h'             ])
    call(['/home/lky/gv/gview.sh', gjf_file])

#===============================================================================
if __name__ == "__main__":
  import sys
  if len(sys.argv) == 1:
    call(['/home/lky/gv/gview.sh'])
  elif len(sys.argv) != 2:
    print("Usage: gview input.[xyz|trj|gjf|log|out|cif|fch|mol]")
    exit(1)
  else:
    input_file = sys.argv[1]
    if input_file.endswith('.xyz'):
      Natom, title, coord, fcoords = load_xyz(input_file, 1)
      print(f'{Nframe} frames found in {input_file}')
      if Nframe == 1:
        single_xyz_edit(input_file)
      elif Nframe > 1:
        multi_xyz_visual(input_file)
      else:
        print("Error: No valid frames found in the xyz file.")
        exit(1)
    elif input_file.endswith('.trj'):
      Natom, title, coord, fcoords = load_xyz(input_file, 1)
      print(f'{Nframe} frames found in {input_file}')
      if Nframe == 1:
        single_xyz_edit(input_file)
      elif Nframe > 1:
        multi_xyz_visual(input_file)
      else:
        print("Error: No valid frames found in the xyz file.")
        exit(1)
    elif input_file.endswith('.gjf'):
      call(['/home/lky/gv/gview.sh', input_file])
    elif input_file.endswith('.log'):
      call(['/home/lky/gv/gview.sh', input_file])
    elif input_file.endswith('.out'):
      out_freq_visual(input_file)
    elif input_file.endswith('.cif'):
      call(['/home/lky/gv/gview.sh', input_file])
    elif input_file.endswith('.fch'):
      call(['/home/lky/gv/gview.sh', input_file])
    elif input_file.endswith('.molden') or input_file.endswith('.molden.input'):
      MOLDEN_freq_visual(input_file)
    elif input_file.endswith('.mol') or input_file.endswith('.mol2'):
      call(['/home/lky/gv/gview.sh', input_file])
    elif input_file.endswith('.pdb'):
      call(['/home/lky/gv/gview.sh', input_file])
    elif input_file.endswith('.cub'):
      call(['/home/lky/gv/gview.sh', input_file])
    elif input_file.endswith('.cdxml'):
      cdxml_visual(input_file)
    else:
      print("unrecognize input file type")
      exit(1)
