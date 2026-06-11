#!/home/lky/miniconda3/envs/polygview/bin/python
'''
GaussView Laucher with support for xyz editing
Author: Dirac4pi
env:polygview
'''

import sys
from os import path, remove, environ
from re import search
from subprocess import call, run, PIPE

if sys.platform.startswith('win'):
  exe_name = 'gview.exe'
else:
  exe_name = 'gview.sh'
GV = environ.get("GV_ROOT", "/opt/gv")
exe_path = path.join(GV, exe_name)

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
  import time
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
  while True:
    # ensure that no gview running when editing an XYZ file
    if sys.platform.startswith('win'):
      check_process = run(
          ['tasklist', '/FI', f'IMAGENAME eq gview.exe', '/NH'], 
          stdout=PIPE,text=True
      )
      if not 'gview.exe' in check_process.stdout.lower():
        break
    else:
      check_process = run(['pgrep', '-f', 'gview.exe'], stdout=PIPE)
      if not check_process.stdout:
        break
    print(f"\nGView seems running, it will trigger Inter-Process-Communication"\
          +" (IPC) issues •_•")
    print(f"Please SAVE AND CLOSE all the GView process to proceed!")
    input("Press ENTER after they're closed ...")
    time.sleep(0.2) # time buffer to close gview
  print(f"calling gview ...")
  try:
    call([exe_path, gjf_file])
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
def visual_orca_inp(inp_file: str) -> None:
  """
  Visualise and edit ORCA input file via GaussView
  --
  inp_file (str): Path to the original ORCA input file.\n
  return: None
  """
  import time
  header_lines = []
  atoms = []
  footer_lines = []
  in_coord_block = False
  coord_block_finished = False
  with open(inp_file, 'r') as f_inp:
    for line in f_inp:
      clean_line = line.strip()
      # Extract everything before and including the '* xyz [C] [M]' line
      if not in_coord_block and not coord_block_finished:
        header_lines.append(line)
        if clean_line.lower().startswith("* xyz"):
          if clean_line.lower().startswith("* xyzfile"):
            print("this ORCA input file has external coordinate reference")
            exit(1)
          in_coord_block = True
      elif in_coord_block:
        if clean_line == "*":
          in_coord_block = False
          coord_block_finished = True
        elif clean_line:
          atoms.append(clean_line)
      elif coord_block_finished:
        footer_lines.append(line)
  xyz_file = path.splitext(inp_file)[0] + "_frominp.xyz"
  with open(xyz_file, 'w') as f_xyz:
    f_xyz.write(f"{len(atoms)}\n")
    f_xyz.write(f"Extracted from {inp_file}\n")
    for atom in atoms:
      f_xyz.write(f"{atom}\n")
  single_xyz_edit(xyz_file)
  time.sleep(1.0) # time buffer to ensure the edited xyz is saved
  modified_atoms = []
  with open(xyz_file, 'r') as f_xyz_out:
    lines = f_xyz_out.readlines()
    if len(lines) >= 3:
      num_atoms = int(lines[0].strip())
      modified_atoms = [l.strip() for l in lines[2:2 + num_atoms]]
  with open(inp_file, 'w') as f_out:
    for h_line in header_lines:
      f_out.write(h_line)
    for modified_atom in modified_atoms:
      f_out.write(f"  {modified_atom}\n")
    f_out.write("*\n")
    for f_line in footer_lines:
      f_out.write(f_line)
  if path.exists(xyz_file): remove(xyz_file)
  print(f"Updated {inp_file} saved")

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
      # at least for ORCA, geomeTRIC and CP2K, last value of title line is energy.
      try:
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
    call([exe_path, log_file])
  except Exception as e:
    print(f"GaussView error: {e}")
    return
  remove(log_file)

#-------------------------------------------------------------------------------
def OUT_visual(outfile:str) -> None:
  '''
  Visualise OPT, FREQ, OPT FREQ jobs from Gaussian/ORCA output file.
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
      call([exe_path, outfile])
    elif found == 'ORCA':
      # Polygview currently does not support the visualization of OUT files
      # generated by ORCA potential energy surface scan tasks.
      with open(outfile, 'r') as f:
        for line in f:
          if 'SURFACE SCAN STEP' in line:
            print("PolyGView (and also OfakeG) currently does not support the "\
            "visualization of OUT files generated by ORCA potential energy "\
            "surface scan tasks, please use the _trj.xyz file generated by "\
            "ORCA instead!")
            exit(1)
      if sys.platform.startswith('win'):
        OfakeG_path = path.join(GV, 'OfakeG.exe')
      else:
        OfakeG_path = path.join(GV, 'OfakeG')
      call([OfakeG_path, outfile])
      outfile_gau = outfile.removesuffix('.out') + "_fake.out"
      # add force on nuclears
      all_gauforces = []
      orcagrad = []
      with open(outfile, 'r') as f:
        while True:
          line = f.readline()
          if not line:
            break
          if line.startswith('CARTESIAN GRADIENT'):
            line = f.readline()
            line = f.readline()
            line = f.readline()
            while ":" in line:
              orcagrad.append(line.strip())
              line = f.readline()
            gauforce = \
'''\n -------------------------------------------------------------------
 Center     Atomic                   Forces (Hartrees/Bohr)
 Number     Number              X              Y              Z
 -------------------------------------------------------------------'''
            for i in range(len(orcagrad)):
              parts = orcagrad[i].split()
              index = int(parts[0])
              atomic_num = elements.get(parts[1], 0)
              fx = -float(parts[3])    # force = - gradient
              fy = -float(parts[4])
              fz = -float(parts[5])
              gauforce += f'\n {index:3d}       {atomic_num:3d}      ' +\
                          f'{fx:14.6f}   {fy:14.6f}   {fz:14.6f}'
            gauforce += \
'\n -------------------------------------------------------------------\n'
            all_gauforces.append(gauforce)
          orcagrad = []
      with open(outfile_gau, 'r') as fo:
        gau_lines = fo.readlines()
      new_gau_lines = []
      oframe = 0
      for lineo in gau_lines:
        new_gau_lines.append(lineo)
        if 'SCF Done:' in lineo:
          if oframe < len(all_gauforces):
            new_gau_lines.append(all_gauforces[oframe])
            oframe += 1
        # ORCA does not output atomic force information during frequency
        # calculations; by default, we use the same data as the final
        # frame of the optimization.
        if oframe >= 1 and 'Harmonic frequencies' in lineo:
          new_gau_lines.append(all_gauforces[oframe-1])
          oframe += 1
      with open(outfile_gau, 'w') as fo:
        fo.writelines(new_gau_lines)
      print('\nForces on nuclears have been added!')
      call([exe_path, outfile_gau])
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
    if sys.platform.startswith('win'):
      MfakeG_path = path.join(GV, 'MfakeG.exe')
    else:
      MfakeG_path = path.join(GV, 'MfakeG')
    call([MfakeG_path, MOLDEN_file])
    outfile_gau = MOLDEN_file.removesuffix('.mol') + "_fake.out"
    call([exe_path, outfile_gau])
    remove(outfile_gau)

#-------------------------------------------------------------------------------
def cdxml_visual_obabel(cdxml_file:str) -> None:
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
    call([exe_path, gjf_file])

#-------------------------------------------------------------------------------
def cdxml_visual_rdkit(cdxml_file: str, output_format='xyz', optimize=True) \
                                                                      -> None:
  '''
  Converts an input .cdxml file to 3D conformation files.
  Automatically adds explicit hydrogens, optimizes molecular geometry, 
  and strictly preserves chirality and E/Z stereochemistry from the 2D drawing.
  --
  :cdxml_file: Path to the input .cdxml file.\n
  :output_format: Desired output format ('sdf','mol',or'xyz').\n
  :optimize: Whether to perform force field optimization.\n
  returns: None.\n
  '''
  from rdkit import Chem
  from rdkit.Chem import AllChem
  try:
    open(cdxml_file, 'r')
  except:
    print(f"Error: File {cdxml_file} can not be opened.")
    exit(1)
  try:
    mols = Chem.MolsFromCDXMLFile(cdxml_file)
  except AttributeError:
    raise RuntimeError(
      "Error: Current RDKit version does not support parsing CDXML directly.\n"
      "Please upgrade using: pip install --upgrade rdkit"
    )
  if not mols:
    raise ValueError(f"Error: No valid molecules parsed from {cdxml_file}.")
  base_name = path.splitext(path.basename(cdxml_file))[0]
  if len(mols) > 1:
    print(f"Warning: Multiple molecules found in {cdxml_file}."+\
          "Processing each separately!")
    exit(1)
  if mols is None:
    print(f"Warning: No molecule found in {cdxml_file}!")
  mol_h = Chem.AddHs(mols[0])
  # ETKDGv3 utilizes distance geometry and experimental torsion-angle
  # knowledge, it strictly penalizes inversions of predefined chiral
  # centers and E/Z double bonds.
  params = AllChem.ETKDGv3()
  # Improves conformational accuracy for small ring systems
  params.useSmallRingTorsions = True
  params.randomSeed = 42
  embed_status = AllChem.EmbedMolecule(mol_h, params)
  if embed_status != 0:
    print(f"Warning: Failed to generate 3D coordinates"\
          +" Check for severe steric clashes.")
    exit(1)
  # Force field optimization
  if optimize:
    try:
      # Prioritize MMFF94, which is highly accurate for organic molecules
      if AllChem.MMFFHasAllMoleculeParams(mol_h):
        AllChem.MMFFOptimizeMolecule(mol_h,maxIters=1000,nonBondedThresh=100.)
      else:
        # Fallback to UFF if MMFF parameters are missing 
        # (e.g., due to metals/heavy heteroatoms)
        print(f"Molecule lacks MMFF parameters."\
              +"Falling back to UFF optimization.")
        if AllChem.UFFHasAllMoleculeParams(mol_h):
          AllChem.UFFOptimizeMolecule(mol_h, maxIters=1000)
        else:
          print(f"Molecule lacks UFF parameters."+\
                "Skipping optimization.")
    except Exception as e:
      print(f"Warning: Optimization failed for molecule: {e}")
  # export 3D structure
  out_name = f"{base_name}_3D.{output_format.lower()}"
  try:
    if output_format.lower() == 'sdf':
      writer = Chem.SDWriter(out_name)
      writer.write(mol_h)
      writer.close()
    elif output_format.lower() == 'mol':
      Chem.MolToMolFile(mol_h, out_name)
    elif output_format.lower() == 'xyz':
      Chem.MolToXYZFile(mol_h, out_name)
    else:
      raise ValueError(f"Unsupported output format: {output_format}")
    print(f"Molecule processed and saved to -> {out_name}")
  except Exception as e:
    print(f"Failed to save molecule: {e}")
  if output_format.lower() == 'xyz':
    single_xyz_edit(out_name)
  elif output_format.lower() in ['sdf', 'mol']:
    call([exe_path, out_name])

#===============================================================================
if __name__ == "__main__":
  import sys
  if len(sys.argv) == 1:
    call([exe_path])
  elif len(sys.argv) != 2:
    print("Usage: gview input.[xyz|inp|trj|gjf|log|out|cif|fch|mol]")
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
    elif input_file.endswith('.inp'):
      visual_orca_inp(input_file)
    elif input_file.endswith('.gjf'):
      call([exe_path, input_file])
    elif input_file.endswith('.log'):
      call([exe_path, input_file])
    elif input_file.endswith('.out'):
      OUT_visual(input_file)
    elif input_file.endswith('.cif'):
      call([exe_path, input_file])
    elif input_file.endswith('.fch'):
      call([exe_path, input_file])
    elif input_file.endswith('.molden') or input_file.endswith('.molden.input'):
      MOLDEN_freq_visual(input_file)
    elif input_file.endswith('.mol') or input_file.endswith('.mol2'):
      call([exe_path, input_file])
    elif input_file.endswith('.pdb'):
      call([exe_path, input_file])
    elif input_file.endswith('.pdb1'):
      call([exe_path, input_file])
    elif input_file.endswith('.ml2'):
      call([exe_path, input_file])
    elif input_file.endswith('.gmmx'):
      call([exe_path, input_file])
    elif input_file.endswith('.com'):
      call([exe_path, input_file])
    elif input_file.endswith('.cub') or input_file.endswith('.cube'):
      print("It's NOT recommended to use GView to open cube files; using VMD"+\
            " or vis2c for visualization instead.")
      if input_file.endswith('.cube'):
        exit(1)
      call([exe_path, input_file])
    elif input_file.endswith('.cdxml'):
      cdxml_visual_rdkit(input_file, output_format='xyz', optimize=True)
      # cdxml_visual_obabel(input_file)  # alternative method using Open Babel
    elif input_file.endswith('.sdf'):
      call([exe_path, input_file])
    else:
      print("unrecognize input file type")
      exit(1)
