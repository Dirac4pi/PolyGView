#!/home/lky/miniconda3/envs/polygview/bin/python
'''
Convert .gjf to ORCA .inp
Author: Dirac4pi
env:polygview
'''

from gjf2xyz import gjf2xyz
from pathlib import Path

def is_two_integers(s: str) -> bool:
  """
  Determine whether the string is 2 integers separated by a single space or not
  --
  :line: input string\n
  return: true/false
  """
  parts = s.split()
  if len(parts) != 2:
    return False
  try:
    int(parts[0])
    int(parts[1])
  except ValueError:
    return False
  else:
    return True


#-------------------------------------------------------------------------------
def gjf2inp(gjf:str) -> None:
  '''
  Convert gjf files to ORCA inp files
  --
  :gjf: gjf file to be converted.
  '''
  csfound = False
  with open(gjf, 'r', encoding='utf-8') as rgjf:
    for line in rgjf:
      s = line.split()
      if is_two_integers(line):
        charge = int(s[0])
        spin = int(s[1])
        csfound = True
        break
  if not csfound:
    print(f"Charge and spin multiplicity not found in {gjf}, set 0 1")
    charge = 0
    spin = 1
  inp = gjf.removesuffix('.gjf') + '.inp'
  xyz = gjf.removesuffix('.gjf') + '.xyz'
  gjf2xyz(gjf)
  print(f"successfully generate {xyz}")
  # this is the default keywords for ORCA inp files, change it as you like
  content = f'''! wB97X-D3 def2-SVP OPT RIJCOSX CPCM(Water) def2/J
! TightSCF NoAutoStart MiniPrint NoPop
%MAXCORE 3000
%PAL NPROCS 64 END
%geom
  maxstep 0.03
  maxiter 100
  Scan
    B 0 8 = 1.8, 3.0, 10
  end
end
* XYZFILE {charge} {spin} {Path(xyz).resolve()}
'''
  try:
    with open(inp, 'w', encoding='utf-8') as file:
      file.write(content)
    print(f"successfully generate {inp}")
  except Exception as e:
    print(f"error: {e}")


#===============================================================================
if __name__ == "__main__":
  import sys
  arguments = sys.argv[1:]
  if not arguments:
    print("NO arguments provided!")
    print("usage: gjf2inp.py 1.gjf 2.gjf 3.gjf ...")
    sys.exit(1)
  for arg in arguments:
    gjf2inp(arg)
