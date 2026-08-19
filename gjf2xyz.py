#!/home/lky/miniconda3/envs/polygview/bin/python
'''
Convert .gjf to .xyz
Author: Dirac4pi
env:polygview
'''

#-------------------------------------------------------------------------------
def isfloat(string:str) -> bool:
  '''
  Determine whether the string is a float or not
  --
  :string: input string\n
  return: true/false
  '''
  try:
    float(string)
  except:
    return False
  else:
    return True

#-------------------------------------------------------------------------------
def iscood(line:str) -> bool:
  '''
  Determine whether the string is a coordinate or not
  --
  :line: input string\n
  return: true/false
  '''
  s = line.split()
  if len(s) != 4:
    return False
  elif s[0][0] == '!':
    return False
  elif not s[0][0].isalpha():
    return False
  elif isfloat(s[1]) and isfloat(s[2]) and isfloat(s[3]):
    return True
  else:
    return False

#-------------------------------------------------------------------------------
def gjf2xyz(gjf:str) -> None:
  '''
  Convert gjf files to xyz siles
  --
  :gjf: gjf file to be converted.
  '''
  with open(gjf, 'r', encoding='utf-8') as rgjf:
    xyz = gjf.removesuffix('.gjf') + '.xyz'
    wxyz = open(xyz, 'w', encoding='utf-8')
    cood = []
    lattice = []
    natom = 0
    for line in rgjf:
      if line.startswith('!'):
        continue
      elif iscood(line):
        if line.find('Tv') != -1:   # lattice infomation in extended xyz
          s = line.split()
          lattice.append(s[1])
          lattice.append(s[2])
          lattice.append(s[3])
        else:
          i = line.find('(')
          if i != -1:
            j = line.find(')')
            line = line[0:i]+line[j+1:]
          cood.append(line)
          natom += 1
    if lattice == []:
      wxyz.write(str(natom)+'\n')
      wxyz.write('\n')
      wxyz.writelines(cood)
      wxyz.write('\n')
      wxyz.close()
    else:
      wxyz.write(str(natom)+'\n')
      wxyz.write('Lattice=\"')
      for i in range(8):
        wxyz.write(lattice[i]+' ')
      wxyz.write(lattice[8])
      wxyz.write('\" Properties=species:S:1:pos:R:3\n')
      wxyz.writelines(cood)
      wxyz.write('\n')
      wxyz.close()

#===============================================================================
if __name__ == "__main__":
  import sys
  arguments = sys.argv[1:]
  if not arguments:
    print("NO arguments provided!")
    print("usage: gjf2xyz.py file1.gjf file2.gjf file3.gjf ...")
    sys.exit(1)
  for arg in arguments:
    gjf2xyz(arg)
