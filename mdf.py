'''
Calculate the Mean Distance Fluctuation (MDF) for each atom in a given
multi-frame XYZ file.
Author: Dirac4pi
env:base
'''

import numpy as np
import os
import sys

#-------------------------------------------------------------------------------
def calc_mdf(xyzfile: str) -> np.ndarray:
  """
  Calculate the Mean Distance Fluctuation (MDF) for each atom in a given
  multi-frame XYZ file.
  --
  :xyz_file: input file name\n
  return: tuple of atomic symbols and coordinates
  Returns: MDF values for the N atoms, T frames.
  """
  if not os.path.exists(xyzfile):
    raise FileNotFoundError(f"File not found: {xyzfile}")
  sum_d = None
  sum_d2 = None
  T = 0
  N = None
  with open(xyzfile, 'r') as f:
    while True:
      line = f.readline()
      if not line: # end of file
        break
      line = line.strip()
      if not line: # empty lines
        continue
      try:
        current_n = int(line)
      except ValueError:
        raise ValueError(f"Unable to read atom count at frame {T+1}.")
      if N is None: # first frame
        N = current_n
        if N <= 1:
          return np.zeros(N)  # A single atom has no distance fluctuation
        sum_d = np.zeros((N, N), dtype=np.float64)
        sum_d2 = np.zeros((N, N), dtype=np.float64)
      elif current_n != N:
        raise ValueError(\
      f"Atom count mismatch: expected {N}, but got {current_n} at frame {T+1}.")
      f.readline()   # empty line or comment line
      coords = np.zeros((N, 3), dtype=np.float64)
      for i in range(N):
        parts = f.readline().split()
        # Standard XYZ format: Atom_symbol X Y Z
        coords[i, 0] = float(parts[1])
        coords[i, 1] = float(parts[2])
        coords[i, 2] = float(parts[3])
      # Calculate the N x N distance matrix for the current frame
      # diff shape: (N, N, 3)
      diff = coords[:, np.newaxis, :] - coords[np.newaxis, :, :]
      dist_matrix = np.sqrt(np.sum(diff**2, axis=-1))
      sum_d += dist_matrix
      sum_d2 += dist_matrix**2
      T += 1
  if T == 0:
    raise ValueError("No valid trajectory frames were read from the file.")
  mean_d = sum_d / T
  mean_d2 = sum_d2 / T
  variance = np.maximum(mean_d2 - mean_d**2, 0.0) # sqrt need non-negative value
  DF = np.sqrt(variance)
  MDF = np.sum(DF, axis=1) / (N - 1)  # N-1 non-zero distances per atom
  return T, N, MDF

#===============================================================================
if __name__ == "__main__":
  if len(sys.argv) == 2:
    input_file = sys.argv[1]
    T, N, vmdf = calc_mdf(input_file)
    print(f"Processed {T} frames with {N} atoms.")
    with open(input_file, 'r') as f:
      line = f.readline()
      line = f.readline()
      print(f"{'Index':>5s}  {'Symbol':>6s}  {'MDF(Angstrom)':>10s}")
      for i in range(N):
        parts = f.readline().split()
        print(f"{i+1:3d}      {parts[0]:<2.2}      {vmdf[i]:.6f}")
  else:
    exit('Usage: python mdf.py traj.xyz')