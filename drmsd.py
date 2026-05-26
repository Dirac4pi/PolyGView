'''
Provide the Distance RMSD for two different conformations of the same molecule.
Author: Dirac4pi
env:base
'''

import sys
import numpy as np

#-------------------------------------------------------------------------------
def read_xyz(filepath:str) -> tuple:
  """
  Reads an XYZ file and returns atomic symbols and coordinates.
  --
  :xyz_file: input file name\n
  return: tuple of atomic symbols and coordinates
  """
  with open(filepath, 'r') as f:
    lines = f.readlines()
  if not lines:
    raise ValueError(f"File {filepath} is empty.")
  # Get the number of atoms from the first line
  num_atoms = int(lines[0].strip())
  coords = np.zeros((num_atoms, 3), dtype=float)
  symbols = []
  # Parse coordinates starting from the third line
  for i in range(num_atoms):
    parts = lines[i + 2].split()
    symbols.append(parts[0])
    coords[i] = [float(parts[1]), float(parts[2]), float(parts[3])]
  return symbols, coords

#-------------------------------------------------------------------------------
def drmsd(coords_a:np.ndarray, coords_b:np.ndarray) -> float:
  """
  Calculates the Distance Root Mean Square Deviation (dRMSD)
  between two sets of coordinates using purely numpy operations.
  --
  :coords_a: first set of coordinates\n
  :coords_b: second set of coordinates\n
  return: dRMSD value
  """
  num_atoms = coords_a.shape[0]
  if num_atoms < 2:
    return 0.0
  # Compute pairwise differences using numpy broadcasting
  # coords shape: (N, 3) -> expanded diff shape: (N, N, 3)
  dist_a = coords_a[:, np.newaxis, :] - coords_a[np.newaxis, :, :]
  dist_matrix_a = np.sqrt(np.sum(dist_a ** 2, axis=-1))
  dist_b = coords_b[:, np.newaxis, :] - coords_b[np.newaxis, :, :]
  dist_matrix_b = np.sqrt(np.sum(dist_b ** 2, axis=-1))
  # Extract only the upper triangular part (i < j)
  # This avoids double-counting pairs and excludes the zero-diagonal
  upper_tri_idx = np.triu_indices(num_atoms, k=1)
  pdist_a = dist_matrix_a[upper_tri_idx]
  pdist_b = dist_matrix_b[upper_tri_idx]
  # Calculate the final dRMSD
  #! Main diagonal does not contribute to the dRMSD, so it's N(N-1)/2
  drmsd = np.sqrt(np.mean((pdist_a - pdist_b) ** 2))
  return drmsd

#-------------------------------------------------------------------------------
def dmaxd(coords_a:np.ndarray, coords_b:np.ndarray) -> float:
  """
  Calculates the Distance Maximum Deviation (dmaxd)
  between two sets of coordinates using purely numpy operations.
  --
  :coords_a: first set of coordinates\n
  :coords_b: second set of coordinates\n
  return: dmaxd value
  """
  num_atoms = coords_a.shape[0]
  if num_atoms < 2:
    return 0.0
  # Compute pairwise differences using numpy broadcasting
  dist_a = coords_a[:, np.newaxis, :] - coords_a[np.newaxis, :, :]
  dist_matrix_a = np.sqrt(np.sum(dist_a ** 2, axis=-1))
  dist_b = coords_b[:, np.newaxis, :] - coords_b[np.newaxis, :, :]
  dist_matrix_b = np.sqrt(np.sum(dist_b ** 2, axis=-1))
  # Extract only the upper triangular part (i < j)
  upper_tri_idx = np.triu_indices(num_atoms, k=1)
  pdist_a = dist_matrix_a[upper_tri_idx]
  pdist_b = dist_matrix_b[upper_tri_idx]
  # Calculate the final dmaxd
  dmaxd = np.max(np.abs(pdist_a - pdist_b))
  return dmaxd

#===============================================================================
if __name__ == "__main__":
  # Require exactly two arguments (the two xyz files)
  if len(sys.argv) != 3:
    print("Usage: python drmsd.py <conf_A.xyz> <conf_B.xyz>")
    sys.exit(1)
  xyz_file_a = sys.argv[1]
  xyz_file_b = sys.argv[2]
  try:
    sym_a, coords_a = read_xyz(xyz_file_a)
    sym_b, coords_b = read_xyz(xyz_file_b)
    # Verify atom counts match between the two conformations
    if len(sym_a) != len(sym_b):
      print("Error: Atom counts differ between the two files!")
      sys.exit(1)
    drmsd_val = drmsd(coords_a, coords_b)
    dmaxd_val = dmaxd(coords_a, coords_b)
    # Print the result to the terminal
    print(f"Distance MAXD between '{xyz_file_a}' and '{xyz_file_b}': {dmaxd_val:.6f} Angstroms")
    print(f"Distance RMSD between '{xyz_file_a}' and '{xyz_file_b}': {drmsd_val:.6f} Angstroms")
  except Exception as err:
    print(f"Error processing files: {err}")
    sys.exit(1)