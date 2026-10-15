#!/home/lky/miniconda3/envs/polygview/bin/python
"""Analyze one fixed conformer or compare two conformers without molecular SCF.

Two XYZ files, paired 1-based selections, English reports and optional CSV/JSON.
Atomic density preparation is a separate, explicitly requested atomic HF step.
"""
import argparse
import csv
import itertools
import json
import math
import os
from pathlib import Path
import sys
import time
import numpy as np
from scipy.linalg import eigh, solve
from scipy.integrate import trapezoid
from scipy.interpolate import CubicHermiteSpline

BOHR = 0.529177210903
KCAL = 627.5094740631
VERSION = '0.3.0'
THOLE_SOURCE = 'van Duijnen and Swart, J. Phys. Chem. A 1998, 102, 2399; DOI 10.1021/jp980221f; Table 7 empirical linear model, Eq. 9c'
THOLE = dict(H=3.5020,C=10.1756,N=7.6048,O=6.3940,F=2.9413,S=19.7422,Cl=16.1145,Br=22.6671,I=34.2434)
ELEMENTS = ('X H He Li Be B C N O F Ne Na Mg Al Si P S Cl Ar K Ca Sc Ti V Cr Mn Fe Co Ni Cu Zn Ga Ge As Se Br Kr Rb Sr Y Zr Nb Mo Tc Ru Rh Pd Ag Cd In Sn Sb Te I Xe Cs Ba La Ce Pr Nd Pm Sm Eu Gd Tb Dy Ho Er Tm Yb Lu Hf Ta W Re Os Ir Pt Au Hg Tl Pb Bi Po At Rn').split()


def read_xyz(path):
  lines = Path(path).read_text().splitlines()
  if len(lines)<2:
    raise ValueError('XYZ requires an atom count and a comment line')
  n = int(lines[0].strip())
  if n<1 or len(lines)<n+2:
    raise ValueError('Invalid XYZ atom count')
  symbols, coords = [], []
  for line in lines[2:n+2]:
    row = line.split()
    if len(row)<4:
      raise ValueError('Malformed XYZ atom line')
    symbol = row[0].capitalize()
    if symbol not in ELEMENTS[1:]:
      raise ValueError('Unsupported element: '+symbol)
    symbols.append(symbol)
    coords.append([float(x) for x in row[1:4]])
  xyz = np.asarray(coords)
  if not np.isfinite(xyz).all():
    raise ValueError('Nonfinite coordinates')
  if n>1:
    d = np.linalg.norm(xyz[:,None]-xyz[None,:],axis=2)
    if np.min(d[np.triu_indices(n,1)])<1e-5:
      raise ValueError('Coincident nuclei')
  return symbols, xyz


def indices(text,n):
  if text.strip().lower() in ('c','all'):
    return list(range(n))
  out=[]
  for token in text.split(','):
    parts=token.strip().split('-')
    if len(parts)==1:
      out.append(int(parts[0])-1)
    elif len(parts)==2:
      a,b=map(int,parts)
      if b<a:
        raise ValueError('Descending selection range')
      out.extend(range(a-1,b))
    else:
      raise ValueError('Invalid atom selection')
  if not out or len(set(out))!=len(out) or min(out)<0 or max(out)>=n:
    raise ValueError('Empty, duplicate or out-of-range atom selection')
  return out


def align(a,b,sel=None):
  sel=np.arange(len(a)) if sel is None else np.asarray(sel)
  ac,bc=a[sel].mean(0),b[sel].mean(0)
  u,s,vt=np.linalg.svd((b[sel]-bc).T@(a[sel]-ac))
  rot=u@np.diag([1.,1.,np.linalg.det(u@vt)])@vt
  aligned=(b-bc)@rot+ac
  rms=float(np.sqrt(np.mean(np.sum((aligned[sel]-a[sel])**2,axis=1))))
  rank=int(np.sum(s>max(s[0],1.)*1e-10))
  return aligned,dict(rotation=rot,translation=ac-bc@rot,rmsd_A=rms,alignment_rank=rank,
    warning='Axis-dependent comparisons are ambiguous for a collinear/single-atom alignment.' if rank<2 else '')


def fragments(path,n,ref_ids):
  if path is None:
    return {f'atom_{ref_ids[i]+1}':[i] for i in range(n)}
  raw=json.loads(Path(path).read_text())
  if not isinstance(raw,dict) or not raw:
    raise ValueError('Fragments must be a nonempty JSON object')
  lookup={j+1:i for i,j in enumerate(ref_ids)}
  out={}
  used=set()
  for name,items in raw.items():
    if not isinstance(items,list) or not items:
      raise ValueError('Fragment entries must be nonempty atom-index lists')
    if any(type(x) is not int or x not in lookup for x in items):
      raise ValueError('Fragments must use selected original REF atom numbers')
    group=[lookup[x] for x in items]
    if len(set(group))!=len(group) or used.intersection(group):
      raise ValueError('Fragments must be disjoint without duplicate atoms')
    used.update(group)
    out[str(name)]=group
  rest=sorted(set(range(n))-used)
  if rest:
    if 'unassigned' in out:
      raise ValueError('Reserved fragment name: unassigned')
    out['unassigned']=rest
  return out


def build_integrals(symbols,xyz,basis,ecp,max_ao):
  from pyscf import gto
  import pyscf
  z=sum(ELEMENTS.index(s) for s in symbols)
  mol=gto.M(atom=list(zip(symbols,xyz)),unit='Angstrom',basis=basis,
    ecp=ecp or {},charge=0,spin=z%2,cart=False,verbose=0)
  if mol.nao_nr()>max_ao:
    raise ValueError(f'{mol.nao_nr()} AOs exceed --max-ao {max_ao}')
  s=mol.intor('int1e_ovlp')
  t=mol.intor('int1e_kin')
  d=1/np.sqrt(np.diag(s))
  s=d[:,None]*s*d[None,:]
  t=d[:,None]*t*d[None,:]
  loc=mol.ao_loc_nr()
  meta=[]
  for shell in range(mol.nbas):
    atom=mol.bas_atom(shell)
    l=mol.bas_angular(shell)
    ex=mol.bas_exp(shell)
    meta.append(dict(atom=atom,shell=shell,l=l,start=int(loc[shell]),stop=int(loc[shell+1]),
      min_exponent=float(np.min(ex)),max_exponent=float(np.max(ex))))
  return s,t,meta,dict(pyscf_version=pyscf.__version__,basis=basis,ecp=ecp,nao=mol.nao_nr(),
    electron_state='No molecular electronic state is solved; neutral/parity spin used only to initialize Mole.')


def health(s,t,meta,thresholds):
  vals,vec=eigh(s)
  if vals[0]<-1e-10:
    raise ValueError('Overlap has a significantly negative eigenvalue')
  rows=[]
  kinetic=[]
  for tau in thresholds:
    low=vals<tau
    weights=np.sum(vec[:,low]**2,axis=1)
    for sh in meta:
      weight=float(weights[sh['start']:sh['stop']].sum())
      if weight>1e-12:
        rows.append(dict(threshold=tau,weight=weight,**sh))
    keep=vals>tau
    if np.any(keep):
      x=vec[:,keep]/np.sqrt(vals[keep])
      eig=eigh(x.T@t@x,eigvals_only=True)
      kinetic.append(dict(threshold=tau,rank=int(keep.sum()),max_Ha=float(eig[-1]),
        p95_Ha=float(np.quantile(eig,.95)),min_Ha=float(eig[0])))
  return dict(status='OK',minimum=float(vals[0]),maximum=float(vals[-1]),
    condition=float(vals[-1]/vals[0]) if vals[0]>0 else None,
    low_counts={str(tau):int(np.sum(vals<tau)) for tau in thresholds},
    shell_weights=rows,kinetic=kinetic),vals


def metric(sa,sb,tau,project=False):
  ea,ua=eigh(sa)
  eb=eigh(sb,eigvals_only=True)
  n=len(ea)
  if min(ea[0],eb[0])<=tau:
    if not project:
      return dict(status='ILL_CONDITIONED',conclusion='Full-space logarithmic metric withheld: overlap below the selected threshold. Use --metric-project for a common retained-space diagnostic.')
    ev,u=eigh((sa+sb)*.5)
    q=u[:,ev>tau]
    if q.shape[1]==0:
      raise ValueError('No common retained metric space')
    aa,bb=q.T@sa@q,q.T@sb@q
    if min(eigh(aa,eigvals_only=True)[0],eigh(bb,eigvals_only=True)[0])<=tau:
      return dict(status='ILL_CONDITIONED',conclusion='Common projection still fails the overlap threshold; metric withheld.')
    projected=True
  else:
    aa,bb=sa,sb
    projected=False
  gamma=eigh(bb,aa,eigvals_only=True)
  logs=np.log(gamma)
  eta=float(np.max(np.abs(gamma-1)))
  rms=float(np.linalg.norm(logs)/np.sqrt(len(logs)))
  if np.max(np.abs(logs))<1e-7:
    text='No resolved AO-metric change in the common coordinate frame.'
  else:
    text=f'With unchanged AO coefficients, squared norms can change by factors {gamma[0]:.4g} to {gamma[-1]:.4g}; this is not a prediction of SCF failure.'
  return dict(status='OK',projected=projected,dimension=len(gamma),original_dimension=n,
    gamma=gamma,log_gamma=logs,log_distance=float(np.linalg.norm(logs)),log_rms=rms,
    max_abs_log=float(np.max(np.abs(logs))),norm_w_minus_i=eta,conclusion=text)


def fingerprint(sa,sb,meta,a,b,centers,radius):
  ea=eigh(sa,eigvals_only=True)
  eb=eigh(sb,eigvals_only=True)
  rows=[]
  for atom in centers:
    da=np.linalg.norm(a-a[atom],axis=1)
    db=np.linalg.norm(b-b[atom],axis=1)
    neighbors=np.where(np.minimum(da,db)<=radius)[0]
    ids=[i for sh in meta if sh['atom'] in neighbors for i in range(sh['start'],sh['stop'])]
    la=eigh(sa[np.ix_(ids,ids)],eigvals_only=True)
    lb=eigh(sb[np.ix_(ids,ids)],eigvals_only=True)
    rows.append(dict(center=atom,neighbors=neighbors.tolist(),nao=len(ids),
      distance=float(np.linalg.norm(lb-la)),rms=float(np.linalg.norm(lb-la)/np.sqrt(len(ids)))))
  rows.sort(key=lambda r:r['rms'],reverse=True)
  delta=eb-ea
  text='No resolved overlap-spectrum difference.' if np.linalg.norm(delta)<1e-8 else 'Overlap environments differ; local rankings locate structural changes, not electronic-state changes.'
  return dict(status='OK',distance=float(np.linalg.norm(delta)),rms=float(np.linalg.norm(delta)/np.sqrt(len(delta))),
    spectrum_ref=ea,spectrum_probe=eb,local=rows,conclusion=text)


def pair_records(pa,pb,frags,ref_ids,probe_ids):
  n=len(pa)
  rows=[]
  for i in range(n):
    for j in range(i+1,n):
      rows.append(dict(ref_atoms=f'{ref_ids[i]+1}-{ref_ids[j]+1}',probe_atoms=f'{probe_ids[i]+1}-{probe_ids[j]+1}',
        ref_kcal=float(pa[i,j]+pa[j,i])*KCAL,probe_kcal=float(pb[i,j]+pb[j,i])*KCAL,
        delta_kcal=float(pb[i,j]+pb[j,i]-pa[i,j]-pa[j,i])*KCAL))
  grouped=[]
  for f,g in itertools.combinations_with_replacement(frags,2):
    ii,jj=frags[f],frags[g]
    if f==g:
      va=float(pa[np.ix_(ii,ii)].sum())
      vb=float(pb[np.ix_(ii,ii)].sum())
    else:
      va=float(pa[np.ix_(ii,jj)].sum()+pa[np.ix_(jj,ii)].sum())
      vb=float(pb[np.ix_(ii,jj)].sum()+pb[np.ix_(jj,ii)].sum())
    grouped.append(dict(fragment_A=f,fragment_B=g,ref_kcal=va*KCAL,probe_kcal=vb*KCAL,delta_kcal=(vb-va)*KCAL))
  return sorted(rows,key=lambda r:abs(r['delta_kcal']),reverse=True),sorted(grouped,key=lambda r:abs(r['delta_kcal']),reverse=True)


def dispersion(symbols,a,b,charge,method,frags,ref_ids,probe_ids):
  from dftd4.interface import DispersionModel,DampingParam
  from importlib.metadata import version
  numbers=np.asarray([ELEMENTS.index(s) for s in symbols],dtype=np.int32)
  param=DampingParam(method=method)
  data=[]
  for coords in (a,b):
    model=DispersionModel(numbers,coords/BOHR,charge=float(charge))
    energy=float(np.asarray(model.get_dispersion(param,grad=False)['energy']))
    pair=model.get_pairwise_dispersion(param)
    additive=np.asarray(pair['additive pairwise energy'])
    nonadd=np.asarray(pair['non-additive pairwise energy'])
    if additive.shape!=(len(a),len(a)) or nonadd.shape!=additive.shape:
      raise ValueError('Unexpected D4 pair matrix shape')
    total=float(additive.sum()+nonadd.sum())
    if abs(total-energy)>1e-8:
      raise ValueError('D4 pairwise sum does not match total energy; backend convention not supported')
    data.append((energy,additive,nonadd))
  rows,groups=pair_records(data[0][1],data[1][1],frags,ref_ids,probe_ids)
  delta=(data[1][0]-data[0][0])*KCAL
  if abs(delta)<1e-7:
    text='No resolved D4 dispersion-energy difference.'
  else:
    text=f'{"PROBE" if delta<0 else "REF"} has more favorable model dispersion by {abs(delta):.4f} kcal/mol; this is not a total conformer-energy ranking.'
  return dict(status='OK',method=method,charge=charge,version=version('dftd4'),ref_kcal=data[0][0]*KCAL,
    probe_kcal=data[1][0]*KCAL,delta_kcal=delta,pairs=rows,fragments=groups,
    additive_delta_kcal=float((data[1][1].sum()-data[0][1].sum())*KCAL),
    three_body_delta_kcal=float((data[1][2].sum()-data[0][2].sum())*KCAL),conclusion=text)


def polar_params(path,symbols,ref_ids):
  if path:
    p=json.loads(Path(path).read_text())
    if p.get('model')!='thole_linear' or p.get('units')!='bohr^3' or not p.get('source'):
      raise ValueError('Polar JSON requires model=thole_linear, units=bohr^3 and a source')
    table=p.get('elements',{})
    atom=p.get('atoms',{})
    alpha=[atom.get(str(i+1),table.get(s)) for i,s in zip(ref_ids,symbols)]
    if any(x is None for x in alpha):
      raise ValueError('Polar parameters do not cover every selected atom')
    damp=float(p['a'])
    source=p['source']
  else:
    missing=sorted(set(symbols)-set(THOLE))
    if missing:
      raise ValueError('Published default Thole set has no parameters for '+','.join(missing)+'; supply --polar-params, no metal fallback is used')
    alpha=[THOLE[s] for s in symbols]
    damp=1.7278
    source=THOLE_SOURCE
  alpha=np.asarray(alpha,dtype=float)
  if not np.isfinite(alpha).all() or np.any(alpha<=0) or not np.isfinite(damp) or damp<=0:
    raise ValueError('Polarizabilities and damping must be positive finite numbers')
  return alpha,damp,source


def polar_tensor(xyz,alpha,damping):
  r=xyz/BOHR
  n=len(r)
  relay=np.diag(np.repeat(1/alpha,3))
  for i in range(n):
    for j in range(i):
      v=r[i]-r[j]
      dist=np.linalg.norm(v)
      if dist<1e-8:
        raise ValueError('Coincident polarizable sites')
      u=dist/(damping*(alpha[i]*alpha[j])**(1/6))
      f3,f5=(4*u**3-3*u**4,u**4) if u<1 else (1.,1.)
      block=np.eye(3)*f3/dist**3-3*f5*np.outer(v,v)/dist**5
      relay[3*i:3*i+3,3*j:3*j+3]=block
      relay[3*j:3*j+3,3*i:3*i+3]=block.T
  eig=eigh(relay,eigvals_only=True)
  if eig[0]<=1e-10*eig[-1]:
    raise ValueError('Polarization response matrix is unstable or nearly singular')
  field=np.tile(np.eye(3),(n,1))
  response=solve(relay,field,assume_a='pos').reshape(n,3,3)
  tensor=response.sum(axis=0)
  tensor=(tensor+tensor.T)/2
  return tensor,response,float(eig[-1]/eig[0])


def polarization(symbols,a,b,frags,params,ref_ids,axis):
  al,damp,source=polar_params(params,symbols,ref_ids)
  tensors=[]
  responses=[]
  conditions=[]
  for xyz in (a,b):
    t,r,c=polar_tensor(xyz,al,damp)
    tensors.append(t)
    responses.append(r)
    conditions.append(c)
  def summary(t):
    vals,vec=eigh(t)
    anis=float(np.sqrt(((vals[0]-vals[1])**2+(vals[1]-vals[2])**2+(vals[2]-vals[0])**2)/2))
    return dict(mean_A3=float(np.trace(t)/3*BOHR**3),anisotropy_A3=anis*BOHR**3,
      eigenvalues_A3=vals*BOHR**3,eigenvectors=vec,tensor_A3=t*BOHR**3)
  aa,bb=map(summary,tensors)
  rows=[]
  for f,ids in frags.items():
    ta=responses[0][ids].sum(0)*BOHR**3
    tb=responses[1][ids].sum(0)*BOHR**3
    rows.append(dict(fragment=f,delta_trace_mean_A3=float(np.trace(tb-ta)/3),
      response_difference_F_A3=float(np.linalg.norm(tb-ta))))
  rows.sort(key=lambda r:r['response_difference_F_A3'],reverse=True)
  out=dict(status='OK',source=source,damping=damp,ref=aa,probe=bb,conditions=conditions,fragments=rows,
    note='Effective isotropic atomic parameters; frozen classical response, not spin-resolved quantum polarizability.')
  norm=float(np.linalg.norm(tensors[1]-tensors[0])*BOHR**3)
  delta=bb['mean_A3']-aa['mean_A3']
  out['conclusion']='No resolved static-response difference.' if norm<1e-7 else f'Model mean response changes by {delta:+.4f} A^3; tensor-change norm is {norm:.4f} A^3. No Raman intensity or reactivity is inferred.'
  if axis is not None:
    i,j=axis
    if np.linalg.norm(a[j]-a[i])<1e-8:
      raise ValueError('Zero response axis')
    direction=(a[j]-a[i])/np.linalg.norm(a[j]-a[i])
    out['common_axis_ref_atoms']=[ref_ids[i]+1,ref_ids[j]+1]
    out['axis_ref_A3']=float(direction@tensors[0]@direction*BOHR**3)
    out['axis_probe_A3']=float(direction@tensors[1]@direction*BOHR**3)
  return out


def prepare_density(elements,basis,directory,overwrite=False):
  from pyscf import gto
  from pyscf.scf import atom_hf
  from pyscf.dft import numint
  import pyscf
  root=Path(directory)
  root.mkdir(parents=True,exist_ok=True)
  for symbol in elements:
    z=ELEMENTS.index(symbol)
    path=root/(symbol+'.json')
    if path.exists() and not overwrite:
      raise FileExistsError(str(path))
    mol=gto.M(atom=f'{symbol} 0 0 0',basis=basis,spin=z%2,verbose=0,cart=False)
    if mol.nelectron == 1:
      mf = atom_hf.AtomHF1e(mol)
    else:
      mf = atom_hf.AtomSphAverageRHF(mol)
    mf.conv_tol=1e-10
    mf.kernel()
    if not mf.converged:
      raise ValueError('Atomic reference did not converge: '+symbol)
    dm = np.asarray(mf.make_rdm1())
    # Convert spin-resolved density matrices to the total density.
    if dm.ndim == 3 and dm.shape[0] == 2:
      dm = dm.sum(axis=0)
    nao = mol.nao_nr()
    if dm.ndim != 2 or dm.shape != (nao, nao):
      raise ValueError(
        f'{symbol}: unexpected density-matrix shape {dm.shape}; '
        f'expected ({nao}, {nao})'
      )
    if not np.isfinite(dm).all():
      raise ValueError(f'{symbol}: nonfinite atomic density matrix')
    rad=np.concatenate(([0.],np.geomspace(1e-6,40.,2400)))
    points=np.zeros((len(rad),3))
    points[:,0]=rad
    ao=numint.eval_ao(mol,points,deriv=1)
    rho=np.einsum('pi,ij,pj->p',ao[0],dm,ao[0])
    grad=2*np.einsum('pi,ij,pj->p',ao[1],dm,ao[0])
    charge=float(trapezoid(4*np.pi*rad**2*rho,x=rad))
    if abs(charge-z)>max(.02,z*.003):
      raise ValueError('Atomic density radial normalization failed: '+symbol)
    obj=dict(element=symbol,r_bohr=rad.tolist(),rho=rho.tolist(),drho=grad.tolist(),electrons=charge,
      source=f'PySCF {pyscf.__version__} spherical-average neutral atomic HF; basis={basis}; generated by quantmol',
      units='atomic',basis=basis)
    path.write_text(json.dumps(obj,indent=2))
  return root


def density_tables(directory,symbols):
  if directory is None:
    raise ValueError('IGM needs --density-library DIR; use the separate --prepare-density command once. Atomic SCF is never run implicitly.')
  tables={}
  for symbol in set(symbols):
    p=json.loads((Path(directory)/(symbol+'.json')).read_text())
    if p.get('element')!=symbol or p.get('units')!='atomic' or not p.get('source'):
      raise ValueError('Invalid density-table metadata: '+symbol)
    r=np.asarray(p['r_bohr'],float)
    rho=np.asarray(p['rho'],float)
    dr=np.asarray(p['drho'],float)
    if len(r)<4 or rho.shape!=r.shape or dr.shape!=r.shape or r[0]!=0 or np.any(np.diff(r)<=0):
      raise ValueError('Invalid radial density grid: '+symbol)
    if not all(np.isfinite(x).all() for x in (r,rho,dr)) or np.any(rho<0):
      raise ValueError('Invalid radial density values')
    if abs(rho[-1])>1e-10 or abs(dr[-1])>1e-9:
      raise ValueError('Density grid does not extend into a negligible tail')
    spline = CubicHermiteSpline(r, rho, dr, extrapolate=False)
    tables[symbol] = (
      spline,
      spline.derivative(),
      float(r[-1]),
      p['source']
    )
  return tables


def atom_gradient(points,center,table):
  vec=points-center
  r=np.linalg.norm(vec,axis=1)
  interp,deriv,end,_=table
  d=np.zeros(len(points))
  mask=r<=end
  d[mask]=deriv(r[mask])
  return vec*(d/np.maximum(r,1e-15))[:,None]


# Covalent radii in Angstrom; Cordero-style common-element values.
IGM_RADII = dict(H=.31,B=.84,C=.76,N=.71,O=.66,F=.57,Si=1.11,P=1.07,S=1.05,
  Cl=1.02,Ge=1.20,As=1.19,Se=1.20,Br=1.20,I=1.39)
IGM_NONMETALS = set('H He B C N O F Ne Si P S Cl Ar Ge As Se Br Kr Sb Te I Xe At Rn'.split())


def igm_topology(symbols,a,b,bond_file=None,ref_ids=None,scale=1.15):
  n=len(symbols)
  metal={i for i,x in enumerate(symbols) if x not in IGM_NONMETALS}
  def infer(xyz):
    edges=set()
    for i in range(n):
      for j in range(i):
        if i in metal or j in metal:
          continue
        if symbols[i] not in IGM_RADII or symbols[j] not in IGM_RADII:
          raise ValueError('No IGM covalent radius; provide --igm-bonds for this element')
        if np.linalg.norm(xyz[i]-xyz[j])<=scale*(IGM_RADII[symbols[i]]+IGM_RADII[symbols[j]]):
          edges.add((j,i))
    return edges
  if bond_file:
    raw=json.loads(Path(bond_file).read_text())
    lookup={j+1:i for i,j in enumerate(ref_ids if ref_ids is not None else range(n))}
    edges=set()
    if not isinstance(raw,list):
      raise ValueError('--igm-bonds requires a JSON list of original REF index pairs')
    for pair in raw:
      if not isinstance(pair,list) or len(pair)!=2 or any(type(x) is not int or x not in lookup for x in pair):
        raise ValueError('Invalid IGM bond pair')
      i,j=sorted(lookup[x] for x in pair)
      if i==j:
        raise ValueError('Self bond in IGM topology')
      if i not in metal and j not in metal:
        edges.add((i,j))
    ea=eb=edges
    source='Explicit common nonmetal covalent graph; metal edges omitted from exclusion paths.'
  else:
    ea,eb=infer(a),infer(b)
    source='Union of distance-inferred nonmetal graphs in REF/PROBE; metal edges omitted.'
  graph=[set() for _ in range(n)]
  for i,j in ea|eb:
    graph[i].add(j)
    graph[j].add(i)
  near=set(ea|eb)
  for i in range(n):
    for j in graph[i]:
      for k in graph[j]:
        if k!=i:
          near.add(tuple(sorted((i,k))))
  return dict(metals=metal,edges=ea|eb,near=near,
    changed_edges=[list(x) for x in sorted(ea^eb)],source=source)


def igm_box(a,b,ids,step,padding):
  centers=(a[ids] if b is None else np.vstack((a[ids],b[ids])))/BOHR
  h=step/BOHR
  lo=centers.min(0)-padding/BOHR
  shape=np.maximum(1,np.ceil((centers.max(0)+padding/BOHR-lo)/h).astype(int))
  return lo,shape,h


def igm_pair(a,b,symbols,ii,jj,tables,step,padding,thresholds,max_points,
    box=None,level=0,expand=0):
  # The refined grids partition exactly the same coarse-grid physical box.
  lo,shape,h=igm_box(a,b,ii+jj,step,padding) if box is None else box
  lo=np.array(lo,copy=True)
  shape=np.array(shape,copy=True)
  if expand:
    layers=int(np.ceil(expand/BOHR/h))
    lo-=layers*h
    shape+=2*layers
  shape*=2**level
  h=h/2**level
  size=int(np.prod(shape.astype(object)))
  if size>max_points:
    raise ValueError(f'Pair grid has {size} cells, over --igm-max-points {max_points}')
  conformers=(a,) if b is None else (a,b)
  totals=np.zeros(len(conformers))
  vols=np.zeros((len(conformers),len(thresholds)))
  edge=np.zeros(len(conformers))
  dv=h**3
  for start in range(0,size,8192):
    flat=np.arange(start,min(start+8192,size))
    ix=np.stack(np.unravel_index(flat,tuple(shape)),axis=1)
    points=lo+(ix+.5)*h
    boundary=np.any((ix==0)|(ix==shape-1),axis=1)
    for k,coords in enumerate(conformers):
      sums=[]
      for group in (ii,jj):
        grad=np.zeros_like(points)
        for i in group:
          grad+=atom_gradient(points,coords[i]/BOHR,tables[symbols[i]])
        sums.append(grad)
      dg=np.maximum(0,np.linalg.norm(sums[0],axis=1)+np.linalg.norm(sums[1],axis=1)-np.linalg.norm(sums[0]+sums[1],axis=1))
      totals[k]+=dg.sum()*dv
      vols[k]+=[np.count_nonzero(dg>tau)*dv*BOHR**3 for tau in thresholds]
      if boundary.any():
        edge[k]=max(edge[k],float(dg[boundary].max()))
  return totals,vols,dict(cells=size,step_A=h*BOHR,boundary_max_delta_g=edge.tolist(),
    box_low_bohr=lo.tolist(),box_high_bohr=(lo+shape*h).tolist())


def igm_candidates(symbols,a,b,frags,cutoff,exclude_depth,include_coordination,topology):
  candidates=[]
  filtered=[]
  metadata={}
  for f,g in itertools.combinations(frags,2):
    ii,jj=frags[f],frags[g]
    atom_pair=len(ii)==len(jj)==1
    cross=[tuple(sorted((i,j))) for i in ii for j in jj]
    metal_pair=any(i in topology['metals'] for i in ii+jj)
    boundary=any(pair in topology['edges'] for pair in cross)
    category='COORDINATION_OR_METAL' if metal_pair else 'COVALENT_BOUNDARY' if boundary else 'NONBONDED'
    reason=None
    if atom_pair:
      pair=cross[0]
      if metal_pair and not include_coordination:
        reason='Metal-site pair excluded from default nonbonded ranking'
      elif exclude_depth>=1 and pair in topology['edges']:
        reason='1-2 covalent neighbor'
      elif exclude_depth>=2 and pair in topology['near']:
        reason='1-3 covalent neighbor'
    if reason:
      filtered.append(dict(fragment_A=f,fragment_B=g,reason=reason))
      continue
    distance=min(float(np.linalg.norm(a[ii,None]-a[None,jj],axis=2).min()),
      float(np.linalg.norm(b[ii,None]-b[None,jj],axis=2).min()))
    if distance<=cutoff:
      candidates.append((f,g))
      metadata[(f,g)]=dict(category=category,min_distance_A=distance,
        crosses_covalent_boundary=boundary)
  candidates.sort()
  return candidates,metadata,filtered


def igm(symbols,a,b,frags,directory,step,padding,thresholds,cutoff,max_pairs,max_points,refine,
    mode='fast',refine_top=12,exclude_depth=2,include_coordination=False,bond_file=None,
    ref_ids=None,rtol=.02,atol=1e-5,delta_atol=1e-4,pad_extra=1.,safety=3.):
  tables=density_tables(directory,symbols)
  topology=igm_topology(symbols,a,b,bond_file,ref_ids)
  candidates,metadata,filtered=igm_candidates(symbols,a,b,frags,cutoff,exclude_depth,include_coordination,topology)
  selected=candidates if max_pairs==0 else candidates[:max_pairs]
  rows=[]
  excluded=[]
  def store(record,ints,vols,info):
    record.update(ref_integral=float(ints[0]),probe_integral=float(ints[1]),
      delta_integral=float(ints[1]-ints[0]),**info)
    record['volumes']=[dict(threshold=tau,ref_A3=float(vols[0,i]),probe_A3=float(vols[1,i]),
      delta_A3=float(vols[1,i]-vols[0,i])) for i,tau in enumerate(thresholds)]
  for f,g in selected:
    try:
      ints,vols,info=igm_pair(a,b,symbols,frags[f],frags[g],tables,step,padding,thresholds,max_points)
      row=dict(fragment_A=f,fragment_B=g,integral_unit='bohr^-1',
        resolution='NOT_CHECKED',**metadata[(f,g)])
      store(row,ints,vols,info)
      row['coarse_integrals']=ints.tolist()
      row['coarse_delta']=float(ints[1]-ints[0])
      row['coarse_volumes_A3']=vols.tolist()
      rows.append(row)
    except ValueError as exc:
      excluded.append(dict(fragment_A=f,fragment_B=g,reason=str(exc),stage='coarse'))
  rows.sort(key=lambda r:abs(r['coarse_delta']),reverse=True)
  count=len(rows) if mode=='complete' or refine else min(refine_top,len(rows))
  for row in rows[:count]:
    f,g=row['fragment_A'],row['fragment_B']
    ii,jj=frags[f],frags[g]
    box=igm_box(a,b,ii+jj,step,padding)
    try:
      fine,vfine,info=igm_pair(a,b,symbols,ii,jj,tables,step,padding,thresholds,max_points,box=box,level=1)
      large,vlarge,largeinfo=igm_pair(a,b,symbols,ii,jj,tables,step,padding,thresholds,max_points,
        box=box,level=1,expand=pad_extra)
      coarse=np.array(row['coarse_integrals'])
      delta_fine=float(fine[1]-fine[0])
      delta_large=float(large[1]-large[0])
      grid_err=np.abs(fine-coarse)
      pad_err=np.abs(large-fine)
      grid_delta=abs(delta_fine-row['coarse_delta'])
      pad_delta=abs(delta_large-delta_fine)
      budget=grid_delta+pad_delta
      baseline_ok=bool(np.all(grid_err<=atol+rtol*np.abs(fine)) and
        np.all(pad_err<=atol+rtol*np.abs(large)))
      boundary_ok=max(largeinfo['boundary_max_delta_g'])<min(thresholds)
      sign_ok=row['coarse_delta']*delta_fine>0 and delta_fine*delta_large>0
      resolved=baseline_ok and boundary_ok and sign_ok and abs(delta_large)>max(delta_atol,safety*budget)
      row.update(refined=True,grid_abs_errors=grid_err.tolist(),padding_abs_errors=pad_err.tolist(),
        refinement_relative_errors=(grid_err/np.maximum(abs(fine),1e-12)).tolist(),
        fine_delta=delta_fine,grid_delta_sensitivity=grid_delta,padding_delta_sensitivity=pad_delta,
        sensitivity_budget=budget,baseline_checks_pass=baseline_ok,boundary_check_pass=boundary_ok,
        sign_check_pass=sign_ok,grid_volume_changes_A3=(vfine-np.array(row['coarse_volumes_A3'])).tolist(),
        padding_volume_changes_A3=(vlarge-vfine).tolist(),
        resolution=('RESOLVED_INCREASE' if delta_large>0 else 'RESOLVED_DECREASE') if resolved else 'NUMERICALLY_UNRESOLVED')
      store(row,large,vlarge,largeinfo)
      row['volume_note']='Threshold volumes are reported with numerical sensitivity; no new/lost-contact classification is assigned.'
    except ValueError as exc:
      row.update(resolution='RESOURCE_LIMIT',refinement_reason=str(exc),refined=False)
  rows.sort(key=lambda r:abs(r['delta_integral']),reverse=True)
  checked=[r for r in rows if r['resolution'].startswith('RESOLVED')]
  nonbonded=[r for r in checked if r['category']=='NONBONDED']
  coverage=len(rows)==len(candidates)
  text=f'Coarse scan covers {len(rows)}/{len(candidates)} distance-screened candidates; {len(checked)} tested changes are numerically resolved.'
  if nonbonded:
    lead=max(nonbonded,key=lambda r:abs(r['delta_integral']))
    text+=f" Largest resolved nonbonded descriptor change AMONG CHECKED PAIRS: {lead['fragment_A']} / {lead['fragment_B']} ({lead['delta_integral']:+.4g} bohr^-1)."
  else:
    text+=' No numerically resolved nonbonded increase/decrease is established among checked pairs.'
  if any(r['resolution']=='NOT_CHECKED' for r in rows):
    text+=' Unrefined pairs are not certified; this is not a globally validated ranking.'
  text+=' Pair-IGM changes are reference-density descriptors, not interaction energies.'
  return dict(status='OK' if coverage else 'PARTIAL',pairs=rows,excluded=excluded,filtered=filtered,
    candidate_pairs=len(candidates),evaluated_pairs=len(rows),truncated=len(selected)<len(candidates),
    coarse_coverage_complete=coverage,refined_pairs=sum(r.get('refined',False) for r in rows),
    resolved_increases=sum(r['resolution']=='RESOLVED_INCREASE' for r in rows),
    resolved_decreases=sum(r['resolution']=='RESOLVED_DECREASE' for r in rows),
    thresholds_au=thresholds,step_A=step,padding_A=padding,candidate_cutoff_A=cutoff,mode=mode,
    source={x:tables[x][3] for x in tables},conclusion=text,topology=dict(source=topology['source'],
      changed_edges=topology['changed_edges'],metal_indices=sorted(topology['metals'])),
    numerical_tolerances=dict(rtol=rtol,atol=atol,delta_atol=delta_atol,padding_extra_A=pad_extra,safety=safety),
    note='1-2/1-3 exclusions apply to singleton pairs. Named fragment boundaries are classified, not subtracted. Numerical checks are sensitivity tests, not rigorous error bounds.')


def igm_single(symbols,a,frags,directory,step,padding,thresholds,cutoff,max_pairs,max_points,refine,
    mode='fast',refine_top=12,exclude_depth=2,include_coordination=False,bond_file=None,
    ref_ids=None,rtol=.02,atol=1e-5,pad_extra=1.):
  tables=density_tables(directory,symbols)
  topology=igm_topology(symbols,a,a,bond_file,ref_ids)
  candidates,metadata,filtered=igm_candidates(symbols,a,a,frags,cutoff,exclude_depth,include_coordination,topology)
  selected=candidates if max_pairs==0 else candidates[:max_pairs]
  rows=[]
  excluded=[]
  def store(row,ints,vols,info):
    row.update(integral=float(ints[0]),**info)
    row['volumes']=[dict(threshold=tau,volume_A3=float(vols[0,i])) for i,tau in enumerate(thresholds)]
  for f,g in selected:
    try:
      vals,vols,info=igm_pair(a,None,symbols,frags[f],frags[g],tables,step,padding,thresholds,max_points)
      row=dict(fragment_A=f,fragment_B=g,integral_unit='bohr^-1',resolution='NOT_CHECKED',**metadata[(f,g)])
      store(row,vals,vols,info)
      row['coarse_integral']=float(vals[0])
      row['coarse_volumes_A3']=vols[0].tolist()
      rows.append(row)
    except ValueError as exc:
      excluded.append(dict(fragment_A=f,fragment_B=g,stage='coarse',reason=str(exc)))
  rows.sort(key=lambda row:row['coarse_integral'],reverse=True)
  count=len(rows) if mode=='complete' or refine else min(refine_top,len(rows))
  for row in rows[:count]:
    ii,jj=frags[row['fragment_A']],frags[row['fragment_B']]
    box=igm_box(a,None,ii+jj,step,padding)
    try:
      fine,vfine,info=igm_pair(a,None,symbols,ii,jj,tables,step,padding,thresholds,max_points,box=box,level=1)
      large,vlarge,largeinfo=igm_pair(a,None,symbols,ii,jj,tables,step,padding,thresholds,max_points,
        box=box,level=1,expand=pad_extra)
      grid=abs(float(fine[0])-row['coarse_integral'])
      pad=abs(float(large[0]-fine[0]))
      grid_ok=grid<=atol+rtol*abs(float(fine[0]))
      pad_ok=pad<=atol+rtol*abs(float(large[0]))
      bound=max(largeinfo['boundary_max_delta_g'])<min(thresholds)
      state='GRID_SENSITIVE' if not grid_ok else 'BOX_SENSITIVE' if not pad_ok or not bound else 'CHECKED_CONTACT'
      row.update(refined=True,fine_integral=float(fine[0]),grid_abs_sensitivity=grid,padding_abs_sensitivity=pad,
        grid_check_pass=bool(grid_ok),padding_check_pass=bool(pad_ok),boundary_check_pass=bool(bound),
        resolution=state,grid_volume_changes_A3=(vfine[0]-np.array(row['coarse_volumes_A3'])).tolist(),
        padding_volume_changes_A3=(vlarge[0]-vfine[0]).tolist())
      store(row,large,vlarge,largeinfo)
    except ValueError as exc:
      row.update(refined=False,resolution='RESOURCE_LIMIT',refinement_reason=str(exc))
  rows.sort(key=lambda row:row['integral'],reverse=True)
  checked=[r for r in rows if r['resolution']=='CHECKED_CONTACT' and r['category']=='NONBONDED']
  coverage=len(rows)==len(candidates)
  text=f'Absolute reference-density contacts: {len(rows)}/{len(candidates)} screened candidates evaluated.'
  if checked:
    lead=max(checked,key=lambda row:row['integral'])
    text+=f" Largest checked NONBONDED descriptor among checked pairs: {lead['fragment_A']} / {lead['fragment_B']} ({lead['integral']:.4g} bohr^-1)."
  else:
    text+=' No numerically checked nonbonded contact is established.'
  if any(r['resolution']=='NOT_CHECKED' for r in rows):
    text+=' Remaining coarse ranks are unverified.'
  text+=' Descriptor size is not attraction, steric repulsion or interaction energy.'
  return dict(status='OK' if coverage else 'PARTIAL',single=True,pairs=rows,filtered=filtered,excluded=excluded,
    candidate_pairs=len(candidates),evaluated_pairs=len(rows),coarse_coverage_complete=coverage,
    truncated=len(selected)<len(candidates),refined_pairs=sum(r.get('refined',False) for r in rows),
    checked_pairs=sum(r['resolution']=='CHECKED_CONTACT' for r in rows),mode=mode,
    thresholds_au=thresholds,step_A=step,padding_A=padding,candidate_cutoff_A=cutoff,
    source={x:tables[x][3] for x in tables},topology=dict(source=topology['source'],
      metal_indices=sorted(topology['metals']),changed_edges=topology['changed_edges']),
    numerical_tolerances=dict(rtol=rtol,atol=atol,padding_extra_A=pad_extra),conclusion=text,
    note='Absolute-integral checks, not self-comparison cancellation. Pair descriptors are nonadditive; volumes are threshold-dependent.')


def dispersion_single(symbols,a,charge,method,frags,ref_ids,bond_file=None):
  from dftd4.interface import DispersionModel,DampingParam
  from importlib.metadata import version
  model=DispersionModel(np.array([ELEMENTS.index(x) for x in symbols],dtype=np.int32),a/BOHR,charge=float(charge))
  param=DampingParam(method=method)
  energy=float(np.asarray(model.get_dispersion(param,grad=False)['energy']))
  pair=model.get_pairwise_dispersion(param)
  additive=np.asarray(pair['additive pairwise energy'])
  nonadd=np.asarray(pair['non-additive pairwise energy'])
  if additive.shape!=(len(a),len(a)) or nonadd.shape!=additive.shape:
    raise ValueError('Unexpected D4 pair matrix shape')
  if abs(float(additive.sum()+nonadd.sum())-energy)>1e-8:
    raise ValueError('D4 pairwise sum does not match total energy')
  topology=igm_topology(symbols,a,a,bond_file,ref_ids)
  owner={i:name for name,group in frags.items() for i in group}
  rows=[]
  for i in range(len(a)):
    for j in range(i+1,len(a)):
      pairid=(i,j)
      if i in topology['metals'] or j in topology['metals']:
        category='METAL_INVOLVING'
      elif pairid in topology['edges']:
        category='COVALENT_12'
      elif pairid in topology['near']:
        category='COVALENT_13'
      else:
        category='NONBONDED'
      rows.append(dict(atoms=f'{ref_ids[i]+1}-{ref_ids[j]+1}',element_pair=f'{symbols[i]}-{symbols[j]}',
        distance_A=float(np.linalg.norm(a[i]-a[j])),energy_kcal=float(additive[i,j]+additive[j,i])*KCAL,
        category=category,fragment_A=owner[i],fragment_B=owner[j],same_fragment=owner[i]==owner[j]))
  rows.sort(key=lambda row:row['energy_kcal'])
  groups=[]
  for f,g in itertools.combinations_with_replacement(frags,2):
    ii,jj=frags[f],frags[g]
    val=additive[np.ix_(ii,jj)].sum()
    if f!=g:
      val+=additive[np.ix_(jj,ii)].sum()
    groups.append(dict(fragment_A=f,fragment_B=g,energy_kcal=float(val)*KCAL,same_fragment=f==g))
  groups.sort(key=lambda row:row['energy_kcal'])
  summaries=[dict(category=cat,energy_kcal=sum(r['energy_kcal'] for r in rows if r['category']==cat),
    pair_count=sum(r['category']==cat for r in rows)) for cat in sorted(set(r['category'] for r in rows))]
  nb=[r for r in rows if r['category']=='NONBONDED']
  text='Absolute D4 correction is not total stability, unfolding energy or a fragment binding energy.'
  if nb:
    text+=f" Largest-magnitude attractive NONBONDED two-body contribution in this ranking: {nb[0]['atoms']} ({nb[0]['energy_kcal']:.4g} kcal/mol)." if nb[0]['energy_kcal']<0 else ' No negative nonbonded two-body contribution found.'
  return dict(status='OK',single=True,method=method,charge=charge,version=version('dftd4'),
    total_kcal=energy*KCAL,additive_kcal=float(additive.sum())*KCAL,three_body_kcal=float(nonadd.sum())*KCAL,
    pairs=rows,fragments=groups,categories=summaries,conclusion=text,
    note='Topology categories are model labels; all pair terms remain included in total energy. No bond-cutting reference is defined.')


def polarization_single(symbols,a,frags,params,ref_ids,axis):
  alpha,damping,source=polar_params(params,symbols,ref_ids)
  tensor,response,condition=polar_tensor(a,alpha,damping)
  vals,vec=eigh(tensor)
  anis=float(np.sqrt(((vals[0]-vals[1])**2+(vals[1]-vals[2])**2+(vals[2]-vals[0])**2)/2))*BOHR**3
  rows=[]
  for name,group in frags.items():
    contribution=response[group].sum(0)*BOHR**3
    rows.append(dict(fragment=name,trace_mean_A3=float(np.trace(contribution)/3),tensor_A3=contribution))
  rows.sort(key=lambda row:abs(row['trace_mean_A3']),reverse=True)
  out=dict(status='OK',single=True,source=source,damping=damping,condition=condition,
    tensor_A3=tensor*BOHR**3,eigenvalues_A3=vals*BOHR**3,eigenvectors=vec,
    mean_A3=float(np.trace(tensor)/3*BOHR**3),anisotropy_A3=anis,fragments=rows,
    conclusion='Classical static-response tensor under the stated parameters; not quantum spin-resolved response or Raman activity. Fragment terms are coupled response contributions, not isolated-fragment polarizabilities.')
  if axis is not None:
    i,j=axis
    direction=a[j]-a[i]
    if np.linalg.norm(direction)<1e-8:
      raise ValueError('Zero response axis')
    direction/=np.linalg.norm(direction)
    out['axis_A3']=float(direction@tensor@direction*BOHR**3)
  return out


def export_single(prefix,result,overwrite):
  root=Path(prefix)
  datasets={}
  modules=result['modules']
  h=modules.get('basis_health',{})
  if h.get('status')=='OK':
    datasets['health_shells']=h['value']['shell_weights']
    datasets['kinetic']=h['value']['kinetic']
  fp=modules.get('fingerprint',{})
  if fp.get('status')=='OK':
    datasets['om_spectrum']=[dict(mode=i+1,eigenvalue=x) for i,x in enumerate(fp['spectrum'])]
  for key,tag in [('dispersion','d4'),('polarization','polar')]:
    item=modules.get(key,{})
    if item.get('status')=='OK':
      datasets[tag+'_fragments']=item['fragments']
      if key=='dispersion':
        datasets['d4_pairs']=item['pairs']
        datasets['d4_categories']=item['categories']
      else:
        datasets['polar_tensor']=[dict(i=i,j=j,value_A3=item['tensor_A3'][i,j]) for i in range(3) for j in range(3)]
  ig=modules.get('igm',{})
  for field in ('pairs','filtered','excluded'):
    datasets['igm_'+field]=ig.get(field,[])
  datasets['igm_volumes']=[dict(fragment_A=r['fragment_A'],fragment_B=r['fragment_B'],resolution=r['resolution'],**v)
    for r in ig.get('pairs',[]) for v in r['volumes']]
  paths=[Path(str(root)+'_quantmol.json')]+[Path(str(root)+'_'+name+'.csv') for name,rows in datasets.items() if rows]
  if not overwrite and any(path.exists() for path in paths):
    raise FileExistsError('Output already exists; use --overwrite')
  root.parent.mkdir(parents=True,exist_ok=True)
  paths[0].write_text(json.dumps(clean(result),indent=2,allow_nan=False))
  for name,rows in datasets.items():
    write_csv(str(root)+'_'+name+'.csv',rows)


def report_single(name,data,top,ref_ids):
  print('\n'+name.upper().center(74,'='))
  print('Status:',data['status'])
  if data['status'] not in ('OK','PARTIAL'):
    print(data.get('reason',data.get('conclusion','')))
    return
  if name=='basis_health':
    h=data['value']
    cond='undefined' if h['condition'] is None else f"{h['condition']:.3e}"
    print(f"lambda_min={h['minimum']:.3e}; condition={cond}; low modes={h['low_counts']}")
    tau=max(float(x) for x in h['low_counts'])
    rows=sorted([r for r in h['shell_weights'] if r['threshold']==tau],key=lambda r:r['weight'],reverse=True)
    for r in rows[:top]:
      print(f"  Atom {ref_ids[r['atom']]+1}; shell {r['shell']}; l={r['l']}; low-space weight={r['weight']:.4g}")
    print('Kinetic maxima (Ha): '+', '.join(f"{r['threshold']:.0e}:{r['max_Ha']:.5g}" for r in h['kinetic']))
  elif name=='fingerprint':
    print(f"OM signature saved: {len(data['spectrum'])} eigenvalues; analysis basis={data['basis']}")
  elif name=='dispersion':
    print(f"D4 {data['method']}: total={data['total_kcal']:.6f}; two-body={data['additive_kcal']:.6f}; three-body={data['three_body_kcal']:.6f} kcal/mol")
    for row in data['categories']:
      print(f"  {row['category']}: {row['energy_kcal']:.5f} kcal/mol; pairs={row['pair_count']}")
    for row in [r for r in data['pairs'] if r['category']=='NONBONDED'][:top]:
      print(f"  {row['atoms']} ({row['element_pair']}): {row['distance_A']:.3f} A; {row['energy_kcal']:.5f} kcal/mol")
    if data['named_fragments']:
      print('  Leading cross-fragment two-body sums (not binding energies):')
      for row in [r for r in data['fragments'] if not r['same_fragment']][:top]:
        print(f"    {row['fragment_A']} / {row['fragment_B']}: {row['energy_kcal']:.5f} kcal/mol")
  elif name=='polarization':
    print(f"Mean={data['mean_A3']:.5f}; anisotropy={data['anisotropy_A3']:.5f} A^3; response condition={data['condition']:.3g}")
    if 'axis_A3' in data:
      print(f"Specified axis response={data['axis_A3']:.5f} A^3")
  elif name=='igm':
    print(f"Candidates={data['candidate_pairs']}; evaluated={data['evaluated_pairs']}; coarse coverage={data['coarse_coverage_complete']}")
    print(f"Refined={data['refined_pairs']}; checked contacts={data['checked_pairs']}; mode={data['mode']}")
    # Keep absolute-integral order; do not silently re-rank by verification status.
    for row in data['pairs'][:top]:
      print(f"  {row['fragment_A']} / {row['fragment_B']} [{row['category']}]: J={row['integral']:.5g} bohr^-1; {row['resolution']}")
      if 'grid_abs_sensitivity' in row:
        print(f"    Absolute J grid/padding sensitivity={row['grid_abs_sensitivity']:.2g}/{row['padding_abs_sensitivity']:.2g}")
      elif 'refinement_reason' in row:
        print('    '+row['refinement_reason'])
    print(f"Topology-filtered={len(data['filtered'])}; failed coarse grids={len(data['excluded'])}")
  print('Conclusion:',data.get('conclusion',''))
  print(f"Elapsed: {data.get('elapsed_s',0):.3f} s")


def main_single(args,thresholds,it):
  if args.probe_atoms or args.align_atoms:
    raise ValueError('--probe-atoms and --align-atoms require two-conformer mode')
  symbols0,coords=read_xyz(args.ref_xyz)
  text='c' if args.all else args.ref_atoms or input('Enter atoms (1-based; C for all): ')
  ids=indices(text,len(symbols0))
  symbols=[symbols0[i] for i in ids]
  a=coords[ids]
  lookup={original:i for i,original in enumerate(ids)}
  fs=fragments(args.fragments,len(a),ids)
  axis=None
  if args.response_axis:
    raw=args.response_axis.split(':')
    if len(raw)!=2:
      raise ValueError('--response-axis requires i:j')
    picked=indices(','.join(raw),len(symbols0))
    if len(picked)!=2 or any(i not in lookup for i in picked):
      raise ValueError('Response axis must use two selected original atom indices')
    axis=[lookup[i] for i in picked]
  print('Mode: SINGLE CONFORMER (fixed geometry)')
  print('Selected original atoms:',', '.join(str(i+1) for i in ids))
  print('Selected atoms form the complete modeled subsystem; no caps or omitted environment.')
  print(f'D4 damping={args.d4_method}; input charge={args.charge}; no comparison/metric is computed.')
  result=dict(version=VERSION,analysis_mode='single',settings=vars(args),
    mapping=dict(atoms=[i+1 for i in ids]),fragments=fs,coordinate_frame='input XYZ, no rigid alignment',modules={})
  modules=result['modules']
  def run(name,func):
    start=time.perf_counter()
    try:
      item=func()
    except (ImportError,ValueError,OSError,KeyError,RuntimeError,np.linalg.LinAlgError) as exc:
      item=dict(status='UNAVAILABLE',single=True,reason=f'{type(exc).__name__}: {exc}')
    item['elapsed_s']=time.perf_counter()-start
    modules[name]=item
    report_single(name,item,args.top,ids)
    return item
  if not args.skip_integrals:
    try:
      started=time.perf_counter()
      sm,tm,meta,backend=build_integrals(symbols,a,args.basis,args.ecp,args.max_ao)
      result['integrals']=dict(**backend,elapsed_s=time.perf_counter()-started)
      def hs():
        val,_=health(sm,tm,meta,thresholds)
        n=val['low_counts'][str(max(thresholds))]
        text=f'{n} near-redundant modes below {max(thresholds):g}; inspect shell participation before changing the basis.' if n else f'No overlap eigenvalues below {max(thresholds):g}; no near-redundancy warning at this threshold.'
        return dict(status='OK',single=True,value=val,conclusion=text+' This is not a general SCF convergence guarantee.')
      run('basis_health',hs)
      def fps():
        matrix=sm
        if args.fingerprint_basis and args.fingerprint_basis!=args.basis:
          matrix,_,_,_=build_integrals(symbols,a,args.fingerprint_basis,args.ecp,args.max_ao)
        return dict(status='OK',single=True,spectrum=eigh(matrix,eigvals_only=True),
          basis=args.fingerprint_basis or args.basis,
          conclusion='Structural overlap signature for future comparisons; no electronic-state or anomaly label is assigned without a reference.')
      run('fingerprint',fps)
    except (ImportError,ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
      for name in ('basis_health','fingerprint'):
        modules[name]=dict(status='UNAVAILABLE',single=True,reason=f'{type(exc).__name__}: {exc}')
        report_single(name,modules[name],args.top,ids)
  modules['metric']=dict(status='NOT_APPLICABLE',single=True,conclusion='Metric distortion requires two conformers; intentionally not calculated.')
  if not args.skip_d4:
    def ds():
      item=dispersion_single(symbols,a,args.charge,args.d4_method,fs,ids,args.igm_bonds)
      item['named_fragments']=bool(args.fragments)
      return item
    run('dispersion',ds)
  if not args.skip_polar:
    run('polarization',lambda:polarization_single(symbols,a,fs,args.polar_params,ids,axis))
  if not args.skip_igm:
    run('igm',lambda:igm_single(symbols,a,fs,args.density_library,args.igm_step,args.igm_padding,it,
      args.igm_cutoff,args.igm_max_pairs,args.igm_max_points,args.igm_refine,mode=args.igm_mode,
      refine_top=args.igm_refine_top,exclude_depth=args.igm_exclude_depth,
      include_coordination=args.igm_include_coordination,bond_file=args.igm_bonds,ref_ids=ids,
      rtol=args.igm_rtol,atol=args.igm_atol,pad_extra=args.igm_padding_extra))
  bad=[name for name,item in modules.items() if item['status'] not in ('OK','NOT_APPLICABLE')]
  result['run_status']='PARTIAL' if bad else 'COMPLETE'
  result['unavailable_or_withheld']=bad
  print('\n'+'SUMMARY'.center(74,'='))
  print('Run status:',result['run_status'],'; absolute diagnostics only, no stability/reactivity score.')
  if bad:
    print('Unavailable/partial:',', '.join(bad))
  if args.csv_prefix:
    export_single(args.csv_prefix,result,args.overwrite)
    print('Single-conformer JSON and CSV exported with prefix:',args.csv_prefix)
  return 3 if bad and args.strict else 0


def clean(obj):
  if isinstance(obj,np.ndarray):
    return clean(obj.tolist())
  if isinstance(obj,np.generic):
    return clean(obj.item())
  if isinstance(obj,float) and not math.isfinite(obj):
    return None
  if isinstance(obj,dict):
    return {str(k):clean(v) for k,v in obj.items()}
  if isinstance(obj,(list,tuple)):
    return [clean(x) for x in obj]
  return obj


def write_csv(path,rows):
  if not rows:
    return
  fields=list(dict.fromkeys(k for r in rows for k in r))
  with open(path,'w',newline='') as f:
    writer=csv.DictWriter(f,fieldnames=fields)
    writer.writeheader()
    for r in rows:
      writer.writerow({k:json.dumps(clean(v)) if isinstance(v,(list,dict,np.ndarray)) else clean(v) for k,v in r.items()})


def export(prefix,result,overwrite):
  if result.get('analysis_mode')=='single':
    return export_single(prefix,result,overwrite)
  root=Path(prefix)
  paths=[Path(str(root)+'_quantmol.json')]
  datasets={}
  modules=result['modules']
  h=modules.get('basis_health',{})
  if h.get('status')=='OK':
    datasets['health_shells']=[dict(conformer=c,**r) for c in ('ref','probe') for r in h[c]['shell_weights']]
    datasets['kinetic']=[dict(conformer=c,**r) for c in ('ref','probe') for r in h[c]['kinetic']]
  fp=modules.get('fingerprint',{})
  if fp.get('status')=='OK':
    datasets['om_local']=fp['local']
    datasets['om_spectrum']=[dict(mode=i+1,ref=x,probe=y,delta=y-x) for i,(x,y) in enumerate(zip(fp['spectrum_ref'],fp['spectrum_probe']))]
  met=modules.get('metric',{})
  if met.get('status')=='OK':
    datasets['metric_spectrum']=[dict(mode=i+1,gamma=x,log_gamma=y) for i,(x,y) in enumerate(zip(met['gamma'],met['log_gamma']))]
  dd=modules.get('dispersion',{})
  if dd.get('status')=='OK':
    datasets['d4_pairs']=dd['pairs']
    datasets['d4_fragments']=dd['fragments']
  pp=modules.get('polarization',{})
  if pp.get('status')=='OK':
    datasets['polar_fragments']=pp['fragments']
    datasets['polar_tensor']=[dict(i=i,j=j,ref_A3=pp['ref']['tensor_A3'][i,j],probe_A3=pp['probe']['tensor_A3'][i,j]) for i in range(3) for j in range(3)]
  gg=modules.get('igm',{})
  datasets['igm_pairs']=gg.get('pairs',[])
  datasets['igm_filtered']=gg.get('filtered',[])
  datasets['igm_excluded']=gg.get('excluded',[])
  datasets['igm_volumes']=[dict(fragment_A=r['fragment_A'],fragment_B=r['fragment_B'],resolution=r['resolution'],**v) for r in gg.get('pairs',[]) for v in r['volumes']]
  paths += [Path(str(root)+'_'+name+'.csv') for name,rows in datasets.items() if rows]
  if not overwrite and any(p.exists() for p in paths):
    raise FileExistsError('Output already exists; use --overwrite')
  root.parent.mkdir(parents=True,exist_ok=True)
  paths[0].write_text(json.dumps(clean(result),indent=2,allow_nan=False))
  for name,rows in datasets.items():
    write_csv(str(root)+'_'+name+'.csv',rows)


def report(name,data,top,ref_ids):
  print('\n'+name.upper().center(74,'='))
  status=data.get('status','UNKNOWN')
  print('Status:',status)
  if status not in ('OK','PARTIAL'):
    print(data.get('conclusion',data.get('reason','No result.')))
    return
  if name=='basis_health':
    for c in ('ref','probe'):
      h=data[c]
      cond=f"{h['condition']:.3e}" if h['condition'] is not None else 'undefined'
      print(f"{c.upper():5s}: lambda_min={h['minimum']:.3e}; condition={cond}; low modes={h['low_counts']}")
      if h['shell_weights']:
        tau=max(float(x) for x in h['low_counts'])
        rows=sorted([r for r in h['shell_weights'] if r['threshold']==tau],key=lambda r:r['weight'],reverse=True)
        for r in rows[:top]:
          print(f"  REF-index {ref_ids[r['atom']]+1}; shell {r['shell']}; l={r['l']}; exponent_min={r['min_exponent']:.3g}; low-space weight={r['weight']:.3g}")
      if h['kinetic']:
        print('  Retained-space kinetic maximum (Ha): '+', '.join(f"{r['threshold']:.0e}:{r['max_Ha']:.4g}" for r in h['kinetic']))
  elif name=='fingerprint':
    print(f"Global OM distance={data['distance']:.6g}; RMS={data['rms']:.6g}")
    for r in data['local'][:top]:
      print(f"  REF center {ref_ids[r['center']]+1}: local RMS={r['rms']:.5g}; common neighborhood atoms={len(r['neighbors'])}")
  elif name=='metric':
    print(f"Log distance={data['log_distance']:.6g}; log RMS={data['log_rms']:.6g}; max |log gamma|={data['max_abs_log']:.6g}")
    print(f"Squared-norm ratio range=[{data['gamma'][0]:.6g}, {data['gamma'][-1]:.6g}]; dimension={data['dimension']}; projected={data['projected']}")
  elif name=='dispersion':
    print(f"D4 {data['method']}: REF={data['ref_kcal']:.6f}; PROBE={data['probe_kcal']:.6f}; delta={data['delta_kcal']:+.6f} kcal/mol")
    print(f"Two-body delta={data['additive_delta_kcal']:+.6f}; three-body delta={data['three_body_delta_kcal']:+.6f}")
    for r in data['fragments'][:top]:
      print(f"  {r['fragment_A']} / {r['fragment_B']}: {r['delta_kcal']:+.5f} kcal/mol (two-body)")
  elif name=='polarization':
    for c in ('ref','probe'):
      print(f"{c.upper()}: mean={data[c]['mean_A3']:.5f}; anisotropy={data[c]['anisotropy_A3']:.5f} A^3")
    if 'axis_ref_A3' in data:
      print(f"Common REF axis: {data['axis_ref_A3']:.5f} -> {data['axis_probe_A3']:.5f} A^3")
    for r in data['fragments'][:top]:
      print(f"  {r['fragment']}: response tensor change={r['response_difference_F_A3']:.5f} A^3")
    print('Model:',data['source'])
  elif name=='igm':
    print(f"Candidates={data['candidate_pairs']}; coarse-evaluated={data['evaluated_pairs']}; complete-coverage={data['coarse_coverage_complete']}")
    print(f"Refined={data['refined_pairs']}; resolved + / - = {data['resolved_increases']} / {data['resolved_decreases']}; mode={data['mode']}")
    ranked=sorted(data['pairs'],key=lambda r:(not r['resolution'].startswith('RESOLVED'),-abs(r['delta_integral'])))
    for r in ranked[:top]:
      print(f"  {r['fragment_A']} / {r['fragment_B']} [{r['category']}]: {r['delta_integral']:+.5g} bohr^-1; {r['resolution']}")
      if 'sensitivity_budget' in r:
        print(f"    Grid/padding delta sensitivity={r['grid_delta_sensitivity']:.2g}/{r['padding_delta_sensitivity']:.2g}; baseline/boundary pass={r['baseline_checks_pass']}/{r['boundary_check_pass']}")
      elif 'refinement_reason' in r:
        print('    '+r['refinement_reason'])
    print(f"Excluded topology pairs={len(data['filtered'])}; failed coarse grids={len(data['excluded'])}")
    if data['topology']['changed_edges']:
      print('WARNING: inferred covalent graph changes; inspect topology or supply --igm-bonds.')
  print('Conclusion:',data.get('conclusion',''))
  print(f"Elapsed: {data.get('elapsed_s',0):.3f} s")


def parser():
  p=argparse.ArgumentParser(description=__doc__)
  p.add_argument('ref_xyz',nargs='?')
  p.add_argument('probe_xyz',nargs='?')
  p.add_argument('--ref-atoms')
  p.add_argument('--probe-atoms')
  p.add_argument('--all',action='store_true')
  p.add_argument('--align-atoms',help='Selected original REF atom numbers')
  p.add_argument('--basis',default='def2-svp')
  p.add_argument('--fingerprint-basis',help='Optional separate PySCF basis name')
  p.add_argument('--ecp',help='Optional PySCF ECP name; always recorded')
  p.add_argument('--charge',type=int,default=0,help='Charge of the complete analyzed molecule/fragment')
  p.add_argument('--fragments',help='Disjoint named fragments, original REF indices, JSON')
  p.add_argument('--top',type=int,default=5)
  p.add_argument('--thresholds',default='1e-6,1e-8,1e-10')
  p.add_argument('--metric-threshold',type=float,default=1e-8)
  p.add_argument('--metric-project',action='store_true')
  p.add_argument('--local-centers',help='Original REF indices; default all selected atoms')
  p.add_argument('--local-radius',type=float,default=4.)
  p.add_argument('--max-ao',type=int,default=1500)
  p.add_argument('--d4-method',default='pbe0',help='Damping parameter set, not a molecular DFT calculation')
  p.add_argument('--polar-params',help='Custom published/user-calibrated linear Thole JSON')
  p.add_argument('--response-axis',help='Original REF atom indices, e.g. 1:25')
  p.add_argument(
    '--density-library',
    default='/home/lky/gv/polygview/atomic_density',
    help=(
      'Directory of radial atomic-density JSON files '
      '(default: /home/lky/gv/polygview/atomic_density)'
    )
  )
  p.add_argument('--igm-step',type=float,default=.25)
  p.add_argument('--igm-padding',type=float,default=2.)
  p.add_argument('--igm-thresholds',default='1e-4,1e-3,1e-2')
  p.add_argument('--igm-cutoff',type=float,default=4.)
  p.add_argument('--igm-max-pairs',type=int,default=0,help='Explicit coarse-scan cap; 0 scans all candidates')
  p.add_argument('--igm-max-points',type=int,default=1000000)
  p.add_argument('--igm-refine',action='store_true')
  p.add_argument('--igm-mode',choices=('fast','complete'),default='fast')
  p.add_argument('--igm-refine-top',type=int,default=12)
  p.add_argument('--igm-exclude-depth',type=int,choices=(0,1,2),default=2)
  p.add_argument('--igm-include-coordination',action='store_true')
  p.add_argument('--igm-bonds',help='JSON nonmetal covalent bond list; original REF indices')
  p.add_argument('--igm-rtol',type=float,default=.02)
  p.add_argument('--igm-atol',type=float,default=1e-5)
  p.add_argument('--igm-delta-atol',type=float,default=1e-4)
  p.add_argument('--igm-padding-extra',type=float,default=1.)
  p.add_argument('--igm-safety',type=float,default=3.)
  p.add_argument('--prepare-density',help='Separate atomic HF cache preparation directory')
  p.add_argument('--density-elements',default='H,C,N,O')
  p.add_argument('--density-basis',default='def2-tzvp')
  p.add_argument('--csv-prefix',help='Write JSON and CSV; no export if omitted')
  p.add_argument('--overwrite',action='store_true')
  p.add_argument('--strict',action='store_true',help='Exit 3 if any requested module cannot run')
  for name in ('integrals','d4','polar','igm'):
    p.add_argument('--skip-'+name,action='store_true')
  p.add_argument('--version',action='version',version=VERSION)
  return p


def main(argv=None):
  p=parser()
  args=p.parse_args(argv)
  if args.prepare_density:
    elements=[x.strip().capitalize() for x in args.density_elements.split(',')]
    if any(x not in ELEMENTS[1:] for x in elements):
      p.error('Invalid --density-elements')
    root=prepare_density(elements,args.density_basis,args.prepare_density,args.overwrite)
    print('Prepared explicitly requested atomic HF references:',root)
    return 0
  if not args.ref_xyz:
    p.error('Provide one XYZ for single mode, or REF and PROBE XYZ for comparison')
  thresholds=[float(x) for x in args.thresholds.split(',')]
  it=[float(x) for x in args.igm_thresholds.split(',')]
  positives=thresholds+it+[args.metric_threshold,args.local_radius,args.igm_step,args.igm_padding,args.igm_cutoff,
    args.igm_rtol,args.igm_atol,args.igm_delta_atol,args.igm_padding_extra,args.igm_safety]
  if any(not np.isfinite(x) or x<=0 for x in positives):
    p.error('Thresholds, lengths and steps must be positive finite numbers')
  if args.igm_max_pairs<0 or min(args.top,args.max_ao,args.igm_max_points,args.igm_refine_top)<1:
    p.error('Counts must be positive')
  if args.all and (args.ref_atoms or args.probe_atoms):
    p.error('--all conflicts with selections')
  if args.probe_xyz and bool(args.ref_atoms)!=bool(args.probe_atoms):
    p.error('Provide both paired selections')
  if not args.probe_xyz:
    return main_single(args,thresholds,it)
  sa,a0=read_xyz(args.ref_xyz)
  sb,b0=read_xyz(args.probe_xyz)
  if args.all:
    ra,pb=list(range(len(sa))),list(range(len(sb)))
  else:
    rt=args.ref_atoms if args.ref_atoms else input('Enter REF atoms (1-based; C for all): ')
    pt=args.probe_atoms if args.probe_atoms else input('Enter PROBE atoms (paired order; C for all): ')
    ra,pb=indices(rt,len(sa)),indices(pt,len(sb))
  if len(ra)!=len(pb) or any(sa[i]!=sb[j] for i,j in zip(ra,pb)):
    raise ValueError('Paired selections must have equal size and matching elements')
  symbols=[sa[i] for i in ra]
  a,b=a0[ra],b0[pb]
  lookup={j:i for i,j in enumerate(ra)}
  def selected(text):
    ids=indices(text,len(sa))
    if any(i not in lookup for i in ids):
      raise ValueError('Requested atom is not in the analyzed selection')
    return [lookup[i] for i in ids]
  als=selected(args.align_atoms) if args.align_atoms else None
  b,alignment=align(a,b,als)
  fs=fragments(args.fragments,len(a),ra)
  axis=None
  if args.response_axis:
    parts=args.response_axis.split(':')
    if len(parts)!=2:
      p.error('--response-axis requires i:j')
    axis=selected(','.join(parts))
    if len(axis)!=2:
      p.error('Axis requires two different atoms')
  print('Atom mapping (REF:PROBE):',', '.join(f'{i+1}:{j+1}' for i,j in zip(ra,pb)))
  print('Selected atoms form the COMPLETE modeled subsystem; cut covalent fragments are not capped.')
  print(f"Fixed geometry; PROBE aligned to REF; alignment RMSD={alignment['rmsd_A']:.6f} A")
  print('Defaults: D4 damping='+args.d4_method+'; neutral charge='+str(args.charge))
  if alignment['warning']:
    print('WARNING:',alignment['warning'])
  result=dict(version=VERSION,settings=vars(args),mapping=dict(ref=[x+1 for x in ra],probe=[x+1 for x in pb]),
    alignment=alignment,fragments=fs,modules={})
  modules=result['modules']
  def run(name,func):
    start=time.perf_counter()
    try:
      data=func()
    except (ImportError,ValueError,OSError,KeyError,RuntimeError,np.linalg.LinAlgError) as exc:
      data=dict(status='UNAVAILABLE',reason=f'{type(exc).__name__}: {exc}')
    data['elapsed_s']=time.perf_counter()-start
    modules[name]=data
    report(name,data,args.top,ra)
    return data
  if not args.skip_integrals:
    start=time.perf_counter()
    try:
      s1,t1,meta,back=build_integrals(symbols,a,args.basis,args.ecp,args.max_ao)
      s2,t2,meta2,_=build_integrals(symbols,b,args.basis,args.ecp,args.max_ao)
      if meta!=meta2:
        raise ValueError('AO shell metadata mismatch')
      result['integrals']=dict(**back,elapsed_s=time.perf_counter()-start)
      def health_pair():
        h1,_=health(s1,t1,meta,thresholds)
        h2,_=health(s2,t2,meta,thresholds)
        tau=max(thresholds)
        c1,c2=h1['low_counts'][str(tau)],h2['low_counts'][str(tau)]
        if c1!=c2:
          text=f'Low-overlap mode count changes {c1}->{c2} at {tau:g}; effective-rank handling deserves review.'
        elif c1:
          text=f'Both structures have {c1} near-redundant modes at {tau:g}; shell weights locate participation, not unique causation.'
        else:
          text=f'Neither structure has overlap eigenvalues below {tau:g}; no near-redundancy warning at the specified diagnostic threshold.'
        return dict(status='OK',ref=h1,probe=h2,conclusion=text+' Kinetic spectra are auxiliary diagnostics, not SCF convergence predictions.')
      run('basis_health',health_pair)
      centers=selected(args.local_centers) if args.local_centers else list(range(len(a)))
      def fp():
        if args.fingerprint_basis and args.fingerprint_basis!=args.basis:
          f1,_,fm,_=build_integrals(symbols,a,args.fingerprint_basis,args.ecp,args.max_ao)
          f2,_,_,_=build_integrals(symbols,b,args.fingerprint_basis,args.ecp,args.max_ao)
        else:
          f1,f2,fm=s1,s2,meta
        return fingerprint(f1,f2,fm,a,b,centers,args.local_radius)
      run('fingerprint',fp)
      run('metric',lambda:metric(s1,s2,args.metric_threshold,args.metric_project))
    except (ImportError,ValueError,RuntimeError,np.linalg.LinAlgError) as exc:
      for name in ('basis_health','fingerprint','metric'):
        modules[name]=dict(status='UNAVAILABLE',reason=f'{type(exc).__name__}: {exc}',elapsed_s=time.perf_counter()-start)
        report(name,modules[name],args.top,ra)
  if not args.skip_d4:
    run('dispersion',lambda:dispersion(symbols,a,b,args.charge,args.d4_method,fs,ra,pb))
  if not args.skip_polar:
    run('polarization',lambda:polarization(symbols,a,b,fs,args.polar_params,ra,axis))
  if not args.skip_igm:
    run('igm',lambda:igm(symbols,a,b,fs,args.density_library,args.igm_step,args.igm_padding,it,
      args.igm_cutoff,args.igm_max_pairs,args.igm_max_points,args.igm_refine,
      mode=args.igm_mode,refine_top=args.igm_refine_top,exclude_depth=args.igm_exclude_depth,
      include_coordination=args.igm_include_coordination,bond_file=args.igm_bonds,ref_ids=ra,
      rtol=args.igm_rtol,atol=args.igm_atol,delta_atol=args.igm_delta_atol,
      pad_extra=args.igm_padding_extra,safety=args.igm_safety))
  bad=[k for k,v in modules.items() if v['status']!='OK']
  result['run_status']='PARTIAL' if bad else 'COMPLETE'
  result['unavailable_or_withheld']=bad
  print('\n'+'SUMMARY'.center(74,'='))
  print('Run status:',result['run_status'],'; no combined stability/reactivity score is assigned.')
  if bad:
    print('Unavailable/withheld:',', '.join(bad))
  if args.csv_prefix:
    export(args.csv_prefix,result,args.overwrite)
    print('JSON and populated CSV tables exported with prefix:',args.csv_prefix)
  return 3 if bad and args.strict else 0


if __name__=='__main__':
  try:
    sys.exit(main())
  except (ValueError,OSError,ImportError,EOFError) as exc:
    print(f'quantmol.py: {exc}',file=sys.stderr)
    sys.exit(2)
