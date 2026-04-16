#!/usr/bin/env python3
"""
Check for self-interaction in KR or RKR simulations (per-molecule count).

Works in the UT (LAMMPS) frame using unwrapped coordinates from dump.

For each molecule:
  - d_unwrap = |xu_j - xu_i| (unwrapped intra-molecular distance)
  - d_minimg = minimum image distance using the triclinic box
  - wrapping: any intra-molecular pair with d_unwrap != d_minimg
  - self-interaction: any intra-molecular pair with d_minimg < rcut
                      AND d_unwrap != d_minimg (wrapped)

Counts are per molecule: a molecule with multiple self-interacting pairs
is counted once.

Options:
  --rcut RCUT         cutoff distance (default 2^(1/6))
  --min-ree REE       only consider molecules with end-to-end distance >= REE
                      (unwrapped, using first and last atom by id)
                      Useful to filter out coiled molecules when studying
                      extensional self-interaction.
  --ref-data FILE     reference equilibrium data file (LAMMPS data format).
                      If given, min-ree is set to 2 * mean(Ree) computed
                      from the reference structure, unless --min-ree is
                      also specified.

Usage:
  python check_self_interaction_mol.py <dump_file> [--rcut 1.122] \
      [--min-ree 30] [--ref-data equilibrium.data]
"""

import sys
import argparse
import numpy as np


def read_dump_frame(f):
    line = f.readline()
    if not line:
        return None
    timestep = int(f.readline().strip())
    f.readline()
    natoms = int(f.readline().strip())
    f.readline()
    box_lines = []
    for _ in range(3):
        box_lines.append(f.readline().strip().split())

    xlo_bound = float(box_lines[0][0])
    xhi_bound = float(box_lines[0][1])
    xy = float(box_lines[0][2]) if len(box_lines[0]) > 2 else 0.0
    ylo_bound = float(box_lines[1][0])
    yhi_bound = float(box_lines[1][1])
    xz = float(box_lines[1][2]) if len(box_lines[1]) > 2 else 0.0
    zlo_bound = float(box_lines[2][0])
    zhi_bound = float(box_lines[2][1])
    yz = float(box_lines[2][2]) if len(box_lines[2]) > 2 else 0.0

    xlo = xlo_bound - min(0.0, xy, xz, xy + xz)
    xhi = xhi_bound - max(0.0, xy, xz, xy + xz)
    ylo = ylo_bound - min(0.0, yz)
    yhi = yhi_bound - max(0.0, yz)
    lx = xhi - xlo
    ly = yhi - ylo
    lz = zhi_bound - zlo_bound

    H = np.array([[lx, xy, xz],
                  [0,  ly, yz],
                  [0,  0,  lz]])

    header = f.readline().strip().split()
    cols = header[2:]
    idx = {c: i for i, c in enumerate(cols)}

    atoms = [f.readline().strip().split() for _ in range(natoms)]

    ids  = np.array([int(a[idx['id']])  for a in atoms])
    mols = np.array([int(a[idx['mol']]) for a in atoms])
    xu   = np.array([float(a[idx['xu']]) for a in atoms])
    yu   = np.array([float(a[idx['yu']]) for a in atoms])
    zu   = np.array([float(a[idx['zu']]) for a in atoms])
    pos = np.column_stack([xu, yu, zu])

    return {'timestep': timestep, 'H': H, 'ids': ids, 'mols': mols, 'pos': pos}


def check_frame(frame, rcut, min_ree=0.0):
    """Per-molecule flags: wrapping / self-interacting.
    If min_ree > 0, only molecules with unwrapped end-to-end distance
    (first atom id vs last atom id) >= min_ree are considered.
    """
    pos = frame['pos']
    mols = frame['mols']
    ids = frame['ids']
    H = frame['H']
    Hinv = np.linalg.inv(H)

    mol_ids = np.unique(mols)
    n_wrap = 0
    n_self = 0
    n_considered = 0
    self_mols = []

    for mid in mol_ids:
        aidx = np.where(mols == mid)[0]
        if len(aidx) < 2:
            continue

        # sort by atom id to identify chain endpoints
        order = np.argsort(ids[aidx])
        mol_pos = pos[aidx][order]

        # end-to-end (unwrapped) filter
        if min_ree > 0.0:
            ree = np.linalg.norm(mol_pos[-1] - mol_pos[0])
            if ree < min_ree:
                continue
        n_considered += 1

        is_wrap = False
        is_self = False

        for i in range(len(mol_pos)):
            for j in range(i + 1, len(mol_pos)):
                dr = mol_pos[j] - mol_pos[i]
                d_unwrap = np.linalg.norm(dr)
                s = Hinv @ dr
                s -= np.round(s)
                d_min = np.linalg.norm(H @ s)

                if abs(d_unwrap - d_min) > 0.1:
                    is_wrap = True
                    if d_min < rcut:
                        is_self = True
                        break
            if is_self:
                break

        if is_wrap:
            n_wrap += 1
        if is_self:
            n_self += 1
            self_mols.append(int(mid))

    return n_considered, n_wrap, n_self, self_mols


def read_lammps_data(path):
    """Minimal LAMMPS data reader for Atoms section.
    Returns (ids, mols, xyz). Supports atom_style that places mol-id in col 2
    (e.g., angle, bond, full). Image flags are applied if present.
    """
    ids = []; mols = []; xyz = []
    with open(path) as f:
        lines = f.read().splitlines()
    # find "Atoms" section
    i = 0
    while i < len(lines) and not lines[i].strip().startswith('Atoms'):
        i += 1
    if i >= len(lines):
        raise RuntimeError(f"'Atoms' section not found in {path}")
    i += 1
    # skip blank lines and optional style comment
    while i < len(lines) and lines[i].strip() == '':
        i += 1
    # box sizes needed for image unwrap
    bx = by = bz = None; xy = xz = yz = 0.0
    for line in lines:
        t = line.strip()
        if t.endswith('xlo xhi'):
            lo, hi = map(float, t.split()[:2]); bx = hi - lo
        elif t.endswith('ylo yhi'):
            lo, hi = map(float, t.split()[:2]); by = hi - lo
        elif t.endswith('zlo zhi'):
            lo, hi = map(float, t.split()[:2]); bz = hi - lo
        elif t.endswith('xy xz yz'):
            xy, xz, yz = map(float, t.split()[:3])
    # parse Atoms lines (until next section header or EOF)
    section_headers = ('Velocities', 'Bonds', 'Angles', 'Dihedrals',
                       'Impropers', 'Masses')
    while i < len(lines):
        ln = lines[i].strip()
        if ln == '' or ln.startswith(section_headers):
            break
        toks = ln.split()
        # atom_style angle: id mol type x y z [ix iy iz]
        aid = int(toks[0]); mid = int(toks[1])
        x = float(toks[3]); y = float(toks[4]); z = float(toks[5])
        ix = int(toks[6]) if len(toks) >= 9 else 0
        iy = int(toks[7]) if len(toks) >= 9 else 0
        iz = int(toks[8]) if len(toks) >= 9 else 0
        # unwrap using image flags (orthorhombic triclinic ok for Ree purposes)
        xu = x + ix*bx + iy*xy + iz*xz
        yu = y + iy*by + iz*yz
        zu = z + iz*bz
        ids.append(aid); mols.append(mid); xyz.append([xu, yu, zu])
        i += 1
    return (np.array(ids), np.array(mols), np.array(xyz))


def mean_ree_from_data(path):
    ids, mols, xyz = read_lammps_data(path)
    rees = []
    for mid in np.unique(mols):
        idx = np.where(mols == mid)[0]
        if len(idx) < 2: continue
        order = np.argsort(ids[idx])
        p = xyz[idx][order]
        rees.append(np.linalg.norm(p[-1] - p[0]))
    if not rees:
        raise RuntimeError(f"No molecules with >= 2 atoms in {path}")
    return float(np.mean(rees))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('dump')
    ap.add_argument('--rcut', type=float, default=2.0**(1.0/6.0))
    ap.add_argument('--min-ree', type=float, default=None,
                    help='only count molecules with Ree >= min-ree (unwrapped)')
    ap.add_argument('--ref-data', default=None,
                    help='equilibrium data file; sets min-ree = 2 * mean(Ree)')
    args = ap.parse_args()

    dump_file = args.dump
    rcut = args.rcut

    if args.min_ree is not None:
        min_ree = args.min_ree
        ref_note = "(user-supplied)"
    elif args.ref_data is not None:
        ref_ree = mean_ree_from_data(args.ref_data)
        min_ree = 2.0 * ref_ree
        ref_note = f"(2 x mean Ree = 2 x {ref_ree:.3f} from {args.ref_data})"
    else:
        min_ree = 0.0
        ref_note = "(no filter)"

    print(f"Dump file: {dump_file}")
    print(f"Cutoff distance: {rcut:.6f}")
    if min_ree > 0:
        print(f"Filter: only molecules with Ree >= {min_ree:.3f} {ref_note}")
    else:
        print(f"Filter: none {ref_note}")
    print()
    print(f"{'Frame':>5} {'Step':>10} {'#considered':>11} {'#wrap':>7} "
          f"{'#self':>7}  self_mols (first 10)")

    totals = {'considered': 0, 'wrap': 0, 'self': 0, 'frames': 0,
              'nmols_max': 0}

    with open(dump_file, 'r') as f:
        frame_idx = 0
        while True:
            frame = read_dump_frame(f)
            if frame is None:
                break
            n_con, n_wrap, n_self, self_mols = check_frame(frame, rcut, min_ree)
            nmols_total = len(np.unique(frame['mols']))
            totals['nmols_max'] = max(totals['nmols_max'], nmols_total)
            print(f"{frame_idx:>5} {frame['timestep']:>10} {n_con:>11} "
                  f"{n_wrap:>7} {n_self:>7}  {self_mols[:10]}")
            totals['considered'] += n_con
            totals['wrap'] += n_wrap
            totals['self'] += n_self
            totals['frames'] += 1
            frame_idx += 1

    print()
    print(f"Total frames: {totals['frames']}")
    print(f"Total molecules: {totals['nmols_max']}")
    if totals['frames'] > 0:
        print(f"Sum over frames: considered = {totals['considered']}, "
              f"wrapping = {totals['wrap']}, self-interacting = {totals['self']}")
        print(f"Average per frame: considered = "
              f"{totals['considered']/totals['frames']:.1f}, "
              f"wrapping = {totals['wrap']/totals['frames']:.1f}, "
              f"self-interacting = {totals['self']/totals['frames']:.2f}")

    if totals['self'] == 0:
        print("No self-interaction detected" +
              (f" (among Ree>={min_ree:.1f} molecules)" if min_ree > 0 else "") +
              ".")
    else:
        print("WARNING: Self-interaction detected!")


if __name__ == '__main__':
    main()
