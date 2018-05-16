import numpy as np
import scipy.optimize
from pymol import cmd
from .setbetahelix import helicize_beta_peptide
from ..utils import get_atom_parameters
from typing import Union, Iterable
import itertools


def get_hbdists(angles, resmin, resmax, resistep, selection='all', helicize=True):
    if helicize:
        helicize_beta_peptide(angles, selection)
    dists = []
    for r, rstep in zip(range(resmin, resmax + 1), itertools.cycle(resistep)):
        if rstep == 0: continue
        if r<resmin or r+rstep<resmin or r>resmax or r+rstep>resmax: continue
        selections = ['({}) and resi {} and name HN'.format(selection, r),
                      '({}) and resi {} and name O'.format(selection, r + rstep)]
        if not all([cmd.count_atoms(sel) for sel in selections]):
            continue
        dists.append(cmd.get_distance(*selections))
    return np.array(dists)

def get_dadists(angles, resmin, resmax, resistep, selection='all', helicize=True):
    if helicize:
        helicize_beta_peptide(angles, selection)
    dists = []
    for r, rstep in zip(range(resmin, resmax + 1), itertools.cycle(resistep)):
        if rstep == 0: continue
        if r<resmin or r+rstep<resmin or r>resmax or r+rstep>resmax: continue
        selections = ['({}) and resi {} and name N'.format(selection, r), # donor
                      '({}) and resi {} and name O'.format(selection, r + rstep)] # acceptor
        if not all([cmd.count_atoms(sel) for sel in selections]):
            continue
        dists.append(cmd.get_distance(*selections))
    return np.array(dists)


def get_non_collinearity(angles, resmin, resmax, resistep, selection='all', helicize=True):
    if helicize:
        helicize_beta_peptide(angles, selection)
    angles = []
    for r, rstep in zip(range(resmin, resmax + 1), itertools.cycle(resistep)):
        if rstep==0: continue
        if r<resmin or r+rstep<resmin or r>resmax or r+rstep>resmax: continue
        selections = ['({}) and resi {} and name N'.format(selection, r),
                      '({}) and resi {} and name HN'.format(selection, r),
                      '({}) and resi {} and name O'.format(selection, r + rstep),
                      '({}) and resi {} and name C'.format(selection, r + rstep)]
        if any([cmd.count_atoms(sel) == 0 for sel in selections]):
            continue
        coords = [cmd.get_coords(sel) for sel in selections]
        vec1 = coords[1] - coords[0]  # N->HN vector
        vec2 = coords[2] - coords[3]  # C->O vector
        cosphi = (vec1 * vec2).sum() / (vec1 ** 2).sum() ** 0.5 / (vec2 ** 2).sum() ** 0.5
        angles.append(np.arccos(cosphi) * 180 / np.pi)
    return angles


def get_hbangles(angles, resmin, resmax, resistep, selection='all', helicize=True):
    if helicize:
        helicize_beta_peptide(angles, selection)
    angles = []
    for r, rstep in zip(range(resmin, resmax + 1), itertools.cycle(resistep)):
        if rstep == 0: continue
        if r<resmin or r+rstep<resmin or r>resmax or r+rstep>resmax: continue
        selections = ['({}) and resi {} and name N'.format(selection, r),
                      '({}) and resi {} and name HN'.format(selection, r),
                      '({}) and resi {} and name O'.format(selection, r + rstep)]
        if any([cmd.count_atoms(sel) == 0 for sel in selections]):
            continue
        angles.append(cmd.get_angle(*selections))
    return angles


def targetfunc(angles, resmin, resmax, resistep, hbondlength, selection):
    angles_ = [(angles[3*i+0], angles[3*i+1], angles[3*i+2]) for i in range(len(angles)//3)]
    dists = get_hbdists(angles_, resmin, resmax, resistep, selection, helicize=True)
    hbangles = get_hbangles(angles_, resmin, resmax, resistep, selection, helicize=False)
    if any([a<120 for a in hbangles]):
        return 1000
    val = np.sqrt(np.sum(dists - hbondlength) ** 2) ** 0.5
    #    val = np.std(dists)
    #    print(val)
    val = np.max((dists - hbondlength) ** 2) ** 0.5
    return val


def targetfunc_hbondlen(hbondlen, init, resmin, resmax, resistep, selection):
    result = scipy.optimize.minimize(targetfunc, np.array(init), args=(resmin, resmax, resistep, hbondlen[0], selection),
                                     method='Nelder-Mead')
    return targetfunc(result.x, resmin, resmax, resistep, hbondlen[0], selection)


def optimize_beta_helix(selection, initphi:float, inittheta:float, initpsi:float,
                        resistep:Union[int,Iterable[int]], hbondlen:float):
    """
    DESCRIPTION

        Try to find optimum torsion angles of a beta-peptide helix

    USAGE

        optimize_beta_helix selection, initphi, inittheta, initpsi, resistep, hbondlen

    ARGUMENTS

        selection: the selection to operate on. Typically the beta backbone, including the
            amide oxygen and hydrogens

        initphi: the initial value of phi (degrees)

        inittheta: the initial value of theta (degrees)

        initpsi: the initial value of psi (degrees)

        resistep: the number of residues a hydrogen bond crosses. I.e. the hydrogen of the
            i-th residue is bonded to the oxygen of the i+<resistep>-th residue. Can be
            negative of course. Alternatively, it can be a list of integers, which will be
            applied consecutively for the amino acid components, cycling if exhausted.

        hbondlen: the desired hydrogen bond length.
    """
    def print_summary(values, label):
        print('   {}: {}'.format(label, ', '.join(['{:.2f}'.format(v) for v in values])))
        print('       {:.2f} \xb1 {:.2f} ({:.2f} to {:.2f})'.format(np.mean(values), np.std(values), np.min(values),
                                                                    np.max(values)))
    atompars = get_atom_parameters(selection)
    resimin = min([a['resv'] for a in atompars])
    resimax = max([a['resv'] for a in atompars])
    if isinstance(resistep, str):
        resistep=resistep.strip()
        if (resistep.startswith('(') and resistep.endswith(')')) or (resistep.startswith('[') and resistep.endswith(']')):
            resistep = [int(x.strip()) for x in resistep[1:-1].split(',')]
        else:
            resistep = [int(resistep)]
    elif isinstance(resistep, int):
        resistep = [resistep]
    init = [float(initphi), float(inittheta), float(initpsi)]*len(resistep)
    hbondlen = float(hbondlen)
    result = scipy.optimize.minimize(targetfunc, np.array(init), args=(resimin, resimax, resistep, hbondlen, selection),
                                     method='Nelder-Mead')
    badness = targetfunc(result.x, resimin, resimax, resistep, hbondlen, selection)
    angles_ = [(result.x[3*i+0], result.x[3*i+1], result.x[3*i+2]) for i in range(len(result.x)//3)]
    hbangles = get_hbangles(angles_, resimin, resimax, resistep, selection, helicize=False)
    dists = get_hbdists(angles_, resimin, resimax, resistep, selection, helicize=False)
    print('L_hb={:.3f} A: badness: {:.3f} Angles: {}'.format(
        hbondlen, badness,
        ', '.join(['({:.2f}, {:.2f}, {:.2f})'.format(*a) for a in angles_])))
    print_summary(dists, 'Distances')
    print_summary(hbangles, 'Hydrogen bond angles')
    return angles_
