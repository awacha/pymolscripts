import numpy as np
import scipy.optimize
from pymol import cmd
from .setbetahelix import helicize_beta_peptide
from ..utils import get_atom_parameters


def get_hbdists(angles, resmin, resmax, resistep, selection='all', helicize=True):
    if helicize:
        helicize_beta_peptide('({}, {}, {})'.format(*angles), selection)
    dists = []
    for r in range(resmin, resmax + 1):
        selections = ['({}) and resi {} and name HN'.format(selection, r),
                      '({}) and resi {} and name O'.format(selection, r + resistep)]
        if not all([cmd.count_atoms(sel) for sel in selections]):
            continue
        dists.append(cmd.get_distance(*selections))
    return np.array(dists)

def get_dadists(angles, resmin, resmax, resistep, selection='all', helicize=True):
    if helicize:
        helicize_beta_peptide('({}, {}, {})'.format(*angles), selection)
    dists = []
    for r in range(resmin, resmax + 1):
        selections = ['({}) and resi {} and name N'.format(selection, r), # donor
                      '({}) and resi {} and name O'.format(selection, r + resistep)] # acceptor
        if not all([cmd.count_atoms(sel) for sel in selections]):
            continue
        dists.append(cmd.get_distance(*selections))
    return np.array(dists)


def get_non_collinearity(angles, resmin, resmax, resistep, selection='all', helicize=True):
    if helicize:
        helicize_beta_peptide('({}, {}, {})'.format(*angles), selection)
    angles = []
    for r in range(resmin, resmax + 1):
        selections = ['({}) and resi {} and name N'.format(selection, r),
                      '({}) and resi {} and name HN'.format(selection, r),
                      '({}) and resi {} and name O'.format(selection, r + resistep),
                      '({}) and resi {} and name C'.format(selection, r + resistep)]
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
        helicize_beta_peptide('({}, {}, {})'.format(*angles), selection)
    angles = []
    for r in range(resmin, resmax + 1):
        selections = ['({}) and resi {} and name N'.format(selection, r),
                      '({}) and resi {} and name HN'.format(selection, r),
                      '({}) and resi {} and name O'.format(selection, r + resistep)]
        if any([cmd.count_atoms(sel) == 0 for sel in selections]):
            continue
        angles.append(cmd.get_angle(*selections))
    return angles


def targetfunc(angles, resmin, resmax, resistep, hbondlength, selection):
    #    print('Targetfunc: {}, {}, {}'.format(angles, resmin, resmax))
    dists = get_hbdists(angles, resmin, resmax, resistep, selection, helicize=True)
    hbangles = get_hbangles(angles, resmin, resmax, resistep, selection, helicize=False)
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


def optimize_beta_helix(selection, initphi:float, inittheta:float, initpsi:float, resistep:int, hbondlen:float):
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
            negative of course.

        hbondlen: the desired hydrogen bond length.
    """
    def print_summary(values, label):
        print('   {}: {}'.format(label, ', '.join(['{:.2f}'.format(v) for v in values])))
        print('       {:.2f} \xb1 {:.2f} ({:.2f} to {:.2f})'.format(np.mean(values), np.std(values), np.min(values),
                                                                    np.max(values)))
    init = [float(initphi), float(inittheta), float(initpsi)]
    atompars = get_atom_parameters(selection)
    resimin = min([a['resv'] for a in atompars])
    resimax = max([a['resv'] for a in atompars])
    resistep = int(resistep)
    hbondlen = float(hbondlen)
    result = scipy.optimize.minimize(targetfunc, np.array(init), args=(resimin, resimax, resistep, hbondlen, selection),
                                     method='Nelder-Mead')
    badness = targetfunc(result.x, resimin, resimax, resistep, hbondlen, selection)
    hbangles = get_hbangles(result.x, resimin, resimax, resistep, selection, helicize=False)
    dists = get_hbdists(result.x, resimin, resimax, resistep, selection, helicize=False)
    print('L_hb={:.3f} A: badness: {:.3f} Angles: {:.2f}, {:.2f}, {:.2f}'.format(hbondlen, badness, *result.x))
    print_summary(dists, 'Distances')
    print_summary(hbangles, 'Hydrogen bond angles')
    return list(result.x)

if __name__ == '__main__':
    init = [-140.3, 66.5, -136.8]
    resistep = 2
    # init = [-134.457, 50.916, -136.7757]
    # targetfunc(init, 0,9)
    # result = scipy.optimize.minimize(targetfunc, np.array(init), args=(0, 9, 1.5), method='Nelder-Mead')
    # print(result.x)

    results = []
    for hblen in np.linspace(1.5, 3, 15):
        result = scipy.optimize.minimize(targetfunc, init, args=(0, 9, resistep, hblen), method='Nelder-Mead')
        dists = get_hbdists(result.x, 0, 9, resistep)
        badness = targetfunc(result.x, 0, 9, resistep, hblen)
        coll_angles = get_non_collinearity(result.x, 0, 9, resistep)
        hbangles = get_hbangles(result.x, 0, 9, resistep)
        results.append((hblen, result, dists, badness, coll_angles, hbangles))




    for hblen, result, dists, badness, coll_angles, hbangles in results:
        print('L_hb={:.3f} A: badness: {:.3f} Angles: {:.2f}, {:.2f}, {:.2f}'.format(hblen, badness, *result.x))
        print_summary(dists, 'Distances')
        #    print_summary(coll_angles, 'Collinearity')
        print_summary(hbangles, 'Hydrogen bond angles')

    # result= scipy.optimize.minimize(targetfunc_hbondlen, [2.0], args=(init, 0, 9, resistep))
    # hbondlen = result.x[0]
    # print('Optimum hydrogen bond length:', hbondlen)

    # result = scipy.optimize.minimize(targetfunc, init, args=(0, 9, resistep, hbondlen))
    # print('Optimum torsion angles: ',result.x)
