import numpy as np
import scipy.optimize
from pymol import cmd
from pymolscripts.setbetahelix import helicize_beta_peptide


def get_dists(angles, resmin, resmax, resistep):
    helicize_beta_peptide('({}, {}, {})'.format(*angles))
    dists = []
    for r in range(resmin, resmax + 1):
        selections = ['resi {} and name HN'.format(r),
                      'resi {} and name O'.format(r + resistep)]
        if not all([cmd.count_atoms(sel) for sel in selections]):
            continue
        dists.append(cmd.get_distance(*selections))
    return np.array(dists)


def get_non_collinearity(angles, resmin, resmax, resistep):
    helicize_beta_peptide('({}, {}, {})'.format(*angles))
    angles = []
    for r in range(resmin, resmax + 1):
        selections = ['resi {} and name N'.format(r),
                      'resi {} and name HN'.format(r),
                      'resi {} and name O'.format(r + resistep),
                      'resi {} and name C'.format(r + resistep)]
        if any([cmd.count_atoms(sel) == 0 for sel in selections]):
            continue
        coords = [cmd.get_coords(sel) for sel in selections]
        vec1 = coords[1] - coords[0]  # N->HN vector
        vec2 = coords[2] - coords[3]  # C->O vector
        cosphi = (vec1 * vec2).sum() / (vec1 ** 2).sum() ** 0.5 / (vec2 ** 2).sum() ** 0.5
        angles.append(np.arccos(cosphi) * 180 / np.pi)
    return angles


def get_hbangles(angles, resmin, resmax, resistep):
    helicize_beta_peptide('({}, {}, {})'.format(*angles))
    angles = []
    for r in range(resmin, resmax + 1):
        selections = ['resi {} and name N'.format(r),
                      'resi {} and name HN'.format(r),
                      'resi {} and name O'.format(r + resistep)]
        if any([cmd.count_atoms(sel) == 0 for sel in selections]):
            continue
        angles.append(cmd.get_angle(*selections))
    return angles


def targetfunc(angles, resmin, resmax, resistep, hbondlength):
    #    print('Targetfunc: {}, {}, {}'.format(angles, resmin, resmax))
    dists = get_dists(angles, resmin, resmax, resistep)
    val = np.sqrt(np.sum(dists - hbondlength) ** 2) ** 0.5
    #    val = np.std(dists)
    #    print(val)
    val = np.max((dists - hbondlength) ** 2) ** 0.5
    return val


def targetfunc_hbondlen(hbondlen, init, resmin, resmax, resistep):
    result = scipy.optimize.minimize(targetfunc, np.array(init), args=(resmin, resmax, resistep, hbondlen[0]),
                                     method='Nelder-Mead')
    return targetfunc(result.x, resmin, resmax, resistep, hbondlen[0])


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
        dists = get_dists(result.x, 0, 9, resistep)
        badness = targetfunc(result.x, 0, 9, resistep, hblen)
        coll_angles = get_non_collinearity(result.x, 0, 9, resistep)
        hbangles = get_hbangles(result.x, 0, 9, resistep)
        results.append((hblen, result, dists, badness, coll_angles, hbangles))


    def print_summary(values, label):
        print('   {}: {}'.format(label, ', '.join(['{:.2f}'.format(v) for v in values])))
        print('       {:.2f} \xb1 {:.2f} ({:.2f} to {:.2f})'.format(np.mean(values), np.std(values), np.min(values),
                                                                    np.max(values)))


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
