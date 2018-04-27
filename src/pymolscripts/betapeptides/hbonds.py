from pymol import cmd

from ..utils import iterate_indices

def find_hbonds(selection_donor_and_hydrogen, selection_acceptor=None, dmin:float=1., dmax:float=2.5, anglemin:float=135):
    dmin = float(dmin)
    dmax = float(dmax)
    anglemin = float(anglemin)
    selection_hydrogen = '({}) and (e. H)'.format(selection_donor_and_hydrogen)
    found_bonds = []
    if selection_acceptor is None:
        selection_acceptor = selection_donor_and_hydrogen
    for acceptor in iterate_indices('({}) and (e. N+O)'.format(selection_acceptor)):
        for hydrogen in iterate_indices('({}) and (e. H)'.format(selection_hydrogen)):
            donor = list(iterate_indices('neighbor ((idx {}) and ({}))'.format(hydrogen, selection_hydrogen)))[0]
            dist = cmd.get_distance('(idx {}) and ({})'.format(hydrogen, selection_hydrogen),
                                    '(idx {}) and ({})'.format(acceptor, selection_acceptor))
            angle = cmd.get_angle('(idx {}) and ({})'.format(donor, selection_donor_and_hydrogen),
                                  '(idx {}) and ({})'.format(hydrogen, selection_hydrogen),
                                  '(idx {}) and ({})'.format(acceptor, selection_acceptor))
            if dist>=dmin and dist<=dmax and angle>=anglemin:
                found_bonds.append((donor, hydrogen, acceptor, dist, angle))

    for donor, hydrogen, acceptor, dist, angle in found_bonds:
        yield (acceptor, hydrogen, dist)

def generate_hbond_constraints(selection, filename, dmin=1, dmax=1, anglemin=135):
    """
    DESCRIPTION

        Generate distance constraints for hydrogen bonds

    USAGE

        generate_hbond_constraints selection, filename [, dmin [, dmax [, anglemin ]]]

    ARGUMENTS

        selection: the selection to operate on, containing the donors, the acceptors and the hydrogens

        filename: the file name to write the constraints to (a GROMACS .itp file)

        dmin: minimum hydrogen-acceptor distance to consider

        dmax: maximum hydrogen-acceptor distance to consider

        anglemin: minimum donor-hydrogen-acceptor angle to consider
    """
    dmin = float(dmin)
    dmax = float(dmax)
    anglemin = float(anglemin)
    with open(filename, 'wt') as f:
        f.write('[ constraints ]\n')
        for idx1, idx2, dist in find_hbonds(selection, selection, dmin, dmax, anglemin):
            f.write('{:10d}{:10d} 2 {:10.4f}\n'.format(idx1[1], idx2[1], dist / 10.))

def generate_hbond_restraints(selection, filename, strength=1000, mindist=1.7, maxdist=2.3, anglemin=135):
    """
    DESCRIPTION

        Generate distance constraints for hydrogen bonds

    USAGE

        generate_hbond_restraints selection, filename [, strength [, mindist [, maxdist [, anglemin ]]]]]

    ARGUMENTS

        selection: the selection to operate on, containing the donors, the acceptors and the hydrogens

        filename: the file name to write the constraints to (a GROMACS .itp file)

        strength: bond strength (kJ mol-1 nm-2)

        mindist: minimum hydrogen-acceptor distance to consider

        maxdist: maximum hydrogen-acceptor distance to consider

        anglemin: minimum donor-hydrogen-acceptor angle to consider

    NOTES

        A GROMACS .itp file will be written. The restraint potential has the form:

                 / 1/2 strength * (r-mindist)^2                                       if r < mindist
                |
                |  0                                                                  if mindist <= r < maxdist
        V(r) = <
                |  1/2 strength * (r-maxdist)^2                                       if maxdist <= r < maxdist1
                |
                 \ 1/2 strength * (maxdist1 - maxdist) * (2*r - maxdist1 - maxdist)   if maxdist1 <= r

        maxdist1 = maxdist + (maxdist-mindist)*0.5 is used.
    """
    mindist=float(mindist)
    maxdist=float(maxdist)
    maxdist1=maxdist+(maxdist-mindist)*0.5
    anglemin = float(anglemin)
    strength = float(strength)
    with open(filename, 'wt') as f:
        f.write('[ bonds ]\n')
        for idx1, idx2, dist in find_hbonds(selection, selection, mindist, maxdist, anglemin):
            f.write('{:10d}{:10d} 10 {:10.4f} {:10.4f} {:10.4f} {:.6f}\n'.format(
                idx1[1], idx2[1], mindist / 10., maxdist/10., maxdist1/10., strength))

def beta_hbonds(selection_hydrogen, selection_acceptor, dmin=1, dmax=3):
    i = 0
    for o, h, d in find_hbonds(selection_hydrogen, selection_acceptor, dmin, dmax):
        cmd.distance('dist{:04d}'.format(i), '(idx {}) and ({})'.format(o, selection_acceptor),
                     '(idx {}) and ({})'.format(h, selection_hydrogen), mode=0)
        cmd.group('hbonds', 'dist{:04d}'.format(i))
        i += 1
