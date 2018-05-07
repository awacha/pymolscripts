from pymol import cmd

from ..utils import iterate_indices

def unbond_close_hydrogen_bonds(selection='all'):
    """
    DESCRIPTION

        Removes spurious chemical bonds between an amide hydrogen and an amide
        oxygen in a peptide/protein

    USAGE

        unbond_close_hydrogen_bonds [selection]

    ARGUMENTS

        selection: the selection to operate on, containing the acceptors and
        the hydrogens (defaults to 'all')

    NOTES

        Some molecular formats do not save bond information. In this case Pymol
        tries to guess bonds based on inter-atomic distances, which sometimes
        results in spurious bonds between atoms which lie too close. This
        function removes these bonds between a hydrogen bond acceptor and a
        hydrogen of a peptide bond.
    """

    for oxygen in iterate_indices('({}) and (e. O) and name O'.format(selection)):
        for hydrogen in iterate_indices('({}) and (neighbor idx {}) and name HN and e. H'.format(selection, oxygen)):
            cmd.unbond('({}) and idx {}'.format(selection, oxygen), '({}) and idx {}'.format(selection, hydrogen))

def find_hbonds(selection_donor_and_hydrogen, selection_acceptor=None, dmin:float=1., dmax:float=2.5, anglemin:float=135):
    dmin = float(dmin)
    dmax = float(dmax)
    anglemin = float(anglemin)
    selection_hydrogen = '({}) and (e. H)'.format(selection_donor_and_hydrogen)
    found_bonds = []
    if selection_acceptor is None:
        selection_acceptor = selection_donor_and_hydrogen
    for acceptor in iterate_indices('({}) and (e. N+O)'.format(selection_acceptor)):
        print('Trying acceptor idx {}'.format(acceptor))
        for hydrogen in iterate_indices('({}) and (e. H) and not (neighbor (idx {} and ({})))'.format(selection_hydrogen, acceptor, selection_acceptor)):
            donors = list(iterate_indices('neighbor ((idx {}) and ({})) and not ((idx {}) and ({}))'.format(hydrogen, selection_hydrogen, acceptor, selection_acceptor)))
            donor=donors[0]
            dist = cmd.get_distance('(idx {}) and ({})'.format(hydrogen, selection_hydrogen),
                                    '(idx {}) and ({})'.format(acceptor, selection_acceptor))
            angle = cmd.get_angle('(idx {}) and ({})'.format(donor, selection_donor_and_hydrogen),
                                  '(idx {}) and ({})'.format(hydrogen, selection_hydrogen),
                                  '(idx {}) and ({})'.format(acceptor, selection_acceptor))
            if dist>=dmin and dist<=dmax and angle>=anglemin:
                found_bonds.append((donor, hydrogen, acceptor, dist, angle))

    for donor, hydrogen, acceptor, dist, angle in found_bonds:
        yield (acceptor, hydrogen, dist)

def generate_hbond_constraints(selection, filename, dmin=1, dmax=2.5, anglemin=135):
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
            f.write('{:10d}{:10d} 2 {:10.4f}\n'.format(idx1, idx2, dist / 10.))

def generate_hbond_restraints_piecewise(selection, filename, strength=1000, mindist=1.1, maxdist=2.5, anglemin=130):
    """
    DESCRIPTION

        Generate piecewise distance restraints for hydrogen bonds

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
            f.write('{:10d}{:10d} 10 {:10.4f} {:10.4f} {:10.4f} {:.6f} ; original distance: {:.6f}\n'.format(
                idx1, idx2, mindist / 10., maxdist/10., maxdist1/10., strength, dist/10.))

def generate_hbond_restraints_harmonic(selection, filename, strength=1000, distance=1.8, mindist_detect=1.1, maxdist_detect=2.5, anglemin_detect=130):
    """
    DESCRIPTION

        Generate harmonic distance restraints for hydrogen bonds

    USAGE

        generate_hbond_restraints_harmonic selection, filename [, strength [, distance [, mindist_detect [, maxdist_detect [, anglemin_detect ]]]]]]

    ARGUMENTS

        selection: the selection to operate on, containing the donors, the acceptors and the hydrogens

        filename: the file name to write the constraints to (a GROMACS .itp file)

        strength: bond strength (kJ mol-1 nm-2) (default: 1000)

        distance: the hydrogen-acceptor distance to restrain to (default: 1.8 A)

        mindist_detect: minimum hydrogen-acceptor distance to consider (default: 1.1 A)

        maxdist_detect: maximum hydrogen-acceptor distance to consider (default: 2.5 A)

        anglemin_detect: minimum donor-hydrogen-acceptor angle to consider (default: 100°)

    NOTES

        A GROMACS .itp file will be written. The restraint potential has the form:

        V(r) = strength * ( r - distance )^2

    """
    mindist=float(mindist_detect)
    maxdist=float(maxdist_detect)
    anglemin = float(anglemin_detect)
    strength = float(strength)
    distance = float(distance)
    with open(filename, 'wt') as f:
        f.write('[ bonds ]\n')
        for idx1, idx2, dist in find_hbonds(selection, selection, mindist, maxdist, anglemin):
            f.write('{:10d}{:10d} 6 {:10.4f} {:.6f} ; original distance: {:.6f}\n'.format(
                idx1, idx2, distance/10., strength, dist/10.))


def beta_hbonds(selection_hydrogen, selection_acceptor, dmin=1, dmax=3):
    i = 0
    for o, h, d in find_hbonds(selection_hydrogen, selection_acceptor, dmin, dmax):
        cmd.distance('dist{:04d}'.format(i), '(idx {}) and ({})'.format(o, selection_acceptor),
                     '(idx {}) and ({})'.format(h, selection_hydrogen), mode=0)
        cmd.group('hbonds', 'dist{:04d}'.format(i))
        i += 1
