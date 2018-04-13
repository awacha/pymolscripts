# coding: utf-8

import math

from pymol import cmd


def save_gro(filename, selection='(all)', bx=None, by=None, bz=None):
    """
    DESCRIPTION
    
        "save_gro" saves content to a .gro file in a similar fashion as "save"

    USAGE

        save_gro filename [, selection [, bx [, by [, bz ]]]]

    ARGUMENTS
        
        filename = string: file path to be written

        selection = string: atoms to save {default: (all)}
        
        bx = float: rectangular box size X dimension in Angströms
        
        by = float: rectangular box size Y dimension in Angströms

        bz = float: rectangular box size Z dimension in Angströms
        
    NOTES

        In contrast to "save", only the current state can be saved
        as the GRO format is a single-structure format.
        
        If any of bx, by and bz is not supplied, the box will be 
        determined from the get_symmetry() command.
        
        Atoms are written with increasing ID.

    SEE ALSO

        load, save

    """
    with open(filename, 'wt') as f:
        f.write('{}\n'.format(selection))
        f.write('{:>5d}\n'.format(cmd.count_atoms(selection)))
        space = {'l': []}
        cmd.iterate_state(-1, selection, 'l.append((ID, resi,resn,name,index,x,y,z))', space=space)
        for ID, resi, resn, name, index, x, y, z in sorted(space['l'], key=lambda a: (a[1], a[4])):
            f.write('{:>5d}{:<5s}{:>5s}{:>5d}{:>8.3f}{:>8.3f}{:>8.3f}\n'.format(int(resi), resn, name, int(ID) + 1,
                                                                                float(x) * 0.1, float(y) * 0.1,
                                                                                float(z) * 0.1))
        if bx is None and by is None and bz is None:
            try:
                a, b, c, alpha, beta, gamma, spacegroup = cmd.get_symmetry(selection)
                alpha *= math.pi / 180.
                beta *= math.pi / 180.
                gamma *= math.pi / 180.
                a /= 10
                b /= 10
                c /= 10
                v1x = a
                v2x = b * math.cos(gamma)
                v2y = b * math.sin(gamma)
                v3x = c * math.cos(beta)
                v3y = c * (math.cos(alpha) - math.cos(gamma) * math.cos(beta)) / (math.sin(gamma))
                v3z = (c ** 2 - v3x ** 2 - v3y ** 2) ** 0.5
                f.write(
                    ' {:>9.5f} {:>9.5f} {:>9.5f} {:>9.5f} {:>9.5f} {:>9.5f} {:>9.5f} {:>9.5f} {:>9.5f}\n'.format(v1x,
                                                                                                                 v2y,
                                                                                                                 v3z, 0,
                                                                                                                 0, v2x,
                                                                                                                 0, v3x,
                                                                                                                 v3y))
            except TypeError:
                f.write(' 0 0 0\n')
        else:
            f.write(' {} {} {}\n'.format(bx, by, bz))
    print('Saved {} atoms.'.format(len(space['l'])))
