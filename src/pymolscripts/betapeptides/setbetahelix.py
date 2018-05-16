from __future__ import print_function
import logging
from pymol import cmd
import collections, numbers, itertools

logger=logging.getLogger(__name__)

# torsion angles for various helices taken from:
# Beke et al., Journal of Computational Chemistry 2006, vol. 27, no. 1 pp. 20-38, DOI 10.1002/jcc.20299

helixtypes = [
    ('Z6M', 'g+[g+]a', (126.7, 62.6, 152.7), 'B3LYP/6-311++G(d,p)'),
    ('Z6P', 'g-[g-]a', (-126.7, -62.6, -152.7), 'B3LYP/6-311++G(d,p)'),
    ('Z8M', 'g+[g+]g-', (47.5, 53.5, -104.3), 'RHF/3-21G'),
    ('Z8P', 'g-[g-]g+', (-47.5, -53.5, 104.3), 'RHF/3-21G'),
    ('H8M', 'g+[a]g+', (76.8, -120.6, 52.7), 'B3LYP/6-311++G(d,p)'),
    ('H8P', 'g-[a]g-', (-76.8, 120.6, -52.7), 'B3LYP/6-311++G(d,p)'),
    ('H10M', 'g-[g-]g-', (-77.5, -51.8, -75.1), 'B3LYP/6-311++G(d,p)'),
    ('H10P', 'g+[g+]g+', (77.5, 51.8, 75.1), 'B3LYP/6-311++G(d,p)'),
    ('H12M', 'g+[g-]g+', (92.3, -90.0, 104.6), 'B3LYP/6-311++G(d,p)'),
    ('H12P', 'g-[g+]g-', (-92.3, 90.0, -104.6), 'B3LYP/6-311++G(d,p)'),
    ('H14M', 'a[g-]a', (-140.3, 66.5, -136.8), 'B3LYP/6-311++G(d,p)'),
    ('H14P', 'a[g+]a', (140.3, -66.5, 136.8), 'B3LYP/6-311++G(d,p)'),
    ('SM', 'g+[a]a', (70.5, 176.2, 168.9), 'RHF/3-21G'),
    ('SP', 'g-[a]a', (-70.5, -176.2, -168.9), 'RHF/3-21G'),
]


def set_beta_helix(prevC, N, CB, CA, C, nextN, helixtype, selection='all'):
    """Set torsion angles of a beta peptide backbone

    Inputs:
       prevC: selection of the amide carbon in the previous amino acid
       N: selection of the amide Nitrogen in this amino acid
       CB: selection of the beta carbon
       CA: selection of the alpha carbon
       C: selection of the amide carbon
       nextN: selection of the amide nitrogen of the next amino acid
       helixtype: the type of the helix or a tuple of 3 floats: the 3 angles
       selection: the main selection in which we operate. Defaults to "*"
"""
    try:
        helixparam = [h for h in helixtypes if h[0] == helixtype or h[1] == helixtype][0]
        angles = helixparam[2]
    except IndexError:
        if isinstance(helixtype, str):
            # try to interpret helixtype as a tuple or a list
            helixtype=helixtype.strip()
            try:
                if (helixtype.startswith('(') and helixtype.endswith(')')) or (
                        helixtype.startswith('[') and helixtype.endswith(']')):
                    helixtype = helixtype[1:-1]
                    angles = [float(x.strip()) for x in helixtype.split(',')]
            except:
                raise ValueError('Unknown helix type: {}'.format(helixtype))
        elif (isinstance(helixtype, collections.Sequence) and
            all([isinstance(x, numbers.Real) for x in helixtype])):
            angles = helixtype
        else:
            raise ValueError('Unknown helix type: {}'.format(helixtype))

    print('Helixtype in set_beta_helix: {} (type: {})'.format(helixtype, type(helixtype)))
    atoms = ['({}) and ({})'.format(selection, atomidx) for atomidx in [prevC, N, CB, CA, C, nextN]]
    for i, angle in enumerate(angles):
        cmd.set_dihedral(*(atoms[i:i + 4] + [angle]))
    cmd.delete('pk1')
    cmd.delete('pk2')
    cmd.delete('pkbond')
    cmd.delete('pkmol')


helixdocs = "\n    Known helix types:\n"
helixdocs += '      Perczel short name | IUPAC name |  phi   | theta  |   psi  |    theory level     \n'
helixdocs += '      -------------------+------------+--------+--------+--------+---------------------\n'
for ht in helixtypes:
    helixdocs += '      {:<19s}| {:<10s} | {:=6.1f} | {:=6.1f} | {:=6.1f} | {:<11s}\n'.format(ht[0], ht[1], ht[2][0],
                                                                                              ht[2][1], ht[2][2], ht[3])
set_beta_helix.__doc__ += helixdocs


def helicize_beta_peptide(helixtype, selection='all'):
    """
    DESCRIPTION

    Adjust the torsion angles of a beta-peptide to different helical conformations

    USAGE

    helicize_beta_peptide helixtype [, selection]
    
    ARGUMENTS
    
    helixtype = the type of the helix (either short or IUPAC name),
        or a tuple of 3 floats, representing three torsional angles, or
        a list of tuples / short names / IUPAC names.
    
    selection = the selection to operate on. Must be a single peptide chain with 
        unique residue IDs (default: all)
    
    NOTES"""

    if isinstance(helixtype, str):
        for perczelname, iupacname, angles, theorylevel in helixtypes:
            helixtype = helixtype.replace(iupacname, perczelname).replace(perczelname, '({}, {}, {})'.format(*angles))
        helixtype=helixtype.strip()
        if not all([h in '0123456789.,()[] -+efg' for h in helixtype]):
            raise ValueError('Helixtype parameter contains an invalid character (only numbers, parentheses, brackets, space and commas are accepted)')
        helixtype = eval(helixtype, {}, {})
    assert isinstance(helixtype, collections.Iterable)
    if all([isinstance(x, numbers.Real) for x in helixtype]) and len(helixtype) == 3:
        helixtype = [helixtype]
    print('Helixtypes: {}'.format(helixtype))
    space = {'lis': []}
    cmd.iterate(selection, 'lis.append(resv)', space=space)
    residues = sorted(set(space['lis']))
    for r, ht in zip(sorted(residues), itertools.cycle(helixtype)):
        if len(ht)!=3 and not all([isinstance(x, numbers.Real) for x in ht]):
            raise ValueError('Invalid helixtype: {}'.format(ht))
        calpha = '({}) and (name CA) and (resi {})'.format(selection, r)
        cbeta = '({}) and (name CB+CB1) and (resi {})'.format(selection, r)
        c = '({}) and (name C) and (resi {})'.format(selection, r)
        n = '({}) and (name N) and (resi {})'.format(selection, r)
        prevc = '(neighbor ({})) and (name C)'.format(n)
        nextn = '(neighbor ({})) and (name N)'.format(c)
        prevo = '(neighbor ({})) and (name O)'.format(prevc)
        hn = '(neighbor ({})) and (name H+HN)'.format(n)
        o = '(neighbor ({})) and (name O)'.format(c)
        nexthn = '(neighbor ({})) and (name H+HN)'.format(nextn)
        for name, sel in [('CA', calpha), ('CB', cbeta), ('C', c), ('N', n), ('prevC', prevc), ('nextN', nextn)]:
            cnt = cmd.count_atoms(sel)
            if cnt != 1:
                #logger.warning('Error in residue {}: number of {} atoms found is {}'.format(r, name, cnt))
                break
        else:
            set_beta_helix(prevc, n, cbeta, calpha, c, nextn, ht, selection)
        for n_, h_, c_, o_ in [
            (n, hn, prevc, prevo),
            (nextn, nexthn, c, o)
        ]:
            if cmd.count_atoms(n_) + cmd.count_atoms(h_) + cmd.count_atoms(c_) + cmd.count_atoms(o_) == 4:
                cmd.set_dihedral(h_,n_,c_,o_,180.)
    cmd.orient(selection)


helicize_beta_peptide.__doc__ += helixdocs
