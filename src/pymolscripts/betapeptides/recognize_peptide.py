from __future__ import print_function

import logging
import os
import sys

import networkx as nx
from chempy import Atom, Bond
from pymol import cmd

from ..utils import iterate_indices, iterate_neighbours

INVALID_RESID = 0

logger = logging.getLogger(__name__)


def find_peptide_bonds(selection):
    logger.debug('Finding peptide bonds in selection {}:'.format(selection))
    for idx in iterate_indices('({}) and (e. N)'.format(selection)):
        logger.debug('Candidate N: {}'.format(idx))
        # this nitrogen should have:
        # - at least one hydrogen neighbour
        # - one carbon neighbour such as it has exactly one oxygen
        #      neighbour which has no other neighbours
        hydrogens = list(iterate_indices('(neighbor (idx {})) and (e. H) and ({})'.format(idx, selection)))
        if not hydrogens:
            # check if idx is part of a proline ring
            if not ((cmd.count_atoms('(e. C) and (byring idx {}) and ({})'.format(idx, selection)) == 4) and
                    (cmd.count_atoms('(byring idx {}) and ({})'.format(idx, selection)) == 5)):
                # idx is not part of a 5-ring with 4 other carbon atoms
                continue
            logger.debug('Idx {} is a proline nitrogen'.format(idx))
            hydrogens = list(
                iterate_indices('(neighbor (idx {0})) and (byring idx {0}) and ({1})'.format(idx, selection)))
            # now hydrogens contains two CARBON atoms!!!
        carbon = None
        oxygen = None
        for c in iterate_indices('(neighbor (idx {})) and (e. C) and ({})'.format(idx, selection)):
            logger.debug('Candidate C: {}'.format(c))
            for o in iterate_indices('(neighbor (idx {})) and (e. O) and ({})'.format(c, selection)):
                logger.debug('Candidate O: {}'.format(o))
                oneighbours = list(iterate_neighbours(o, selection))
                if len(oneighbours) == 1:
                    logger.debug('Carbon is {}, oxygen is {}'.format(c, o))
                    carbon = c
                    oxygen = o
                    break
            else:
                logger.debug('No good oxygens for this C')
                continue
            break
        else:
            logger.debug('No good carbons for this N.')
            continue
        if carbon is None or oxygen is None:
            logger.debug('No carbon or no oxygen -> no N')
            continue
        for h in hydrogens:
            dih = cmd.get_dihedral('(idx {}) and ({})'.format(h, selection),
                                   '(idx {}) and ({})'.format(idx, selection),
                                   '(idx {}) and ({})'.format(carbon, selection),
                                   '(idx {}) and ({})'.format(oxygen, selection))
            if abs(dih) > 140 or abs(dih) < 40:
                hydrogen = h
                break
        else:
            logger.debug('No appropriate (planar) hydrogens for this nitrogen, oxygen and carbon.')
            continue
        if (len(hydrogens) > 1) and not cmd.count_atoms('(e. C) and (idx {}) and ({})'.format(hydrogen, selection)):
            # if the nitrogen has more than one hydrogens and the trans hydrogen is not a carbon
            # (i.e. this peptide bond does not belong to a proline)
            continue
        logger.debug('Found peptide bond: {}, {}, {}, {}'.format(hydrogen, idx, carbon, oxygen))
        yield (hydrogen, idx, carbon, oxygen)


def get_param(idx, param, selection='all'):
    space = {'lis': []}
    cmd.iterate('(idx {}) and ({})'.format(idx, selection), 'lis.append({})'.format(param), space=space)
    return space['lis'][0]


def get_resv(idx, selection='all'):
    return get_param(idx, 'resv', selection)


def get_chain(idx, selection='all'):
    return get_param(idx, 'chain', selection)


def get_elem(idx, selection='all'):
    return get_param(idx, 'elem', selection)


def number_residues(selection):
    """
    DESCRIPTION

    Number residues in a peptid chain

    USAGE

    number_residues selection

    ARGUMENTS

    selection = a selection containing the peptide chain
    """
    logger.debug('number_residues in selection "{}"'.format(selection))
    cmd.alter(selection, 'resv={}'.format(0))
    peptide_bonds = list(find_peptide_bonds(selection))
    logger.debug('Number of peptide bonds found: {}'.format(len(peptide_bonds)))
    for i, n in enumerate([n for h, n, c, o in peptide_bonds]):
        cmd.alter('(idx {}) and ({})'.format(n, selection), 'resv={}'.format(i + 1))
    for i, n in enumerate([n for h, n, c, o in peptide_bonds]):
        flood_fill_resi(n, get_resv(n, selection), [c for h, n, c, o in peptide_bonds], selection)
    # now set the residue index of the amide bond C=O-s
    for h, n, c, o in peptide_bonds:
        for idx in iterate_neighbours(c, selection):
            if not cmd.count_atoms('(idx {}) and ({}) and (idx {} or idx {})'.format(idx, selection, n, o)):
                # this is the alpha carbon of the current residue
                cmd.alter('(idx {} or idx {}) and ({})'.format(c, o, selection),
                          'resv={}'.format(get_resv(idx, selection)))
    cmd.sort(selection)
    # now residue "0" will be the n-terminal residue. Adjust the numbers to be consecutive
    residue_order = []
    pairs = [(get_resv(c, selection), get_resv(n, selection)) for h, n, c, o in peptide_bonds]
    residue_order = list(pairs[0])
    while True:
        # try to add the next number at the end
        try:
            nextpair = [p for p in pairs if p[0] == residue_order[-1]][0]
            residue_order.append(nextpair[1])
            continue
        except IndexError:
            pass
        try:
            prevpair = [p for p in pairs if p[1] == residue_order[0]][0]
            residue_order.insert(0, prevpair[0])
            continue
        except IndexError:
            pass
        break

    indices_for_residues = [
        list(iterate_indices('({}) and (resi {})'.format(selection, resi)))
        for resi in range(len(peptide_bonds) + 1)]
    for newresi, oldresi in enumerate(residue_order):
        sel = '({}) and '.format(selection) + '(idx ' + '+'.join([str(i) for i in indices_for_residues[oldresi]]) + ')'
        logger.debug('SEL:   {}'.format(sel))
        cmd.alter(sel, 'resv={}'.format(newresi))
    cmd.sort(selection)
    return list(range(len(peptide_bonds) + 1))


def flood_fill_resi(idx, resi, untouchable=None, selection='all'):
    if untouchable is None:
        untouchable = []
    for neighbour in iterate_neighbours(idx, selection):
        if neighbour in untouchable:
            continue
        if get_resv(neighbour, selection) == INVALID_RESID:
            cmd.alter('({}) and (idx {})'.format(selection, neighbour), 'resv={}'.format(resi))
            flood_fill_resi(neighbour, resi, untouchable, selection)
        else:
            pass


def number_chains(selection='all', numbers='ABCDEFGHIJKLMNOPQRSTUVWXYZ'):
    """
    DESCRIPTION

    Find chains and number them consecutively

    USAGE

    number_chains [selection [, numbers]]

    ARGUMENTS

    selection = the selection in which this function operates

    numbers = list of symbols to be applied. Defaults to capital letters of the ABC
    """
    cmd.alter(selection, 'chain=""')
    chains = iter(numbers)
    foundchains = []
    logger.debug('Finding chains among {} atoms'.format(cmd.count_atoms(selection)))
    for idx in iterate_indices(selection):
        if get_chain(idx, selection):
            # if this atom already has a chain ID set, continue with the next atom.
            continue
        foundchains.append(next(chains))
        cmd.alter('(bymol idx {}) and ({})'.format(idx, selection), 'chain="{}"'.format(foundchains[-1]))
    cmd.sort(selection)
    logger.debug('Found chains: ' + ', '.join([str(c) for c in foundchains]))
    return foundchains


def recognize_peptide(rtpfile, selection='all'):
    """
    DESCRIPTION

    Analyze a model, assign chain and residue IDs and name atoms according to
    a GROMACS forcefield

    USAGE

    recognize_peptide rtpfile [, selection]

    ARGUMENTS

    rtpfile = a GROMACS residue database (.rtp file) or a GROMACS forcefield directory

    selection = selected group of atoms to operate on. Defaults to 'all'.

    NOTES

    First the selection is divided into chains according to connectivity. Each
    chain is subsequently divided into residues by finding peptide bonds. The
    residues are numbered consecutively from 0, starting at the N terminus.
    Lastly, every residue of every chain is compared to the entries in the rtp
    file and if matching is found, the residue name and atom names are updated.

    CAVEATS

    The residue matching is done using graph isomorphism (VF2 algorithm, as
    implemented by the networkx package), taking only bonds and element symbols
    into account. This has the following consequences:

    1) There is no distinction between ligands of the same element. E.g. the
        naming of the two hydrogens bound to the nitrogen of the amide group at
        the end of the sidechain of Gln and Asn is undefined, although they
        have slightly different partial charges in the Charmm36m forcefield.
        This has to be fixed by the user.

    2) If the RTP file contains isomorphous residues, the first one will be
        matched. A typical example for this is stereoisomeric amino acids:
        e.g. ALA and DALA, GLN and DGLN in the Charmm36m forcefield. This also
        needs to be fixed by the user.
    """
    chains = number_chains(selection)
    rtpdata = list(graphs_from_rtp(rtpfile))
    for ch in chains:
        print('Looking at chain {}'.format(ch))
        residues = number_residues('({}) and (chain {})'.format(selection, ch))
        print('  Found residues {} to {}'.format(min(residues), max(residues)))
        for r in residues:
            # analyze each residue
            resn = match_amino_acid('({}) and (chain {}) and (resi {})'.format(selection, ch, r), rtpdata)
            if resn is not None:
                print('  Residue {}/{}/ {}'.format(ch, r, resn))
            else:
                print('  Residue {}/{}/ not matched'.format(ch, r))
            if resn in ['B3Q', 'DB3Q']:
                fix_gln_hydrogens('({}) and (chain {}) and (resi {})'.format(selection, ch, r))
    cmd.sort()


def fix_gln_hydrogens(selection):
    """Fix the naming of cis and trans hydrogens in glutamine"""
    # HZ21 is the cis, HZ22 is the trans
    for hindex in iterate_indices('({}) and (name HZ21+HZ22)'.format(selection)):
        if abs(cmd.get_dihedral('(idx {}) and ({})'.format(hindex, selection),
                                '(name NZ2) and ({})'.format(selection),
                                '(name CE) and ({})'.format(selection),
                                '(name OZ1) and ({})'.format(selection))) < 90:
            cmd.alter('(idx {}) and ({})'.format(hindex, selection),
                      'name="HZ21"')
        else:
            cmd.alter('(idx {}) and ({})'.format(hindex, selection),
                      'name="HZ22"')


def select_peptide_bonds(selection:str='all', newselectionprefix:str='pb_'):
    """
    DESCRIPTION

        Select the peptide bonds

    USAGE

        select_peptide_bonds [selection [, newselectionprefix]]

    ARGUMENTS

        selection: selection to operate on, typically a peptide/protein. Defaults to 'all'.

        newselectionprefix: the prefix used for the new selection. Defaults to 'pb_'

    NOTES

        A new selection will be created for each peptide bond in the molecule, named as
        <newselectionprefix><number>, by default pb_0, pb_1, pb_2 etc.
    """
    i = 0
    for h, n, c, o in find_peptide_bonds(selection):
        cmd.select('{}{}'.format(newselectionprefix, i), 'idx {}+{}+{}+{}'.format(h, n, c, o))
        i += 1


def selection_to_graph(selection:str):
    """Make a graph from a selection"""
    G = nx.Graph()
    model = cmd.get_model(selection)
    for atom in model.atom:
        assert isinstance(atom, Atom)
        G.add_node(atom.index, Z=atom.get_number(), idx=atom.index, element=atom.symbol)
        # print('Adding atom {}'.format(atom.index))
    for bond in model.bond:
        assert isinstance(bond, Bond)
        idx1 = model.atom[bond.index[0]].index
        idx2 = model.atom[bond.index[1]].index
        # print('Adding edge between {} and {}'.format(idx1, idx2))
        G.add_edge(idx1, idx2, order=bond.order, idx1=idx1, idx2=idx2)
    return G


def iterate_rtpfiles(rtpfile):
    if sys.version_info.major == 3:
        str_type = str
    elif sys.version_info.major == 2:
        str_type = basestring
    else:
        raise ValueError('Unsupported Python version')
    if os.path.isdir(rtpfile):
        for f in os.listdir(rtpfile):
            if f.endswith('.rtp'):
                yield os.path.join(rtpfile, f)
    else:
        yield rtpfile


def iterate_residues(rtpfile):
    for fn in iterate_rtpfiles(rtpfile):
        with open(fn, 'rt') as f:
            lastresidue = None
            lastheader = None
            atoms = None
            bonds = None
            for l in f:
                if l.strip().startswith(';'):
                    continue
                l = l.strip()
                if not l:
                    continue
                try:
                    l = l[:l.index(';')].strip()
                except ValueError:
                    pass
                if l.startswith('[') and l.endswith(']'):
                    resname = l[1:-1].strip()
                    if resname not in ['bondedtypes', 'atoms', 'bonds', 'impropers', 'cmap']:
                        if lastresidue is not None:
                            yield lastresidue, atoms, bonds
                        lastresidue = resname
                        atoms = []
                        bonds = []
                    else:
                        lastheader = resname
                elif lastheader == 'atoms':
                    atoms.append(l.split()[0])
                elif lastheader == 'bonds':
                    bonds.append(l.split()[:2])
            if lastresidue is not None:
                yield lastresidue, atoms, bonds


def graphs_from_rtp(rtpfile):
    if isinstance(rtpfile, list) and all([isinstance(e, tuple) for e in rtpfile]):
        for resname, graph in rtpfile:
            yield resname, graph
        return

    def generate_graph(atoms, bonds):
        G = nx.Graph()
        for i, a in enumerate(atoms):
            G.add_node(a, idx=i, element=a[0])
        for b in bonds:
            if b[0] in atoms and b[1] in atoms:
                G.add_edge(b[0], b[1], idx1=atoms.index(b[0]), idx2=atoms.index(b[1]))
        return G

    for resname, atoms, bonds in iterate_residues(rtpfile):
        yield resname, generate_graph(atoms, bonds)


def match_amino_acid(selection, rtpfile):
    def nodematch(attr1, attr2):
        return attr1['element'] == attr2['element']

    Gselection = selection_to_graph(selection)
    for resname, Grtp in graphs_from_rtp(rtpfile):
        gm = nx.algorithms.isomorphism.GraphMatcher(Gselection, Grtp, nodematch)
        if gm.is_isomorphic():
            for m in gm.mapping:
                cmd.alter('idx {}'.format(m), 'name="{}"'.format(gm.mapping[m]))
                cmd.alter('idx {}'.format(m), 'resn="{}"'.format(resname))
            return resname
    return None


def select_beta_backbone(selectionname='bbone', originalselection='all'):
    """
    DESCRIPTION

    Select the backbone of beta-peptides

    USAGE

    select_beta_backbone name [, originalselection]

    ARGUMENTS

    name = name of the new selection

    originalselection = superset in which the backbone will be searched
    """
    cmd.select(selectionname, '({}) and name CA+CB+CB1+CC+C+O+N+HN'.format(originalselection))


def order_atoms_in_peptide(rtpfile, selection='all'):
    """
    DESCRIPTION

    Order the atoms in the peptide according to the Gromacs .rtp database

    USAGE

    order_atoms_in_peptide rtpfile [, selection]

    ARGUMENTS

    rtpfile: either a .rtp file (Gromacs residue topology database) or a Gromacs forcefield directory

    selection: defaults to 'all'
    """
    residuetypes = list(iterate_residues(rtpfile))
    residues = {(a.resi_number, a.resn) for a in cmd.get_model(selection).atom}
    i = 0
    for resi, resn in sorted(residues, key=lambda x: x[0]):
        resn, atomnames, bondtypes = [r for r in residuetypes if r[0] == resn][0]
        for a in atomnames:
            print(i, a)
            cmd.alter('({}) and (resi {}) and (resn {}) and (name {})'.format(selection, resi, resn, a),
                      'ID={}'.format(i))
            cmd.alter('({}) and (resi {}) and (resn {}) and (name {})'.format(selection, resi, resn, a),
                      'rank={}'.format(i))
            i += 1
