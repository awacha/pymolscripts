from __future__ import print_function
from .utils import iterate_indices, iterate_neighbours
from pymol import cmd
from chempy import Atom, Bond
import sys
import networkx as nx
INVALID_RESID=0

def find_peptide_bonds(selection):
    for idx in iterate_indices('({}) and (e. N)'.format(selection)):
        # this nitrogen should have:
        # - at least one hydrogen neighbour
        # - one carbon neighbour such as it has exactly one oxygen
        #      neighbour which has no other neighbours
        hydrogens = list(iterate_indices('(neighbor (idx {})) and (e. H)'.format(idx)))
        if not hydrogens:
            # check if idx is part of a proline ring
            if not ((cmd.count_atoms('(e. C) and (byring idx {})'.format(idx)) == 4) and
                (cmd.count_atoms('byring idx {}'.format(idx))==5)):
                # idx is not part of a 5-ring with 4 other carbon atoms
                continue
            #print('Idx {} is a proline nitrogen'.format(idx))
            hydrogens = list(iterate_indices('(neighbor (idx {0})) and (byring idx {0})'.format(idx)))
            # now hydrogens contains two CARBON atoms!!!
        carbon = None
        oxygen = None
        for c in iterate_indices('(neighbor (idx {})) and (e. C)'.format(idx)):
            for o in iterate_indices('(neighbor (idx {})) and (e. O)'.format(c)):
                oneighbours = list(iterate_neighbours(o))
                if oneighbours == [c]:
                    carbon = c
                    oxygen = o
                    break
            else:
                continue
            break
        else:
            continue
        if carbon is None or oxygen is None:
            continue
        for h in hydrogens:
            dih=cmd.get_dihedral('idx {}'.format(h),'idx {}'.format(idx),'idx {}'.format(carbon),'idx {}'.format(oxygen))
            if abs(dih)>160:
                hydrogen = h
                break
        else:
            continue
        if (len(hydrogens) >1) and not cmd.count_atoms('(e. C) and (idx {})'.format(hydrogen)):
            # if the nitrogen has more than one hydrogens and the trans hydrogen is not a carbon
            # (i.e. this peptide bond does not belong to a proline)
            continue
        yield (hydrogen,idx,carbon,oxygen)

def get_param(idx, param):
    space = {'lis':[]}
    cmd.iterate('idx {}'.format(idx), 'lis.append({})'.format(param), space=space)
    return space['lis'][0]

def get_resv(idx):
    return get_param(idx, 'resv')

def get_chain(idx):
    return get_param(idx, 'chain')

def get_elem(idx):
    return get_param(idx, 'elem')

def number_residues(selection):
    """
    DESCRIPTION

    Number residues in a peptid chain

    USAGE

    number_residues selection

    ARGUMENTS

    selection = a selection containing the peptide chain
    """
    cmd.alter(selection, 'resv={}'.format(0))
    peptide_bonds=list(find_peptide_bonds(selection))
    for i, n in enumerate([n for h,n,c,o in peptide_bonds]):
        cmd.alter('idx {}'.format(n),'resv={}'.format(i+1))
    for i,n in enumerate([n for h,n,c,o in peptide_bonds]):
        flood_fill_resi(n, get_resv(n), [c for h,n,c,o in peptide_bonds])
    cmd.sort(selection)
    # now residue "0" will be the n-terminal residue. Adjust the numbers to be consecutive
    residue_order = []
    pairs = [(get_resv(c),get_resv(n)) for h,n,c,o in peptide_bonds]
    residue_order = list(pairs[0])
    while True:
        # try to add the next number at the end
        try:
            nextpair = [p for p in pairs if p[0]==residue_order[-1]][0]
            residue_order.append(nextpair[1])
            continue
        except IndexError:
            pass
        try:
            prevpair = [p for p in pairs if p[1]==residue_order[0]][0]
            residue_order.insert(0, prevpair[0])
            continue
        except IndexError:
            pass
        break

    indices_for_residues = [
        list(iterate_indices('({}) and (resi {})'.format(selection, resi)))
        for resi in range(len(peptide_bonds)+1)]
    for newresi, oldresi in enumerate(residue_order):
        cmd.alter('idx '+'+'.join([str(i) for i in indices_for_residues[oldresi]]), 'resv={}'.format(newresi))
    cmd.sort(selection)
    return list(range(len(peptide_bonds)+1))


def flood_fill_resi(idx, resi, untouchable=None):
    if untouchable is None:
        untouchable = []
    for neighbour in iterate_neighbours(idx):
        if neighbour in untouchable:
            continue
        if get_resv(neighbour) == INVALID_RESID:
            cmd.alter('idx {}'.format(neighbour), 'resv={}'.format(resi))
            flood_fill_resi(neighbour, resi)
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
    foundchains=[]
    for idx in iterate_indices(selection):
        if get_chain(idx):
            # if this atom already has a chain ID set, continue with the next atom.
            continue
        foundchains.append(next(chains))
        cmd.alter('(bymol idx {}) and ({})'.format(idx, selection), 'chain="{}"'.format(foundchains[-1]))
    cmd.sort(selection)
    print('Found chains: '+', '.join([str(c) for c in foundchains]))
    return foundchains

def recognize_peptide(rtpfile, selection='all'):
    """
    DESCRIPTION

    Analyze a model, assign chain and residue IDs and name atoms according to
    a GROMACS forcefield

    USAGE

    recognize_peptide rtpfile [, selection]

    ARGUMENTS

    rtpfile = a GROMACS residue database (.rtp file)

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
    for ch in chains:
        print('Looking at chain {}')
        residues = number_residues('({}) and (chain {})'.format(selection, ch))
        print('Found residues {} to {}'.format(min(residues),max(residues)))
        for r in residues:
            # analyze each residue
            print ('Matching chain {}, residue {}'.format(ch, r))
            match_amino_acid('({}) and (chain {}) and (resi {})'.format(selection, ch, r), rtpfile)
    cmd.sort()


#def selection_from_indices(lis):
#    return 'idx '+'+'.join([str(x) for x in lis])
#
# def name_amino_acid(selection, n, c, o):
#     #print('Naming amino acid in selection: ({}). N: {}, C: {}, O: {}'.format(selection,n,c,o))
#     cmd.alter('idx {}'.format(n),'name="N"')
#     cmd.alter('idx {}'.format(c), 'name="C"')
#     cmd.alter('idx {}'.format(o), 'name="O"')
#     named = [n,c,o]
#     calpha=list(iterate_indices('(neighbor idx {}) and ({}) and not ({}) and (e. C)'.format(
#         c,selection,selection_from_indices(named))))
#     if not len(calpha)==1:
#         raise ValueError('No CA or more C atoms bound to C')
#     calpha=calpha[0]
#     cmd.alter('idx {}'.format(calpha), 'name="CA"')
#     named.append(calpha)
#     distance_found=[n,c,o,calpha]
#     distances_from_calpha={calpha:0}
#     previous_atom_in_chain={}
#     while True:
#         for idx in list(distances_from_calpha.keys()):
#             for nbr in iterate_neighbours(idx):
#                 if nbr in distance_found:
#                     continue
#                 if get_elem(nbr) == 'H':
#                     continue
#                 distances_from_calpha[nbr]=distances_from_calpha[idx]+1
#                 previous_atom_in_chain[nbr]=idx
#                 distance_found.append(nbr)
#         if not cmd.count_atoms('({}) and not ({}) and not (e. H)'.format(selection, selection_from_indices(distance_found))):
#             break
#     for idx in distances_from_calpha:
#         symbols='ABGDEZHTIKLMNPO'
#         same_level_atoms = [i for i in distances_from_calpha if distances_from_calpha[i]==distances_from_calpha[idx]]
#         if len(same_level_atoms) == 1:
#             cmd.alter('idx {}'.format(idx), 'name="{}{}"'.format(get_elem(idx),symbols[distances_from_calpha[idx]]))
#         else: # len(same_level_atoms) cannot be 0, the list contains at least `idx`
#             pass
#             #ToDo
#             # more atoms are at the same distance from CA
#             # Rule 1): check the previous atom
#
#     # cmd.sort() must not be called, it messes up the indices!

def select_peptide_bonds(selection, newselectionprefix='pb_'):
    i=0
    for h,n,c,o in find_peptide_bonds(selection):
        cmd.select('{}{}'.format(newselectionprefix,i), 'idx {}+{}+{}+{}'.format(h,n,c,o))
        i+=1


def selection_to_graph(selection):
    G=nx.Graph()
    model = cmd.get_model(selection)
    for atom in model.atom:
        assert isinstance(atom, Atom)
        G.add_node(atom.index, {'Z':atom.get_number(), 'idx':atom.index, 'element':atom.symbol})
        #print('Adding atom {}'.format(atom.index))
    for bond in model.bond:
        assert isinstance(bond, Bond)
        idx1=model.atom[bond.index[0]].index
        idx2=model.atom[bond.index[1]].index
        #print('Adding edge between {} and {}'.format(idx1, idx2))
        G.add_edge(idx1, idx2, {'order':bond.order, 'idx1':idx1, 'idx2':idx2})
    return G

def graphs_from_rtp(rtpfile):
    def generate_graph(atoms, bonds):
        G=nx.Graph()
        for i,a in enumerate(atoms):
            G.add_node(a, {'idx':i,'element':a[0]})
        for b in bonds:
            if b[0] in atoms and b[1] in atoms:
                G.add_edge(b[0], b[1], {'idx1':atoms.index(b[0]), 'idx2':atoms.index(b[1])})
        return G
    if sys.version_info.major==3:
        str_type=str
    elif sys.version_info.major == 2:
        str_type=basestring
    else:
        raise ValueError('Unsupported Python version')
    if isinstance(rtpfile, str_type):
        rtpfile = [rtpfile]
    for fn in rtpfile:
        with open(fn, 'rt') as f:
            lastresidue = None
            lastheader = None
            atoms = None
            bonds = None
            for l in f:
                if l.strip().startswith(';'):
                    continue
                l=l.strip()
                if not l:
                    continue
                try:
                    l=l[:l.index(';')]
                except ValueError:
                    pass
                if l.startswith('[') and l.endswith(']'):
                    resname = l[1:-1].strip()
                    if resname not in ['bondedtypes', 'atoms', 'bonds', 'impropers', 'cmap']:
                        if lastresidue is not None:
                            yield lastresidue, generate_graph(atoms, bonds)
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
                yield lastresidue, generate_graph(atoms, bonds)

def match_amino_acid(selection, rtpfile):
    def nodematch(attr1, attr2):
        return attr1['element'] == attr2['element']

    Gselection = selection_to_graph(selection)
    #print([n[1] for n in Gselection.nodes(data=True)])
    #labels = {n[0]:n[1]['element'] for n in Gselection.nodes(data=True)}
    #print('Labels:',labels)
    #nx.draw(Gselection, ax=plt.gca(), pos = nx.spring_layout(
    #            Gselection, k=0.01,
    #            iterations=500), with_labels=True, labels = labels)
    #plt.gcf().canvas.draw()
    #plt.show()
    for resname, Grtp in graphs_from_rtp(rtpfile):
        gm = nx.algorithms.isomorphism.GraphMatcher(Gselection, Grtp, nodematch)
        if gm.is_isomorphic():
            print('Matched:',resname)
            for m in gm.mapping:
                #print('  {} -> {}'.format(m,gm.mapping[m]))
                cmd.alter('idx {}'.format(m), 'name="{}"'.format(gm.mapping[m]))
                cmd.alter('idx {}'.format(m), 'resn="{}"'.format(resname))
            return resname
    print('No match.')
    return None


cmd.extend('number_residues', number_residues)
cmd.extend('number_chains', number_chains)
cmd.extend('recognize_peptide', recognize_peptide)
cmd.extend('select_peptide_bonds', select_peptide_bonds)
cmd.extend('match_amino_acid', match_amino_acid)
