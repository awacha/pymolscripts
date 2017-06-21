from .utils import iterate_indices, iterate_neighbours, Atom
from pymol import cmd
from functools import cmp_to_key

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
            print('Idx {} is a proline nitrogen'.format(idx))
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
        flood_fill_resi(n, Atom(n).resv, [c for h,n,c,o in peptide_bonds])
    cmd.sort(selection)
    # now residue "0" will be the n-terminal residue. Adjust the numbers to be consecutive
    residue_order = []
    pairs = [(Atom(c).resv,Atom(n).resv) for h,n,c,o in peptide_bonds]
    print('Pairs: ')
    for p in pairs:
        print('({}, {})'.format(p[0],p[1]))
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

    print('Residue order: {}'.format(residue_order))
    indices_for_residues = [
        list(iterate_indices('({}) and (resi {})'.format(selection, resi)))
        for resi in range(len(peptide_bonds)+1)]
    for indices, resi in zip(indices_for_residues, residue_order):
        cmd.alter('idx '+'+'.join([str(i) for i in indices]), 'resv={}'.format(resi))
    cmd.sort(selection)
    return list(range(len(peptide_bonds)+1))

cmd.extend('number_residues', number_residues)

def flood_fill_resi(idx, resi, untouchable=None):
    if untouchable is None:
        untouchable = []
    for neighbour in iterate_neighbours(idx):
        if neighbour in untouchable:
            continue
        if Atom(neighbour).resv == INVALID_RESID:
            cmd.alter('idx {}'.format(neighbour), 'resv={}'.format(resi))
            #print('RESI of {} is now {}'.format(neighbour, Atom(neighbour).resv))
            flood_fill_resi(neighbour, resi)
        else:
            #print('RESI of atom {} is {}, which is not invalid ({})'.format(neighbour, Atom(neighbour).resv, INVALID_RESID))
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
        if Atom(idx).chain:
            # if this atom already has a chain ID set, continue with the next atom.
            continue
        foundchains.append(next(chains))
        cmd.alter('(bymol idx {}) and ({})'.format(idx, selection), 'chain="{}"'.format(foundchains[-1]))
    cmd.sort(selection)
    print('Found chains: '+', '.join([str(c) for c in foundchains]))
    return foundchains

cmd.extend('number_chains', number_chains)

def recognize_peptide(selection='all'):
    chains = number_chains(selection)
    for c in chains:
        residues = number_residues('{} and (chain {})'.format(selection, c))
        cmd.alter('{} and (chain {})'.format(selection, c),'name=""')
        peptide_bonds = find_peptide_bonds('{} and (chain {})'.format(selection, c))
        for h, n,c,o in peptide_bonds:
            cmd.alter('idx {}'.format(n),'name="N"')
            cmd.alter('idx {}'.format(c),'name="C"')
            cmd.alter('idx {}'.format(o),'name="O"')
        for r in residues:
            # analyze each residue
            analyze_residue(selection)


cmd.extend('recognize_peptide', recognize_peptide)

def analyze_residue(selection):
    residue