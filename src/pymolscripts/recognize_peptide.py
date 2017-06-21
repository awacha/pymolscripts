from .utils import iterate_indices, iterate_neighbours, Atom
from pymol import cmd

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
    cmd.alter(selection, 'resv={}'.format(INVALID_RESID))
    peptide_bonds=list(find_peptide_bonds(selection))
    for i, n in enumerate([n for h,n,c,o in peptide_bonds]):
        cmd.alter('idx {}'.format(n),'resv={}'.format(i+1))
    #for c,o in [(c,o) for h,n,c,o in peptide_bonds]:
    #    cmd.alter('idx {} or idx {}'.format(c,o),'resv={}'.format(INVALID_RESID+100000000))
    for i,n in enumerate([n for h,n,c,o in peptide_bonds]):
        #print(n)
        flood_fill_resi(n, Atom(n).resv, [c for h,n,c,o in peptide_bonds])
    #cmd.alter('resi {}'.format(INVALID_RESID), 'resi=0')

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

def sort_resids(selection):
    pass