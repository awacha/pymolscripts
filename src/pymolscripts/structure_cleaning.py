from pymol import cmd
from .utils import iterate_indices

def select_wrong_bond_numbers(element, maxnbonds, selectionname='sele'):
    """
    DESCRIPTION

        "select_wrong_bond_numbers" selects atoms with superfluous bonds

    USAGE

        select_wrong_bond_numbers element, maxnbonds [, selectionname]

    ARGUMENTS

        element = chemical symbol of the element (C, N, S, etc.)

        maxnbonds = maximum number of bonds formed by this element

        selectionname = the name of the selection to save

    """
    wrong_indices=[]
    for idx in iterate_indices('e. {}'.format(element)):
        neighbours = list(iterate_indices('neighbor idx {}'.format(idx)))
        if len(neighbours) > int(maxnbonds):
            wrong_indices.append(idx)
    cmd.select(selectionname, 'idx '+'+'.join([str(i) for i in wrong_indices]))
    print('Selected {} atoms.'.format(len(wrong_indices)))

cmd.extend('select_wrong_bond_numbers', select_wrong_bond_numbers)