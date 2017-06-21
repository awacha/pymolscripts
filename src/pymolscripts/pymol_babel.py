from __future__ import print_function
from pymol import cmd
import openbabel

def get_atom_parameters(selection):
    namespace = {'lis':[]}
    exposed_variables = ['model', 'name', 'resn', 'resi', 'resv', 'chain',
                         'alt','elem', 'q', 'b', 'segi', 'type',
                         'formal_charge', 'partial_charge', 'numeric_type',
                         'text_type', 'stereo', 'ID', 'rank', 'index', 'vdw',
                         'ss', 'color', 'reps', 'protons', 's'] # not supported: p
    command='lis.append({'+', '.join(['"{0}": {0}'.format(ev) for ev in exposed_variables])+'})'
    cmd.iterate('({})'.format(selection),command, space=namespace)
    return namespace['lis']


def pymol2ob(selection: str) -> openbabel.OBMol:
    atoms = get_atom_parameters(selection)
    mol = openbabel.OBMol()
    for at in atoms:
        atm = openbabel.OBAtom()
        atm.SetId(at['ID'])
        atm.Set
    
    pass
