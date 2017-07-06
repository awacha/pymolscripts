from pymol import cmd

def select_atoms(selection):
    for idx in iterate_indices(selection):
        yield Atom(idx)

class Atom(object):
    _exposed_variables = ['model', 'name', 'resn', 'resi', 'resv', 'chain',
                         'alt','elem', 'q', 'b', 'segi', 'type',
                         'formal_charge', 'partial_charge', 'numeric_type',
                         'text_type', 'stereo', 'ID', 'rank', 'index', 'vdw',
                         'ss', 'color', 'reps', 'protons', 's'] # not supported: p
    def __init__(self, index_or_selection):
        if isinstance(index_or_selection, (int, float)):
            self._index = int(index_or_selection)
        else:
            idx = list(iterate_indices(index_or_selection))
            if len(idx)>1:
                raise ValueError('Selection contains more than one atoms')
            elif not idx:
                raise ValueError('Empty selection')
            else:
                self._index=idx[0]
        self._cache = {}
        #self.refresh()

    def refresh(self):
        namespace = {'lis':[]}
        command='lis.append({'+', '.join(['"{0}": {0}'.format(ev) for ev in self._exposed_variables])+'})'
        cmd.iterate('(idx {})'.format(self._index),command, space=namespace)
        if not namespace['lis']:
            raise ValueError('Atom not found with index {}'.format(self._index))
        for ev in sorted(namespace['lis'][0]):
            self._cache[ev]=namespace['lis'][0][ev]

    def __str__(self):
        s="Atom:\n"
        for ev in sorted(self._exposed_variables):
            s+='  {}: {} (type {})\n'.format(ev,getattr(self,ev),type(getattr(self,ev)))
        return s

    def __getattr__(self, name):
        if name in self._exposed_variables:
            space = {'lis':[]}
            cmd.iterate('(idx {})'.format(self._index), 'lis.append({})'.format(name),space=space)
            self._cache[name]=space['lis'][0]
            return self._cache[name]
        else:
            return super(Atom, self).__getattribute__(name)

    def __setattr__(self, name, value):
        if name in self._exposed_variables:
            cmd.alter('idx {}'.format(self._index), '{}={}'.format(name, value))
        else:
            super(Atom, self).__setattr__(name, value)


def iterate_indices(selection):
    space = {'lis':[]}
    cmd.iterate(selection, 'lis.append(index)', space=space)
    for l in space['lis']:
        yield l

def iterate_neighbours(selection):
    if isinstance(selection, int):
        selection = '(idx {})'.format(selection)
    for idx in iterate_indices('neighbor ({})'.format(selection)):
        yield idx

def get_atom_parameters(selection):
    namespace={'lis':[]}
    exposed_variables=['model','name','resn','resi','resv','chain','alt','elem','q','b','segi','type','formal_charge','partial_charge','numeric_type','text_type','stereo','ID','rank','index','vdw','ss','color','reps','protons','s'] # not supported: p
    command='lis.append({'+', '.join(['"{0}": {0}'.format(ev) for ev in exposed_variables])+'})'
    cmd.iterate('({})'.format(selection),command, space=namespace)
    return namespace['lis']



cmd.extend('get_atom_parameters',get_atom_parameters)
