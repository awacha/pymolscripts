from pymol import cmd


def save_pdb_ordered(filename, selection='all', state=-1):
    """
    DESCRIPTION

        Save a PDB file while preserving the atom order.

    USAGE

        save_pdb_ordered filename [, selection [, state]]

    ARGUMENTS

        filename: the file name to save to

        selection: selected atoms to save to. Defaults to 'all'

        state: the state from where the coordinates are taken.
            Special values are:
                 0 : saves all states in a multi-model PDB file
                -1 : saves the current state (this is the default)
    """
    state = int(state)
    with open(filename, 'wt') as f:
        if state == -1:
            states = [cmd.get_state()]
        elif state == 0:
            states = list(range(cmd.count_states(selection)))
        else:
            states = [state]
        for i, s in enumerate(states):
            if len(states) > 1:
                f.write('MODEL{:>9d}\n'.format(i + 1))
            for a in sorted(cmd.get_model(selection, s).atom, key=lambda a: a.id):
                f.write('ATOM  {:5d} {:>4s} {:3s} {}{:4d}    {:8.3f}{:8.3f}{:8.3f}  1.00  0.00           {:s}\n'.format(
                    a.index, a.name, a.resn, a.chain, int(a.resi), a.coord[0], a.coord[1], a.coord[2], a.symbol
                ))
            if len(states) > 1:
                f.write('ENDMDL\n')
