from pymol import cmd
import numpy as np

def select_intralayer_waters(sele_headgroup, sele_water, selectionname):
    """
    DESCRIPTION

        "select_intralayer_waters" selects the water residues inside a bilayer

    USAGE
   
        select_intralayer_waters sele_headgroup, sele_water, selectionname

    ARGUMENTS
   
        sele_headgroup = string: a selection containing the head group atoms 
            (e.g. "resn P")
        
        sele_water = string: a selection of the water residues (e.g. "resn SOL")

        selectionname = string: the name of the selected water molecules

    NOTES
 
        The bilayer must lie in the xy plane. The z coordinates of the headgroup
        atoms are collected and the mean z of the upper and lower leaflet is 
        found. All waters where either atom is between these two limits are
        selected.
    """
    space={'lis':[]}
    cmd.iterate_state(-1, sele_headgroup, 'lis.append(z)',space=space)
    zmean=np.mean(space['lis'])
    zmin = np.mean([l for l in space['lis'] if l<zmean])
    zmax = np.mean([l for l in space['lis'] if l>zmean])
    print('zmin: {}'.format(zmin))
    print('zmax: {}'.format(zmax))
    sel='byres (({}) and (z>{}) and (z<{}))'.format(sele_water,zmin,zmax)
    print('selection:',sel)
    cmd.select(selectionname, sel)
    print('{} atoms selected.'.format(cmd.count_atoms(sel)))

cmd.extend('select_intralayer_waters',select_intralayer_waters)
