from pymol import cmd
import math

def save_gro(filename, selection='(all)'):
    """
    DESCRIPTION
    
        "save_gro" saves content to a .gro file in a similar fashion as "save"

    USAGE

        save_gro filename [, selection]

    ARGUMENTS
        
        filename = string: file path to be written

        selection = string: atoms to save {default: (all)}

    NOTES

        In contrast to "save", only the current state can be saved.

    SEE ALSO

        load, save

    """
    with open(filename, 'wt') as f:
        f.write('{}\n'.format(selection))
        f.write('{:>5d}\n'.format(cmd.count_atoms(selection)))
        space={'l':[]}
        cmd.iterate_state(-1,selection,'l.append((resi,resn,name,index,x,y,z))',space=space)
        for resi, resn, name, index, x,y,z in space['l']:
            f.write('{:>5d}{:<5s}{:>5s}{:>5d}{:>8.3f}{:>8.3f}{:>8.3f}\n'.format(int(resi),resn,name,int(index),float(x)*0.1,float(y)*0.1,float(z)*0.1))
        try:
            a,b,c,alpha,beta,gamma,spacegroup=cmd.get_symmetry()
            alpha*=math.pi/180.
            beta*=math.pi/180.
            gamma*=math.pi/180.
            a/=10
            b/=10
            c/=10
            v1x = a
            v2x = b*math.cos(gamma)
            v2y = b*math.sin(gamma)
            v3x = c*math.cos(beta)
            v3y = c*(math.cos(alpha)-math.cos(gamma)*math.cos(beta))/(math.sin(gamma))
            v3z = (c**2-v3x**2-v3y**2)**0.5
            f.write(' {:>9.5f} {:>9.5f} {:>9.5f} {:>9.5f} {:>9.5f} {:>9.5f} {:>9.5f} {:>9.5f} {:>9.5f}\n'.format(v1x,v2y,v3z,0,0,v2x,0,v3x,v3y))
        except TypeError:
            f.write(' 0 0 0\n')
            

cmd.extend('save_gro',save_gro)

