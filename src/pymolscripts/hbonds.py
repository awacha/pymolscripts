from pymol import cmd

def iterate_hbonds(selection, name='__hbonds'):
    # idea stolen from http://pymolwiki.org/index.php/get_raw_distances, Takanori Nakane and Thomas Holder
    cmd.dist(name, selection, selection, mode=2,label=0)
    raw_objects = cmd.get_session(name, 1,1,0,0)['names']
    space={'xyz2idx':{}}
    state=cmd.get_state()
    cmd.iterate_state(state, selection, 'xyz2idx[x,y,z] = (model, index)', space=space)
    points = raw_objects[0][5][2][state-1][1]
    for i in range(0,len(points),6):
        x1,y1,z1=tuple(points[i:i+3])
        x2,y2,z2=tuple(points[i+3:i+6])
        yield space['xyz2idx'][x1,y1,z1], space['xyz2idx'][x2,y2,z2], ((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)**0.5
    cmd.delete(name)
    return

def generate_hbond_constraints(selection, filename, name='__hbonds'):
    with open(filename, 'wt') as f:
        f.write('[ constraints ]\n')
        for idx1, idx2, dist in iterate_hbonds(selection, name):
            f.write('{:10d}{:10d}{:10.4f}\n'.format(idx1[1],idx2[1],dist))

cmd.extend('generate_hbond_constraints', generate_hbond_constraints)
