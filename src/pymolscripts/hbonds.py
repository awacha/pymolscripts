from pymol import cmd
from .utils import iterate_indices

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
            f.write('{:10d}{:10d} 2 {:10.4f}\n'.format(idx1[1],idx2[1],dist/10.))

def generate_hbond_restraints(selection, filename, strength=1000, name='__hbonds'):
    strength=float(strength)
    with open(filename, 'wt') as f:
        f.write('[ bonds ]\n')
        for idx1, idx2, dist in iterate_hbonds(selection, name):
            f.write('{:10d}{:10d} 6 {:10.4f} {:.6f}\n'.format(idx1[1],idx2[1],dist/10.,strength))


def find_hbonds(selection_hydrogen, selection_acceptor, dmin=1, dmax=3):
    dmin = float(dmin)
    dmax=float(dmax)
    for oxygen in iterate_indices(selection_acceptor):
        hydrogendists = [(h, cmd.get_distance('idx {}'.format(oxygen),
                                              'idx {}'.format(h))) for h in iterate_indices(selection_hydrogen)]
        hydrogendists = [(h,d) for h,d in hydrogendists if d>=dmin and d<=dmax]
        if not hydrogendists:
            continue
        hydrogen, dist = sorted(hydrogendists, key=lambda x:x[1])[0]
        # this hydrogen is the nearest one to this oxygen. Try it the other way round
        oxygendists = [(o, cmd.get_distance('idx {}'.format(hydrogen),
                                            'idx {}'.format(o))) for o in iterate_indices(selection_acceptor)]
        oxygendists = [(o, d) for o,d in oxygendists if d >=dmin and d<=dmax]
        assert oxygendists
        o1, dist1 = sorted(oxygendists, key = lambda x:x[1])[0]
        if o1==oxygen:
            yield (oxygen, hydrogen, dist)

def beta_hbonds(selection_hydrogen, selection_acceptor, dmin=1, dmax=3):
    i=0
    for o, h, d in find_hbonds(selection_hydrogen, selection_acceptor, dmin, dmax):
        cmd.distance('dist{:04d}'.format(i), 'idx {}'.format(o), 'idx {}'.format(h), mode=0)
        cmd.group('hbonds', 'dist{:04d}'.format(i))
        i+=1

cmd.extend('beta_hbonds', beta_hbonds)

cmd.extend('generate_hbond_constraints', generate_hbond_constraints)

cmd.extend('generate_hbond_restraints', generate_hbond_restraints)