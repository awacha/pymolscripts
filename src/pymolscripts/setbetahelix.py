from __future__ import print_function
from pymol import cmd
# torsion angles for various helices taken from:
# Beke et al., Journal of Computational Chemistry 2006, vol. 27, no. 1 pp. 20-38, DOI 10.1002/jcc.20299

helixtypes=[
    ('Z6M','g+[g+]a',(126.7,62.6,152.7),'B3LYP/6-311++G(d,p)'),
    ('Z6P','g-[g-]a',(-126.7,-62.6,-152.7),'B3LYP/6-311++G(d,p)'),
    ('Z8M','g+[g+]g-',(47.5,53.5,-104.3),'RHF/3-21G'),
    ('Z8P','g-[g-]g+',(-47.5,-53.5,104.3),'RHF/3-21G'),
    ('H8M','g+[a]g+',(76.8,-120.6,52.7),'B3LYP/6-311++G(d,p)'),
    ('H8P','g-[a]g-',(-76.8,120.6,-52.7),'B3LYP/6-311++G(d,p)'),
    ('H10M','g-[g-]g-',(-77.5,-51.8,-75.1),'B3LYP/6-311++G(d,p)'),
    ('H10P','g+[g+]g+',(77.5,51.8,75.1),'B3LYP/6-311++G(d,p)'),
    ('H12M','g+[g-]g+',(92.3,-90.0,104.6),'B3LYP/6-311++G(d,p)'),
    ('H12P','g-[g+]g-',(-92.3,90.0,-104.6),'B3LYP/6-311++G(d,p)'),
    ('H14M','a[g-]a',(-140.3,66.5,-136.8),'B3LYP/6-311++G(d,p)'),
    ('H14P','a[g+]a',(140.3,-66.5,136.8),'B3LYP/6-311++G(d,p)'),
    ('SM','g+[a]a',(70.5,176.2,168.9),'RHF/3-21G'),
    ('SP','g-[a]a',(-70.5,-176.2,-168.9),'RHF/3-21G'),
]

def set_beta_helix(prevC, N, CB, CA, C, nextN, helixtype, selection = 'all'):
    """Set torsion angles of a beta peptide backbone

    Inputs:
       prevC: selection of the amide carbon in the previous amino acid
       N: selection of the amide Nitrogen in this amino acid
       CB: selection of the beta carbon
       CA: selection of the alpha carbon
       C: selection of the amide carbon
       nextN: selection of the amide nitrogen of the next amino acid
       helixtype: the type of the helix
       selection: the main selection in which we operate. Defaults to "*"
"""
    try:
        helixparam = [h for h in helixtypes if h[0]==helixtype or h[1]==helixtype][0]
    except IndexError:
        raise ValueError('Unknown helix type: {}'.format(helixtype))
    atoms = ['({}) and ({})'.format(selection, atomidx) for atomidx in [prevC, N, CB, CA, C, nextN]]
    for i, angle in enumerate(helixparam[2]):
        cmd.set_dihedral(*(atoms[i:i+4]+[angle]))
    cmd.delete('pk1')
    cmd.delete('pk2')
    cmd.delete('pkbond')
    cmd.delete('pkmol')
    
helixdocs="\n    Known helix types:\n"
helixdocs+='      Perczel short name | IUPAC name |  phi   | theta  |   psi  |    theory level     \n'
helixdocs+='      -------------------+------------+--------+--------+--------+---------------------\n'
for ht in helixtypes:
    helixdocs+='      {:<19s}| {:<10s} | {:=6.1f} | {:=6.1f} | {:=6.1f} | {:<11s}\n'.format(ht[0],ht[1],ht[2][0],ht[2][1],ht[2][2],ht[3])
set_beta_helix.__doc__+=helixdocs


cmd.extend('set_beta_helix',set_beta_helix)

def get_atom_parameters(selection):
    namespace={'lis':[]}
    exposed_variables=['model','name','resn','resi','resv','chain','alt','elem','q','b','segi','type','formal_charge','partial_charge','numeric_type','text_type','stereo','ID','rank','index','vdw','ss','color','reps','protons','s'] # not supported: p
    command='lis.append({'+', '.join(['"{0}": {0}'.format(ev) for ev in exposed_variables])+'})'
    cmd.iterate('({})'.format(selection),command, space=namespace)
    return namespace['lis']



cmd.extend('get_atom_parameters',get_atom_parameters)

def select_by_index(sel, name):
    namespace={'lis':[]}
    cmd.iterate(sel, 'lis.append(index)', space=namespace)
    cmd.select(name, ' or '.join(['(index {})'.format(x) for x in namespace['lis']]))
    

def select_beta_backbone(firstnitrogen, selectionname='backbone', selection_n = None, selection_ca = None, selection_cb = None, selection_c=None):
    """Find and select backbone atoms in a beta peptide

    Inputs:
        firstnitrogen: selection of the first nitrogen
        selectionname: name of the final selection
        selection_n: name of the selection containing the nitrogens
        selection_ca: name of the selection containing the alpha carbons
        selection_cb: name of the selection containing the beta carbons
        selection_c: name of the selection containing the amide carbons
"""
    # first find all nitrogens
    cmd.select('_nitrogens',firstnitrogen)
    lastlength=0
    while len(get_atom_parameters('_nitrogens'))>lastlength:
        lastlength=len(get_atom_parameters('_nitrogens'))
        cmd.select('_nitrogens', '(_nitrogens extend 4) and (symbol N)')
    # nitrogens are now in the _nitrogens selection
    # find amide carbons. They are neighbours of nitrogens and they also have an oxygen neighbour
    namespace={'lis':[]}
    cmd.iterate('(_nitrogens extend 1) and (symbol C)', 'lis.append(index)',space=namespace)
    # namespace['lis'] now has the indices of carbon neighbours of N
    cbs=[]
    cs=[]
    for idx in namespace['lis']:
        elements=[a['elem'] for a in get_atom_parameters('neighbor (index {})'.format(idx))]
        if elements.count('O')==1 and elements.count('N')==1 and elements.count('C')==1 and len(elements)==3:
            cs.append(idx)
        else:
            cbs.append(idx)
    cmd.select('_c',' or '.join(['(index {})'.format(x) for x in cs]))
    cmd.select('_cb',' or '.join(['(index {})'.format(x) for x in cbs]))
    cmd.select('_ca','(neighbor _cb) and (neighbor _c) and (symbol C)')

    select_by_index('(_c) or (_cb) or (_ca) or (_nitrogens)', selectionname)
    if selection_c is not None:
        select_by_index('(_c)', selection_c)
    if selection_ca is not None:
        select_by_index('(_ca)', selection_ca)
    if selection_cb is not None:
        select_by_index('(_cb)', selection_cb)
    if selection_n is not None:
        select_by_index('(_nitrogens)', selection_n)
    cmd.delete('(_c)')
    cmd.delete('(_cb)')
    cmd.delete('(_ca)')
    cmd.delete('(_nitrogens)')

cmd.extend('select_beta_backbone',select_beta_backbone)
    
def name_peptide_atoms(firstnitrogen):
    """Name the atoms in a beta peptide.

    Inputs:
       firstnitrogen: a selection of a single atom, the first nitrogen in the
          peptide backbone
    """
    select_beta_backbone(firstnitrogen,'backbone','_n','_ca','_cb','_c')
    cmd.alter('(_c)', 'name="C"')
    cmd.alter('(_ca)', 'name="CA"')
    cmd.alter('(_cb)', 'name="CB"')
    cmd.alter('(_n)', 'name="N"')
    cmd.alter('((neighbor _c) and (symbol O))', 'name="O"')
    cmd.alter('((neighbor _n) and (symbol H))', 'name="H"')
    cmd.alter('((neighbor _cb) and (symbol H))', 'name="HB"')
    cmd.select('alreadynamed','name C+CA+CB+N+O+H+HB')
    # name CA hydrogens
    space={'lis':[]}
    cmd.iterate('(_ca)','lis.append(index)', space=space)
    for calpha in space['lis']:
        spc={'lis':[]}
        cmd.iterate('(neighbor idx. {}) and (symbol H)'.format(calpha),'lis.append(index)', space=spc)
        for i,idx in enumerate(spc['lis']):
            cmd.alter('(idx. {})'.format(idx),'name="HA{}"'.format(i+1))
    cmd.select('alreadynamed','(alreadynamed) or (name HA*)')
    space['lis']=[]
    cmd.iterate('(bychain _n)','lis.append((resi, resn))', space=space)
    for resi,resn in sorted(set(space['lis'])):
        print(resi, resn)
        if resn == 'B3V':
            name_val_sidechain('resi {}'.format(resi))
        elif resn == 'B3A':
            name_ala_sidechain('resi {}'.format(resi))
        elif resn == 'B3L':
            name_leu_sidechain('resi {}'.format(resi))
        elif resn == 'ACE':
            name_ace_sidechain('resi {}'.format(resi))
        else:
            print('Unknown resn: {}. Resi: {}'.format(resn, resi))
        
    
def name_val_sidechain(selection):
    """Name the atoms in a valine sidechain

    Inputs:
       selection: a selection containing a SINGLE valine residue
    """
    cmd.select('notnamed', '({}) and not (name C+CA+CB+N+O+H+HA*+HB)'.format(selection))
    # name CG atom
    cmd.alter('(neighbor (name CB)) and (symbol C) and (notnamed)', 'name="CG"')
    cmd.select('notnamed','(notnamed) and not (name CG)'.format(selection))
    # name Hydrogen of CG
    cmd.alter('(neighbor (name CG)) and (symbol H) and (notnamed)', 'name="HG"')
    cmd.select('notnamed','(notnamed) and not(name HG)'.format(selection))
    # name CD1 and CD2
    space={'l':[],'l1':[]}
    cmd.iterate('(neighbor (name CG)) and (symbol C) and (notnamed)', 'l.append(index)',space=space)
    print(len(space['l']))
    for i,idx in enumerate(space['l']):
        print('CD{}'.format(i+1))
        cmd.alter('index {}'.format(idx), 'name="CD{}"'.format(i+1))
        space['l1']=[]
        cmd.iterate('(neighbor (index {})) and (symbol H) and (notnamed)'.format(idx), 'l1.append(index)', space=space)
        for j, idx1 in enumerate(space['l1']):
            cmd.alter('index {}'.format(idx1), 'name="HD{}{}"'.format(i+1,j+1))
    
def name_ala_sidechain(selection):
    """Name the atoms in a alanine sidechain

    Inputs:
       selection: a selection containing a SINGLE alanine residue
    """
    cmd.select('notnamed', '({}) and not (name C+CA+CB+N+O+H+HA*+HB)'.format(selection))
    # name CG atom
    cmd.alter('(neighbor (name CB)) and (symbol C) and (notnamed)', 'name="CG"')
    cmd.select('notnamed','(notnamed) and not (name CG)'.format(selection))
    # name Hydrogens of CG
    space={'l':[]}
    cmd.iterate('(neighbor (name CG)) and (symbol H) and (notnamed)', 'l.append(index)', space=space)
    for i,idx in enumerate(space['l']):
        cmd.alter('index {}'.format(idx), 'name="HG{}"'.format(i+1))

def name_leu_sidechain(selection):
    """Name the atoms in a leucine sidechain

    Inputs:
       selection: a selection containing a SINGLE leucine residue
    """
    cmd.select('notnamed', '({}) and not (name C+CA+CB+N+O+H+HA*+HB)'.format(selection))
    # name CG atom
    cmd.alter('(neighbor (name CB)) and (symbol C) and (notnamed)', 'name="CG"')
    cmd.select('notnamed','(notnamed) and not (name CG)'.format(selection))
    # name Hydrogens of CG
    space={'l':[]}
    cmd.iterate('(neighbor (name CG)) and (symbol H) and (notnamed)', 'l.append(index)', space=space)
    if len(space['l'])==1:
        cmd.alter('index {}'.format(space['l'][0]), 'name="HG"')
    else:
        for i,idx in enumerate(space['l']):
            cmd.alter('index {}'.format(idx), 'name="HG{}"'.format(i+1))
    cmd.select('notnamed','(notnamed) and not (name HG*)')
    # name CD atom
    cmd.alter('(neighbor (name CG)) and (symbol C) and (notnamed)', 'name="CD"')
    cmd.select('notnamed','(notnamed) and not (name CD)'.format(selection))
    # name Hydrogen of CD
    space={'l':[]}
    cmd.iterate('(neighbor (name CD)) and (symbol H) and (notnamed)', 'l.append(index)', space=space)
    if len(space['l'])==1:
        cmd.alter('index {}'.format(space['l'][0]), 'name="HD"')
    else:
        for i,idx in enumerate(space['l']):
            cmd.alter('index {}'.format(idx), 'name="HD{}"'.format(i+1))
    cmd.select('notnamed','(notnamed) and not (name HD*)')
    # name CE1 and CE2
    space={'l':[],'l1':[]}
    cmd.iterate('(neighbor (name CD)) and (symbol C) and (notnamed)', 'l.append(index)',space=space)
    for i,idx in enumerate(space['l']):
        cmd.alter('index {}'.format(idx), 'name="CE{}"'.format(i+1))
        space['l1']=[]
        cmd.iterate('(neighbor (index {})) and (symbol H) and (notnamed)'.format(idx), 'l1.append(index)', space=space)
        for j, idx1 in enumerate(space['l1']):
            cmd.alter('index {}'.format(idx1), 'name="HE{}{}"'.format(i+1,j+1))

def name_ace_sidechain(selection):
    cmd.alter('({}) and (symbol O)'.format(selection),'name="O"')
    cmd.alter('(neighbor (({}) and (symbol O))) and (symbol C)'.format(selection), 'name="C"')
    cmd.alter('(neighbor (neighbor (({}) and (symbol O))) and (symbol C)) and (symbol C)'.format(selection), 'name="CH3"')
    space={'l':[]}
    cmd.iterate('(neighbor (({}) and (name CH3))) and (symbol H)'.format(selection), 'l.append(index)', space=space)
    for i, idx in enumerate(space['l']):
        cmd.alter('index {}'.format(idx),'name="HH3{}"'.format(i+1))

cmd.extend('name_peptide_atoms',name_peptide_atoms)
    
def helicize_peptide(firstnitrogen, helixtype):
    select_beta_backbone(firstnitrogen, '_backbone', '_n', '_ca', '_cb','_c')
    space={'l':[]}
    cmd.iterate(firstnitrogen, 'l.append(index)',space=space)
    if len(space['l'])!=1:
        raise ValueError('The selection "{}" must contain only one atom, a nitrogen'.format(firstnitrogen))
    nitrogen = space['l'][0]
    while True:
        space['l']=[]
        cmd.iterate('index {}'.format(nitrogen),'l.append(resi)',space=space)
        try:
            space['l1']=[]
            cmd.iterate('(resi {}) and (name N)'.format(int(space['l'][0])+1),
                        'l1.append(index)', space=space)
            nextnitrogen = space['l1'][0]
        except IndexError:
            break
        set_beta_helix('(neighbor (index {})) and (name C)'.format(nitrogen),
                       '(index {})'.format(nitrogen),
                       '(neighbor (index {})) and (name CB)'.format(nitrogen),
                       '(byres (index {})) and (name CA)'.format(nitrogen),
                       '(byres (index {})) and (name C)'.format(nitrogen),
                       '(index {})'.format(nextnitrogen), helixtype
        )
        nitrogen=nextnitrogen 
    cmd.orient('(bychain ({}))'.format(firstnitrogen))
        
cmd.extend('helicize_peptide',helicize_peptide)

    
