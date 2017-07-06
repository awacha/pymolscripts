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

def helicize_beta_peptide(helixtype, selection='all'):
    """
    DESCRIPTION

    Adjust the torsion angles of a beta-peptide to different helical conformations

    USAGE

    set_beta_helix helixtype [, selection]
    
    ARGUMENTS
    
    helixtype = the type of the helix (either short or IUPAC name)
    
    selection = the selection to operate on. Must be a single peptide chain with 
        unique residue IDs (default: all)
    
    NOTES"""
    space = {'lis':[]}
    cmd.iterate(selection, 'lis.append(resv)', space=space)
    residues = set(space['lis'])
    for r in residues:
        calpha = '({}) and (name CA) and (resi {})'.format(selection, r)
        cbeta = '({}) and (name CB) and (resi {})'.format(selection,r)
        c = '({}) and (name C) and (resi {})'.format(selection,r)
        n = '({}) and (name N) and (resi {})'.format(selection,r)
        prevc = '(neighbor ({})) and (name C)'.format(n)
        nextn = '(neighbor ({})) and (name N)'.format(c)
        for name, selection in [('CA',calpha),('CB',cbeta),('C',c),('N',n),('prevC',prevc), ('nextN',nextn)]:
            cnt = cmd.count_atoms(selection)
            if cnt!=1:
                print('Error in residue {}: number of {} atoms found is {}'.format(r, name, cnt))
                break
        else:
            set_beta_helix(prevc, n, cbeta, calpha, c, nextn, helixtype)
    cmd.orient(selection)

helicize_beta_peptide.__doc__+=helixdocs

cmd.extend('helicize_beta_peptide',helicize_beta_peptide)

    
