##
## list of the simbolic name of the lines.
keys = ['E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'G1', 'G2', 'G3', \
        'G4', 'G5', 'G6', 'F12', 'F13', 'F14', 'F15', 'F16', 'F23', \
	'F24', 'F25', 'F26', 'F34', 'F35', 'F36', 'F45', 'F46', 'F56']

## given a permutation prm (of the numbers 1, 2,..., 27) and a line rt
## expressed by a symbolic name, like "E1", "F23", etc.)
## returns a line (again of the form "E1", "G1", ...) which is the image
## under the permutation of the line rt.
## example:
## prm = Permutation('(1, 2)(4, 5, 7)')
## image_simbolic_line(prm, 'E1')
## sage: 'E2'

def image_simbolic_line(prm, rt):
    ind = keys.index(rt)
    Fdi_ind = prm(ind+1)
    return keys[Fdi_ind-1]

## returns a list of tuples of labels of lines. Each tuple represents a cycle,
## i.e. if a permutation is given by the cicles (a, b, c)(d, e)(f, g, h),
## where a, b, c, ... are integer in the interval [1, 27],
## returns a list of the form [(R1, R2, R3), (R4, R5), (R6, R7, R8)]
## where R1 is the simbolic name, precisely is keys[a-1], 
## R2 is keys[b-1] etc.
## example:
## prm = Permutation('(1, 2)(4, 5, 7)')
## perm_to_labels_lines(prm)
## sage: [('E1', 'E2'), ('E3',), ('E4', 'E5', 'G1'), ('E6',)]

def perm_to_labels_lines(prm):
    cics = prm.cycle_tuples()
    risp = []
    for cl in cics:
        risp.append(tuple(keys[i-1] for i in cl))
    return(risp)

## given a number i in range(1, 28), returns the line of position
## i-1 in keys:
## example:
## simbolic_line_from_index(5)
## sage: 'E5'

def simbolic_line_from_index(i):
    return(keys[i-1])
    

## find the lines (given with simbolic name) that remain fixed under
## the permutation prm
## example:
## prm = Permutation('(1,19,7,22)(2,3,9,8)(4,13,10,14)(5,20,11,23)(6,21,12,24)(15,27)(16,17)(25,26)')
## fixed_lines(prm)
## sage: ['F23']
## warning: is required that prm.orbit() is available.

def fixed_lines(prm):
    risp = []
    for i in range(27):
        if len(prm.orbit(i+1)) == 1:
            risp.append(i)
    return([keys[j] for j in risp])

## define the 45 tritangent planes 
## Their position is in accordance with table 1 of the paper.

tritangentPlanes = \
[['E1', 'G2', 'F12'], ['E1', 'G3', 'F13'], ['E1', 'G4', 'F14'],\
 ['E1', 'G5', 'F15'], ['E1', 'G6', 'F16'], ['E2', 'G1', 'F12'],\
 ['E2', 'G3', 'F23'],  ['E2', 'G4', 'F24'], ['E2', 'G5', 'F25'],\
 ['E2', 'G6', 'F26'], ['E3', 'G1', 'F13'], ['E3', 'G2', 'F23'],\
 ['E3', 'G4', 'F34'], ['E3', 'G5', 'F35'], ['E3', 'G6', 'F36'],\
 ['E4', 'G1', 'F14'], ['E4', 'G2', 'F24'], ['E4', 'G3', 'F34'],\
 ['E4', 'G5', 'F45'], ['E4', 'G6', 'F46'], ['E5', 'G1', 'F15'],\
 ['E5', 'G2', 'F25'], ['E5', 'G3', 'F35'], ['E5', 'G4', 'F45'],\
 ['E5', 'G6', 'F56'], ['E6', 'G1', 'F16'], ['E6', 'G2', 'F26'],\
 ['E6', 'G3', 'F36'], ['E6', 'G4', 'F46'], ['E6', 'G5', 'F56'],\
 ['F12', 'F34', 'F56'], ['F12', 'F35', 'F46'], ['F12', 'F36', 'F45'],\
 ['F13', 'F24', 'F56'], ['F13', 'F25', 'F46'], ['F13', 'F26', 'F45'],\
 ['F14', 'F23', 'F56'], ['F14', 'F25', 'F36'], ['F14', 'F26', 'F35'],\
 ['F15', 'F23', 'F46'], ['F15', 'F24', 'F36'], ['F15', 'F26', 'F34'],\
 ['F16', 'F23', 'F45'], ['F16', 'F24', 'F35'], ['F16', 'F25', 'F34']]

## defines 45 names: 't1', 't2', ..., 't45'
## label_tritantent_plane = ["t"+str(i) for i in range(1, 46)]

### given a tritangent plane (as a triplet of names of lines), returns
### its label "t_n" according to the table 1 of the paper.
## example:
## name_tritangent_plane(['E1', 'G3', 'F13'])
## sage: 't2'

def name_tritangent_plane(pl):
    n = tritangentPlanes.index(pl)
    return('t'+str(n+1))


## given a name of a tritangent plane, like "t3" or "t41", returns
## the triplet of names of lines corresponding to that plane
## in accordance with the table 1 of the paper.
## example:
## tritangent_plane_from_name('t7')
## sage: ['E2', 'G3', 'F23']

def tritangent_plane_from_name(nm):
    n = int(nm[1:]) 
    return(tritangentPlanes[n-1])
 
### given an index in the set {1, 2, ..., 45}, which corresponds to one
### of the 45 tritangent planes, returns the simbolic name of that plane
### according to table 1 of the paper.
## example:
## simbolic_plane_from_index(3)
## sage: 
def simbolic_plane_from_index(n):
    return('t'+str(n))



## given a permutation (of the numbers 1, 2, ..., 27) and a tritangent
## plane expressed as a triplet
## of simbolic lines, returns the image of the plane under the permutation
## example:
## prm = Permutation('(1,19,7,22)(2,3,9,8)(4,13,10,14)(5,20,11,23)(6,21,12,24)(15,27)(16,17)(25,26)')
## image_tritangent_plane(prm, ['E1', 'G4', 'F14'])
## sage: ['F13', 'F24', 'F56']

def image_tritangent_plane(prm, plane):
    conversion1 = [prm(keys.index(rt)+1) for rt in plane]
    conversion2 = list(map(lambda uu: keys[uu-1], conversion1))
    for pl in tritangentPlanes:
        if Set(conversion2) == Set(pl):
            return(pl)
    return("Impossible!")


## find the planes (given as triplets of simbolic lanes) that remain fixed under the permutation prm of the 45 planes.
## example:
## prm = Permutation('(1,16,8)(2,37,13)(4,38,24)(5,39,29)(7,31,11)(9,33,21)(10,32,26)(12,18,34)(14,28,25)(15,23,30)(19,41,22)(20,44,27)(35,43,42)(36,40,45)')
## fixed_planes(prm)
## sage: [['E1', 'G4', 'F14'], ['E2', 'G1', 'F12'], ['E4', 'G2', 'F24']]

def fixed_planes(prm):
    risp = []
    for i in range(45):
        if len(prm.orbit(i+1)) == 1:
            risp.append(i)
    return([tritangentPlanes[j] for j in risp])



## Define the group of permutation of the numbers 1, ..., 45
## which is the group of permutation of the 45 planes for the family SE
def perm_group_planes(SE, simm_SE):
    grpPlanes = []
    for simm in simm_SE:
        simm = simm[0]
        prm = Permutation(SE.apply_proj_to_lines(simm)).to_permutation_group_element()
        prmExt = Permutation([tritangentPlanes.index(image_tritangent_plane(prm, pln))+1 for pln in tritangentPlanes])
        grpPlanes.append(prmExt.to_permutation_group_element())
    return(grpPlanes)

## Define the group of permutation of the numbers 1, ..., 27
## which is the group of permutation of the 27 lines for the family SE
def perm_group_lines(SE, simm_SE):
    grpLines = []
    for simm in simm_SE:
        simm = simm[0]
        prm = Permutation(SE.apply_proj_to_lines(simm)).to_permutation_group_element()
        grpLines.append(prm)
    return(grpLines)

## returns a list of tuples of planes (expressed as triplets of labels of
## lines). Each tuple represent a cycle,
## i.e. if a permutation is given by the cicles (a, b, c)(d, e)(f, g, h)
## return a list of the form [(R1, R2, R3), (R4, R5), (R6, R7, R8)]
## where R1 is the triplet [x, y, z] (x, y, z are simbolic names of
## lines), precisely is tritangentPlanes[a-1], 
## R2 is tritangentPlanes[b-1] etc.
def perm_to_labels_planes(prm):
    cics = prm.cycle_tuples()
    risp = []
    for cl in cics:
        risp.append(tuple(tritangentPlanes[i-1] for i in cl))
    return(risp)
