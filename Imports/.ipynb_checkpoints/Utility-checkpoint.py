def find_conic_discriminant(conic):
    # do not consider l, m as variables for the conic
    P = conic.parent()
    v = [variable for variable in conic.variables() if variable in P.gens()[0:4]]
    conic_discriminant = matrix([
        [conic.coefficient(v[0] ^ 2),
         (1 / 2) * conic.coefficient(v[0]).coefficient(v[1]),
         (1 / 2) * conic.coefficient(v[0]).coefficient(v[2])],
        [(1 / 2) * conic.coefficient(v[1]).coefficient(v[0]),
         conic.coefficient(v[1] ^ 2),
         (1 / 2) * conic.coefficient(v[1]).coefficient(v[2])],
        [(1 / 2) * conic.coefficient(v[2]).coefficient(v[0]),
         (1 / 2) * conic.coefficient(v[2]).coefficient(v[1]),
         conic.coefficient(v[2] ^ 2)]
    ])
    return conic_discriminant.det().factor()

# L-sets ------------------------------------------------------------------------------------------


def find_all_L_sets(lines_dict):
    five_tuples = []
    keys = lines_dict.keys()
    for key1 in keys:
        for key2 in get_incident_keys(key1):
            for key3 in set(get_non_incident_keys(key1)).intersection(get_incident_keys(key2)):
                for key4 in set(get_incident_keys(key1)).intersection(
                    get_non_incident_keys(key2)).intersection(
                    get_incident_keys(key3)):
                    
                    for key5 in set(get_non_incident_keys(key1)).intersection(
                        get_incident_keys(key2)).intersection(
                        get_non_incident_keys(key3)).intersection(
                        get_non_incident_keys(key4)):
                        
                        five_tuples.append([key1, key2, key3, key4, key5])
    return five_tuples


def get_non_incident_keys(key):
    not_incident_lines = []
    if key[0] == 'E':
        index = int(key[1])
        not_incident_lines += ['E' + str(i + 1) for i in range(6) if i + 1 != index]
        not_incident_lines.append('G' + str(index))
        not_incident_lines += ['F' + str(i + 1) + str(j + 1) for i in range(6) for j in range(i + 1, 6)
                               if i + 1 != index and j + 1 != index]
    elif key[0] == 'G':
        index = int(key[1])
        not_incident_lines += ['G' + str(i + 1) for i in range(6) if i + 1 != index]
        not_incident_lines.append('E' + str(index))
        not_incident_lines += ['F' + str(i + 1) + str(j + 1) for i in range(6) for j in range(i + 1, 6)
                               if i + 1 != index and j + 1 != index]
    else:
        index1 = int(key[1])
        index2 = int(key[2])
        not_incident_lines += ['E' + str(i + 1) for i in range(6) if i + 1 != index1 and i + 1 != index2]
        not_incident_lines += ['G' + str(i + 1) for i in range(6) if i + 1 != index1 and i + 1 != index2]
        not_incident_lines += ['F' + str(i + 1) + str(j + 1) for i in range(6) for j in range(i + 1, 6)
                               if len(set([i + 1, j + 1, index1, index2])) < 4]
        not_incident_lines.remove('F' + str(index1) + str(index2))
    return not_incident_lines


def get_incident_keys(key):
    incident_lines = []
    if key[0] == 'E':
        index = int(key[1])
        incident_lines += ['G' + str(i + 1) for i in range(6) if i + 1 != index]
        incident_lines += ['F' + str(i + 1) + str(j + 1) for i in range(6) for j in range(i + 1, 6)
                           if i + 1 == index or j + 1 == index]
    elif key[0] == 'G':
        index = int(key[1])
        incident_lines += ['E' + str(i + 1) for i in range(6) if i + 1 != index]
        incident_lines += ['F' + str(i + 1) + str(j + 1) for i in range(6) for j in range(i + 1, 6)
                           if i + 1 == index or j + 1 == index]
    else:
        index1 = int(key[1])
        index2 = int(key[2])
        incident_lines += ['E' + str(i + 1) for i in range(6) if i + 1 == index1 or i + 1 == index2]
        incident_lines += ['G' + str(i + 1) for i in range(6) if i + 1 == index1 or i + 1 == index2]
        incident_lines += ['F' + str(i + 1) + str(j + 1) for i in range(6) for j in range(i + 1, 6)
                           if len(set([i + 1, j + 1, index1, index2])) == 4]
    return incident_lines


#from the classification as E, F, G returns the actual lines in plucker coordinates
def get_L_set_in_plucker(classified_lines, L_set):
    return tuple(map(lambda uu: classified_lines[uu], L_set))

        
# A = l1 ∩ l2, B = l1 ∩ l4, C = l3 ∩ l4, D = l2 ∩ l3, E = l2 ∩ l5
# P := l4 ∩ <l2 + l5>, Q := <P + D> ∩ l5.
def get_five_points_in_general_position(L_set):
    A = L_set[0].intersection_point(L_set[1])
    B = L_set[0].intersection_point(L_set[3])
    C = L_set[2].intersection_point(L_set[3])
    D = L_set[1].intersection_point(L_set[2])
    E = L_set[1].intersection_point(L_set[4])
    
    plane_l2_l5 = get_plane_containing_two_incident_lines(L_set[1], L_set[4])
    plane_coeff = plane_coefficients(plane_l2_l5)
    a = vector(plane_coeff[1:4])
    pl = L_set[3].plucker
    d = vector(pl[0:3])
    m = vector([pl[5], -pl[4], pl[3]])
    P = [a.dot_product(d)]+list(a.cross_product(m) - plane_coeff[0]*d) 
    
    plane_l3_l4 = get_plane_containing_two_incident_lines(L_set[2], L_set[3])
    if plane_l3_l4 is None:
        print(L_set)
    line_P_D = Line([plane_l2_l5, plane_l3_l4])
    Q = line_P_D.intersection_point(L_set[4])
    return A,B,C,E,Q
    
def get_plane_containing_two_incident_lines(line1, line2):
    P1, P2 = line1.points
    PL2 = line2.points
    vrs = line1.P.gens()[0:4]
    for point in PL2:
        #check if point is on line1
        if matrix([P1, P2, point]).rank() > 2:
            #take determinant of 4x4 matrix to get equation
            plane_factored = matrix([P1, P2, point, vrs]).det().factor()
            return [fct[0] for fct in plane_factored if [v in fct[0].variables() for v in vrs] != [False for i in range(4)]][0]
                        

def solve_linear_system2(eqns, variables, param):
    A = matrix([[eqn.coefficient(var) for var in variables] for eqn in eqns])
    b = matrix([sum([-eqn.coefficient(par)*par for par in param]) for eqn in eqns]).T
    return A.solve_right(b)


def find_projectivity(base_five_points, L_set2):
    P = L_set2[0].P
    S = PolynomialRing(P.base_ring(), 21, 'v')
    SS = PolynomialRing(P.base_ring(), P.gens()+S.gens())
    vrs = SS.gens()[-21:]
    points2 = get_five_points_in_general_position(L_set2)
    M = matrix([[var for var in vrs[i:i+4]] for i in range(0, 15, 4)])
    system = matrix([base_five_points[i] for i in range(5)])*M
    b = diagonal_matrix(vrs[-5:])*matrix(points2[i] for i in range(5))
    eqn = system - b
    sol = solve_linear_system2(eqn.list(), vrs[0:20], [vrs[-1]])
    sol = sol*sol.denominator()
    sol = [SS(el) for el in sol.list()]
    sost = {vrs[i] : sol[i] for i in range(16)}
    proj = M.subs(sost).subs({SS.gens()[-1] : 1})
    proj = (proj/gcd(proj.list()) ).list()
    proj = matrix([[P(el) for el in proj[i:i+4]] for i in range(0, 15, 4)])
    proj = proj/gcd(gcd(proj.coefficients()[i].coefficients()) for i in range(len(proj.coefficients())))
    return proj*proj.denominator()
    
def find_all_projectivities(classified_lines, L_set, base_five_points):
    L2 = get_L_set_in_plucker(classified_lines, L_set)
    M = find_projectivity(base_five_points, L2)
    return M  

def find_proj_parallel(args):
    return find_all_projectivities(*args)

def find_all_proj_parallel(cl_lines, all_L_sets):
    if os.name == "nt":
        freeze_support()
    pool = mp.Pool(mp.cpu_count()-1)
    L_set_base = get_L_set_in_plucker(cl_lines, ['E1', 'G4', 'E2', 'G3', 'E3'])
    base_five_points = get_five_points_in_general_position(L_set_base)
    all_param = ((cl_lines, L_set, base_five_points) for L_set in all_L_sets)
    result = pool.map(find_proj_parallel, all_param) 
    pool.close()
    return result


#Cubics-----------------------------------------------------------------------------------------------------


def find_simmetry(cubic, proj):
    P = cubic.parent()
    vrs = vector(P.gens()[0:4])
    change_coord = vector(vrs)*proj 
    sost = {vrs[i] : change_coord[i] for i in range(4)}
    new_cubic = cubic.subs(sost)
    if cubic.coefficient(vrs[0]^2*vrs[1])*new_cubic - new_cubic.coefficient(vrs[0]^2*vrs[1])*cubic == 0:
        return proj
    else:
        return None

def find_simmetries_wrapper(args):
    return find_simmetry(*args)

def find_simmetries_parallel(cubic, all_projectivities):
    if os.name == "nt":
        freeze_support()
    pool = mp.Pool(mp.cpu_count()-1)
    all_param = ((cubic, proj) for proj in all_projectivities)
    result = [el for el in pool.map(find_simmetries_wrapper, all_param) if el is not None]
    pool.close()
    return result
    
    
# Utility -----------------------------------------------------------------------------

        
def plane_coefficients(plane):
    return [plane.coefficient(vr) for vr in plane.parent().gens()[0:4]]


def are_vectors_proportional(elem1, elem2):
    e1 = elem1.nonzero_positions()
    e2 = elem2.nonzero_positions()
    if len(e1)!=len(e2):
        return False
    for i in range(len(e1)):
        if e1[i] != e2[i]:
            return False
    lin_comb = elem2[e2[0]]*elem1 - elem1[e1[0]]*elem2
    if lin_comb == vector([0 for i in range(len(elem1))]):
        return True
    else:
        return False
    
def are_matrices_equal(m1, m2):
    elem1 = vector(m1.list())
    elem2 = vector(m2.list())
    return are_vectors_proportional(elem1, elem2)
    
    
def solve_linear_system(eqns, variables, param):
    A = matrix([[eqn.coefficient(var) for var in variables] for eqn in eqns])
    b = matrix([sum([-eqn.coefficient(par)*par for par in param]) for eqn in eqns]).T
    sol = A.adjugate()*b
    return [sol[i,0] for i in range(len(variables))] + [det(A)*par for par in param]

            
def find_all_tritangent_planes(cl_lines):
    all_triplets = find_all_triplets_of_coplanar_lines()
    planes = []
    for triplet in all_triplets:
        line1 = cl_lines.get(triplet[0])
        line2 = cl_lines.get(triplet[1])        
        plane = get_plane_containing_two_incident_lines(line1, line2)
        lines_dict = {k:cl_lines.get(k) for k in triplet}
        planes.append(tritangent_plane(plane, lines_dict))
    return planes


def find_all_triplets_of_coplanar_lines():
    all_triplets = []
    for i in range(1,7):
        for j in range(1,7):
            if i < j:
                all_triplets.append(['E'+str(i), 'G'+str(j), 'F'+str(i)+str(j)])
            elif i>j:
                all_triplets.append(['E'+str(i), 'G'+str(j), 'F'+str(j)+str(i)])                
    for j in range(2,7):
        for k in range(2,6):
            for l in range(3,7):
                for m in range(2,6):
                    for n in range(3,7):
                        if k<l and m<n and k<n and set((j,k,l,m,n)) == set((2,3,4,5,6)):
                            if ['F'+str(1)+str(j), 'F'+str(m)+str(n), 'F'+str(k)+str(l)] not in all_triplets:
                                all_triplets.append(['F'+str(1)+str(j), 'F'+str(k)+str(l), 'F'+str(m)+str(n)])
    return all_triplets    
    
def remove_sing_factors(poly, sing):
    sing_fact = [el[0] for el in list(sing)]+[-el[0] for el in list(sing)]
    poly_f = poly.factor()
    poly_fact = [el[0] for el in list(poly_f)]
    poly_powers = [el[1] for el in list(poly_f)]
    product = 1
    for fact in sing_fact:
        if fact in poly_fact:
            product = product*fact^poly_powers[poly_fact.index(fact)]
    return poly//product
    
    
def change_coord(proj):
    vrs = proj.base_ring().gens()[0:4]
    change_coord = vector(vrs)*proj
    return {vrs[i] : change_coord[i] for i in range(4)}
    
    
def find_all_permutations(keys):
    all_perm = []
    for key1 in keys:
        non_in_1 = set(get_non_incident_keys(key1))
        in_1 = set(get_incident_keys(key1))
        for key2 in non_in_1:
            non_in_2 = set(get_non_incident_keys(key2))
            in_2 = set(get_incident_keys(key2))
            for key3 in non_in_1.intersection(non_in_2):
                non_in_3 = set(get_non_incident_keys(key3))
                in_3 = set(get_incident_keys(key3))
                for key4 in non_in_1.intersection(non_in_2).intersection(non_in_3):
                    non_in_4 = set(get_non_incident_keys(key4))
                    in_4 = set(get_incident_keys(key4))
                    for key5 in non_in_1.intersection(non_in_2).intersection(non_in_3).intersection(non_in_4):
                        non_in_5 = set(get_non_incident_keys(key5))
                        in_5 = set(get_incident_keys(key5))
                        for key6 in non_in_1.intersection(non_in_2).intersection(non_in_3).intersection(non_in_4).intersection(non_in_5):
                            E = [key1, key2, key3, key4, key5, key6]
                            non_in_6 = set(get_non_incident_keys(key6))
                            in_6 = set(get_incident_keys(key6))
                            not_in = [non_in_1, non_in_2, non_in_3, non_in_4, non_in_5, non_in_6]
                            inc = [in_1, in_2, in_3, in_4, in_5, in_6]
                            G = []
                            for i in range(6):
                                s = set([key for key in keys if key not in E])
                                for j in range(6):
                                    if j!=i:
                                        s.intersection_update(inc[j])
                                    else:
                                        s.intersection_update(not_in[j])
                                G.append(list(s)[0])
                            F=[]
                            for i in range(5):
                                for j in range(i+1, 6):
                                    s = set([key for key in keys if key not in E and key not in G])
                                    for k in range(6):
                                        if i == k or j == k:
                                            s.intersection_update(inc[k])
                                        else:
                                            s.intersection_update(not_in[k])
                                    F.append(list(s)[0])
                            all_perm.append(E+G+F)                           
    return all_perm
    
def find_conditions_for_subfamilies(cubic, projectivities, simmetries):
    mon = ((x+y+z+t)^3).monomials()
    conditions = []
    sing_cubics_factored = cubic.sing_cubic.factor()
    for M in [proj for proj in projectivities if proj not in simmetries]:
        sost = change_coord(M)
        new_cubic = remove_sing_factors(cubic.eqn.subs(sost), sing_cubics_factored)    
        minor = matrix([[new_cubic.coefficient(mn) for mn in mon], [cubic.eqn.coefficient(mn) for mn in mon]]).minors(2)
        minor = [remove_sing_factors(el, sing_cubics_factored) for el in minor if el !=0]
        prim_deco = cubic.P.ideal(minor).radical().primary_decomposition()
        for ideale in prim_deco:
            if is_ideal_valid(ideale, cubic, sing_cubics_factored):
                conditions.append(ideale.gens())                        
    return list(set(conditions))


def is_ideal_valid(ideal, cubic, sing_cubics):
    if sing_cubics in ideal:
        return False
    for poly in list(set([pl.condition for pl in cubic.tritangent_planes if pl.condition != 0])):
        if poly in ideal:
            return False
    return True    
    
def apply_proj_to_eck(proj, eck):
    new_indices = []
    for i in range(len(eck)):
        new_indices.append(eck.index(eck[i]*proj)+1)
    return new_indices
    
def get_dual_coordinates(pl):
    d = list(pl[0:3])
    m = [pl[5], -pl[4], pl[3]]
    return m+d
        
def get_planes(pl):
    #p01, p02, p03, p23, p31, p12
    dpl = get_dual_coordinates(pl)
    vrs = vector([x,y,z,t])
    eqns = [vector([0, -dpl[0], -dpl[1], -dpl[2]]), 
            vector([dpl[0], 0, dpl[5], -dpl[4]]), 
            vector([dpl[1], dpl[5], 0, -dpl[3]]), 
            vector([dpl[2], -dpl[4], dpl[3], 0])]
    planes = [eqn.dot_product(vrs) for eqn in eqns if eqn.dot_product(vrs)!=0]
    for plane in planes[1:]:
        coeff1 = vector(plane_coefficients(planes[0]))
        coeff2 = vector(plane_coefficients(plane))
        if not are_vectors_proportional(coeff1, coeff2):
            return [planes[0], plane]
    return None
    
    
def apply_proj_to_lines(proj, lines):
    new_indices = []
    for i in range(len(lines)):
        new_indices.append(lines.index(lines[i].apply_proj(proj))+1)
    return new_indices
    
def from_perm_to_labels(perm):
    keys = ['E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'F12', 'F13', 'F14', 'F15', 'F16', 'F23', 'F24', 'F25', 'F26', 'F34', 'F35', 'F36', 'F45', 'F46', 'F56']
    return [keys[perm.dict()[key]-1] for key in perm.dict()]