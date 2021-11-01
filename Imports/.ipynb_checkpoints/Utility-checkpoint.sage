def find_conic_discriminant(conic):
    v = [variable for variable in conic.variables() if variable in conic.parent().gens()[0:4]]
    discriminant_matrix = matrix([
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
    return discriminant_matrix.det().factor()


def find_all_L_sets():
    all_L_sets = []
    keys = ['E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'F12', 'F13', 'F14', 'F15', 'F16', 'F23', 'F24', 'F25', 'F26', 'F34', 'F35', 'F36', 'F45', 'F46', 'F56']
    for key1 in keys:
        incident_1 = set(get_incident_keys(key1))
        non_incident_1 = set(get_non_incident_keys(key1))
        for key2 in incident_1:
            incident_2 = set(get_incident_keys(key2))
            non_incident_2 = set(get_non_incident_keys(key2))
            for key3 in non_incident_1.intersection(incident_2):
                incident_3 = set(get_incident_keys(key3))
                non_incident_3 = set(get_non_incident_keys(key3))
                for key4 in incident_1.intersection(non_incident_2).intersection(incident_3):   
                    non_incident_4 = set(get_non_incident_keys(key4))
                    for key5 in non_incident_1.intersection(incident_2).intersection(non_incident_3).intersection(non_incident_4):                       
                        all_L_sets.append((key1, key2, key3, key4, key5))
    return all_L_sets

# E_i does not meet G_j if i=j
# E_i, G_i do not meet F_jk if i!=j and i!=k
# F_ij does not meet F_kl if i,j,k,l are not all different
def get_non_incident_keys(key):
    non_incident_keys = []
    if key[0] == 'E':
        index = int(key[1])
        non_incident_keys += ['E' + str(i + 1) for i in range(6) if i + 1 != index]
        non_incident_keys.append('G' + str(index))
        non_incident_keys += ['F' + str(i + 1) + str(j + 1) for i in range(6) for j in range(i + 1, 6)
                              if i + 1 != index and j + 1 != index]
    elif key[0] == 'G':
        index = int(key[1])
        non_incident_keys += ['G' + str(i + 1) for i in range(6) if i + 1 != index]
        non_incident_keys.append('E' + str(index))
        non_incident_keys += ['F' + str(i + 1) + str(j + 1) for i in range(6) for j in range(i + 1, 6)
                              if i + 1 != index and j + 1 != index]
    else:
        index1 = int(key[1])
        index2 = int(key[2])
        non_incident_keys += ['E' + str(i + 1) for i in range(6) if i + 1 != index1 and i + 1 != index2]
        non_incident_keys += ['G' + str(i + 1) for i in range(6) if i + 1 != index1 and i + 1 != index2]
        non_incident_keys += ['F' + str(i + 1) + str(j + 1) for i in range(6) for j in range(i + 1, 6)
                              if len({i + 1, j + 1, index1, index2}) < 4]
        non_incident_keys.remove('F' + str(index1) + str(index2))
    return non_incident_keys


# E_i meets G_j if i!=j
# E_i, G_i meets F_jk if i=j or i=k
# F_ij meets F_kl if i,j,k,l are all different
def get_incident_keys(key):
    incident_keys = []
    if key[0] == 'E':
        index = int(key[1])
        incident_keys += ['G' + str(i + 1) for i in range(6) if i + 1 != index]
        incident_keys += ['F' + str(i + 1) + str(j + 1) for i in range(6) for j in range(i + 1, 6)
                          if i + 1 == index or j + 1 == index]
    elif key[0] == 'G':
        index = int(key[1])
        incident_keys += ['E' + str(i + 1) for i in range(6) if i + 1 != index]
        incident_keys += ['F' + str(i + 1) + str(j + 1) for i in range(6) for j in range(i + 1, 6)
                          if i + 1 == index or j + 1 == index]
    else:
        index1 = int(key[1])
        index2 = int(key[2])
        incident_keys += ['E' + str(i + 1) for i in range(6) if i + 1 == index1 or i + 1 == index2]
        incident_keys += ['G' + str(i + 1) for i in range(6) if i + 1 == index1 or i + 1 == index2]
        incident_keys += ['F' + str(i + 1) + str(j + 1) for i in range(6) for j in range(i + 1, 6)
                          if len({i + 1, j + 1, index1, index2}) == 4]
    return incident_keys


def find_all_permutations():
    all_perm = []
    keys = ['E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'F12', 'F13', 'F14', 'F15', 'F16',
            'F23', 'F24', 'F25', 'F26', 'F34', 'F35', 'F36', 'F45', 'F46', 'F56']
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
                        for key6 in non_in_1.intersection(non_in_2).intersection(non_in_3).intersection(
                                non_in_4).intersection(non_in_5):
                            non_in_6 = set(get_non_incident_keys(key6))
                            in_6 = set(get_incident_keys(key6))
                            E = [key1, key2, key3, key4, key5, key6]
                            non_incident = [non_in_1, non_in_2, non_in_3, non_in_4, non_in_5, non_in_6]
                            incident = [in_1, in_2, in_3, in_4, in_5, in_6]
                            G = []
                            for i in range(6):
                                s = set([key for key in keys if key not in E])
                                for j in range(6):
                                    if j != i:
                                        s.intersection_update(incident[j])
                                    else:
                                        s.intersection_update(non_incident[j])
                                G.append(list(s)[0])
                            F = []
                            for i in range(5):
                                for j in range(i + 1, 6):
                                    s = set([key for key in keys if key not in E and key not in G])
                                    for k in range(6):
                                        if i == k or j == k:
                                            s.intersection_update(incident[k])
                                        else:
                                            s.intersection_update(non_incident[k])
                                    F.append(list(s)[0])
                            all_perm.append(E + G + F)
    return all_perm
        
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

#TOFIX            
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
    line_P_D = Line([plane_l2_l5, plane_l3_l4])
    Q = line_P_D.intersection_point(L_set[4])
    return A,B,C,E,Q
    
#TOFIX            
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
                        

def solve_linear_system_in_fraction_field(eqns, variables, param):
    A = matrix([[eqn.coefficient(var) for var in variables] for eqn in eqns])
    b = matrix([sum([-eqn.coefficient(par)*par for par in param]) for eqn in eqns]).T
    return A.solve_right(b)


#TOFIX            
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
    sol = solve_linear_system_in_fraction_field(eqn.list(), vrs[0:20], [vrs[-1]])
    sol = sol*sol.denominator()
    sol = [SS(el) for el in sol.list()]
    sost = {vrs[i] : sol[i] for i in range(16)}
    proj = M.subs(sost).subs({SS.gens()[-1] : 1})
    proj = (proj/gcd(proj.list()) ).list()
    proj = matrix([[P(el) for el in proj[i:i+4]] for i in range(0, 15, 4)])
    proj = proj/gcd(gcd(proj.coefficients()[i].coefficients()) for i in range(len(proj.coefficients())))
    return proj*proj.denominator()
    
# Utility -----------------------------------------------------------------------------

        
def plane_coefficients(plane):
    return [plane.coefficient(vr) for vr in plane.parent().gens()[0:4]]


def are_vectors_proportional(vec1, vec2):
    nz1 = vec1.nonzero_positions()
    nz2 = vec2.nonzero_positions()
    #check if the vectors have same number of     
    #nonzero elements in the same positions
    if len(nz1)!=len(nz2):
        return False
    for i in range(len(nz1)):
        if nz1[i] != nz2[i]:
            return False
    #check if they are proportional
    lin_comb = vec2[nz2[0]]*vec1 - vec1[nz1[0]]*vec2
    if lin_comb == vector([0 for i in range(len(vec1))]):
        return True
    else:
        return False
    
def are_matrices_equal(m1, m2):
    vec1 = vector(m1.list())
    vec2 = vector(m2.list())
    return are_vectors_proportional(vec1, vec2)
    
#solve a linear system by computing the adjugate matrix   
def solve_linear_system(eqns, variables, param):
    A = matrix([[eqn.coefficient(var) for var in variables] for eqn in eqns])
    b = matrix([sum([-eqn.coefficient(par)*par for par in param]) for eqn in eqns]).T
    sol = A.adjugate()*b
    return [sol[i,0] for i in range(len(variables))] + [det(A)*par for par in param]

#TOFIX            
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
 
    
def remove_sing_factors(poly, sing_locus):
    sing_locus_factors = [el[0] for el in list(sing_locus)]+[-el[0] for el in list(sing_locus)]
    poly_factorization = poly.factor()
    poly_factors = [el[0] for el in list(poly_factorization)]
    poly_powers = [el[1] for el in list(poly_factorization)]
    product = 1
    for fact in sing_locus_factors:
        if fact in poly_factors:
            product = product*fact^poly_powers[poly_factors.index(fact)]
    return poly//product
    
    
def change_coord(proj):
    vrs = proj.base_ring().gens()[0:4]
    change_coord = vector(vrs)*proj
    return {vrs[i] : change_coord[i] for i in range(4)}

    
def find_conditions_for_subfamilies(cubic, projectivities, simmetries):
    mon = ((x+y+z+t)^3).monomials()
    conditions = []
    for M in [proj for proj in projectivities if proj not in simmetries]:
        sost = change_coord(M)
        new_cubic = remove_sing_factors(cubic.eqn.subs(sost), cubic.sing_cubic)    
        minor = list(set(matrix([[new_cubic.coefficient(mn) for mn in mon], [cubic.eqn.coefficient(mn) for mn in mon]]).minors(2)))
        minor = [remove_sing_factors(el, cubic.sing_cubic) for el in minor if el !=0]
        prim_deco = cubic.P.ideal(minor).radical().primary_decomposition()
        for ideale in prim_deco:
            if is_ideal_valid(ideale, cubic, cubic.sing_cubic):
                conditions.append(ideale.gens())                        
    return list(set(conditions))


def is_ideal_valid(ideal, cubic, sing_cubics):
    if sing_cubics in ideal:
        return False
    for poly in list(set([pl.conditions for pl in cubic.tritangent_planes if pl.conditions !=0])):
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