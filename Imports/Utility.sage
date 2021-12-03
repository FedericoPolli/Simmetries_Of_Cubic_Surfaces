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

# An L_set is a quintuple (l1, l2, l3, l4, l5) where 
# l1, l3, l5 are skew, l2 meets l1, l3, l5 while l4 meets l1 and l3.
def find_all_L_sets():
    all_L_sets = []
    keys = ['E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'F12', 'F13', 'F14', 'F15', 'F16',
            'F23', 'F24', 'F25', 'F26', 'F34', 'F35', 'F36', 'F45', 'F46', 'F56']
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
                    for key5 in non_incident_1.intersection(incident_2).intersection(non_incident_3).intersection(
                            non_incident_4):
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

# six mutually skew lines are enough to determine a permutation. It finds all the possible
# sets of 6 mutually skew lines E_1, ..., E_6 and then computes the correspondnig G_i and F_ij
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
                            non_incident = [non_in_1, non_in_2, non_in_3, non_in_4, non_in_5, non_in_6]
                            incident = [in_1, in_2, in_3, in_4, in_5, in_6]
                            E = [key1, key2, key3, key4, key5, key6]
                            
                            #G_i is the unique line meeting all E_j except E_i
                            G = []
                            for i in range(6):
                                s = set([key for key in keys if key not in E])
                                for j in range(6):
                                    if j != i:
                                        s.intersection_update(incident[j])
                                    else:
                                        s.intersection_update(non_incident[j])
                                G.append(list(s)[0])
                                
                            #F_ij is the unique line meeting E_i, E_j but not the other E_k
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
    for i in range(1, 7):
        for j in range(1, 7):
            if i < j:
                all_triplets.append(['E' + str(i), 'G' + str(j), 'F' + str(i) + str(j)])
            elif i > j:
                all_triplets.append(['E' + str(i), 'G' + str(j), 'F' + str(j) + str(i)])
    for j in range(2, 7):
        for k in range(2, 6):
            for l in range(3, 7):
                for m in range(2, 6):
                    for n in range(3, 7):
                        if k < l and m < n and k < n and {j, k, l, m, n} == {2, 3, 4, 5, 6}:
                            if ['F' + str(1) + str(j), 'F' + str(m) + str(n),
                                'F' + str(k) + str(l)] not in all_triplets:
                                all_triplets.append(
                                    ['F' + str(1) + str(j), 'F' + str(k) + str(l), 'F' + str(m) + str(n)])
    return all_triplets


# A = l1 ∩ l2, B = l1 ∩ l4, C = l3 ∩ l4, D = l2 ∩ l3, E = l2 ∩ l5
# P := l4 ∩ <l2 + l5>, Q := <P + D> ∩ l5.
def get_five_points_in_general_position(L_set):
    A = L_set[0].intersection_point(L_set[1])
    B = L_set[0].intersection_point(L_set[3])
    C = L_set[2].intersection_point(L_set[3])
    E = L_set[1].intersection_point(L_set[4])
    plane_l2_l5 = L_set[1].get_plane_containing_another_incident_line(L_set[4])
    plane_l3_l4 = L_set[2].get_plane_containing_another_incident_line(L_set[3])
    line_P_D = Line([plane_l2_l5, plane_l3_l4])
    Q = line_P_D.intersection_point(L_set[4])
    return A, B, C, E, Q


def get_two_planes_containing_line(point1, point2):
    P = point1.parent()
    vrs = P.gens()[0:4]
    possible_planes = [el for el in matrix([point1, point2, vrs]).minors(3) if el !=0]
    plane1 = possible_planes[0]
    coeff1 = vector(plane_coefficients(plane1))
    for plane in possible_planes[1:]:
        coeff2 = vector(plane_coefficients(plane))
        if not are_vectors_proportional(coeff1, coeff2):
            return [plane1, plane]
    return None


def solve_linear_system_in_fraction_field(eqns, variables, param):
    A = matrix([[eqn.coefficient(var) for var in variables] for eqn in eqns])
    b = matrix([sum([-eqn.coefficient(par) * par for par in param]) for eqn in eqns]).T
    return A.solve_right(b)


# finds the matrix sending the five points of the first L-set to the five points obtained from the other L-set
# maybe should be moved to cubic.sage
def find_projectivity(L_set1, L_set2):
    P = L_set2[0].P
    S = PolynomialRing(P.base_ring(), 21, 'v')  #build new polynomial ring with the variables needed
    SS = PolynomialRing(P.base_ring(), P.gens() + S.gens())
    vrs = SS.gens()[-21:]
    points1 = get_five_points_in_general_position(L_set1)
    points2 = get_five_points_in_general_position(L_set2)
    M = matrix([[var for var in vrs[i:i + 4]] for i in range(0, 15, 4)])
    
    # for each point P_i in the base five points, we want P_i*M = P_i'
    system = matrix([points1[i] for i in range(5)]) * M
    b = diagonal_matrix(vrs[-5:]) * matrix(points2[i] for i in range(5))
    eqn = system - b
    sol = solve_linear_system_in_fraction_field(eqn.list(), vrs[0:20], [vrs[-1]])
    
    #recast solution in polynomial ring
    sol = sol * sol.denominator()
    sol = [SS(el) for el in sol.list()]
    sost = {vrs[i]: sol[i] for i in range(16)}
    
    # get actual matrix without the added variables 
    proj = M.subs(sost).subs({SS.gens()[-1]: 1})
    
    #simplify matrix
    proj = (proj / gcd(proj.list())).list()
    proj = matrix([[P(el) for el in proj[i:i + 4]] for i in range(0, 15, 4)])
    proj = proj / gcd(gcd(proj.coefficients()[i].coefficients()) for i in range(len(proj.coefficients())))
    return proj * proj.denominator()


def plane_coefficients(plane):
    return [plane.coefficient(vr) for vr in plane.parent().gens()[0:4]]


def are_vectors_proportional(vec1, vec2):
    nz1 = vec1.nonzero_positions()
    nz2 = vec2.nonzero_positions()
    # check if the vectors have same number of nonzero elements
    if len(nz1) != len(nz2):
        return False
    # check if the nonzero elements are in the same positions
    for i in range(len(nz1)):
        if nz1[i] != nz2[i]:
            return False
    # check if vectors are proportional
    lin_comb = vec2[nz2[0]] * vec1 - vec1[nz1[0]] * vec2
    if lin_comb == vector([0 for _ in range(len(vec1))]):
        return True
    else:
        return False


def are_matrices_equal(m1, m2):
    vec1 = vector(m1.list())
    vec2 = vector(m2.list())
    return are_vectors_proportional(vec1, vec2)


# solve a linear system by computing the adjugate matrix
def solve_linear_system(eqns, variables, param):
    A = matrix([[eqn.coefficient(var) for var in variables] for eqn in eqns])
    b = matrix([sum([-eqn.coefficient(par) * par for par in param]) for eqn in eqns]).T
    sol = A.adjugate() * b
    return [sol[i, 0] for i in range(len(variables))] + [det(A) * par for par in param]

# removes singular factors from polynomial by checking each factor 
# of the singular locus against the polynomial
# sing_locus is a Factorization object
def remove_sing_factors(poly, sing_locus):
    sing_locus_factors = [el[0] for el in list(sing_locus)] + [-el[0] for el in list(sing_locus)]
    poly_factorization = poly.factor()
    poly_factors = [el[0] for el in list(poly_factorization)]
    poly_powers = [el[1] for el in list(poly_factorization)]
    product = 1
    for factor in sing_locus_factors:
        if factor in poly_factors:
            product = product * factor ^ poly_powers[poly_factors.index(factor)]
    return poly // product

#return a dictionary with the coordinate change associated to the projectivity
def change_coord(proj):
    vrs = proj.base_ring().gens()[0:4]
    coordinate_change = vector(vrs) * proj
    return {vrs[i]: coordinate_change[i] for i in range(4)}

#given a permutation of the 27 lines, returns the associated permuted labels
def from_perm_to_labels(perm):
    labels = ['E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'F12', 'F13', 'F14', 'F15', 'F16',
            'F23', 'F24', 'F25', 'F26', 'F34', 'F35', 'F36', 'F45', 'F46', 'F56']
    return [labels[perm.dict()[key] - 1] for key in perm.dict()]
  
#given the permuted labels, returns a permutation group element
def from_labels_to_perm(labels):
    keys = ['E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'F12', 'F13', 'F14', 'F15', 'F16', 'F23', 'F24', 'F25', 'F26', 'F34', 'F35', 'F36', 'F45', 'F46', 'F56']
    return Permutation([labels.index(label)+1 for label in keys]).to_permutation_group_element()

def are_cubics_same(cubic1, cubic2):
    mon = (sum(cubic1.P.gens()[0:4]) ^ 3).monomials()
    coeffs = matrix([[cubic1.coefficient(mn) for mn in mon], [cubic2.coefficient(mn) for mn in mon]]).minors(2)
    return list(set(coeffs)) == [0]

def get_permuted_L_set(perm):
    keys = ['E1', 'E2', 'E3', 'E4', 'E5', 'E6', 'G1', 'G2', 'G3', 'G4', 'G5', 'G6', 'F12', 'F13', 'F14', 'F15', 'F16', 'F23', 'F24', 'F25', 'F26', 'F34', 'F35', 'F36', 'F45', 'F46', 'F56']
    labels = from_perm_to_labels(perm)
    return [labels[keys.index(label)] for label in ['E1', 'G4', 'E2', 'G3', 'E3']]