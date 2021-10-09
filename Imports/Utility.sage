# Finding lines------------------------------------------------------------------------


def find_all_lines_on_cubic_surface(line, surface):
    all_lines = [line]
    #First line was already done
    for i in range(3):
        all_lines += get_new_lines_on_cubic_surface(all_lines[i], surface, all_lines)
        if len(all_lines) == 27:
            break
    return all_lines


#Starting from the given lines, it tries to find new lines
def get_new_lines_on_cubic_surface(line, surface, starting_lines):
    possible_new_lines = []
    possible_new_lines += find_lines_when_first_parameter_is_nonzero(line, surface)
    possible_new_lines += find_lines_when_first_parameter_is_zero(line, surface)
    new_lines = []
    for l in possible_new_lines:
        flag = True
        for l2 in starting_lines:
            if l == l2:
                flag = False
                break
        if flag:
            new_lines.append(l)
    return new_lines

def find_lines_when_first_parameter_is_nonzero(line, surface):
    P, conic, plane = define_equations(line, surface, 0)
    vr = P.gens()[0:4]
    conic_discriminant = find_conic_discriminant(conic)
    new_lines = []
    for eqn in [factor[0] for factor in conic_discriminant]:
        # exclude case when the eliminated variable is zero
        if eqn.variables() == (P.gens()[-2],):
            continue
        # avoid factors without any parameter
        if P.gens()[-2] and P.gens()[-1] not in eqn.variables():
            continue
        first_param_root = eqn.coefficient(P.gens()[-1])
        second_param_root = -eqn.coefficient(P.gens()[-2])
        factors_of_conic = conic.subs({P.gens()[-2]: first_param_root, P.gens()[-1]: second_param_root}).factor()
        factors = [factor[0] for factor in factors_of_conic if [v in factor[0].variables() for v in vr] != [False for i in range(4)]]
        non_param_plane = plane.subs({P.gens()[-2]: first_param_root, P.gens()[-1]: second_param_root})
        new_lines.append(Line([non_param_plane, factors[0]]))
        new_lines.append(Line([non_param_plane, factors[1]]))
    return new_lines
        
        
def find_lines_when_first_parameter_is_zero(line, surface):
    P, conic, plane = define_equations(line, surface, 1)
    vr = P.gens()[0:4]
    conic_discriminant = find_conic_discriminant(conic)
    new_lines = []
    for eqn in [factor[0] for factor in conic_discriminant]:
        if eqn.variables() == (P.gens()[-2],):
            factors_of_conic = conic.subs({P.gens()[-2]: 0, P.gens()[-1]: 1}).factor()
            non_param_plane = plane.subs({P.gens()[-2]: 0, P.gens()[-1]: 1})
            factors = [f[0] for f in factors_of_conic if [v in f[0].variables() for v in vr] != [False for i in range(4)]]
            new_lines.append(Line([non_param_plane, factors[0]]))
            new_lines.append(Line([non_param_plane, factors[1]]))
        else:
            continue
    return new_lines


def define_equations(line, surface, k):
    P = surface.parent()
    plane = line.planes[0] * P.gens()[-2] + line.planes[1] * P.gens()[-1]
    var = line.planes[k].variables()[0]
    param = [vr for vr in P.gens()[0:4] if vr != var]
    sol = solve_linear_system([plane], [var], param)
    sost = {param[i-1]: sol[i] for i in range(1,4)}
    sost[var]= sol[0]
    for fact in surface.subs(sost).factor():
        for deg in fact[0].degrees()[0:4]:
            if deg>1:
                return P, fact[0], plane
  
  
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


# Classify lines ---------------------------------------------------------------------------------


def classify_lines_on_cubic_surface(all_lines):
    E = find_E(all_lines)
    G = find_G(E, [line for line in all_lines if line not in E])
    F = find_F(E, [line for line in all_lines if line not in E and line not in G])
    lines_dict = dict(E1=E[0], E2=E[1], E3=E[2], E4=E[3], E5=E[4], E6=E[5],
                  G1=G[0], G2=G[1], G3=G[2], G4=G[3], G5=G[4], G6=G[5],
                  F12=F[0][1], F13=F[0][2], F14=F[0][3], F15=F[0][4], F16=F[0][5],
                  F23=F[1][2], F24=F[1][3], F25=F[1][4], F26=F[1][5], F34=F[2][3], 
                  F35=F[2][4], F36=F[2][5], F45=F[3][4], F46=F[3][5], F56=F[4][5])
    return lines_dict


def find_E(all_lines):
    P = all_lines[0].P
    E1 = Line([y,z])
    E2 = Line([x,t])
    E3 = Line([x-y, z+t])
    G3 = Line([x-t, y-z])
    G4 = Line([x,y])
    E4=[]
    for line in all_lines:
        incidence_relations = [line.are_incident(e) for e in [E1, E2, E3, G4]]
        if incidence_relations == [False, False, False, False] and line.are_incident(G3):
            E4 = line
            break
    E = get_other_skew_lines([E1, E2, E3, E4], all_lines)
    return E


def get_other_skew_lines(skew_lines, all_lines):
    for line in all_lines:
        is_line_skew = true
        for e in skew_lines:
            if line.are_incident(e):
                is_line_skew = false
        if is_line_skew:       
            skew_lines.append(line)
        if len(skew_lines) == 6:
            break
    return skew_lines


def find_G(E, possible_lines):
    G = [0 for i in range(6)]
    for line in possible_lines:
        incidence_relations = [line.are_incident(e) for e in E]
        if incidence_relations.count(true) == 5:
            G[incidence_relations.index(false)] = line
            continue 
    return G
 
    
def find_F(E, possible_lines):
    F = [[None for i in range(6)] for j in range(6)]
    for line in possible_lines:
        incidence_relations = [line.are_incident(e) for e in E]
        first_index = incidence_relations.index(true)
        second_index = incidence_relations.index(true, first_index+1)
        F[first_index][second_index] = line
    return F


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
    line_P_D = Line([plane_l2_l5, plane_l3_l4])
    Q = line_P_D.intersection_point(L_set[4])
    return A,B,C,E,Q


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
                        
            
def find_all_tritangent_planes(cl_lines):
    keys = cl_lines.keys()
    all_triplets = find_all_triplets_of_coplanar_lines(keys)
    planes = []
    for triplet in all_triplets:
        line1 = cl_lines.get(triplet[0])
        line2 = cl_lines.get(triplet[1])        
        plane = get_plane_containing_two_incident_lines(line1, line2)
        lines_dict = {k:cl_lines.get(k) for k in triplet}
        planes.append(tritangent_plane(plane, lines_dict))
    return planes


def find_all_triplets_of_coplanar_lines(keys):
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
    fact = poly.factor()
    factors = [el[0] for el in list(fact)]
    mcd = [el[0] for el in list(gcd(fact, sing))]
    return prod([el for el in factors if el not in mcd])
    
    
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