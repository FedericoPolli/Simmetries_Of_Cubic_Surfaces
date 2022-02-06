import os
import pickle
import multiprocessing as mp


class Cubic:
    def __init__(self, eqn, line, sing_locus, lines=None, cl_lines=None, tritangent_planes=None, eck_points=None):
        self.eqn = eqn
        self.P = eqn.parent()

        # want the sing_locus to be already factorized to speedup calculations
        if isinstance(sing_locus, Factorization):
            self.sing_locus = sing_locus
        else:
            self.sing_locus = sing_locus.factor()

        # if no parameters passed calculate everything
        if lines is None or cl_lines is None or tritangent_planes is None:
            self.lines = self.find_all_lines(line)
            self.cl_lines = self.classify_lines()
            self.lines = list(self.cl_lines.values())  # redefine self.lines in accordance with subs method
            self.tritangent_planes = self.find_tritangent_planes()
        else:
            self.lines = lines
            self.cl_lines = cl_lines
            self.tritangent_planes = tritangent_planes

        if eck_points is not None:
            self.eckardt_points = eck_points
        else:
            self.eckardt_points = [pl.find_eckardt_point() for pl in self.tritangent_planes if pl.has_eckardt_point()]
        self.eckardt_points_labels = [pl.labels for pl in self.tritangent_planes if pl.has_eckardt_point()]
        self.L_set_base = ('E1', 'G4', 'E2', 'G3', 'E3')

    def __str__(self):
        return self.eqn.__str__()

    def __repr__(self):
        return self.eqn.__repr__()

    def __repr_html__(self):
        return self.eqn.__repr_html__()

    def reduce(self, ideal):
        sing_locus = ideal.reduce(self.sing_locus.value())
        if sing_locus == 0:
            raise ValueError('Cubic is singular')
        sing_locus = sing_locus.factor()
        eqn = remove_sing_factors(ideal.reduce(self.eqn), sing_locus)
        eqn = eqn / (eqn.factor().unit())
        cl_lines = {key: line.reduce(ideal) for key, line in self.cl_lines.items()}
        lines = list(cl_lines.values())
        tritangent_planes = [pl.reduce(ideal, cl_lines, sing_locus) for pl in self.tritangent_planes]
        eck_points = ["Not Implemented"]
        return Cubic(eqn, lines[0], sing_locus, lines, cl_lines, tritangent_planes, eck_points)

        # update with sostitution all members of the class

    def subs(self, sost):
        sing_subs = self.sing_locus.value().subs(sost)
        if sing_subs == 0:
            raise ValueError('Cubic is singular')    
        # multiply by denominator^2 to avoid missing singular factors
        sing_locus = (self.P(sing_subs * (sing_subs.denominator()) ^ 2)).factor().radical()
        eqn = remove_sing_factors(self.P(self.eqn.subs(sost).numerator()), sing_locus)
        eqn = eqn / (eqn.factor().unit())
        cl_lines = {key: value.subs(sost) for key, value in self.cl_lines.items()}  # calls Line.subs()
        lines = list(cl_lines.values())  # to improve performance avoid resubbing
        tritangent_planes = [pl.update_with_sostitution(sost, cl_lines, sing_locus) for pl in self.tritangent_planes]
        return Cubic(eqn, lines[0], sing_locus, lines, cl_lines, tritangent_planes)

    def factor(self):
        return self.eqn.factor()

    # Finding lines ------------------------------------------------------------------------------------

    def find_all_lines(self, line):
        all_lines = [line]
        for i in range(3):
            all_lines += self._get_new_lines_on_cubic_surface(all_lines[i], all_lines)
        return all_lines

    # Gets possible new lines and filter to get only the actually new ones
    def _get_new_lines_on_cubic_surface(self, line, starting_lines):
        possible_new_lines = []
        possible_new_lines += self._find_lines_when_first_parameter_is_nonzero(line)
        possible_new_lines += self._find_lines_when_first_parameter_is_zero(line)
        new_lines = [line for line in possible_new_lines if line not in starting_lines]
        return new_lines

    # TBD, it is used to intersect cubic with pencil of planes
    def _find_lines_when_first_parameter_is_nonzero(self, line):
        conic, plane = self._define_equations(line, 0)
        P = self.P
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
            factors = [factor[0] for factor in factors_of_conic if
                       [v in factor[0].variables() for v in vr] != [False for _ in range(4)]]
            non_param_plane = plane.subs({P.gens()[-2]: first_param_root, P.gens()[-1]: second_param_root})
            new_lines.append(Line([non_param_plane, factors[0]]))
            new_lines.append(Line([non_param_plane, factors[1]]))
        return new_lines

    # TBD, it is used to intersect cubic with pencil of planes
    def _find_lines_when_first_parameter_is_zero(self, line):
        conic, plane = self._define_equations(line, 1)
        P = self.P
        vr = P.gens()[0:4]
        conic_discriminant = find_conic_discriminant(conic)
        new_lines = []
        for eqn in [factor[0] for factor in conic_discriminant]:
            if eqn.variables() == (P.gens()[-2],):
                factors_of_conic = conic.subs({P.gens()[-2]: 0, P.gens()[-1]: 1}).factor()
                non_param_plane = plane.subs({P.gens()[-2]: 0, P.gens()[-1]: 1})
                factors = [f[0] for f in factors_of_conic if
                           [v in f[0].variables() for v in vr] != [False for _ in range(4)]]
                new_lines.append(Line([non_param_plane, factors[0]]))
                new_lines.append(Line([non_param_plane, factors[1]]))
            else:
                continue
        return new_lines

    # TBD, it is used to intersect cubic with pencil of planes
    def _define_equations(self, line, k):
        plane = line.planes[0] * self.P.gens()[-2] + line.planes[1] * self.P.gens()[-1]
        if k==0:
            var = line.planes[0].variables()[0]
        else:
            var = [var for var in line.planes[k].variables() if var != line.planes[0].variables()[0]][0]
        param = [vr for vr in self.P.gens()[0:4] if vr != var]
        sol = solve_linear_system([plane], [var], param)
        sost = {param[i - 1]: sol[i] for i in range(1, 4)}
        sost[var] = sol[0]
        for fact in self.eqn.subs(sost).factor():
            for deg in fact[0].degrees()[0:4]:
                if deg > 1:
                    return fact[0], plane

    # Classify lines ------------------------------------------------------------------------------------

    # classifies the lines as E, G, F
    def classify_lines(self):
        E = self._find_E()
        G = self._find_G(E, [line for line in self.lines if line not in E])
        F = self._find_F(E, [line for line in self.lines if line not in E and line not in G])
        lines_dict = {'E1': E[0], 'E2': E[1], 'E3': E[2], 'E4': E[3], 'E5': E[4], 'E6': E[5], 'G1': G[0], 'G2': G[1],
                      'G3': G[2], 'G4': G[3], 'G5': G[4], 'G6': G[5], 'F12': F[0][1], 'F13': F[0][2], 'F14': F[0][3],
                      'F15': F[0][4], 'F16': F[0][5], 'F23': F[1][2], 'F24': F[1][3], 'F25': F[1][4], 'F26': F[1][5],
                      'F34': F[2][3], 'F35': F[2][4], 'F36': F[2][5], 'F45': F[3][4], 'F46': F[3][5], 'F56': F[4][5]}
        return lines_dict

    def _find_E(self):
        vrs = self.P.gens()[0:4]
        # Base L-set
        E1 = Line([vrs[1], vrs[2]])
        E2 = Line([vrs[0], vrs[3]])
        E3 = Line([vrs[0] - vrs[1], vrs[2] + vrs[3]])
        G3 = Line([vrs[0] - vrs[3], vrs[1] - vrs[2]])
        G4 = Line([vrs[0], vrs[1]])
        E4 = []
        # E4 is the only line not incident to E1, E2, E3, G4 intersecting G3
        for line in self.lines:
            incidence_relations = [line.are_incident(e) for e in [E1, E2, E3, G4]]
            if incidence_relations == [False, False, False, False] and line.are_incident(G3):
                E4 = line
                break
        E = self._get_other_skew_lines([E1, E2, E3, E4])
        # we want E5 to be the line whose last plucker coordinate is divisible by f
        if self.P.gens()[8].divides(E[-1].plucker[-1]):
            E[-2], E[-1] = E[-1], E[-2]
        return E

    def _get_other_skew_lines(self, skew_lines):
        for line in self.lines:
            if True not in [line.are_incident(e) for e in skew_lines]:
                skew_lines.append(line)
            if len(skew_lines) == 6:  # there are at most 6 skew lines
                break
        return skew_lines

    def _find_G(self, E, possible_lines):
        G = [0 for _ in range(6)]
        # G_i is the only line meeting all E_j except E_i
        for line in possible_lines:
            incidence_relations = [line.are_incident(e) for e in E]
            if incidence_relations.count(True) == 5:
                G[incidence_relations.index(False)] = line
                continue
        return G

    def _find_F(self, E, possible_lines):
        F = [[None for _ in range(6)] for _ in range(6)]
        # Fij is the only line meeting E_i, E_j but not other E_k
        for line in possible_lines:
            incidence_relations = [line.are_incident(e) for e in E]
            first_index = incidence_relations.index(True)
            second_index = incidence_relations.index(True, first_index + 1)
            F[first_index][second_index] = line
        return F

    # creates and populates all the 45 tritangent planes
    def find_tritangent_planes(self):
        all_triplets = find_all_triplets_of_coplanar_lines()
        planes = []
        for triplet in all_triplets:
            line1 = self.cl_lines[triplet[0]]
            line2 = self.cl_lines[triplet[1]]
            plane = line1.get_plane_containing_another_incident_line(line2)
            lines_dict = {k: self.cl_lines[k] for k in triplet}
            planes.append(TritangentPlane(plane, lines_dict, self.sing_locus))
        return planes

    # Find projectivities -----------------------------------------------------------------------------

    # find the permutations which send the set of this cubic's Eckardt points to itself
    def find_admissible_permutations(self):
        with open('all_permutations.pickle', 'rb') as fil:
            all_permutations = pickle.load(fil)
        adm_perm = []
        keys = list(self.cl_lines.keys())
        labels = [sorted(label) for label in self.eckardt_points_labels]
        for perm_label in all_permutations:
            flag = True
            # check if each permuted label is still in Eckardt points labels
            for triple in labels:
                new_triple = [perm_label[keys.index(label)] for label in triple]
                if sorted(new_triple) not in labels:
                    flag = False
                    break
            if flag is True:
                adm_perm.append(from_labels_to_perm(perm_label))
        return adm_perm

    # given a list of permutations of the 27 lines, it finds the permuted base L-sets
    # and finds the associated projectivities
    def find_admissible_projectivities(self, adm_perm=None):
        if adm_perm is None:
            adm_perm = self.find_admissible_permutations()
        resulting_L_sets = list(set([get_permuted_L_set(perm) for perm in adm_perm]))
        return self.find_all_proj_parallel(resulting_L_sets)

    # from the classification as E, F, G returns the actual lines in plucker coordinates
    def get_L_set_in_plucker(self, L_set):
        return tuple(map(lambda uu: self.cl_lines[uu], L_set))

    # finds all the projectivities from the given L-sets in parallel
    def find_all_proj_parallel(self, all_L_sets):
        if os.name == "nt":
            mp.freeze_support()
        pool = mp.Pool(mp.cpu_count() - 1)
        all_param = ((self.L_set_base, L_set) for L_set in all_L_sets)
        result = pool.map(self.find_proj_parallel_wrapper, all_param)
        pool.close()
        return result

    def find_proj_parallel_wrapper(self, args):
        return self.find_projectivity(*args)

    def find_projectivity(self, L_set1, L_set2):
        L1 = self.get_L_set_in_plucker(L_set1)
        L2 = self.get_L_set_in_plucker(L_set2)
        M = find_projectivity(L1, L2)
        return M

    # Find simmetries -------------------------------------------------------------------------------------

    # finds which projectivities is also a simmetry of this cubic
    def find_simmetries(self, projectivities=None):
        if projectivities is None:
            projectivities = self.find_admissible_projectivities()
        return self.find_simmetries_parallel(projectivities)

    # find simmetries in parallel
    def find_simmetries_parallel(self, all_projectivities):
        if os.name == "nt":
            mp.freeze_support()
        pool = mp.Pool(mp.cpu_count() - 1)
        all_param = ((proj,) for proj in all_projectivities)
        result = [el for el in pool.map(self.find_simmetries_wrapper, all_param) if el is not None]
        pool.close()
        return result

    def find_simmetries_wrapper(self, args):
        return self.find_simmetry(*args)

    # check if the list of coefficients of this cubic and the new one are proportional
    def find_simmetry(self, proj):
        sost = change_coordinates(proj)
        new_cubic = self.eqn.subs(sost)
        if self.are_cubics_same(new_cubic):
            return proj
        else:
            return None

    def are_cubics_same(self, other):
        mon = (sum(self.P.gens()[0:4]) ^ 3).monomials()
        coeffs = matrix([[self.coefficient(mn) for mn in mon], [other.coefficient(mn) for mn in mon]]).minors(2)
        for el in coeffs:
            if el != 0:
                return False
        return True

    # --------------------------------------------------------------------------------------------------

    def find_conditions_on_L_sets(self, possible_L_sets, filter_ideals = "", factors = None):
        L_set_ideals = []
        all_ideals = []
        for L_set in possible_L_sets:
            proj = self.find_projectivity(L_set_base, L_set)
            ideals = self.find_decomposed_conditions_on_cubic(proj, factors)
            if filter_ideals == "eck":
                ideals = [ide for ide in ideals if self.is_ideal_valid(ide)]
            else:
                ideals = [ide for ide in ideals if not self.is_ideal_singular(ide)]
            all_ideals += ideals
            L_set_ideals.append([L_set, ideals])
        return L_set_ideals, all_ideals
    
    def find_conditions_on_cubic(self, proj):
        #change coordinates, obtain new cubic and find conditions on the parameters
        vrs = self.P.gens()[0:4]
        mon = (sum(vrs) ^ 3).monomials()
        sost = change_coordinates(proj)
        new_cubic = self.eqn.subs(sost)
        current_conds = matrix([[self.eqn.coefficient(mn) for mn in mon], [new_cubic.coefficient(mn) for mn in mon]]).minors(2)
        no_sing_conds = [remove_sing_factors(el, self.sing_locus) for el in list(set(current_conds)) if el !=0]
        return self.P.ideal(no_sing_conds)

    def find_decomposed_conditions_on_cubic(self, proj, factors = None):
        if factors is None:
            factors = self.P.gens()[5]
        ideal = self.find_conditions_on_cubic(proj)
        return ideal.saturation(factors)[0].radical().primary_decomposition('gtz')
    
    
    # studies the projectivities which do not send this cubic to itself in general to try to find
    # values for the parameters for which some of these projectivities become simmetries
    def find_conditions_for_subfamilies(self, projectivities, simmetries=[]):
        projs = [proj for proj in projectivities if proj not in simmetries]
        return self.find_conditions_parallel(projs)
    
    
    def find_conditions_parallel(self, projs):
        if os.name == "nt":
            mp.freeze_support()
        pool = mp.Pool(mp.cpu_count() - 1)
        all_param = ((proj,) for proj in projs)
        result = [el for el in pool.map(self.find_conditions_wrapper, all_param) if el is not None]
        pool.close()
        return list(set([item for sublist in result for item in sublist])) 
    

    def find_conditions_wrapper(self, args):
        return self.find_conditions(*args)
    
    
    def find_conditions(self, proj):
        prim_deco = self.find_decomposed_conditions_on_cubic(proj)
        return [ideale for ideale in prim_deco if self.is_ideal_valid(ideale)]


    # check if ideal does not contain singular locus and  if it does not contain
    # polynomials which would cause this cubic to have more eckardt points
    def is_ideal_valid(self, ideal):
        if self.is_ideal_singular(ideal):
            return False
        if self.does_ideal_give_more_eck_points(ideal):
            return False
        return True

    def is_ideal_singular(self, ideal):
        return self.sing_locus.value() in ideal

    def does_ideal_give_more_eck_points(self, ideal):
        for poly in list(set([pl.condition for pl in self.tritangent_planes if pl.condition != 0])):
            if poly in ideal:
                return True
        return False

    # applies simmetry to eckardt points and returns the associated permutation
    def apply_proj_to_eck(self, proj):
        new_indices = []
        eck = self.eckardt_points
        for i in range(len(eck)):
            new_indices.append(eck.index(eck[i] * proj) + 1)
        return new_indices

    # applies simmetriy to lines and returns the associated permutation
    def apply_proj_to_lines(self, proj):
        new_indices = []
        for i in range(len(self.lines)):
            new_indices.append(self.lines.index(self.lines[i].apply_proj(proj)) + 1)
        return new_indices

    def coefficients(self):
        return self.eqn.coefficients()

    def coefficient(self, var):
        return self.eqn.coefficient(var)
