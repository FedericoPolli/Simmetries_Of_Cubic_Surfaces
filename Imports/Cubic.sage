import os
import pickle
import multiprocessing as mp


class Cubic:
    def __init__(self, eqn, line, sing_cubic, lines=None, cl_lines=None, tritangent_planes=None):
        self.eqn = eqn
        self.P = eqn.parent()
        if isinstance(sing_cubic, Factorization):
            self.sing_cubic = sing_cubic
        else:
            self.sing_cubic = sing_cubic.factor()
        if lines is None or cl_lines is None or tritangent_planes is None:
            self.lines = self.find_all_lines_on_cubic_surface(line)
            self.cl_lines = self.classify_lines_on_cubic_surface()
            self.tritangent_planes = self.find_tritangent_planes()
        else:
            self.lines = lines
            self.cl_lines = cl_lines
            self.tritangent_planes = tritangent_planes
        self.eckardt_points = [pl.find_eckardt_point() for pl in self.tritangent_planes if pl.conditions == 0]
        self.eckardt_points_labels = [pl.labels for pl in self.tritangent_planes if pl.conditions == 0]

    def __str__(self):
        return self.eqn.__str__()

    def __repr__(self):
        return self.eqn.__repr__()

    def __repr_html__(self):
        return self.eqn.__repr_html__()

    def subs(self, sost):
        sing_subs = self.sing_cubic.value().subs(sost)
        sing_cubic = (self.P(sing_subs * (sing_subs.denominator()) ^ 2)).factor()
        eqn = remove_sing_factors(self.P(self.eqn.subs(sost).numerator()), sing_cubic)
        lines = [line.subs(sost) for line in self.lines]
        cl_lines = {key: value.subs(sost) for key, value in self.cl_lines.items()}
        for pl in self.tritangent_planes:
            pl.set_sing_cubic(sing_cubic)
        tritangent_planes = [pl.subs(sost) for pl in self.tritangent_planes]
        return Cubic(eqn, None, sing_cubic, lines, cl_lines, tritangent_planes)

    def factor(self):
        return self.eqn.factor()

    def find_all_lines_on_cubic_surface(self, line):
        all_lines = [line]
        # First line was already done
        for i in range(3):
            all_lines += self._get_new_lines_on_cubic_surface(all_lines[i], all_lines)
            if len(all_lines) == 27:
                break
        return all_lines

    # Starting from the given lines, it tries to find new lines
    def _get_new_lines_on_cubic_surface(self, line, starting_lines):
        possible_new_lines = []
        possible_new_lines += self._find_lines_when_first_parameter_is_nonzero(line)
        possible_new_lines += self._find_lines_when_first_parameter_is_zero(line)
        new_lines = []
        for line1 in possible_new_lines:
            flag = True
            for line2 in starting_lines:
                if line1 == line2:
                    flag = False
                    break
            if flag:
                new_lines.append(line1)
        return new_lines

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

    def _define_equations(self, line, k):
        plane = line.planes[0] * self.P.gens()[-2] + line.planes[1] * self.P.gens()[-1]
        var = line.planes[k].variables()[0]
        param = [vr for vr in self.P.gens()[0:4] if vr != var]
        sol = solve_linear_system([plane], [var], param)
        sost = {param[i - 1]: sol[i] for i in range(1, 4)}
        sost[var] = sol[0]
        for fact in self.eqn.subs(sost).factor():
            for deg in fact[0].degrees()[0:4]:
                if deg > 1:
                    return fact[0], plane

    # Classify lines ---------------------------------------------------------------------------------

    def classify_lines_on_cubic_surface(self):
        E = self._find_E()
        G = self._find_G(E, [line for line in self.lines if line not in E])
        F = self._find_F(E, [line for line in self.lines if line not in E and line not in G])
        lines_dict = dict(E1=E[0], E2=E[1], E3=E[2], E4=E[3], E5=E[4], E6=E[5],
                          G1=G[0], G2=G[1], G3=G[2], G4=G[3], G5=G[4], G6=G[5],
                          F12=F[0][1], F13=F[0][2], F14=F[0][3], F15=F[0][4], F16=F[0][5],
                          F23=F[1][2], F24=F[1][3], F25=F[1][4], F26=F[1][5], F34=F[2][3],
                          F35=F[2][4], F36=F[2][5], F45=F[3][4], F46=F[3][5], F56=F[4][5])
        return lines_dict

    def _find_E(self):
        vrs = self.P.gens()[0:4]
        E1 = Line([vrs[1], vrs[2]])
        E2 = Line([vrs[0], vrs[3]])
        E3 = Line([vrs[0] - vrs[1], vrs[2] + vrs[3]])
        G3 = Line([vrs[0] - vrs[3], vrs[1] - vrs[2]])
        G4 = Line([vrs[0], vrs[1]])
        E4 = []
        for line in self.lines:
            incidence_relations = [line.are_incident(e) for e in [E1, E2, E3, G4]]
            if incidence_relations == [False, False, False, False] and line.are_incident(G3):
                E4 = line
                break
        E = self._get_other_skew_lines([E1, E2, E3, E4])
        if self.P.gens()[-3].divides(E[-1].plucker[-1]):
            E[-2], E[-1] = E[-1], E[-2]
        return E

    def _get_other_skew_lines(self, skew_lines):
        for line in self.lines:
            is_line_skew = True
            for e in skew_lines:
                if line.are_incident(e):
                    is_line_skew = False
            if is_line_skew:
                skew_lines.append(line)
            if len(skew_lines) == 6:
                break
        return skew_lines

    def _find_G(self, E, possible_lines):
        G = [0 for _ in range(6)]
        for line in possible_lines:
            incidence_relations = [line.are_incident(e) for e in E]
            if incidence_relations.count(True) == 5:
                G[incidence_relations.index(False)] = line
                continue
        return G

    def _find_F(self, E, possible_lines):
        F = [[None for _ in range(6)] for _ in range(6)]
        for line in possible_lines:
            incidence_relations = [line.are_incident(e) for e in E]
            first_index = incidence_relations.index(True)
            second_index = incidence_relations.index(True, first_index + 1)
            F[first_index][second_index] = line
        return F

    def find_tritangent_planes(self):
        all_triplets = find_all_triplets_of_coplanar_lines()
        planes = []
        for triplet in all_triplets:
            line1 = self.cl_lines.get(triplet[0])
            line2 = self.cl_lines.get(triplet[1])
            plane = get_plane_containing_two_incident_lines(line1, line2)
            lines_dict = {k: self.cl_lines.get(k) for k in triplet}
            planes.append(TritangentPlane(plane, lines_dict, self.sing_cubic))
        return planes

    def find_simmetries(self, projectivities=None):
        if projectivities is None:
            projectivities = self.find_admissible_projectivities()
        return self.find_simmetries_parallel(projectivities)

    def find_simmetry(self, proj):
        vrs = self.P.gens()[0:4]
        mon = (sum(vrs) ^ 3).monomials()
        change_coord = vector(vrs) * proj
        sost = {vrs[i]: change_coord[i] for i in range(4)}
        new_cubic = self.eqn.subs(sost)
        coeffs = matrix([[self.eqn.coefficient(mn) for mn in mon], [new_cubic.coefficient(mn) for mn in mon]]).minors(2)
        if coeffs == [0 for _ in range(len(coeffs))]:
            return proj
        else:
            return None

    def find_simmetries_wrapper(self, args):
        return self.find_simmetry(*args)

    def find_simmetries_parallel(self, all_projectivities):
        if os.name == "nt":
            mp.freeze_support()
        pool = mp.Pool(mp.cpu_count() - 1)
        all_param = ((proj,) for proj in all_projectivities)
        result = [el for el in pool.map(self.find_simmetries_wrapper, all_param) if el is not None]
        pool.close()
        return result

    def find_admissible_projectivities(self, adm_perm=None):
        if adm_perm is None:
            adm_perm = self.find_admissible_permutations()
        resulting_L_sets = []
        for perm in adm_perm:
            perm_L_set = [perm[list(self.cl_lines.keys()).index(label)] for label in ['E1', 'G4', 'E2', 'G3', 'E3']]
            if perm_L_set not in resulting_L_sets:
                resulting_L_sets.append(perm_L_set)
        return self.find_all_proj_parallel(resulting_L_sets)

    # from the classification as E, F, G returns the actual lines in plucker coordinates
    def get_L_set_in_plucker(self, L_set):
        return tuple(map(lambda uu: self.cl_lines[uu], L_set))

    def find_all_projectivities(self, L_set, base_five_points):
        L2 = self.get_L_set_in_plucker(L_set)
        M = find_projectivity(base_five_points, L2)
        return M

    def find_proj_parallel(self, args):
        return self.find_all_projectivities(*args)

    def find_all_proj_parallel(self, all_L_sets):
        if os.name == "nt":
            mp.freeze_support()
        pool = mp.Pool(mp.cpu_count() - 1)
        L_set_base = self.get_L_set_in_plucker(['E1', 'G4', 'E2', 'G3', 'E3'])
        base_five_points = get_five_points_in_general_position(L_set_base)
        all_param = ((L_set, base_five_points) for L_set in all_L_sets)
        result = pool.map(self.find_proj_parallel, all_param)
        pool.close()
        return result

    def find_admissible_permutations(self):
        with open('all_permutations.pickle', 'rb') as fil:
            all_permutations = pickle.load(fil)
        adm_perm = []
        keys = list(self.cl_lines.keys())
        labels = [sorted(el.labels) for el in self.tritangent_planes if el.conditions == 0]
        for perm in all_permutations:
            flag = True
            for triple in labels:
                new_triple = []
                for label in triple:
                    new_triple.append(perm[keys.index(label)])
                if sorted(new_triple) not in labels:
                    flag = False
                    break
            if flag is True:
                adm_perm.append(perm)
        return adm_perm

    def find_conditions_for_subfamilies(self, projectivities, simmetries):
        vrs = self.P.gens()[0:4]
        mon = (sum(vrs) ^ 3).monomials()
        conditions = []
        for M in [proj for proj in projectivities if proj not in simmetries]:
            sost = change_coord(M)
            new_cubic = remove_sing_factors(self.eqn.subs(sost), self.sing_cubic)
            minor = list(set(
                matrix([[new_cubic.coefficient(mn) for mn in mon], [self.eqn.coefficient(mn) for mn in mon]]).minors(2)))
            minor = [remove_sing_factors(el, self.sing_cubic) for el in minor if el != 0]
            prim_deco = self.P.ideal(minor).radical().primary_decomposition()
            for ideale in prim_deco:
                if self.is_ideal_valid(ideale):
                    conditions.append(ideale.gens())
        return list(set(conditions))

    def is_ideal_valid(self, ideal):
        if self.sing_cubic.value() in ideal:
            return False
        for poly in list(set([pl.conditions for pl in self.tritangent_planes if pl.conditions != 0])):
            if poly in ideal:
                return False
        return True
