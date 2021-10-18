class Cubic():
    def __init__(self, eqn, line, sing_cubic, lines = None, cl_lines = None, tritangent_planes = None):
        self.eqn = eqn
        self.P = eqn.parent()
        self.sing_cubic = sing_cubic
        if lines is None or cl_lines is None or tritangent_planes is None:
            self.lines = self.find_all_lines_on_cubic_surface(line)
            self.cl_lines = self.classify_lines_on_cubic_surface()
            self.tritangent_planes = self.find_tritangent_planes(self.sing_cubic.factor())
        else:
            self.lines = lines
            self.cl_lines = cl_lines
            self.tritangent_planes = tritangent_planes
        self.eckardt_points = [pl.find_eckardt_point() for pl in self.tritangent_planes if pl.conditions == 0]            
        
    def __str__(self):
        return self.eqn.__str__()
        
    def __repr__(self):
        return self.eqn.__repr__()
        
    def __repr_html__(self):
        return self.eqn.__repr_html__()

    def subs(self, sost):
        sing_cubic = self.P(self.sing_cubic.subs(sost).numerator())
        eqn = remove_sing_factors(self.P(self.eqn.subs(sost).numerator()), sing_cubic.factor())
        lines = [line.subs(sost) for line in self.lines]
        cl_lines = {key:value.subs(sost) for key, value in self.cl_lines.items()}
        for pl in self.tritangent_planes:
            pl.set_sing_cubic(sing_cubic.factor())
        tritangent_planes = [pl.subs(sost) for pl in self.tritangent_planes]
        return Cubic(eqn, None, sing_cubic, lines, cl_lines, tritangent_planes)

    def factor(self):
        return self.eqn.factor()
    
    def find_all_lines_on_cubic_surface(self, line):
        all_lines = [line]
        #First line was already done
        for i in range(3):
            all_lines += self._get_new_lines_on_cubic_surface(all_lines[i], all_lines)
            if len(all_lines) == 27:
                break
        return all_lines

    #Starting from the given lines, it tries to find new lines
    def _get_new_lines_on_cubic_surface(self, line, starting_lines):
        possible_new_lines = []
        possible_new_lines += self._find_lines_when_first_parameter_is_nonzero(line)
        possible_new_lines += self._find_lines_when_first_parameter_is_zero(line)
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
            factors = [factor[0] for factor in factors_of_conic if [v in factor[0].variables() for v in vr] != [False for i in range(4)]]
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
                factors = [f[0] for f in factors_of_conic if [v in f[0].variables() for v in vr] != [False for i in range(4)]]
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
        sost = {param[i-1]: sol[i] for i in range(1,4)}
        sost[var]= sol[0]
        for fact in self.eqn.subs(sost).factor():
            for deg in fact[0].degrees()[0:4]:
                if deg>1:
                    return fact[0], plane
                
# Classify lines ---------------------------------------------------------------------------------

    def classify_lines_on_cubic_surface(self):
        E = self._find_E()
        G = self._find_G(E, [line for line in self.lines if line not in E])
        if (e^2).divides(G[-1].plucker[-1]):
            E[-2], E[-1] = E[-1], E[-2]
            G[-2], G[-1] = G[-1], G[-2]
        F = self._find_F(E, [line for line in self.lines if line not in E and line not in G])
        lines_dict = dict(E1=E[0], E2=E[1], E3=E[2], E4=E[3], E5=E[4], E6=E[5],
                      G1=G[0], G2=G[1], G3=G[2], G4=G[3], G5=G[4], G6=G[5],
                      F12=F[0][1], F13=F[0][2], F14=F[0][3], F15=F[0][4], F16=F[0][5],
                      F23=F[1][2], F24=F[1][3], F25=F[1][4], F26=F[1][5], F34=F[2][3], 
                      F35=F[2][4], F36=F[2][5], F45=F[3][4], F46=F[3][5], F56=F[4][5])
        return lines_dict

    def _find_E(self):
        vrs = self.P.gens()[0:4]
        E1 = Line([vrs[1],vrs[2]])
        E2 = Line([vrs[0],vrs[3]])
        E3 = Line([vrs[0]-vrs[1], vrs[2]+vrs[3]])
        G3 = Line([vrs[0]-vrs[3], vrs[1]-vrs[2]])
        G4 = Line([vrs[0],vrs[1]])
        E4=[]
        for line in self.lines:
            incidence_relations = [line.are_incident(e) for e in [E1, E2, E3, G4]]
            if incidence_relations == [False, False, False, False] and line.are_incident(G3):
                E4 = line
                break
        E = self._get_other_skew_lines([E1, E2, E3, E4])
        #testing
        if f.divides(E[-1].plucker[-1]):
            E[-2], E[-1] = E[-1], E[-2]
        return E

    def _get_other_skew_lines(self, skew_lines):
        for line in self.lines:
            is_line_skew = true
            for e in skew_lines:
                if line.are_incident(e):
                    is_line_skew = false
            if is_line_skew:       
                skew_lines.append(line)
            if len(skew_lines) == 6:
                break
        return skew_lines

    def _find_G(self, E, possible_lines):
        G = [0 for i in range(6)]
        for line in possible_lines:
            incidence_relations = [line.are_incident(e) for e in E]
            if incidence_relations.count(true) == 5:
                G[incidence_relations.index(false)] = line
                continue 
        return G

    def _find_F(self, E, possible_lines):
        F = [[None for i in range(6)] for j in range(6)]
        for line in possible_lines:
            incidence_relations = [line.are_incident(e) for e in E]
            first_index = incidence_relations.index(true)
            second_index = incidence_relations.index(true, first_index+1)
            F[first_index][second_index] = line
        return F
    
    def find_tritangent_planes(self, sing_cubic):
        keys = self.cl_lines.keys()
        all_triplets = find_all_triplets_of_coplanar_lines(keys)
        planes = []
        for triplet in all_triplets:
            line1 = self.cl_lines.get(triplet[0])
            line2 = self.cl_lines.get(triplet[1])        
            plane = get_plane_containing_two_incident_lines(line1, line2)
            lines_dict = {k:self.cl_lines.get(k) for k in triplet}
            planes.append(tritangent_plane(plane, lines_dict, sing_cubic))
        return planes
        
    def find_simmetries(self, projectivities):
        return find_simmetries_parallel(self.eqn, projectivities)

    def find_simmetries_np(self, projectivities):
        simm = []
        vrs = vector(self.P.gens()[0:4])
        for proj in projectivities:
            change_coord = vector(vrs)*proj 
            sost = {vrs[i] : change_coord[i] for i in range(4)}
            new_cubic = self.eqn.subs(sost)
            if self.eqn.coefficient(vrs[0]^2*vrs[1])*new_cubic - new_cubic.coefficient(vrs[0]^2*vrs[1])*self.eqn == 0:
                simm.append(proj)
        return simm
        
    def find_admissible_projectivities(self):
        adm_perm = self.find_admissible_permutations() 
        resulting_L_sets = []
        for perm in adm_perm:
            perm_L_set = [perm[list(self.cl_lines.keys()).index(label)] for label in ['E1', 'G4', 'E2', 'G3', 'E3']]
            if perm_L_set not in resulting_L_sets:
                resulting_L_sets.append(perm_L_set)
        return find_all_proj_parallel(self.cl_lines, resulting_L_sets) 


    def find_admissible_permutations(self):
        with open('all_permutations.pickle', 'rb') as ff:
            all_permutations = pickle.load(ff) 
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