class TritangentPlane:
    def __init__(self, plane, lines_dict, sing_locus, conditions=None):
        self.sing_locus = sing_locus
        self.plane = remove_sing_factors(plane, self.sing_locus)
        self.poly_ring = plane.parent()
        self.lines_dict = lines_dict
        self.labels = list(lines_dict.keys())
        self.lines = list(lines_dict.values())
        if conditions is None:
            self.conditions = self.find_conditions_for_eckardt_point()
        else:
            self.conditions = conditions

    def update_with_sostitution(self, sost, classified_lines, sing_locus):
        plane = self.poly_ring(self.plane.subs(sost).numerator())
        lines_dict = {key:classified_lines[key] for key in self.labels}
        conditions = self.conditions.subs(sost)
        conditions = self.poly_ring(conditions * conditions.denominator())
        if conditions != 0:
            conditions = remove_sing_factors(conditions, sing_locus)
        return TritangentPlane(plane, lines_dict, sing_locus, conditions)        

    def __eq__(self, other):
        if isinstance(other, TritangentPlane):
            return sorted(self.labels) == sorted(other.labels)
        return False

    def find_conditions_for_eckardt_point(self):
        P1 = self.lines[0].intersection_point(self.lines[1])
        P2 = self.lines[1].intersection_point(self.lines[2])
        mins = matrix([P1, P2]).minors(2)
        conditions = list(set([cond for cond in mins if cond != 0]))
        conditions1 = []
        for cond in conditions:
            conditions1.append(remove_sing_factors(prod([el[0] for el in list(cond.factor())]), self.sing_locus))        
        #print(gcd(conditions1))
        #print(P.ideal(conditions1).radical().primary_decomposition())
        return gcd(conditions1)

    def find_eckardt_point(self):
        if self.conditions == 0:
            return self.lines[1].intersection_point(self.lines[2])
        else:
            raise ValueError('Conditions are not 0')

    def __repr__(self):
        return self.labels.__repr__() +', '+ self.plane.__repr__()
