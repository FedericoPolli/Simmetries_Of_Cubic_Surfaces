class tritangent_plane():
    def __init__(self, plane, lines_dict, sing_cubic, conditions = None):
        self.sing_cubic = sing_cubic
        self.plane = remove_sing_factors(plane, self.sing_cubic)
        self.P = plane.parent()
        self.lines_dict = lines_dict
        self.labels = list(lines_dict.keys())
        self.lines = list(lines_dict.values())
        if conditions is None:
            self.conditions = self.find_conditions_for_eckardt_point()
        else:
            self.conditions = conditions

    def set_sing_cubic(self, sing_cubic):
        self.sing_cubic = sing_cubic
    
    def subs(self, sost):
        plane = self.P(self.plane.subs(sost).numerator())
        lines_dict = {key:value.subs(sost) for key, value in self.lines_dict.items()}
        conditions = self.conditions.subs(sost)
        conditions = self.P(conditions*conditions.denominator())
        if conditions != 0:
            conditions = remove_sing_factors(conditions, self.sing_cubic)
        return tritangent_plane(plane, lines_dict, self.sing_cubic, conditions)
    
    def __eq__(self, other):
        if isinstance(other, tritangent_plane):
            return self.plane == other.plane
        return False
    
    def find_conditions_for_eckardt_point(self):
        P1 = self.lines[0].intersection_point(self.lines[1])
        P2 = self.lines[1].intersection_point(self.lines[2])
        mins = matrix([P1, P2]).minors(2)
        conditions = [cond for cond in mins if cond != 0]
        
        #remove higher powers
        conditions1 = []
        for cond in conditions:
            conditions1.append(prod([el[0] for el in list(cond.factor())]))
        condition = gcd(conditions1)
        
        #remove singular factors
        if condition == 0:
            return condition
        else:
            return remove_sing_factors(condition, self.sing_cubic)         
        
    def find_eckardt_point(self):
        if self.conditions == 0:
            return self.lines[1].intersection_point(self.lines[2])
        else:
            raise ValueError('Conditions are not 0')
    
    def __repr__(self):
        return self.labels.__repr__()+self.plane.__repr__()