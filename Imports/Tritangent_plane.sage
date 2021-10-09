class tritangent_plane():
    def __init__(self, plane, lines_dict, sing_cubic):
        self.plane = plane
        self.P = plane.parent()
        self.lines_dict = lines_dict
        self.labels = list(lines_dict.keys())
        self.lines = list(lines_dict.values())
        self.conditions = self.find_conditions_for_eckardt_point(sing_cubic)

    def __eq__(self, other):
        if isinstance(other, tritangent_plane):
            return self.plane == other.plane
        return False
    
    def find_conditions_for_eckardt_point(self, sing_cubic):
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
            return remove_sing_factors(condition, sing_cubic)         
        
    def find_eckardt_point(self):
        if self.conditions == 0:
            return self.lines[1].intersection_point(self.lines[2])
        else:
            raise ValueError('Conditions are not 0')
    
    def __repr__(self):
        return self.plane.__repr__()