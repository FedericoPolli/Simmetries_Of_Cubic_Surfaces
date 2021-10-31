class Point():  
    def __init__(self, components):
        self.P = components[0].parent()
        mcd = gcd(components)
        temp = vector([comp//mcd for comp in components])
        temp = temp/gcd([gcd(el.coefficients()) for el in temp])
        if temp[temp.nonzero_positions()[0]] < 0:
            self.components = -temp
        else:
            self.components = temp
    
    def subs(self, sost):
        components = vector([comp.subs(sost)for comp in self.components])
        return Point([self.P(el) for el in components*components.denominator()])
    
    def __str__(self):
        return self.components.__str__()
        
    def __repr__(self):
        return self.components.__repr__()
        
    def __repr_html__(self):
        return self.components.__repr_html__()
     
    def __iter__(self):
        return (el for el in self.components)
        
    def __len__(self):
        return len(self.components)
        
    def __getitem__(self, item):
        return self.components.__getitem__(item)
    
    def __eq__(self, other_point):
        if isinstance(other_point, Point):
            return are_vectors_proportional(self.components, other_point.components)
        return False
        
    def __mul__(self, other):
        return Point(self.components * other)
        
    def __rmul__(self, other):
        return Point(other * self.components)
        
    def parent(self):
        return self.P