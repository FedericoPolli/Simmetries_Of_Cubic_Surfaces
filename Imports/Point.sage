class Point:   
    """
    Class that represents a point in any projective space, 
    i.e. a vector defined up to multiplication with a nonzero scalar.
    """

    def __init__(self, components):
        self._poly_ring = components[0].parent()
        mcd = gcd(components)
        #divide each component by gcd with exact division
        temp = vector([comp//mcd for comp in components])
        
        # calculate the gcd of the coefficients of each component
        #  as gcd of polynomials returns 1 as coefficient
        temp = temp/gcd([gcd(el.coefficients()) for el in temp])
        
        #make it so that the first element is positive 
        # (helps in determining the plucker coordinates of
        # of the base L_set E_1, G_4, E_2, G_3, E_3)
        if temp[temp.nonzero_positions()[0]] < 0:
            self.components = -temp
        else:
            self.components = temp
    
    #returns a new Point where each component has the new substitution
    def subs(self, sost):
        components = vector([comp.subs(sost)for comp in self.components])
        #do not want to work in fraction fields
        return Point([self._poly_ring(el) for el in components * components.denominator()])
    
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
        return self._poly_ring
