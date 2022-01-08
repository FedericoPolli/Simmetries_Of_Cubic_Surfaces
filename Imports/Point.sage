class Point:   
    """
    Class that represents a point in any projective space, 
    i.e. a vector defined up to multiplication with a nonzero scalar.
    """

    def __init__(self, components):
        if list(set(components)) == [0]:
            raise ValueError(str(components)+" is not a projective point")
            
        self._poly_ring = components[0].parent()
        mcd = gcd(components)
        # divide each component by gcd with exact division
        temp = vector([comp//mcd for comp in components])
        
        # calculate the gcd of the coefficients of each component
        #  as gcd of polynomials returns 1 as coefficient
        try:
            temp = temp/gcd([gcd(el.coefficients()) for el in temp])
        except Exception:
            pass
        
        # make it so that the first element is positive
        # (helps in determining the plucker coordinates of
        # of the base L_set E_1, G_4, E_2, G_3, E_3)
        if temp[temp.nonzero_positions()[0]] < 0:
            self._components = -temp
        else:
            self._components = temp
    
    def reduce(self, ideal):
        return Point([ideal.reduce(el) for el in self._components])

    # returns a new Point where each component has the new substitution
    def subs(self, substitution):
        components = vector([comp.subs(substitution) for comp in self._components])
        # do not want to work in fraction fields
        return Point([self._poly_ring(el) for el in components * components.denominator()])
    
    def __str__(self):
        return self._components.__str__()
        
    def __repr__(self):
        return self._components.__repr__()
        
    def __repr_html__(self):
        return self._components.__repr_html__()
     
    def __iter__(self):
        return (el for el in self._components)
        
    def __len__(self):
        return len(self._components)
        
    def __getitem__(self, item):
        return self._components.__getitem__(item)
    
    def __eq__(self, other_point):
        if isinstance(other_point, Point):
            return are_vectors_proportional(self._components, other_point._components)
        return False
    
    #uses mul, rmul of underlying vector
    def __mul__(self, other):
        return Point(self._components * other)
        
    def __rmul__(self, other):
        return Point(other * self._components)
    
    def __div__(self, other):
        return Point(self._components / other)
    
    def __add__(self, other):
        return Point(self._components + other.get_components())
    
    def __sub__(self, other):
        return Point(self._components - other.get_components())
        
    def parent(self):
        return self._poly_ring

    def get_components(self):
        return self._components
