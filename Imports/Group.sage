class Group:
    """
    Class that represent a matrix group
    """
    def __init__(self, elements):
        self.elements = elements
        self.order = len(self.elements)
        self.base_ring = self.elements[0].base_ring()
        self.identity = matrix.identity(4).change_ring(self.base_ring)

    def __getitem__(self, item):
        return self.elements.__getitem__(item)

    def __iter__(self):
        return (el for el in self.elements)

    def __str__(self):
        if self.is_abelian():
            return 'Abelian group of order ' + str(self.order)
        else:
            return 'Non abelian group of order ' + str(self.order)

    def __repr__(self):
        return self.__str__()

    def __repr_html__(self):
        return self.__str__()

    def is_in_group(self, g):
        for mat in self.elements:
            if are_matrices_equal(g, mat):
                return True
        return False

    def is_group(self):
        return False not in [self.is_in_group(mat1 * mat2) for mat1 in self.elements for mat2 in self.elements]

    def is_abelian(self):
        for mat1 in self.elements:
            for mat2 in self.elements:
                if not are_matrices_equal(mat1 * mat2, mat2 * mat1):
                    return False
        return True

    def get_order(self):
        return self.order

    def get_divisors_of_order(self):
        divisors = []
        for i in range(1, self.order // 2 + 1):
            if self.order % i == 0:
                divisors.append(i)
        divisors.append(self.order)
        return divisors

    def get_order_of_elements(self):
        orders = []
        for el in self.elements:
            for n in self.get_divisors_of_order():
                if are_matrices_equal(el ^ n, self.identity):
                    orders.append(n)
                    break
        return orders

    def apply_to(self, v):
        return [v * el for el in self.elements]
