class Line:
    def __init__(self, planes, points=None, plucker=None):
        self.P = planes[0].parent()
        self.planes = planes
        if points is None or plucker is None:
            self.points = self.get_two_points_on_line()
            self.plucker = Point(vector(matrix([p for p in self.points]).minors(2)))
        else:
            self.points = points
            self.plucker = plucker

    def __str__(self):
        return self.plucker.__str__() + ' ' + self.planes.__str__()

    def __repr__(self):
        return self.plucker.__repr__() + ' ' + self.planes.__repr__()

    def __eq__(self, other_line):
        if isinstance(other_line, Line):
            return self.plucker == other_line.plucker
        return False

    def subs(self, sost):
        plucker = vector([pl.subs(sost) for pl in self.plucker])
        plucker = Point([self.P(el) for el in plucker * plucker.denominator()])
        planes = [self.P(pl.subs(sost).numerator()) for pl in self.planes]
        if matrix([plane_coefficients(plane) for plane in planes]).minors(2) == [0 for _ in range(6)]:
            planes = get_planes(plucker)
        points = [pl.subs(sost) for pl in self.points]
        if points[0] == points[1]:
            return Line(planes)
        return Line(planes, points, plucker)

    def get_two_points_on_line(self):
        vr = self.P.gens()
        M = matrix([plane_coefficients(plane) for plane in self.planes])
        minors = M.minors(2)
        index = next(i[0] for i in enumerate(minors) if i[1] != 0)
        possible_columns = {0: (0, 1), 1: (0, 2), 2: (0, 3), 3: (1, 2), 4: (1, 3), 5: (2, 3)}
        columns = possible_columns.get(index)
        other_columns = tuple({0, 1, 2, 3} - set(columns))
        sol = solve_linear_system(self.planes, [vr[columns[0]], vr[columns[1]]],
                                  [vr[other_columns[0]], vr[other_columns[1]]])
        sost = {vr[columns[0]]: sol[0], vr[columns[1]]: sol[1], vr[other_columns[0]]: sol[2],
                vr[other_columns[1]]: sol[3]}
        point1_sost = {vr[other_columns[0]]: 0, vr[other_columns[1]]: 1}
        point2_sost = {vr[other_columns[0]]: 1, vr[other_columns[1]]: 0}

        not_ordered_vars = [vr[columns[0]], vr[columns[1]], vr[other_columns[0]], vr[other_columns[1]]]
        ordered_vars = []
        for i in range(4):
            for el in not_ordered_vars:
                if el == vr[i]:
                    ordered_vars.append(el)
                    continue
        point1 = vector([el.subs(sost).subs(point1_sost) for el in ordered_vars])
        point2 = vector([el.subs(sost).subs(point2_sost) for el in ordered_vars])
        return Point(point1), Point(point2)

    def are_incident(self, other_line):
        d = [self.plucker[i] * other_line.plucker[5 - i] for i in range(6)]
        return d[0] - d[1] + d[2] + d[3] - d[4] + d[5] == 0

    def intersection_point(self, other_line):
        vrs = self.P.gens()
        plane1, plane2 = self.planes
        plane3, plane4 = other_line.planes
        if matrix([plane_coefficients(plane) for plane in [plane1, plane2, plane3]]).rank() == 3:
            sol = solve_linear_system([plane1, plane2, plane3], vrs[0:3], [vrs[3]])
            return Point(vector([s.subs({vrs[3]: 1}) for s in sol]))
        else:
            sol = solve_linear_system([plane1, plane2, plane4], vrs[0:3], [vrs[3]])
            return Point(vector([s.subs({vrs[3]: 1}) for s in sol]))

    def parent(self):
        return self.P

    def apply_proj(self, proj):
        sost = change_coord(proj)
        new_planes = [plane.subs(sost) for plane in self.planes]
        return Line(new_planes)
