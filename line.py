from decimal import Decimal, getcontext

from vector import Vector

getcontext().prec = 30

class Line:

    NO_NONZERO_ELTS_FOUND_MSG = 'No nonzero elements found'
    NO_INTERSECTION_PARALLEL_LINES_MSG = 'No interesection. Lines are parallel'
    NO_INTERSECTION_EQUAL_LINES_MSG = 'No intersection. Lines are equal'


    def __init__(self, normal_vector=None, constant_term=None):
        self.dimension = 2

        if not normal_vector:
            all_zeros = ['0']*self.dimension
            normal_vector = Vector(all_zeros)
        self.normal_vector = normal_vector

        if not constant_term:
            constant_term = Decimal('0')
        self.constant_term = Decimal(constant_term)

        self.set_basepoint()

    def set_basepoint(self):
        a = self.normal_vector.coordinates[0]
        b = self.normal_vector.coordinates[1]

        if a != 0:
            self.basepoint = [self.constant_term / a, 0]
        elif b != 0:
            self.basepoint = [0, self.constant_term / b]
        else:
            raise Exception(self.NO_NONZERO_ELTS_FOUND_MSG)

    def parallel_to(self, l):
        return self.normal_vector.is_parallel(l.normal_vector)
    
    def __eq__(self, l):
        bp1 = Vector(self.basepoint)
        bp2 = Vector(l.basepoint)
        v = bp1.minus(bp2)
        return v.is_orthogonal(self.normal_vector) and v.is_orthogonal(l.normal_vector)

    def intersection_with(self, l):
        try:
            if self == l:
                raise Exception(self.NO_INTERSECTION_EQUAL_LINES_MSG)
            elif self.parallel_to(l):
                raise Exception(self.NO_INTERSECTION_PARALLEL_LINES_MSG)
            else:
                a = self.normal_vector.coordinates[0]
                b = self.normal_vector.coordinates[1]
                c = l.normal_vector.coordinates[0]
                d = l.normal_vector.coordinates[1]
                k1 = self.constant_term
                k2 = l.constant_term
                x = (d*k1 - b*k2) / (a*d - b*c)
                y = (-1*c*k1 + a*k2) / (a*d - b*c)
            return [x, y]            

        except Exception as e:
            if str(e) == self.NO_INTERSECTION_EQUAL_LINES_MSG:
                print(str(e))
            elif str(e) == self.NO_INTERSECTION_PARALLEL_LINES_MSG:
                print(str(e))
            else:
                raise e


l1 = Line(Vector([4.046, 2.836]), 1.21)
l2 = Line(Vector([10.115, 7.09]), 3.025)
print(l1.intersection_with(l2))

l1 = Line(Vector([7.204, 3.182]), 8.68)
l2 = Line(Vector([8.172, 4.114]), 9.883)
print(l1.intersection_with(l2))

l1 = Line(Vector([1.182, 5.562]), 6.744)
l2 = Line(Vector([1.773, 8.343]), 9.525)
print(l1.intersection_with(l2))

