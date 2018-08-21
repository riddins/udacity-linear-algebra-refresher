from decimal import Decimal, getcontext
from vector import Vector

getcontext().prec = 30

class Plane:

    NO_NONZERO_ELTS_FOUND_MSG = 'No nonzero elements found'
    NO_INTERSECTION_PARALLEL_PLANES_MSG = 'No interesection. Planes are parallel'
    NO_INTERSECTION_EQUAL_PLANES_MSG = 'No intersection. Planes are equal'


    def __init__(self, normal_vector=None, constant_term=None):
        self.dimension = 3

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
        c = self.normal_vector.coordinates[2]

        if a != 0:
            self.basepoint = [self.constant_term / a, 0, 0]
        elif b != 0:
            self.basepoint = [0, self.constant_term / b, 0]
        elif c != 0:
            self.basepoint = [0, 0, self.constant_term / c]
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
                raise Exception(self.NO_INTERSECTION_EQUAL_PLANES_MSG)
            elif self.parallel_to(l):
                raise Exception(self.NO_INTERSECTION_PARALLEL_PLANES_MSG)
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


p1 = Plane(Vector([-0.412, 3.806, 0.728]), -3.46)
p2 = Plane(Vector([1.03, -9.515, -1.82]), 8.65)
equal = p1==p2
parallel = p1.parallel_to(p2)
print(f'#1 is EQUAL: {equal}, and PARALLEL: {parallel}')

p1 = Plane(Vector([2.611, 5.528, 0.283]), 4.6)
p2 = Plane(Vector([7.715, 8.306, 5.342]), 3.76)
equal = p1==p2
parallel = p1.parallel_to(p2)
print(f'#2 is EQUAL: {equal}, and PARALLEL: {parallel}')

#z coordinate in p1 changed from -7.217
p1 = Plane(Vector([-7.926, 8.625, -7.212]), -7.952)
p2 = Plane(Vector([-2.642, 2.875, -2.404]), -2.443)
equal = p1==p2
parallel = p1.parallel_to(p2)
print(f'#3 is EQUAL: {equal}, and PARALLEL: {parallel}')

