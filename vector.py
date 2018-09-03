import math
from decimal import Decimal, getcontext

getcontext().prec = 30

class Vector():

    CANNOT_NORMALIZE_ZERO_VECTOR_MSG = 'Cannot normalize the zero vector'
    NO_UNIQUE_PARALLEL_COMPONENT_MSG = 'No unique parallel component'

    def __init__(self, coordinates):
        try:
            if not coordinates:
                raise ValueError
            self.coordinates = tuple(coordinates)
            self.dimension = len(coordinates)

        except ValueError:
            raise ValueError('The coordinates must be nonempty')

        except TypeError:
            raise TypeError('The coordinates must be an iterable')


    def __str__(self):
        return 'Vector: {}'.format(self.coordinates)


    def __eq__(self, v):
        return self.coordinates == v.coordinates
  
    def plus(self, v):
        newCoordinates = []
        for x, y in zip(self.coordinates, v.coordinates):
            newCoordinates.append(x + y)

        return Vector(newCoordinates) 

    def minus(self, v):
        newCoordinates = []
        for x, y in zip(self.coordinates, v.coordinates):
            newCoordinates.append(x - y)

        return Vector(newCoordinates) 

    def times_scalar(self, s):
        newCoordinates = []
        for x in self.coordinates:
            newCoordinates.append(x * Decimal(s))

        return Vector(newCoordinates)

    def magnitude(self):
        coordinates_squared = [x**2 for x in self.coordinates]
        return Decimal(math.sqrt(sum(coordinates_squared)))

    def normalized(self):
        try:
            return self.times_scalar(Decimal('1.0')/self.magnitude())

        except ZeroDivisionError:
            raise Exception('Cannot normalize the zero vector')

    def dot(self, v):
        return sum([x * y for x, y in zip(self.coordinates, v.coordinates)])

    def angle_with(self, v, in_degrees=False):
        try:
            u1 = self.normalized()
            u2 = v.normalized()
            angle_in_radians = math.acos(round(u1.dot(u2), 10))

            if in_degrees:
                return math.degrees(angle_in_radians)
            else:
                return angle_in_radians

        except Exception as e:
            if str(e) == self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG:
                raise Exception('Cannot compute an angle with the zero vector')
            else:
                raise e

    def is_parallel(self, v):
        return self.is_zero() or v.is_zero() or self.angle_with(v) == 0 or self.angle_with(v) == math.pi

    def is_zero(self, tolerance=1e-10):
        return self.magnitude() < tolerance

    def is_orthogonal(self, v, tolerance=1e-10):
        return abs(self.dot(v)) < tolerance 

    def component_parallel_to(self, basis):
        try:
            u = basis.normalized()
            weight = self.dot(u)
            return u.times_scalar(weight)

        except Exception as e:
            if str(e) == self.CANNOT_NORMALIZE_ZERO_VECTOR_MSG:
                raise Exception(self.NO_UNIQUE_PARALLEL_COMPONENT_MSG)
            else:
                raise e

    def component_orthogonal_to(self, basis):
        try:
            projection = self.component_parallel_to(basis)
            return self.minus(projection)

        except Exception as e:
            if str(e) == self.NO_UNIQUE_PARALLEL_COMPONENT_MSG:
                raise Exception(self.NO_UNIQUE_PARALLEL_COMPONENT_MSG)
            else:
                raise e

    def cross_product(self, v):
        new_coordinates = []
        x1 = self.coordinates[0]
        y1 = self.coordinates[1]
        z1 = self.coordinates[2]
        x2 = v.coordinates[0]
        y2 = v.coordinates[1]
        z2 = v.coordinates[2]

        new_coordinates.append(y1*z2 - y2*z1)
        new_coordinates.append(-1 * (x1*z2 - x2*z1))
        new_coordinates.append(x1*y2 - x2*y1)

        return Vector(new_coordinates)

    def area_of_parallelogram(self, v):
        if self.dimension == 3 and v.dimension == 3:
            return self.cross_product(v).magnitude()

    def area_of_triangle(self, v):
        if self.dimension == 3 and v.dimension == 3:
            return self.area_of_parallelogram(v) / 2

def test_vector_plus():
    u = Vector([1,2,3])
    v = Vector([1,1,1])
    return u.plus(v) == Vector([x + y for x, y in zip(u.coordinates, v.coordinates)])

def test_vector_minus():
    u = Vector([1,2,3])
    v = Vector([1,1,1])
    return u.minus(v) == Vector([x - y for x, y in zip(u.coordinates, v.coordinates)])

def test_vector_times_scalar():
    u = Vector([1,2,3])
    s = 2
    return u.times_scalar(s) == Vector([x * s for x in u.coordinates])

def test_vector_magnitude():
    u = Vector([1,2,3])
    return u.magnitude() == math.sqrt(sum([x**2 for x in u.coordinates]))

def test_vector_normalized():
    u = Vector([1,2,3])
    return u.normalized() == u.times_scalar(1.0 / u.magnitude())

def test_vector_dot():
    u = Vector([1,2,3])
    v = Vector([3,4,5])
    return u.dot(v) == sum([x * y for x, y in zip(u.coordinates, v.coordinates)])

def test_vector_angle_with():
    u = Vector([1,2,3])
    v = Vector([4,5,6])
    return u.angle_with(v) == math.acos(u.dot(v) / (u.magnitude() * v.magnitude()))

def test_vector_is_parallel():
    u = Vector([1,2,3])
    v = u.normalized()
    return u.is_parallel(v) and u.is_parallel(v.times_scalar(-1))

def test_vector_is_parallel_negative_case():
    u = Vector([1,2,3])
    v = Vector([1,2,4])
    return u.is_parallel(v) == False 

def test_vector_is_orthogonal():
    u = Vector([1,2,3])
    v = Vector([1,1,-1])
    return u.dot(v) == 0 and u.is_orthogonal(v)

