from decimal import Decimal, getcontext
from copy import deepcopy

from vector import Vector
from plane import Plane

getcontext().prec = 30


class LinearSystem(object):

    ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG = 'All planes in the system should live in the same dimension'
    NO_SOLUTIONS_MSG = 'No solutions'
    INF_SOLUTIONS_MSG = 'Infinitely many solutions'

    def __init__(self, planes):
        try:
            d = planes[0].dimension
            for p in planes:
                assert p.dimension == d

            self.planes = planes
            self.dimension = d

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)


    def swap_rows(self, row1, row2):
        tmp = self.planes[row1]
        self.planes[row1] = self.planes[row2]
        self.planes[row2] = tmp


    def multiply_coefficient_and_row(self, coefficient, row):
        p = self.planes[row]
        n = p.normal_vector
        new_coordinates = tuple(x * coefficient for x in n.coordinates)
        new_constant = p.constant_term * coefficient
        self.planes[row] = Plane(Vector(new_coordinates), new_constant)


    def add_multiple_times_row_to_row(self, coefficient, row_to_add, row_to_be_added_to):
        plane_to_add = self.planes[row_to_add]
        plane_to_add_to = self.planes[row_to_be_added_to]
        n = plane_to_add.normal_vector
        new_n = n.times_scalar(coefficient)
        new_n = new_n.plus(plane_to_add_to.normal_vector)
        constant = plane_to_add.constant_term * coefficient
        constant = constant + plane_to_add_to.constant_term
        self.planes[row_to_be_added_to] = Plane(new_n, constant)



    def indices_of_first_nonzero_terms_in_each_row(self):
        num_equations = len(self)
        num_variables = self.dimension

        indices = [-1] * num_equations

        for i,p in enumerate(self.planes):
            try:
                indices[i] = p.first_nonzero_index(p.normal_vector)
            except Exception as e:
                if str(e) == Plane.NO_NONZERO_ELTS_FOUND_MSG:
                    continue
                else:
                    raise e

        return indices


    def compute_triangular_form(self):
        system = deepcopy(self)
        non_zero_indices = system.indices_of_first_nonzero_terms_in_each_row()
        i = 0
        while (i < len(system.planes)):
            swapped = False
            j = 1
            while (j < len(system.planes) - i):
                if system.planes[i] == system.planes[i+j]:
                    system.planes[i+j] = Plane()
                elif non_zero_indices[i] == non_zero_indices[i+j]:
                    index_of_coefficient = non_zero_indices[i]
                    coefficient = -1 * system.planes[i+j].normal_vector[index_of_coefficient] / system[i].normal_vector[index_of_coefficient]
                    system.add_multiple_times_row_to_row(coefficient, i, i+j)
                    non_zero_indices = system.indices_of_first_nonzero_terms_in_each_row()
                elif non_zero_indices[i] > non_zero_indices[i+j]:
                    system.swap_rows(i, i+j)
                    non_zero_indices = system.indices_of_first_nonzero_terms_in_each_row()
                    swapped = True
                    break
                j += 1
            if not swapped:
                i += 1
        return system


    def compute_rref(self):
        tf = self.compute_triangular_form()
        non_zero_indices = tf.indices_of_first_nonzero_terms_in_each_row()
        i = len(tf) - 1
        while (i >= 0):
            nv = tf[i].normal_vector
            if not nv.is_zero():
                indicie_of_coefficient = non_zero_indices[i]
                coefficient = nv[indicie_of_coefficient]
                if coefficient != 1:
                    tf.multiply_coefficient_and_row(1 / coefficient, i)
                    non_zero_indices = tf.indices_of_first_nonzero_terms_in_each_row()
                j = 1
                while (j <=  i):
                    coefficient = tf.planes[i-j].normal_vector[indicie_of_coefficient]
                    if coefficient != 0:
                        tf.add_multiple_times_row_to_row(-1 * coefficient, i, i-j)
                        non_zero_indices = tf.indices_of_first_nonzero_terms_in_each_row()
                    j += 1
            i -= 1
        return tf

    def do_gaussian_elimination_and_extract_solution(self):
        rref = self.compute_rref()

        non_zero_coefficients = rref.indices_of_first_nonzero_terms_in_each_row()
        num_pivots = sum([1 if x >=0 else 0 for x in non_zero_coefficients])
        solution_coordinates = [''] * self.dimension
        print(solution_coordinates)
        for i, p in enumerate(rref):
            if p.normal_vector.is_zero() and p.constant_term != 0:
                raise Exception(self.NO_SOLUTIONS_MSG)
            if non_zero_coefficients[i] >= 0:
                indicie_of_variable = non_zero_coefficients[i]
                solution_coordinates[indicie_of_variable] = p.constant_term
        if num_pivots < rref.dimension:
            raise Exception(self.INF_SOLUTIONS_MSG)
        return Vector(solution_coordinates)

    def __len__(self):
        return len(self.planes)


    def __getitem__(self, i):
        return self.planes[i]


    def __setitem__(self, i, x):
        try:
            assert x.dimension == self.dimension
            self.planes[i] = x

        except AssertionError:
            raise Exception(self.ALL_PLANES_MUST_BE_IN_SAME_DIM_MSG)


    def __str__(self):
        ret = 'Linear System:\n'
        temp = ['Equation {}: {}'.format(i+1,p) for i,p in enumerate(self.planes)]
        ret += '\n'.join(temp)
        return ret


class MyDecimal(Decimal):
    def is_near_zero(self, eps=1e-10):
        return abs(self) < eps


p0 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
p1 = Plane(normal_vector=Vector(['0','1','0']), constant_term='2')
p2 = Plane(normal_vector=Vector(['1','1','-1']), constant_term='3')
p3 = Plane(normal_vector=Vector(['1','0','-2']), constant_term='2')

s = LinearSystem([p0,p1,p2,p3])

print(s.indices_of_first_nonzero_terms_in_each_row())
print('{},{},{},{}'.format(s[0],s[1],s[2],s[3]))
print(len(s))
print(s)

s[0] = p1
print(s)

print(MyDecimal('1e-9').is_near_zero())
print(MyDecimal('1e-11').is_near_zero())

#19.Quiz Coding Row Operations
p0 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
p1 = Plane(normal_vector=Vector(['0','1','0']), constant_term='2')
p2 = Plane(normal_vector=Vector(['1','1','-1']), constant_term='3')
p3 = Plane(normal_vector=Vector(['1','0','-2']), constant_term='2')

s = LinearSystem([p0,p1,p2,p3])
s.swap_rows(0,1)
if not (s[0] == p1 and s[1] == p0 and s[2] == p2 and s[3] == p3):
    print('test case 1 failed')

s.swap_rows(1,3)
if not (s[0] == p1 and s[1] == p3 and s[2] == p2 and s[3] == p0):
    print('test case 2 failed')

s.swap_rows(3,1)
if not (s[0] == p1 and s[1] == p0 and s[2] == p2 and s[3] == p3):
    print('test case 3 failed')

s.multiply_coefficient_and_row(1,0)
if not (s[0] == p1 and s[1] == p0 and s[2] == p2 and s[3] == p3):
    print('test case 4 failed')

s.multiply_coefficient_and_row(-1,2)
if not (s[0] == p1 and
        s[1] == p0 and
        s[2] == Plane(normal_vector=Vector(['-1','-1','1']), constant_term='-3') and
        s[3] == p3):
    print('test case 5 failed')

s.multiply_coefficient_and_row(10,1)
if not (s[0] == p1 and
        s[1] == Plane(normal_vector=Vector(['10','10','10']), constant_term='10') and
        s[2] == Plane(normal_vector=Vector(['-1','-1','1']), constant_term='-3') and
        s[3] == p3):
    print('test case 6 failed')

s.add_multiple_times_row_to_row(0,0,1)
if not (s[0] == p1 and
        s[1] == Plane(normal_vector=Vector(['10','10','10']), constant_term='10') and
        s[2] == Plane(normal_vector=Vector(['-1','-1','1']), constant_term='-3') and
        s[3] == p3):
    print('test case 7 failed')

s.add_multiple_times_row_to_row(1,0,1)
if not (s[0] == p1 and
        s[1] == Plane(normal_vector=Vector(['10','11','10']), constant_term='12') and
        s[2] == Plane(normal_vector=Vector(['-1','-1','1']), constant_term='-3') and
        s[3] == p3):
    print('test case 8 failed')

s.add_multiple_times_row_to_row(-1,1,0)
if not (s[0] == Plane(normal_vector=Vector(['-10','-10','-10']), constant_term='-10') and
        s[1] == Plane(normal_vector=Vector(['10','11','10']), constant_term='12') and
        s[2] == Plane(normal_vector=Vector(['-1','-1','1']), constant_term='-3') and
        s[3] == p3):
    print('test case 9 failed')

# 20 Quiz: Coding Triagular Form
p1 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['0','1','1']), constant_term='2')
s = LinearSystem([p1,p2])
t = s.compute_triangular_form()
if not (t[0] == p1 and
        t[1] == p2):
    print('Triangular Form: test case 1 failed')
#print('TF Test1:\n{}'.format(t))

p1 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['1','1','1']), constant_term='2')
s = LinearSystem([p1,p2])
t = s.compute_triangular_form()
if not (t[0] == p1 and
        t[1] == Plane(constant_term='1')):
    print('Triangular Form: test case 2 failed')
#print('TF Test2:\n{}'.format(t))

p1 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['0','1','0']), constant_term='2')
p3 = Plane(normal_vector=Vector(['1','1','-1']), constant_term='3')
p4 = Plane(normal_vector=Vector(['1','0','-2']), constant_term='2')
s = LinearSystem([p1,p2,p3,p4])
t = s.compute_triangular_form()
if not (t[0] == p1 and
        t[1] == p2 and
        t[2] == Plane(normal_vector=Vector(['0','0','-2']), constant_term='2') and
        t[3] == Plane()):
    print('Triangular Form: test case 3 failed')
#print('TF Test3:\n{}'.format(t))

p1 = Plane(normal_vector=Vector(['0','1','1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['1','-1','1']), constant_term='2')
p3 = Plane(normal_vector=Vector(['1','2','-5']), constant_term='3')
s = LinearSystem([p1,p2,p3])
t = s.compute_triangular_form()
if not (t[0] == Plane(normal_vector=Vector(['1','-1','1']), constant_term='2') and
        t[1] == Plane(normal_vector=Vector(['0','1','1']), constant_term='1') and
        t[2] == Plane(normal_vector=Vector(['0','0','-9']), constant_term='-2')):
    print('Triangular Form: test case 4 failed')
#print('TF Test4:\n{}'.format(t))

# 21 Quiz: Coding RREF
p1 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['0','1','1']), constant_term='2')
s = LinearSystem([p1,p2])
r = s.compute_rref()
if not (r[0] == Plane(normal_vector=Vector(['1','0','0']), constant_term='-1') and
        r[1] == p2):
    print('RREF: test case 1 failed')

p1 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['1','1','1']), constant_term='2')
s = LinearSystem([p1,p2])
r = s.compute_rref()
if not (r[0] == p1 and
        r[1] == Plane(constant_term='1')):
    print('RREF: test case 2 failed')

p1 = Plane(normal_vector=Vector(['1','1','1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['0','1','0']), constant_term='2')
p3 = Plane(normal_vector=Vector(['1','1','-1']), constant_term='3')
p4 = Plane(normal_vector=Vector(['1','0','-2']), constant_term='2')
s = LinearSystem([p1,p2,p3,p4])
r = s.compute_rref()
if not (r[0] == Plane(normal_vector=Vector(['1','0','0']), constant_term='0') and
        r[1] == p2 and
        r[2] == Plane(normal_vector=Vector(['0','0','-2']), constant_term='2') and
        r[3] == Plane()):
    print('RREF: test case 3 failed')

p1 = Plane(normal_vector=Vector(['0','1','1']), constant_term='1')
p2 = Plane(normal_vector=Vector(['1','-1','1']), constant_term='2')
p3 = Plane(normal_vector=Vector(['1','2','-5']), constant_term='3')
s = LinearSystem([p1,p2,p3])
r = s.compute_rref()
if not (r[0] == Plane(normal_vector=Vector(['1','0','0']), constant_term=Decimal('23')/Decimal('9')) and
        r[1] == Plane(normal_vector=Vector(['0','1','0']), constant_term=Decimal('7')/Decimal('9')) and
        r[2] == Plane(normal_vector=Vector(['0','0','1']), constant_term=Decimal('2')/Decimal('9'))):
    print('RREF: test case 4 failed')

p1 = Plane(Vector(['5.862', '1.178', '-10.366']), '-8.15')
p2 = Plane(Vector(['-2.931', '-0.589', '5.183']), '-4.075')
s = LinearSystem([p1,p2])
r = s.compute_rref()
print('***first test:***')
#print('TF:\n{}'.format(s.compute_triangular_form()))
try:
    solution = r.do_gaussian_elimination_and_extract_solution()
    print('RREF:\n{}'.format(Plane(solution)))
except Exception as e:
    if str(e) == LinearSystem.NO_SOLUTIONS_MSG or str(e) == LinearSystem.INF_SOLUTIONS_MSG:
        print(str(e))
    else: raise e

p1 = Plane(Vector([8.631, 5.112, -1.816]), -5.113)
p2 = Plane(Vector([4.315, 11.132, -5.27]), -6.775)
p3 = Plane(Vector([-2.158, 3.01, -1.727]), -0.831)
s = LinearSystem([p1,p2,p3])
r = s.compute_rref()
print('***second test:***')
#print('TF:\n{}'.format(s.compute_triangular_form()))
try:
    solution = r.do_gaussian_elimination_and_extract_solution()
    print('RREF:\n{}'.format(Plane(solution)))
except Exception as e:
    if str(e) == LinearSystem.NO_SOLUTIONS_MSG or str(e) == LinearSystem.INF_SOLUTIONS_MSG:
        print(e)
    else: raise e

p1 = Plane(Vector([5.262, 2.739, -9.878]), -3.441)
p2 = Plane(Vector([5.111, 6.358, 7.638]), -2.152)
p3 = Plane(Vector([2.016, -9.924, -1.367]), -9.278)
p4 = Plane(Vector([2.167, -13.543, -18.883]), -10.567)
s = LinearSystem([p1,p2,p3,p4])
r = s.compute_rref()
print('***third test:****')
#print('TF:\n{}'.format(s.compute_triangular_form()))
try:
    solution = r.do_gaussian_elimination_and_extract_solution()
    print('RREF:\n{}'.format(Plane(solution)))
except Exception as e:
    if str(e) == LinearSystem.NO_SOLUTIONS_MSG or str(e) == LinearSystem.INF_SOLUTIONS_MSG:
        print(e)
    else: raise e

