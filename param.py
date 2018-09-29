class Parametrization:

    BASEPT_AND_DIR_VECTORS_MUST_BE_IN_SAME_DIM = (
        'The basepoint and direction vectors should all live in the same '
        'dimension')

    def __init__(self, basepoint, direction_vectors):

        self.basepoint = basepoint
        self.direction_vectors = direction_vectors
        self.dimension = self.basepoint.dimension

        try:
            for v in direction_vectors:
                assert v.dimension == self.dimension

        except AssertionError:
            raise Exception(self.BASEPT_AND_DIR_VECTORS_MUST_BE_IN_SAME_DIM)


    def __str__(self):

        system = []
        for i in range(self.dimension):
            my_string = 'x_{} = {}'.format(i+1, round(self.basepoint[i],3))
            for j, v in enumerate(self.direction_vectors):
                n = round(v[i], 3)
                my_string += ' + {}t_{}'.format(n, j+1)  
            system.append(my_string)
        my_string = 'System\n'
        my_string += '\n'.join(system)
        return my_string

