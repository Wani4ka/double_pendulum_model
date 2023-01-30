from math import sqrt


class Matrix:
    def __init__(self, a):
        self.vectors = [Vector([j for j in i]) for i in a]

    def __len__(self):
        return len(self.vectors)

    def __sub__(self, other):
        if isinstance(other, Matrix) and len(other.vectors) == len(self.vectors):
            return Matrix([self.vectors[i] - other.vectors[i] for i in range(len(self.vectors))])

    def __str__(self):
        return '\n'.join(map(str, self.vectors))

    def __mul__(self, other):
        if not isinstance(other, Matrix):
            return Matrix([v * other for v in self.vectors])
        return None

    def __rmul__(self, other):
        return Matrix([v * other for v in self.vectors])

    def __truediv__(self, other):
        if not isinstance(other, Vector):
            return None
        result = self.vectors[:]
        for i in range(len(result)):
            for j in range(len(result[i])):
                result[i][j] /= other[i]
        return Matrix(result)

    def column(self, i):
        return Vector([v[i] for v in self.vectors])

    def __getitem__(self, item):
        return self.vectors[item]

    def __add__(self, other):
        result = self.vectors[:]
        for i in range(len(result)):
            for j in range(len(result[i])):
                result[i][j] += other.vectors[i][j]
        return Matrix(result)

    def append(self, vector):
        self.vectors.append(vector)


class Vector:
    def __init__(self, els):
        self.components = list(els)[:]

    def __add__(self, other):
        comps = self.components[:]
        if not hasattr(other, '__iter__'):
            for i in range(len(comps)):
                comps[i] += other
            return Vector(comps)

        if len(comps) != len(other):
            raise Exception("Vectors must be in one dimension")
        for i in range(len(comps)):
            comps[i] += other[i]
        return Vector(comps)

    def __iter__(self):
        return iter(self.components)

    def __sub__(self, other):
        comps = self.components[:]
        if len(comps) != len(other):
            raise Exception("Vectors must be in one dimension")
        for i in range(len(comps)):
            comps[i] -= other.components[i]
        return Vector(comps)

    def __mul__(self, other):
        if isinstance(other, Vector):
            return sum([self.components[i] * other.components[i] for i in range(len(self.components))])
        else:
            return Vector([comp * other for comp in self.components])

    def __rmul__(self, other):
        return self.__mul__(other)

    def __truediv__(self, other):
        if not isinstance(other, Vector):
            return self.__mul__(1 / other)
        result = self.components[:]
        for i in range(len(result)):
            if other.components[i] != 0:
                result[i] /= other.components[i]
            else:
                result[i] = -(10**10)
        return Vector(result)

    def __len__(self):
        return len(self.components)

    def __str__(self):
        return str(self.components)

    def append(self, el):
        self.components.append(el)

    def __getitem__(self, item):
        return self.components[item]

    def __setitem__(self, item, value):
        self.components[item] = value

	# Евклидова норма
    def norm(self):
        return sqrt(sum([x ** 2 for x in self.components]))

    def abs(self):
        return Vector([abs(comp) for comp in self.components])


def arrange(start, stop, step):
    current = start + step
    result = []
    while current <= stop:
        result.append(current)
        current += step
    return result
