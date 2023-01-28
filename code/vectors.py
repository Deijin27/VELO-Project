"""

Solution to Worksheet 1
with additional functionality added to FourVector.

"""

from numpy import arctanh, sinh, arccos, cos, sin, cosh


class Vector:
    """
    The 'Vector' class represents a vector of length 4.
    """

    def __str__(self):
        """
        Return a string to print of this vector.
        """
        return "%r %r %r %r" % tuple(self.v)

    def __repr__(self):
        """
        Return the representation of this vector.
        """
        return "%s(%r, %r, %r, %r)" % tuple([self.__class__.__name__] + self.v)

    def __init__(self, v0, v1, v2, v3):
        """
        Initialise the class with its four components.
        """
        self.v = [v0, v1, v2, v3]

    def __getitem__(self, i):
        """
        Return the 'i' component of the vector. Note, this can be a
        slice.
        """
        return self.v[i]

    def __setitem__(self, i, s):
        """
        Set the 'i' component of the vector with the scalar
        's'. Alternatively, 'i' can be a slice and 's' a sub-vector.
        """
        self.v[i] = s

    def __len__(self):
        """
        Return the length of the vector, e.g. 4.
        """
        return 4

    def __pos__(self):
        """
        Return a copy of this vector with the '+' operator applied to
        each element.
        """
        return self.__class__(+self[0], +self[1], +self[2], +self[3])

    def __neg__(self):
        """
        Return a copy of this vector with the '-' operator applied to
        each element.
        """
        return self.__class__(-self[0], -self[1], -self[2], -self[3])

    def __iadd__(self, v):
        """
        Augmented assignment '+=' for adding a vector 'v' to this vector.
        """
        for i in range(0, 4): self[i] += v[i]
        return self

    def __isub__(self, v):
        """
        Augmented assignment '-=' for subtracting a vector 'v' from
        this vector.
        """
        for i in range(0, 4): self[i] -= v[i]
        return self

    def __add__(self, v):
        """
        Return the addition of this vector with the vector 'v'.
        """
        u = +self
        u += v
        return u

    def __sub__(self, v):
        """
        Return the subtraction of this vector with the vector 'v'.
        """
        u = +self
        u -= v
        return u

    def __invert__(self):
        """
        Return the complex transpose of this vector.
        """
        v = +self
        for i in range(0, 4):
            try: v[i] = self[i].conjugate()
            except: v[i] = self[i]
        return v

    def __imul__(self, x):
        """
        Augmented assignment '*=' for multiplying this vector with a
        vector, matrix, or scalar 'x'.
        """
        # The vector case.
        try:
            x[0]
            self = sum(self[i]*x[i] for i in range(0, 4))
        except:
            # The matrix case.
            try:
                x[0, 0]
                u = +self
                for i in range(0, 4):
                    self[i] = sum([u[j]*x[j, i] for j in range(0, 4)])
            # The scalar case.
            except:
                for i in range(0, 4): self[i] *= x
        return self

    def __mul__(self, x):
        """
        Return the multiplication of this vector with a vector, matrix, or
        scalar 'x'.
        """
        u = +self
        u *= x
        return u

    def __rmul__(self, x):
        """
        Return the multiplication of a vector, matrix, or scalar 'x'
        with this vector. The operation x*v where x is a vector or
        matrix is not used.
        """
        return self*x

    def __itruediv__(self, s):
        """
        Augmented assignment '/=' for dividing this vector with a
        scalar 's'.
        """
        for i in range(0, 4): self[i] /= s
        return self

    def __truediv__(self, s):
        """
        Return the division of this vector by a scalar 's'. The
        reflected operator, 's/v', cannot be implemented since this is
        not a defined mathematical operation.
        """
        u = +self
        u /= s
        return u

    def __ipow__(self, i):
        """
        Augmented assignment '**=' for raising this vector to the
        integer power 'i'. For even 'i' this is a scalar and odd 'i' a
        vector.
        """
        if i < 0: raise ValueError('power must be positive')
        u = (~self)*self
        if i == 2: self = u
        # The even case.
        elif i % 2 == 0: self = u**int(i/2)
        # The odd case.
        elif i % 2 == 1: self *= u**int((i - 1)/2)
        return self

    def __pow__(self, i):
        """
        Return this vector raised to the integer power 'i'. For even
        'i' this is a scalar and odd 'i' a vector.
        """
        u = +self
        u **= i
        return u

    def __abs__(self):
        """
        Return the norm of the vector.
        """
        from math import sqrt
        return sqrt(self**2)


class FourVector(Vector):
    """
    The 'FourVector' class represents a physics four-vector.
    """

    def __init__(self, *args):
        """
        Constructs the four-vector from its base Vector class.
        """
        Vector.__init__(self, *args)

    def __invert__(self):
        """
        Return the lowered index of this vector, e.g. p_mu.
        """
        v = +self
        for i in range(1, 4): v[i] = -v[i]
        return v

    @classmethod
    def from_pseudo(cls, pseudorapidity, angle_of_azimuth, momentum_magnitude, energy=None):
        """ Generate cartesian four vector from pseudorapidity, azimuth etc."""
        eta = pseudorapidity
        phi = angle_of_azimuth
        p = momentum_magnitude

        pt = p / cosh(eta)

        px = pt * cos(phi)
        py = pt * sin(phi)
        pz = pt * sinh(eta)

        return cls(energy, px, py, pz)

    def momentum_magnitude(self):
        """ Return the magnitude of the momentum. """
        E, px, py, pz = self.v
        return (px**2 + py**2 + pz**2)**0.5

    def pseudorapidity(self):
        """ Return the pseudorapidity. """
        pz = self[3]
        p = self.momentum_magnitude()
        return -arctanh(pz/p)

    def transverse_momentum(self):
        """ Return the transverse momentum, i.e. perpendicular to z-axis."""
        pz = self[3]
        eta = self.pseudorapidity()
        return pz / sinh(eta)

    def angle_of_azimuth(self):
        """ Return the angle of azimuth. """
        px = self[1]
        pt = self.transverse_momentum()
        return arccos(px / pt)
