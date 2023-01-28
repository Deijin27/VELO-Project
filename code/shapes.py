"""
Shapes in Euclidean space.
"""

import numpy as np


class Point(list):
    """ A point in Euclidean space.

    Has standard functionality of a list except:

    Public Methods
    --------------
    x, y, z :
        Give the x, y, or z component of the point respectively.
    gaussian_smear :
        Randomly gaussian smear point in the x and y axes.

    Static Methods
    --------------
    scalar_separation :
        Find the scalar separation of two Point objects.
    avg :
        Get the average point of a list of pints.

    Magic Methods
    -------------
    __init__ :
        Initialise the 3 components of a new point.
    __add__ :
        Add the components of two points together.
    __iadd__ :
        Use += to add the components of two points together.
    __neg__ :
        Switch the signs of the components of the point.
    __sub__ :
        Subtract the components of the first point from the second.
    __isub__ :
        Use -= to subtract the components of the first point from the second.
    __truediv__ :
        Divide the point by a number.
    """

    def __init__(self, x, y, z=0):
        """ Create Point

        Parameters
        ----------
        x, y, z : float or int
            Components of the point position.
        """
        super().__init__([x, y, z])

    def x(self):
        """ Give the x-component of the point. """
        return self[0]

    def y(self):
        """ Give the y-component of the point. """
        return self[1]

    def z(self):
        """ Give the z-component of the point. """
        return self[2]

    def gaussian_smear(self, xsigma, ysigma):
        """ Randomly gaussian smear point in the x and y axes. """
        x, y, z = self
        x2 = np.random.normal(x, xsigma)
        y2 = np.random.normal(y, ysigma)
        return self.__class__(x2, y2, z)

    @staticmethod
    def scalar_separation(point_a, point_b, dimentions=3):
        """ Find the scalar separation of two Point objects """
        if dimentions not in (1, 2, 3):
            raise ValueError("dimentions must be 1, 2 or 3")
        sqrsum = (point_a[0] - point_b[0])**2
        if dimentions > 1:
            sqrsum += (point_a[1] - point_b[1])**2
        if dimentions == 3:
            sqrsum += (point_a[2] - point_b[2])**2
        return (sqrsum)**0.5

    def __add__(self, other):
        """ Add the components of two points together. """
        return self.__class__(*[a+b for a, b in zip(self, other)])

    def __iadd__(self, other):
        """ Use += to add the components of two points together. """
        return self + other

    def __neg__(self):
        """ Switch the signs of the components of the point. """
        return self.__class__(*[-a for a in self])

    def __sub__(self, other):
        """ Subtract the components of the first point from the second. """
        return self.__class__(*[a-b for a, b in zip(self, other)])

    def __isub__(self, other):
        """ Use -= to subtract the components of the first point from the second. """
        return self - other

    def __truediv__(self, num):
        """ Divide the point by a number. """
        return self.__class__(*[a/num for a in self])

    @staticmethod
    def avg(list_of_points):
        """ Get the average point of a list of pints. """
        tot = list_of_points[0]
        for i in list_of_points[1:]:
            tot += i
        return tot / len(list_of_points)


class Rectangle:
    """ A rectangle in a Euclidean xy plane.

    Attributes
    ----------
    xmin, xmax, ymin, ymax : float or int
        Bounding limits of the rectangle.

    Magic Methods
    -------------
    __init__ :
        Create a new rectangle instance.
    __contains__ :
        Return bool of whether the Point is within the rectangle.
    """

    def __init__(self, xmin, xmax, ymin, ymax):
        """ Create Rectangle instance.

        Parameters
        ----------
        xmin, xmax, ymin, ymax : float or int
            Bounding limits of the rectangle.
        """
        self.xmin = xmin
        self.xmax = xmax
        self.ymin = ymin
        self.ymax = ymax

    def __contains__(self, point):
        """ Return bool of whether the Point is within the rectangle. """
        if self.xmin <= point.x() <= self.xmax and self.ymin <= point.y() <= self.ymax:
            return True
        return False


class CompositeShape:
    """ A composite shape in a Euclidean xy plane.

    Attributes
    ----------
    sub_shapes : list of Rectangle and/or CompositeShape
        List of shapes that make up the composite shape.

    Magic Methods
    -------------
    __init__ :
        Create a new composite shape instance.
    __contains__ :
        Return bool of whether the Point is within the shape.
    """

    def __init__(self, *sub_shapes):
        """ Create a new composite shape from an arbitrary
        number of sub shape parameters."""
        self.sub_shapes = list(sub_shapes)

    def __contains__(self, point):
        """ Return bool of whether the Point is within the shape. """
        return any(point in shape for shape in self.sub_shapes)


class Prism:
    """ A prism in Euclidean space. It's face is in the xy plane.

    Attributes
    ----------
    xy_face : Rectangle or CompositeShape
        Shape in the xy plane that is the face of prism.
    zmin, zmax : int or float
        z-limits of the prism.

    Magic Methods
    -------------
    __init__ :
        Create new Prism instance.
    __contains__ :
        Return bool of whether the Point is within the prism.
    """

    def __init__(self, xy_face, zmin, zmax):
        """ Create a new prism instance.

        Parameters
        ----------
        xy_face : Rectangle or CompositeShape
            Shape in the xy plane that is the face of the sensor's prism shape.
        zmin, zmax : int or float
            z-limits of the prism.
        """
        self.xy_face = xy_face
        self.zmin = zmin
        self.zmax = zmax

    def __contains__(self, point):
        """ Return bool of whether the Point is within the prism. """
        in_xy = point in self.xy_face
        in_z = self.zmin <= point.z() <= self.zmax
        return all((in_xy, in_z))


if __name__ == "__main__":

    p = Point(1, 2, 3)

    r = Rectangle(0, 2, 0, 2)
    r2 = Rectangle(3, 4, 3, 4)
    cs = CompositeShape(r, r2)
    pr = Prism(cs, 2, 4)

    print("p in r:", p in r)
    print("p in cs:", p in cs)
    print("p in pr:", p in pr)
