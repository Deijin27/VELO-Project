"""

Vertex Locator

As part of the coding project:

"Computer Modelling of a Vertex Locator"

by Mia Boulter

---

Dependencies:

- numpy
- matplotlib (including mpl_toolkits, which should be there by default)

---

Other files that are part of the package (should be in the same folder):

- shapes.py :
    Shapes used to construct a VELOSensor.

- vectors.py :
    FourVector for Particle

- progress.py :
    Shell/Command-line progress bar used for one of the processes which takes a while.
    Note: some shells such as IDLE's one do not allow rewriting of lines, so this progress
    bar will not work. To disable there is an option in the call of the any function which uses it.

"""

import numpy as np
import matplotlib.pyplot as plt
import mpl_toolkits.mplot3d as m3d

from progress import ProgressBar
from shapes import Point, Rectangle, CompositeShape, Prism
from vectors import FourVector


class Particle:
    """ A point particle in Eulcidian space.

    Attributes
    ----------
    four_vector : FourVector
        Momentum four-vector of the particle.
    init_point : Point
        Initial point in euclidean space.
    point : Point
        Current point in euclidean space.
    """

    def __init__(self, four_vector, init_point=Point(0, 0, 0)):
        """ Create Particle

        Parameters
        ----------
        four_vector : FourVector
            Momentum four-vector of the particle.
        init_point : Point, optional
            Initial point in euclidean space, defaults to Point(0, 0, 0).
        """
        self.four_vector = four_vector
        self.init_point = init_point
        self.point = init_point

    def reset(self):
        """ Change the particle's `point` back to it's `init_point` value. """
        self.point = self.init_point

    def advance(self, z_distance):
        """ Advance particle to a z position which is the sum of it's current `point`
        z component and the given `z_distance`, then apply the corresponding changes
        it's x and y positions dependent of the proportions of it's momenta.

        Care should be taken since this allows movement backwards in time from it's `init_point`.
        """
        E, px, py, pz = self.four_vector.v

        time = z_distance / pz

        dpoint = Point(px*time, py*time, z_distance)

        self.point += dpoint


class VELOSensor(Prism):
    """ A class describing a VELO sensor.

    Instance Attributes
    -------------------
    z_position : float or int
        Position in the z-direction where the sensor is.
    xy_face : Rectangle or CompositeShape
        Shape in the xy plane that is the face of the sensor's prism shape.
    zmin, zmax : int or float
        z-limits of the sensor.
    direction : str ('left', 'right', or other)
        Whether the sensor is left or right, or some other thing.
    sid : str, optional
        Sensor identifier.

    Static Attributes
    -----------------
    std_left_shape : CompositeShape
        Standard left shape of VELO sensor.
    std_right_shape : CompositeShape
        Standard right shape of VELO sensor.
    
    Public Methods
    --------------
    draw :
        Draw the sensor on an m3d.Axes3D
    details :
        Return a string containing details of the sensor.

    Class Methods
    -------------
    left_sensor :
        Return a VELOSensor to the standard left sensor specifications.
    right_sensor :
        Return a VELOSensor to the standard right sensor specifications.

    Static Methods
    --------------
    generate_sid :
        Genearate a sensor id using the given parameters.
    """
    std_left_shape = CompositeShape(Rectangle(-30.71, -2.55, -33.51, 33.51), Rectangle(-39.96, 2.55, 2.55, 33.51))
    std_right_shape = CompositeShape(Rectangle(2.55, 30.71, -33.51, 33.51), Rectangle(-2.55, 39.96, -33.51, -2.55))

    def __init__(self, xy_face, z_position, direction, sid=None):
        """ Initialise a VELOSensor

        Parameters
        ----------
        xy_face : Rectangle or CompositeShape
            Shape in the xy plane that is the face of the sensor's prism shape.
        z_position : int or float
            Position in the z-direction where the sensor is.
        direction : str ('left', 'right', or other)
            Whether the sensor is left or right, or some other thing.
        sid : str, optional
            Sensor identifier. Generates automatically using generate_sid if not given.
        """
        self.z_position = z_position
        z_thickness = 0.2
        self.direction = direction
        self.sid = sid
        if sid is None:
            self.sid = self.generate_sid(z_position, direction)
        super().__init__(xy_face, self.z_position - z_thickness/2, self.z_position + z_thickness/2)

    @classmethod
    def left_sensor(cls, z_position):
        """ Return a VELOSensor to the standard left sensor specifications. """
        direction = "left"
        sid = cls.generate_sid(z_position, direction)
        shape = cls.std_left_shape
        return cls(shape, z_position, direction, sid)

    @classmethod
    def right_sensor(cls, z_position):
        """ Return a VELOSensor to the standard right sensor specifications. """
        direction = "right"
        shape = cls.std_right_shape
        sid = cls.generate_sid(z_position, direction)
        return cls(shape, z_position, direction, sid)

    @staticmethod
    def generate_sid(z_position, direction):
        """ Genearate a sensor id using the given parameters.

        Parameters
        ----------
        z_position : int or float
            Position in the z-direction where the sensor is.
        direction : str ('left', 'right', or other)
            Whether the sensor is left or right, or some other thing.
        """

        if direction == 'left':
            dmark = "L"
        elif direction == 'right':
            dmark = "R"
        else:
            dmark = "U"
        if z_position < 0:
            sign = "N"
        else:
            sign = "P"
        return "{}-{}{}".format(dmark, sign, str(abs(z_position)).zfill(3))

    def draw(self, ax, alpha=0.2, color=None):
        """ Draw the sensor on an m3d.Axes3D """
        if color is None:
            color = np.random.rand(3,)

        z = self.z_position
        zarr = np.array([[z, z], [z, z]])

        if self.direction == 'right':
            # bottom portion
            x = np.array([[-2.55, -2.55], [39.96, 39.96]])
            y = np.array([[-33.51, -2.55], [-33.51, -2.55]])
            ax.plot_surface(x, y, zarr, alpha=alpha, color=color)
            # top portion
            x = np.array([[2.55, 2.55], [30.71, 30.71]])
            y = np.array([[-2.55, 33.51], [-2.55, 33.51]])
            ax.plot_surface(x, y, zarr, alpha=alpha, color=color)
        elif self.direction == 'left':
            # bottom portion
            x = np.array([[-30.71, -30.71], [-2.55, -2.55]])
            y = np.array([[-33.51, 2.55], [-33.51, 2.55]])
            ax.plot_surface(x, y, zarr, alpha=alpha, color=color)
            # top portion
            x = np.array([[-39.96, -39.96], [2.55, 2.55]])
            y = np.array([[2.55, 33.51], [2.55, 33.51]])
            ax.plot_surface(x, y, zarr, alpha=alpha, color=color)
        else:
            raise ValueError("Invalid direction")

    def details(self):
        """ Return a string containing details of the sensor. """
        lines = [
            "VELOSensor:",
            f"sid = {self.sid}",
            f"z_position = {self.z_position}",
            f"z_limits = {self.zmin} - {self.zmax}"
            f"direction = {self.direction}"
        ]
        return "\n   ".join(lines)


class VELO:
    """ Class representing a VELO

    Instance Attributes
    -------------------
    sensors : dict
        Dictionary where a key is an sid and a value is the corresponding sensor.

    Static Attributes
    -----------------
    std_left_z_positions, std_right_z_positions : list
        Standard positions of sensors in a VELO for left and right sensors respectively.

    Public Methods
    --------------
    hits :
        Find all possible hits of a given particle in the VELO.
    efficiency_from_pseudorapidity :
        Find the efficiency given a pseudorapidity.
    efficiency_by_pseudorapidity :
        Given a range of pseudorapidities, get all the corresponding efficencies.
    draw :
        Draw all sensors on an m3d.Axes3D

    Class Methods
    -------------
    standard :
        Create a VELO according to standard specifications.

    """
    std_left_z_positions = [-277, -252, -227, -202, -132, -62, -37, -12, 13, 38, 63,
                            88, 113, 138, 163, 188, 213, 238, 263, 325, 402, 497,
                            616, 661, 706, 751]

    std_right_z_positions = [-289, -264, -239, -214, -144, -74, -49, -24, 1, 26, 51,
                             76, 101, 126, 151, 176, 201, 226, 251, 313, 390, 485,
                             604, 649, 694, 739]

    def __init__(self, left_z_positions, right_z_positions):
        """ Create new VELO instance
        
        Parameters
        ----------
        left_z_positions : list of float and/or int
            Positions of each left sensor of the VELO in the z-direction.
        right_z_positions : list of float and/or int
            Positions of each right sensor of the VELO in the z-direction.
        """
        left_sensors = [VELOSensor.left_sensor(pos) for pos in left_z_positions]
        right_sensors = [VELOSensor.right_sensor(pos) for pos in right_z_positions]
        left_sensors = {sensor.sid: sensor for sensor in left_sensors}
        right_sensors = {sensor.sid: sensor for sensor in right_sensors}

        self.sensors = {}
        self.sensors.update(left_sensors)
        self.sensors.update(right_sensors)

    @classmethod
    def standard(cls):
        """ Create a VELO according to standard specifications. """
        return cls(cls.std_left_z_positions, cls.std_right_z_positions)

    def hits(self, particle, percentage_hit_efficiency=100, gaussian_smear_sigma=0):
        """ For a particle starting at Point(0, 0, 0), find all possible hits in the VELO.

        Parameters
        ----------
        particle : Particle
            The particle to find the hits for.
        percentage_hit_efficiency : int or float, 0 <= val <= 100, optional
            Percentage change of a possible hit actually being detected. Default is 100.

        Returns
        -------
        hits : Hits
            All possible hits in the VELO given the input particle.
        """

        # ENSURE CORRECT START POINT

        zero = Point(0, 0, 0)
        if particle.init_point != zero:
            particle.init_point = zero
            print("Warning: VELO.hits method forces particle init_point to be Point(0, 0, 0)")
        particle.reset()

        # FIND SENSORS IN THE CORRECT DIRECTION

        if particle.four_vector[3] > 0: # move in positive z direction
            # don't care about order, but need to select the sensors in the correct
            # z direction so the particle doesn't move backwards in time from it's initial position
            sensors = [sensor for sensor in self.sensors.values() if sensor.z_position > 0]

        elif particle.four_vector[3] < 0:
            sensors = [sensor for sensor in self.sensors.values() if sensor.z_position < 0]

        else:
            sensors = []

        # FIND HITS

        hits = Hits()
        for sensor in sensors:
            dist = sensor.z_position - particle.point.z()
            particle.advance(dist)
            if particle.point in sensor:
                hits[sensor.sid] = particle.point
        particle.reset()

        if percentage_hit_efficiency != 100:
            hits.apply_hit_efficiency(percentage_hit_efficiency)

        if gaussian_smear_sigma > 0:
            hits.gaussian_smear(gaussian_smear_sigma, gaussian_smear_sigma)

        return hits

    def efficiency_from_pseudorapidity(self, pseudorapidity, percentage_hit_efficiency=100,
                                       angle_of_azimuth_range=np.arange(0, 360, 1)):
        """ Find the percentage reconstructibility efficiency given a pseudorapidity.

        Parameters
        ----------
        pseudorapidity : int or float
            Pseudorapidity value to calculate for.
        percentage_hit_efficiency : int or float, 0 <= val <= 100, optional
            Percentage change of a possible hit actually being detected. Default is 100.
        angle_of_azimuth_range : iterable of int or float values, optional
            The range of angles of azimuth to use to calcualte the efficiency. If you
            change the value of this it is likely to be to adjust the accuracy, by say
            reducing the step in a range.

        Returns
        -------
        percentage_efficiency : int or float
            Percentage reconstructibility efficiency.
        """

        total_count = 0
        reconstructable_count = 0
        momentum_magnitude = 10 # makes no difference since momentum proportions are what matter

        for angle_of_azimuth in angle_of_azimuth_range:

            four_vector = FourVector.from_pseudo(pseudorapidity, angle_of_azimuth,
                                                 momentum_magnitude)
            particle = Particle(four_vector)
            hits = self.hits(particle, percentage_hit_efficiency)
            if len(hits) >= 3:
                reconstructable_count += 1
            total_count += 1

        percentage_efficiency = reconstructable_count / total_count * 100

        return percentage_efficiency

    def efficiency_by_pseudorapidity(self, pseudorapidity_range,
                                     percentage_hit_efficiency=100,
                                     angle_of_azimuth_range=np.arange(0, 360, 1),
                                     progress_bar_enabled=True):
        """ Find the percentage reconstructibility efficiency for each
        pseudorapidity in the given iterable.

        Parameters
        ----------
        pseudorapidity_range : iterable of int or float values
            Pseudorapidity value to calculate for.
        percentage_hit_efficiency : int or float, 0 <= val <= 100, optional
            Percentage change of a possible hit actually being detected. Default is 100.
        angle_of_azimuth_range : iterable of int or float values, optional
            The range of angles of azimuth to use to calcualte the efficiency. If you
            change the value of this it is likely to be to adjust the accuracy, by say
            reducing the step in a range.

        Returns
        -------
        efficiencies : list of int or float values
            Percentage reconstructibility efficiencies corresponding to pseudorapidities.
        """
        
        if progress_bar_enabled:
            progress_count = 0
            progress = ProgressBar(len(pseudorapidity_range), prefix="Progress")

        efficiencies = []

        for pseudorapidity in pseudorapidity_range:

            efficiencies.append(self.efficiency_from_pseudorapidity(pseudorapidity,
                                                                    percentage_hit_efficiency,
                                                                    angle_of_azimuth_range))
            if progress_bar_enabled:
                progress_count += 1
                progress.update(progress_count)

        return np.array(efficiencies)

    def draw(self, ax, alpha=0.2):
        """ Draw all sensors on an m3d.Axes3D """
        for sensor in self.sensors.values():
            sensor.draw(ax, alpha)


class Hits(dict):
    """ Dictionary where sid is key and point of hit is value. """

    def apply_hit_efficiency(self, percentage_hit_efficiency):
        """ Apply a random chance of a hit not being detected such that
        on average the percentage of hits that continue to exist with this dictionary
        afterwards is equal to the percentage_hit_efficiency. """

        for sid in list(self.keys()):
            rand = np.random.random_sample()
            if rand > percentage_hit_efficiency/100:
                del self[sid]

    def gaussian_smear(self, xsigma, ysigma):
        """ gaussian smear all points """
        for key, val in self.items():
            self[key] = val.gaussian_smear(xsigma, ysigma)

    def as_point_array(self):
        """ Return all hit points as an array. """
        return np.array([point for point in self.values()])

    def fit(self):
        """ Fit points to a straight line.

        Returns
        -------
        unit_vector : array
            Unit vector of fit, [x, y, z].
        data_mean : array
            Point on line of fit.
        line_point : lambda function
            Returns points on the best fit line when passed a parameter for
            some scalar distance from the data mean.
        """
        data = self.as_point_array()
        data_mean = np.mean(data, axis=0)

        _, _, vh = np.linalg.svd(data - data_mean)
        unit_vector = vh[0]

        line_point = lambda x, uv=unit_vector, dm=data_mean: uv * x + dm

        return unit_vector, data_mean, line_point


def example_efficiency_comparison(progress_bar_enabled=True):
    """ Comparison of how efficiencies vary with pseudorapidity for
    100% and 98% hit efficiency. """

    # create a VELO sensor
    my_velo = VELO.standard()

    # pseudorapidity range to sample over
    pr_range = np.arange(-7, 7, 0.01)

    # FIND EFFICIENCIES ###############################

    print(r"Finding efficiencies for 100% hit efficiency...")

    efficiency = my_velo.efficiency_by_pseudorapidity(pr_range,
                                                      progress_bar_enabled=progress_bar_enabled)

    print(r"Finding efficiencies for 98% hit efficiency...")

    efficiency2 = my_velo.efficiency_by_pseudorapidity(pr_range, percentage_hit_efficiency=98, 
                                                       progress_bar_enabled=progress_bar_enabled)

    # PLOT RESULTS #####################################

    plt.subplot(2, 1, 1)
    plt.grid(True)
    plt.title("Comparison of how efficiency is dependent on pseudorapidity for two hit efficencies")
    plt.xlabel(r"Pseudorapidity $\eta$")
    plt.ylabel(r"% Efficiency")
    plt.plot(pr_range, efficiency, label="100% Hit Efficiency")
    plt.plot(pr_range, efficiency2, color="red", label="98% Hit Efficiency")
    plt.legend(loc='best')

    plt.subplot(2, 1, 2)
    plt.grid(True)
    plt.title(r"The difference in efficiency between 100% and 98%")
    plt.xlabel(r"Pseudorapidity $\eta$")
    plt.ylabel(r"Difference in % Efficiency")
    plt.plot(pr_range, efficiency - efficiency2)

    plt.tight_layout()
    plt.show()


def example_fit_plot():
    """ Example of how one might fit and plot the hits on a velo. """

    # GENERATE GAUSSIAN SMEARED HITS ######################################

    # create a VELO sensor
    my_velo = VELO.standard()

    # create a particle with a particular four vector
    part = Particle(FourVector(1, 0.9, 1.1, 8))

    # find the hits that would occur with the aformentioned particle
    hits = my_velo.hits(part)

    # smear the hits by a certain standard deviation in x and y
    hits.gaussian_smear(0.012, 0.012)

    # get hits as lists of points for plotting
    data = hits.as_point_array()

    # FIT AND GENERATE POINTS #############################################

    # fit and return lambda function for point on the line
    unit_vector, data_mean, line_point = hits.fit()
    print("unit vector:", unit_vector)
    print("data mean:", data_mean)

    # judge a reasonable length for the line
    line_half_len = Point.scalar_separation(data[0], data[-1]) / 2

    # generate line points
    line_points = np.array([line_point(-line_half_len), line_point(line_half_len)])

    # PLOT ################################################################

    # create a 3D axis
    ax = m3d.Axes3D(plt.figure())

    # plot sensors that were hit
    for sid in hits.keys():
        sens = my_velo.sensors[sid]
        sens.draw(ax, alpha=0.2)

    # plot points
    ax.scatter3D(*data.T)

    # plot fit line
    ax.plot3D(*line_points.T)

    # plot origin
    ax.scatter3D(0, 0, 0, color='red')

    # extra formatting
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel("z")
    plt.gca().set_aspect('equal', adjustable='box')

    # show plot
    plt.show()


if __name__ == "__main__":

    #example_efficiency_comparison(progress_bar_enabled=True)
    example_fit_plot()
