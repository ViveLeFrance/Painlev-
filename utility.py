from typing import List, Dict
from sage.all import var, solve, CC, simplify, Expression, point3d, line3d
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np



class LefschetzFibration:

    def __init__(self, variables, domain_equation, fibration_equation) -> None:
        self.variables = variables
        self.domain = domain_equation
        self.fibration = fibration_equation

    def __call__(self, argument):
        return self.fibration.subs(argument)

    def get_critical_points(self):
        """Solves for when the gradient of the domain
        is parallel to the differential of the fibration."""
        G = self.domain
        f = self.fibration
        a = var('a', domain=CC) # additional variable to solve for parallelity


        constraints = [G==0]
        constraints.extend([G.diff(variable) == a*f.diff(variable) for variable in self.variables])

        points = solve(constraints, self.variables + [a])

        return points

    def get_critical_values(self):
        crit_points = self.get_critical_points()
        return list(set([self.__call__(x) for x in crit_points]))

    def get_fibre(self, point, variable=None):
        """Solves for the fibre of a given point in the specified variable.
        If no variable is specified, it will solve for the first variable in the list.
        Note that we assume the fibration to be linear in all variables."""

        if variable is None:
            variable = self.variables[0]

        f = self.fibration
        G = self.domain

        fib_solution = solve(f == point, variable) # f linear, so we expect a single equation

        return G.subs(fib_solution).simplify()
    
    def get_matching_path(self, rho_eq, origin_fibre=None, target_fibre=None, steps=70, solvefor=None, path=None):

        if (origin_fibre is None or target_fibre is None) and path is None:
            raise ValueError("Please provide a path or origin and target fibres.")

        t = var('t', domain=CC)

        if solvefor is None:
            solvefor = self.variables[0]

        fibre_t = self.get_fibre(t, solvefor)
        rho_eq_t = rho_eq.subs(solvefor == fibre_t)

        variables = self.variables
        # variables.remove(solvefor)

        matching_path = {}
        if not path:
            
            for s in np.linspace(0,1,steps):
                fibre_s = fibre_t.subs(t==(1-s)*origin_fibre + s*target_fibre)
                rho_eq_s = rho_eq_t.subs(t==(1-s)*origin_fibre + s*target_fibre)
                rho_s = LefschetzFibration(variables, fibre_s, rho_eq_s)
                matching_path[s] = rho_s.get_critical_values()
        else:
            for index, point in enumerate(path):
                fibre_s = fibre_t.subs(t=point)
                rho_eq_s = rho_eq_t.subs(t=point)
                rho_s = LefschetzFibration(variables, fibre_s, rho_eq_s)
                matching_path[index] = rho_s.get_critical_values()

        return matching_path

    
    
    
def NumericalRoots(expr):
    """Returns the numerical roots of the polynomial 'expr'."""
    coeffs = expr.coefficients(sparse=False)
    coeffs = [complex(coefficient) for coefficient in coeffs]
    return np.polynomial.polynomial.polyroots(coeffs)

def sort_by_angle(points: List[complex], origin_fibre: complex = 0, anticlockwise: bool = False):
    """Sorts a list of points by their argument, starting from the negative real axis."""
    points = [complex(point) for point in points]

    points = sorted(points, key=lambda point: (-np.pi+np.angle(point-origin_fibre))%(2*np.pi))
    if not anticlockwise:
        points.reverse()
    return points

def plot_points_ordered(points: List[complex], title: str = None, fig=None, ax=None, origin_fibre = 0, anticlockwise=False):
    """Plots a list of points on the complex plane, ordered anticlockwise by argument."""
    # Sage complex type is not compatible with python's, but can be coerced
    points = [complex(point) for point in points]

    # Sort points by argument
    points = sort_by_angle(points, origin_fibre=origin_fibre, anticlockwise=anticlockwise)

    real = [point.real for point in points]
    imag = [point.imag for point in points]
    
    if ax is None:
        fig, ax = plt.subplots()

    ax.spines['left'].set_position(('data', 0))
    # Move bottom spine to y=0
    ax.spines['bottom'].set_position(('data', 0))

    # Remove the top and right spines
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')

    for index, point in enumerate(points):
        ax.text(point.real+0.05, point.imag, str(index), fontsize=12, color='blue')
    
    ax.set_title(title)
    
    ax.plot(real, imag, 'ro', markersize=5)
    ax.grid(True)

    return fig, ax
    # plt.show()


def plot_path(path: Dict[complex, List[complex]], title: str = None, origin_fibre=0, anticlockwise=False):
    plot_points_origin = path[0]
    plot_points_target = path[1]

    fig, ax = plot_points_ordered(plot_points_origin, title=title, origin_fibre=origin_fibre, anticlockwise=anticlockwise)

    ax.plot([point.real() for point in plot_points_target], [point.imag() for point in plot_points_target], 'co', markersize=5)

    path_points_x = []
    path_points_y = []
    for step in path.values():

        path_points_x.extend([point.real() for point in step])
        
        path_points_y.extend([point.imag() for point in step])

    ax.plot(path_points_x, path_points_y, 'bo', markersize=2)

    return fig, ax

def plot_path_3d(path: Dict[complex, List[complex]], title: str = None, origin_fibre=0, anticlockwise=False):
   

    points = []

    for index, step in enumerate(path.values()):
        
        for point in step:


            points.append((point.real(), point.imag(), index/len(path)))

    init_points  =[(point.real(), point.imag(), 0) for point in path[0]]
    final_points = [(point.real(), point.imag(), 1) for point in path[1]]
    init_plot = point3d(init_points, size=20, color='red')
    target_plot = point3d(final_points, size=20, color='green')
    trajectory_plot = point3d(points, size=10, color='blue')

    axis_length = 2

    # Create axes (x, y, z) at the origin
    x_axis = line3d([(0,0,0), (axis_length,0,0)], color='black', thickness=1.5)
    y_axis = line3d([(0,0,0), (0,axis_length,0)], color='black', thickness=1.5)
    z_axis = line3d([(0,0,0), (0,0,axis_length)], color='black', thickness=1.5)

    axes = x_axis + y_axis + z_axis

    plot = trajectory_plot + init_plot + target_plot + axes

    plot.show()

    



def pl_path(points: List[complex], steps=70):
    path = []

    for i in range(0, len(points)-1):
        for s in np.linspace(0,1,steps):
            path.append((1-s)*points[i] + s*points[(i+1)])
    return path

def pl_path_1(origin_fibre, target_fibre, offset = None, steps=70, above=True):
    if offset is None:
        offset = np.abs(target_fibre - origin_fibre)/4

    theta = np.arctan(1/2)
    hyp = offset / np.sin(theta)

    if above is True:
        intermediate_point = origin_fibre + hyp*np.exp(1j*theta)
    else:
        intermediate_point = origin_fibre - hyp*np.exp(1j*theta)

    return pl_path([origin_fibre, intermediate_point, target_fibre], steps=steps)


def trace_preimage(rho: LefschetzFibration, t, path: List[complex], title=None, solvefor=None):
    if solvefor is None:
        solvefor = rho.variables[0]

    fibre_rho_t = rho.get_fibre(t, solvefor)
    

    fibres = []

    for s in path:
        fibre_rho_s = fibre_rho_t.subs(t=s)
        fibres.append(NumericalRoots(fibre_rho_s))

    # fibres = np.array(fibres)

    plot_points_real = []
    plot_points_imag = []
    for preimage in fibres:
        plot_points_real.append([value.real for value in preimage])
        plot_points_imag.append([value.imag for value in preimage])

    fig, ax = plt.subplots()    

    ax.spines['left'].set_position(('data', 0))
    # Move bottom spine to y=0
    ax.spines['bottom'].set_position(('data', 0))

    # Remove the top and right spines
    ax.spines['right'].set_color('none')
    ax.spines['top'].set_color('none')

    ax.set_title(title)

    init_sols = NumericalRoots(fibre_rho_t.subs(t==path[0]))
    init_real = [value.real for value in init_sols]
    init_imag = [value.imag for value in init_sols]
    
    final_sols = NumericalRoots(fibre_rho_t.subs(t==path[-1]))
    final_real = [value.real for value in final_sols]
    final_imag = [value.imag for value in final_sols]

    regular_sols = NumericalRoots(fibre_rho_t.subs(t==path[len(path)//2]))
    regular_real = [value.real for value in regular_sols]
    regular_imag = [value.imag for value in regular_sols]


    
    
    ax.plot(init_real, init_imag, 'ro', markersize=5)
    ax.plot(final_real, final_imag, 'co', markersize=5)
    ax.plot(regular_real, regular_imag, 'o', color='purple', markersize=5)


    ax.grid(True)

    ax.plot(plot_points_real, plot_points_imag, 'bo', markersize=2)

    return fig, ax






# x,y,z = var('x, y, z', domain=CC)



# # Equation of affine surface
# G_eq = x**3 + x*y**2 + z**2 -1
# fibration_eq = x+y

# f = LefschetzFibration([x,y,z], G_eq, fibration_eq)

# rho_eq = x

# crit_points_f = f.get_critical_points()
# crit_values_f = f.get_critical_values()

# origin_fibre = 0
# target_fibre = crit_values_f[0]

# matching_path = f.get_matching_path(rho_eq, 0, crit_values_f[0], solvefor=y)

# fig, ax = plot_path(matching_path)

# plt.show()



# # x,y,z,t = var('x, y, z, t', domain=CC)

# # # Equation of affine surface
# # G_eq = x**3 + x*y**2 + z*2 -1
# # fibration_eq = x+y

# # f = LefschetzFibration([x,y,z], G_eq, fibration_eq)

# # print(f.get_fibre(t, y))

