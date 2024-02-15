from typing import List, Dict
from sage.all import var, solve, CC, simplify, Expression
import matplotlib.pyplot as plt
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

        variables = self.variables

        constraints = [G==0]
        constraints.extend([G.diff(variable) == a*f.diff(variable) for variable in variables])

        points = solve(constraints, variables + [a])

        return points

    def get_critical_values(self):
        crit_points = self.get_critical_points()
        return [self.__call__(x) for x in crit_points]

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
    
    def get_matching_path(self, rho_eq, origin_fibre, target_fibre, steps=70, solvefor=None):

        t = var('t', domain=CC)

        fibre_t = self.get_fibre(t, solvefor)
        rho_eq_t = rho_eq.subs(solvefor == fibre_t)

        matching_path = {}

        for s in np.linspace(0,1,steps):
            fibre_s = fibre_t.subs(t==(1-s)*origin_fibre + s*target_fibre)
            rho_eq_s = rho_eq_t.subs(t==(1-s)*origin_fibre + s*target_fibre)
            rho_s = LefschetzFibration(self.variables, fibre_s, rho_eq_s)
            matching_path[s] = rho_s.get_critical_values()

        return matching_path
    
    
def NumericalRoots(expr):
    """Returns the numerical roots of the polynomial 'expr'."""
    coeffs = expr.coefficients(sparse=False)
    coeffs = [complex(coefficient) for coefficient in coeffs]
    return np.polynomial.polynomial.polyroots(coeffs)

def sort_by_angle(points: List[complex], origin_fibre: complex = 0, anticlockwise: bool = False):
    """Sorts a list of points by their argument, anticlockwise."""
    points = [complex(point) for point in points]

    points = sorted(points, key=lambda point: np.pi+np.angle(point-origin_fibre)%(2*np.pi))
    if not anticlockwise:
        points.reverse()
    return points

def plot_points_ordered(points: List[complex], title: str = None, fig=None, ax=None, origin_fibre = 0):
    """Plots a list of points on the complex plane, ordered anticlockwise by argument."""
    # Sage complex type is not compatible with python's, but can be coerced
    points = [complex(point) for point in points]

    # Sort points by argument
    points = sort_by_angle(points, origin_fibre=origin_fibre, anticlockwise=False)

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


def plot_path(path: Dict[complex, List[complex]], title: str = None):
    plot_points_origin = path[0]
    plot_points_target = path[1]

    fig, ax = plot_points_ordered(plot_points_origin, title=title)

    ax.plot([point.real() for point in plot_points_target], [point.imag() for point in plot_points_target], 'co', markersize=5)

    path_points_x = []
    path_points_y = []
    for step in path.values():

        path_points_x.extend([point.real() for point in step])
        
        path_points_y.extend([point.imag() for point in step])

    ax.plot(path_points_x, path_points_y, 'bo', markersize=2)

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

