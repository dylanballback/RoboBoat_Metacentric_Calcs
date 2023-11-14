import math

def calculate_modified_ellipse_moment_of_area(a, b, width):
    """
    Calculate the second moment of area about the horizontal axis for a modified ellipse.

    :param a: Semi-major axis of the original ellipse (half of the major diameter).
    :param b: Semi-minor axis of the original ellipse (half of the minor diameter).
    :param width: The width of the rectangle added between the two halves of the ellipse.
    :return: The second moment of area about the horizontal axis.
    """
    # Moment of inertia for one semi-ellipse about its own axis
    I_semi_ellipse = (math.pi / 8) * a * b**3

    # Area of one semi-ellipse
    A_semi_ellipse = (math.pi / 2) * a * b

    # Distance from the centroid of the semi-ellipse to the horizontal axis
    d = (4 * b) / (3 * math.pi)

    # Adjusted moment of inertia for one semi-ellipse about the horizontal axis
    I_adjusted_semi_ellipse = I_semi_ellipse + A_semi_ellipse * d**2

    # Moment of inertia for the rectangle about its own axis
    I_rectangle = (1 / 12) * width * (2 * b)**3

    # Total second moment of area
    I_total = 2 * I_adjusted_semi_ellipse + I_rectangle

    return I_total

# Example usage
a = 50  # Semi-major axis length in mm
b = 30  # Semi-minor axis length in mm
width = 20  # Width of the rectangle in mm

I_total = calculate_modified_ellipse_moment_of_area(a, b, width)
print(f"The second moment of area about the horizontal axis: {I_total} mm^4")





def calculate_ellipse_moment_of_area(a, b):
    """
    Calculate the second moment of area about the major axis for an ellipse.

    :param a: Semi-major axis of the ellipse (half of the major diameter).
    :param b: Semi-minor axis of the ellipse (half of the minor diameter).
    :return: The second moment of area about the major axis.
    """
    I_x = (math.pi / 4) * a * b**3
    return I_x

# Example usage
a = 50  # Semi-major axis length in mm
b = 30  # Semi-minor axis length in mm

I_x = calculate_ellipse_moment_of_area(a, b)
print(f"The second moment of area about the major axis: {I_x} mm^4")




#Get displacement from CAD.

def calculate_metacenter_for_ellipse(a, b, displacement):
    """
    Calculate the metacenter height above the center of buoyancy for an elliptical cross-section.

    :param a: Semi-major axis of the ellipse (half of the major diameter) in meters.
    :param b: Semi-minor axis of the ellipse (half of the minor diameter) in meters.
    :param displacement: Volume of water displaced by the hull (displacement) in cubic meters.
    :return: Height of the metacenter above the center of buoyancy in meters.
    """
    # Calculate the second moment of area for the ellipse
    I = (math.pi / 4) * a * b**3

    # Calculate the height of the metacenter above the center of buoyancy
    BM = I / displacement

    return BM

# Example usage
a = 5  # Semi-major axis in meters
b = 2  # Semi-minor axis in meters
displacement = 300  # Displacement in cubic meters

BM = calculate_metacenter_for_ellipse(a, b, displacement)
print(f"Height of metacenter above center of buoyancy: {BM} meters")




def calculate_metacentric_height(BM, KG):
    """
    Calculate the metacentric height of a vessel.

    :param BM: Height of the metacenter above the center of buoyancy in meters.
    :param KG: Height of the center of gravity above the keel in meters.
    :return: Metacentric height in meters.
    """
    GM = BM - KG
    return GM

# Example usage
BM = 1.5  # Example value in meters
KG = 1.0  # Example value in meters

GM = calculate_metacentric_height(BM, KG)
print(f"Metacentric height (GM): {GM} meters")