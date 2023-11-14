

def calculate_center_of_gravity(components, origin):
    """
    Calculate the center of gravity of a boat.

    :param components: A list of tuples, each representing a component.
                       Each tuple contains (name, mass, (x, y, z)) where name is a string,
                       mass is in kg, and (x, y, z) are the coordinates of the component's
                       center of mass in mm relative to the origin.
    :param origin: The origin (0,0,0) point of the boat in mm.
    :return: The center of gravity of the boat as a tuple (x, y, z).
    """
    total_mass = 0
    total_mass_times_x = 0
    total_mass_times_y = 0
    total_mass_times_z = 0

    for name, mass, (x, y, z) in components:
        # Adjust coordinates relative to the origin
        adjusted_x = x - origin[0]
        adjusted_y = y - origin[1]
        adjusted_z = z - origin[2]

        total_mass += mass
        total_mass_times_x += mass * adjusted_x
        total_mass_times_y += mass * adjusted_y
        total_mass_times_z += mass * adjusted_z

    if total_mass == 0:
        raise ValueError("Total mass cannot be zero")

    center_of_gravity = (
        total_mass_times_x / total_mass,
        total_mass_times_y / total_mass,
        total_mass_times_z / total_mass
    )

    return center_of_gravity



components = [
    ("Hull", 0.156, (0, 0, -164)),
    ("T200_1", 0.156, (0, 230.6, -345.5)),
    ("T200_2", 0.156, (0, -250, -345.5)),   # Name, Mass in kg, coordinates in mm
    ("Electronics_Enclosure", 6.3, (1, 2.8, -198)),
    ("Velodyne_Puck", 1.65, (0.58, 436, 46)),
    ("SteroLabs_Camera", 0.23, (-2, 432.5, 98)),
    #need to add vector nav pucks
    #rocket device 
    #rocket antenna
    #e-stop
    #top plate
]


origin = (225, 500, 341)  # Define your origin here
center_of_gravity = calculate_center_of_gravity(components, origin)
print("Center of Gravity (x, y, z):", center_of_gravity, "mm")


def calculate_metacentric_heights(displacement, center_of_gravity, beam, length):
    """
    Calculate the transverse and longitudinal metacentric heights of a boat.

    :param displacement: Volume of water displaced by the boat (V), in cubic meters
    :param center_of_gravity: Height of the boat's center of gravity (G) above the baseline, in meters
    :param beam: Width of the boat at its widest point (B), in meters
    :param length: Length of the boat (L), in meters
    :return: Transverse Metacentric Height (GMt) and Longitudinal Metacentric Height (GMl) in meters
    """
    # Calculate Transverse Metacentric Radius
    BMt = (beam ** 2) / (12 * displacement)

    # Calculate Longitudinal Metacentric Radius
    BMl = (length ** 2) / (12 * displacement)

    # Calculate Metacentric Heights
    GMt = BMt - center_of_gravity
    GMl = BMl - center_of_gravity

    return GMt, GMl

# Example usage
displacement = 0.03322  # in cubic meters (from spreadsheet)
center_of_gravity_z_mm = center_of_gravity[2] # in mm
center_of_gravity_z_m = center_of_gravity_z_mm / 1000 #now in m
beam = 450 /1000  # in meters
length = 1  # in meters


GMt, GMl = calculate_metacentric_heights(displacement, center_of_gravity_z_m, beam, length)
print("Transverse Metacentric Height (GMt):", GMt, "meters")
print("Longitudinal Metacentric Height (GMl):", GMl, "meters")

#Boat Hull overall dimensions 1000 mm in y, 450mm in x, and 341 mm in z

