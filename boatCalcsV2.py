import pprint


"""
Offset and Volumes:
Offset: 10 mm, Volume: 965.3563220419986 cubic mm
Offset: 20 mm, Volume: 2317.0685378752714 cubic mm
Offset: 30 mm, Volume: 3909.4908298255123 cubic mm
Offset: 40 mm, Volume: 5695.029668585284 cubic mm
Offset: 50 mm, Volume: 7674.351555923939 cubic mm
Offset: 60 mm, Volume: 9857.59499545448 cubic mm
Offset: 70 mm, Volume: 12253.627883523002 cubic mm
Offset: 80 mm, Volume: 14872.065670262795 cubic mm
Offset: 90 mm, Volume: 17721.794233459615 cubic mm
Offset: 100 mm, Volume: 20804.429657454522 cubic mm
Offset: 110 mm, Volume: 24073645.383786465 cubic mm
Offset: 120 mm, Volume: 27475.49246326336 cubic mm
Offset: 130 mm, Volume: 30963.05134522178 cubic mm
Offset: 140 mm, Volume: 34488.10201405947 cubic mm
Offset: 150 mm, Volume: 38022.3072296451 cubic mm
Offset: 160 mm, Volume: 41556.60025863865 cubic mm
Offset: 170 mm, Volume: 45091.0022283055 cubic mm
Offset: 180 mm, Volume: 48625.43093226241 cubic mm
Offset: 190 mm, Volume: 52159.814078463794 cubic mm
Offset: 200 mm, Volume: 55693.554378497094 cubic mm
Offset: 210 mm, Volume: 59227.818694045454 cubic mm
Offset: 220 mm, Volume: 62762.083009535374 cubic mm
"""


"""
Offset, Volumes (mm³), Lengths (mm), Widths (mm), Center of Gravity (mm):
Offset: 10 mm, Volume: 965356.3220419986 mm³, Length: 689.7726398685952 mm, Width: 224.77680640344795 mm, Center of Gravity (0.30007935611115477, 0.30001207417669096, -21.651948260848137) (mm):
Offset: 20 mm, Volume: 2317068.5378752714 mm³, Length: 725.9259259319572 mm, Width: 260.97407050926057 mm, Center of Gravity (0.30042431255213514, 0.30050598860047734, -21.086362447499717) (mm):
Offset: 30 mm, Volume: 3909490.829825512 mm³, Length: 755.5555555604138 mm, Width: 284.50813147265524 mm, Center of Gravity (0.3000584904480588, 0.2998450345251152, -20.5172935032466) (mm):
Offset: 40 mm, Volume: 5695029.668585284 mm³, Length: 785.1851851888704 mm, Width: 305.00000000378697 mm, Center of Gravity (0.30004024470440893, 0.29993129704538674, -19.944763855038808) (mm):
Offset: 50 mm, Volume: 7674351.555923939 mm³, Length: 814.8148148173268 mm, Width: 325.0000000032855 mm, Center of Gravity (0.30002998760671734, 0.29991691418795036, -19.363593338890762) (mm):
Offset: 60 mm, Volume: 9857594.995454479 mm³, Length: 844.4444444457836 mm, Width: 345.0000000024301 mm, Center of Gravity (0.30002339272671785, 0.2999504422053776, -18.77191124382883) (mm):
Offset: 70 mm, Volume: 12253627.883523002 mm³, Length: 874.0740740742401 mm, Width: 365.0000000014551 mm, Center of Gravity (0.30001887374173164, 0.2999455686385841, -18.169748920602032) (mm):
Offset: 80 mm, Volume: 14872065.670262795 mm³, Length: 903.7037037026967 mm, Width: 385.0000000005909 mm, Center of Gravity (0.3000156954043761, 0.2999753731497804, -17.5575826262228) (mm):
Offset: 90 mm, Volume: 17721794.233459614 mm³, Length: 933.3333333311532 mm, Width: 405.0000000000182 mm, Center of Gravity (0.30002684511816674, 0.2999171406814513, -16.936108792008042) (mm):
Offset: 100 mm, Volume: 20804429.657454524 mm³, Length: 958.1672533291883 mm, Width: 423.206400854453 mm, Center of Gravity (0.30005534701559133, 0.29996774037973534, -16.307567608627924) (mm):
Offset: 110 mm, Volume: 24073645.383786466 mm³, Length: 975.5961873600427 mm, Width: 436.071065000763 mm, Center of Gravity (0.3000536251203804, 0.300605236740445, -15.681303068732873) (mm):
Offset: 120 mm, Volume: 27475492.46326336 mm³, Length: 987.6938930415545 mm, Width: 444.4790240969344 mm, Center of Gravity (0.30002898677195033, 0.30004366561395945, -15.064239369287108) (mm):
Offset: 130 mm, Volume: 30963051.345221777 mm³, Length: 995.4528167686499 mm, Width: 448.9959570815114 mm, Center of Gravity (0.30033973347825826, 0.3015662460963317, -14.459942683789302) (mm):
Offset: 140 mm, Volume: 34488102.01405947 mm³, Length: 999.3798396948256 mm, Width: 450.0006894684381 mm, Center of Gravity (0.2999844151178106, 0.29996033751281487, -13.87106565364634) (mm):
Offset: 150 mm, Volume: 38022307.229645096 mm³, Length: 1000.0014605249592 mm, Width: 450.0006894684381 mm, Center of Gravity (0.29998290432125474, 0.29997285708362176, -13.297458007547673) (mm):
Offset: 160 mm, Volume: 41556600.25863865 mm³, Length: 1000.0014605249592 mm, Width: 450.0006894684381 mm, Center of Gravity (0.29994748697047574, 0.29980204002222416, -12.73635866271114) (mm):
Offset: 170 mm, Volume: 45091002.228305504 mm³, Length: 1000.0014605249592 mm, Width: 450.0006894684381 mm, Center of Gravity (0.2999406075908592, 0.29971970734930936, -12.18482347320143) (mm):
Offset: 180 mm, Volume: 48625430.93226241 mm³, Length: 1000.0014605249592 mm, Width: 450.0006894684381 mm, Center of Gravity (0.29993503590954723, 0.299698939147548, -11.640776937173493) (mm):
Offset: 190 mm, Volume: 52159814.07846379 mm³, Length: 1000.0014605249592 mm, Width: 450.0006894684381 mm, Center of Gravity (0.29992442745707876, 0.29964583328127853, -11.102704693845583) (mm):
Offset: 200 mm, Volume: 55693554.378497094 mm³, Length: 1000.0014605249592 mm, Width: 450.0006894684381 mm, Center of Gravity (0.299989056197164, 0.30007448126364084, -10.569530003562644) (mm):
Offset: 210 mm, Volume: 59227818.694045454 mm³, Length: 1000.0014605249592 mm, Width: 450.0006894684381 mm, Center of Gravity (0.2999893339156391, 0.3000878705387151, -10.040264126828731) (mm):
Offset: 220 mm, Volume: 62762083.00953537 mm³, Length: 1000.0014605249592 mm, Width: 450.0006894684381 mm, Center of Gravity (0.29998958035467427, 0.3000997519139576, -9.514294295690211) (mm):


"""


# first value is mm from bottom (water line), volume in cm ^3 
#40mm wider and 75mm shorter 
Offset_and_Volumes_40_75={
    10: (965356.3220419986, 689.7726398685952, 224.77680640344795, 0.30007935611115477, 0.30001207417669096, -21.651948260848137), 
    20: (2317068.5378752714, 725.9259259319572, 260.97407050926057, 0.30042431255213514, 0.30050598860047734, -21.086362447499717), 
    30: (3909490.829825512, 755.5555555604138, 284.50813147265524, 0.3000584904480588, 0.2998450345251152, -20.5172935032466), 
    40: (5695029.668585284, 785.1851851888704, 305.00000000378697, 0.30004024470440893, 0.29993129704538674, -19.944763855038808), 
    50: (7674351.555923939, 814.8148148173268, 325.0000000032855, 0.30002998760671734, 0.29991691418795036, -19.363593338890762), 
    60: (9857594.995454479, 844.4444444457836, 345.0000000024301, 0.30002339272671785, 0.2999504422053776, -18.77191124382883), 
    70: (12253627.883523002, 874.0740740742401, 365.0000000014551, 0.30001887374173164, 0.2999455686385841, -18.169748920602032), 
    80: (14872065.670262795, 903.7037037026967, 385.0000000005909, 0.3000156954043761, 0.2999753731497804, -17.5575826262228), 
    90: (17721794.233459614, 933.3333333311532, 405.0000000000182, 0.30002684511816674, 0.2999171406814513, -16.936108792008042), 
    100: (20804429.657454524, 958.1672533291883, 423.206400854453, 0.30005534701559133, 0.29996774037973534, -16.307567608627924), 
    110: (24073645.383786466, 975.5961873600427, 436.071065000763, 0.3000536251203804, 0.300605236740445, -15.681303068732873), 
    120: (27475492.46326336, 987.6938930415545, 444.4790240969344, 0.30002898677195033, 0.30004366561395945, -15.064239369287108), 
    130: (30963051.345221777, 995.4528167686499, 448.9959570815114, 0.30033973347825826, 0.3015662460963317, -14.459942683789302), 
    140: (34488102.01405947, 999.3798396948256, 450.0006894684381, 0.2999844151178106, 0.29996033751281487, -13.87106565364634), 
    150: (38022307.229645096, 1000.0014605249592, 450.0006894684381, 0.29998290432125474, 0.29997285708362176, -13.297458007547673), 
    160: (41556600.25863865, 1000.0014605249592, 450.0006894684381, 0.29994748697047574, 0.29980204002222416, -12.73635866271114), 
    170: (45091002.228305504, 1000.0014605249592, 450.0006894684381, 0.2999406075908592, 0.29971970734930936, -12.18482347320143), 
    180: (48625430.93226241, 1000.0014605249592, 450.0006894684381, 0.29993503590954723, 0.299698939147548, -11.640776937173493), 
    190: (52159814.07846379, 1000.0014605249592, 450.0006894684381, 0.29992442745707876, 0.29964583328127853, -11.102704693845583), 
    200: (55693554.378497094, 1000.0014605249592, 450.0006894684381, 0.299989056197164, 0.30007448126364084, -10.569530003562644), 
    210: (59227818.694045454, 1000.0014605249592, 450.0006894684381, 0.2999893339156391, 0.3000878705387151, -10.040264126828731), 
    220: (62762083.00953537, 1000.0014605249592, 450.0006894684381, 0.29998958035467427, 0.3000997519139576, -9.514294295690211)
    }



def calculate_and_append_mass_displaced(volume_dict):
    """
    Calculate the mass displaced for each volume in the dictionary and append it to the dictionary.
    Assumes the volumes are in mm^3.

    :param volume_dict: Dictionary with offsets as keys and tuples of (volume, y_length, x_width, cg_x, cg_y, cg_z) as values.
    :return: Dictionary with offsets as keys and tuples of (volume, mass displaced, y_length, x_width, cg_x, cg_y, cg_z) as values.
    """
    density_of_water = 1000  # kg/m^3
    conversion_factor = 1e-9  # Convert mm^3 to m^3

    # Calculate and append mass displaced in kg for each volume displaced
    for offset, (volume, y_length, x_width, cg_x, cg_y, cg_z) in volume_dict.items():
        mass_displaced = volume * conversion_factor * density_of_water
        volume_dict[offset] = (volume, mass_displaced, y_length, x_width, cg_x, cg_y, cg_z)

    return volume_dict



Calculated_40_75 = calculate_and_append_mass_displaced(Offset_and_Volumes_40_75)
# water height from bottom of hull, volume of water displaced, mass displaced in Kg)
pp = pprint.PrettyPrinter(depth=7)
pp.pprint(Calculated_40_75)


def calculate_metacentric_heights_for_all_offsets(offset_volumes_dict, components, origin):
    metacentric_heights_results = {}

    for offset, (volume_mm3, mass_kg, y_length_mm, x_width_mm, cg_x_mm, cg_y_mm, cg_z_mm) in offset_volumes_dict.items():
        # Update Hull's CG and Calculate Center of Gravity
        hull_cg = (cg_x_mm, cg_y_mm, cg_z_mm)
        center_of_gravity = calculate_center_of_gravity(components, origin, hull_cg)

        # Convert displacement volume from mm³ to m³
        displacement_m3 = volume_mm3 * 1e-9

        # Convert beam and length from mm to m
        beam_m = x_width_mm / 1000
        length_m = y_length_mm / 1000

        # Calculate Metacentric Heights
        GMt, GMl = calculate_metacentric_heights(displacement_m3, center_of_gravity[2] / 1000, beam_m, length_m)
        
        # Store the results
        metacentric_heights_results[offset] = (GMt, GMl)

    return metacentric_heights_results


def calculate_center_of_gravity(components, origin, hull_cg):
    
    """
    Calculate the center of gravity of a boat.

    :param components: A list of tuples, each representing a component.
                       Each tuple contains (name, mass, (x, y, z)) where name is a string,
                       mass is in kg, and (x, y, z) are the coordinates of the component's
                       center of mass in mm relative to the origin.
    :param origin: The origin (0,0,0) point of the boat in mm.
    :return: The center of gravity of the boat as a tuple (x, y, z).
    """
    # Update the center of gravity of the hull component
    for i, (name, mass, _) in enumerate(components):
        if name == "Hull":
            components[i] = (name, mass, hull_cg)

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



def calculate_all_metacentric_heights(offset_volumes_dict):
    """
    Calculate the metacentric heights for each offset height in the dictionary.

    :param offset_volumes_dict: Dictionary with offset heights and tuples of parameters.
    :return: Dictionary with offset heights and metacentric heights.
    """
    metacentric_heights_dict = {}

    # Initialize components list for each iteration to avoid side effects
    for offset, (volume_mm3, mass_kg, y_length_mm, x_width_mm, cg_x_mm, cg_y_mm, cg_z_mm) in offset_volumes_dict.items():
        components = [
            ("T200_1", 0.156, (0, 230.6, -49)),
            ("T200_2", 0.156, (0, -250, -49)),
            ("Electronics_Enclosure", 6.3, (1, 2.8, 88)),
            ("Velodyne_Puck", 1.65, (0.58, 436, 267)),
            ("SteroLabs_Camera", 0.23, (-2, 432.5, 319)),
            ("water_pump", 2, (0, 0, 100))
            # Add other components here
        ]

        hull_cg = (cg_x_mm, cg_y_mm, cg_z_mm)
        center_of_gravity = calculate_center_of_gravity(components, origin, hull_cg)# Convert displacement volume from mm³ to m³
        displacement_m3 = volume_mm3 * 1e-9

        # Convert CG's z-coordinate from mm to m
        cg_z_m = cg_z_mm / 1000

        # Convert beam and length from mm to m
        beam_m = x_width_mm / 1000
        length_m = y_length_mm / 1000

        # Calculate Metacentric Heights
        GMt, GMl = calculate_metacentric_heights(displacement_m3, cg_z_m, beam_m, length_m)
        
        # Store the results in the dictionary
        metacentric_heights_dict[offset] = (GMt, GMl)

    return metacentric_heights_dict


"""
components = [
    ("T200_1", 0.156, (0, 230.6, -49)),
    ("T200_2", 0.156, (0, -250, -49)),   # Name, Mass in kg, coordinates in mm
    ("Electronics_Enclosure", 6.3, (1, 2.8, 88)),
    ("Velodyne_Puck", 1.65, (0.58, 436, 267)),
    ("SteroLabs_Camera", 0.23, (-2, 432.5, 319)),
    ("water_pump", 2, (0, 0, 100))
    #need to add vector nav pucks
    #rocket device 
    #rocket antenna
    #e-stop
    #top plate
]
"""

origin = (0, 0, 0)  # Define your origin here

# Calculate Metacentric Heights for each offset
metacentric_heights_results = calculate_metacentric_heights_for_all_offsets(Offset_and_Volumes_40_75, components, origin)

# Printing the results
for offset, (GMt, GMl) in metacentric_heights_results.items():
    print(f"Offset: {offset} mm, GMt: {GMt} m, GMl: {GMl} m")