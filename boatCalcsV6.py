import pprint
import math
import matplotlib.pyplot as plt


"""
11/14/23

I fixed the script and got everthing working well. 
Now I am going to try to use Matplotlib to plot the CG, CB, and GM  
"""






'''
Offset, Volumes (mm³), Lengths (mm), Widths (mm), Center of Gravity (mm):
Offset: 10 mm, Volume: 1227618.7275407966 mm³, Length: 689.7726398680547 mm, Width: 264.7768064034481 mm, Center of Gravity (0.0007455613112383741, -0.0020769311218571396, 5.401104627582356) (mm):
Offset: 20 mm, Volume: 2863007.7478833557 mm³, Length: 725.9259259316937 mm, Width: 300.97407050926057 mm, Center of Gravity (0.003642049026305605, -0.0017365549276188164, 10.97118496853593) (mm):
Offset: 30 mm, Volume: 4751677.9055522075 mm³, Length: 755.5555555623795 mm, Width: 324.5081314726553 mm, Center of Gravity (0.0004402732916969354, -0.0012872871220881432, 16.585421326249197) (mm):
Offset: 40 mm, Volume: 6845471.649743572 mm³, Length: 785.1851851902209 mm, Width: 345.00000000378697 mm, Center of Gravity (0.0002713390311137598, -0.0006231197070011296, 22.24213936579634) (mm):
Offset: 50 mm, Volume: 9144732.941907631 mm³, Length: 814.814814819635 mm, Width: 365.00000000328544 mm, Center of Gravity (0.00015588681306966112, -0.0007466678736051113, 27.983969458238818) (mm):
Offset: 60 mm, Volume: 11659800.778181601 mm³, Length: 844.4444444469736 mm, Width: 385.0000000024301 mm, Center of Gravity (0.00034348257339278376, 0.0008841706576423469, 33.82712039453754) (mm):
Offset: 70 mm, Volume: 14399619.405949375 mm³, Length: 874.0740740742929 mm, Width: 405.0000000014551 mm, Center of Gravity (0.0006347659895206137, -0.0007930205699890891, 39.77171884419292) (mm):
Offset: 80 mm, Volume: 17373587.875892773 mm³, Length: 903.7037037021095 mm, Width: 425.0000000005909 mm, Center of Gravity (2.2765957091719358e-05, -0.00040760721615462625, 45.81344096345305) (mm):
Offset: 90 mm, Volume: 20590721.978321414 mm³, Length: 933.3333333308726 mm, Width: 445.00000000001813 mm, Center of Gravity (-9.02801847424417e-05, -0.0007583381565617886, 51.94637281212028) (mm):
Offset: 100 mm, Volume: 24053619.77637813 mm³, Length: 958.167253329188 mm, Width: 463.206400854453 mm, Center of Gravity (0.006850542641622104, -0.0004841940763949948, 58.152473307082175) (mm):
Offset: 110 mm, Volume: 27708457.14891749 mm³, Length: 975.5961873600427 mm, Width: 476.0710650007629 mm, Center of Gravity (-0.00028157406882167724, 0.00046411391929512824, 64.33738043911725) (mm):
Offset: 120 mm, Volume: 31502815.83840659 mm³, Length: 987.6938930415545 mm, Width: 484.4790240969344 mm, Center of Gravity (-3.336666802411514e-05, 0.0004117935778169546, 70.44241628858252) (mm):
Offset: 130 mm, Volume: 35386969.66961597 mm³, Length: 995.4528167686499 mm, Width: 488.99595708151134 mm, Center of Gravity (0.0038850345564561684, 0.02114648528543317, 76.43171009037643) (mm):
Offset: 140 mm, Volume: 39311247.09433885 mm³, Length: 999.3798396948255 mm, Width: 490.0006894684381 mm, Center of Gravity (-0.0005104078823858771, -0.00040308776796338375, 82.27959522194935) (mm):
Offset: 150 mm, Volume: 43245352.11047468 mm³, Length: 1000.0014605249592 mm, Width: 490.0006894684381 mm, Center of Gravity (-0.000638247926454949, -0.0008325387840374399, 87.98539733020343) (mm):
Offset: 160 mm, Volume: 47179711.188874096 mm³, Length: 1000.0014605249592 mm, Width: 490.0006894684381 mm, Center of Gravity (-0.0009215148033181464, -0.002050243627852999, 93.57380449382013) (mm):
Offset: 170 mm, Volume: 51114130.06277282 mm³, Length: 1000.0014605249592 mm, Width: 490.0006894684381 mm, Center of Gravity (-0.0010079836818649435, -0.002917694105769586, 99.0716986262689) (mm):
Offset: 180 mm, Volume: 55048581.823114686 mm³, Length: 1000.0014605249592 mm, Width: 490.0006894684381 mm, Center of Gravity (-0.0010686478178678627, -0.003156368496587958, 104.49845689839452) (mm):
Offset: 190 mm, Volume: 58982911.310619205 mm³, Length: 1000.0014605249592 mm, Width: 490.0006894684381 mm, Center of Gravity (-0.00111565135249736, -0.0042268688977542546, 109.8681611091994) (mm):
Offset: 200 mm, Volume: 62917294.878949516 mm³, Length: 1000.0014605249592 mm, Width: 490.0006894684381 mm, Center of Gravity (-0.0012098928783949026, -0.004796256089857898, 115.19167841436136) (mm):
Offset: 210 mm, Volume: 66850795.06888761 mm³, Length: 1000.0014605249592 mm, Width: 490.0006894684381 mm, Center of Gravity (-0.00031589449527247293, 0.00014216838284064615, 120.476356689737) (mm):
Offset: 220 mm, Volume: 70785040.81976028 mm³, Length: 1000.0014605249592 mm, Width: 490.0006894684381 mm, Center of Gravity (-0.00030013773999859605, 0.0001948827355430227, 125.72999838812248) (mm):

'''

# Data representing offset heights from the bottom, volume, dimensions, and center of gravity for each offset
#40mm wider and 75mm shorter 
Offset_and_Volumes_40_75={
    # Each entry in the dictionary represents data for a specific offset from the waterline
    # The values are: volume in mm³, length in mm, width in mm, and center of gravity coordinates in mm (x, y, z)
    # Example: 10: (volume_mm3, y_length_mm, x_width_mm, cg_x_mm, cg_y_mm, cg_z_mm)
    10: (1227618.7275407966, 689.7726398680547, 264.7768064034481, 0.0007455613112383741, -0.0020769311218571396, 5.401104627582356), 
    20: (2863007.7478833557, 725.9259259316937, 300.97407050926057, 0.003642049026305605, -0.0017365549276188164, 10.97118496853593), 
    30: (4751677.9055522075, 755.5555555623795, 324.5081314726553, 0.0004402732916969354, -0.0012872871220881432, 16.585421326249197), 
    40: (6845471.649743572, 785.1851851902209, 345.00000000378697, 0.0002713390311137598, -0.0006231197070011296, 22.24213936579634), 
    50: (9144732.941907631, 814.814814819635, 365.00000000328544, 0.00015588681306966112, -0.0007466678736051113, 27.983969458238818), 
    60: (11659800.778181601, 844.4444444469736, 385.0000000024301, 0.00034348257339278376, 0.0008841706576423469, 33.82712039453754), 
    70: (14399619.405949375, 874.0740740742929, 405.0000000014551, 0.0006347659895206137, -0.0007930205699890891, 39.77171884419292), 
    80: (17373587.875892773, 903.7037037021095, 425.0000000005909, 2.2765957091719358e-05, -0.00040760721615462625, 45.81344096345305), 
    90: (20590721.978321414, 933.3333333308726, 445.00000000001813, -9.02801847424417e-05, -0.0007583381565617886, 51.94637281212028), 
    100: (24053619.77637813, 958.167253329188, 463.206400854453, 0.006850542641622104, -0.0004841940763949948, 58.152473307082175), 
    110: (27708457.14891749, 975.5961873600427, 476.0710650007629, -0.00028157406882167724, 0.00046411391929512824, 64.33738043911725), 
    120: (31502815.83840659, 987.6938930415545, 484.4790240969344, -3.336666802411514e-05, 0.0004117935778169546, 70.44241628858252), 
    130: (35386969.66961597, 995.4528167686499, 488.99595708151134, 0.0038850345564561684, 0.02114648528543317, 76.43171009037643), 
    140: (39311247.09433885, 999.3798396948255, 490.0006894684381, -0.0005104078823858771, -0.00040308776796338375, 82.27959522194935), 
    150: (43245352.11047468, 1000.0014605249592, 490.0006894684381, -0.000638247926454949, -0.0008325387840374399, 87.98539733020343), 
    160: (47179711.188874096, 1000.0014605249592, 490.0006894684381, -0.0009215148033181464, -0.002050243627852999, 93.57380449382013), 
    170: (51114130.06277282, 1000.0014605249592, 490.0006894684381, -0.0010079836818649435, -0.002917694105769586, 99.0716986262689), 
    180: (55048581.823114686, 1000.0014605249592, 490.0006894684381, -0.0010686478178678627, -0.003156368496587958, 104.49845689839452), 
    190: (58982911.310619205, 1000.0014605249592, 490.0006894684381, -0.00111565135249736, -0.0042268688977542546, 109.8681611091994), 
    200: (62917294.878949516, 1000.0014605249592, 490.0006894684381, -0.0012098928783949026, -0.004796256089857898, 115.19167841436136), 
    210: (66850795.06888761, 1000.0014605249592, 490.0006894684381, -0.00031589449527247293, 0.00014216838284064615, 120.476356689737), 
    220: (70785040.81976028, 1000.0014605249592, 490.0006894684381, -0.00030013773999859605, 0.0001948827355430227, 125.72999838812248)
    }


# Function to calculate and append the mass displaced for each volume in the dictionary
def calculate_and_append_mass_displaced(volume_dict):
    density_of_water = 1000  # kg/m^3
    conversion_factor = 1e-9  # Convert mm^3 to m^3

    # Loop through each entry in the dictionary
    for offset, (volume, y_length, x_width, cg_x, cg_y, cg_z) in volume_dict.items():
        # Calculate mass displaced using the formula: volume * conversion_factor * density_of_water
        mass_displaced = volume * conversion_factor * density_of_water
        # Update the dictionary entry with the calculated mass
        volume_dict[offset] = (volume, mass_displaced, y_length, x_width, cg_x, cg_y, cg_z)

    return volume_dict

# Apply the function to the dictionary
Calculated_40_75 = calculate_and_append_mass_displaced(Offset_and_Volumes_40_75)

# Pretty print the updated dictionary
pp = pprint.PrettyPrinter(depth=7)
pp.pprint(Calculated_40_75)


# Function to calculate the center of gravity of a boat
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
    # Initialize total mass and moment accumulators
    total_mass = 0
    total_mass_times_x = 0
    total_mass_times_y = 0
    total_mass_times_z = 0

    # Loop through each component to calculate the total mass and moments
    for name, mass, (x, y, z) in components:
        adjusted_x = x - origin[0]
        adjusted_y = y - origin[1]
        adjusted_z = z - origin[2]

        total_mass += mass
        total_mass_times_x += mass * adjusted_x
        total_mass_times_y += mass * adjusted_y
        total_mass_times_z += mass * adjusted_z

    # Check if total mass is zero to avoid division by zero error
    if total_mass == 0:
        raise ValueError("Total mass cannot be zero")

    # Calculate the center of gravity
    center_of_gravity = (
        total_mass_times_x / total_mass,
        total_mass_times_y / total_mass,
        total_mass_times_z / total_mass
    )

    return center_of_gravity




def ellipse_moment_of_area(a, b, added_width):
    """
    Calculate the second moment of area for a modified ellipse.

    :param a: Semi-major axis of the original ellipse in mm.
    :param b: Semi-minor axis of the original ellipse in mm.
    :param added_width: The width added to the ellipse in mm.
    :return: The second moment of area in mm^4.
    """
    # Second moment of area for the original ellipse
    I_ellipse = (math.pi / 4) * a * (b ** 3)
    #print(f"I_ellipse = {I_ellipse} for A= {a}, B = {b}, and added width= {added_width}")

    # Adjust a for the added width
    a += added_width / 2

    # Second moment of area for the added rectangle
    height = 2 * a  # Full length of the ellipse
    I_rectangle = (1 / 12) * height * added_width**3
    #print(f"I_rectangle = {I_rectangle} for added_width = {added_width}, height = {height}")
    # Total second moment of area
    I_total =  I_ellipse + I_rectangle
    return I_ellipse, I_rectangle, I_total

"""
#This was used for testing the ellipse_moment_of_area function 

def calculate_all_moments_of_area(offset_volumes_dict, added_width):
    moments_of_area_dict = {}

    for offset, (_, _, y_length_mm, x_width_mm, _, _, _) in offset_volumes_dict.items():
        # a and b values for the ellipse
        a = y_length_mm / 2
        b = (x_width_mm - added_width) / 2

        # Calculate the moment of area
        I_ellipse, I_rectangle, I_total = ellipse_moment_of_area(a, b, added_width)
        moments_of_area_dict[offset] = I_total

        # Store the individual and total moments of area in the dictionary
        moments_of_area_dict[offset] = {
            "I_ellipse": I_ellipse,
            "I_rectangle": I_rectangle,
            "I_total": I_total
        }

    return moments_of_area_dict

# Assuming 40mm added width
added_width = 40
moments_of_area_results = calculate_all_moments_of_area(Calculated_40_75, added_width)

# Print the results
for offset, I in moments_of_area_results.items():
    print(f"Offset: {offset} mm, Second Moment of Area: {I} mm^4")

"""



# Function to calculate the metacentric heights
def calculate_metacentric_heights(displacement, center_of_gravity_z, center_of_buoyancy_z, I, beam, length):
    """
    Calculate the transverse and longitudinal metacentric heights of a boat.

    :param displacement: Volume of water displaced by the boat (V), in cubic meters.
    :param center_of_gravity_z: Height of the boat's center of gravity (G) above the baseline, in meters.
    :param center_of_buoyancy_z: Height of the boat's center of buoyancy (B) above the baseline, in meters.
    :param I: Second moment of area of the waterline plane in mm^4.
    :param beam: Width of the boat at its widest point (B), in meters.
    :param length: Length of the boat (L), in meters.
    :return: Transverse Metacentric Height (GMt) and Longitudinal Metacentric Height (GMl) in meters.
    """
    # Convert I from mm^4 to m^4
    I_m4 = I * 1e-12

    # Calculate Transverse and Longitudinal Metacentric Radii
    BMt = I_m4 / (beam * displacement)
    BMl = I_m4 / (length * displacement)

    # Calculate Metacentric Heights
    GMt_m = BMt - (center_of_gravity_z - center_of_buoyancy_z)
    GMl_m = BMl - (center_of_gravity_z - center_of_buoyancy_z)

    # Convert Metacentric Heights from m to mm
    GMt_mm = GMt_m * 1000
    GMl_mm = GMl_m * 1000

    return GMt_mm, GMl_mm
    
    
# Loop through each offset in the dictionary


# Function to calculate metacentric heights for all offsets
def calculate_all_metacentric_heights(offset_volumes_dict):
    metacentric_heights_dict = {}
    moments_of_area_dict = {}
    cg_dict = {}  # Dictionary to store the center of gravity
    # Loop through each offset in the dictionary
    for offset, (volume_mm3, mass_kg, y_length_mm, x_width_mm, cb_x_mm, cb_y_mm, cb_z_mm) in offset_volumes_dict.items():
        # Define the components of the boat including their mass and location
        components = [
            ("Hull", 5, (0, 0, 128.3)),
            ("T200_1", 0.156, (0, 230.6, -49)),
            ("T200_2", 0.156, (0, -250, -49)),
            ("Electronics_Enclosure", 6.3, (1, 2.8, 88)),
            ("Velodyne_Puck", 1.65, (0.58, 436, 267)),
            ("SteroLabs_Camera", 0.23, (-2, 432.5, 319)),
            ("water_pump", 2, (0, 0, 100))
            # Add other components here
        ]

        origin = (0, 0, 0)
        constant_cg = calculate_center_of_gravity(components, origin)
        # Store the CG in the dictionary
        cg_dict[offset] = constant_cg

        # Center of buoyancy Z coordinate in meters
        center_of_buoyancy_z_m = cb_z_mm / 1000

        # Assuming 40mm added width
        added_width = 40

        # Calculate the moment of area for the current offset
        a = y_length_mm / 2
        b = x_width_mm / 2
        I_ellipse, I_rectangle, I_total = ellipse_moment_of_area(a, b, added_width)

        # Store the individual and total moments of area in the dictionary
        moments_of_area_dict[offset] = {
            "I_ellipse": I_ellipse,
            "I_rectangle": I_rectangle,
            "I_total": I_total
        }

        # Convert displacement volume from mm³ to m³
        displacement_m3 = volume_mm3 * 1e-9
        # Convert beam and length from mm to m 
        beam_m = x_width_mm / 1000
        length_m = y_length_mm / 1000

        # Calculate metacentric heights for the current offset
        GMt_mm, GMl_mm = calculate_metacentric_heights(displacement_m3, constant_cg[2] / 1000, center_of_buoyancy_z_m, I_total, beam_m, length_m)

        # Store the results in the dictionary
        metacentric_heights_dict[offset] = (GMt_mm, GMl_mm)

    return metacentric_heights_dict, moments_of_area_dict, cg_dict


# Calculate metacentric heights and CG for all offsets
metacentric_heights_results, moments_of_area_results, cg_results = calculate_all_metacentric_heights(Calculated_40_75)


"""
# Print the results
for offset, (GMt, GMl) in metacentric_heights_results[0].items():
    print(f"Offset: {offset} mm, GMt: {GMt} m, GMl: {GMl} m")

# Print the results
for offset, moment_data in metacentric_heights_results[1].items():
    # This loop will go through each key-value pair in the dictionary moments_of_area_results.
    # 'offset' will be the key (e.g., 10, 20, 30, etc.).
    # 'moment_data' will be the value, which is itself a dictionary containing 'I_ellipse', 'I_rectangle', and 'I_total'.

    I_ellipse = moment_data['I_ellipse'] # Access the 'I_ellipse' value from the inner dictionary.
    I_rectangle = moment_data['I_rectangle'] # Access the 'I_rectangle' value from the inner dictionary.
    I_total = moment_data['I_total']  # Access the 'I_total' value from the inner dictionary.
    # Print the extracted values along with the offset.
    print(f"Offset: {offset} mm, Ellipse: {I_ellipse}, Rectangle: {I_rectangle}, Total Second Moment of Area: {I_total} mm^4")
"""


def print_all_data(offset_volumes_dict, combined_results):
    metacentric_heights_dict, moments_of_area_dict, cg_dict = combined_results

    for offset in offset_volumes_dict:
        # Retrieve data from dictionaries
        volume_data = offset_volumes_dict[offset]
        metacentric_data = metacentric_heights_dict.get(offset, (None, None))
        moment_data = moments_of_area_dict.get(offset, {})
        cg_data = cg_dict.get(offset, (0, 0, 0))  # Get CG data, default to (0,0,0)

        # Unpack data
        volume_mm3, mass_kg, y_length_mm, x_width_mm, cb_x_mm, cb_y_mm, cb_z_mm = volume_data
        GMt_mm, GMl_mm = metacentric_data
        I_ellipse = moment_data.get('I_ellipse', 0)
        I_rectangle = moment_data.get('I_rectangle', 0)
        I_total = moment_data.get('I_total', 0)

        # Print data for each offset
        print(f"Offset: {offset} mm")
        print(f"Volume: {volume_mm3:.2f} mm³, Mass: {mass_kg:.2f} kg")
        print(f"Dimensions (Length x Width): {y_length_mm:.2f} mm x {x_width_mm:.2f} mm")
        print(f"Center of Buoyancy (X, Y, Z): {cb_x_mm:.2f} mm, {cb_y_mm:.2f} mm, {cb_z_mm:.2f} mm")
        print(f"Center of Gravity (X, Y, Z): {cg_data[0]:.2f} mm, {cg_data[1]:.2f} mm, {cg_data[2]:.2f} mm")
        print(f"Metacentric Heights - Transverse Metacentric Height (GMt): {GMt_mm:.2f} mm, Longitudinal Metacentric Height (GMl): {GMl_mm:.2f} mm")
        print(f"Moments of Area - Ellipse: {I_ellipse:.2f} mm⁴, Rectangle: {I_rectangle:.2f} mm⁴, Total: {I_total:.2f} mm⁴")
        print("-" * 50)

# Example usage
print_all_data(Calculated_40_75, (metacentric_heights_results, moments_of_area_results, cg_results))


#just plot CG CB and Metacentric Height GMt
def plot_stability_data(offset_volumes_dict, combined_results):
    """
    Plots the center of gravity, center of buoyancy, and metacentric height for each offset.

    :param offset_volumes_dict: Dictionary with volume and dimension data.
    :param combined_results: Tuple containing dictionaries of metacentric heights, moments of area, and center of gravity.
    """
    metacentric_heights_dict, _, cg_dict = combined_results

    offsets = list(offset_volumes_dict.keys())
    center_of_gravity_z = [cg_dict[offset][2] for offset in offsets]  # Z component in mm
    center_of_buoyancy_z = [offset_volumes_dict[offset][6] for offset in offsets]  # Z component in mm
    GMt = [metacentric_heights_dict[offset][0] for offset in offsets]  
    plt.figure(figsize=(10, 6))

    plt.plot(offsets, center_of_gravity_z, label='Center of Gravity (Z)', color='blue')
    plt.plot(offsets, center_of_buoyancy_z, label='Center of Buoyancy (Z)', color='green')
    plt.plot(offsets, GMt, label='Metacentric Height (GMt)', color='red')

    plt.xlabel('Offset from bottom (mm)')
    plt.ylabel('Height (mm)')
    plt.title('Boat Stability Data (40mm wider & 75mm shorter)')
    plt.legend()
    plt.grid(True)
    plt.show()


plot_stability_data(Calculated_40_75, (metacentric_heights_results, moments_of_area_results, cg_results))



def plot_stability_data(offset_volumes_dict, combined_results, cutoff_offset=0):
    """
    Plots the center of gravity, center of buoyancy, metacentric height, and mass for each offset above a certain cutoff.

    :param offset_volumes_dict: Dictionary with volume and dimension data.
    :param combined_results: Tuple containing dictionaries of metacentric heights, moments of area, and center of gravity.
    :param cutoff_offset: The offset value from which to start plotting (default is 0).
    """
    metacentric_heights_dict, _, cg_dict = combined_results

    # Filter data based on the cutoff_offset
    filtered_offsets = [offset for offset in offset_volumes_dict.keys() if offset >= cutoff_offset]
    
    center_of_gravity_z = [cg_dict[offset][2] for offset in filtered_offsets]
    center_of_buoyancy_z = [offset_volumes_dict[offset][6] for offset in filtered_offsets]
    GMt = [metacentric_heights_dict[offset][0] for offset in filtered_offsets] 
    mass = [offset_volumes_dict[offset][1] for offset in filtered_offsets]

    fig, ax1 = plt.subplots(figsize=(10, 6))

    ax1.plot(filtered_offsets, center_of_gravity_z, label='Center of Gravity (Z)', color='blue')
    ax1.plot(filtered_offsets, center_of_buoyancy_z, label='Center of Buoyancy (Z)', color='green')
    ax1.plot(filtered_offsets, GMt, label='Transverse Metacentric Height (GMt)', color='red')
    ax1.set_xlabel('Offset from bottom (mm)')
    ax1.set_ylabel('Height (mm)', color='black')

    ax2 = ax1.twinx()
    ax2.plot(filtered_offsets, mass, label='Mass (kg)', color='purple', linestyle='dashed')
    ax2.set_ylabel('Mass (kg)', color='purple')

    ax1.set_title('Boat Stability and Mass Data (40mm wider & 75mm shorter)')
    ax1.legend(loc='upper left')
    ax2.legend(loc='upper right')
    ax1.grid(True)

    plt.show()

# Example usage with cutoff
plot_stability_data(Calculated_40_75, (metacentric_heights_results, moments_of_area_results, cg_results), cutoff_offset=50)
