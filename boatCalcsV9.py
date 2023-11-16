import pprint
import math
import matplotlib.pyplot as plt
from hullData import Offset_and_Volumes_40_75, Offset_and_Volumes_80_75, Offset_and_Volumes_80_105, Offset_and_Volumes_40_105
import re

"""
11/15/23

I think that I am calulating the CG for each of the different hulls wrong, because if 
I change the overall height of the hull then the CG should change but I currently done 
represent that in the program below. 
"""


# Function to calculate and append the mass displaced for each volume in the dictionary
def calculate_and_append_mass_displaced(volume_dict):
    density_of_water = 1000  # kg/m^3
    conversion_factor = 1e-9  # Convert mm^3 to m^3

     # Loop through each entry in the dictionary
    for offset, data in volume_dict.items():
        # Unpack the data. Adjust the number of variables here to match your actual data structure
        volume, y_length, x_width, cg_x, cg_y, cg_z, *rest = data
        mass_displaced = volume * conversion_factor * density_of_water
        # Update the dictionary entry with the calculated mass
        volume_dict[offset] = (volume, mass_displaced, y_length, x_width, cg_x, cg_y, cg_z, *rest)

    return volume_dict



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



# Function to calculate metacentric heights for all offsets
def calculate_all_metacentric_heights(offset_volumes_dict, dataset_name):
    metacentric_heights_dict = {}
    moments_of_area_dict = {}
    cg_dict = {}  # Dictionary to store the center of gravity

    # Extract added_width and subtracted_height from dataset name
    match = re.search(r"Hull (\d+)mm wider and (\d+)mm shorter", dataset_name)
    if match:
        added_width = int(match.group(1))
        subtracted_height = int(match.group(2))
    else:
        added_width = 40  # Default value if not found
        subtracted_height = 75  # Default value if not found

    #print(added_width)
    #print(subtracted_height)

    #These are from fusion360 properities of the whole boat hull body.
    hull_cg_40_75 = ("Hull", 5, (0, 0, 128.347804)) #mm
    hull_cg_80_75 = ("Hull", 5, (0, 0, 127.446557)) #mm
    hull_cg_80_105 = ("Hull", 5, (0, 0, 110.454906)) #mm
    hull_cg_40_105 = ("Hull", 5, (0, 0, 111.23766)) #mm

    # Loop through each offset in the dictionary
    for offset, data in offset_volumes_dict.items():
        # Adjust the number of variables here to match your actual data structure
        volume_mm3, mass_kg, y_length_mm, x_width_mm, cb_x_mm, cb_y_mm, cb_z_mm, *rest = data


        # Define the components of the boat including their mass and location
        #This will set the correct hull CG for the different hulls 
        if added_width == 40 and subtracted_height == 75:
            components = [
            ("T200_1", 0.156, (0, 230.6, -49)),
            ("T200_2", 0.156, (0, -250, -49)),
            ("Electronics_Enclosure", 5, (0, 0, 116.5)),
            ("Velodyne_Puck", 1.65, (0.58, 436, 267)),
            ("SteroLabs_Camera", 0.23, (-2, 432.5, 319)),
            ("water_pump", 2, (0, 0, 80)),
            ("Hull", 5, (0, 0, 128.347804))
            ] 
        elif added_width == 80 and subtracted_height == 75:
            components = [
            ("T200_1", 0.156, (0, 230.6, -49)),
            ("T200_2", 0.156, (0, -250, -49)),
            ("Electronics_Enclosure", 5, (0, 0, 116.5)),
            ("Velodyne_Puck", 1.65, (0.58, 436, 267)),
            ("SteroLabs_Camera", 0.23, (-2, 432.5, 319)),
            ("water_pump", 2, (0, 0, 80)),
            ("Hull", 5, (0, 0, 127.446557))
            ] 
        elif added_width == 80 and subtracted_height == 105:
            components = [
            ("T200_1", 0.156, (0, 230.6, -49)),
            ("T200_2", 0.156, (0, -250, -49)),
            ("Electronics_Enclosure", 5, (0, 0, 88)),
            ("Velodyne_Puck", 1.65, (0.58, 436, 237)),
            ("SteroLabs_Camera", 0.23, (-2, 432.5, 289)),
            ("water_pump", 2, (0, 0, 80)),
            ("Hull", 5, (0, 0, 110.454906))
            ]
        elif added_width == 40 and subtracted_height == 105:
            components = [
            ("T200_1", 0.156, (0, 230.6, -49)),
            ("T200_2", 0.156, (0, -250, -49)),
            ("Electronics_Enclosure", 5, (0, 0, 88)),
            ("Velodyne_Puck", 1.65, (0.58, 436, 237)),
            ("SteroLabs_Camera", 0.23, (-2, 432.5, 289)),
            ("water_pump", 2, (0, 0, 80)),
            ("Hull", 5, (0, 0, 111.23766))
            ]
        else:
            print("ERROR: Bad Hull CG Data")

        #print(components)
        origin = (0, 0, 0)
        constant_cg = calculate_center_of_gravity(components, origin)
        # Store the CG in the dictionary
        cg_dict[offset] = constant_cg

        # Center of buoyancy Z coordinate in meters
        center_of_buoyancy_z_m = cb_z_mm / 1000

        

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


def calculate_roll_times(offset_volumes_dict, metacentric_heights_dict):
    """
    Calculate the roll time for each offset.

    :param offset_volumes_dict: Dictionary with volume and dimension data.
    :param metacentric_heights_dict: Dictionary with calculated metacentric heights.
    :return: Dictionary of roll times for each offset.
    """
    roll_times_dict = {}

    for offset, metacentric_data in metacentric_heights_dict.items():
        GMt_mm = metacentric_data[0]  # Transverse Metacentric Height in millimeters
        if GMt_mm <= 0:
            continue  # Skip if GMt is zero or negative

        GMt_m = GMt_mm / 1000  # Convert GMt from mm to meters
        beam_m = offset_volumes_dict[offset][3] / 1000  # Convert width from mm to meters
        g = 9.81  # Acceleration due to gravity in m/s^2

        # Calculate the roll time using the formula
        T = 2 * math.pi * math.sqrt(beam_m / (g * GMt_m))
        # Calculate the stability by T / Beam 
        stability = T / beam_m
        roll_times_dict[offset] = (T, stability)


    return roll_times_dict



def print_all_data(offset_volumes_dict, combined_results):
    metacentric_heights_dict, moments_of_area_dict, cg_dict, roll_times_dict = combined_results

    for offset in offset_volumes_dict:
        # Retrieve data from dictionaries
        volume_data = offset_volumes_dict[offset]
        metacentric_data = metacentric_heights_dict.get(offset, (None, None))
        moment_data = moments_of_area_dict.get(offset, {})
        cg_data = cg_dict.get(offset, (0, 0, 0))  # Default to (0,0,0) if not present
        roll_time_data = roll_times_dict.get(offset, (None, None))  # Default to (None, None) if not present

        # Unpack data
        volume_mm3, mass_kg, y_length_mm, x_width_mm, cb_x_mm, cb_y_mm, cb_z_mm = volume_data
        GMt_mm, GMl_mm = metacentric_data
        I_ellipse = moment_data.get('I_ellipse', 0)
        I_rectangle = moment_data.get('I_rectangle', 0)
        I_total = moment_data.get('I_total', 0)
        roll_time, stability = roll_time_data

        # Print data for each offset
        print(f"Offset: {offset} mm")
        print(f"Volume: {volume_mm3:.2f} mm³, Mass: {mass_kg:.2f} kg")
        print(f"Dimensions (Length x Width): {y_length_mm:.2f} mm x {x_width_mm:.2f} mm")
        print(f"Center of Buoyancy (X, Y, Z): {cb_x_mm:.2f} mm, {cb_y_mm:.2f} mm, {cb_z_mm:.2f} mm")
        print(f"Center of Gravity (X, Y, Z): {cg_data[0]:.2f} mm, {cg_data[1]:.2f} mm, {cg_data[2]:.2f} mm")
        print(f"Metacentric Heights - Transverse Metacentric Height (GMt): {GMt_mm:.2f} mm, Longitudinal Metacentric Height (GMl): {GMl_mm:.2f} mm")
        print(f"Moments of Area - Ellipse: {I_ellipse:.2f} mm⁴, Rectangle: {I_rectangle:.2f} mm⁴, Total: {I_total:.2f} mm⁴")
        if roll_time and stability:
            print(f"Roll Time: {roll_time:.2f} seconds, Stability: {stability:.2f}")
        print("-" * 50)

# Example usage
#print_all_data(Calculated_40_75, (metacentric_heights_results, moments_of_area_results, cg_results, roll_times_results))

def print_dictionary_structure(dict_data):
    for key, value in dict_data.items():
        print(f"Key: {key}, Number of Elements: {len(value)}, Value: {value}")


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


#plot_stability_data(Calculated_40_75, (metacentric_heights_results, moments_of_area_results, cg_results))



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
#plot_stability_data(Calculated_40_75, (metacentric_heights_results, moments_of_area_results, cg_results), cutoff_offset=50)


def plot_combined_roll_times(hulls_data):
    """
    Plots the roll times for all hulls on a single graph.

    :param hulls_data: Dictionary with keys as hull names and values as hull data dictionaries.
    """
    plt.figure(figsize=(12, 8))

    # Iterate over each hull's data and plot
    for hull_name, hull_data in hulls_data.items():
        # Perform necessary calculations for each hull
        calculated_data = calculate_and_append_mass_displaced(hull_data)
        metacentric_heights_results, _, _ = calculate_all_metacentric_heights(calculated_data)
        roll_times_results = calculate_roll_times(calculated_data, metacentric_heights_results)

        # Prepare data for plotting
        offsets = list(roll_times_results.keys())
        roll_times = [rt[0] for rt in roll_times_results.values()]  # Assuming first element in tuple is roll time
        
        # Plotting each hull's roll time
        plt.plot(offsets, roll_times, label=hull_name)

    # Customizing the plot
    plt.xlabel('Offset from Bottom (mm)')
    plt.ylabel('Roll Time (seconds)')
    plt.title('Roll Times for Different Hulls')
    plt.legend()
    plt.grid(True)

    plt.show()


def plot_stability_data(offset_volumes_dict, combined_results, subplot_position, cutoff_offset=0):
    metacentric_heights_dict, _, cg_dict = combined_results

    # Filter data based on the cutoff_offset
    filtered_offsets = [offset for offset in offset_volumes_dict.keys() if offset >= cutoff_offset]
    
    center_of_gravity_z = [cg_dict[offset][2] for offset in filtered_offsets]
    center_of_buoyancy_z = [offset_volumes_dict[offset][6] for offset in filtered_offsets]
    GMt = [metacentric_heights_dict[offset][0] for offset in filtered_offsets] 
    mass = [offset_volumes_dict[offset][1] for offset in filtered_offsets]

    ax1 = plt.subplot(subplot_position)
    ax1.plot(filtered_offsets, center_of_gravity_z, label='Center of Gravity (Z)', color='blue')
    ax1.plot(filtered_offsets, center_of_buoyancy_z, label='Center of Buoyancy (Z)', color='green')
    ax1.plot(filtered_offsets, GMt, label='Transverse Metacentric Height (GMt)', color='red')
    ax1.set_xlabel('Offset from bottom (mm)')
    ax1.set_ylabel('Height (mm)', color='black')

    ax2 = ax1.twinx()
    ax2.plot(filtered_offsets, mass, label='Mass (kg)', color='purple', linestyle='dashed')
    ax2.set_ylabel('Mass (kg)', color='purple')

    ax1.legend(loc='upper left')
    ax2.legend(loc='upper right')
    ax1.grid(True)


#only displays 50mm ofset to 150mm 
def plot_combined_stability_data(hulls_processed_data):
    plt.figure(figsize=(12, 8))

    # Define colors for different hulls
    colors = ['blue', 'green', 'red', 'purple']
    hull_names = list(hulls_processed_data.keys())

    for hull_name, color in zip(hull_names, colors):
        data = hulls_processed_data[hull_name]
        offset_volumes_dict = data["calculated_data"]
        metacentric_heights_dict = data["metacentric_heights"]
        cg_dict = data["cg"]

        # Filter offsets between 50mm and 150mm
        filtered_offsets = [offset for offset in offset_volumes_dict.keys() if 50 <= offset <= 150]
        center_of_gravity_z = [cg_dict[offset][2] for offset in filtered_offsets]
        center_of_buoyancy_z = [offset_volumes_dict[offset][6] for offset in filtered_offsets]
        GMt = [metacentric_heights_dict[offset][0] for offset in filtered_offsets] 

        plt.plot(filtered_offsets, center_of_gravity_z, label=f'{hull_name} - CG (Z)', color=color, linestyle='dotted')
        plt.plot(filtered_offsets, center_of_buoyancy_z, label=f'{hull_name} - CB (Z)', color=color, linestyle='dashdot')
        plt.plot(filtered_offsets, GMt, label=f'{hull_name} - GMt', color=color, linestyle='solid')

    plt.xlabel('Offset from bottom (mm)')
    plt.ylabel('Height (mm)')
    plt.title('Combined Boat Stability Data (50mm to 150mm Offset)')
    plt.legend(loc='upper left')
    plt.grid(True)
    plt.show()


def plot_combined_roll_times_and_stability(hulls_processed_data):
    plt.figure(figsize=(12, 8))

    # Define colors for different hulls
    colors = ['blue', 'green', 'red', 'purple']
    hull_names = list(hulls_processed_data.keys())

    for hull_name, color in zip(hull_names, colors):
        data = hulls_processed_data[hull_name]
        offset_volumes_dict = data["calculated_data"]
        roll_times_data = data["roll_times"]

        filtered_offsets = [offset for offset in offset_volumes_dict.keys() if 50 <= offset <= 150]
        roll_times = [roll_times_data[offset][0] for offset in filtered_offsets if offset in roll_times_data]
        stability = [roll_times_data[offset][1] for offset in filtered_offsets if offset in roll_times_data]

        #plt.plot(filtered_offsets, roll_times, label=f'{hull_name} - Roll Time', color=color, linestyle='solid')
        plt.plot(filtered_offsets, stability, label=f'{hull_name} - Stability', color=color, linestyle='dashed')

    plt.xlabel('Offset from bottom (mm)')
    plt.ylabel('Roll Time and Stability')
    plt.title('Combined Roll Times and Stability Data (50mm to 150mm Offset)')
    plt.legend(loc='upper left')
    plt.grid(True)
    plt.show()



#plot_combined_roll_times(hulls)
hulls = {
        "Hull 40mm wider and 75mm shorter": Offset_and_Volumes_40_75,
        "Hull 80mm wider and 75mm shorter": Offset_and_Volumes_80_75,
        "Hull 80mm wider and 105mm shorter": Offset_and_Volumes_80_105,
        "Hull 40mm wider and 105mm shorter": Offset_and_Volumes_40_105
    }



def main():
    hulls_processed_data = {}

    for hull_name, hull_data in hulls.items():
        calculated_data = calculate_and_append_mass_displaced(hull_data)
        metacentric_heights_results, moments_of_area_results, cg_results = calculate_all_metacentric_heights(calculated_data, hull_name)
        roll_times_results = calculate_roll_times(calculated_data, metacentric_heights_results)

        hulls_processed_data[hull_name] = {
            "calculated_data": calculated_data,
            "metacentric_heights": metacentric_heights_results,
            "moments_of_area": moments_of_area_results,
            "cg": cg_results,
            "roll_times": roll_times_results
        }

    # Now, use the processed data for each hull
    for hull_name, results in hulls_processed_data.items():
        print(f"\nProcessing data for {hull_name}")
        print_all_data(results["calculated_data"], (results["metacentric_heights"], results["moments_of_area"], results["cg"], results["roll_times"]))
        # Additional plotting functions can be called here as needed


    plot_combined_stability_data(hulls_processed_data)
    plot_combined_roll_times_and_stability(hulls_processed_data)
    # Create a figure for the subplots
    plt.figure(figsize=(20, 12))  # Adjust size as needed

    subplot_positions = [221, 222, 223, 224]  # Positions in a 2x2 grid
    for (hull_name, results), subplot_pos in zip(hulls_processed_data.items(), subplot_positions):
        plot_stability_data(results["calculated_data"], 
                            (results["metacentric_heights"], results["moments_of_area"], results["cg"]),
                            subplot_pos)
        plt.subplot(subplot_pos).set_title(hull_name)

    plt.tight_layout()
    plt.show()

if __name__ == "__main__":
    main()