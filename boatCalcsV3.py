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




'''Offset, Volumes (mm³), Lengths (mm), Widths (mm), Center of Gravity (mm):
Offset: 10 mm, Volume: 965356.3220419986 mm³, Length: 689.7726398685952 mm, Width: 224.77680640344795 mm, Center of Gravity (0.0006387627264570206, 3.077821067121178e-05, 5.48051739151866) (mm):
Offset: 20 mm, Volume: 2317068.5378752714 mm³, Length: 725.9259259319572 mm, Width: 260.97407050926057 mm, Center of Gravity (0.004088327136260728, 0.004969922448535002, 11.136375525002862) (mm):
Offset: 30 mm, Volume: 3909490.829825512 mm³, Length: 755.5555555604138 mm, Width: 284.50813147265524 mm, Center of Gravity (0.00043010609549754, -0.001639618305086521, 16.827064967534042) (mm):
Offset: 40 mm, Volume: 5695029.668585284 mm³, Length: 785.1851851888704 mm, Width: 305.00000000378697 mm, Center of Gravity (0.00024764865899862265, -0.0007769931023710042, 22.55236144961195) (mm):
Offset: 50 mm, Volume: 7674351.555923939 mm³, Length: 814.8148148173268 mm, Width: 325.0000000032855 mm, Center of Gravity (0.0001450776820827171, -0.0009208216767347954, 28.364066611092404) (mm):
Offset: 60 mm, Volume: 9857594.995454479 mm³, Length: 844.4444444457836 mm, Width: 345.0000000024301 mm, Center of Gravity (7.912888208783109e-05, -0.0005855415024624833, 34.280887561711744) (mm):
Offset: 70 mm, Volume: 12253627.883523002 mm³, Length: 874.0740740742401 mm, Width: 365.0000000014551 mm, Center of Gravity (3.3939032225704935e-05, -0.0006342771703971017, 40.30251079397971) (mm):
Offset: 80 mm, Volume: 14872065.670262795 mm³, Length: 903.7037037026967 mm, Width: 385.0000000005909 mm, Center of Gravity (2.155658670455196e-06, -0.0003362320584343603, 46.42417373777203) (mm):
Offset: 90 mm, Volume: 17721794.233459614 mm³, Length: 933.3333333311532 mm, Width: 405.0000000000182 mm, Center of Gravity (0.00011365279657671312, -0.0009185567417252516, 52.638912079919606) (mm):
Offset: 100 mm, Volume: 20804429.657454524 mm³, Length: 958.1672533291883 mm, Width: 423.206400854453 mm, Center of Gravity (0.00039867177082264504, -0.0004125597588849894, 58.92432391372079) (mm):
Offset: 110 mm, Volume: 24073645.383786466 mm³, Length: 975.5961873600427 mm, Width: 436.071065000763 mm, Center of Gravity (0.00038145281871349024, 0.0059624038482114505, 65.1869693126713) (mm):
Offset: 120 mm, Volume: 27475492.46326336 mm³, Length: 987.6938930415545 mm, Width: 444.4790240969344 mm, Center of Gravity (0.0001350693344126519, 0.0003466925833561785, 71.35760630712895) (mm):
Offset: 130 mm, Volume: 30963051.345221777 mm³, Length: 995.4528167686499 mm, Width: 448.9959570815114 mm, Center of Gravity (0.0032425363974919463, 0.015572497407078667, 77.40057316210701) (mm):
Offset: 140 mm, Volume: 34488102.01405947 mm³, Length: 999.3798396948256 mm, Width: 450.0006894684381 mm, Center of Gravity (-0.0003106472069847044, -0.00048658842808968483, 83.28934346353662) (mm):
Offset: 150 mm, Volume: 38022307.229645096 mm³, Length: 1000.0014605249592 mm, Width: 450.0006894684381 mm, Center of Gravity (-0.0003257551725432961, -0.00036139272002078116, 89.02541992452329) (mm):
Offset: 160 mm, Volume: 41556600.25863865 mm³, Length: 1000.0014605249592 mm, Width: 450.0006894684381 mm, Center of Gravity (-0.000679928680333286, -0.002069563333996771, 94.63641337288863) (mm):
Offset: 170 mm, Volume: 45091002.228305504 mm³, Length: 1000.0014605249592 mm, Width: 450.0006894684381 mm, Center of Gravity (-0.0007487224764984868, -0.0028928900631447485, 100.15176526798572) (mm):
Offset: 180 mm, Volume: 48625430.93226241 mm³, Length: 1000.0014605249592 mm, Width: 450.0006894684381 mm, Center of Gravity (-0.000804439289618375, -0.0031005720807580905, 105.5922306282651) (mm):
Offset: 190 mm, Volume: 52159814.07846379 mm³, Length: 1000.0014605249592 mm, Width: 450.0006894684381 mm, Center of Gravity (-0.0009105238143031213, -0.0036316307434530737, 110.9729530615442) (mm):
Offset: 200 mm, Volume: 55693554.378497094 mm³, Length: 1000.0014605249592 mm, Width: 450.0006894684381 mm, Center of Gravity (-0.0002642364134508579, 0.0006548490801699947, 116.30469996437358) (mm):
Offset: 210 mm, Volume: 59227818.694045454 mm³, Length: 1000.0014605249592 mm, Width: 450.0006894684381 mm, Center of Gravity (-0.00026145922869980076, 0.0007887418309127092, 121.59735873171272) (mm):
Offset: 220 mm, Volume: 62762083.00953537 mm³, Length: 1000.0014605249592 mm, Width: 450.0006894684381 mm, Center of Gravity (-0.0002589948383480145, 0.0009075555833376603, 126.85705704309791) (mm):

'''

# Data representing offset heights from the bottom, volume, dimensions, and center of gravity for each offset
#40mm wider and 75mm shorter 
Offset_and_Volumes_40_75={
    # Each entry in the dictionary represents data for a specific offset from the waterline
    # The values are: volume in mm³, length in mm, width in mm, and center of gravity coordinates in mm (x, y, z)
    # Example: 10: (volume_mm3, y_length_mm, x_width_mm, cg_x_mm, cg_y_mm, cg_z_mm)
    10: (965356.3220419986, 689.7726398685952, 224.77680640344795, 0.0006387627264570206, 3.077821067121178e-05, 5.48051739151866), 
    20: (2317068.5378752714, 725.9259259319572, 260.97407050926057, 0.004088327136260728, 0.004969922448535002, 11.136375525002862), 
    30: (3909490.829825512, 755.5555555604138, 284.50813147265524, 0.00043010609549754, -0.001639618305086521, 16.827064967534042), 
    40: (5695029.668585284, 785.1851851888704, 305.00000000378697, 0.00024764865899862265, -0.0007769931023710042, 22.55236144961195), 
    50: (7674351.555923939, 814.8148148173268, 325.0000000032855, 0.0001450776820827171, -0.0009208216767347954, 28.364066611092404), 
    60: (9857594.995454479, 844.4444444457836, 345.0000000024301, 7.912888208783109e-05, -0.0005855415024624833, 34.280887561711744), 
    70: (12253627.883523002, 874.0740740742401, 365.0000000014551, 3.3939032225704935e-05, -0.0006342771703971017, 40.30251079397971), 
    80: (14872065.670262795, 903.7037037026967, 385.0000000005909, 2.155658670455196e-06, -0.0003362320584343603, 46.42417373777203), 
    90: (17721794.233459614, 933.3333333311532, 405.0000000000182, 0.00011365279657671312, -0.0009185567417252516, 52.638912079919606), 
    100: (20804429.657454524, 958.1672533291883, 423.206400854453, 0.00039867177082264504, -0.0004125597588849894, 58.92432391372079), 
    110: (24073645.383786466, 975.5961873600427, 436.071065000763, 0.00038145281871349024, 0.0059624038482114505, 65.1869693126713), 
    120: (27475492.46326336, 987.6938930415545, 444.4790240969344, 0.0001350693344126519, 0.0003466925833561785, 71.35760630712895), 
    130: (30963051.345221777, 995.4528167686499, 448.9959570815114, 0.0032425363974919463, 0.015572497407078667, 77.40057316210701), 
    140: (34488102.01405947, 999.3798396948256, 450.0006894684381, -0.0003106472069847044, -0.00048658842808968483, 83.28934346353662), 
    150: (38022307.229645096, 1000.0014605249592, 450.0006894684381, -0.0003257551725432961, -0.00036139272002078116, 89.02541992452329), 
    160: (41556600.25863865, 1000.0014605249592, 450.0006894684381, -0.000679928680333286, -0.002069563333996771, 94.63641337288863), 
    170: (45091002.228305504, 1000.0014605249592, 450.0006894684381, -0.0007487224764984868, -0.0028928900631447485, 100.15176526798572), 
    180: (48625430.93226241, 1000.0014605249592, 450.0006894684381, -0.000804439289618375, -0.0031005720807580905, 105.5922306282651), 
    190: (52159814.07846379, 1000.0014605249592, 450.0006894684381, -0.0009105238143031213, -0.0036316307434530737, 110.9729530615442), 
    200: (55693554.378497094, 1000.0014605249592, 450.0006894684381, -0.0002642364134508579, 0.0006548490801699947, 116.30469996437358), 
    210: (59227818.694045454, 1000.0014605249592, 450.0006894684381, -0.00026145922869980076, 0.0007887418309127092, 121.59735873171272), 
    220: (62762083.00953537, 1000.0014605249592, 450.0006894684381, -0.0002589948383480145, 0.0009075555833376603, 126.85705704309791)
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
def calculate_center_of_gravity(components, origin, hull_cg):
    # Update the hull component with the provided center of gravity coordinates
    for i, (name, mass, _) in enumerate(components):
        if name == "Hull":
            components[i] = (name, mass, hull_cg)

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



# Function to calculate the metacentric heights
def calculate_metacentric_heights(displacement, center_of_gravity, beam, length):
    # Calculate transverse and longitudinal metacentric radii
    BMt = (beam ** 2) / (12 * displacement)
    BMl = (length ** 2) / (12 * displacement)

    # Calculate metacentric heights
    GMt = BMt - center_of_gravity
    GMl = BMl - center_of_gravity
    return GMt, GMl



# Function to calculate metacentric heights for all offsets
def calculate_all_metacentric_heights(offset_volumes_dict):
    metacentric_heights_dict = {}

    # Loop through each offset in the dictionary
    for offset, (volume_mm3, mass_kg, y_length_mm, x_width_mm, cg_x_mm, cg_y_mm, cg_z_mm) in offset_volumes_dict.items():
        # Define the components of the boat including their mass and location
        components = [
            ("T200_1", 0.156, (0, 230.6, -49)),
            ("T200_2", 0.156, (0, -250, -49)),
            ("Electronics_Enclosure", 6.3, (1, 2.8, 88)),
            ("Velodyne_Puck", 1.65, (0.58, 436, 267)),
            ("SteroLabs_Camera", 0.23, (-2, 432.5, 319)),
            ("water_pump", 2, (0, 0, 100))
            # Add other components here
        ]

        # Use the center of gravity from the offset data
        hull_cg = (cg_x_mm, cg_y_mm, cg_z_mm)
        center_of_gravity = calculate_center_of_gravity(components, (0, 0, 0), hull_cg)

        # Convert displacement volume from mm³ to m³
        displacement_m3 = volume_mm3 * 1e-9
        # Convert beam and length from mm to m 
        beam_m = x_width_mm / 1000
        length_m = y_length_mm / 1000

        # Calculate metacentric heights for the current offset
        GMt, GMl = calculate_metacentric_heights(displacement_m3, center_of_gravity[2] / 1000, beam_m, length_m)

        # Store the results in the dictionary
        metacentric_heights_dict[offset] = (GMt, GMl)

    return metacentric_heights_dict


# Calculate metacentric heights for all offsets
metacentric_heights_results = calculate_all_metacentric_heights(Calculated_40_75)

# Print the results
for offset, (GMt, GMl) in metacentric_heights_results.items():
    print(f"Offset: {offset} mm, GMt: {GMt} m, GMl: {GMl} m")