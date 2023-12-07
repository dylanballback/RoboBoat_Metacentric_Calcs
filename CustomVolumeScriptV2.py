import adsk.core, adsk.fusion, traceback
from CustomVolumeDatabase import insert_data  # Importing from your database script

# --------------------------------------------------------------------------------------------------

def create_offset_plane(plane_type, distance_mm):
    app = adsk.core.Application.get()
    ui = app.userInterface
    design = adsk.fusion.Design.cast(app.activeProduct)
    rootComp = design.rootComponent

    # Find the origin plane based on the input
    if plane_type.lower() == 'xy':
        originPlane = rootComp.xYConstructionPlane
    elif plane_type.lower() == 'xz':
        originPlane = rootComp.xZConstructionPlane
    elif plane_type.lower() == 'yz':
        originPlane = rootComp.yZConstructionPlane
    else:
        ui.messageBox('Invalid plane type. Please specify XY, XZ, or YZ.')
        return None

    # Create an offset construction plane
    constructionPlanes = rootComp.constructionPlanes
    planeInput = constructionPlanes.createInput()
    distanceValue = adsk.core.ValueInput.createByReal(distance_mm / 10.0)  # Convert mm to cm
    planeInput.setByOffset(originPlane, distanceValue)
    offsetPlane = constructionPlanes.add(planeInput)

    return offsetPlane

# Example usage:
# plane = create_offset_plane('XY', 50)  # Creates an offset plane from the XY plane at 50mm



# --------------------------------------------------------------------------------------------------



def get_center_of_gravity(body_name):
    app = adsk.core.Application.get()
    ui = app.userInterface
    design = adsk.fusion.Design.cast(app.activeProduct)
    rootComp = design.rootComponent

    # Find the body by name
    body = rootComp.bRepBodies.itemByName(body_name)
    if body is None: 
        ui.messageBox(f'Body "{body_name}" not found')
        return None

    # Get the physical properties of the body
    physicalProperties = body.getPhysicalProperties(adsk.fusion.CalculationAccuracy.HighCalculationAccuracy)
    centerOfGravity = physicalProperties.centerOfMass

    # Convert the center of gravity point to the unit system used in the design
    unitsMgr = design.fusionUnitsManager
    x = unitsMgr.convert(centerOfGravity.x, 'cm', unitsMgr.defaultLengthUnits)
    y = unitsMgr.convert(centerOfGravity.y, 'cm', unitsMgr.defaultLengthUnits)
    z = unitsMgr.convert(centerOfGravity.z, 'cm', unitsMgr.defaultLengthUnits)

    return adsk.core.Point3D.create(x, y, z)

# Example usage:
# cog = get_center_of_gravity('Body1')
# if cog:
#     app = adsk.core.Application.get()
#     ui = app.userInterface
#     ui.messageBox(f'Center of Gravity: X = {cog.x}, Y = {cog.y}, Z = {cog.z}')



# --------------------------------------------------------------------------------------------------



def move_body(body_name, axis, distance_mm):
    app = adsk.core.Application.get()
    design = adsk.fusion.Design.cast(app.activeProduct)
    rootComp = design.rootComponent

    # Find the body by name
    body = rootComp.bRepBodies.itemByName(body_name)
    if body is None:
        app.userInterface.messageBox(f'Body "{body_name}" not found')
        return False

    # Create a transform to move the body
    vector = adsk.core.Vector3D.create(0, 0, 0)
    if axis.upper() == 'X':
        vector.x = distance_mm / 10.0  # Converting mm to cm
    elif axis.upper() == 'Y':
        vector.y = distance_mm / 10.0  # Converting mm to cm
    elif axis.upper() == 'Z':
        vector.z = distance_mm / 10.0  # Converting mm to cm
    else:
        app.userInterface.messageBox('Invalid axis. Please specify X, Y, or Z.')
        return False

    transform = adsk.core.Matrix3D.create()
    transform.translation = vector

    # Move the body
    body.transform = transform
    return True

# Example usage:
# move_body('Body1', 'X', 50)



# --------------------------------------------------------------------------------------------------



def non_uniform_scale(body_name, axis, scale_factor):
    app = adsk.core.Application.get()
    ui = app.userInterface
    design = adsk.fusion.Design.cast(app.activeProduct)
    rootComp = design.rootComponent

    # Find the body by name
    body = rootComp.bRepBodies.itemByName(body_name)
    if body is None:
        ui.messageBox(f'Body "{body_name}" not found')
        return False

    if not (0 <= scale_factor <= 1):
        ui.messageBox('Scale factor must be between 0 and 1')
        return False

    # Define the scaling factors for each axis
    x_factor, y_factor, z_factor = 1.0, 1.0, 1.0
    if axis.lower() == 'x':
        x_factor = scale_factor
    elif axis.lower() == 'y':
        y_factor = scale_factor
    elif axis.lower() == 'z':
        z_factor = scale_factor
    else:
        ui.messageBox('Invalid axis. Please specify X, Y, or Z.')
        return False

    # Create a non-uniform scaling transform
    scalingMatrix = adsk.core.Matrix3D.create()
    scalingMatrix.setCell(0, 0, x_factor)
    scalingMatrix.setCell(1, 1, y_factor)
    scalingMatrix.setCell(2, 2, z_factor)

    # Apply the scaling transform to the body
    entities = adsk.core.ObjectCollection.create()
    entities.add(body)
    transform = adsk.core.Matrix3D.create()
    transform.translation = adsk.core.Vector3D.create(0, 0, 0)
    rootComp.features.moveFeatures.add(entities, transform, scalingMatrix, None, None, None)

    return True

# Example usage:
# non_uniform_scale('Body1', 'X', 0.5)



# --------------------------------------------------------------------------------------------------



def get_body_names_from_component(component):
    body_names = []
    for body in component.bRepBodies:
        body_names.append(body.name)
    return body_names


# Example usage
"""
def run(context):
    app = adsk.core.Application.get()
    ui = app.userInterface
    try:
        design = adsk.fusion.Design.cast(app.activeProduct)
        rootComp = design.rootComponent
        # Assuming you want to get the body names from the root component
        # Replace 'rootComp' with your specific component if necessary
        body_names = get_body_names_from_component(rootComp)
        ui.messageBox("Body names:\n" + "\n".join(body_names))
    except:
        if ui:
            ui.messageBox('Failed:\n{}'.format(traceback.format_exc()))

"""

# --------------------------------------------------------------------------------------------------



def split_hull_at_yz(app, rootComp, your_hull_body):
    # Create a YZ construction plane at the origin
    yzPlane = rootComp.yZConstructionPlane
    splitFeatures = rootComp.features.splitBodyFeatures
    splitInput = splitFeatures.createInput(your_hull_body, yzPlane, True)
    splitFeature = splitFeatures.add(splitInput)
    return splitFeature.bodies

# --------------------------------------------------------------------------------------------------


def manipulate_and_join_bodies(app, rootComp, left_body, right_body, increment):
    # Move the bodies
    transform_left = adsk.core.Matrix3D.create()
    transform_left.translation = adsk.core.Vector3D.create(-increment / 2, 0, 0)
    left_body.transform = transform_left

    transform_right = adsk.core.Matrix3D.create()
    transform_right.translation = adsk.core.Vector3D.create(increment / 2, 0, 0)
    right_body.transform = transform_right

    # Extrude logic to join the bodies
    # This will depend on the specific geometry of the hull and may require more complex handling.
    # ...


# --------------------------------------------------------------------------------------------------


def calculate_and_store_data(lowerBody, transform, offset):
    volume_mm3 = lowerBody.volume * 1000  # Convert cm³ to mm³

    # ... [Center of gravity calculation using the transform] ...

    y_length_mm = (lowerBody.boundingBox.maxPoint.y - lowerBody.boundingBox.minPoint.y) * 10
    x_width_mm = (lowerBody.boundingBox.maxPoint.x - lowerBody.boundingBox.minPoint.x) * 10

    # Store data in a dictionary or directly insert into the database
    data = {'offset': offset, 'volume': volume_mm3, 'y_length': y_length_mm, 'x_width': x_width_mm, 'cg_x': cg_x, 'cg_y': cg_y, 'cg_z': cg_z}
    # insert_data(...) or return data
    return data


# --------------------------------------------------------------------------------------------------

def create_offset_plane(rootComp, offset):
    constructionPlanes = rootComp.constructionPlanes
    planeInput = constructionPlanes.createInput()
    offsetValue = adsk.core.ValueInput.createByReal(offset / 10.0)  # Convert mm to cm
    xyPlane = rootComp.xYConstructionPlane
    planeInput.setByOffset(xyPlane, offsetValue)
    offsetPlane = constructionPlanes.add(planeInput)
    return offsetPlane

# --------------------------------------------------------------------------------------------------

def split_body(boatHullOcc, originalBody, offsetPlane):
    splitBodyFeats = boatHullOcc.component.features.splitBodyFeatures
    splitBodyInput = splitBodyFeats.createInput(originalBody, offsetPlane, True)
    splitBodyFeature = splitBodyFeats.add(splitBodyInput)
    return splitBodyFeature

# --------------------------------------------------------------------------------------------------

def calculate_volume_and_cg(lowerBody, transform):
    volume_mm3 = lowerBody.volume * 1000  # Convert cm³ to mm³
    physicalProperties = lowerBody.physicalProperties
    centerOfGravity = physicalProperties.centerOfMass

    cg_x, cg_y, cg_z = apply_transform_to_point(centerOfGravity, transform)
    return volume_mm3, cg_x, cg_y, cg_z

# --------------------------------------------------------------------------------------------------




# --------------------------------------------------------------------------------------------------

def perform_calculations_and_store_data(rootComp, boatHullOcc, originalBody):
    offset_volume_data = {}

    for offset in range(10, 196, 10):
        # Create an offset construction plane
        constructionPlanes = rootComp.constructionPlanes
        planeInput = constructionPlanes.createInput()
        offsetValue = adsk.core.ValueInput.createByReal(offset / 10.0)  # Convert mm to cm
        xyPlane = rootComp.xYConstructionPlane
        planeInput.setByOffset(xyPlane, offsetValue)
        offsetPlane = constructionPlanes.add(planeInput)

        # Create SplitBodyFeatureInput
        splitBodyFeats = boatHullOcc.component.features.splitBodyFeatures
        splitBodyInput = splitBodyFeats.createInput(originalBody, offsetPlane, True)
        splitBodyFeature = splitBodyFeats.add(splitBodyInput)

        # Continue if no new body is created
        if splitBodyFeature.bodies.count < 2:
            continue

        # Determine the lower body
        body1, body2 = splitBodyFeature.bodies.item(0), splitBodyFeature.bodies.item(1)
        centroid1 = get_centroid(body1)
        centroid2 = get_centroid(body2)
        lowerBody = body1 if centroid1.z < centroid2.z else body2

        # Calculate volume and center of gravity
        volume_mm3, cg_x, cg_y, cg_z = calculate_volume_cg_and_transform(lowerBody, boatHullOcc.transform)

        # Measure lengths and widths
        y_length_mm = get_length_in_y(lowerBody)
        x_width_mm = get_width_in_x(lowerBody)

        # Store the calculated data
        offset_volume_data[offset] = (volume_mm3, y_length_mm, x_width_mm, cg_x, cg_y, cg_z)

    return offset_volume_data

def get_centroid(body):
    return adsk.core.Point3D.create(
        (body.boundingBox.minPoint.x + body.boundingBox.maxPoint.x) / 2,
        (body.boundingBox.minPoint.y + body.boundingBox.maxPoint.y) / 2,
        (body.boundingBox.minPoint.z + body.boundingBox.maxPoint.z) / 2
    )

def calculate_volume_cg_and_transform(body, transform):
    volume_mm3 = body.volume * 1000  # Convert cm³ to mm³
    physicalProperties = body.physicalProperties
    centerOfGravity = physicalProperties.centerOfMass
    cg_x, cg_y, cg_z = apply_transform_to_point(centerOfGravity, transform)
    return volume_mm3, cg_x, cg_y, cg_z

def apply_transform_to_point(centerOfGravity, transform):
    # Manually apply the transformation matrix to the center of gravity point
    # Extract the elements of the transformation matrix
    a11 = transform.getCell(0, 0)
    a12 = transform.getCell(0, 1)
    a13 = transform.getCell(0, 2)
    a14 = transform.getCell(0, 3)
    a21 = transform.getCell(1, 0)
    a22 = transform.getCell(1, 1)
    a23 = transform.getCell(1, 2)
    a24 = transform.getCell(1, 3)
    a31 = transform.getCell(2, 0)
    a32 = transform.getCell(2, 1)
    a33 = transform.getCell(2, 2)
    a34 = transform.getCell(2, 3)

    # Apply transformation to the center of gravity
    cg_x = (a11 * centerOfGravity.x + a12 * centerOfGravity.y + a13 * centerOfGravity.z + a14) * 10
    cg_y = (a21 * centerOfGravity.x + a22 * centerOfGravity.y + a23 * centerOfGravity.z + a24) * 10
    cg_z = (a31 * centerOfGravity.x + a32 * centerOfGravity.y + a33 * centerOfGravity.z + a34) * 10

    return cg_x, cg_y, cg_z

def get_length_in_y(body):
    return (body.boundingBox.maxPoint.y - body.boundingBox.minPoint.y) * 10

def get_width_in_x(body):
    return (body.boundingBox.maxPoint.x - body.boundingBox.minPoint.x) * 10

def split_body(rootComp, body, plane_type, distance_mm):
    offset_plane = create_offset_plane(plane_type, distance_mm, rootComp)
    splitFeatures = rootComp.features.splitBodyFeatures
    splitInput = splitFeatures.createInput(body, offset_plane, True)
    return splitFeatures.add(splitInput)


def recombine_bodies(rootComp, bodies):
    combineFeatures = rootComp.features.combineFeatures
    toolBodies = adsk.core.ObjectCollection.create()
    for body in bodies:
        toolBodies.add(body)

    combineInput = combineFeatures.createInput(bodies.item(0), toolBodies)
    combineInput.operation = adsk.fusion.FeatureOperations.JoinFeatureOperation
    return combineFeatures.add(combineInput)


def manipulate_bodies(splitBodies, increment):
    left_body, right_body = splitBodies.bodies.item(0), splitBodies.bodies.item(1)

    # Moving left and right bodies
    move_body(left_body.name, 'X', -increment / 2)
    move_body(right_body.name, 'X', increment / 2)

    # Extrude logic here (implementation depends on the specific geometry)
    # You need to create a feature to extrude and join these bodies
# --------------------------------------------------------------------------------------------------




def run(context):
    try:
        app = adsk.core.Application.get()
        ui = app.userInterface
        design = adsk.fusion.Design.cast(app.activeProduct)
        rootComp = design.rootComponent

        # Find the "Boat Hull" occurrence.
        boatHullOcc = None
        for occ in rootComp.occurrences:
            if occ.component.name == "Boat Hull":
                boatHullOcc = occ
                break

        if not boatHullOcc:
            ui.messageBox('Boat Hull occurrence not found')
            return

        # Get the body named "Body1" from the "Boat Hull" component.
        originalBody = boatHullOcc.component.bRepBodies.itemByName("Body1")
        if not originalBody:
            ui.messageBox('Body1 not found in Boat Hull')
            return

        # Perform initial calculations on the original body
        # Placeholder for your calculation logic
        # initial_data = perform_initial_calculations(originalBody)
        # insert_data(**initial_data)  # Store initial data in the database

        # Splitting, Manipulating, and Recombining the Body
        total_increment = 200  # Total value for increments
        increment_value = 10  # Increment value

        for i in range(0, total_increment, increment_value):
            splitBodies = split_body(rootComp, originalBody, 'XY', i)
            manipulate_bodies(splitBodies, i)
            recombinedBody = recombine_bodies(rootComp, splitBodies.bodies)
            manipulated_data = perform_calculations(recombinedBody)
            insert_data(**manipulated_data)

    except:
        if ui:
            ui.messageBox('Failed:\n{}'.format(traceback.format_exc()))


"""


def run(context):
    ui = None
    try:
        # Get the application and user interface.
        app = adsk.core.Application.get()
        ui = app.userInterface

        # Get the active design.
        product = app.activeProduct
        design = adsk.fusion.Design.cast(product)
        rootComp = design.rootComponent

        # Find the "Boat Hull:1" occurrence. Adjust the name if necessary.
        boatHullOcc = None
        for occ in rootComp.occurrences:
            if occ.component.name == "Boat Hull":
                boatHullOcc = occ
                break

        if not boatHullOcc:
            ui.messageBox('Boat Hull:1 occurrence not found')
            return

        # Get "Body1" from the "Boat Hull:1" component.
        originalBody = boatHullOcc.component.bRepBodies.itemByName("Body1")
        if not originalBody:
            ui.messageBox('Body1 not found in Boat Hull:1')
            return

        # Split the hull at the YZ plane and get the bodies
        split_bodies = split_hull_at_yz(rootComp, originalBody)

        # Assuming split_bodies[0] is left and split_bodies[1] is right
        left_body = split_bodies[0]
        right_body = split_bodies[1]

        total_increment = 200  # Total value to increment
        increment_value = 10  # Individual increment value

        for i in range(0, total_increment, increment_value):
            manipulate_and_join_bodies(rootComp, left_body, right_body, increment_value)

            # Run calculations on modified hull
            run_calculations_and_store_data()

            # Reset or recreate the bodies for the next iteration if necessary

    except:
        if ui:
            ui.messageBox('Failed:\n{}'.format(traceback.format_exc()))
"""


            
"""
def run(context):
    ui = None
    try:
        # Get the application and user interface.
        app = adsk.core.Application.get()
        ui = app.userInterface

        # Get the active design.
        product = app.activeProduct
        design = adsk.fusion.Design.cast(product)
        rootComp = design.rootComponent

        # Find the "Boat Hull:1" occurrence. Adjust the name if necessary.
        boatHullOcc = None
        for occ in rootComp.occurrences:
            if occ.component.name == "Boat Hull":
                boatHullOcc = occ
                break

        if not boatHullOcc:
            ui.messageBox('Boat Hull:1 occurrence not found')
            return

        # Get "Body1" from the "Boat Hull:1" component.
        originalBody = boatHullOcc.component.bRepBodies.itemByName("Body1")
        if not originalBody:
            ui.messageBox('Body1 not found in Boat Hull:1')
            return

        # Dictionary to store offset heights and volumes.
        offset_volume_data = {}

        # Loop from 10mm to 225mm, stepping 10mm each time.
        for offset in range(10, 196, 10):
            # Create an offset construction plane from the root component's XY plane.
            constructionPlanes = rootComp.constructionPlanes
            planeInput = constructionPlanes.createInput()
            offsetValue = adsk.core.ValueInput.createByReal(offset / 10.0)  # Convert mm to cm if necessary
            xyPlane = rootComp.xYConstructionPlane
            planeInput.setByOffset(xyPlane, offsetValue)
            offsetPlane = constructionPlanes.add(planeInput)

            # Create SplitBodyFeatureInput using the construction plane as the splitting tool.
            splitBodyFeats = boatHullOcc.component.features.splitBodyFeatures
            splitBodyInput = splitBodyFeats.createInput(originalBody, offsetPlane, True)
            
            # Try to create the split body feature.
            splitBodyFeature = splitBodyFeats.add(splitBodyInput)

            # Check if new body was created, continue to the next iteration if not.
            if splitBodyFeature.bodies.count < 2:
                continue

            # Determine the lower body by comparing the centroid Z values.
            body1 = splitBodyFeature.bodies.item(0)
            body2 = splitBodyFeature.bodies.item(1)

            # Get centroid of each body
            centroid1 = adsk.core.Point3D.create(
                (body1.boundingBox.minPoint.x + body1.boundingBox.maxPoint.x) / 2,
                (body1.boundingBox.minPoint.y + body1.boundingBox.maxPoint.y) / 2,
                (body1.boundingBox.minPoint.z + body1.boundingBox.maxPoint.z) / 2
            )
            centroid2 = adsk.core.Point3D.create(
                (body2.boundingBox.minPoint.x + body2.boundingBox.maxPoint.x) / 2,
                (body2.boundingBox.minPoint.y + body2.boundingBox.maxPoint.y) / 2,
                (body2.boundingBox.minPoint.z + body2.boundingBox.maxPoint.z) / 2
            )

            # Assume lower body has a smaller Z value of centroid
            lowerBody = body1 if centroid1.z < centroid2.z else body2

            # Measure the volume of the lower body in cm³ and convert to mm³
            volume_mm3 = lowerBody.volume * 1000  # Convert cm³ to mm³

            # Get the physical properties of the lower body to find the center of gravity.
            physicalProperties = lowerBody.physicalProperties
            centerOfGravity = physicalProperties.centerOfMass

            # Get the transformation matrix of the boat hull occurrence
            transform = boatHullOcc.transform

            # Manually apply the transformation matrix to the center of gravity point
            # Extract the elements of the transformation matrix
            a11 = transform.getCell(0, 0)
            a12 = transform.getCell(0, 1)
            a13 = transform.getCell(0, 2)
            a14 = transform.getCell(0, 3)
            a21 = transform.getCell(1, 0)
            a22 = transform.getCell(1, 1)
            a23 = transform.getCell(1, 2)
            a24 = transform.getCell(1, 3)
            a31 = transform.getCell(2, 0)
            a32 = transform.getCell(2, 1)
            a33 = transform.getCell(2, 2)
            a34 = transform.getCell(2, 3)

            # Apply transformation to the center of gravity
            cg_x = (a11 * centerOfGravity.x + a12 * centerOfGravity.y + a13 * centerOfGravity.z + a14) * 10
            cg_y = (a21 * centerOfGravity.x + a22 * centerOfGravity.y + a23 * centerOfGravity.z + a24) * 10
            cg_z = (a31 * centerOfGravity.x + a32 * centerOfGravity.y + a33 * centerOfGravity.z + a34) * 10

            # Measure the longest length in Y direction and convert to mm
            y_length_mm = (lowerBody.boundingBox.maxPoint.y - lowerBody.boundingBox.minPoint.y) * 10

            # Measure the widest section in X direction and convert to mm
            x_width_mm = (lowerBody.boundingBox.maxPoint.x - lowerBody.boundingBox.minPoint.x) * 10

            # Store the offset height, volume, Y length, and X width in the dictionary.
            offset_volume_data[offset] = (volume_mm3, y_length_mm, x_width_mm, cg_x, cg_y, cg_z)

            # Combine the split bodies back together.
            combineFeatures = boatHullOcc.component.features.combineFeatures
            toolBodies = adsk.core.ObjectCollection.create()
            for body in splitBodyFeature.bodies:
                if body != lowerBody:
                    toolBodies.add(body)
            combineInput = combineFeatures.createInput(originalBody, toolBodies)
            combineInput.operation = adsk.fusion.FeatureOperations.JoinFeatureOperation
            combineFeatures.add(combineInput)

        # Create a string to display all volumes in mm³, lengths, and widths in mm
        result_message = "Offset, Volumes (mm³), Lengths (mm), Widths (mm), Center of Gravity (mm):\n"
        for offset, (volume_mm3, y_length_mm, x_width_mm, cg_x, cg_y, cg_z) in offset_volume_data.items():
            result_message += f"Offset: {offset} mm, Volume: {volume_mm3} mm³, Length: {y_length_mm} mm, Width: {x_width_mm} mm, Center of Gravity {cg_x, cg_y, cg_z} (mm):\n"

        # Display the result in a single message box
        ui.messageBox(result_message)

         # Convert the offset_volume_data to a string representation of a Python dictionary
        result_message = str(offset_volume_data)

        # Display the result in a single message box
        ui.messageBox(result_message)
       
    except:
        if ui:
            ui.messageBox('Failed:\n{}'.format(traceback.format_exc()))




"""