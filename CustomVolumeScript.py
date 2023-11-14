import adsk.core, adsk.fusion, traceback

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
        for offset in range(10, 226, 10):
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
