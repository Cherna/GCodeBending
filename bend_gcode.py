import numpy as np
import math
from scipy.interpolate import CubicSpline
import matplotlib.pyplot as plt
import re
from collections import namedtuple

Point2D = namedtuple('Point2D', 'x y')
GCodeLine = namedtuple('GCodeLine', 'x y z e f')


#################   USER INPUT PARAMETERS   #########################
GCODE_FOLDER = "gcode/"  # Folder where the GCode files are located
INPUT_FILE_NAME = GCODE_FOLDER + "Cylinder.gcode"
OUTPUT_FILE_NAME = GCODE_FOLDER + "BENT_" + INPUT_FILE_NAME.split("/")[-1]
LAYER_HEIGHT = 0.25 #Layer height of the sliced gcode
WARNING_ANGLE = 30 #Maximum Angle printable with your setup
MINIMUM_EXTRUSION = 0.0001  # Minimum extrusion so the E motor does not go BRrRrRrRrR

# Setup the arcing spline
SPLINE_X = [125, 125.035, 125.085, 125.5]
SPLINE_Z = [0, 3.5, 6, 16]
SPLINE_ANGLES_DEGREES = [0, 10, 4, 5]  # degrees

#################   USER INPUT PARAMETERS END  #########################

# Convert angles from degrees to radians
SPLINE_ANGLES = [np.radians(a) for a in SPLINE_ANGLES_DEGREES]

# Define boundary conditions using slope (tan(angle)) at the start and end
bc_start = (1, np.tan(SPLINE_ANGLES[0]))
bc_end   = (1, np.tan(SPLINE_ANGLES[-1]))

# Define spline with boundary conditions
SPLINE = CubicSpline(SPLINE_Z, SPLINE_X, bc_type=(bc_start, bc_end))

DISCRETIZATION_LENGTH = 0.01  # Discretization length for the spline length lookup table

# Optional: create a placeholder spline length lookup table
SplineLookupTable = [0.0]

# Plotting the spline in 2D
xs = np.linspace(0, SPLINE_Z[-1], 200)

# --- 3D Visualization of the spline (XZ plane, Y=0) ---
fig3d = plt.figure(figsize=(8, 6))
ax3d = fig3d.add_subplot(111, projection='3d')
ax3d.plot(SPLINE(xs), np.zeros_like(xs), xs, label='Spline (XZ)')
ax3d.scatter(SPLINE_X, np.zeros_like(SPLINE_X), SPLINE_Z, color='red', label='Control Points')

# Calculate bounds and apply even range to both axes
x_min = min(SPLINE_X)
x_max = max(SPLINE_X)
z_min = min(SPLINE_Z)
z_max = max(SPLINE_Z)

x_center = (x_min + x_max) / 2
z_center = (z_min + z_max) / 2

half_range = max((x_max - x_min), (z_max - z_min)) / 2 + 1  # add padding

ax3d.set_xlim(x_center - half_range, x_center + half_range)
ax3d.set_ylim(-half_range, half_range)  # Y axis (centered at 0)
ax3d.set_zlim(z_center - half_range, z_center + half_range)

# Set equal aspect ratio
ax3d.set_box_aspect([np.ptp([x_center - half_range, x_center + half_range]),
                     np.ptp([-half_range, half_range]),
                     np.ptp([z_center - half_range, z_center + half_range])])

ax3d.set_xlabel('X')
ax3d.set_ylabel('Y')
ax3d.set_zlabel('Z')
ax3d.legend()
plt.show()


def getNormalPoint(currentPoint: Point2D, derivative: float, distance: float) -> Point2D: #claculates the normal of a point on the spline
    angle = np.arctan(derivative) + math.pi /2
    return Point2D(currentPoint.x + distance * np.cos(angle), currentPoint.y + distance * np.sin(angle))

def parseGCode(currentLine: str) -> GCodeLine: #parse a G-Code line
    thisLine = re.compile(r'(?i)^[gG][0-3](?:\s+x(?P<x>-?[0-9.]{1,15})|\s+y(?P<y>-?[0-9.]{1,15})|\s+z(?P<z>-?[0-9.]{1,15})|\s+e(?P<e>-?[0-9.]{1,15})|\s+f(?P<f>-?[0-9.]{1,15}))*')
    lineEntries = thisLine.match(currentLine)
    if lineEntries:
        return GCodeLine(lineEntries.group('x'), lineEntries.group('y'), lineEntries.group('z'), lineEntries.group('e'), lineEntries.group('f'))

def writeLine(G, X, Y, Z, F = None, E = None): #write a line to the output file
    outputSting = "G" + str(int(G)) + " X" + str(round(X,5)) + " Y" + str(round(Y,5)) + " Z" + str(round(Z,3))
    if E is not None:
        outputSting = outputSting + " E" + str(round(float(E),5))
    if F is not None:
        outputSting = outputSting + " F" + str(int(float(F)))
    outputFile.write(outputSting + "\n")

def onSplineLength(Zheight) -> float: #calculates a new z height if the spline is followed
    for i in range(len(SplineLookupTable)):
        height = SplineLookupTable[i]
        if height >= Zheight:
            return i * DISCRETIZATION_LENGTH
    print("Error! Spline not defined high enough!")

def createSplineLookupTable():
    heightSteps = np.arange(DISCRETIZATION_LENGTH, SPLINE_Z[-1], DISCRETIZATION_LENGTH)
    for i in range(len(heightSteps)):
        height = heightSteps[i]
        SplineLookupTable.append(SplineLookupTable[i] + np.sqrt((SPLINE(height)-SPLINE(height-DISCRETIZATION_LENGTH))**2 + DISCRETIZATION_LENGTH**2))

lastPosition = Point2D(0, 0)
currentZ = 0.0
lastZ = 0.0
currentLayer = 0
relativeMode = False
createSplineLookupTable()

with open(INPUT_FILE_NAME, "r") as gcodeFile, open(OUTPUT_FILE_NAME, "w+") as outputFile:
    for currentLine in gcodeFile:
        if currentLine[0] == ";":  # If it's a comment line
            outputFile.write(currentLine)
            continue
        if currentLine.find("G91 ") != -1:  # Filter relative commands (skip G91)
            continue  # Skip this line, no need to write it
        if currentLine.find("G90 ") != -1:  # Set absolute mode (skip G90)
            continue  # Skip this line, no need to write it
        if relativeMode:  # If in relative mode don't do anything
            outputFile.write(currentLine)
            continue

        currentLineCommands = parseGCode(currentLine)
        if currentLineCommands is not None:  # If current command is valid G-code
            if currentLineCommands.z is not None:  # If there is a Z height in the command
                currentZ = float(currentLineCommands.z)

            if currentLineCommands.x is None or currentLineCommands.y is None:  # If no X/Y movement
                if currentLineCommands.z is not None:  # Only Z movement (e.g., Z-hop)
                    outputFile.write("G1 ")
                    if currentLineCommands.f is not None:
                        outputFile.write(" F" + str(currentLineCommands.f))
                    outputFile.write("\n")
                    lastZ = currentZ
                    continue
                outputFile.write(currentLine)
                continue

            # Parse X and Y values
            x = float(currentLineCommands.x)
            y = float(currentLineCommands.y)
            currentPosition = Point2D(x, y)

            # Calculate midpoint
            midpointX = lastPosition.x + (currentPosition.x - lastPosition.x) / 2  # Look for midpoint
            distToSpline = midpointX - SPLINE_X[0]

            # Correct the Z-height if the spline gets followed
            correctedZHeight = onSplineLength(currentZ)

            angleSplineThisLayer = np.arctan(SPLINE(correctedZHeight, 1))  # Inclination angle this layer
            angleLastLayer = np.arctan(SPLINE(correctedZHeight - LAYER_HEIGHT, 1))  # Inclination angle previous layer
            heightDifference = np.sin(angleSplineThisLayer - angleLastLayer) * distToSpline * -1  # Layer height difference

            transformedGCode = getNormalPoint(
                Point2D(correctedZHeight, SPLINE(correctedZHeight)),
                SPLINE(correctedZHeight, 1),
                currentPosition.x - SPLINE_X[0]
            )

            # Check if a move is below Z = 0
            if float(transformedGCode.x) <= 0.0: 
                print("Warning! Movement below build platform. Check your spline!")

            # Detect implausible moves
            if transformedGCode.x < 0 or np.abs(transformedGCode.x - currentZ) > 50:
                print("Warning! Possibly implausible move detected at height " + str(currentZ) + " mm!")
                outputFile.write(currentLine)
                continue    
            # Check for self-intersection
            if (LAYER_HEIGHT + heightDifference) < 0:
                print("ERROR! Self-intersection at height " + str(currentZ) + " mm! Check your spline!")

            # Check the angle of the printed layer and warn if it's above the machine limit
            if angleSplineThisLayer > (WARNING_ANGLE * np.pi / 180.):
                print("Warning! Spline angle is", (angleSplineThisLayer * 180. / np.pi), "at height", str(currentZ), "mm! Check your spline!")

            if currentLineCommands.e is not None:  # If extrusion is present
                extrusionAmount = float(currentLineCommands.e) * ((LAYER_HEIGHT + heightDifference) / LAYER_HEIGHT)
                if extrusionAmount < MINIMUM_EXTRUSION:
                    extrusionAmount = MINIMUM_EXTRUSION
            else:
                extrusionAmount = None

            feedrate = float(currentLineCommands.f) if currentLineCommands.f is not None else None

            writeLine(1, transformedGCode.y, currentPosition.y, transformedGCode.x, feedrate, extrusionAmount)
            lastPosition = currentPosition
            lastZ = currentZ
        else:
            outputFile.write(currentLine)
    
print("GCode bending finished!")