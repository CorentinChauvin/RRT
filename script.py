"""
"""
import numpy as np
import cv2


# Definition of constants
VOID = 0
OBSTACLE = 1
START = 2
GOAL = 3


def loadImageFromImage(path):
    """ Load an image and return its pixel matrix
    """
    image = cv2.imread(path, cv2.IMREAD_COLOR)

    if image is None:
        print("Can't open image")

    return image


def getHeight(image):
    return len(image)


def getWidth(image):
    return len(image[0])


def imageToTextMap(image):
    """ Convert the map from its RGB image to a text representation
    """
    global VOID, OBSTACLE, START, GOAL
    height = getHeight(image)
    width = getWidth(image)

    textMap = [[0]*width for k in range(height)]

    for i in range(height):
        for j in range(width):
            pixel = image[i][j]

            if np.array_equal(pixel, [255, 255, 255]):
                textMap[i][j] = VOID
            elif np.array_equal(pixel, [0, 0, 0]):
                textMap[i][j] = OBSTACLE
            elif np.array_equal(pixel, [255, 0, 0]):    # GRB
                textMap[i][j] = START
            elif np.array_equal(pixel, [0, 255, 0]):
                textMap[i][j] = GOAL

    return textMap


def saveMapAsText(path, textMap):
    """ Save the map in a text file
    """
    height = getHeight(textMap)
    width = getWidth(textMap)

    mapFile = open(path, 'w')

    if mapFile is None:
        print("Map file couldn't be opened")
        return None

    for i in range(height):
        string = ""
        for j in range(width):
            string += str(textMap[i][j])
        string += "\n"
        mapFile.write(string)

    mapFile.close()


def loadMapFromText(path):
    """ Load an image and return its text representation
    """
    mapFile = open(path, 'r')

    if mapFile is None:
        print("Map file couldn't be opened")
        return None

    textMap = []

    for line in mapFile:
        textMap.append([])
        for c in line:
            if c != '\n':
                textMap[-1].append(int(c))

    return textMap



"""imagePath = "images/map1.bmp"
image = loadImageFromImage(imagePath)
textMap = imageToTextMap(image)
saveMapAsText("images/textMap.txt", textMap)
textMap2 = loadMapFromText("images/textMap.txt")
saveMapAsText("images/textMap2.txt", textMap2)"""











