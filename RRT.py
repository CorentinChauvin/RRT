#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
    Alexis Dupuis and Corentin Chauvin-Hameau
    SYMOU - 2018

    Adaptation of the RRT algorithm for a (1, 1) mobile robot
"""

import random as rd
from math import sqrt, cos, sin, tan, pi
from time import clock
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from copy import deepcopy
from a_star import Node, a_star
from map_utilitary import loadMapFromText


# Global variables
step_size = 0.5         # number of seconds we apply the command
merge_threshold = 0.5   # cartesian distance needed between two nodes needed to merge
Ne = 100    # number of integration steps for the command
Nc = 10     # number of steps for the collision detection

Vmax = 1.0  # maximum velocity of the robot
L = 1.5     # half the distance between the two fixed wheels
a = 0.5     # distance between the centers of two fixed wheels and the steering wheel
d = 0.1     # distance between the steering wheel and the "sugar" (controlled point)

textMap = loadMapFromText("images/textMap.txt")
Xmax = 50               # dimensions of the map in meters
Ymax = 50
Xpix = len(textMap[0])  # dimensions of the map in pixels
Ypix = len(textMap)

discriminant = 0.15         # minimal cartesian distance between nodes
discriminantRefusals = 0    # number of refusals due to discriminant
maxMergesNumber = 20        # maximum of merges before stopping the algorithm


def XYspace(q):
    """ Extract the 2D coordinates of the configuration
    """
    return q[:2]


def PIXspace(q):
    """ Convert the 2D coordinates in pixel coordinates
    """
    q = XYspace(q)
    q = (int(q[0]*Xpix/Xmax), int(q[1]*Ypix/Ymax))
    return q


def cart_dist(q1, q2):
    """ Cartesian distance between two vectors of same dimensions
    """
    if len(q1) != len(q2):
        print("ERROR: computing cartesian distance of vectors of different dimensions")
        return -1

    squares = 0
    for k in range(len(q1)):
        squares += (q1[k]-q2[k])**2

    return sqrt(squares)


def checkcollision(q1, q2):
    """ Is there an obstacle in the strait line joining 2 configs ?
    """
    for i in range(1, Nc+1):
        (xi, yi) = (q1[0]+i*(q2[0]-q1[0])/Nc, q1[1]+i*(q2[1]-q1[1])/Nc)
        (xi, yi) = PIXspace((xi, yi))
        if textMap[xi][yi] != 0:
            return True
    return False


def sugarPosition(q, a, d):
    """ Compute the position of the controlled point (sugar)

        Parameters :
            - q : configuration of the robot$
            - a : distance between the steerable wheel and the robot origin
            - d : distance between the sugar and the steerable wheel
    """
    x = q[0]
    y = q[1]
    theta = q[2]
    beta = q[3]

    x_sugar = x + a*cos(theta) + d*cos(theta + beta)
    y_sugar = y + a*sin(theta) + d*sin(theta + beta)

    return (x_sugar, y_sugar)


def randomNode(T):
    """ Choose a random command and apply it to a node
    """

    global discriminantRefusals

    # Pick a random config
    x_rand = rd.random()*Xmax
    y_rand = rd.random()*Ymax
    theta_rand = rd.random()*2*pi
    beta_rand = rd.random()*2*pi

    q_rand = [x_rand, y_rand, theta_rand, beta_rand]

    # Find the closest node of this random config in T
    bestDistance = 100000
    bestNode = None

    for node in T:
        distance = cart_dist(node.q, q_rand)
        if distance < bestDistance:
            bestDistance = distance
            bestNode = node

    q_near = bestNode.q
    node_qnear = bestNode

    # Compute the desired velocity of the sugar
    sugarPos = sugarPosition(q_near, a, d)
    desiredSugarPos = sugarPosition(q_rand, a, d)

    dx = desiredSugarPos[0] - sugarPos[0]
    dy = desiredSugarPos[1] - sugarPos[1]

    vx = Vmax * dx / sqrt(dx**2 + dy**2)
    vy = Vmax * dy / sqrt(dx**2 + dy**2)

    # Applying the sugar tracking command
    [x, y, theta, beta] = q_near
    success = True  # to check whether this configuration is ok

    for i in range(Ne):
        xm_dot = beta_dot = 0

        if beta != 0:
            A = L/tan(beta)*cos(theta) - a*sin(theta) - d*sin(theta+beta)
            B = -d*sin(theta+beta)
            C = L/tan(beta)*sin(theta) + a*cos(theta) + d*cos(theta+beta)
            D = d*cos(theta+beta)

            invDet = 1 / (A*D - B*C)

            xm_dot = L/tan(beta) * invDet * (D*vx - B*vy)
            beta_dot = invDet * (A*vy - C*vx)
        else:
            xm_dot = Vmax
            beta_dot = 1/d* (cos(theta)*vy - sin(theta)*vx)

        beta += step_size/Ne*beta_dot
        theta += step_size/Ne*tan(beta)/L
        x += step_size/Ne*cos(theta)*xm_dot
        y += step_size/Ne*sin(theta)*xm_dot

        # Check for collision
        if textMap[int(x*Xpix/Xmax)][int(y*Ypix/Ymax)] != 0:
            success = False  # We found at least one collision
            break

        beta = beta % (2*pi)
        theta = theta % (2*pi)

    q = [x, y, theta, beta]

    if success:
        # Check whether the resulting node is not too near from an existing one
        for node in T:
            if cart_dist(node.q, q) < discriminant:
                success = False
                discriminantRefusals = discriminantRefusals + 1
                break

    if success:
        # Create the resulting node
        node_qnew = Node(q)
        return node_qnear, node_qnew

    return None, None


def extendRRT(T, node_qnear, node_qnew):
    """ Extend the tree T given the nodes of qnear and qnew
        Add also a link between these two nodes
    """

    T.append(node_qnew)

    node_qnear.neighbors.append(node_qnew)
    node_qnew.neighbors.append(node_qnear)


def buildRRT(n, q0, qgoal):
    """ Perform the RRT algorithm
    """

    global discriminantRefusals

    # Initialisation of trees
    startNode = Node(q0)
    endNode = Node(qgoal)

    T = [[startNode], [endNode]]

    currentTree = 0  # Tree to grow
    otherTree = 1

    mergesNumber = 0  # Number of merges achieved

    # Main loop
    i = 0
    while i < n and mergesNumber < maxMergesNumber:
        i = i + 1
        (node_qnear, node_qnew) = randomNode(T[currentTree])

        if not (node_qnew is None):
            extendRRT(T[currentTree], node_qnear, node_qnew)

            # Try to merge qnew to the other tree:
            for node_q in T[otherTree]:
                if (cart_dist(node_q.q[:2], node_qnew.q[:2]) < merge_threshold
                    and not checkcollision(node_q.q, node_qnew.q)):

                    node_qnew.neighbors.append(node_q)
                    node_q.neighbors.append(node_qnew)

                    mergesNumber = mergesNumber + 1
                    break

        # If you're here merging hasn't succeeded: swap trees and continue
        currentTree = (currentTree + 1) % 2
        otherTree = (otherTree + 1) % 2


    print("Number of merges : {}".format(mergesNumber))
    print("Discriminant refusals : {}".format(discriminantRefusals))
    return T


def displayTree(T1, T2, path=None):
    """ Display the two trees
        If path is given, display also the path generated by A*
    """

    global textMap

    # Initialise map matrix
    mapMatrix = deepcopy(textMap)

    # Fill the matrix with the nodes of the trees
    for node in T1:
        [x, y] = PIXspace(node.q)
        mapMatrix[x][y] = 5

    for node in T2:
        [x, y] = PIXspace(node.q)
        mapMatrix[x][y] = 6

    # Display path on an other image
    pathMatrix = deepcopy(mapMatrix)
    if path is not None:
        for node in path:
            [x, y] = PIXspace(node.q)
            pathMatrix[x][y] = 7

    # Fill start and goal
    [x, y] = PIXspace(T1[0].q)
    mapMatrix[x][y] = 2
    pathMatrix[x][y] = 2
    [x, y] = PIXspace(T2[0].q)
    mapMatrix[x][y] = 3
    pathMatrix[x][y] = 3

    # Display image without path
    fig = plt.figure(1)
    cmap1 = ListedColormap(['w', 'k', 'r', 'g', 'r', 'y', 'b'])
    ax1 = fig.add_subplot(121)
    ax1.matshow(mapMatrix, cmap=cmap1)

    # Display image with path
    cmap2 = ListedColormap(['w', 'k', 'r', 'g', 'r', 'y', 'b', 'grey'])
    ax2 = fig.add_subplot(122)
    ax2.matshow(pathMatrix, cmap=cmap2)



startTime = clock()
[T1, T2] = buildRRT(10000, [35, 30, 0, 0], [1, 1, 1, 0])
RRTTime = clock()

T = T1 + T2
path = a_star(T, [35, 30, 0, 0], [1, 1, 1, 0])
a_starTime = clock()

displayTree(T1, T2, path)

print("RRT execution time : {}".format(RRTTime - startTime))
print("A* execution time : {}".format(a_starTime - RRTTime))
