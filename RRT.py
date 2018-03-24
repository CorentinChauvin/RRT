# -*- coding: utf-8 -*-
"""
Adaptation du RRT à un robot (1,1)

Plutôt que de verifier si une configuration tirée est atteignable, on va
fabriquer une configuration atteignable, en tirant une commande et en
l'appliquant à un noeud random de l'arbre.

Puis on testera si le point XY atteint n'entre pas en collision avec un
obstacle...
"""
# (x,y,theta,beta)

import random as rd
from math import sqrt, cos, sin, tan, pi, atan2
import matplotlib.pyplot as plt
from matplotlib.colors import ListedColormap
from copy import deepcopy



#######

class Node:
    def __init__(self, q=None):
        if q is None:
            self.q = [0, 0, 0, 0]
        else:
            self.q = q
        self.neighbors = []
        self.commands = []  # For each neighbor, the command to reach it


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


######


# Global variables:
VOID = 0
OBSTACLE = 1
START = 2
GOAL = 3

step_size = 0.5   # nbr de secondes durant lesquelles on applique la consigne
merge_threshold = 0.5
Ne = 100  # nbr de pas d'integration pour l'application de la commande
Nc = 10  # nbr de pas pour la detection de collision

Vmax = 1.0
BETA_POINTmax = 1.0/step_size
L = 1.5
a = 0.5
d = 0.1

Xmax = 50
Ymax = 50

textMap = loadMapFromText("images/textMap.txt")
Xpix = len(textMap[0])
Ypix = len(textMap)

discriminant = 2.0  # minimal distance between points

###


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
    node_qnear = Node(q_near)

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


#    if globalSuccess:
#        # Check whether the resulting node is not too near from an existing one
#        for node in T:
#            if cart_dist(node.q, q) < discriminant:
#                globalSuccess = False

    if success:
        # Create the resulting node
        q = [x, y, theta, beta]
        node_qnew = Node(q)
        return node_qnear, node_qnew

    print(".")
    return None, None


def extendRRT(T, node_qnear, node_qnew):
    T.append(node_qnew)

    node_qnear.neighbors.append(node_qnew)
    node_qnear.commands.append(0)  # <<<<<<<<<<<<< To complete
    node_qnew.neighbors.append(node_qnear)
    node_qnew.commands.append(0)    # <<<<<<<<<<<<<

    return T


def buildRRT(n, q0, qgoal):
    growtree = 1  # which tree to grow

    # Initialisation of trees
    startNode = Node(q0)
    endNode = Node(qgoal)

    T1 = [startNode]
    T2 = [endNode]

    T = [T1, T2]

    # Main loop
    for i in range(n):
#        if i%(n//100)==0:
#            print(i)

        if growtree == 1:
            (node_qnear, node_qnew) = randomNode(T1)
            if not (node_qnew is None):
                T1 = extendRRT(T1, node_qnear, node_qnew)

                # Try merging qnew to T2:
                for node_q in T2:
                    if (cart_dist(node_q.q[:2], node_qnew.q[:2]) < merge_threshold
                        and not checkcollision(node_q.q, node_qnew.q)):
                        node_qnew.neighbors.append(node_q)
                        node_qnew.commands.append(0)  # <<<<<<<<<<<<<<<<<
                        node_q.neighbors.append(node_qnew)
                        node_q.commands.append(0)  # <<<<<<<<<<<<<<<<<

                        return [T1, T2]

                # If you're here merging hasnt succeeded: swap trees and continue
                growtree = 2    # <<<<<<<<<<<<<<<

        else:
            (node_qnear, node_qnew) = randomNode(T2)
            if not (node_qnew is None):
                T2 = extendRRT(T2, node_qnear, node_qnew)

                # Try merging qnew to T1:
                for node_q in T1:
                    if (cart_dist(node_q.q[:2], node_qnew.q[:2]) < merge_threshold
                        and not checkcollision(node_q.q, node_qnew.q)):
                        node_qnew.neighbors.append(node_q)
                        node_qnew.commands.append(0)  # <<<<<<<<<<<<<<<<<
                        node_q.neighbors.append(node_qnew)
                        node_q.commands.append(0)  # <<<<<<<<<<<<<<<<<

                        return [T1, T2]

                # If you're here merging hasnt succeeded: swap trees and continue
                growtree = 1    # <<<<<<<<<<<<<<<

    print("Nope------------")
    return [T1, T2]


def displayTree(T1, T2):
    global textMap

    """fig = plt.figure()
    ax = fig.add_subplot(111)
    x = [[0,0,0,0,0,0],[0,0,0,0,0,0], [0,1,1,2,1,1], [0,2,3,0,0,1], [0,1,1,1,1,1]]

    # Display matrix
    cax = ax.matshow(x,cmap=cmap)"""

    # Initialise map matrix
    mapMatrix = deepcopy(textMap)

    # Fill the matrix with the nodes of the trees
    for node in T1:
        [x, y] = PIXspace(node.q)
        mapMatrix[x][y] = 5

    for node in T2:
        [x, y] = PIXspace(node.q)
        mapMatrix[x][y] = 6

    # Fill start and goal
    [x, y] = PIXspace(T1[0].q)
    mapMatrix[x][y] = 2
    [x, y] = PIXspace(T2[0].q)
    mapMatrix[x][y] = 3


    # Display matrix
    cmap = ListedColormap(['w', 'k', 'r', 'g', 'r', 'y', 'b'])
    fig = plt.figure()
    ax = fig.add_subplot(111)
    ax.matshow(mapMatrix,cmap=cmap)





[T1, T2] = buildRRT(4000, (35, 30, 0, 0), (1, 1, 1, 0))

displayTree(T1, T2)


