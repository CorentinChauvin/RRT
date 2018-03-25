#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
    Alexis Dupuis and Corentin Chauvin-Hameau
    SYMOU - 2018

    A* algorithm for (1, 1) robot
"""

from math import sqrt


class Node:
    def __init__(self, q=None):
        if q is None:
            self.q = [0, 0, 0, 0]
        else:
            self.q = q
        self.neighbors = []


def h(node1, node2):
    """ Compute the heuristic cost between two nodes
        (cartesian distance between configurations)
    """

    squares = 0
    for k in range(4):
        squares += (node1.q[k]-node2.q[k])**2

    return sqrt(squares)


def reconstructPath(node):
    """ Return the path taken to arrive to this node
    """

    path = [node]

    while node.cameFrom is not None:
        node = node.cameFrom
        path.insert(0, node)

    return path


def a_star(tree, q_start, q_goal):
    """ Apply the A* algorithm to a tree between start and goal

        Arguments :
            - tree : the list of every nodes
            - q_start : starting configuration
            - q_goal : end configuration
            The nodes of the tree must store their neighbors
        Return :
            - a list of passing nodes if a path is found
            - None if no path is found
    """

    print("Entering a_star...")

    # Finding the nodes corresponding to q_start and q_goal in tree
    start = None
    goal = None
    for node in tree:
        if node.q == q_start:
            start = node
        elif node.q == q_goal:
            goal = node

    if start is None:
        print("Error: start specified for a* not found in the tree")
        return None
    if goal is None:
        print("Error: goal specified for a* not found in the tree")
        return None

    # Initialisation
    closedSet = []
    openSet = [start]

    for node in tree:
        node.g = 100000         # heuristic cost
        node.f = 100000         # total cost
        node.cameFrom = None    # the previous node to reach this one

    start.g = 0
    start.f = h(start, goal)

    # Main loop
    while openSet != []:
        # Find the node with the lowest total score in openSet
        current = openSet[0]
        indexCurrent = 0   # index of the current node in openSet

        for k in range(1, len(openSet)):
            node = openSet[k]

            if node.f < current.f:
                current = node
                indexCurrent = k

        if current == goal:
            print("Success of a_star")
            return reconstructPath(goal)

        del openSet[indexCurrent]
        closedSet.append(current)

        for neighbor in current.neighbors:
            if neighbor in closedSet:
                continue

            if neighbor not in openSet:
                openSet.append(neighbor)

            # Check if the path from this node is not better
            g = current.g + h(current, neighbor)
            if g < neighbor.g:
                neighbor.cameFrom = current
                neighbor.g = g
                neighbor.f = g + h(neighbor, goal)

    # No path found
    print("Failure of a_star")
    return None











