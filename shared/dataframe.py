import numpy as np
import random


"""
# ---------------------------------------------------------------------------------------------------
# FUNCTIONS for PD.DATAFRAME/LIST/UMANAGER/NUMBER
# ---------------------------------------------------------------------------------------------------
List related:

    list_unwrap
        FUNCTION: unwrap one layer of a list with first element per sublist
        SYNTAX:   list_unwrap(lst: list)
    
    str_to_float
        FUNCTION: transform a string into a list of floats
        SYNTAX:   str_to_float(string: str)
"""


def list_unwrap(lst: list):
    """
    Unwrap one layer of a list with first element per sublist

    Examples:
    input list: [[a1,a2],[b1,b2],[c1,c2],[d1,d2]]
    output list: [a1,b1,c1,d1]
    :param lst: list, the list to be unwrapped
    :return: out: list
    """
    out = list()
    for i in range(len(lst)):
        out.append(lst[i][0])
    return out


def str_to_float(string: str):
    """
    Transform a string into a list of floats

    Examples:
    input string: (24 characters)
    [5.55, 6.53, 7.35, 8.91]
    output list: (4 elements)
    [5.55, 6.53, 7.35, 8.91]
    :param string: str, string to be converted
    :return: out: list
    """
    out = [float(i) for i in string[1:-1].split(', ')]
    return out
