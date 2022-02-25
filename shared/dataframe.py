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
    
    find_pos
        FUNCTION: find the position of the first value in linearly increased list that is larger 
                  than or equal to the given value
        SYNTAX:   find_pos(value: int or float, increase_lst: list)
    
    mean_list
        FUNCTION: calculate mean list from a list of lists, excluding 0 numbers
        SYNTAX:   mean_list(lst: list)
        
    list_exclude_zero
        FUNCTION: exclude zero from list y and corresponding index in x
        SYNTAX:   list_exclude_zero(x: list, y: list)
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


def find_pos(value: int or float, increase_lst: list):
    """
    Find the position of the first value in linearly increased list that is larger than or equal to the
    given value

    Note: if value > max(increase_lst), which means such position does not exist in the given list. This
        function will return len(increase_lst), which is the last position of the list + 1

    Usage examples:
    1) used to look for bleach frame
       > bleach_frame = find_pos(bleach_time, time_tseries)

    :param value: int or float
    :param increase_lst: list, has to be a linearly increased list
    :return: out: position of the first value in linearly increased list that is larger than or equal to
                the given value, start from 0

    """
    out = len(increase_lst)-1
    i = 0
    while i < len(increase_lst):
        if value <= increase_lst[i]:
            out = i
            break
        else:
            i += 1

    return out


def mean_list(lst: list):
    """
    Calculate mean list from a list of lists, excluding 0 numbers

    :param lst: list, input list
                [[...], [...], [...], [...],...[...]]
    :return: out: list
    """
    out = []
    for i in range(len(lst[0])):
        n = 0
        s = 0
        for j in range(len(lst)):
            if lst[j][i] != 0:
                s = s + lst[j][i]
                n = n + 1
        out.append(s*1.0/n)

    return out


def list_exclude_zero(x: list, y: list):
    """
    Exclude zero from list y and corresponding index in x

    :param x: list
    :param y: list
    :return:
    """
    x_out = []
    y_out = []
    for i in range(len(y)):
        if y[i] != 0:
            x_out.append(x[i])
            y_out.append(y[i])
    return x_out, y_out
