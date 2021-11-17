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
