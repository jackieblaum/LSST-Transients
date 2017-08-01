import numpy as np


def shift_without_replacement(v, shiftx, shifty):
    """
    Shift the vector v without replacing elements that falls off the edges (and inserting zeros to keep the shape on
    the other side)

    :param v:
    :param shiftx:
    :param shifty:
    :return:
    """
    return my_roll_x(my_roll_y(v, shifty), shiftx)


def my_roll_x(x, shift):

    if shift > 0:

        return np.pad(x, ((0, 0), (shift, 0)), mode='constant')[:, :-shift]

    elif shift < 0:

        return np.pad(x, ((0, 0), (0, -shift)), mode='constant')[:, -shift:]

    else:

        return x


def my_roll_y(y, shift):
    if shift > 0:

        return np.pad(y, ((shift, 0), (0, 0)), mode='constant')[:-shift, :]

    elif shift < 0:

        return np.pad(y, ((0, -shift), (0, 0)), mode='constant')[-shift:, :]

    else:

        return y
