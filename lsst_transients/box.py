class Box(object):
    '''
    A box object is a square region, many of which form a grid.
    '''

    def __init__(self, x, y, side):
        '''
        Constructs a Box object.

        :param x: The horizontal location of the center of the box
        :param y: The vertical location of the center of the box
        :param side: The length of the sides of the box
        '''

        self._x = x
        self._y = y
        self._side = side

    @property
    def x(self):

        return self._x

    @property
    def y(self):

        return self._y

    @property
    def side(self):

        return self._side

    def __str__(self):
        '''
        Overrides the string function for a Box object.

        :return: a string in the form "box(_x, _y, width, height, angle)" to create a box region in DS9
        '''

        return "box(" + str(self.x) + ", " + str(self._y) + ", " + str(self._side) + ", " + str(self._side) + ", 0)"
