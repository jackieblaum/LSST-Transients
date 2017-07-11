class Region(object):
    '''
    A Region object is a square or circular region, many of which form a grid.
    '''

    def __init__(self, x, y, d, shape):
        '''
        Constructs a Region object.

        :param x: The horizontal location of the center of the region
        :param y: The vertical location of the center of the region
        :param d: The length of the sides of the region if it is a square, otherwise the diameter
		  of the region if it is a circle
	:param shape: A string that specifies the shape of the region, either "square" or "circle"
        '''

        self._x = x
        self._y = y
        self._d = d
        self._shape = shape

    @property
    def x(self):

        return self._x

    @property
    def y(self):

        return self._y

    @property
    def d(self):

        return self._d


    def __str__(self):
        '''
        Overrides the string function for a Box object.

        :return: a string in the form "box(_x, _y, width, height, angle)" to create a box region in DS9
        '''
        if self._shape == "square":
            return "box(" + str(self.x) + ", " + str(self._y) + ", " + str(self._d) + ", " + str(self._d) + ", 0)"

        elif self._shape == "circle":
            return "circle(" + str(self.x) + ", " + str(self._y) + ", " + str(self._d) +")"

