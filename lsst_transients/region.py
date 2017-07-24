import numpy as np

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


    def get_boundingbox(self, maxx, maxy):
        '''

        :return corners: An array of sky coordinate pairs for the four corners of the bounding box with a padding of 5 pixels
        '''

        padding = 5

        corner1 = [int(np.floor(float(self.x) - (float(self.d) / 2.0) - padding)),
                   int(np.floor(float(self.y) - (float(self.d) / 2.0) - padding))]
        corner2 = [int(np.floor(float(self.x) - (float(self.d) / 2.0) - padding)),
                   int(np.ceil(float(self.y) + (float(self.d) / 2.0) + padding))]
        corner3 = [int(np.ceil(float(self.x) + (float(self.d) / 2.0) + padding)),
                   int(np.floor(float(self.y) - (float(self.d) / 2.0) - padding))]
        corner4 = [int(np.ceil(float(self.x) + (float(self.d) / 2.0) + padding)),
                   int(np.ceil(float(self.y) + (float(self.d) / 2.0) + padding))]

        if corner1[0] < 1:
            corner1[0] = 1
        if corner1[1] < 1:
            corner1[1] = 1
        if corner2[0] < 1:
            corner2[0] = 1
        if corner2[1] > int(maxy):
            corner2[1] = int(maxy)
        if corner3[0] > int(maxx):
            corner3[0] = int(maxx)
        if corner3[1] < 1:
            corner3[1] = 1
        if corner4[0] > int(maxx):
            corner4[0] = int(maxx)
        if corner4[1] > int(maxy):
            corner4[1] = int(maxy)

        return [corner1, corner2, corner3, corner4]


