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

        self.x = x
        self.y = y
        self.side = side

    def __str__(self):
        '''
        Overrides the string function for a Box object.

        :return: a string in the form "box(x, y, width, height, angle)" to create a box region in DS9
        '''

        return "box(" + str(self.x) + ", " + str(self.y) + ", " + str(self.side) + ", " + str(self.side) + ", 0)"


class Grid(object):
    '''
    The constructor takes the x and y coordinates of the center of the image (xcenter, ycenter)
    and the maximum and minimum x and y coordinates of the image (xcenter, ycenter).

    It also sets the overlap factor (overlapfactor) and side length (side) of the boxes in the grid,
    and then calculates the overlap by multiplying these two values.

    An empty array of Box objects is also instantiated.
    '''

    def __init__(self, xcenter, ycenter, xymax, xymin, overlapfactor, side):
        '''
        Constructs a grid object. A box list is instantiated, and the overlap length is set to be the product of the overlap factor and side length.

        :param xcenter: The x-coordinate of the center of the image
        :param ycenter: The y-coordinate of the center of the image
        :param xymax: The maximum x- and y-coordinates of the image
        :param xymin: The minimum x- and y-coordinates of the image
        :param overlapfactor: The fraction of the length of the box by which the boxes overlap one another on the grid
        :param side: The length of the sides of the boxes in the grid
        '''

        self.xcenter = xcenter
        self.ycenter = ycenter
        self.xymax = xymax
        self.xymin = xymin
        self.boxes = []
        self.overlapfactor = overlapfactor
        self.side = side
        self.overlap = overlapfactor * side

    def _out_of_range(self, x, y):
        '''
        A helper method used by do_steps() in order to determine if the Box is in range of the region of
        interest and should be added to the Box array.

        Returns true if out of range, returns false if in range.

        *Not yet fully implemented to account for rotation angle.

        :param x: The x-coordinate of the center of the current Box
        :param y: The y-coordinate of the center of the current Box
        :return: False if the box is out of the region of interest, True otherwise
        '''

        if x > self.xymin and x < self.xymax and y > self.xymin and y < self.xymax:

            return False

        else:

            return True


    def _do_steps(self, steps, x, y, xory, addorsub):
        '''
        Finds the boxes in the current row or column and adds them to the Box array if they are in range.

        :param steps: The number of boxes needed to be added in the row or column
        :param x: The x-coordinate of the current Box.
        :param y: The y-coordinate of the current Box.
        :param xory: Either "x" or "y"; tells whether boxes are being added to the row or column
        :param addorsub: Either "add" or "sub"; tells which direction to move along the row or column
        :return: True if at least one of the boxes is in the region of interest, False otherwise
        '''

        some_in_range = False

        # Needs to be changed to calculate angular distance later
        for i in range(steps):

            if addorsub == "add":

                if xory == "y":
                    y = y + (self.side - self.overlap)

            elif xory == "x":
                    x = x + (self.side - self.overlap)

            elif addorsub == "sub":
                if xory == "y":
                    y = y - (self.side - self.overlap)

                elif xory == "x":
                    x = x - (self.side - self.overlap)

            if not self._out_of_range(x, y):
                self.boxes.append(Box(x, y, self.side))
                some_in_range = True

        return some_in_range


    def _do_round(self, x, y, steps):
        '''
         A helper method used by gen_grid() to add the next round of Box objects to the Box array.

         The some_in_range variables keep track of whether any boxes are within the region of interest. They remain false until a Box is added to the array in the current round.

        :param x: The x-coordinate of the current Box.
        :param y: The y-coordinate of the current Box.
        :param steps: The number of boxes needed to be added in the row or column
        :return: True if at least one Box is in the field of view, False if all boxes in the round
        are out of range.
        '''

        some_in_range_1, some_in_range_2, some_in_range_3, some_in_range_4, some_in_range_5 = False

        if steps == 0:

            self.boxes.append(Box(x, y, self.side))
            some_in_range_1 = True


        else:

            some_in_range_2 = self._do_steps(steps, x, y, y, "sub")
            some_in_range_3 = self._do_steps(steps, x, y, x, "sub")
            some_in_range_4 = self._do_steps(steps, x, y, y, "add")
            some_in_range_5 = self._do_steps(steps, x, y, x, "add")

        return some_in_range_1 or some_in_range_2 or some_in_range_3 or some_in_range_4 or some_in_range_5


    def _gen_grid(self, x, y, steps):
        '''
        Recursion is used to fill the Box array. The base case is when all the Box objects
        created by do_round() fall outside of the field of view. This round of Box objects
        is still added to the array in order to make sure that the edges of the region are
        accounted for. Otherwise, do_round() is called in order to add Box objects to the
        Box array.

        :param x: The x-coordinate of the current Box.
        :param y: The y-coordinate of the current Box.
        :param steps: The number of boxes needed to be added in the row or column; increases by 2 for each round
        :return: No return value
        '''

        if self._do_round(x, y, steps) == 1:

            return None


        else:

            #Needs to be changed to calculate angular distance later
            self._gen_grid(x + self.side - self.overlap, y + self.side - self.overlap, steps + 2)


    def get_grid(self):
        '''
        Calls gen_grid() using the location of the center of the region of interest and
        the starting step value 0.

        Returns the array of Box objects.

        :return:
        '''

        self._gen_grid(self.xcenter, self.ycenter, 0)

        return self.boxes


#Needs to be changed to take user input later
grid = Grid(2000, 2000, 4100, 0, 0.5, 400)
boxes = grid.get_grid()

# Write the grid to a DS9 region file

with open("grid.reg", "w+") as f:
    f.write("physical\n")

    for box in boxes:
        f.write("%s\n" % str(box))
