from box import Box
from utils.cartesian_product import cartesian_product

import numpy as np


class Grid(object):
    '''
    The constructor takes the _x and _y coordinates of the center of the image (xcenter, ycenter)
    and the maximum and minimum _x and _y coordinates of the image (xcenter, ycenter).

    It also sets the overlap factor (overlapfactor) and _side length (_side) of the boxes in the grid,
    and then calculates the overlap by multiplying these two values.

    An empty array of Box objects is also instantiated.
    '''

    def __init__(self, xcenter, ycenter, xmin, xmax, ymin, ymax, overlapfactor, side):
        '''
        Constructs a grid object. A box list is instantiated, and the overlap length is set to be the product of the overlap factor and _side length.

        :param xcenter: The _x-coordinate of the center of the image
        :param ycenter: The _y-coordinate of the center of the image
        :param xymax: The maximum _x- and _y-coordinates of the image
        :param xymin: The minimum _x- and _y-coordinates of the image
        :param overlapfactor: The fraction of the length of the box by which the boxes overlap one another on the grid
        :param side: The length of the sides of the boxes in the grid
        '''

        self.xcenter = xcenter
        self.ycenter = ycenter

        self._xmin = xmin
        self._xmax = xmax
        self._ymin = ymin
        self._ymax = ymax

        self.overlapfactor = overlapfactor
        self.side = side
        self.overlap = overlapfactor * side

        self.debug = False

    def _in_range(self, a_box):
        '''
        A helper method used by do_steps() in order to determine if the Box is in range of the region of
        interest and should be added to the Box array.

        Returns true if out of range, returns false if in range.

        *Not yet fully implemented to account for rotation angle.

        :param a_box: the box to be checked
        :return: False if the box is out of the region of interest, True otherwise
        '''

        x = a_box.x
        y = a_box.y

        if self._xmin < x < self._xmax and self._ymin < y < self._ymax:

            return True

        else:

            return False

    # def _do_steps(self, steps, _x, _y, xory, addorsub):
    #     '''
    #     Finds the boxes in the current row or column and adds them to the Box array if they are in range.
    #
    #     :param steps: The number of boxes needed to be added in the row or column
    #     :param _x: The _x-coordinate of the current Box.
    #     :param _y: The _y-coordinate of the current Box.
    #     :param xory: Either "_x" or "_y"; tells whether boxes are being added to the row or column
    #     :param addorsub: Either "add" or "sub"; tells which direction to move along the row or column
    #     :return: True if at least one of the boxes is in the region of interest, False otherwise
    #     '''
    #
    #     some_in_range = False
    #
    #     # Needs to be changed to calculate angular distance later
    #     for i in range(steps):
    #
    #         if addorsub == "add":
    #
    #             if xory == "_y":
    #                 _y = _y + (self._side - self.overlap)
    #
    #         elif xory == "_x":
    #             _x = _x + (self._side - self.overlap)
    #
    #         elif addorsub == "sub":
    #             if xory == "_y":
    #                 _y = _y - (self._side - self.overlap)
    #
    #             elif xory == "_x":
    #                 _x = _x - (self._side - self.overlap)
    #
    #         if not self._out_of_range(_x, _y):
    #             self.boxes.append(Box(_x, _y, self._side))
    #             some_in_range = True
    #
    #     return some_in_range
    #
    # def _do_round(self, _x, _y, steps):
    #     '''
    #      A helper method used by gen_grid() to add the next round of Box objects to the Box array.
    #
    #      The some_in_range variables keep track of whether any boxes are within the region of interest. They remain
    #      false until a Box is added to the array in the current round.
    #
    #     :param _x: The _x-coordinate of the current Box.
    #     :param _y: The _y-coordinate of the current Box.
    #     :param steps: The number of boxes needed to be added in the row or column
    #     :return: True if at least one Box is in the field of view, False if all boxes in the round
    #     are out of range.
    #     '''
    #
    #     some_in_range_1, some_in_range_2, some_in_range_3, some_in_range_4, some_in_range_5 = False
    #
    #     if steps == 0:
    #
    #         self.boxes.append(Box(_x, _y, self._side))
    #         some_in_range_1 = True
    #
    #
    #     else:
    #
    #         some_in_range_2 = self._do_steps(steps, _x, _y, _y, "sub")
    #         some_in_range_3 = self._do_steps(steps, _x, _y, _x, "sub")
    #         some_in_range_4 = self._do_steps(steps, _x, _y, _y, "add")
    #         some_in_range_5 = self._do_steps(steps, _x, _y, _x, "add")
    #
    #     return some_in_range_1 or some_in_range_2 or some_in_range_3 or some_in_range_4 or some_in_range_5
    #
    # def _gen_grid(self, _x, _y, steps):
    #     '''
    #     Recursion is used to fill the Box array. The base case is when all the Box objects
    #     created by do_round() fall outside of the field of view. This round of Box objects
    #     is still added to the array in order to make sure that the edges of the region are
    #     accounted for. Otherwise, do_round() is called in order to add Box objects to the
    #     Box array.
    #
    #     :param _x: The _x-coordinate of the current Box.
    #     :param _y: The _y-coordinate of the current Box.
    #     :param steps: The number of boxes needed to be added in the row or column; increases by 2 for each round
    #     :return: No return value
    #     '''
    #
    #     if self._do_round(_x, _y, steps) == 1:
    #
    #         return None
    #
    #
    #     else:
    #
    #         # Needs to be changed to calculate angular distance later
    #         self._gen_grid(_x + self._side - self.overlap, _y + self._side - self.overlap, steps + 2)

    def _do_round(self, delta):

        # Let's figure out the iteration we are in
        iteration = delta / (self.side - self.overlap)

        # Let's compute the start point of the round
        start_x = self.xcenter + delta
        start_y = self.ycenter + delta

        # How many steps do we need to complete one _side of the round?
        steps = iteration * 2

        xs = np.linspace(start_x - steps * (self.side - self.overlap), start_x, steps + 1)
        ys = np.linspace(start_y - steps * (self.side - self.overlap), start_y, steps + 1)

        this_boxes = map(lambda (x, y): Box(x, y, self.side), cartesian_product([xs,ys]))

        return this_boxes

    @staticmethod
    def distance(x1, y1, x2, y2):

        return np.sqrt((x1-x2)**2 + (y1-y2)**2)

    def get_grid(self):
        '''
        Calls gen_grid() using the location of the center of the region of interest and
        the starting step value 0.

        Returns the array of Box objects.

        :return:
        '''

        # Figure out maximum diagonal
        d1 = self.distance(self._xmin, self._ymin, self.xcenter, self.ycenter)
        d2 = self.distance(self._xmin, self._ymax, self.xcenter, self.ycenter)

        longest_diagonal = max(d1, d2)

        all_boxes = self._do_round(longest_diagonal)

        survived_boxes = filter(self._in_range, all_boxes)

        return survived_boxes

        # iteration = 1
        #
        # boxes = [Box(self.xcenter, self.ycenter, self.side)]
        #
        # # Infinite loop, will exit with break
        # while True:
        #
        #     # Compute the delta in coordinates
        #     delta = iteration * (self.side - self.overlap)
        #
        #     print(delta)
        #
        #     # This will contain the boxes generated in this round
        #     boxes_in_this_round = self._do_round(delta)
        #
        #     # print("Generated %s boxes" % (len(boxes_in_this_round)))
        #     #
        #     # min_x = min(map(lambda this_box:this_box.x, boxes_in_this_round))
        #     # min_y = min(map(lambda this_box: this_box.y, boxes_in_this_round))
        #     #
        #     # if min_x < 0:
        #     #
        #     #     self.debug = True
        #     #
        #     #     #import pdb;pdb.set_trace()
        #     #
        #     # print("Minimum _x: %s, Minimum _y: %s" % (min_x, min_y))
        #
        #     # Filter out boxes that are outside of the field of view
        #     #survived_boxes = filter(self._in_range, boxes_in_this_round)
        #
        #     survived_boxes = []
        #
        #     for box in boxes_in_this_round:
        #
        #         if self._in_range(box):
        #
        #             #print("Box %s succeeded" % box)
        #
        #             survived_boxes.append(box)
        #
        #         else:
        #             pass
        #             #print("Box %s failed" % box)
        #
        #     print("%s boxes survived" % (len(survived_boxes)))
        #
        #     if len(survived_boxes) == 0 or len(boxes_in_this_round) > 5000:
        #
        #         # We are finished
        #         break
        #
        #     else:
        #
        #         # Let's add these boxes and go on with the next round
        #
        #         boxes.extend(survived_boxes)
        #
        #         iteration += 1

        return boxes


if __name__ == "__main__":

    # Needs to be changed to take user input later
    grid = Grid(2000, 2000, 0, 4100, 0, 4100, 0.5, 400)
    boxes = grid.get_grid()

    # Write the grid to a DS9 region file

    with open("grid.reg", "w+") as f:
        f.write("physical\n")

        for box in boxes:
            f.write("%s\n" % str(box))