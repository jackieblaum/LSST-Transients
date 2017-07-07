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
    

    def pix_to_wcs(self, filename, boxes):
        '''
        Converts the center coordinates for each box in the array from pixels to WCS and
        stores them in a new array.

        :param filename: The name of the file from which the HDUlist will be loaded
        :param boxes: The array of boxes with locations in pixel coordinates
        :return: The array of boxes with locations in WCS
        '''

        arr = []
        for box in boxes:

            arr.append([box.x, box.y])


        boxes_wcs = pix2world(filename, arr)

        return boxes_wcs
    

if __name__ == "__main__":
    
    
    # Collect input from the user
    parser = argparse.ArgumentParser(description="Grid Region Generator")
    parser.add_argument('-i', '--input', type=str, help='Input filename', required=True)
    parser.add_argument('-o', '--output', type=str, help='Output filename', required=True)
    parser.add_argument('-s', '--side', type=float, help='Length of the side of the boxes in pixels', required=True)
    parser.add_argument('-f', '--fraction', type=float, help='Fraction of the box for which they overlap', required=True)

    args = parser.parse_args()


    # Open the header from the fits file inputted by the user
    header = pyfits.getheader(args.input, 0)

    # Set variables by reading from the header
    center_pix_x = header['CRPIX1']
    center_pix_y = header['CRPIX2']
    max_pix_x = header['NAXIS1']
    max_pix_y = header['NAXIS2']
    rotation_angle = header['ROTANG']

    # Generate the grid
    grid = Grid(2000, 2000, 0, 4100, 0, 4100, 0.5, 400)
    boxes = grid.get_grid()

    # Convert the locations of the boxes to WCS
    boxes_wcs = grid.pix_to_wcs(args.input, boxes)

    # Convert the side length of the boxes to WCS
    side_wcs = grid.pix_to_wcs(args.input, [args.side])

    # Write the grid to a DS9 region file
    with open(args.output, "w+") as f:
        f.write("physical\n")

        for box in np.nditer(boxes_wcs):
            f.write("box(%f, %f, %f, %f, 0)" % (box[0], box[1], side_wcs[0], side_wcs[0]))

            
    print("Input file: %s" % args.input)
    print("Output file: %s" % args.output)
    print("Side length of boxes in pixels/WCS: %f/ %f" % (args.side, side_wcs[0]))
    print("Overlap (fraction of box): %f" % args.fraction)
