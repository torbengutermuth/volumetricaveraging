import math
import numpy as np
from grid_TG import GRID
import logging


class GRIDDATA(GRID):
    def __init__(self, min_x, min_y, min_z, max_x, max_y, max_z, gridsize, id_string):
        """
        This function initiates the grid_Data class
        :param min_x: minimal x value that shall be used
        :param min_y: minimal y value
        :param min_z: ...
        :param max_x: ...
        :param max_y: ...
        :param max_z: ...
        :param gridsize: gridisize used
        :param id_string: string to identify the grid
        """
        self.id_string = id_string
        self.loggingname = "griddata_" + self.id_string
        self.logger = logging.getLogger("SCOODER." + self.loggingname)
        self.logger.info(self.loggingname + " of SCOODER was started")
        self.size_gridspace = gridsize
        self.min_x = min_x
        self.min_y = min_y
        self.min_z = min_z
        self.max_x = (
            max_x + 2 * self.size_gridspace
        )  # This addition of a gridspace is so that the number of spaces is good
        self.max_y = (
            max_y + 2 * self.size_gridspace
        )  # Otherwise the float operation might round down
        self.max_z = (
            max_z + 2 * self.size_gridspace
        )  # In addition, one more is needed in case edge cases are filled
        # Because if not, due to the fact that we always look as 000 first and then advance in the positive
        # We would get an axis error if we have a gridvoxel filled on a positive edge

        self.len_x = self.max_x - self.min_x
        self.len_y = self.max_y - self.min_y
        self.len_z = self.max_z - self.min_z

        self.num_spaces_x = int(math.ceil(self.len_x / self.size_gridspace))
        self.num_spaces_y = int(math.ceil(self.len_y / self.size_gridspace))
        self.num_spaces_z = int(math.ceil(self.len_z / self.size_gridspace))

        grid_start = []  # list of lists of lists is generated
        for i in range(self.num_spaces_x):
            column = []
            for i in range(self.num_spaces_y):
                row = []
                for i in range(self.num_spaces_z):
                    row.append(0)
                column.append(row)
            grid_start.append(column)
        grid_start = np.array(grid_start)  # lists are converted to np array
        grid_start = grid_start.astype("float64")
        self.data = grid_start

    def find_equal_spot(self, np_arr, config):
        """
        This function figures out how we have to move a grid to fit the grid we have
        :param np_arr: The numpy array we want to add to the grid
        :param config: The config of the grid we want to add
        :return: results of the change we have to do, vector of two lists with the grid shift and how much of that needs
            to be added to the grid
        """
        min_x_arr = config["min_x"]
        min_y_arr = config["min_y"]
        min_z_arr = config["min_z"]
        start_x = self.min_x
        start_y = self.min_y
        start_z = self.min_z
        a = 0
        b = 0
        c = 0
        # This calculates how much we need to adjust our minimal x,y and z for the new grid that we add
        while start_x < (min_x_arr - self.size_gridspace):
            start_x += self.size_gridspace
            a += 1
        while start_y < (min_y_arr - self.size_gridspace):
            start_y += self.size_gridspace
            b += 1
        while start_z < (min_z_arr - self.size_gridspace):
            c += 1
            start_z += self.size_gridspace
        results = self.calc_distribution(np_arr, config, a, b, c)
        return results

    def calculate_common_volume(self, arr1, arr2):
        ### calculate euclidian between two lists
        if len(arr1) != len(arr2):
            raise BaseException
        xapart = self.size_gridspace - abs(arr1[0] - arr2[0])
        yapart = self.size_gridspace - abs(arr1[1] - arr2[1])
        zapart = self.size_gridspace - abs(arr1[2] - arr2[2])
        volume_different = xapart * yapart * zapart
        effective = volume_different / self.size_gridspace ** 3
        return effective

    def make_edges(self, a, b, c, change):
        e1 = [a, b, c]
        e2 = [a, b, c + change]
        e3 = [a, b + change, c]
        e4 = [a + change, b, c]
        e5 = [a + change, b + change, c]
        e6 = [a + change, b, c + change]
        e7 = [a, b + change, c + change]
        e8 = [a + change, b + change, c + change]
        edges = [e1, e2, e3, e4, e5, e6, e7, e8]
        return edges

    def calc_distribution(self, np_arr, config, x, y, z):
        change = self.size_gridspace
        a = self.min_x + x * self.size_gridspace
        b = self.min_y + y * self.size_gridspace
        c = self.min_z + z * self.size_gridspace
        edges = self.make_edges(a, b, c, change)

        array_oi = [config["min_x"], config["min_y"], config["min_z"]]
        edges_normal = self.make_edges(x, y, z, 1)
        ident_flag = False

        d_ges = 0
        results = []
        i = 0
        for edge in edges:
            volume_together = self.calculate_common_volume(array_oi, edge)
            edgeoi = edges_normal[i]
            d_ges = d_ges + volume_together
            results.append([edgeoi, volume_together])
            i += 1
        for result in results:
            result[1] = result[1] / d_ges
        return results

    def make_changes(self, results, value, x, y, z):
        for result in results:
            a = x + result[0][0]
            b = y + result[0][1]
            c = z + result[0][2]
            self.data[a][b][c] = self.data[a][b][c] + value * result[1]

    def add_np_to_grid(self, np_arr, config):
        shape_arr = np_arr.shape
        shape_self = self.data.shape
        for i in range(3):
            if (
                shape_arr[i] > shape_self[i]
            ):  # This checks if the grid we want to add is bigger than ours, which shouldnt happen
                raise BaseException
        equal = self.find_equal_spot(np_arr, config)
        for x in range(shape_arr[0]):
            for y in range(shape_arr[1]):
                for z in range(shape_arr[2]):
                    value = np_arr[x][y][z]
                    if value > 0:
                        self.make_changes(equal, value, x, y, z)
