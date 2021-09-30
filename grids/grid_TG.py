import itertools
import logging
import math
import os

import numpy as np
from gridData import Grid as grid_alt


class Grid:
    def __init__(self, id_string, data_grid_pre, size_gridspace, stauch):
        """
        initializes grid with necessary information
        :param id_string: name of the grid, used for saving
        :param data_grid_pre: data to initiate the grid on
        :param size_gridspace: size of the gridvoxels in each dimension
        :param stauch: stauch factor, which results in certain shineout
        """
        data_grid = data_grid_pre.copy()
        self.id_string = str(id_string)
        self.loggingname = "grid_" + self.id_string
        self.logger = logging.getLogger("grids." + self.loggingname)
        self.logger.info(self.loggingname + " of grids was started")
        blub = 6 + 12 * (1 / math.sqrt(2)) + 8 * (1 / math.sqrt(3))
        if stauch < 0 or stauch > 1:
            self.logger.error("Stauch factor not between 0 and 1 but:" + str(stauch))
            raise RuntimeError
        self.stauch = stauch
        self.shineout = blub * self.stauch / (blub * self.stauch + 1)
        self.value_alter = 1 / (1 + blub * self.stauch)  # 1x +  15.2368xt = 1
        if size_gridspace <= 0:
            self.logger.error(
                "Size gridspace below or equal to 0 :" + str(size_gridspace)
            )
            raise RuntimeError
        self.size_gridspace = size_gridspace
        self.min_x = min(data_grid["X"]) - self.size_gridspace
        self.min_y = min(data_grid["Y"]) - self.size_gridspace
        self.min_z = min(data_grid["Z"]) - self.size_gridspace
        self.max_x = max(data_grid["X"]) + self.size_gridspace
        self.max_y = max(data_grid["Y"]) + self.size_gridspace
        self.max_z = max(data_grid["Z"]) + self.size_gridspace

        self.x_change = (self.max_x - self.min_x) % self.size_gridspace
        self.max_x = self.max_x + (self.size_gridspace - self.x_change)
        self.len_x = self.max_x - self.min_x

        self.y_change = (self.max_y - self.min_y) % self.size_gridspace
        self.max_y = self.max_y + (self.size_gridspace - self.y_change)
        self.len_y = self.max_y - self.min_y

        self.z_change = (self.max_z - self.min_z) % self.size_gridspace
        self.max_z = self.max_z + (self.size_gridspace - self.z_change)
        self.len_z = self.max_z - self.min_z

        data_grid["X"] = data_grid["X"] - self.min_x
        data_grid["Y"] = data_grid["Y"] - self.min_y
        data_grid["Z"] = data_grid["Z"] - self.min_z

        self.num_spaces_x = int(math.ceil(self.len_x / self.size_gridspace)) + 1
        self.num_spaces_y = int(math.ceil(self.len_y / self.size_gridspace)) + 1
        self.num_spaces_z = int(math.ceil(self.len_z / self.size_gridspace)) + 1
        self.dim = self.num_spaces_x * self.num_spaces_y * self.num_spaces_z
        self.logger.debug(
            "Num spaces x"
            + str(self.num_spaces_x)
            + "Num spaces y"
            + str(self.num_spaces_y)
            + "Num spaces z"
            + str(self.num_spaces_z)
            + "min x"
            + str(self.min_x)
            + "max x"
            + str(self.max_x)
            + "min y"
            + str(self.min_y)
            + "max y"
            + str(self.max_y)
            + "min z"
            + str(self.min_z)
            + "max z"
            + str(self.max_z)
        )

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

    def alter_grid(self, coord1, coord2, coord3):
        self.logger.debug(
            "alter grid function called : "
            + str(coord1)
            + " "
            + str(coord2)
            + " "
            + str(coord3)
        )
        """
        This function will alter the grid at a certain point and this change will "shine outside"
         Now we want to let that hit "shine outside" of its coordinate
         There are three different kinds of cubes outside of the core, middle, side edge and edge
         The shine out effect has a factor that softens it (stauch) as well as factors depending on the geometry of the cube towards the middle (1/length_to_core)
         core = 1x 6 middles * 1 = 6xt 12 sides * 0.7071 = 8.4852xt 8 corners * 0.5773 = 4.6184xt = 1x + 15.2368xt
         1x +  15.2368xt = 1
        :param coord1: x coord
        :param coord2: y coord
        :param coord3: z coord
        :return:
        """
        for a, b, c in itertools.product(range(3), range(3), range(3)):
            a = a - 1
            b = b - 1
            c = c - 1
            if abs(a) + abs(b) + abs(c) == 0:
                self.data[coord1][coord2][coord3] = (
                    self.data[coord1][coord2][coord3] + self.value_alter
                )
            elif abs(a) + abs(b) + abs(c) == 1:
                if (coord1 + a) < 0 or (coord2 + b) < 0 or (coord3 + c) < 0:
                    self.logger.warning(
                        "Gridnumber too low "
                        + str(coord1)
                        + "_"
                        + str(coord2)
                        + "_"
                        + str(coord3)
                        + "_"
                        + str(a)
                        + "_"
                        + str(b)
                        + "_"
                        + str(c)
                    )
                    continue
                elif (
                    (coord1 + a) >= self.num_spaces_x
                    or (coord2 + b) >= self.num_spaces_y
                    or (coord3 + c) >= self.num_spaces_z
                ):
                    self.logger.warning(
                        "Gridnumber too high "
                        + str(coord1)
                        + "_"
                        + str(coord2)
                        + "_"
                        + str(coord3)
                        + "_"
                        + str(a)
                        + "_"
                        + str(b)
                        + "_"
                        + str(c)
                    )
                    continue
                else:
                    self.data[coord1 + a][coord2 + b][coord3 + c] = (
                        self.data[coord1 + a][coord2 + b][coord3 + c]
                        + self.value_alter * self.stauch * 1
                    )
            elif abs(a) + abs(b) + abs(c) == 2:
                if (coord1 + a) < 0 or (coord2 + b) < 0 or (coord3 + c) < 0:
                    self.logger.warning(
                        "Gridnumber too low "
                        + str(coord1)
                        + "_"
                        + str(coord2)
                        + "_"
                        + str(coord3)
                        + "_"
                        + str(a)
                        + "_"
                        + str(b)
                        + "_"
                        + str(c)
                    )
                    continue
                elif (
                    (coord1 + a) >= self.num_spaces_x
                    or (coord2 + b) >= self.num_spaces_y
                    or (coord3 + c) >= self.num_spaces_z
                ):
                    self.logger.warning(
                        "Gridnumber too high "
                        + str(coord1)
                        + "_"
                        + str(coord2)
                        + "_"
                        + str(coord3)
                        + "_"
                        + str(a)
                        + "_"
                        + str(b)
                        + "_"
                        + str(c)
                    )
                    continue
                else:
                    self.data[coord1 + a][coord2 + b][coord3 + c] = self.data[
                        coord1 + a
                    ][coord2 + b][coord3 + c] + self.value_alter * self.stauch * (
                        1 / math.sqrt(2)
                    )
            elif abs(a) + abs(b) + abs(c) == 3:
                if (coord1 + a) < 0 or (coord2 + b) < 0 or (coord3 + c) < 0:
                    self.logger.warning(
                        "Gridnumber too low "
                        + str(coord1)
                        + "_"
                        + str(coord2)
                        + "_"
                        + str(coord3)
                        + "_"
                        + str(a)
                        + "_"
                        + str(b)
                        + "_"
                        + str(c)
                    )
                    continue
                elif (
                    (coord1 + a) >= self.num_spaces_x
                    or (coord2 + b) >= self.num_spaces_y
                    or (coord3 + c) >= self.num_spaces_z
                ):
                    self.logger.warning(
                        "Gridnumber too high "
                        + str(coord1)
                        + "_"
                        + str(coord2)
                        + "_"
                        + str(coord3)
                        + "_"
                        + str(a)
                        + "_"
                        + str(b)
                        + "_"
                        + str(c)
                    )
                    continue
                else:
                    self.data[coord1 + a][coord2 + b][coord3 + c] = self.data[
                        coord1 + a
                    ][coord2 + b][coord3 + c] + self.value_alter * self.stauch * (
                        1 / math.sqrt(3)
                    )

    def sum_grid(self):
        """
        sums up the grid
        :return: sum of the grid
        """
        return np.sum(a=self.data.flatten())

    def max_grid(self):
        """
        calculates the maximum value in the grid
        :return: max value of grid
        """
        return np.amax(self.data.flatten())

    def save_grid_density(self, extra_string=""):
        self.logger.debug("save grid density function called")
        origin = (self.min_x, self.min_y, self.min_z)
        delta = (self.size_gridspace, self.size_gridspace, self.size_gridspace)
        g = grid_alt(self.data, origin=origin, delta=delta)
        g.export(self.id_string + "-" + extra_string + ".dx", type="double")

    def save_PDB(self, cutoff_perc):
        """
        This converts the grid to a PDB file
        :param cutoff_perc:
        :return:
        """
        self.logger.debug("save pdb function called percentage: " + str(cutoff_perc))
        grid_temp = self.data.copy()
        summe_del = 0
        summe_ges = np.sum(grid_temp)
        count_gridspaces = 0
        list_coordinates = []
        array_sorted = np.sort(grid_temp, axis=None)
        while summe_del < (cutoff_perc * summe_ges):
            count_gridspaces = count_gridspaces + 1
            summe_del = summe_del + array_sorted[-count_gridspaces]
            hit = array_sorted[-count_gridspaces]
            coordinate = np.where(grid_temp == hit)
            x = (
                coordinate[0] * self.size_gridspace
                + self.min_x
                + 0.5 * self.size_gridspace
            )  # Here 0.5 Gridspaces are added to adjust for
            y = (
                coordinate[1] * self.size_gridspace
                + self.min_y
                + 0.5 * self.size_gridspace
            )
            z = (
                coordinate[2] * self.size_gridspace
                + self.min_z
                + 0.5 * self.size_gridspace
            )  # To correct from the start of the space to the middle
            coordinates_trans = [x, y, z]
            list_coordinates.append(coordinates_trans)
        pdb = " "
        counter = 1
        final_pdb = ""
        mark = "HETATM"
        atomname = "TB"
        astring = " "
        resid = 1
        occ = 1
        temp = 1
        elname = "Mg"
        for coordinates in list_coordinates:
            for ccount in range(len(coordinates[0])):
                a = mark.ljust(6)  # atom#6s
                b = str(counter).rjust(5)  # aomnum#5d
                c = atomname.center(4)  # atomname$#4s
                d = atomname.ljust(3)  # resname#1s
                e = astring.rjust(1)  # Astring
                f = str(resid).rjust(4)  # resnum
                g = str("%8.3f" % (float(coordinates[0][ccount]))).rjust(8)  # x
                h = str("%8.3f" % (float(coordinates[1][ccount]))).rjust(8)  # y
                i = str("%8.3f" % (float(coordinates[2][ccount]))).rjust(8)  # z\
                j = str("%6.2f" % (float(occ))).rjust(6)  # occ
                k = str("%6.2f" % (float(temp))).ljust(6)  # temp
                l = elname.rjust(12)  # elname
                # pdb_line = "%s%s %s %s %s%s    %s%s%s%s%s%s\n"% a,b,c,d,e,f,g,h,i,j,k,l,l,l,l,l,l,l
                pdb_line = "{}{} {} {} {}{}    {}{}{}{}{}{}\n".format(
                    a, b, c, d, e, f, g, h, i, j, k, l
                )
                counter = counter + 1
                final_pdb = final_pdb + pdb_line
        tosave = self.id_string + str(int(cutoff_perc * 100)) + ".pdb"
        if os.path.exists(tosave):
            os.remove(tosave)
        f = open(tosave, "w")
        f.write(final_pdb)
        f.close()
        self.logger.info(
            "PDB of "
            + self.loggingname
            + " saved to "
            + tosave
            + " from "
            + os.getcwd()
        )

    def subtract_grid(self, other_grid):
        """
        This subtracts two identically setup grids from each other
        This is only done after a check that they are identical
        If not identical, errors or plain wrong results are to be expected
        :param other_grid:
        :return:
        """
        self.logger.debug("subtract grid function called")
        if (
            self.min_x == other_grid.min_x
            and self.max_x == other_grid.max_x
            and self.max_y == other_grid.max_y
            and self.min_y == other_grid.min_y
            and self.max_z == other_grid.max_z
            and self.min_z == other_grid.min_z
        ):
            sum_a = self.sum_grid()
            sum_b = other_grid.sum_grid()
            print("Grid A sum: " + str(sum_a))
            print("Grid B sum: " + str(sum_b))
            frac = sum_a / sum_b
            print("Fraction a/b: " + str(frac))
            self.data = np.absolute(self.data - other_grid.data)
            sum_post = self.sum_grid()
            print("Subtractive grid sum: " + str(sum_post))
            frac_sub = sum_post / (sum_a + sum_b)
            print("Subtractive grid fraction: " + str(frac_sub))
            frac_min = np.absolute(sum_a - sum_b) / (sum_a + sum_b)
            print("Minimal post grid fraction: " + str(frac_min))
            return [sum_a, sum_b, frac, sum_post, frac_sub, frac_min]
        else:
            self.logger.error("The grids are unequal, subtraction not done!")
            raise RuntimeError("The grids are unequal, subtraction not done!")

    def euclidian_type_distance(self):
        """
        Calculates a euclidian type distance for a subtractive grid
        Using this calculation is only rational if the grid investigated is a subtractive grid
        :return:
        """
        data_raw = self.data.copy()
        euclidian = math.sqrt(np.sum(np.square(data_raw)))
        return euclidian

    def volume_used(self):
        """
        Returns the used volume in the grid
        :return:
        """
        nonzero_count = np.count_nonzero(self.data)
        nonzero_volume = nonzero_count * (self.size_gridspace ** 3)
        return nonzero_volume

    def volume_cube(self):
        """
        Returns the total volume of the grid
        :return:
        """
        nonzero_volume_cube = self.dim * (self.size_gridspace ** 3)
        return nonzero_volume_cube

    def indeces_nonzero(self):
        """
        eturns indices for all nonzero voxels
        :return:
        """
        indeces = np.transpose(np.nonzero(self.data))
        return indeces

    def multiply_by_value(self, value):
        self.data = self.data * value
