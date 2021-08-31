from rdkit import Chem
import pandas as pd
import numpy as np
import os
import argparse
import logging
import math
from grid_TG import GRID
import json
import pathlib
import random
from grid_datareadin import GRIDDATA


class SCOODER:
    def __init__(self, args):
        self.logger = logging.getLogger("SCOODER")
        self.arguments = args
        self.modus = self.arguments["modus"]
        self.path = pathlib.Path(__file__).parent.absolute()
        self.start_dir = os.getcwd()
        self.allperc = [0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 0.99]
        if self.modus == "nosc" or self.modus == "query":
            self.name = os.path.basename(self.arguments["dockingposes"]).split(".")[0]
            os.mkdir(self.name)
            os.chdir(self.name)
            self.factor_cube = 6 + 12 * (1 / math.sqrt(2)) + 8 * (1 / math.sqrt(3))
            self.stauch = 5 * (
                1 / self.factor_cube
            )  # between 0 and 1 1 no stauching and close to zero no shining outside
            self.size_gridspace = self.arguments["gridsize"]  # angstroem
            self.shineout = (
                self.factor_cube * self.stauch / (self.factor_cube * self.stauch + 1)
            )
            self.num_build = 0
            self.dic_features = {}
            self.read_features()
        if self.modus == "query":
            self.print_pdbs = self.arguments["pdb"]
        if self.modus == "analysis":
            self.max_values = {}
        self.init_logs()
        self.logger.info("Used arguments: " + str(self.arguments))

    def init_logs(self):
        """
        Initiates the logs
        Nothing much to see here
        :return:
        """
        ### Initiate logging using multiple logging files
        self.logger.setLevel(logging.DEBUG)
        formatter = logging.Formatter(
            "%(asctime)s - %(name)s - %(levelname)s - %(message)s", "%d-%m-%Y %H:%M:%S"
        )
        ch = logging.StreamHandler()
        ch.setLevel(logging.ERROR)
        ch.setFormatter(formatter)
        self.logger.addHandler(ch)
        fh = logging.FileHandler(os.path.join(os.getcwd(), "scooder_debug.log"))
        fh2 = logging.FileHandler(os.path.join(os.getcwd(), "scooder_info.log"))
        fh.setLevel(logging.DEBUG)
        fh.setFormatter(formatter)
        self.logger.addHandler(fh)
        fh2.setLevel(logging.INFO)
        fh2.setFormatter(formatter)
        self.logger.addHandler(fh2)

    def read_features(self):
        """
        This function reads the features from feature_dictionaries.json
        It also checks if the command line parameter exists in the json and if not raises a RuntimeError
        :return:
        """
        self.logger.debug("read_features function called")
        f = open(os.path.join(self.path, "feature_dictionaries.json"))
        data = json.load(f)
        f.close()
        features_used = {}
        if self.arguments["features"] in data:
            features_prelim = data[self.arguments["features"]]
        else:
            self.logger.error("Feature was not found in features_dictionaries.json")
            raise RuntimeError("Feature was not found in features_dictionaries.json")
        for bla in features_prelim:
            features_used[bla] = Chem.MolFromSmarts(features_prelim[bla])
        self.dic_features = features_used
        return True

    def create_parameter_dictionary(self, grid_dictionary):
        """
        This function creates a dictionary with important parameters that all grids in the scooder run must have in common
        It chooses a random grid because I dont want to hardcode one and all grids MUST have these in common
        Writes the parameter dictionary in a json file
        :return:
        """
        self.logger.debug("create_parameter_dictionary function called")
        self.parameter_dictionary = {}
        self.parameter_dictionary["basename"] = self.name
        key, val = random.choice(list(grid_dictionary.items()))
        self.parameter_dictionary["gridsize"] = grid_dictionary[key].size_gridspace

        self.parameter_dictionary["min_x"] = grid_dictionary[key].min_x
        self.parameter_dictionary["max_x"] = grid_dictionary[key].max_x
        self.parameter_dictionary["min_y"] = grid_dictionary[key].min_y
        self.parameter_dictionary["max_y"] = grid_dictionary[key].max_y
        self.parameter_dictionary["min_z"] = grid_dictionary[key].min_z
        self.parameter_dictionary["max_z"] = grid_dictionary[key].max_z
        self.parameter_dictionary["features"] = self.arguments["features"]

        with open("config.json", "w") as json_file:
            json.dump(self.parameter_dictionary, json_file)

    def output_data_nosc(self, grid_dictionary):
        """
        This outputs the data inside the grids as numpy arrays
        This is only used in the nosc runs to save the grids
        :return:
        """
        self.logger.debug("output_data_nosc function called")
        for grid in grid_dictionary:
            np.save(grid, grid_dictionary[grid].data)

    def read_mols(self, file):
        """
        This function uses the SDMolSupplier from rdkit to read the molecules from the dockingposes file
        Therefore, a SD File is required and no other format will be accepted
        Also checks if more than 0 molecules have been succesfully read
        Continues if mol is None, but doesnt filter mols in any other way, later there will be a heavy atom filter!
        :return:
        """
        self.logger.debug("read_mols function called")
        os.chdir("..")
        self.counter_full = 0
        suppl = Chem.SDMolSupplier(file)
        for mol in suppl:
            if mol is None:
                continue
            self.counter_full += 1
        if self.counter_full <= 0:
            logging.warning("No Molecules read")
            raise RuntimeError
        os.chdir(self.name)
        return suppl

    def transform_mols(self, suppl):
        """
        This function transforms the raw molecules into the raw feature data that we want
        Here a heavy atom filter is employed of 5 heavy atoms to filter very small molecules
        If your docking is correct, this should not bother you, if you dock very small fragments this might be a problem
        :return:
        """
        self.logger.debug("transform_mols function called")
        feature_data = []
        i = 0
        j = 0
        for mol in suppl:
            if mol is None or mol.GetNumAtoms() < 5:
                j = j + 1
                continue
            k = 0
            for feature in self.dic_features:
                if mol.HasSubstructMatch(self.dic_features[feature]):
                    b = mol.GetSubstructMatches(self.dic_features[feature])
                    c = [d[0] for d in b]
                    e = mol.GetProp("_Name")
                    feature_data.append([i, feature, c, e, j])
                k = k + 1
            i = i + 1
            j = j + 1
        return pd.DataFrame(feature_data)

    def extract_information(self, features_raw, suppl):
        """
        This function converts the raw feature information into more usable features
        Structure : molecule_id, interaction_id, atomnumber, x, y, z
        :return:
        """
        self.logger.debug("extract_information function called")
        dataframe_full = []

        for tup in features_raw.itertuples():
            atomnumbers = tup[3]
            molecule_idoi = tup[1]
            molecule_idoi_suppl = tup[5]
            interaction_type = tup[2]
            for atomnumber in atomnumbers:
                pos = (
                    suppl[molecule_idoi_suppl]
                    .GetConformer()
                    .GetAtomPosition(atomnumber)
                )
                dataframe_full.append(
                    [
                        molecule_idoi + 1,
                        interaction_type,
                        atomnumber + 1,
                        pos.x,
                        pos.y,
                        pos.z,
                    ]
                )
        features_processed = pd.DataFrame(
            dataframe_full,
            columns=["Molecule_ID", "Interaction_ID", "Atomnumber", "X", "Y", "Z"],
        )
        num_build = features_processed["Molecule_ID"].max()
        return [num_build, features_processed]

    def color_dataframe(self, features_processed):
        """
        This function converts the processed data into the data necessary to introduce to the grid
        For that, we need to convert the xyz coordinates into the used indices in the numpy array
        First we divide by the gridsize, then we have to normalise this by the minimum in each dimension
        Due to the shift we do in the gridclass, we have to resemble this shift here (+1)
        :return:
        """
        self.logger.debug("color_dataframe function called")
        color_df = []

        for index, row in features_processed.iterrows():
            x_vector = int(row["X"] / self.size_gridspace)
            y_vector = int(row["Y"] / self.size_gridspace)
            z_vector = int(row["Z"] / self.size_gridspace)
            color_df.append(
                [
                    row["Interaction_ID"],
                    x_vector,
                    y_vector,
                    z_vector,
                    int(row["Molecule_ID"]),
                ]
            )
        color_df = pd.DataFrame(
            color_df,
            columns=[
                "Interaction Type",
                "x_vector",
                "y_vector",
                "z_vector",
                "Molecule ID",
            ],
        )
        color_df["x_vector"] = color_df["x_vector"] - min(color_df["x_vector"]) + 1
        color_df["y_vector"] = color_df["y_vector"] - min(color_df["y_vector"]) + 1
        color_df["z_vector"] = color_df["z_vector"] - min(color_df["z_vector"]) + 1
        return color_df

    def make_grids(self, features_processed, color_df):
        """
        This function actually creates the grid objects from the gathered data
        We create one grid per SMARTS pattern we use in our features list
        To initiate our grid we need to feed it the processed features
        After that, we iterate over the color_df and alter the grid each time with the respective change
        If the modus is traditional, we also calculate the necessary information for later on the fly
        :return:
        """
        self.logger.debug("make_grids function called")
        grid_dictionary = {}

        i = 0
        for j in self.dic_features:

            # Now we color it
            molecule_number = 1

            name_grid = str(j) + "grid"

            grid_conversion = GRID(
                name_grid, features_processed, self.size_gridspace, self.stauch
            )
            for index, row in color_df.iterrows():
                if row[0] == j:  # This checks for the type of interaction
                    grid_conversion.alter_grid(
                        row[1], row[2], row[3]
                    )  # This colors the grid

            grid_dictionary[j] = grid_conversion
            i += 1
        return grid_dictionary

    def calc_all_nonzero_indices(self, grid_dictionary):
        """
        This function calculates the volume that is used by all features combined
        This is necessary in the trad mode because we use this number for normalisation
        :return:
        """
        self.logger.debug("calc_all_nonzero_indices function called")
        all_indices_nonzero = []
        for grid in grid_dictionary:
            new_arrays = grid_dictionary[grid].indeces_nonzero()
            for indice in new_arrays:
                all_indices_nonzero.append(list(indice))

        unique_indices_nonzero = [
            list(x) for x in set(tuple(x) for x in all_indices_nonzero)
        ]
        self.volume_used_allmethods = (
            len(unique_indices_nonzero) * self.size_gridspace ** 3
        )
        return True

    def write_pdbs(self, grid_dictionary, percentage=0):
        """
        This function actually just writes the pdbs of the given percentage
        In the trad mode, this percentage will be the user-provided percentage
        In other modes we iterate over all percentages in the self.allperc list
        :param percentage: either given or all percentages in self.allperc used
        :return:
        """
        self.logger.debug("write_pdbs function called")
        if not os.path.isdir("grids_pdbs"):
            os.mkdir("grids_pdbs")
        os.chdir("grids_pdbs")
        if percentage == 0:
            if self.modus == "analysis" or self.modus == "query":
                percentage = self.allperc
            else:
                percentage = [self.allperc[0]]
        for grid in grid_dictionary:
            for perc in percentage:
                grid_dictionary[grid].save_PDB(perc)
        os.chdir("..")

    def write_density(self, grid_dictionary):
        """
        This function writes out the densities of all grids
        Due to the fact that this always exports ALL information, no percentage necessary
        :return:
        """
        self.logger.debug("write_density function called")
        if not os.path.isdir("grids_densities"):
            os.mkdir("grids_densities")
        os.chdir("grids_densities")
        for grid in grid_dictionary:
            grid_dictionary[grid].save_grid_density()
            if self.arguments["norm_max"]:
                maximum = grid_dictionary[grid].max_grid()
                grid_dictionary[grid].multiply_by_value(1 / maximum)
                grid_dictionary[grid].save_grid_density("maximum")
                grid_dictionary[grid].multiply_by_value(maximum)
            if self.arguments["norm_sum"]:
                summe = grid_dictionary[grid].sum_grid()
                grid_dictionary[grid].multiply_by_value(1 / summe)
                grid_dictionary[grid].save_grid_density("summe")
                grid_dictionary[grid].multiply_by_value(summe)
            if self.arguments["norm_manual"] != None:
                grid_dictionary[grid].multiply_by_value(
                    1 / self.arguments["norm_manual"]
                )
                grid_dictionary[grid].save_grid_density("manual")
                grid_dictionary[grid].multiply_by_value(self.arguments["norm_manual"])

        os.chdir("..")

    def read_configs(self):
        """
        This function is used in the analysis mode and reads the config file of all folders given to analyse
        It calculates the minmimal and maximal points of the grid that needs to be constructed
        In addition it makes some sanity checks if all grids have been done with the identical features and gridsize
        :return:
        """
        self.logger.debug("read_configs function called")
        i = 0
        for arg in self.arguments["folders"]:
            config_path = os.path.join(arg, "config.json")
            with open(config_path, "r") as f:
                parameters = json.load(f)
                f.close()
            if i == 0:
                self.max_values["min_x"] = parameters["min_x"]
                self.max_values["min_y"] = parameters["min_y"]
                self.max_values["min_z"] = parameters["min_z"]

                self.max_values["max_x"] = parameters["max_x"]
                self.max_values["max_y"] = parameters["max_y"]
                self.max_values["max_z"] = parameters["max_z"]
                self.max_values["gridsize"] = parameters["gridsize"]
                self.max_values["features"] = parameters["features"]
            else:
                if parameters["min_x"] < self.max_values["min_x"]:
                    self.max_values["min_x"] = parameters["min_x"]
                if parameters["min_y"] < self.max_values["min_y"]:
                    self.max_values["min_y"] = parameters["min_y"]
                if parameters["min_z"] < self.max_values["min_z"]:
                    self.max_values["min_z"] = parameters["min_z"]
                if parameters["max_x"] > self.max_values["max_x"]:
                    self.max_values["max_x"] = parameters["max_x"]
                if parameters["max_y"] > self.max_values["max_y"]:
                    self.max_values["max_y"] = parameters["max_y"]
                if parameters["max_z"] > self.max_values["max_z"]:
                    self.max_values["max_z"] = parameters["max_z"]
                if (
                    self.max_values["gridsize"] != parameters["gridsize"]
                    or self.max_values["features"] != parameters["features"]
                ):
                    self.logger.error("Grids with different gridsizes/features given")
                    self.logger.error("This was not intended. Script exiting")
                    raise BaseException
            i += 1
        self.arguments["features"] = self.max_values["features"]
        return True

    def combine_dockings(self):
        """
        This function is used in the analysis mode
        It now reads all the grids in the different folders and combines them into one
        :return:
        """
        self.logger.debug("combine_dockings function called")
        grid_dictionary = {}
        for feature in self.dic_features:
            grid_dictionary[feature] = GRIDDATA(
                self.max_values["min_x"],
                self.max_values["min_y"],
                self.max_values["min_z"],
                self.max_values["max_x"],
                self.max_values["max_y"],
                self.max_values["max_z"],
                self.max_values["gridsize"],
                str(feature),
            )
        for arg in self.arguments["folders"]:
            config_path = os.path.join(arg, "config.json")
            with open(config_path, "r") as f:
                parameters = json.load(f)
            f.close()
            for feature in self.dic_features:
                filename = feature + ".npy"
                path_to_file = os.path.join(arg, filename)
                herbert = np.load(path_to_file)
                grid_dictionary[feature].add_np_to_grid(herbert, parameters)
        return grid_dictionary

    def run(self):
        """
        Run function
        Run function runs
        Runs well true
        Runs bad false
        Runny thing
        :return:
        """
        self.logger.debug("run function called")
        if self.modus == "nosc":
            suppl = self.read_mols(self.arguments["dockingposes"])
            feature_data = self.transform_mols(suppl)
            infos = self.extract_information(feature_data, suppl)
            num_build = infos[0]
            features_processed = infos[1]
            color_df = self.color_dataframe(features_processed)
            grid_dictionary = self.make_grids(features_processed, color_df)
            self.create_parameter_dictionary(grid_dictionary)
            self.output_data_nosc(grid_dictionary)
            return True
        elif self.modus == "analysis":
            if self.read_configs():
                if self.read_features():
                    grid_dictionary = self.combine_dockings()
                    if self.arguments["pdb"]:
                        self.write_pdbs(grid_dictionary)
                    if self.arguments["density"]:
                        self.write_density(grid_dictionary)
                    return True
        elif self.arguments["modus"] == "query":
            suppl = self.read_mols(self.arguments["dockingposes"])
            feature_data = self.transform_mols(suppl)
            infos = self.extract_information(feature_data, suppl)
            num_build = infos[0]
            features_processed = infos[1]
            color_df = self.color_dataframe(features_processed)
            grid_dictionary = self.make_grids(features_processed, color_df)
            if self.arguments["pdb"]:
                self.write_pdbs(grid_dictionary)
            if self.arguments["density"]:
                self.write_density(grid_dictionary)
            return True
        else:
            raise RuntimeError


def main():
    args = parser.parse_args()
    args_dict = vars(args)
    Herbert = SCOODER(args_dict)
    Herbert.run()


if __name__ == "__main__":
    # Parser and subparser added
    parser = argparse.ArgumentParser(description="Provide settings for SCOODER")
    subparsers = parser.add_subparsers(
        dest="modus", help="Different Modes of operation for SCOODER"
    )
    # Nosc Mode
    parser_nosc = subparsers.add_parser(
        "nosc", help="Mode for saving of grids to disk to analyse them at a later point"
    )
    parser_nosc.add_argument(
        "--dockingposes", "-f", help="File of the docking poses as SDF", required=True
    )
    parser_nosc.add_argument(
        "--gridsize",
        "-g",
        type=float,
        help="Gridsize as float in Angstroem",
        required=True,
    )
    parser_nosc.add_argument(
        "--features",
        help="What feature collection found in feature_dictionaries.json do you want to use? \n"
        "E.g. normal ,onlyrings, onlyhb, extended",
    )

    # Analysis Mode
    parser_analysis = subparsers.add_parser(
        "analysis", help="Analysis Mode to analyse one or multiple nosc runs"
    )
    parser_analysis.add_argument(
        "--folders",
        nargs="+",
        help="The names of the folders on which analysis shall be performed",
        required=True,
    )
    parser_analysis.add_argument("--pdb", action="store_true")
    parser_analysis.add_argument("--density", action="store_true")
    parser_analysis.add_argument("--norm_max", action="store_true")
    parser_analysis.add_argument("--norm_sum", action="store_true")
    parser_analysis.add_argument(
        "--norm_manual",
        type=float,
        help="Introduce a custom value by which the grid is divided before saving the densities",
    )

    # Query Mode
    parser_query = subparsers.add_parser(
        "query", help="Query Mode to analyse a single docking"
    )
    parser_query.add_argument(
        "--dockingposes", "-f", help="File of the docking poses as SDF", required=True
    )
    parser_query.add_argument(
        "--features",
        help="What feature collection found in feature_dictionaries.json do you want to use? \n"
        "E.g. normal ,onlyrings, onlyhb, extended",
    )
    parser_query.add_argument(
        "--gridsize",
        "-g",
        type=float,
        help="Gridsize as float in Angstroem",
        required=True,
    )
    parser_query.add_argument("--pdb", action="store_true")
    parser_query.add_argument("--density", action="store_true")
    parser_query.add_argument("--norm_max", action="store_true")
    parser_query.add_argument("--norm_sum", action="store_true")
    parser_query.add_argument(
        "--norm_manual",
        type=float,
        help="Introduce a custom value by which the grid is divided before saving the densities",
        default=None,
    )
    main()
