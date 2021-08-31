from unittest import TestCase
from SCOODER import SCOODER
import os
import rdkit
import pandas as pd
from ast import literal_eval
import math
import shutil
import pathlib


class TestScooder(TestCase):
    def setUp(self):
        if os.path.isdir("Testrun"):
            shutil.rmtree("Testrun")
        os.mkdir("Testrun")
        os.chdir("Testrun")
        self.path = pathlib.Path(__file__).parent.absolute()
        self.args_nosc = {
            "modus": "nosc",
            "dockingposes": "./../testfiles/05_moe_prepared.sdf",
            "gridsize": 0.5,
            "features": "normal",
        }
        self.args_analysis = {
            "modus": "analysis",
            "folders": ["adrb2_fred", "adrb2_hybrid"],
            "pdb": True,
            "density": True,
            "norm_sum": False,
            "norm_max": False,
            "norm_manual": None,
        }
        self.args_analysis2 = {
            "modus": "analysis",
            "folders": ["adrb2_fred"],
            "pdb": True,
            "density": True,
            "norm_sum": False,
            "norm_max": False,
            "norm_manual": None,
        }
        self.args_query = {
            "modus": "query",
            "dockingposes": "./../testfiles/05_moe_prepared.sdf",
            "gridsize": 0.5,
            "pdb": True,
            "features": "normal",
            "density": True,
            "norm_sum": False,
            "norm_max": False,
            "norm_manual": None,
        }
        self.args_analysis3 = {
            "modus": "analysis",
            "folders": ["adrb2_fred"],
            "pdb": True,
            "density": True,
            "norm_sum": True,
            "norm_max": True,
            "norm_manual": None,
        }
        self.args_analysis4 = {
            "modus": "analysis",
            "folders": ["adrb2_fred"],
            "pdb": True,
            "density": True,
            "norm_sum": True,
            "norm_max": True,
            "norm_manual": 400,
        }

    def tearDown(self):
        os.chdir(self.path)
        if os.path.isdir("Testrun"):
            shutil.rmtree("Testrun")

    def test_init(self):
        herbert = SCOODER(self.args_nosc)
        self.assertEqual("05_moe_prepared", herbert.name)
        self.assertEqual(True, os.path.isfile("scooder_debug.log"))
        self.assertEqual(True, os.path.isfile("scooder_info.log"))

    def test_read_mols(self):
        herbert = SCOODER(self.args_nosc)
        suppl = herbert.read_mols(self.args_nosc["dockingposes"])
        self.assertIsInstance(suppl, rdkit.Chem.rdmolfiles.SDMolSupplier)
        self.assertEqual(herbert.counter_full, 510)

    def test_features_readin(self):
        herbert = SCOODER(self.args_nosc)
        herbert.read_features()
        self.assertEqual(7, len(herbert.dic_features))
        self.args_nosc["features"] = "onlyhb"
        herbert = SCOODER(self.args_nosc)
        herbert.read_features()
        self.assertEqual(2, len(herbert.dic_features))
        self.args_nosc["features"] = "onlyrings"
        herbert = SCOODER(self.args_nosc)
        herbert.read_features()
        self.assertEqual(2, len(herbert.dic_features))
        self.args_nosc["features"] = "extended"
        herbert = SCOODER(self.args_nosc)
        herbert.read_features()
        self.assertEqual(8, len(herbert.dic_features))

    def test_transform_mols(self):
        herbert = SCOODER(self.args_nosc)
        suppl = herbert.read_mols(self.args_nosc["dockingposes"])
        features_raw = herbert.transform_mols(suppl)
        to_compare = pd.read_csv(
            "./../../testfiles/features_raw_moeprep.csv",
            index_col=0,
            header=0,
            names=["A", "B", "C", "D", "E"],
            converters={"C": literal_eval},
        )
        features_raw.columns = ["A", "B", "C", "D", "E"]
        to_compare_extra = to_compare["C"].to_list()
        standard_extra = features_raw["C"].to_list()
        to_compare = to_compare.drop(["C"], axis=1)
        standard = features_raw.drop(["C"], axis=1)
        self.assertEqual(True, to_compare.compare(standard).empty)
        for i in range(len(standard_extra)):
            interesting_list = list(standard_extra[i])
            interesting_compare = list(to_compare_extra[i])
            for j in range(len(interesting_list)):
                self.assertEqual(interesting_list[j], interesting_compare[j])

    def test_extract_information(self):
        herbert = SCOODER(self.args_nosc)
        suppl = herbert.read_mols(self.args_nosc["dockingposes"])
        infos = herbert.extract_information(herbert.transform_mols(suppl), suppl)
        num_build = infos[0]
        features_processed = infos[1]
        datatypes = {
            "Molecule_ID": int,
            "Interaction_ID": str,
            "Atomnumber": int,
            "X": float,
            "Y": float,
            "Z": float,
        }
        to_compare = pd.read_csv(
            os.path.join("..", "..", "testfiles", "features_processed_moreprep.csv"),
            header=0,
            index_col=0,
            dtype=datatypes,
        )
        colnames = ["Molecule_ID", "Interaction_ID", "Atomnumber", "X", "Y", "Z"]
        self.assertEqual(True, self.compare_dataframes(features_processed, to_compare))

    def compare_dataframes(self, dataframe_a, dataframe_b):
        colnames = dataframe_a.columns.values
        print(colnames)
        for bla in colnames:
            series_a = dataframe_a[bla].squeeze()
            series_b = dataframe_b[bla].squeeze()
            comparison = series_a.compare(series_b)
            if not comparison.empty:
                for i in range(len(series_a)):
                    if not math.isclose(series_a[i], series_b[i]):
                        return False
        return True

    def test_color_dataframe(self):
        herbert = SCOODER(self.args_nosc)
        suppl = herbert.read_mols(self.args_nosc["dockingposes"])
        infos = herbert.extract_information(herbert.transform_mols(suppl), suppl)
        num_build = infos[0]
        features_processed = infos[1]
        color_df = herbert.color_dataframe(features_processed)
        to_compare = pd.read_csv(
            filepath_or_buffer=os.path.join(
                "..", "..", "testfiles", "color_df_moeprep.csv"
            ),
            header=0,
            index_col=0,
        )
        self.assertEqual(True, color_df.compare(to_compare).empty)

    def test_make_grids(self):
        herbert = SCOODER(self.args_nosc)
        suppl = herbert.read_mols(self.args_nosc["dockingposes"])
        infos = herbert.extract_information(herbert.transform_mols(suppl), suppl)
        num_build = infos[0]
        features_processed = infos[1]
        color_df = herbert.color_dataframe(features_processed)
        grid_dictionary = herbert.make_grids(features_processed, color_df)
        for feature in herbert.dic_features:
            sum = grid_dictionary[feature].sum_grid()
            if feature == "donor":
                self.assertAlmostEqual(sum, 1029)
            elif feature == "acceptor":
                self.assertAlmostEqual(sum, 1715)
            elif feature == "aromatic":
                self.assertAlmostEqual(sum, 7014)
            elif feature == "halogen":
                self.assertAlmostEqual(sum, 167)
            elif feature == "basic":
                self.assertAlmostEqual(sum, 724)
            elif feature == "acidic":
                self.assertAlmostEqual(sum, 0)
            elif feature == "aliphatic_ring":
                self.assertAlmostEqual(sum, 959)

    def test_calc_all_nonzero_indices(self):
        herbert = SCOODER(self.args_nosc)
        suppl = herbert.read_mols(self.args_nosc["dockingposes"])
        infos = herbert.extract_information(herbert.transform_mols(suppl), suppl)
        num_build = infos[0]
        features_processed = infos[1]
        color_df = herbert.color_dataframe(features_processed)
        grid_dictionary = herbert.make_grids(features_processed, color_df)
        herbert.calc_all_nonzero_indices(grid_dictionary)
        self.assertEqual(831.5, herbert.volume_used_allmethods)

    def test_write_pdbs(self):
        herbert = SCOODER(self.args_nosc)
        suppl = herbert.read_mols(self.args_nosc["dockingposes"])
        infos = herbert.extract_information(herbert.transform_mols(suppl), suppl)
        num_build = infos[0]
        features_processed = infos[1]
        color_df = herbert.color_dataframe(features_processed)
        grid_dictionary = herbert.make_grids(features_processed, color_df)
        herbert.calc_all_nonzero_indices(grid_dictionary)
        herbert.write_pdbs(grid_dictionary)
        count_pdbs = 0
        for file in os.listdir("grids_pdbs"):
            if file.endswith(".pdb"):
                count_pdbs += 1
        self.assertEqual(count_pdbs, len(herbert.dic_features))

    def test_write_density(self):
        herbert = SCOODER(self.args_query)
        suppl = herbert.read_mols(self.args_query["dockingposes"])
        infos = herbert.extract_information(herbert.transform_mols(suppl), suppl)
        num_build = infos[0]
        features_processed = infos[1]
        color_df = herbert.color_dataframe(features_processed)
        grid_dictionary = herbert.make_grids(features_processed, color_df)
        herbert.calc_all_nonzero_indices(grid_dictionary)
        herbert.write_density(grid_dictionary)
        count_pdbs = 0
        for file in os.listdir("grids_densities"):
            if file.endswith(".dx"):
                count_pdbs += 1
        self.assertEqual(count_pdbs, len(herbert.dic_features))

    def test_create_parameter_dictionary(self):
        herbert = SCOODER(self.args_nosc)
        suppl = herbert.read_mols(self.args_nosc["dockingposes"])
        infos = herbert.extract_information(herbert.transform_mols(suppl), suppl)
        num_build = infos[0]
        features_processed = infos[1]
        color_df = herbert.color_dataframe(features_processed)
        grid_dictionary = herbert.make_grids(features_processed, color_df)
        herbert.create_parameter_dictionary(grid_dictionary)
        self.assertEqual(True, os.path.isfile("config.json"))

    def test_output_data_nosc(self):
        herbert = SCOODER(self.args_nosc)
        suppl = herbert.read_mols(self.args_nosc["dockingposes"])
        infos = herbert.extract_information(herbert.transform_mols(suppl), suppl)
        num_build = infos[0]
        features_processed = infos[1]
        color_df = herbert.color_dataframe(features_processed)
        grid_dictionary = herbert.make_grids(features_processed, color_df)
        herbert.output_data_nosc(grid_dictionary)
        self.assertEqual(9, len(os.listdir()))

    def test_nosc_run(self):
        herbert = SCOODER(self.args_nosc)
        self.assertEqual(True, herbert.run())

    def test_read_configs(self):
        shutil.copytree(
            os.path.join(self.path, "testfiles", "adrb2_fred"), "adrb2_fred"
        )
        shutil.copytree(
            os.path.join(self.path, "testfiles", "adrb2_hybrid"), "adrb2_hybrid"
        )
        herbert = SCOODER(self.args_analysis)
        herbert.read_configs()
        self.assertEqual(-3.5137, herbert.max_values["min_x"])
        self.assertEqual(12.093300000000001, herbert.max_values["max_x"])
        self.assertEqual(-3.4751, herbert.max_values["min_y"])
        self.assertEqual(18.078, herbert.max_values["max_y"])
        self.assertEqual(44.2456, herbert.max_values["min_z"])
        self.assertEqual(65.30019999999999, herbert.max_values["max_z"])
        self.assertEqual(0.5, herbert.max_values["gridsize"])

    def test_combine_dockings(self):
        shutil.copytree(
            os.path.join(self.path, "testfiles", "adrb2_fred"), "adrb2_fred"
        )
        shutil.copytree(
            os.path.join(self.path, "testfiles", "adrb2_hybrid"), "adrb2_hybrid"
        )
        herbert = SCOODER(self.args_analysis)
        herbert.read_configs()
        herbert.read_features()
        grid_dictionary = herbert.combine_dockings()
        self.assertAlmostEqual(10373.999999999993, grid_dictionary["donor"].sum_grid())
        self.assertAlmostEqual(
            17395.999999999996, grid_dictionary["acceptor"].sum_grid()
        )
        self.assertAlmostEqual(
            48909.999999999956, grid_dictionary["aromatic"].sum_grid()
        )
        self.assertAlmostEqual(3067.0, grid_dictionary["halogen"].sum_grid())
        self.assertAlmostEqual(4254.999999999999, grid_dictionary["basic"].sum_grid())
        self.assertAlmostEqual(99.0, grid_dictionary["acidic"].sum_grid())
        self.assertAlmostEqual(
            15651.999999999998, grid_dictionary["aliphatic_ring"].sum_grid()
        )
        self.assertAlmostEqual(
            110172.99999999993, grid_dictionary["everything"].sum_grid()
        )

    def test_writestuff_analysis(self):
        shutil.copytree(
            os.path.join(self.path, "testfiles", "adrb2_fred"), "adrb2_fred"
        )
        shutil.copytree(
            os.path.join(self.path, "testfiles", "adrb2_hybrid"), "adrb2_hybrid"
        )
        herbert = SCOODER(self.args_analysis)
        herbert.read_configs()
        herbert.read_features()
        grid_dictionary = herbert.combine_dockings()
        herbert.write_pdbs(grid_dictionary)
        herbert.write_density(grid_dictionary)
        self.assertEqual(80, len(os.listdir("grids_pdbs")))
        self.assertEqual(8, len(os.listdir("grids_densities")))

    def test_run_analysis(self):
        shutil.copytree(
            os.path.join(self.path, "testfiles", "adrb2_fred"), "adrb2_fred"
        )
        shutil.copytree(
            os.path.join(self.path, "testfiles", "adrb2_hybrid"), "adrb2_hybrid"
        )
        herbert = SCOODER(self.args_analysis)
        self.assertEqual(True, herbert.run())
        self.assertEqual(80, len(os.listdir("grids_pdbs")))
        self.assertEqual(8, len(os.listdir("grids_densities")))

    def test_run_query(self):
        herbert = SCOODER(self.args_query)
        self.assertEqual(True, herbert.run())
        self.assertEqual(70, len(os.listdir("grids_pdbs")))
        self.assertEqual(7, len(os.listdir("grids_densities")))

    def test_run_analysis2(self):
        shutil.copytree(
            os.path.join(self.path, "testfiles", "adrb2_fred"), "adrb2_fred"
        )
        herbert = SCOODER(self.args_analysis2)
        self.assertEqual(True, herbert.run())
        self.assertEqual(80, len(os.listdir("grids_pdbs")))
        self.assertEqual(8, len(os.listdir("grids_densities")))

    def test_run_analysis3(self):
        shutil.copytree(
            os.path.join(self.path, "testfiles", "adrb2_fred"), "adrb2_fred"
        )
        herbert = SCOODER(self.args_analysis3)
        self.assertEqual(True, herbert.run())
        self.assertEqual(80, len(os.listdir("grids_pdbs")))
        self.assertEqual(24, len(os.listdir("grids_densities")))

    def test_run_analysis4(self):
        shutil.copytree(
            os.path.join(self.path, "testfiles", "adrb2_fred"), "adrb2_fred"
        )
        herbert = SCOODER(self.args_analysis4)
        self.assertEqual(True, herbert.run())
        self.assertEqual(80, len(os.listdir("grids_pdbs")))
        self.assertEqual(32, len(os.listdir("grids_densities")))
