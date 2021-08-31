from unittest import TestCase
from grid_TG import GRID
import pandas as pd
import numpy as np
import os


class Testgrid(TestCase):
    def test_init(self):
        random_data = pd.read_csv("./testfiles/random_data_raw_input.csv")
        random = GRID("random", random_data, 0.5, 0.2617)
        self.assertAlmostEqual(-20.499698309826396, random.min_x)
        self.assertAlmostEqual(-20.499953189332818, random.min_y)
        self.assertAlmostEqual(-20.49996776362234, random.min_z)
        self.assertAlmostEqual(20.500301690173607, random.max_x)
        self.assertAlmostEqual(20.500046810667182, random.max_y)
        self.assertAlmostEqual(20.50003223637766, random.max_z)
        self.assertAlmostEqual(41.0, random.len_x)
        self.assertAlmostEqual(41.0, random.len_y)
        self.assertAlmostEqual(41.0, random.len_z)
        self.assertEqual(83, random.num_spaces_x)
        self.assertEqual(83, random.num_spaces_y)
        self.assertEqual(83, random.num_spaces_z)
        moe_prepared_data = pd.read_csv("./testfiles/05_moe_prepared_raw_input.csv")
        moe_prepared = GRID("05_moe_prepared", moe_prepared_data, 0.5, 0.2617)
        self.assertAlmostEqual(-1.009, moe_prepared.min_x)
        self.assertAlmostEqual(-3.8489999999999998, moe_prepared.min_y)
        self.assertAlmostEqual(46.37, moe_prepared.min_z)
        self.assertAlmostEqual(9.491, moe_prepared.max_x)
        self.assertAlmostEqual(15.651, moe_prepared.max_y)
        self.assertAlmostEqual(64.37, moe_prepared.max_z)
        self.assertAlmostEqual(10.5, moe_prepared.len_x)
        self.assertAlmostEqual(19.5, moe_prepared.len_y)
        self.assertAlmostEqual(18.000000000000007, moe_prepared.len_z)
        self.assertEqual(22, moe_prepared.num_spaces_x)
        self.assertEqual(40, moe_prepared.num_spaces_y)
        self.assertEqual(38, moe_prepared.num_spaces_z)
        with self.assertRaises(RuntimeError):
            GRID("fail", random_data, 0, 0.2617)
        with self.assertRaises(RuntimeError):
            GRID("fail", random_data, 0.5, -10)
        with self.assertRaises(RuntimeError):
            GRID("fail", random_data, 0.5, 10)

    def test_alter_grid(self):
        random_data = pd.read_csv("./testfiles/random_data_raw_input.csv")
        random = GRID("random", random_data, 0.5, 0.2617)
        random.alter_grid(10, 10, 10)
        self.assertEqual(1, np.sum(random.data))
        random.alter_grid(20, 20, 20)
        self.assertEqual(2, np.sum(random.data))
        random.alter_grid(40, 40, 40)
        self.assertEqual(3, np.sum(random.data))
        random.alter_grid(82, 10, 10)
        self.assertAlmostEqual(3.7322678073143187, np.sum(random.data))
        random.alter_grid(82, 82, 10)
        self.assertAlmostEqual(4.278015495605006, np.sum(random.data))
        random.alter_grid(82, 82, 82)
        self.assertAlmostEqual(4.693271107889538, np.sum(random.data))
        with self.assertRaises(IndexError):
            random.alter_grid(2323, 23123, 12312)

    def test_sum_grid(self):
        random_data = pd.read_csv("./testfiles/random_data_raw_input.csv")
        random = GRID("random", random_data, 0.5, 0.2617)
        self.assertEqual(0, random.sum_grid())
        random.alter_grid(20, 20, 20)
        self.assertAlmostEqual(1, random.sum_grid())
        random.alter_grid(20, 20, 20)
        self.assertAlmostEqual(2, random.sum_grid())
        random.alter_grid(20, 20, 20)
        self.assertAlmostEqual(3, random.sum_grid())
        random.alter_grid(0, 0, 0)
        self.assertAlmostEqual(3.4152556122845326, random.sum_grid())
        random.alter_grid(1, 0, 0)
        self.assertAlmostEqual(3.9610033005752197, random.sum_grid())
        random.alter_grid(1, 1, 0)
        self.assertAlmostEqual(4.693271107889538, random.sum_grid())

    def test_max_grid(self):
        random_data = pd.read_csv("./testfiles/random_data_raw_input.csv")
        random = GRID("random", random_data, 0.5, 0.2617)
        self.assertEqual(0, random.max_grid())
        random.alter_grid(20, 20, 20)
        self.assertAlmostEqual(0.16667948267393465, random.max_grid())
        random.alter_grid(20, 20, 20)
        self.assertAlmostEqual(0.3333589653478693, random.max_grid())
        random.alter_grid(60, 60, 60)
        self.assertAlmostEqual(0.3333589653478693, random.max_grid())
        random.alter_grid(21, 21, 19)
        self.assertAlmostEqual(0.3585429959924404, random.max_grid())
        random.alter_grid(21, 19, 19)
        self.assertAlmostEqual(0.38372702663701147, random.max_grid())
        random.alter_grid(20, 20, 20)
        self.assertAlmostEqual(0.5504065093109461, random.max_grid())

    def test_save_pdb(self):
        random_data = pd.read_csv("./testfiles/random_data_raw_input.csv")
        random = GRID("random", random_data, 0.5, 0.2617)
        random.alter_grid(1, 1, 1)
        for i in range(1, 10):
            random.save_PDB(float(i) / 10)
            name = "random" + str(i * 10) + ".pdb"
            self.assertEqual(True, os.path.isfile(name))
            os.remove(name)

    def test_subtract_grid(self):
        random_data = pd.read_csv("./testfiles/random_data_raw_input.csv")
        random_a = GRID("randoma", random_data, 0.5, 0.2617)
        random_b = GRID("randomb", random_data, 0.5, 0.2617)
        for i in range(1, 60):
            random_a.alter_grid(i, i, i)
            random_b.alter_grid(1, 1, i)
        random_a.subtract_grid(random_b)
        self.assertAlmostEqual(115.05462339663899, random_a.sum_grid())
        self.assertAlmostEqual(0.253919523905472, random_a.max_grid())

    def test_euclidian_type_distance(self):
        random_data = pd.read_csv("./testfiles/random_data_raw_input.csv")
        random = GRID("random", random_data, 0.5, 0.2617)
        self.assertEqual(0, random.euclidian_type_distance())
        random.alter_grid(1, 1, 1)
        self.assertEqual(0.23598391369486413, random.euclidian_type_distance())
        random.alter_grid(1, 1, 1)
        self.assertEqual(0.47196782738972826, random.euclidian_type_distance())
        random.alter_grid(1, 1, 1)
        self.assertEqual(0.7079517410845922, random.euclidian_type_distance())
        random.alter_grid(1, 1, 1)
        self.assertEqual(0.9439356547794565, random.euclidian_type_distance())

    def test_volume_used(self):
        random_data = pd.read_csv("./testfiles/random_data_raw_input.csv")
        random = GRID("random", random_data, 0.5, 0.2617)
        self.assertEqual(0, random.volume_used())
        j = 1
        for i in range(1, 70, 3):
            random.alter_grid(i, i, i)
            self.assertEqual(j * 3.375, random.volume_used())
            j = j + 1

    def test_volume_used_cube(self):
        random_data = pd.read_csv("./testfiles/random_data_raw_input.csv")
        random = GRID("random", random_data, 0.5, 0.2617)
        self.assertEqual(71473.375, random.volume_cube())
        moe_prepared_data = pd.read_csv("./testfiles/05_moe_prepared_raw_input.csv")
        moe_prepared = GRID("05_moe_prepared", moe_prepared_data, 0.5, 0.2617)
        self.assertEqual(4180.0, moe_prepared.volume_cube())

    def test_indeces_nonzero(self):
        random_data = pd.read_csv("./testfiles/random_data_raw_input.csv")
        random = GRID("random", random_data, 0.5, 0.2617)
        self.assertEqual(0, len(random.indeces_nonzero()))
        random.alter_grid(1, 1, 1)
        self.assertEqual(27, len(random.indeces_nonzero()))
        random.alter_grid(1, 1, 1)
        self.assertEqual(27, len(random.indeces_nonzero()))
        random.alter_grid(1, 1, 2)
        self.assertEqual(36, len(random.indeces_nonzero()))

    def test_multiply_by_value(self):
        random_data = pd.read_csv("./testfiles/random_data_raw_input.csv")
        random = GRID("random", random_data, 0.5, 0.2617)
        random.alter_grid(1, 1, 1)
        self.assertEqual(1, random.sum_grid())
        random.multiply_by_value(3)
        self.assertEqual(3, random.sum_grid())
        random.multiply_by_value(1 / 3)
        self.assertAlmostEqual(1, random.sum_grid())
        random.multiply_by_value(0)
        self.assertAlmostEqual(0, random.sum_grid())
