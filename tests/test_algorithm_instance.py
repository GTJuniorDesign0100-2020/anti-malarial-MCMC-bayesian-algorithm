import os
import sys

import numpy as np
import pandas as pd

# Add parent directory to search path, so we can import those files
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from algorithm_instance import AlgorithmInstance


def test_runs_without_error():
    example_file = os.path.join(
        os.path.dirname(__file__),
        '../Angola2017_example.xlsx')
    test = AlgorithmInstance(example_file, [])
