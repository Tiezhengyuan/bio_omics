from copy import deepcopy
from datetime import datetime
import json
import os
import pandas as pd

from ddt import ddt, data, unpack
from unittest import TestCase, mock, skip

# store example data for unit and local testing
DIR_DATA = os.path.join(os.path.dirname(__file__), 'data')
