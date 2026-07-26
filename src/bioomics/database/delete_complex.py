from collections import defaultdict
import os
import pandas as pd

from .delete_db import DeleteDb

class DeleteComplex(DeleteDb):

    def __init__(self, verbose:bool=None):
        super().__init__(verbose)
        self.config['database'] = os.getenv('MYSQL_DB', 'COMPLEX2')