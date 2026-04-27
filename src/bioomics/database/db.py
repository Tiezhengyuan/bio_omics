'''
configurations for mysql.connector.connection
'''

import os

class Db:

    def __init__(self, verbose:bool):
        self.verbose = verbose
        self.config = {
            'host': os.getenv('MYSQL_HOST'),
            'user': os.getenv('MYSQL_USER'),
            'password': os.getenv('MYSQL_PW'),
            'use_pure': True,
            'buffered': True,
            'autocommit': True,    # avoid transaction issues
            'connection_timeout': 60,
            'pool_size': 1,        # Single connection for stability
        }