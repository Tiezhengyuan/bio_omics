'''
pytest -s tests/database/test_crud.py
'''
from ..helper import *
from src.bioomics import CreateDb, DeleteDb, QueryComplex, BuildComplex

DB = 'COMPLEX2'
table_name = 'tmp'
cols = ['name', 'value']

@ddt
class TestQuery(TestCase):

    def setUp(self):
        # drop the leftover table from a previous failed run
        d = DeleteDb(True)
        d.config['database'] = DB
        d.drop_table(table_name)

        c = CreateDb(True)
        c.config['database'] = DB
        cols_str = """
            name varchar(10),
            value INT
        """
        c.create_table(table_name, cols_str)

    def tearDown(self) -> None:
        d = DeleteDb(True)
        d.config['database'] = DB
        d.drop_table(table_name)

    def test_list_data(self):
        query = f'select * from {table_name}'

        records = []
        c = BuildComplex(True)
        c.insert_records(records, table_name, cols)
        res = QueryComplex().list_data(query)
        assert res == []

        records = [('a',1), ('b',1)]
        c = BuildComplex(True)
        c.insert_records(records, table_name, cols)
        res = QueryComplex().list_data(query)
        assert res == [{'name': 'a', 'value': 1},
            {'name': 'b', 'value': 1},]
        res = QueryComplex().list_data(query, True)
        assert list(res) == cols