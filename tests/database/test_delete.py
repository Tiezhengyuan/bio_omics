'''
pytest -s tests/database/test_delete.py
'''
from ..helper import *
from src.bioomics import BuildComplex, CreateDb, DeleteComplex, DeleteDb, \
    QueryComplex

DB = 'COMPLEX2'
table_name = 'tmp'
cols = ['name', 'value']

@ddt
class TestDelete(TestCase):

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

        records = [('a', 1), ('b', 2), ('c', None)]
        b = BuildComplex(True)
        b.insert_records(records, table_name, cols)

    def tearDown(self) -> None:
        d = DeleteDb(True)
        d.config['database'] = DB
        d.drop_table(table_name)

    def test_empty_table(self):
        d = DeleteDb(True)
        d.config['database'] = DB
        res = d.empty_table(table_name)
        assert res

        res = QueryComplex(True).num_rows(table_name)
        assert res == 0

    def test_delete_rows(self):
        d = DeleteDb(True)
        d.config['database'] = DB
        res = d.delete_rows(table_name, 'name', ['a', 'b'])
        assert res

        res = QueryComplex(True).unique_values(table_name, 'name')
        assert res == ['c']

    def test_delete_null_rows(self):
        d = DeleteDb(True)
        d.config['database'] = DB
        res = d.delete_null_rows(table_name, 'value')
        assert res

        res = QueryComplex(True).unique_values(table_name, 'name')
        assert sorted(res) == ['a', 'b']

    def test_remove_null_rows(self):
        # more records than one chunk (default chunk_size = 100)
        records = [(f'r{i}', i) for i in range(250)]
        b = BuildComplex(True)
        b.insert_records(records, table_name, cols)

        dc = DeleteComplex(False)
        del_pool = [f'r{i}' for i in range(150)]
        num_del = dc.remove_null_rows(table_name, 'name', del_pool)
        assert num_del == 150

        res = QueryComplex(True).num_rows(table_name)
        assert res == 103
