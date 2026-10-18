'''
pytest -s tests/database/test_query.py
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

    def test_parse_data(self):
        records = [('a', 1), ('b', 2)]
        c = BuildComplex(True)
        c.insert_records(records, table_name, cols)

        q = QueryComplex(True)
        query = f'select name, value from {table_name} where name in '
        res = q.parse_data(query, ('a', 'b'))
        assert sorted(res, key=lambda i: i['name']) == [
            {'name': 'a', 'value': 1},
            {'name': 'b', 'value': 2},
        ]

    def test_unique_values(self):
        records = [('a', 1), ('b', 2), ('b', 3)]
        c = BuildComplex(True)
        c.insert_records(records, table_name, cols)

        q = QueryComplex(True)
        res = q.unique_values(table_name, 'name')
        assert sorted(res) == ['a', 'b']

    def test_next_id(self):
        q = QueryComplex(True)
        # empty table starts from 1
        res = q.next_id(table_name, 'value')
        assert res == 1

        records = [('a', 1), ('b', 2)]
        c = BuildComplex(True)
        c.insert_records(records, table_name, cols)
        res = q.next_id(table_name, 'value')
        assert res == 3

    def test_max_value(self):
        q = QueryComplex(True)
        res = q.max_value(table_name, 'value')
        assert res is None

        records = [('a', 1), ('b', 2)]
        c = BuildComplex(True)
        c.insert_records(records, table_name, cols)
        res = q.max_value(table_name, 'value')
        assert res == 2

    def test_table_names(self):
        q = QueryComplex(True)
        res = q.table_names()
        assert isinstance(res, tuple)
        assert table_name in res

    def test_combo_ids(self):
        # the query must be valid SQL even when nothing matches
        q = QueryComplex(True)
        res = q.combo_ids(999999, 'AbAg')
        assert res == []

    def test_delete_pdb_ids(self):
        table2 = 'tmp_chain'
        d = DeleteDb(True)
        d.config['database'] = DB
        d.drop_table(table2)

        c = CreateDb(True)
        c.config['database'] = DB
        c.create_table(table2, 'pdb_id varchar(10), chain_type varchar(5)')
        try:
            records = [('1abc', 'H'), ('2abc', 'x'), ('3abc', 'x')]
            b = BuildComplex(True)
            b.insert_records(records, table2, ['pdb_id', 'chain_type'])

            q = QueryComplex(True)
            res = q.delete_pdb_ids(table2)
            assert tuple(sorted(res)) == ('2abc', '3abc')
        finally:
            d.drop_table(table2)