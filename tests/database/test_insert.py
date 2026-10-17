'''
pytest -s tests/database/test_crud.py
'''
from ..helper import *
from src.bioomics import Db, CreateDb, DeleteDb, InsertDb, \
    QueryComplex, BuildComplex

DB = 'COMPLEX2'
table_name = 'tmp'
cols = ['name', 'value']

@ddt
class TestInsert(TestCase):

    def setUp(self):
        # drop the leftover table from a previous failed run
        d = DeleteDb(True)
        d.config['database'] = DB
        d.drop_table(table_name)

        c = CreateDb(True)
        c.config['database'] = DB
        cols_str = """
            name varchar(10) primary key,
            value INT
        """
        c.create_table(table_name, cols_str)

    def tearDown(self) -> None:
        d = DeleteDb(True)
        d.config['database'] = DB
        d.drop_table(table_name)

    def test_db_execute(self):
        b = Db(True)
        b.config['database'] = DB

        #insert one record
        query = f"""
            insert into {table_name} (name, value)
            values ('a', 10)
        """
        res = b.execute(query)
        assert res

    def test_insert(self):
        b = InsertDb(True)
        b.config['database'] = DB

        #insert one reord
        res = b.insert(['b', 20], table_name, cols)
        assert res

        q = QueryComplex(True)
        res = q.num_rows(table_name)
        assert res == 1

    def test_db_func_insert(self):
        records = [('a', 1), ('b', 2),]
        b = InsertDb(True)
        b.config['database'] = DB

        query = f"""
            insert into {table_name} (name, value)
            values (%s, %s)
        """
        func = lambda rec: [rec,]
        res = b.func_insert(records, query, func)
        assert res == (2, [])
        res = QueryComplex(True).num_rows(table_name)
        assert res == 2

    @data(
        [[('a', 1), ('b', 2),], 2, []],
        # extra column
        [[('a', 1), ('b', 2, 2),], 1, [('b', 2, 2)]],
        #duplicates
        [[('a', 1), ('a', 2),], 1, [('a', 2)]],
    )
    @unpack
    def test_insert_records(self, records, expect_succeed, expect_fail):
        b = InsertDb(True)
        b.config['database'] = DB

        good, fail = b.insert_records(records, table_name, cols)
        assert good == expect_succeed
        assert fail == expect_fail
        res = QueryComplex().num_rows(table_name)
        assert res == expect_succeed
    
    def test_insert_batch_records(self):
        records = [('a', 1), ('b', 2),]
        b = BuildComplex(True)
        res = b.insert_batch_records(records, table_name, cols)
        assert res == (2, [])
        res = QueryComplex().num_rows(table_name)
        assert res == 2