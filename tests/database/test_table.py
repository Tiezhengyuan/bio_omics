'''
pytest -s tests/database/test_table.py
'''
from ..helper import *
from src.bioomics import CreateDb, DeleteDb, QueryDb

DB = 'COMPLEX2'
table_name = 'tmp'


@ddt
class TestTable(TestCase):

    def setUp(self) -> None:
        d = DeleteDb()
        d.config['database'] = DB
        d.drop_table(table_name)
       
    def test_crud(self):
        # create
        cols_str = """
            tmp_id int auto_increment primary key,
            name varchar(10),
            value INT
        """
        c = CreateDb(True)
        c.config['database'] = DB
        res = c.create_table(table_name, cols_str)
        assert res

        # table exists
        q = QueryDb(True)
        q.config['database'] = DB
        res = q.table_exists(table_name)
        assert res

        # rows
        res = q.num_rows(table_name)
        assert res == 0

        # drop table
        d = DeleteDb(True)
        d.config['database'] = DB
        res = d.drop_table(table_name)
        assert res

        # table exists
        res = q.table_exists(table_name)
        assert not res

    def test_invalid_table(self):
        c = CreateDb(True)
        c.config['database'] = DB
        res = c.create_table(table_name, '')
        assert not res

        # table exists
        q = QueryDb(True)
        q.config['database'] = DB
        res = q.table_exists(table_name)
        assert not res

        # rows
        res = q.num_rows(table_name)
        assert res == -1

        # drop table
        d = DeleteDb(True)
        d.config['database'] = DB
        res = d.drop_table(table_name)
        assert not res