'''
pytest -s tests/database/test_update.py
'''
from ..helper import *
from src.bioomics import Db, CreateDb, DeleteDb, UpdateDb, QueryComplex

DB = 'COMPLEX2'
table_name = 'tmp'

@ddt
class TestUpdate(TestCase):

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

        b = Db(True)
        b.config['database'] = DB
        records = [('a', 1), ('b', 2)]
        for name, value in records:
            query = f"""
                insert into {table_name} (name, value)
                values ('{name}', {value})
            """
            assert b.execute(query)

    def tearDown(self) -> None:
        d = DeleteDb(True)
        d.config['database'] = DB
        d.drop_table(table_name)

    @data(
        [{'name': 'a', 'value': 10}, 1],
        # no row matches
        [{'name': 'not_exist', 'value': 10}, 0],
    )
    @unpack
    def test_execute_update(self, record, expect_rowcount):
        u = UpdateDb(True)
        u.config['database'] = DB
        query = f"""
            UPDATE {table_name} SET value = %(value)s
            WHERE name = %(name)s
        """
        res = u.execute_update(query, record)
        assert res == expect_rowcount

        q = QueryComplex(True)
        rows = q.list_data(
            f"select value from {table_name} where name = '{record['name']}'")
        if expect_rowcount:
            assert rows[0]['value'] == record['value']
        else:
            assert rows == []

    def test_execute_update_error(self):
        u = UpdateDb(True)
        u.config['database'] = DB
        query = f"""
            UPDATE {table_name} SET value = %(value)s
            WHERE name = %(name)s
        """
        # a missing key in the record raises an error
        res = u.execute_update(query, {'name': 'a'})
        assert res is None
