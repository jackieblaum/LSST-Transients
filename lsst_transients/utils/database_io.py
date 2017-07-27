import contextlib
import sqlite3
import pandas as pd
from pandas.io.sql import get_schema as pd_get_schema

from file_io import sanitize_filename


# Translation table for array types
# (from numpy dtype to Sqlite3 type)

array_types = {'float32': 'REAL',
               'float64': 'REAL',
               'float16': 'REAL',
               'int64': 'INTEGER',
               'int32': 'INTEGER',
               'object': 'TEXT'}


@contextlib.contextmanager
def bulk_operation(db_instance):
    """
    A context manager for bulk operations.

    Use this if you need speed in a loop where you are writing a lot to the database. This sacrifices solidity against
    crashes in exchange for speed (a lot of speed)

    :param db_instance: the instance of the database (you need to have called .connect() first)
    :return:
    """

    # Deactivate writing on disk at every iteration

    db_instance.query("PRAGMA synchronous=OFF")

    try:

        yield

    except:

        raise

    finally:

        # We need to commit before reactivating the synchronous mode, otherwise sqlite3 will fail

        db_instance.connection.commit()

        db_instance.query("PRAGMA synchronous=ON")


class SqliteDatabase(object):

    def __init__(self, database_file):

        # Sanitize the file name (expand env. variables, relative paths and so on)

        self._database_file = sanitize_filename(database_file)

        # These will keep track of the connection

        self._connection = None
        self._cursor = None

    def connect(self):
        """
        Connect to the database (and create if not existent)

        :return: None
        """

        self._connection = sqlite3.connect(sanitize_filename(self._database_file))

        self._cursor = self._connection.cursor()

    def disconnect(self, commit=True):
        """
        Disconnect from the database

        :param commit: (default: True) commit the changes if any. Use False only if you are absolutely sure that there
         is no change to be written on disk
        :return:
        """


        assert self._cursor is not None, "Connection is already closed"

        if commit:

            self._connection.commit()

        self._connection.close()

        self._connection = None
        self._cursor = None

    @property
    def cursor(self):
        """
        Returns the active cursor (or None if there is no active connection)

        :return:
        """

        return self._cursor

    @property
    def connection(self):
        """
        Returns the active connection (or None if there is no connection)

        :return:
        """

        return self._connection

    def query(self, query, *args, **kwargs):
        """
        Execute query.

        :param query: the query to be executed
        :param args: arguments as in the method .execute of the cursor object in sqlite3
        :param kwargs: arguments as in the method .execute of the cursor object in sqlite3, except:
        :param fetch: whether to fetch all results or not. If False, you need to fetch the results yourself using
        the "cursor" property and one of the fetch methods (fetchall, fetchone, or the iterator interface).
        :param commit: whether to commit immediately or not.
        :return: the result of the query, if fetch=True, or an empty list if fetch=False
        """

        assert self._cursor is not None, "You need to connect before executing a query"

        if 'fetch' in kwargs:

            fetch = bool(kwargs['fetch'])

            kwargs.pop('fetch')

        else:

            # Default

            fetch = False

        if 'commit' in kwargs:

            commit = bool(kwargs['commit'])

            kwargs.pop('commit')

        else:

            # Default

            commit = True

        # Execute query

        self._cursor.execute(query, *args, **kwargs)

        # Get the results (if fetch is True)

        if fetch:

            results = self._cursor.fetchall()

        else:

            results = []

        # Commit if commit is true (slow, use False and commit manually if you are doing many queries one after the
        # other)

        if commit:

            self._connection.commit()

        return results

    def remove_table(self, table_name, commit=True):
        """
        Remove a table from the database

        :param table_name:
        :param commit:
        :return:
        """

        # sqlite3 does not allow to use the ? substitution
        # for table names, so we need to do it manually
        table_name = table_name.split(" ")[0]

        self.query("DROP TABLE %s" % table_name, fetch=False, commit=commit)

    def insert_dataframe(self, dataframe, table_name, commit=True):
        """
        Insert a dataframe in a NEW table. It crashes if the table already exists

        :param dataframe:
        :param table_name:
        :param commit: (default: True) whether to commit immediately or not
        :return:
        """

        # Figure out the appropriate dtype conversion
        conv = {}

        for column in dataframe.columns.values:

            conv[column] = array_types[str(dataframe[column].dtype)]

        # Get table definition from Pandas
        table_definition = pd_get_schema(dataframe, table_name, dtype=conv)

        # Create table (but do not commit yet)

        self.query(table_definition, fetch=False, commit=False)

        # Inserts many records at a time
        self.append_dataframe_to_table(dataframe, table_name, commit)

    def append_dataframe_to_table(self, dataframe, table_name, commit=True):

        # How many columns?
        n_columns = dataframe.shape[1]

        # Make a string of ?,?,? with as many ? as columns
        question_marks = ",".join(['?'] * n_columns)

        my_query = '''INSERT INTO %s VALUES (%s)''' % (table_name, question_marks)

        # Insert all data

        self._cursor.executemany(my_query, dataframe.values)

        # Commit if required
        if commit:
            self._connection.commit()

    def get_table_as_dataframe(self, table_name):
        """
        Returns a pandas Dataframe starting from the provided table

        :param table_name:
        :return:
        """

        my_query = 'SELECT * FROM %s' % table_name

        return pd.read_sql_query(my_query, self._connection)
