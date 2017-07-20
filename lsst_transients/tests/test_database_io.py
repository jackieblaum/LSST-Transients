# You only need this if this is a test (do not import this in other scripts)
import pytest

from lsst_transients.utils import database_io
import pandas as pd
import os
import numpy as np


def test_database_io():

    database_name = "test.db"

    # Remove database if it exists

    try:

        os.remove(database_name)

    except OSError:

        pass

    # Fake data to insert into the database

    df = pd.DataFrame()

    df['pippo'] = np.random.uniform(0, 1, 1)
    df['pluto'] = np.random.uniform(0, 1, 1)

    db = database_io.SqliteDatabase(database_name)

    db.connect()

    db.insert_dataframe(df, 'test')

    df2 = db.get_table_as_dataframe("test")

    assert np.allclose(df.values, df2.values)

    # Insert many tables

    with database_io.bulk_operation(db):

        for i in range(1000):

            if i % 100 == 0:

                print("%i" % i)

            db.insert_dataframe(df, 'test%i' % i, commit=False)

    db.disconnect()

    # Remove the database so we do not leave things around after the test is done

    try:

        os.remove(database_name)

    except OSError:

        pass