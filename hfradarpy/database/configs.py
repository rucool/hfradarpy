import datetime as dt
import pandas as pd
from collections import OrderedDict


def mysql_configs():
    database = {
        'drivername': 'mysql+mysqlconnector',
        'username': 'coolops',
        'password': 'SjR9rM',
        'host': 'mysql1.marine.rutgers.edu',
        'database': 'coolops'}
    return database


def mongodb_configs():
    database = {
        'uri': 'mongodb://michaesm:Mia8ge6soGcm@mongodb.marine.rutgers.edu:27017'}
    return database
