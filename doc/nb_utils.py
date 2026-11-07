'''
shared utilities for the analysis notebooks in doc/

- mysql connection configured from the repo .env (mirrors src/bioomics/database/db.py)
- query helpers: query_df / query_scalar
- export helpers: save_fig -> doc/figures/*.png, save_csv -> doc/stats/*.csv

outputs:
    - doc/figures/*.png
    - doc/stats/*.csv
'''
from pathlib import Path
import os

import mysql.connector
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from dotenv import load_dotenv
from scipy import stats


def _find_repo_root(start: Path = None) -> Path:
    '''walk upwards until the repo root (contains .env) is found'''
    start = start or Path(__file__).resolve().parent
    for p in (start, *start.parents):
        if (p / '.env').is_file():
            return p
    return start


REPO_DIR = _find_repo_root()
DOC_DIR = REPO_DIR / 'doc'
FIG_DIR = DOC_DIR / 'figures'
STAT_DIR = DOC_DIR / 'stats'
FIG_DIR.mkdir(parents=True, exist_ok=True)
STAT_DIR.mkdir(parents=True, exist_ok=True)

# database configuration, same keys as src/bioomics/database/db.py
load_dotenv(DOC_DIR / '.env')
load_dotenv(REPO_DIR / '.env')

DB_CONFIG = {
    'host': os.getenv('MYSQL_HOST'),
    'user': os.getenv('MYSQL_USER'),
    'password': os.getenv('MYSQL_PASSWORD'),
    'port': int(os.getenv('MYSQL_PORT', '3306')),
    'database': os.getenv('MYSQL_DB') or os.getenv('MYSQL_DATABASE'),
    'use_pure': True,
    'buffered': True,
    'autocommit': True,
    'connection_timeout': 60,
}


def query_df(query: str) -> pd.DataFrame:
    '''run a SELECT and return the result as a DataFrame (dictionary cursor)'''
    connection = mysql.connector.connect(**DB_CONFIG)
    try:
        cursor = connection.cursor(dictionary=True)
        cursor.execute(query)
        rows = cursor.fetchall()
        cursor.close()
        return pd.DataFrame(rows)
    finally:
        connection.close()


def query_scalar(query: str):
    '''run a SELECT returning a single value'''
    df = query_df(query)
    if len(df) == 0:
        return None
    return df.iloc[0, 0]


def save_fig(fig, name: str, dpi: int = 200) -> Path:
    '''export a matplotlib figure to doc/figures/<name>.png'''
    path = FIG_DIR / f'{name}.png'
    fig.savefig(path, dpi=dpi, bbox_inches='tight')
    print(f'saved figure: {path}')
    return path


def save_csv(df: pd.DataFrame, name: str, index: bool = False) -> Path:
    '''export a statistical table to doc/stats/<name>.csv'''
    path = STAT_DIR / f'{name}.csv'
    df.to_csv(path, index=index)
    print(f'saved table: {path}  shape={df.shape}')
    return path


def describe_groups(df: pd.DataFrame, group_col: str, value_col: str) -> pd.DataFrame:
    '''descriptive statistics per group, ready for export'''
    summ = df.groupby(group_col)[value_col].agg(
        n='count', mean='mean', std='std', min='min',
        q25=lambda s: s.quantile(0.25),
        median='median',
        q75=lambda s: s.quantile(0.75),
        max='max',
    )
    return summ.reset_index()


def pairwise_mannwhitney(df: pd.DataFrame, group_col: str, value_col: str,
                         min_n: int = 500) -> pd.DataFrame:
    '''pairwise Mann-Whitney U tests between groups, Bonferroni corrected'''
    groups = {}
    for name, g in df.groupby(group_col):
        vals = g[value_col].dropna().to_numpy()
        if len(vals) >= min_n:
            groups[name] = vals
    names = list(groups)
    rows = []
    n_tests = len(names) * (len(names) - 1) // 2
    for i in range(len(names)):
        for j in range(i + 1, len(names)):
            u, p = stats.mannwhitneyu(groups[names[i]], groups[names[j]],
                                      alternative='two-sided')
            rows.append({
                'group_a': names[i], 'group_b': names[j],
                'n_a': len(groups[names[i]]), 'n_b': len(groups[names[j]]),
                'median_a': float(np.median(groups[names[i]])),
                'median_b': float(np.median(groups[names[j]])),
                'U': u, 'p_value': p,
                'p_bonferroni': min(p * n_tests, 1.0),
            })
    return pd.DataFrame(rows)


# plotting style
sns.set_theme(style='whitegrid', context='notebook')
PALETTE = sns.color_palette('colorblind')
plt.rcParams['figure.facecolor'] = 'white'
