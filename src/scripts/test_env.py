#!/usr/bin/python
# python scripts/test_env.py
'''
Environment sanity check:
    1. environmental variables loaded from .env
    2. dependencies declared in pyproject.toml
    3. CUDA and GPU-related packages (use CUDA when available)
    4. database connections (MySQL, Redis)
'''
import importlib
import importlib.metadata
import os
import socket
import sys
from pathlib import Path

# load env
# src/.env takes precedence, the repo root .env is the fallback
# (same convention as src/bioomics/database/db.py)
from dotenv import load_dotenv

SCRIPT_DIR = Path(__file__).resolve()
load_dotenv(SCRIPT_DIR.parents[1] / '.env')
load_dotenv(SCRIPT_DIR.parents[2] / '.env')

# section separators, keep the original style
SEPARATOR = '#############'

# accumulated results: (section, name, status, detail)
RESULTS = []
PASS, FAIL, SKIP = 'PASS', 'FAIL', 'SKIP'


def record(section: str, name: str, status: str, detail: str = ''):
    RESULTS.append((section, name, status, detail))
    print(f'[{status:>4}] {name}' + (f' - {detail}' if detail else ''))


def summary() -> int:
    '''print the summary and return the exit code (0 ok, 1 any FAIL)'''
    failed = [r for r in RESULTS if r[2] == FAIL]
    skipped = [r for r in RESULTS if r[2] == SKIP]
    passed = [r for r in RESULTS if r[2] == PASS]
    print(f'\n{SEPARATOR} summary {SEPARATOR}')
    print(f'pass: {len(passed)}  fail: {len(failed)}  skip: {len(skipped)}')
    for section, name, status, detail in failed:
        print(f'  FAIL {section} :: {name}' + (f' - {detail}' if detail else ''))
    return 1 if failed else 0


def section(title: str):
    print(f'\n\n{SEPARATOR} {title} {SEPARATOR}')


#############
# 1. test environmental variables
#############
def test_env_variables():
    section('1. environmental variables')
    # key: (required, path)
    env_keys = {
        # server
        'DATA_HOST_IP': (False, False),
        'DATA_HOSTNAME': (False, False),
        'DATA_HOST_USER': (False, False),
        'DATA_DIR': (False, True),
        # input data
        'rawdata_dir': (True, True),
        # output data
        'output_dir': (True, True),
        'alphafold_dir': (False, True),
        'tmp_dir': (False, True),
        'colabfold_bin_dir': (False, True),
        'colabfold_db': (False, True),
        # external tools
        'tools_dir': (True, True),
        'igblast_dir': (False, True),
        'muscle_bin': (False, True),
        'mmalign_dir': (False, True),
        # mysql database
        'MYSQL_HOST': (True, False),
        'MYSQL_PORT': (True, False),
        'MYSQL_DATABASE': (True, False),
        'MYSQL_USER': (True, False),
        'MYSQL_PASSWORD': (True, False),
        'MYSQL_DATA_DIR': (False, True),
        # other
        'verbose': (False, False),
    }
    for key, (required, is_path) in env_keys.items():
        value = os.getenv(key)
        if value is None or value == '':
            status = FAIL if required else SKIP
            detail = 'required' if required else 'optional, not set'
            record('env', key, status, detail)
            continue
        # strip quotes often left in .env files
        value = value.strip().strip('"\'')
        if is_path:
            exists = Path(value).exists()
            detail = value if exists else f'{value} (not found)'
            record('env', key, PASS if exists else FAIL, detail)
        else:
            # do not print passwords
            shown = '***' if 'PASSWORD' in key else value
            record('env', key, PASS, shown)


#############
# 2. test dependencies
#############
# pyproject.toml dependencies -> import names (None = meta package, no import)
DEPENDENCIES = {
    'bio': 'Bio',
    'bs4': 'bs4',
    'hmmer': None,
    'jupyter': None,
    'lxml': 'lxml',
    'matplotlib': 'matplotlib',
    'matplotlib-venn': 'matplotlib_venn',
    'mysql-connector-python': 'mysql.connector',
    'netcdf4': 'netCDF4',
    'notebook': None,
    'pandas': 'pandas',
    'plotly': 'plotly',
    'prodigy-prot': None,
    'py3dmol': 'py3Dmol',
    'python-dotenv': 'dotenv',
    'redis': 'redis',
    'requests': 'requests',
    'scikit-learn': 'sklearn',
    'seaborn': 'seaborn',
    'sh': 'sh',
    'torch': 'torch',
    'torchaudio': 'torchaudio',
    'torchvision': 'torchvision',
    'urllib3': 'urllib3',
    'venny4py': 'venny4py',
}
DEV_DEPENDENCIES = {
    'ipython': 'IPython',
    'mypy': 'mypy',
    'pytest': 'pytest',
    'ruff': None,
    'ddt': 'ddt',
}


def test_dependencies(deps: dict, group: str = 'deps'):
    section(f'2. dependencies ({group})')
    for dist_name, import_name in deps.items():
        try:
            version = importlib.metadata.version(dist_name)
        except importlib.metadata.PackageNotFoundError:
            record(group, dist_name, FAIL, 'not installed')
            continue
        if import_name is None:
            record(group, dist_name, PASS, f'{version} (meta package)')
            continue
        try:
            module = importlib.import_module(import_name)
            module_version = getattr(module, '__version__', version)
            record(group, dist_name, PASS, f'{module_version}')
        except Exception as e:
            record(group, dist_name, FAIL, f'{version} installed but import failed: {e}')


#############
# 3. test cuda and related packages
#############
def test_torch():
    import torch
    record('cuda', 'torch', PASS, torch.__version__)
    record('cuda', 'torch.version.cuda', PASS, torch.version.cuda or 'none (cpu build)')

    if not torch.cuda.is_available():
        record('cuda', 'torch.cuda', SKIP, 'CUDA not available')
        return

    count = torch.cuda.device_count()
    record('cuda', 'torch.cuda.device_count', PASS, str(count))
    for i in range(count):
        props = torch.cuda.get_device_properties(i)
        memory = props.total_memory / 1024**3
        record('cuda', f'torch.cuda[{i}]', PASS,
               f'{props.name}, {memory:.1f} GB, capability {props.major}.{props.minor}')

    # use CUDA if there is
    try:
        device = torch.device('cuda')
        a = torch.randn(1024, 1024, device=device)
        b = torch.randn(1024, 1024, device=device)
        c = a @ b
        torch.cuda.synchronize()
        record('cuda', 'torch matmul on GPU', PASS, f'{c.shape}, device={c.device}')
    except Exception as e:
        record('cuda', 'torch matmul on GPU', FAIL, str(e))


def test_jax():
    # avoid preallocating 75% of GPU memory (conflicts with torch and other jobs)
    os.environ.setdefault('XLA_PYTHON_CLIENT_PREALLOCATE', 'false')
    try:
        import jax
    except Exception as e:
        record('cuda', 'jax', FAIL, f'import failed: {e}')
        return
    record('cuda', 'jax', PASS, jax.__version__)
    try:
        devices = jax.devices()
        record('cuda', 'jax.devices', PASS,
               ', '.join(f'{d.device_kind} [{d.platform}]' for d in devices))
        if any(d.platform == 'gpu' for d in devices):
            # use CUDA if there is
            x = jax.numpy.ones((1024, 1024))
            y = (x @ x).block_until_ready()
            record('cuda', 'jax matmul on GPU', PASS, f'{y.shape}')
        else:
            record('cuda', 'jax matmul on GPU', SKIP, 'no gpu device')
    except Exception as e:
        record('cuda', 'jax runtime', FAIL, str(e))


def test_cuda_packages():
    section('3. cuda and related packages')
    test_torch()

    # torchvision / torchaudio must match the installed torch
    try:
        import torchvision
        record('cuda', 'torchvision', PASS, torchvision.__version__)
    except Exception as e:
        record('cuda', 'torchvision', FAIL, str(e))
    try:
        import torchaudio
        record('cuda', 'torchaudio', PASS, torchaudio.__version__)
    except Exception as e:
        record('cuda', 'torchaudio', FAIL, str(e))

    test_jax()

    # PyRosetta is optional (not in pyproject.toml)
    try:
        import pyrosetta
        pyrosetta.init('-mute all')
        record('cuda', 'pyrosetta', PASS, pyrosetta.version())
    except ImportError:
        record('cuda', 'pyrosetta', SKIP, 'not installed (optional)')
    except Exception as e:
        record('cuda', 'pyrosetta', FAIL, str(e))


#############
# 4. test database connections
#############
def test_mysql():
    import mysql.connector
    # same configuration as src/bioomics/database/db.py
    config = {
        'host': os.getenv('MYSQL_HOST'),
        'user': os.getenv('MYSQL_USER'),
        'password': os.getenv('MYSQL_PASSWORD'),
        'port': int(os.getenv('MYSQL_PORT', '3306')),
        'use_pure': True,
        'buffered': True,
        'autocommit': True,
        'connection_timeout': 10,
        'pool_size': 1,
    }
    database = os.getenv('MYSQL_DB') or os.getenv('MYSQL_DATABASE')
    if database:
        config['database'] = database

    connection = None
    cursor = None
    try:
        connection = mysql.connector.connect(**config)
        cursor = connection.cursor(dictionary=True)
        cursor.execute('SELECT 1')
        cursor.fetchall()
        record('db', 'mysql connect', PASS, f"{config['host']}:{config['port']}")
        cursor.execute('SELECT VERSION() AS version')
        version = cursor.fetchone()['version']
        record('db', 'mysql server', PASS, f'server {version}')
        record('db', 'mysql database', PASS, database if database else 'not selected')
    except Exception as e:
        record('db', 'mysql connect', FAIL, f"{config['host']}:{config['port']} {e}")
    finally:
        if cursor is not None:
            cursor.close()
        if connection is not None and connection.is_connected():
            connection.close()


def test_redis():
    import redis
    host = 'localhost'
    port = 6379
    # redis is an optional local service, same as tests/connector/test_conn_redis.py
    try:
        with socket.create_connection((host, port), timeout=2):
            pass
    except OSError:
        record('db', 'redis', SKIP, f'not running on {host}:{port}')
        return
    try:
        conn = redis.StrictRedis(host=host, port=port, db=0, decode_responses=True)
        assert conn.ping()
        record('db', 'redis', PASS, f'{host}:{port} ping ok')
    except Exception as e:
        record('db', 'redis', FAIL, f'{host}:{port} {e}')


def test_database_connections():
    section('4. database connections')
    test_mysql()
    test_redis()


def main():
    test_env_variables()
    test_dependencies(DEPENDENCIES)
    test_dependencies(DEV_DEPENDENCIES, group='deps-dev')
    test_cuda_packages()
    test_database_connections()
    sys.exit(summary())


if __name__ == '__main__':
    main()


