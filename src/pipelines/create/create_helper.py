from dotenv import load_dotenv
import os
import sys
# src/pipelines/create/create_helper.py -> src
src_dir = os.path.dirname(os.path.dirname(os.path.dirname(__file__)))
if src_dir not in sys.path:
    sys.path.append(src_dir)


project_dir = os.path.dirname(src_dir)
env_path = os.path.join(project_dir, '.env')
load_dotenv(env_path)