from dotenv import load_dotenv
import os
import sys
src_dir = os.path.dirname(os.path.dirname(__file__))
if src_dir not in sys.path:
    sys.path.append(src_dir)


from bioomics import BuildComplex

env_path = os.path.join(src_dir, 'pipelines', '.env')
load_dotenv(env_path)