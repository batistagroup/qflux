
import sys
import os

ROOT_DIR = os.path.abspath('..')
print(ROOT_DIR)
if ROOT_DIR not in sys.path:
    sys.path.append(ROOT_DIR)
