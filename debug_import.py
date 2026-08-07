import os
import sys

print("Executable:", sys.executable)
print("Working directory:", os.getcwd())

print("\nsys.path:")
for p in sys.path:
    print(repr(p))

import qflux

print("\nImported qflux from:")
print(qflux.__file__)