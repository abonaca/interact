from distutils.core import setup, Extension
import numpy as np

c_ext = Extension("interact3", ["interact_wrap.c", "interact.c", "nr_rand.c", "utility.c"])

setup(ext_modules = [c_ext], 
      include_dirs = np.get_include())
