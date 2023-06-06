from distutils.core import setup, Extension
import numpy as np

c_ext = Extension("interact3", sources=["interact_wrap.c", "interact.c", "nr_rand.c", "utility.c"], libraries=['hdf5', 'hdf5_hl'], library_dirs=['/usr/lib64'])

setup(ext_modules = [c_ext],
      include_dirs = np.get_include())
