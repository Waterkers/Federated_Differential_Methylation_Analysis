#from distutils.core import setup, Extension
#from setuptools import setup, Extension
#from Cython.Build import cythonize
#import numpy

#linbin = Extension('linbinR', ['linbinR.c'], include_dirs=["C:\\Users\\Silke\\AppData\\Roaming\\Python\\Python37\\site-packages\\numpy\\core\\include\\numpy",
#"C:\\Users\\Silke\\AppData\\Roaming\\Python\\Python37\\site-packages\\numpy\\core\\include\\numpy\\arrayobject.h"])
from setuptools import setup, Extension
from Cython.Build import cythonize
import numpy as np

# setup(
    # name='linbinR',
    # ext_modules=
        # Extension(
            # "linbinR",
            # sources=["linbinR.c"]
        # )
    # )

# old version that worked before handing in thesis 
setup(
    name='linbinR',
    ext_modules=cythonize(
        Extension(
            "linbinR",
            sources=["linbinR.pyx"],
            include_dirs=[np.get_include()]
        )
    ),
    install_requires=["numpy"]
    )