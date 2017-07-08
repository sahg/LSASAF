from setuptools import setup, Extension
from Cython.Distutils import build_ext

setup(
    name='pymsg',
    version='0.1',
    cmdclass={'build_ext': build_ext},
    ext_modules=[Extension('pymsg._nav',
                           ['pymsg/_nav.pyx', 'pymsg/MSG_navigation.c'],
                           libraries=['m'],
                           extra_compile_args=['-O3'])
                 ],
    packages=['pymsg']
)
