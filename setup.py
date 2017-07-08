from setuptools import setup, Extension
from Cython.Distutils import build_ext

setup(
    name='lsasaf',
    version='0.1',
    cmdclass={'build_ext': build_ext},
    ext_modules=[Extension('lsasaf._nav',
                           ['lsasaf/_nav.pyx', 'lsasaf/MSG_navigation.c'],
                           libraries=['m'],
                           extra_compile_args=['-O3'])
                 ],
    packages=['lsasaf']
)
