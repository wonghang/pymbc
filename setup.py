from setuptools import setup, Extension
import numpy as np
    
pymbc_module = Extension(
    'pymbc',
    sources = ['pymbc.c', 'bitset.c', 'bipartite_graph.c', 'lyu.c'],
    include_dirs=[np.get_include()],
    define_macros=[('_GNU_SOURCE',1)],
    #undef_macros=['NDEBUG'],
    extra_compile_args = ['-std=gnu99']
)

long_description = """
###########################
Maximum Biclique Search
###########################

===================

Introduction
------------
| An implementation of
| *Lyu, Bingqing, et al. "Maximum biclique search at billion scale." Proceedings of the VLDB Endowment (2020).*
|
"""

setup(
    name='pymbc',
    version='0.0.1',
    description='Maximum Biclique Search',
    long_description=long_description,
    url='https://github.com/wonghang/pymbc',
    author='wonghang',
    author_email="wonghang@gmail.com",
    setup_requires=['numpy>=1.17.0'],
    install_requires=['numpy>=1.17.0'],
    classifiers=[
        'Development Status :: 3 - Alpha',
        'Intended Audience :: Developers',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: BSD License',
        'Programming Language :: C',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.3',
        'Programming Language :: Python :: 3.4',
        'Programming Language :: Python :: 3.5',
        'Programming Language :: Python :: 3.6',
        'Operating System :: POSIX :: Linux',
        'Topic :: Software Development :: Libraries',
        'Topic :: Scientific/Engineering :: Bio-Informatics',
        'Topic :: Scientific/Engineering :: Artificial Intelligence',
        'Topic :: Scientific/Engineering :: Mathematics',
    ],
    keywords='graph theory,biclique,maximum biclique search',
    platforms = ["Linux"],
    include_package_data=True,
    ext_modules = [pymbc_module],
)
