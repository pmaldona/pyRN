from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

setup(
    name='pyRN',

    version='0.6',

    description='Reaction Network Structure and Simulation Library for Resilience',
    long_description='Reaction Network Structure and Simulation Library for Resilience',

    url='https://github.com/pmaldona/pyRN',

    author='Alejandro Bassi, Fionn Daire Keogh, Pedro Maldonado, Tomas Veloz',
    author_email='pmaldona@sax.cl',

    license='GNU v.3',

    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science',
        'Intended Audience :: Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3) ',
        'Programming Language :: Python :: 3.9',
    ],

    keywords='chemical reaction networks, closed structure, decomposition, simulation, random walk, random network generator',

    packages=find_packages(exclude=['tests']),

    # run-time dependencies that will be installed by pip
    # install_requires=['numpy','pandas','beautifulsoup4', 'bitarray','scipy','networkx','libroadrunner','lxml','pyvis','pypoman','matplotlib', 'joblib']
    install_requires=['numpy','pandas','beautifulsoup4', 'bitarray','scipy','networkx','libroadrunner','lxml','pyvis','matplotlib', 'joblib']
)
