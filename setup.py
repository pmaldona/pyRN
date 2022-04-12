
from setuptools import setup, find_packages
# To use a consistent encoding
from codecs import open
from os import path

setup(
    name='CRNS',

    version='0.6',

    description='Closed Reaction Network Structure Library',
    long_description='Closed Reaction Network Structure Library',

    url='https://github.com/pmaldona/CRNS',

    author='Alejandro Bassi, Pedro Maldonado',
    author_email='pmaldona@sax.cl',

    license='GNU v.3',

    classifiers=[
        'Development Status :: 4 - Beta',
        'Intended Audience :: Science',
        'Intended Audience :: Research',
        'License :: OSI Approved :: GNU General Public License v3 (GPLv3) ',
        'Programming Language :: Python :: 3.9',
    ],

    keywords='chemical reaction networks, closed structure',

    packages=find_packages(exclude=['tests']),

    # run-time dependencies that will be installed by pip
    install_requires=['numpy','pandas','BeautifulSoup4', 'bitarray','scipy','networkx']
)
