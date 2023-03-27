from setuptools import setup

# read the contents of your README file
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(
    name='thunder_ase',
    version='0.2.2',
    packages=['thunder_ase'],
    url='https://github.com/thunder-dft/thunder-ase',
    license='LGPL',
    author='renpj',
    author_email='openrpj@gmail.com',
    description='ASE calculator interface for FIREBALL code.'
    long_description=long_description,
    long_description_content_type='text/markdown'
)
