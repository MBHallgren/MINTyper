from setuptools import setup, find_packages
from pathlib import Path

from src.mintyper import __version__

this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

sys.path = [os.path.join(os.path.dirname(os.path.realpath(__file__)),'..')] + sys.path


setup(
    name='mintyper',
    long_description=long_description,
    long_description_content_type='text/markdown',
    version=__version__,
    packages=find_packages(),
    data_files=[],
    include_package_data=True,
    url='https://https://github.com/MBHallgren/mintyper',
    license='',
    install_requires=(),
    author='Malte B. Hallgren',
    scripts=['bin/mintyper'],
    author_email='malhal@food.dtu.dk',
    description='mintyper: an outbreak-detection method for accurate and rapid SNP'
                ' typing of clonal clusters with noisy long reads.'
)