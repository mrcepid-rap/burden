from setuptools import setup, find_packages
from pathlib import Path


def load_requirements(fname: Path):
    reqs = []
    with fname.open('r') as reqs_file:
        for line in reqs_file:
            reqs.append(line.rstrip())
    return reqs


setup(
    name='burden',
    version='1.1.3',
    packages=find_packages(),
    url='https://github.com/mrcepid-rap/mrcepid-runassociationtesting-burden',
    license='',
    author='Eugene Gardner',
    author_email='eugene.gardner@mrc-epid.cam.ac.uk',
    description='',
    install_requires=load_requirements(Path('requirements.txt'))
)
