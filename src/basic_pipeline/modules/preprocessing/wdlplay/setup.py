#!/usr/bin/env python3

from setuptools import setup, find_packages
import versioneer

with open('README.md') as f:
    readme = f.read()

dependencies = []
with open('requirements.txt', 'r') as f:
    for line in f:
        dependencies.append(line.strip())

setup(
    name='wdlplay',
    version=versioneer.get_version(),
    cmdclass=versioneer.get_cmdclass(),
    description='Python library for submitting wdl scripts using caper.',
    long_description=readme,
    long_description_content_type='text/markdown',
    author='Tarjinder Singh',
    author_email='ts3475@cumc.columbia.edu`',
    url='https://github.com/tarjindersingh/wdlplay',
    license='MIT license',
    python_requires='>=3.7',
    packages=find_packages(exclude=('tests')),
    install_requires=dependencies
)
