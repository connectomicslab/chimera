#!/usr/bin/env python

"""The setup script."""

from setuptools import setup, find_packages

with open('README.rst') as readme_file:
    readme = readme_file.read()

with open('HISTORY.rst') as history_file:
    history = history_file.read()

requirements = [ ]

test_requirements = [ ]

setup(
    author="Yasser Aleman Gomez",
    author_email='yasseraleman@gmail.com',
    python_requires='>=3.6',
    classifiers=[
        'Development Status :: 2 - Pre-Alpha',
        'Intended Audience :: Developers',
        'License :: OSI Approved :: Apache Software License',
        'Natural Language :: English',
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.6',
        'Programming Language :: Python :: 3.7',
        'Programming Language :: Python :: 3.8',
    ],
    description="Atlases and parcellation methods are key to unravel the anatomo-functional organization of the human brain. They subdivide brain's anatomy allowing us to perform region-based analyses and narrowing the hypotheses",
    entry_points={
        'console_scripts': [
            'chimera=chimera.chimera:main',
        ],
    },
    install_requires=requirements,
    license="Apache Software License 2.0",
    long_description=readme + '\n\n' + history,
    include_package_data=True,
    keywords='chimera',
    name='chimera',
    packages=find_packages(include=['chimera', 'chimera.*']),
    test_suite='tests',
    tests_require=test_requirements,
    url='https://github.com/yasseraleman/chimera',
    version='0.1.0',
    zip_safe=False,
)
