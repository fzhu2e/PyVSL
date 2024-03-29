from setuptools import setup, find_packages

__version__ = '2024.3.16'

with open('README.rst', 'r') as fh:
    long_description = fh.read()

setup(
    name='PyVSL',
    version=__version__,
    description='VS-Lite in Python',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Feng Zhu',
    author_email='fengzhu@ucar.edu',
    url='https://github.com/fzhu2e/PyVSL',
    packages=find_packages(),
    license='GPL-3.0 License',
    keywords='Proxy System Modeling, Tree-ring Width',
    classifiers=[
        'Programming Language :: Python :: 3.9',
        'Programming Language :: Python :: 3.10',
        'Programming Language :: Python :: 3.11',
    ],
    include_package_data=True,
    install_requires=[
        'oct2py',
    ],
)
