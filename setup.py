from setuptools import setup, find_packages

__version__ = '0.0.8'

with open('README.rst', 'r') as fh:
    long_description = fh.read()

setup(
    name='PyVSL',
    version=__version__,
    description='VS-Lite in Python',
    long_description=long_description,
    long_description_content_type='text/markdown',
    author='Feng Zhu',
    author_email='fzhu@nuist.edu.cn',
    url='https://github.com/fzhu2e/PyVSL',
    packages=find_packages(),
    license='GPL-3.0 License',
    keywords='Proxy System Modeling, Tree-ring Width',
    classifiers=[
        'Programming Language :: Python :: 3',
        'Programming Language :: Python :: 3.9',
    ],
    install_requires=[
        'oct2py',
    ],
)
