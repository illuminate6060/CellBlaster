from setuptools import setup, find_packages
import os

requires = ['numpy',
            'requests',
            'pandas',
            'numpy',
            'seaborn',
            'matplotlib',
            'tqdm',
            'scanpy',
            ]

setup(
      name='CellBlaster', 
      version='1.0', 
      author='Lin Du',
      author_email='3051065449@qq.com',
      license="MIT",
      description="CellBLASTer: A universal plant scRNA-seq annotation tool inspired by cellular BLAST strategies",
      long_description=open('README.md').read(),
      packages=find_packages(), 
      install_requires = requires,
      python_requires = '>=3.7'
)