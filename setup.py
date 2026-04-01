from setuptools import setup, find_packages
import os

requires = ['numpy',
            'requests==2.31.0',
            'pandas==2.1.3',
            'numpy==1.26.0',
            'seaborn==0.13.2',
            'matplotlib==3.8.0',
            'tqdm==4.65.0',
            'scanpy==1.9.8',
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