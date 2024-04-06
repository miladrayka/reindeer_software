from setuptools import setup, find_packages

setup(
    name='REINDEER',
    version='0.1.0',
    description="""REINDEER is a python package to generate
    feature vectors for a protein-ligand complex.""",
    url='https://github.com/miladrayka/reindeer_software',
    author='Milad Rayka',
    author_email='milad.rayka@yahoo.com',
    license='MIT License',
    #packages=['reindeer'],
    package_dir={"reindeer": "reindeer"},
    packages=find_packages(where="reindeer"),
    python_requires=">=3.8",
    install_requires=['pandas==2.2.1',
                      'numpy==1.26.4',
                      'biopython==1.83',
                      'rdkit-pypi',
                      'scipy==1.13.0',
                      'joblib==1.3.2',
                      'typing-extensions==4.11.0',
                      'tqdm==4.66.2',
                      ],

    classifiers=[
        'Development Status :: Under Development',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: POSIX :: Windows',
        'Programming Language :: Python :: 3.9',
    ],
    )
