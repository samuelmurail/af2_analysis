from setuptools import setup, find_packages

version = "0.0.2"

with open('README.md', encoding='utf-8') as readme_file:
    readme = readme_file.read()

requirements = [
    'pandas>=1.3.4',
    'numpy>=1.21',
    'tqdm>=4.0',
    'seaborn>=0.11',
    'pdb_numpy>=0.0.6',
    'cmcrameri>=1.7',
    'nglview>=3.0',
    'ipywidgets>=7.6,
]

setup(
    name='af2_analysis',
    version=version,
    description=(
        '`AF2 analysis` is a python library allowing analysis of Alphafold results.'
    ),
    long_description=readme,
    long_description_content_type='text/markdown',
    author='Samuel Murail',
    author_email="samuel.murail@u-paris.fr",
    url='https://github.com/samuelmurail/af2_analysis',
    packages=find_packages(where='src'),
    package_dir={'af2_analysis': 'src/af2_analysis', 'af2_analysis.format': 'src/af2_analysis/format'},
    entry_points={'console_scripts': ['af2_analysis = af2_analysis.__main__:main']},
    include_package_data=True,
    python_requires='>=3.7',
    install_requires=requirements,
    license='GNUv2.0',
    zip_safe=False,
    classifiers=[
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Developers",
        "Natural Language :: English",
        "License :: OSI Approved :: BSD License",
        "Programming Language :: Python :: 3 :: Only",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: Implementation :: CPython",
        "Programming Language :: Python :: Implementation :: PyPy",
        "Programming Language :: Python",
        "Topic :: Software Development",
    ],
    keywords=[
        "af2_analysis",
        "Python",
        "AlphaFold2",
        "ColabFold",
    ],
)
