from setuptools import setup, find_packages
from pathlib import Path
this_directory = Path(__file__).parent
long_description = (this_directory / "README.md").read_text()

setup(  name='NEOPRENE',
        packages = find_packages(),
        license = "Apache 2.0",
        version='0.0.10',
        description='🌎 Scripts and information to synthetic generation of precipitation based on Point Processes.',
        long_description=long_description,
        long_description_content_type='text/markdown',
        author='Javier Diez Sierra <javier.diez@unican.es>, Salvador Navas <salvador.navas@unican.es>, Manuel del Jesus <manuel.deljesus@unican.es>',
        author_email='javier.diez@unican.es, salvador.navas@unican.es, manuel.deljesus@unican.es',
        maintainer       = 'Manuel del Jesus',
        maintainer_email = 'manuel.deljesus@unican.es',
        url = 'https://github.com/IHCantabria/NEOPRENE/tree/Crear-ejecutable',

        include_package_data = True,
        python_requires='>=3.7, <4',
        install_requires=[
            'numpy',
            'pandas',
            'scipy',
            'datetime',
            'matplotlib',
            'pyyaml',
        ],
        extras_require={'plotting': ['matplotlib>=2.2.0', 'jupyter','jupyterlab']}
        )