from distutils.core import setup
import glob
import os

scripts = glob.glob(os.path.join('scripts', '*.py'))

# noinspection PyPackageRequirements
setup(

    name='LSST-Transients',

    version='0.1',

    packages=['lsst_transients', 'lsst_transients.utils'],

    url='https://github.com/jackieblaum/LSST-Transients',

    license='BSD-3',

    author='Jackie Blaum, Giacomo Vianello',

    author_email='jackie.blaum@gmail.com',

    description='Finding transients in LSST data with a new technique',

    install_requires=['numpy',
                      'scipy',
                      'astropy',
                      'pandas',
                      'photutils',
                      'imageio',
                      'pyregion',
                      'matplotlib'],

    scripts=scripts,

)
