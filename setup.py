from setuptools import setup

# get the version here
pkg_vars = {}

with open("version.py") as fp:
    exec(fp.read(), pkg_vars)

setup(
    name='sn_fit_lc',
    version=pkg_vars['__version__'],
    description='Light curve fitting for supernovae',
    url='http://github.com/lsstdesc/sn_fit_lc',
    author='Philippe Gris',
    author_email='philippe.gris@clermont.in2p3.fr',
    license='BSD',
    packages=['sn_fit', 'sn_fitter', 'sn_fit_input'],
    # All files from folder sn_fit_input
    package_data={'sn_fit_input': ['*.txt']},
    python_requires='>=3.5',
    zip_safe=False,
    install_requires=[
        'sn_tools>=0.1',
        'iminuit>=1.4'
    ],
)
