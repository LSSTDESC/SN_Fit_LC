from setuptools import setup


setup(
    name='sn_fit_lc',
    version='0.1',
    description='Light curve fitting for supernovae',
    url='http://github.com/lsstdesc/sn_fit_lc',
    author='Philippe Gris',
    author_email='philippe.gris@clermont.in2p3.fr',
    license='BSD',
    packages=['sn_fit', 'sn_fitter'],
    python_requires='>=3.5',
    zip_safe=False,
    install_requires=[
        'sn_tools>=0.1',
        'sn_stackers>=0.1'
    ],
)
