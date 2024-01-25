import setuptools

setuptools.setup(
    name = "PyTIE_D",
    version = "0.0.1",
    description = "PyTIE_D is an open-source Python package. It is designed for obtaining molecular descriptor expressions and its numerical values in constant time. The computation time of PyTIE is very efficiently faster than manual computation",
    license = 'MIT',
    author = 'Sahaya Vijay J',
    author_email = 'sahayavijay.j@gmail.com',
    packages=setuptools.find_packages(),
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python'
    ],
    install_requires = ['numpy','math','sympy'],
)
