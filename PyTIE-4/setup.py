import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()   

setuptools.setup(
    name = "PyTIE_SMS_DSE",
    version = "0.0.1",
    description = "PyTIE_SMS_DSE is an open-source Python package. It is designed for obtaining molecular descriptor expressions in constant time. The computation time of PyTIE is very efficiently faster than manual computation",
    license = 'MIT',
    author = 'Sahaya Vijay J',
    author_email = 'sahayavijay.j@gmail.com',
    classifiers = [
        'Development Status :: 3 - Alpha',
        'Environment :: Console',
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
        'Programming Language :: Python'
    ],
    install_requires = ['numpy','math'],
)
