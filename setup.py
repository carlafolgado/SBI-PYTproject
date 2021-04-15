from setuptools import setup, find_packages

with open("README.md", "r") as fh:
    long_description = fh.read()

setup(name='biocomplexbuilder',
    version='0.1.0',
    description='Generate macrocomplexes superimposing and aligning paired interacting elements.',
    #data_files = [("", ['LICENSE'])],
    license = "MIT",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url='https://github.com/carlafolgado/SBI-PYTproject.git',
    author='Oriol Canal Pujol, Carla Folgado Salido and Arnau Llinàs Bertran',
    author_email='oriol.canal01@estudiant.upf.edu, carla.folgado01@estudiant.upf.edu, arnau.llinas01@estudiant.upf.edu',
    packages=['biocomplexbuilder'],
    py_modules = ['biocomplexbuilder/arguments','biocomplexbuilder/builder','biocomplexbuilder/DNAbased_utilities','biocomplexbuilder/utilities' ],
    classifiers=[
    "Programming Language :: Python :: 3",
    "License :: OSI Approved :: MIT License"
    ],
    install_requires=['biopython', 'modeller'],
    include_package_data=True)
