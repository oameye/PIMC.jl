import pathlib
from setuptools import setup

# The directory containing this file
HERE = pathlib.Path(__file__).parent

# The text of the README file
README = (HERE / "README.md").read_text()

# This call to setup() does all the work
setup(
    name="erroranalysis-py",
    version="0.6.11",
    description="Generate error bars and perform binning analysis using jackknife or bootstrap resampling. Calculate average and error in quantum Monte Carlo data (or other data) and on functions of averages (such as fluctuations, skew, and kurtosis).",
    long_description=README,
    long_description_content_type="text/markdown",
    url="https://github.com/nscottnichols/erroranalysis-py",
    author="Nathan Nichols",
    author_email="Nathan.Nichols@uvm.edu",
    license="MIT",
    classifiers=[
        "License :: OSI Approved :: MIT License",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
    ],
    packages=["erroranalysis"],
    include_package_data=True,
    install_requires=["numpy", "matplotlib"],
    entry_points={
        "console_scripts": [
            "ea=erroranalysis.__main__:main",
        ]
    },
)
