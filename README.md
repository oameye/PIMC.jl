# Path Integral Monte Carlo

PIMC.jl is a Julia package that implements the ab-initio continuous path integral Monte Carlo method, as presented by D. Ceperley (1995) and Boninsegni, et al. (2006).
The method provides an efficient way to compute quantum-mechanical properties of a many-body system at finite temperature.

The continuous path integral Monte Carlo method is a statistical technique that uses random sampling to compute properties of a quantum system. In this method, the quantum partition function is approximated by a classical statistical mechanical system in which the particles are represented by continuous paths in imaginary time. The method takes advantage of the fact that at high temperatures, the quantum paths are relatively smooth and can be well approximated by classical paths. 

In PIMC, the continuous paths are sampled using a Monte Carlo algorithm, which generates a set of random paths that are used to estimate the thermodynamic properties of the system. The key idea is to weight each path according to its quantum mechanical probability, so that the statistical average over all paths approximates the true quantum mechanical average. The method has been shown to be highly efficient and accurate, making it a powerful tool for studying many-body quantum systems.


## Installation

To install PIMC, you can use the Julia package manager. Start the Julia REPL and run the following command:

```julia
] add https://github.com/oameye/PIMC.jl
```

## Usage
To use PIMC, simply import the package in your Julia script:
```julia
using PIMC
```

Look in the `example` folder for example how to use the package.

## Contributing
We welcome contributions to PIMC! If you'd like to contribute, please follow these steps:

Fork the repository on Github
Clone your fork to your local machine
Create a new branch for your changes
Make your changes
Commit your changes and push the branch to your fork on Github
Create a pull request to merge your changes into the main repository
We appreciate any contributions, whether it be fixing bugs, adding new features, or improving the documentation.

## License
PIMC is released under the MIT License. See the LICENSE file for more information
