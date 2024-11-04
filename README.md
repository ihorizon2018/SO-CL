# SO-CL

> ### A hybrid Snake Optimizer with Crisscross Learning (SO-CL) strategy for constrained structural optimization
> #### Note: All source code will be uploaded after the paper has been accepted


## Overview

SO-CL (hybrid Snake Optimizer with Crisscross Learning) is an advanced optimization algorithm designed specifically for constrained structural optimization tasks. This repository implements the SO-CL algorithm and demonstrates its effectiveness on several classic benchmark functions and complex engineering design problems.

## Benchmarks

The repository includes implementations for 23 classic benchmark functions, including all test instances of the CEC benchmark. These benchmarks provide a robust testing ground for the optimizer and highlight its performance across various constrained optimization scenarios.

### Running the Benchmark Tests

To test the algorithm on the benchmark functions, you can run the `algorithm.m` file. This will execute SO-CL across the included test instances.

## Engineering Design Problems

SO-CL is also tested on five well-known engineering design problems in the field of structural optimization. These cases involve truss structure optimization and demonstrate the algorithm's practical applicability in engineering.

### Supported Engineering Problems:

1. **200-Truss Static Optimization**
2. **10-Truss Static Optimization**
3. **25-Truss Static Optimization**
4. **72-Truss Static Optimization**
5. **120-Truss Static Optimization**

The truss analyzer function, originally adapted for these cases, has been modified for compatibility with SO-CL. 

### Running the Engineering Design Problems

To execute these design problems, run the `main.m` file. This will allow you to explore the performance of SO-CL on each of the listed engineering challenges.

## Methodology
More details can be found in the 'Paper.pdf'
[A hybrid snake optimizer with crisscross learning strategy for constrained structural optimization]

## License
GNU General Public License v3.0

---

