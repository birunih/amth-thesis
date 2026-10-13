# Sequential Monte Carlo for Noncontiguous Attendance Zones

This repository walks through one step of a Sequential Monte Carlo (SMC) algorithm for sampling school attendance zone maps, first on Iowa precinct data and a toy 10x10 grid. It is part of my senior thesis in Applied Mathematics at Yale (advisor: Prof. Shiro Kuriwaki).

**Status:** work in progress. The repository currently demonstrates generating a single zone/district. The full pipeline is under development.

## Motivation

Attendance zone boundaries shape which students attend which schools, and residential segregation often carries over into school segregation. Policy research suggests that noncontiguous zones, which can connect geographically separated neighborhoods, may be one way to promote integration. Standard algorithms for sampling redistricting plans, developed for legislative districts, generally require contiguity. This project adapts an SMC sampler to allow noncontiguous zones so that existing maps can be compared against a distribution of plausible alternatives.

## What's here

- `iowanoncontig.R`: demonstration on Iowa precinct data
- `gridnoncontig.R`: the same step applied to a toy 10x10 grid

## What the code does

1. Represents the map as an adjacency graph of precincts.
2. To generate one zone, we randomly sample a spanning tree, then cut a pair of edges. The resulting district (which may be noncontiguous) is shown in green.
3. The output is an array of possibly generated districts (green).

## Running it

```r
# install.packages(c("foreign", "igraph", "redist", "sf"))
```

## Authorship and credits

Code is mostly my own, written in R. It builds on the `redist` package and on the SMC framework of McCartan and Imai [Annals of Applied Statistics, 2023]. 

## Next steps

- Extending from one step to full sampling
- Evaluating existing maps against the sampled distribution
- Optimization of new maps