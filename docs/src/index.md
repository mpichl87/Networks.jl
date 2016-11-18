# Networks.jl Documentation


```@contents
```

## Status

![Travis Status](https://travis-ci.org/mpichl87/Networks.jl.svg)
[![Coverage Status](https://coveralls.io/repos/github/mpichl87/Networks.jl/badge.svg?branch=master)](https://coveralls.io/github/mpichl87/Networks.jl?branch=master)

## Description

One method for S-parameter simulation of linear networks is to connect the ports of smaller S-matrices.
[This method](http://qucs.sourceforge.net/docs/technical/technical.pdf) is used in [qucs](http://qucs.sourceforge.net/).

There are two cases to consider:

1. Connect two ports of the same network

    If two ports of the same circuit $S$ are connected, the new S-parameters are

    $S^{\prime}_{ij} = S_{ij} + \frac{ S_{kj} · S_{il} · ( 1 − S_{lk} ) + S_{lj} · S_{ik} · ( 1 − S_{kl} ) + S_{kj} · S_{ll} · S_{ik} + S_{lj} · S_{kk} · S_{il} } { ( 1 − S_{kl} ) · ( 1 − S_{lk} ) − S_{kk} · S_{ll} }$

2. Connect two ports of two different networks:

    Connecting port $k$ of circuit $S$ with port $l$ of circuit $T$, the new S-parameters are

    $S^{\prime}_{ij} = S_{ij} +  \frac{ S_{kj} · T_{ll} · S_{ik} }{ 1 − S_{kk} · T_{ll} }$

    with $i$ and $j$ both being ports of $S$. Furthermore, it is

    $S^{\prime}_{mj} = \frac{ S_{kj} · T_{ml} }{ 1 − S_{kk} · T_{ll} }$

    with $m$ being a port of the circuit $T$.

S-parameters are "Black Boxes".
No internal state is visible.
To make measurement inside a circuit, the matrices are augmented with ports to output internal signals ( voltage, current ).
Because there is no input from this ports, the resulting matrix gets long but narrow.
This output-only part is held in a different matrix as the normal s-parameters.
Additionally, the lines of this matrix can be labeled with a name, e.g the name of the node + "_V" for voltage.

To simulate a network and measure some values:

1. Construct the networks for the individual components.

2. Augment the matrix for the ports you want to measure with output-only ports.

3. Connect the networks together, until only one port for each source remains.

This results in a network with a small s-parameter matrix, which describes the flow of power between all sources through the circuit.
The additional lines in the measurement matrix describe the measured signals in dependence of the power input from each (real) port.

## API doc

```@autodocs
Modules = [Networks]
```

## Index

```@index
```
