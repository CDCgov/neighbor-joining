neighborjoining is a javascript library for creating phylogenetic trees from
distance matrices. It is based on the [Rapid Neighbour-Joining](http://pure.au.dk/ws/files/19821675/rapidNJ.pdf)
algorithm.

# Installation

To install neighborjoining package with NPM use: `npm install neighborjoining`

# Usage

```javascript
var RNJ = new neighborjoining(D, taxa, copyDistanceMatrix, taxonIdAccessor);
```

Description of arguments used in initialization:
* **D** - distance matrix (two dimensional array of size NxN)
* **taxa** - array of size N with taxa data (such as strings or objects). Element (taxon) at first index corresponds to first row/column of the distance matrix and so on.
* **copyDistanceMatrix** - flag specifying whether D can be modified or must be copied. Copying might minimally increase the initialization time. Default: `false`.
* **taxonIdAccessor** - function for retrieving the name/identifier of a taxon. It is only called during Newick tree creation. Default: `function(t){ return t.name }`.

# Example
```javascript
var D = [
    [0,  5,  9,  9, 8],
    [5,  0, 10, 10, 9],
    [9, 10,  0,  8, 7],
    [9, 10,  8,  0, 3],
    [8,  9,  7,  3, 0]
];
var taxa = [
    {
        name: "A",
        genotype: "g1"
    },
    {
        name: "B",
        genotype: "g2"
    },
    {
        name: "C",
        genotype: "g3"
    },
    {
        name: "D",
        genotype: "g4"
    },
    {
        name: "E",
        genotype: "g5"
    }
];
var RNJ = new neighborjoining(D, taxa);
```

So, what can we do with it? If you're interested in the raw data, you can see the entire tree:

```javascript
RNJ.getAsObject();

{
    "taxon": null,
    "length": null,
    "children": [{
        "taxon": {
            "name": "C",
            "genotype": "g3"
        },
        "length": 2,
        "children": []
    }, {
        "taxon": null,
        "length": 2,
        "children": [{
            "taxon": null,
            "length": 3,
            "children": [{
                "taxon": {
                    "name": "A",
                    "genotype": "g1"
                },
                "length": 2,
                "children": []
            }, {
                "taxon": {
                    "name": "B",
                    "genotype": "g2"
                },
                "length": 3,
                "children": []
            }]
        }, {
            "taxon": null,
            "length": 2,
            "children": [{
                "taxon": {
                    "name": "D",
                    "genotype": "g4"
                },
                "length": 2,
                "children": []
            }, {
                "taxon": {
                    "name": "E",
                    "genotype": "g5"
                },
                "length": 1,
                "children": []
            }]
        }]
    }]
}
```

Alternately, if you just need the Newick string...

```javascript
var distances = false;
var round = false;
var p1 = false;
RNJ.getAsNewick(distances, round, p1);
// (C,((A,B),(D,E)));
```

`getAsNewick` takes three arguments:
 - `distances` - Should the output include distances?
 - `round` - Should distances be rounded to the nearest integers?
 - `p1` - Should distances all have one added to them? (Useful if zeros break downstream normalization efforts)

So...

```javascript
RNJ.getAsNewick(true);
// (C:2,((A:2,B:3):3,(D:2,E:1):2):2);

RNJ.getAsNewick(true, true);
// (C:2,((A:2,B:3):3,(D:2,E:1):2):2);

RNJ.getAsNewick(true, true, true);
// (C:3,((A:3,B:4):4,(D:3,E:2):3):3);
```
