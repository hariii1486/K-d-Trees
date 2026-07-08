# KD-Tree Implementation (C++)

## Features

### Tree Construction
- Randomly built KD-Tree
- Optimal (balanced) KD-Tree using median-based splitting

### Insertion
- Dynamic insertion of k-dimensional points
- Insertion queries on existing trees

### Searching & Queries
- Exact match query
- Partial match query
- Region (range) query
- Nearest Neighbor (NN) query
- Radius query

### Deletion
- Deletion of points from the KD-Tree
- Dimension-aware restructuring after deletion

---

## Supported Operations

| Operation | Description |
|---------|------------|
| Insert | Insert a k-dimensional point |
| Search | Search for a point in the KD-Tree |
| Exact Match | Find an exact point |
| Partial Match | Match on selected dimensions |
| Region Query | Range search within a hyper-rectangle |
| NN Query | Nearest neighbor search |
| Radius Query | Points within a given radius |
| Delete | Remove a point from the tree |

---

## Implementation Details

- Language: **C++**
- Dimensionality: **Generic k dimensions**
- Distance Metric: **Euclidean distance**
- Tree Type: **KD-Tree**
- Construction Methods:
  - Random insertion
  - Median-based optimal construction


