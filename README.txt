Brushfire algorythm
The brushfire algorithm, as described by Verwer., computes static distance maps with a
shortest path search similar to Dijkstra’s algorithm with multiple source.
By using a priority queue that orders the expansion of cells by the distance to their closest
obstacle, the propagation spreads in wavefronts that start at the location of obstacles.
The occupied cells are initialized with zero
distances and then inserted into the priority queue open: the function insert(open, s, d) inserts
s into the queue with distance d, or updates the priority if s is already enqueued.
As long as this queue contains cells, the algorithm iteratively calls pop(open), which returns the
cell s with the lowest enqueued distance and removes it from the queue. It then updates
the cells in the 8-connected neighborhood Adj8 (s) of s: if the distance d from a neighbor n to
the closest obstacle of s as specified by obst(s) is smaller than the current value D(n),
the distance value and closest obstacle of n are updated with the obstacle of s.
Furthermore, each updated neighbor cell n is inserted into the priority queue with its new distance value to continue the propagation 
