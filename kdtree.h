#include <stdio.h>


#define MAX_DIM 4



/* kd_node structure of the tree. */
struct kd_node
{
  double point[MAX_DIM];
  struct kd_node *left, *right;
};




void
swap(struct kd_node *x, struct kd_node *y, size_t dim);

double
find_distance(struct kd_node *a, struct kd_node *b, size_t dim);

struct kd_node*
find_median(struct kd_node *left, struct kd_node *right, size_t current_axis);

struct kd_node*
make_tree(struct kd_node *current_node, size_t total_nodes, size_t current_axis, size_t dim);

void
nearest_neighbour(struct kd_node *root, struct kd_node *point, size_t current_axis, size_t dim,
                  struct kd_node **nearest, double *least_dist);