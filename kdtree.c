#include <stdio.h>
#include <math.h>
#include "kdtree.h"




/* Swap 2 nodes of the tree. */
void
swap(struct kd_node *x, struct kd_node *y, size_t dim)
{
  size_t i;
  double tmp;

  /* For all points in the node. */
  for(i=0; i<dim; ++i)
    {
      tmp=x->point[i];
      x->point[i]=y->point[i];
      y->point[i]=tmp;
    }
}





/* Find the distance between 2 nodes of the tree.

   return : distance(squared) between 2 nodes of the tree.
*/
double
find_distance(struct kd_node *a, struct kd_node *b, size_t dim)
{
  size_t i;
  double t_dis, dis = 0;

  /* For all points in the node. */
  for(i=0; i<dim; ++i)
    {
      t_dis = a->point[i] - b->point[i];

      /* As only the relative magnitude of distance is required,
         it is efficient to return distance square rather than using
         `sqrt()` function to calculate actual distace value. */
      dis += t_dis * t_dis;
    }

  return dis;
}





/* Find the median to seperate the hyperspace. Instead of randomly
   chossing the media-point, we use `quickselect alogorithm` to find
   median in average complexity of O(n). This also makes nodes partially
   sorted w.r.t each axis.
   See `https://en.wikipedia.org/wiki/Quickselect` for pseudocode and
   more information of the algorithm.

   return : median point(node) that splits the hyperspace.
*/
struct kd_node*
find_median(struct kd_node *left, struct kd_node *right, size_t current_axis)
{
  double pivotValue;
  struct kd_node *p=NULL, *storeIndex=NULL, *mid=NULL;

  /* False state, return null. */
  if (right <= left) return NULL;

  /* If the tree contains only one element, return that element. */
  if (right == left + 1) return left;

  /* The middle value(here used as pivot) between left and right. */
  mid=left + (right - left) / 2;

  /* Loop until the median(for the current axis) is returned. */
  while(1)
    {
      storeIndex = left;
      pivotValue = mid->point[current_axis];

      /* Select a pivotIndex between left and right. */
      swap(mid, right - 1, MAX_DIM);

      for (p = left; p < right; ++p)
        if (p->point[current_axis] < pivotValue)
          {
            if (p != storeIndex)
                swap(p, storeIndex, MAX_DIM);

            /* Increase the index of partition(pivot). */
            storeIndex++;
          }

      /* Move pivot to its final place. */
      swap(storeIndex, right - 1, MAX_DIM);

      /* If median is found, return it. */
      if (storeIndex->point[current_axis] == mid->point[current_axis])
        return mid;

      /* The pivot is now in its final sorted position. Assign left
         or right based on its current value w.r.t mid. */
      if (storeIndex > mid) right = storeIndex;
      else                  left  = storeIndex;
    }
}





/* Make a kd-tree from a given set of points.
   For tree construction, a median point is selected for each
   axis and the left and right branches are made by comparing
   points based on that axis.

   return : a balanced kd-tree. */
struct kd_node*
make_tree(struct kd_node *current_node, size_t total_nodes, size_t current_axis, size_t dim)
{
  struct kd_node *median_node=NULL;

  /* If no more nodes are present, stop the branch creation. */
  if (!total_nodes) return NULL;

  /* Find the median node for the current axis. */
  median_node = find_median(current_node, current_node + total_nodes, current_axis);

  /* If median node is present, recursively make left and right
     subtree. */
  if(median_node)
    {
      current_axis = (current_axis + 1) % dim;
      median_node->left  = make_tree(current_node,
                                     median_node - current_node,
                                     current_axis, dim);
      median_node->right = make_tree(median_node + 1,
                                     current_node + total_nodes - (median_node + 1),
                                     current_axis, dim);
    }

  return median_node;
}





/* global variable alert! Remove after testing. */
size_t visited;

/* Find the nearest neighbour of the `point`.
   See `https://en.wikipedia.org/wiki/K-d_tree#Nearest_neighbour_search`
   for more information.

  return:
  nearest : The nearest node to the search point.
  least_dist : The distance from the search point to the nearest point.
  */
void
nearest_neighbour(struct kd_node *root, struct kd_node *point, size_t current_axis, size_t dim,
                  struct kd_node **nearest, double *least_dist)
{
  double d, dx, dx2;

  /* If no tree/subtree present, don't search further. */
  if(!root) return;

  /* The distance between search point to the current nearest.*/
  d = find_distance(root, point, dim);

  /* Distance between the splitting coordinate of the search
  point and current node*/
  dx = root->point[current_axis] - point->point[current_axis];
  dx2 = dx*dx;

  visited++;

  /* Check if the current node is nearer than the previous
     nearest node. */
  if(!*nearest || d < *least_dist)
    {
      *least_dist = d;
      *nearest = root;
    }

  /* If exact match found(least distance 0), return it. */
  if(!*least_dist) return;

  current_axis = (current_axis + 1) % dim;

  /* Recursively search in subtrees. */
  nearest_neighbour(dx > 0 ? root->left : root->right, point, current_axis, dim, nearest, least_dist);

  /* Since the hyperplanes are all axis-aligned, to check
  if there is a node in other branch that is nearer to the
  search node is done by a simple comparison to see whether the
  distance between the splitting coordinate of the search
  point and current node is lesser(i.e on same side of hyperplane)
  than the distance (overall coordinates) from the search point to
  the current nearest. */
  if(dx2 >= *least_dist) return;
  nearest_neighbour(dx > 0 ? root->right : root->left, point, current_axis, dim, nearest, least_dist);
}




int main()
{
  double least_dist;
  struct kd_node *root=NULL, *nearest=NULL;
  struct kd_node testNode = {{9, 6, 1, 4}};

  struct kd_node tree[] = {
      {{2, 3, 1, 4}}, {{5, 4, 2, 2}}, {{9, 6, 5, 3}}, {{4, 7, 2, 2}}, {{8, 1, 2, 4}}, {{7, 2, 1, 4}}
  };

  root = make_tree(tree, sizeof(tree) / sizeof(tree[1]), 0, 4);

  visited = 0;
  nearest_neighbour(root, &testNode, 0, 4, &nearest, &least_dist);

  printf(">> KD-tree\nsearching for (%g, %g, %g, %g)\n"
          "nearest (%g, %g, %g, %g) find_distance %g\nseen %zu nodes\n",
          testNode.point[0], testNode.point[1],testNode.point[2], testNode.point[3],
          nearest->point[0], nearest->point[1],nearest->point[2], nearest->point[3], sqrt(least_dist), visited);


  return 0;
}