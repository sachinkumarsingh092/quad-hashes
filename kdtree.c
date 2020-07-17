#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>


#define MAX_DIM 4



struct kd_node
{
  double point[MAX_DIM];
  struct kd_node *left, *right;
};





void 
swap(struct kd_node *x, struct kd_node *y, size_t dim)
{
  size_t i;
  double tmp;

  for(i=0; i<dim; ++i)
    {
      tmp=x->point[i];
      x->point[i]=y->point[i];
      y->point[i]=tmp;
    }
}





double
distance(struct kd_node *a, struct kd_node *b, size_t dim)
{
  size_t i;
  double t_dis, dis = 0;

  for(i=0; i<dim; ++i)
    {
      t_dis = a->point[i] - b->point[i];
      dis += t_dis * t_dis;
    }

  return dis;
}





/* see quickselect method */
struct kd_node*
find_median(struct kd_node *start, struct kd_node *end, size_t idx)
{
  if (end <= start) 
    return NULL;
  if (end == start + 1)
    return start;

  struct kd_node *p, *store, *md = start + (end - start) / 2;
  double pivot;
  while (1)
    {
      pivot = md->point[idx];

      swap(md, end - 1, MAX_DIM);
      for (store = p = start; p < end; p++) {
          if (p->point[idx] < pivot) {
              if (p != store)
                  swap(p, store, MAX_DIM);
              store++;
          }
      }
      swap(store, end - 1, MAX_DIM);

      /* median has duplicate values */
      if (store->point[idx] == md->point[idx])
        return md;

      if (store > md) end = store;
      else            start = store;
    }
}





struct kd_node*
make_tree(struct kd_node *t, size_t len, size_t i, size_t dim)
{
  struct kd_node *n;

  if (!len) return 0;

  if ((n = find_median(t, t + len, i)))
    {
      i = (i + 1) % dim;
      n->left  = make_tree(t, n - t, i, dim);
      n->right = make_tree(n + 1, t + len - (n + 1), i, dim);
    }

  return n;
}





/* global variable alert! */
size_t visited;

void nearest(struct kd_node *root, struct kd_node *nd, size_t i, size_t dim,
        struct kd_node **best, double *best_dist)
{
  double d, dx, dx2;

  if (!root) return;
  d = distance(root, nd, dim);
  dx = root->point[i] - nd->point[i];
  dx2 = dx * dx;

  visited ++;

  if (!*best || d < *best_dist)
    {
      *best_dist = d;
      *best = root;
    }

  /* if chance of exact match is high */
  if (!*best_dist) return;

  if (++i >= dim) i = 0;

  nearest(dx > 0 ? root->left : root->right, nd, i, dim, best, best_dist);
  if (dx2 >= *best_dist) return;
  nearest(dx > 0 ? root->right : root->left, nd, i, dim, best, best_dist);
}





