#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>

#include <chealpix.h>

#include <gnuastro/qsort.h>
#include <gnuastro/table.h>
#include <gnuastro/polygon.h>
#include <gnuastro/statistics.h>



# define M_PI		    3.14159265358979323846	/* Somehow M_PI is not working. */
# define DEG2RAD(d) (d)*M_PI/180




struct Map
{
  size_t size;
  size_t object_id[5];
};



struct quad_candidate
{
  double x, y;
  double distance;
  size_t times_used;
};




/* Make the healpix map using the RA and DEC table columns. 
   Ensure RA and DEC are in degree. If not convert them. */
void
create_healpixmap(double *ra,
                  double *dec,
                  float  *magnitude,
                  size_t col_size,
                  size_t ntop_objects,
                  long   nside,
                  size_t *npolygon,
                  float *polygon)
{
  float *signal=NULL;
  size_t *ordinds=NULL;
  size_t i=0, j=0, k=0;
  double *temp_polygon=NULL;
  size_t *objectid_arr=NULL;
  struct Map *most_bright=NULL;
  long ipring, maxipring=LONG_MIN, minipring=LONG_MAX;


  /* Allocate the arrays to use the ordinnds arry to store the
     indexes and the temp_polygon array to store the points
     to be sorted while sorting the polygons. */
  ordinds=malloc(col_size*ntop_objects*sizeof(*ordinds));

  temp_polygon=malloc(col_size*2*ntop_objects*sizeof(*temp_polygon));

  /* Allocate the columns in the map whcih is same as
      total no of HEALpixes. */
  most_bright=malloc(col_size*sizeof(*most_bright));

  for(i=0;i<col_size;i++)
    {
      for(j=0;j<ntop_objects;j++)
            most_bright[i].object_id[j]=0;

      most_bright[i].size=0;
    }


  /* Allocate signal array. */
  signal=malloc(col_size*sizeof(*signal));

  /* Allocate object id array and initialize. */
  objectid_arr=malloc(col_size*sizeof(*objectid_arr));
  for(i=0;i<col_size; ++i)
    objectid_arr[i]=i;


  /* Use brightness column as a reference to sort object id array
      in incresing order.  */
  gal_qsort_index_single=magnitude;
  qsort(objectid_arr, col_size, sizeof(size_t), gal_qsort_index_single_float32_d);


  /* Allocate size for most bright array containing
      5 most brightest start brightness where index of array is
      the ring index no. of the HEALPix and the value is the size/index
      of the columns and their value at that size/index */

  for(i=0, k=0; i<col_size; ++i)
    {
      double ptheta = DEG2RAD(90-dec[objectid_arr[i]]);
      double pphi   = DEG2RAD(ra[objectid_arr[i]]);

      ang2pix_ring(nside, ptheta, pphi, &ipring);
      // printf("ring = %ld\n", ipring);

      if(maxipring < ipring) maxipring = ipring;
      if(minipring > ipring) minipring = ipring;


      if(most_bright[ipring].size < ntop_objects)
        {
          most_bright[ipring].object_id[most_bright[ipring].size] = i;
          
          /* If even 1 star is present, use it's id in map formation. */
          if(most_bright[ipring].size == 0)
            signal[k++]=(float)ipring;


          // printf("ring-index = %ld, size = %ld, object-index = %ld \n",
          //                         ipring
          //                         , most_bright[ipring].size+1
          //                         , most_bright[ipring].object_id[most_bright[ipring].size]);

          most_bright[ipring].size++; /* Increase the size of the array. */
        }
    }
  
  // printf("max = %zu , min = %zu\n", maxipring, minipring);

  i=0;
  /* Build a polygon array with first `ntop_objects` as a single polyon,
     the next `ntop_objects` another and so on. */
  for(k=0; signal[k]!=0; k++)
    {
      // i=29126;/* For test only. */
      for(j=0; j<ntop_objects; ++j)
      {
        size_t ii = most_bright[(size_t)signal[k]].object_id[j];
        temp_polygon[j*2]=ra[ii];
        temp_polygon[j*2+1]=dec[ii];
        // printf("temp_polygon = \t%lf %lf\n", temp_polygon[j*2], temp_polygon[j*2+1]);
      }

      gal_polygon_vertices_sort(temp_polygon, ntop_objects, ordinds);

      for(j=0; j<ntop_objects; ++j)
      {
        // printf("ordinds %f\n", signal[k]);
        polygon[i*2]=temp_polygon[ordinds[j]*2];
        polygon[i*2+1]=temp_polygon[ordinds[j]*2+1];
        // printf("sorted polygon = %f, %f\n", polygon[i*2], polygon[i*2+1]);
        i++;
      }
    }
  
  *npolygon=i;

  /* Make/Write a healpix of the magnitude data. Interpolates
      missing values. */
  write_healpix_map(signal, nside, "healpix-test.fits", 0, "C");


  /* Free refernce data. */
  free(objectid_arr);
  free(signal);
  free(most_bright);
  free(temp_polygon);
  free(ordinds);
}





/* Make a ds9 complatiable region file. Internal use only. */
void
make_regionfile(char *filename,
                size_t npolygon,
                float *polygon)
{
  FILE *fileptr;
  size_t i=0,j=0, n=0;

  fileptr = fopen(filename, "w+");
  
  fprintf(fileptr, "# Region file format: DS9 version 4.1\n");
  fprintf(fileptr, "global color=green dashlist=8 3 width=1 "
                    "font=\"helvetica 10 normal roman\" select=1 "
                    "highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 "
                    "include=1 source=1\n");
  fprintf(fileptr, "fk5\n");
  for(i=0; i<npolygon/5; ++i)
    {
      if(polygon[j*2] != 0 && polygon[j*2+1]!=0)
        {
          fprintf(fileptr, "polygon(");
          // i=29126;/* For test only.*/
            for(n=0; n<5; ++n)
            {
              if(j%5==0) fprintf(fileptr, "%lf,%lf", polygon[j*2], polygon[j*2+1]);
              else       fprintf(fileptr, ",%lf,%lf", polygon[j*2], polygon[j*2+1]);
              ++j;
            }
          fprintf(fileptr, ")\n");
        }
    }

  /* Close the file. */
  fclose(fileptr);
}






static int
sort_compare(const void *a, const void *b)
{
  struct quad_candidate *p1 = (struct quad_candidate *)a;
  struct quad_candidate *p2 = (struct quad_candidate *)b;

  return ( p1->distance==p2->distance
           ? 0
           : (p1->distance<p2->distance ? -1 : 1) );
}





void 
make_dist_array(float *polygon,
                size_t npolygon,
                long nside)
{
  double dist;
  size_t i, j, k;
  double theta_pix;
  struct quad_candidate *quad_array=NULL;

  quad_array=malloc(npolygon*sizeof(*quad_array));

  for(i=0;i<npolygon;i++)
    {
      quad_array[i].x=0;
      quad_array[i].y=0;
      quad_array[i].distance=0;
      quad_array[i].times_used=0;
    }

  theta_pix=(180/M_PI)*sqrt(M_PI/3)/nside; /* IN degrees. */
  // printf("theta_pix = %lf\n", theta_pix);

  /*TODO
    The below algorithm will take a lot of time and space, O(n^2).
    Use heaps instead, O(lon(n)).*/

  for(i=0; i<npolygon; ++i)
    {
      for(j=0; j<npolygon; ++j)
        {
          dist=(polygon[2*j]-polygon[2*i])*(polygon[2*j]-polygon[2*i])
                + (polygon[2*j+1]-polygon[2*i+1])*(polygon[2*j+1]-polygon[2*i+1]);
          
          if(dist <= theta_pix && quad_array[i].times_used < 16)
            {
              quad_array[j].x=polygon[2*j];
              quad_array[j].y=polygon[2*j+1];
              quad_array[j].distance=dist;
              quad_array[j].times_used++;
              // printf("%lf, %lf, %lf, %zu\n", quad_array[j].x,
              //                                quad_array[j].y,
              //                                quad_array[j].distance,
              //                                quad_array[j].times_used); 

            }
        }

      qsort(quad_array, npolygon, sizeof(struct quad_candidate), sort_compare);

      for(k=0;k<npolygon;++k)
        if(quad_array[k].distance != 0)
          printf("%lf, %lf, %lf, %zu\n", quad_array[k].x,
                                         quad_array[k].y,
                                         quad_array[k].distance,
                                         quad_array[k].times_used);

      printf("======== next star =======\n");

      for(j=0; j<npolygon; ++j)
        {
          quad_array[j].x=0;
          quad_array[j].y=0;
          quad_array[j].distance=0;
          quad_array[j].times_used=0;
        }
    }

  free(quad_array);

}





int main(){
  long nside=64;
  float *polygon;
  size_t npolygon=0;
  // size_t i=0, k=0;
  // struct Map *most_bright=NULL;
  // float *signal=NULL;
  // size_t *objectid_arr=NULL;


  /* Choose columns to read. */
  gal_list_str_t *cols=NULL;
  gal_list_str_add(&cols, "ra", 0);
  gal_list_str_add(&cols, "dec", 0);
  gal_list_str_add(&cols, "phot_g_mean_mag", 0);
  gal_list_str_reverse(&cols);

  /* Read the columns. */
  gal_data_t *ref=gal_table_read("gaia-in-img.fits", "1",
                                  NULL, cols, GAL_TABLE_SEARCH_NAME,
                                  0, -1, 0, NULL);

  /* Seperate columns. */
  gal_data_t *ra= gal_data_copy_to_new_type (ref, GAL_TYPE_FLOAT64);
  gal_data_t *dec=gal_data_copy_to_new_type (ref->next, GAL_TYPE_FLOAT64);
  gal_data_t *mag=gal_data_copy_to_new_type (ref->next->next, GAL_TYPE_FLOAT32);


  /* Check range. */
  /* gal_data_t *mag_max=gal_statistics_maximum (mag);
  gal_data_t *mag_min=gal_statistics_minimum (mag);

  float *max_c3=mag_max->array;
  float *min_c3=mag_min->array;
  */

  /* Make an array of the columns. */
  double *c1=ra->array;
  double *c2=dec->array;
  float  *c3=mag->array;

  /* Allocate space for the polygon array. */
  polygon=malloc(ref->dsize[0]*sizeof(*polygon));

  /* Make a HEALPix map. OPtimal value of nside is used according to 
     users requirement. */
  create_healpixmap(c1, c2, c3, ref->dsize[0], 5, nside, &npolygon, polygon);

  /* Create a ds9 file for visualisation. */
  make_regionfile("polygon.reg", npolygon, polygon);

  make_dist_array(polygon, npolygon, nside);

  /* Free and return. */
  free(polygon);
  gal_list_data_free (ref);

  return 0;
}