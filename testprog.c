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
                  float *polygon)
{
  size_t *ordinds;
  float *signal=NULL;
  size_t i=0, j=0, k=0;
  double *temp_polygon;
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

  /* Allocate signal array. */
  signal=malloc(col_size*sizeof(*most_bright));

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

      // printf("max = %zu , min = %zu\n", maxipring, minipring);

      if(most_bright[ipring].size < ntop_objects)
        {
          most_bright[ipring].object_id[most_bright[ipring].size] = i;

          if(most_bright[ipring].size == 0)
            signal[k++]=(float)ipring;


          // printf("ring-index = %ld, size = %ld, object-index = %ld \n",
          //                         ipring
          //                         , most_bright[ipring].size+1
          //                         , most_bright[ipring].object_id[most_bright[ipring].size]);

          most_bright[ipring].size++; /* Increase the size of the array. */
        }
    }

  k=0;
  /* Build a polygon array with first `ntop_objects` as a single polyon,
     the next `ntop_objects` another and so on. */
  for(i=minipring; i<=maxipring; ++i)
    {
      // i=29126;/* For test only. */
      for(j=0; j<ntop_objects; ++j)
      {
        size_t ii = most_bright[i].object_id[j];
        temp_polygon[j*2]=ra[ii];
        temp_polygon[j*2+1]=dec[ii];
        // printf("temp_polygon = \t%lf %lf\n", temp_polygon[j*2], temp_polygon[j*2+1]);
      }

      gal_polygon_vertices_sort(temp_polygon, ntop_objects, ordinds);

      for(j=0; j<ntop_objects; ++j)
      {
        // printf("ordinds %zu\n", ordinds[j]);
        polygon[k*2]=temp_polygon[ordinds[j]*2];
        polygon[k*2+1]=temp_polygon[ordinds[j]*2+1];
        // printf("sorted polygon = %f, %f\n", temp_polygon[ordinds[j]*2], temp_polygon[ordinds[j]*2+1]);
        k++;
      }
    }
  

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
  for(i=0; i<npolygon; ++i)
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

  /* Close the file. */
  fclose(fileptr);
}




void 
make_dist_array(float *polygon,
                     long nside, 
                     size_t total_stars)
{
  double dist;
  size_t i, j;
  double theta_pix;
  struct quad_candidate *quad_array;

  quad_array=malloc(total_stars*sizeof(quad_array));

  theta_pix=(180/M_PI)*sqrt(M_PI/3)/nside; /* IN degrees. */
  // printf("theta_pix = %lf\n", theta_pix);


  for(i=0; i<total_stars; ++i)
    for(j=0; j<total_stars; ++j)
      {
        dist=(polygon[2*j]-polygon[2*i])*(polygon[2*j]-polygon[2*i])
              + (polygon[2*j+1]-polygon[2*i+1])*(polygon[2*j+1]-polygon[2*i+1]);
        
        if(dist <= theta_pix && quad_array[i].times_used < 16)
          {
            quad_array[i].x=polygon[2*j];
            quad_array[i].y=polygon[2*j+1];
            quad_array[i].distance=dist;
            quad_array[i].times_used++;
            printf("%lf, %lf, %lf, %zu\n", quad_array[i].x,
                                           quad_array[i].y,
                                           quad_array[i].distance,
                                           quad_array[i].times_used); 
          }
      }
  
    
}





int main(){
  long nside=100;
  float *polygon;
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
  create_healpixmap(c1, c2, c3, ref->dsize[0], 5, nside, polygon);

  /* Create a ds9 file for visualisation. */
  make_regionfile("polytest.reg", 1000, polygon);

  make_dist_array(polygon, 64, 1000);

  /* Free and return. */
  free(polygon);
  gal_list_data_free (ref);

  return 0;
}