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
# define RAD2DEG(d) (d)*180/M_PI




struct Map
{
  size_t size;
  // Use a FAM for object_id
  size_t object_id[5];
};



struct quad_candidate
{
  double x, y;
  double distance;
  float brightness;
  size_t times_used;
};







/* Make a ds9 complatiable region file. Internal use only. 
   Rename it to gal_polygon_to_ds9reg. Give color argument.*/
void
gal_polygon_to_ds9reg(char *filename,
                size_t npolygon,
                float *polygon,
                size_t npoints,
                char *color)
{
  FILE *fileptr;
  size_t i=0,j=0, n=0;

  if (!color) color = "green";

  fileptr = fopen(filename, "w+");

  fprintf(fileptr, "# Region file format: DS9 version 4.1\n");
  fprintf(fileptr, "global color=%s dashlist=8 3 width=1 "
                    "font=\"helvetica 10 normal roman\" select=1 "
                    "highlite=1 dash=0 fixed=0 edit=1 move=1 delete=1 "
                    "include=1 source=1\n", color);
  fprintf(fileptr, "fk5\n");
  for(i=0; i<npolygon; ++i) 
    {
      if(polygon[j*2]!=0 && polygon[j*2+1]!=0)
        {
          fprintf(fileptr, "polygon(");
            for(n=0; n<npoints; ++n)
              {
                if(j%npoints==0) fprintf(fileptr, "%lf,%lf", polygon[j*2], polygon[j*2+1]);
                else             fprintf(fileptr, ",%lf,%lf", polygon[j*2], polygon[j*2+1]);
                ++j;
              }
          fprintf(fileptr, ")\n");
        }
    }

  /* Close the file. */
  fclose(fileptr);
}



/* Make the healpix map using the RA and DEC table columns.
   Ensure RA and DEC are in degree. If not convert them. */
void
create_healpixmap(double *ra,
                  double *dec,
                  float  *magnitude,
                  size_t col_size,
                  size_t ntop_objects,
                  long   nside,
                  struct Map *most_bright,
                  float *healpix_id)
{
  long ipring;
  size_t i=0, k=0;
  size_t *objectid_arr=NULL;
  float *healpix_polygon=NULL;


  /* Allocate object id array and initialize. */
  objectid_arr=malloc(col_size*sizeof(*objectid_arr));
  for(i=0;i<col_size; ++i)
    objectid_arr[i]=i;

  healpix_polygon=malloc(col_size*sizeof(*healpix_polygon));

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

      // make a function to find center of healpix.
      double theta, phi;
      double rc, dc;
      double theta_pix=60*sqrt(3/M_PI)/nside;

      if(most_bright[ipring].size < ntop_objects)
        {
          most_bright[ipring].object_id[most_bright[ipring].size] = i;

          /* If even 1 star is present, use it's id in map formation. */
          if(most_bright[ipring].size == 0)
            {
              healpix_id[k]=(float)ipring;

              pix2ang_ring(nside, ipring, &theta, &phi);

              rc = RAD2DEG(phi);
              dc = 90-RAD2DEG(theta);

              printf("ra = %g, dec = %g\n", rc, dc);

              healpix_polygon[8*k]   = rc;
              healpix_polygon[8*k+1] = dc-theta_pix/2;

              healpix_polygon[8*k+2] = rc-theta_pix/2;
              healpix_polygon[8*k+3] = dc;

              healpix_polygon[8*k+4] = rc;
              healpix_polygon[8*k+5] = dc+theta_pix/2;

              healpix_polygon[8*k+6] = rc+theta_pix/2;
              healpix_polygon[8*k+7] = dc;

              k++;
            }


          // printf("ring-index = %ld, size = %ld, object-index = %ld \n",
          //                         ipring
          //                         , most_bright[ipring].size+1
          //                         , most_bright[ipring].object_id[most_bright[ipring].size]);

          most_bright[ipring].size++; /* Increase the size of the array. */
        
        }
    }

  gal_polygon_to_ds9reg("healpix-diamond.reg", k, healpix_polygon, 4, "red");

  /* Make/Write a healpix of the magnitude data. Interpolates
      missing values. */
  write_healpix_map(healpix_id, nside, "healpix-test.fits", 0, "C");

  /* Free refernce data. */
  free(healpix_polygon);
  free(objectid_arr);
}






/* Sort in ascending order. */
static int
sort_distance(const void *a, const void *b)
{
  struct quad_candidate *p1 = (struct quad_candidate *)a;
  struct quad_candidate *p2 = (struct quad_candidate *)b;

  return ( p1->distance==p2->distance
           ? 0
           : (p1->distance<p2->distance ? -1 : 1) );
}



static int
sort_brightness(const void *a, const void *b)
{
  struct quad_candidate *p1 = (struct quad_candidate *)a;
  struct quad_candidate *p2 = (struct quad_candidate *)b;

  return ( p1->brightness==p2->brightness
           ? 0
           : (p1->brightness<p2->brightness ? -1 : 1) );
}



void
make_dist_array(double *ra,
                double *dec,
                float  *magnitude,
                size_t ntop_objects,
                size_t col_size,
                long   nside,
                float  *healpix_id,
                struct Map *most_bright)
{
  double dist_sq;
  double theta_pix;
  size_t npolygon=0;
  size_t nquad_arr=0;
  size_t nquads_num=0;
  float *brightness_arr_plot=NULL;
  struct quad_candidate *polygon=NULL;
  size_t i=0, j=0, k=0, b_ind=0, ii=0;
  struct quad_candidate *quad_cand=NULL;
  struct quad_candidate brightsorted_quad[4]={0};
  
  brightness_arr_plot=malloc(col_size*sizeof(*brightness_arr_plot));

  /* Allocate space for the polygon array. */
  polygon=malloc(col_size*sizeof(*polygon));
  for(i=0;i<col_size;++i)
    {
      polygon[i].x=0;
      polygon[i].y=0;
      polygon[i].brightness=0;
      polygon[i].distance=0;
      polygon[i].times_used=0;
    }

  theta_pix=60*sqrt(3/M_PI)/nside; /* IN degrees. */
  // printf("theta_pix = %lf\n", theta_pix);

  /* TODO
    The below algorithm will take a lot of time and space, O(n^2).
    Use heaps instead, O(log(n)).*/

  /* Create indexed polygon object. */
  for(i=0, k=0; healpix_id[k]!=0; k++)
    for(j=0; j<ntop_objects; ++j)
      {
        ii=most_bright[(size_t)healpix_id[k]].object_id[j];

        polygon[i].x=ra[ii];
        polygon[i].y=dec[ii];
        polygon[i].brightness=magnitude[ii];

        i++;
      }

  npolygon=i;

  quad_cand=malloc(npolygon*sizeof(*quad_cand));
  for(i=0;i<npolygon;i++)
    {
      quad_cand[i].x=0;
      quad_cand[i].y=0;
      quad_cand[i].distance=0;
      quad_cand[i].brightness=0;
      quad_cand[i].times_used=0;
    }

  printf("npolygon = %zu\n", npolygon);

  for(i=0; i<npolygon; ++i)
    {
      /* Using k starts the array form 0 rather than any value on j
         would satisfy(maybe somewhere in the middle of pollygon array for some j)*/
      k=0;
      nquad_arr=0;

      for(j=0; j<npolygon; ++j)
        {

          if( polygon[j].x <= polygon[i].x + theta_pix/2 &&
              polygon[j].x >= polygon[i].x - theta_pix/2  )
            {
                if( polygon[j].y <= polygon[i].y + theta_pix/2 &&
                    polygon[j].y >= polygon[i].y - theta_pix/2 &&
                    polygon[j].times_used <= 10)
                  {
                    /* distance^2 = x*x+y*y */
                    dist_sq=(polygon[j].x-polygon[i].x)*(polygon[j].x-polygon[i].x )
                            +(polygon[j].y-polygon[i].y)*(polygon[j].y-polygon[i].y);

                    quad_cand[k].x=polygon[j].x;
                    quad_cand[k].y=polygon[j].y;
                    quad_cand[k].brightness=polygon[j].brightness;
                    quad_cand[k].distance=dist_sq;
                    quad_cand[k].times_used=polygon[j].times_used++;

                    printf("%lf, %lf, %f, %lf, %zu\n", quad_cand[k].x,
                                                   quad_cand[k].y,
                                                   quad_cand[k].brightness,
                                                   quad_cand[k].distance,
                                                   quad_cand[k].times_used);

                    k++;
                    nquad_arr++;
                  }
            }
        }

      printf("%ld\n", nquad_arr);

      /* Sort them on basis of distance and take the n, n-2, n-4 stars. */
      qsort(quad_cand, nquad_arr, sizeof(struct quad_candidate), sort_distance);

      if(nquad_arr >= 6)
        {
          brightsorted_quad[0] = quad_cand[nquad_arr-1]; 
          brightsorted_quad[1] = quad_cand[nquad_arr-3];
          brightsorted_quad[2] = quad_cand[nquad_arr-5];
          brightsorted_quad[3] = quad_cand[0];

          nquads_num++;
        }
      

      /* Sort according to brightness. */
      qsort(brightsorted_quad, 4, sizeof(struct quad_candidate), sort_brightness);

      // for(j=0;j<nquad_arr;++j)
      //   if(quad_cand[j].distance != 0)
          // printf("%lf, %lf, %f, %lf, %zu\n", quad_cand[j].x,
          //                                    quad_cand[j].y,
          //                                    quad_cand[j].brightness,
          //                                    quad_cand[j].distance,
          //                                    quad_cand[j].times_used);
      
      for(j=0; j<4 && nquad_arr>=6; ++j)
        {
          brightness_arr_plot[2*b_ind]=brightsorted_quad[j].x;
          brightness_arr_plot[2*b_ind+1]=brightsorted_quad[j].y;

          // printf("%lf, %lf, %f, %lf, %zu\n", brightsorted_quad[j].x,
          //                                    brightsorted_quad[j].y,
          //                                    brightsorted_quad[j].brightness,
          //                                    brightsorted_quad[j].distance,
          //                                    brightsorted_quad[j].times_used);
          
          // printf("brightness_arr_plot = %ld %lf %lf\n",j, brightness_arr_plot[2*b_ind],
          //                                           brightness_arr_plot[2*b_ind+1]);
          b_ind++;
        }

      /* Take n, n-2, n-4 stars and sort them according to brightness. */


      for(j=0; j<npolygon; ++j)
        {
          quad_cand[j].x=0;
          quad_cand[j].y=0;
          quad_cand[j].brightness=0;
          quad_cand[j].distance=0;
          quad_cand[j].times_used=0;
        }
    
    }
  
  /* Total number of quads. */
  printf("number of quad = %zu\n", nquads_num);

  /* If quad needs to be sorted in point orger sort them using polygon sort.*/

  /* Plot for check. */
  gal_polygon_to_ds9reg("final-quads.reg", nquads_num, brightness_arr_plot, 4, "red");

  free(quad_cand);
  free(polygon);
  free(brightness_arr_plot);

}





int main(){
  long nside=64;
  size_t i=0, j=0;
  float *healpix_id=NULL;
  struct Map *most_bright=NULL;
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


  /* Make an array of the columns. */
  double *c1=ra->array;
  double *c2=dec->array;
  float  *c3=mag->array;

  /* Allocate the columns in the map whcih is same as
    total no of HEALpixes. */
  most_bright=malloc(ref->dsize[0]*sizeof(*most_bright));
  for(i=0;i<ref->dsize[0];i++)
  {
    for(j=0;j<5;j++)
          most_bright[i].object_id[j]=0;

    most_bright[i].size=0;
  }


  /* Allocate healpix_id array. */
  healpix_id=malloc(ref->dsize[0]*sizeof(*healpix_id));


  /* Make a HEALPix map. OPtimal value of nside is used according to
     users requirement. */
  create_healpixmap(c1, c2, c3, ref->dsize[0], 5, nside, most_bright, healpix_id);

  make_dist_array(c1, c2, c3, 5, ref->dsize[0], nside, healpix_id, most_bright);

  /* Free and return. */
  free(healpix_id);
  free(most_bright);
  gal_list_data_free (ref);

  return 0;
}