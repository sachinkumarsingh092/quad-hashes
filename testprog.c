#include <stdio.h>
#include <stdlib.h>

#include <chealpix.h>

#include <gnuastro/table.h>
#include <gnuastro/statistics.h>

double
deg2rad(double degree)
{
  double pi = 3.14159; 
  return(degree * (pi/180));
}

struct Map
{
  size_t size;
  size_t object_id[5];
};

int main(){
    size_t i;
    long nside=34, ipring;
    double theta, phi;
    struct Map *most_bright;
    size_t total_objects=50000; /* Total objects in the index catalog for HEALPix map*/

    /* Allocate the columns in the map whcih is same as
       total no of HEALpixes. */
    most_bright=malloc(total_objects*sizeof(most_bright));


    /* Choose columns to read. */
    gal_list_str_t *cols=NULL;
    gal_list_str_add(&cols, "ra", 0);
    gal_list_str_add(&cols, "dec", 0);
    gal_list_str_add(&cols, "phot_g_mean_mag", 0);
    gal_list_str_reverse(&cols);

    /* Read the columns. */
    gal_data_t *ref=gal_table_read("gaia-dr2-near-castor.fits", "1",
                                    NULL, cols, GAL_TABLE_SEARCH_NAME,
                                    0, -1, 0, NULL);

    /* Seperate columns. */
    gal_data_t *ra= gal_data_copy_to_new_type (ref, GAL_TYPE_FLOAT64);
    gal_data_t *dec=gal_data_copy_to_new_type (ref->next, GAL_TYPE_FLOAT64);
    gal_data_t *mag=gal_data_copy_to_new_type (ref->next->next, GAL_TYPE_FLOAT32);


    /* Find range. */
    gal_data_t *ra_min=gal_statistics_minimum (ra);
    gal_data_t *ra_max=gal_statistics_maximum (ra);

    /* Make an array of the columns. */
    double *c1=ra->array;
    double *c2=dec->array;
    float  *c3=mag->array;

    double *min_c1=ra_min->array;
    double *max_c1=ra_max->array;

    /* Allocate size for most bright array containing
       5 most brightest start brightness where index of array is
       the ring index no. of the HEALPix. */

    for(i=0; i<ref->dsize[0]; ++i)
      {
        double ptheta = deg2rad(90-c2[i]);
        double pphi   = deg2rad(c1[i]);

        ang2pix_ring(nside, ptheta, pphi, &ipring);
        // printf("ring = %ld\n", ipring);

        if(most_bright[ipring].size < 5)
          {
            most_bright[ipring].object_id[most_bright[ipring].size] = i;

            printf("ring-index = %ld size = %ld , object-index = %ld \n", ipring, most_bright[ipring].size,
                            most_bright[ipring].object_id[most_bright[ipring].size]);

            most_bright[ipring].size++; /* Increase the size of the array. */
          }
      }

    // for(i=0; i<3921; ++i)
    //   {
    //     printf("ring-index = %ld\t", i);
    //     for(size_t j=0; j<5; ++j)
    //     {
    //       printf("object-ids[%ld] = %ld ", j, most_bright[i].object_id[j]);
    //     }
    //     printf("\n");
    //   }


    // for(i=0; i<12; ++i)
    //   {
    //     pix2ang_ring(2, i+1, &theta, &phi);
    //     printf("theta = %lf, phi = %lf\n", theta, phi);
    //   }

    /* Make a healpix of the magnitude data. */
    // write_healpix_map(c3, 4, "healpix-test.fits", 0, "C");

    // printf("%lf %lf\n", min_c1[0], max_c1[0]);

    printf("%lf %lf %f\n", c1[0], c2[0], c3[0]);


    /* Free refernce data. */
    gal_list_data_free (ref);
    free(most_bright);

    return 0;
}