#include <stdio.h>
#include <stdlib.h>

#include <chealpix.h>

#include <gnuastro/table.h>
#include <gnuastro/statistics.h>


int main(){

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

    /* Free refernce data. */
    gal_list_data_free (ref);

    /* Find range. */
    gal_data_t *ra_min=gal_statistics_minimum (ra);
    gal_data_t *ra_max=gal_statistics_maximum (ra);

    /* Make an array of the columns. */
    double *c1=ra->array;
    double *c2=dec->array;
    float  *c3=mag->array;

    double *min_c1=ra_min->array;
    double *max_c1=ra_max->array;

    /* Make a healpix of the magnitude data. */
    write_healpix_map(c3, 2, "healpix-test.fits", 0, "C");

    printf("%lf %lf\n", min_c1[0], max_c1[0]);

    printf("%lf %lf %f\n", c1[0], c2[0], c3[0]);

    return 0;
}