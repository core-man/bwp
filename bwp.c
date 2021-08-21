/*
 *  Calculate body-wave polarization using RFHV.
 *
 *  Author: Jiayuan Yao @ NTU
 *
 *  Revisions:
 *      2020-08-20  Jiayuan Yao  Initial Coding
 *
 */
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <malloc.h>
#include <math.h>
#include "sacio.h"
#include "const.h"

#define deg2rad(degrees) ((degrees) * M_PI / 180.0)
#define rad2deg(radians) ((radians) * 180.0 / M_PI)

typedef struct initial_vals
{
    char   str[MAX_FNAME];      /* sac file name                  */
    float  pRr;                 /* ray parameter in sec/km        */
    float  tt;                  /* RF begin time (hd.b, negative) */
    float  vs_app;              /* apparent Vs from original RF   */
} IVAL;

void usage(void);
int get_file_line(char *fname);
int read_file(char *fname, int line_num, IVAL *IV);
int read_pers(char *fname, int nbwp, float *Tbwp);
int read_sac_data(IVAL *IV, int len, float **data, int *npts, float *b, float *depmin, float *depmax, float *delta, float *tt, int tmark);
int cal_bwp(IVAL *IV, int len, float **Rf, float *Z_rf, int *npts, float *b, float *depmin, float *depmax, float *dt_Rf, int nbwp, float *Tbwp, float **bwp);

void usage() {
    fprintf(stderr, "calculate body-wave polarization using RF          \n");
    fprintf(stderr, "                                                   \n");
    fprintf(stderr, "Usage:                                             \n");
    fprintf(stderr, "  sacstack [-Ddatalist] [-Zz-rf]                   \n");
    fprintf(stderr, "           [-Tperiodlist] [-Ooutfile]              \n");
    fprintf(stderr, "                                                   \n");
    fprintf(stderr, "Options:                                           \n");
    fprintf(stderr, "  -D: data list                                    \n");
    fprintf(stderr, "  -Z: z-rf                                         \n");
    fprintf(stderr, "  -T: perido file                                  \n");
    fprintf(stderr, "  -O: output file                                  \n");
    fprintf(stderr, "  -h: show usage                                   \n");
}

int main(int argc, char *argv[])
{
    int c;
    int error;
                        /* data list, z-rf file, period list, output file */
    char datalist[MAX_FNAME], zrffile[MAX_FNAME];
    char perlist[MAX_FNAME], outfile[MAX_FNAME];

    SACHEAD hdr;        /* SAC headers */
    int tmark=-3;       /* -5(b), -3(o), -2(a), [0-9] */
    int *npts;
    float *b, *delta, *depmin, *depmax, *tt;

    float **Rf, *Z_rf;  /* R and Z receiver function */
    int   nbwp;         /* period number */
    float *Tbwp;        /* periods */
    float **bwp;        /* body-wave polarization */

    int i, j;
    int len=0;          /* sac file number */
    IVAL *IV;           /* input vars */
    FILE *fp;

    error = 0;
    while ((c = getopt(argc, argv, "D:Z:T:O:h")) != -1) {
        switch(c) {
            case 'D':
                sscanf(optarg, "%s", datalist);
                break;
            case 'Z':
                sscanf(optarg, "%s", zrffile);
                break;
            case 'T':
                sscanf(optarg, "%s", perlist);
                break;
    /*
            case '':
                if (sscanf(optarg, "%d/%f/%f", &tmark, &ts, &tw) != 3)
                    error++;
                break;
                */
            case 'O':
                sscanf(optarg, "%s", outfile);
                break;
            case 'h':
                usage();
                return -1;
            default:
                usage();
                return -1;
        }
    }

    if (argc-optind != 0 || error) {
        usage();
        return -1;
    }


    /*********** Read Data Part ***********/
    len  = get_file_line(datalist); /* get file line number */
    nbwp = get_file_line(perlist);  /* get period number */

    /* set memory for struct vars */
    if ((IV = malloc(sizeof(IVAL)*len)) == NULL) {
        fprintf(stderr, "malloc memory error for IV\n");
        return -1;
    }

    /* read data list into strcut vars */
    len = read_file(datalist, len, IV);

    /* set memory for SAC headers and data */
    if ((Rf = (float**)malloc(sizeof(float*)*len)) == NULL) {
        fprintf(stderr, "Error in allocating memory for Rf\n");
        return -1;
    }
    if ((Tbwp = (float*)malloc(sizeof(float)*nbwp)) == NULL) {
        fprintf(stderr, "Error in allocating memory for npts\n");
        return -1;
    }

    if ((npts = (int*)malloc(sizeof(int)*len)) == NULL) {
        fprintf(stderr, "Error in allocating memory for npts\n");
        return -1;
    }
    if ((b = (float*)malloc(sizeof(float)*len)) == NULL) {
        fprintf(stderr, "Error in allocating memory for b\n");
        return -1;
    }
    if ((delta = (float*)malloc(sizeof(float)*len)) == NULL) {
        fprintf(stderr, "Error in allocating memory for delta\n");
        return -1;
    }
    if ((depmin = (float*)malloc(sizeof(float)*len)) == NULL) {
        fprintf(stderr, "Error in allocating memory for depmin\n");
        return -1;
    }
    if ((depmax = (float*)malloc(sizeof(float)*len)) == NULL) {
        fprintf(stderr, "Error in allocating memory for depmax\n");
        return -1;
    }
    /* tt is not used now */
    if ((tt = (float*)malloc(sizeof(float)*len)) == NULL) {
        fprintf(stderr, "Error in allocating memory for tt\n");
        return -1;
    }

    /* read sac headers and data */
    read_sac_data(IV, len, Rf, npts, b, depmin, depmax, delta, tt, tmark);
    Z_rf = read_sac(zrffile, &hdr);

    /* read periods */
    read_pers(perlist, nbwp, Tbwp);

    /*********** Calculation bwp Part ***********/
    /* set memory for bwp */
    if ((bwp = (float**)malloc(len*sizeof(*bwp))) == NULL) {
        fprintf(stderr, "Error in allocating memory for bwp\n");
        return -1;
    }
    for (i = 0; i < len; i++) {
        if ((bwp[i] = (float*)malloc(nbwp*sizeof(float))) == NULL) {
            fprintf(stderr, "Error in allocating memory for bwp\n");
            return -1;
        }
    }

    /* calculate bwp */
    if (cal_bwp(IV, len, Rf, Z_rf, npts, b, depmin, depmax, delta, nbwp, Tbwp, bwp) < 1) {
        fprintf(stderr, "NO waveforms to be used!!!\n");
    }

    /* output bwp */
    if ((fp = fopen(outfile, "w")) == NULL) {
        fprintf(stderr, "Can not open file %s\n", outfile);
        exit(-1);
    }
    for (i = 0; i < len; i++) {
        fprintf(fp, "> %s %f %f\n", IV[i].str, IV[i].pRr, IV[i].vs_app);
        for (j = 0; j < nbwp; j++) {
           fprintf(fp, "%f %f\n", Tbwp[j], bwp[i][j]);
           //fprintf(stdout, "period: %f  Vs,app: %f\n", Tbwp[j], bwp[i][j]);
        }
    }
    fclose(fp);


    /*** free memory ***/
    free(IV);

    free(npts);
    free(b);      free(delta);
    free(depmin); free(depmax);
    free(tt);

    for (i=0; i<len; i++) {
        free(Rf[i]);
        free(bwp[i]);
    }
    free(Rf);
    free(bwp);

    free(Z_rf);
    free(Tbwp);

    return 0;
}




/* **************************************************** */
/* *****************  subroutines ********************* */
/* **************************************************** */

/*
 * get file line number
 */
int get_file_line(char *fname)
{
    int line_num=0;
    char c;
    FILE *fp;

    if ((fp = fopen(fname, "r")) == NULL) {
        fprintf(stderr, "Can not open file %s in read_file\n", fname);
        exit(-1);
    }

    for (c = getc(fp); c != EOF; c = getc(fp))
        if (c == '\n') line_num++;
    //fprintf(stdout, "line number: %d\n", line_num);

    fclose(fp);

    return line_num;
}


/*
 * read file into variables
 */
int read_file(char *fname, int line_num, IVAL *IV)
{
    int i;
    FILE *fp;

    if ((fp = fopen(fname, "r")) == NULL) {
        fprintf(stderr, "Can not open file %s in read_file\n", fname);
        exit(-1);
    }

    for (i=0; i<line_num; i++)
        fscanf(fp, "%s %f %f\n", IV[i].str, &IV[i].pRr, &IV[i].tt);

    fclose(fp);

    if (i != line_num) {
        fprintf(stderr, "Error in reading file %s\n", fname);
        exit(-1);
    }

    return i;
}


/*
 * read periods
 */
int read_pers(char *fname, int line_num, float *Tbwp)
{
    int i;
    FILE *fp;

    if ((fp = fopen(fname, "r")) == NULL) {
        fprintf(stderr, "Can not open file %s in read_file\n", fname);
        exit(-1);
    }

    for (i=0; i<line_num; i++)
        fscanf(fp, "%f\n", &Tbwp[i]);

    fclose(fp);

    if (i != line_num) {
        fprintf(stderr, "Error in reading file %s\n", fname);
        exit(-1);
    }

    return i;
}



/*
 * read SAC headers and data
 */
int read_sac_data(IVAL  *IV,
                  int   len,
                  float **data,
                  int   *npts,
                  float *b,
                  float *depmin,
                  float *depmax,
                  float *delta,
                  float *tt,
                  int   tmark)
{
    int i;
    SACHEAD hdr;

    for (i=0; i<len; i++) {
        data[i] = read_sac(IV[i].str, &hdr);

        npts[i]   = hdr.npts;
        b[i]      = hdr.b;
        depmin[i] = hdr.depmin;
        depmax[i] = hdr.depmax;
        delta[i]  = hdr.delta;
        tt[i]     = *((float *) &hdr + TMARK + tmark);
    }

    return 0;
}


/*
 * measure body-wave polarization
 */
int cal_bwp(IVAL  *IV,
            int   len,
            float **Rf,
            float *Z_rf,
            int   *npts,
            float *b,
            float *depmin,
            float *depmax,
            float *dt_Rf,
            int   nbwp,
            float *Tbwp,
            float **bwp)
{
    int i,j,k;
    int rf_num=0;

    int    Rf_b, Rf_e, nbwp_int;
    float  dt_bwp_int, tau, apparent_angle;
    double R_Rf_amp0, Z_Rf_amp0;

    for (k=0; k<len; k++) {
        if (isnan(depmin[k]) || isinf(depmax[k])) continue;
        if (dt_Rf[k] != dt_Rf[0]) continue; /* non-equale sampling rates */

        dt_bwp_int = dt_Rf[k];        /* RFHV integral dtau */
                                      /* RF first point index w.r.t. 0 */
                                      /* IV.tt is hd.b (negative)      */
        Rf_b = (int)((IV[k].tt - 0.00001) / dt_Rf[k]);
        Rf_e = Rf_b + npts[k] - 1;    /* RF last point index w.r.t. 0  */
        //fprintf(stdout, "%f %d %d %d %d\n", dt_bwp_int, Rf_b, Rf_e, npts[k], nbwp);

        /* apparent Vs at orignal gauss filter */
        IV[k].vs_app = sin(0.5 * atan(Rf[k][0-Rf_b] / Z_rf[0-Rf_b])) / IV[k].pRr;

        /* loop: periods */
        for (i = 0; i < nbwp; i++) {
            R_Rf_amp0 = 0.0;
            Z_Rf_amp0 = 0.0;
            nbwp_int  = (int)((Tbwp[i]+0.00001) / dt_bwp_int);
                                                        /* integral number */
            for (j = 0-nbwp_int; j <= nbwp_int; j++) {
                if ( (j >= Rf_b) && (j <= Rf_e) ) {
                    /* ignore dtau in the integral for amplitue ratio */
                    tau = M_PI * dt_bwp_int * j / (2.0 * Tbwp[i]);
                    R_Rf_amp0 = R_Rf_amp0 + Rf[k][j-Rf_b] * pow(cos(tau), 2);
                    Z_Rf_amp0 = Z_Rf_amp0 +  Z_rf[j-Rf_b] * pow(cos(tau), 2);
                }
            }

            apparent_angle = atan(R_Rf_amp0 / Z_Rf_amp0);
            bwp[k][i]      = sin(0.5*apparent_angle) / IV[k].pRr;
            /*fprintf(stdout, "%f %f %f %f %f\n", Tbwp[i], bwp[k][i],
                        rad2deg(apparent_angle), R_Rf_amp0, Z_Rf_amp0);*/
        }

        rf_num++;
    }

    return rf_num;
}
