/*
 Created by Sebastiano Vascon on 23/03/20.
 
 26, Alexandru Condurache 880890, Matteo Pagano 880833, Marco Pistollato 880889, Marco Quarta 880789
*/


#include <stdio.h>
#include "ip_lib.h"
#include "bmp.h"

ip_mat * ip_mat_create(unsigned int h, unsigned int w,unsigned  int k, float v){
    unsigned int i,j,z;
    /* inizializzazione ip_mat */
    ip_mat * mat;
    if (w > 0 && h > 0 && k > 0){
        mat = (ip_mat *) malloc (sizeof(ip_mat));
        if (mat == NULL){
            printf("allocazione non riuscita");
            exit(1);
        }
        mat -> w = w;
        mat -> h = h;
        mat -> k = k;
        mat -> stat = (stats *) malloc (sizeof(stats) * k);
        if (mat -> stat == NULL) {
            printf("allocazione non riuscita");
            exit(1);
        }
        mat -> data = (float ***) malloc (sizeof(float **) * h);
        if (mat -> data == NULL){
            printf("allocazione non riuscita");
            exit(1);
        }
        for (i = 0; i < h; i ++){
            mat -> data[i] = (float **) malloc (sizeof(float *) * w);
            if (mat -> data [i] == NULL){
                printf("allocazione non riuscita");
                exit(1);
            }
        }
        mat -> data [0][0] = (float *) malloc (sizeof(float) * w * h * k);
        if (mat -> data [0][0] == NULL){
            printf("allocazione non riuscita");
            exit(1);
        }
        for (i = k; i < w * h * k; i += k) {
            mat -> data [(i/k)/w][(i/k)%w] = mat -> data [0][0] + i;
        }
        /* riempimento matrice */
        for (i = 0; i < k; i ++){
            for (j = 0; j < h; j ++){
                for (z = 0; z < w; z ++){
                    mat -> data [j][z][i] = v;
                }
            }
        }
        return mat;
    } else {
        printf("ip_mat_create non riuscita: almeno una delle dimensioni h, w, k risulta essere nulla o negativa");
        exit(1);
    }
}

void ip_mat_free(ip_mat *a){
    if (a){
        unsigned int i;
        free (a -> data [0][0]);
        for (i = 0; i < a -> h; i ++){
            free (a -> data [i]);
        }
        free (a -> data);
        free (a -> stat);
        free(a);
    }
}

float get_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k){
    if(i<a->h && j<a->w &&k<a->k){  /* j>=0 and k>=0 and i>=0 is non sense*/
        return a->data[i][j][k];
    }else{
        printf("Errore get_val!!!");
        exit(1);
    }
}

void set_val(ip_mat * a, unsigned int i,unsigned int j,unsigned int k, float v){
    if(i<a->h && j<a->w &&k<a->k){
        a->data[i][j][k]=v;
    }else{
        printf("Errore set_val!!!");
        exit(1);
    }
}

void compute_stats(ip_mat * t) {
    if (t){
        unsigned int i, j, z;
        float n;
        float max, min, mean;
        for (i = 0; i < t->k; i++) {
            max = t->data[0][0][i];
            min = max;
            mean = 0;
            for (j = 0; j < t->h; j++) {
                for (z = 0; z < t->w; z++) {
                    n = t->data[j][z][i];
                    mean += n;
                    if (n > max) max = n;
                    if (n < min) min = n;
                }
            }
            mean /= t->w * t->h;
            t-> stat [i].min = min;
            t-> stat [i].max = max;
            t-> stat [i].mean = mean;
        }
    }
    else
        printf("compute_stats non riuscita: matrice passata come parametro vuota");
}

float gaussiana (float mean, float var, float x){
    float dev;
    dev = sqrt(var * 1.);
    return (1.f/(dev * sqrt(2. * PI))) * exp(-(((x - mean)*(x - mean))/(2.f * var)));
}

void ip_mat_init_random(ip_mat * t, float mean, float var){
    if (t){
        unsigned int i, j, z;
        float x;
        for (i = 0; i < t -> k; i ++){
            for (j = 0; j < t -> h; j ++){
                for (z = 0; z < t -> w; z ++){
                    x = get_normal_random(mean, sqrt(var));
                    t -> data [j][z][i] = gaussiana (mean, var, x);
                }
            }
        }
        compute_stats(t);
    }
    else
        printf("ip_mat_init_random non riuscita: matrice passata come parametro vuota");
}


ip_mat * ip_mat_copy(ip_mat * in){
    if (in){
        unsigned int i, j, z;
        ip_mat * out;
        out = ip_mat_create(in -> h, in -> w, in -> k, 0.f);
        for (i = 0; i < in -> k; i ++) {
            out -> stat[i].min = in -> stat[i].min;
            out -> stat[i].max = in -> stat[i].max;
            out -> stat[i].mean = in -> stat[i].mean;
        }
        for (i = 0; i < in -> k; i ++){
            for (j = 0; j < in -> h; j ++){
                for (z = 0; z < in -> w; z ++){
                    set_val (out, j, z, i, get_val(in, j, z, i));
                }
            }
        }
        return out;
    } else {
        printf ("ip_mat_copy non riuscita: matrice passata come parametro vuota");
        exit(1);
    }
}

ip_mat * ip_mat_subset(ip_mat * t, unsigned int row_start, unsigned int row_end, unsigned int col_start, unsigned int col_end){
    if (t && row_end - row_start > 0 && col_end - col_start > 0){
        unsigned int i, j, z;
        ip_mat * out;
        out = ip_mat_create (row_end - row_start + 1, col_end - col_start + 1, t -> k, 0.f);
        for (i = 0; i < t -> k; i ++){
            for (j = row_start; j < row_end; j ++){
                for (z = col_start; z < col_end; z ++){
                    set_val (out, j - row_start, z - col_start, i, get_val(t, j, z, i));
                }
            }
        }
        compute_stats (out);
        return out;
    } else {
        printf ("ip_mat_subset non riuscita: matrice passata come parametro vuota e/o indici errati");
        exit(1);
    }
}

ip_mat * ip_mat_concat(ip_mat * a, ip_mat * b, int dimensione){
    if (a && b){
        ip_mat * out;
        unsigned int i, j, z;
        out = NULL;
        switch(dimensione){
            case 0: /* concatenazione su h */
                if (a -> w == b -> w && a -> k == b -> k){
                    out = ip_mat_create(a -> h + b -> h, a -> w, a -> k, 0.f);
                    for (i = 0; i < a -> k; i ++){
                        for (j = 0; j < a -> h + b -> h; j ++){
                            for (z = 0; z < a -> w; z ++){
                                if (j < a -> h){
                                    set_val (out, j, z, i, get_val(a, j, z, i));
                                } else {
                                    set_val (out, j, z, i, get_val(b, j - a -> h, z, i));
                                }
                            }
                        }
                    }
                } else {
                    printf ("ip_mat_concat non riuscita: le due matrici hanno dimensioni diverse");
                    exit(1);
                }
                break;
                
            case 1: /* concatenazione su w */
                if (a -> h == b -> h && a -> k == b -> k){
                    out = ip_mat_create(a -> h, a -> w + b -> w, a -> k, 0.f);
                    for (i = 0; i < a -> k; i ++){
                        for (j = 0; j < a -> h; j ++){
                            for (z = 0; z < a -> w + b -> w; z ++){
                                if (z < a -> w){
                                    set_val (out, j, z, i, get_val(a, j, z, i));
                                } else {
                                    set_val (out, j, z, i, get_val(b, j , z - a -> w, i));
                                }
                            }
                        }
                    }
                } else {
                    printf ("ip_mat_concat non riuscita: le due matrici hanno dimensioni diverse");
                    exit(1);
                }
                break;
                
            case 2: /* concatenazione su k */
                if (a -> h == b -> h && a -> w == b -> w){
                    out = ip_mat_create(a -> h, a -> w, a -> k + b -> k, 0.f);
                    for (i = 0; i < a -> k + b -> k; i ++){
                        for (j = 0; j < a -> h; j ++){
                            for (z = 0; z < a -> w; z ++){
                                if (i < a -> k){
                                    set_val (out, j, z, i, get_val(a, j, z, i));
                                } else {
                                    set_val (out, j, z, i, get_val(b, j , z, i - a -> k));
                                }
                            }
                        }
                    }
                } else {
                    printf ("ip_mat_concat non riuscita: le due matrici hanno dimensioni diverse");
                    exit(1);
                }
                break;
                
            default :
                printf ("dimensione non compresa tra 0 e 2");
                exit(1);
        }
        compute_stats(out);
        return out;
    } else {
        printf ("ip_mat_concat non riuscita: almeno una delle matrici passate come parametro è vuota");
        exit(1);
    }
}

ip_mat * ip_mat_sum(ip_mat * a, ip_mat * b){
    if (a && b && a -> w == b -> w && a -> h == b -> h && a -> k == b -> k){
        unsigned int i, j, z;
        ip_mat * out;
        out = ip_mat_create (a->h, a->w, a->k, 0.f);
        for (i = 0; i < a -> k; i ++){
            for (j = 0; j < a -> h; j ++){
                for (z = 0; z < a -> w; z ++){
                    set_val (out, j, z, i, get_val(a, j, z, i) + get_val(b, j, z, i));
                }
            }
        }
        compute_stats(out);
        return out;
    } else {
        printf ("ip_mat_sum non riuscita: almeno una delle matrici passate come parametro è vuota o matrici differiscono per almeno una dimensione");
        exit (1);
    }
}

ip_mat * ip_mat_sub(ip_mat * a, ip_mat * b){
    if (a && b && a -> w == b -> w && a -> h == b -> h && a -> k == b -> k){
        unsigned int i, j, z;
        ip_mat * out;
        out = ip_mat_create (a->h, a->w, a->k, 0.f);
        for (i = 0; i < a -> k; i ++){
            for (j = 0; j < a -> h; j ++){
                for (z = 0; z < a -> w; z ++){
                    set_val (out, j, z, i, get_val(a, j, z, i) - get_val(b, j, z, i));
                }
            }
        }
        compute_stats(out);
        return out;
    } else {
        printf ("ip_mat_sub non riuscita: almeno una delle matrici passate come parametro è vuota o matrici differiscono per almeno una dimensione");
        exit (1);
    }
}

ip_mat * ip_mat_mul_scalar(ip_mat *a, float c){
    if (a){
        unsigned int i, j, z;
        ip_mat * out;
        out = ip_mat_create (a->h, a->w, a->k, 0.f);
        for (i = 0; i < a -> k; i ++){
            for (j = 0; j < a -> h; j ++){
                for (z = 0; z < a -> w; z ++){
                    set_val (out, j, z, i, get_val(a, j, z, i) * c);
                }
            }
        }
        compute_stats(out);
        return out;
    } else {
        printf ("ip_mat_mul_scalar non riuscita: matrice passata come parametro vuota");
        exit (1);
    }
}
    
ip_mat *  ip_mat_add_scalar(ip_mat *a, float c){
    if (a){
        unsigned int i, j, z;
        ip_mat * out;
        out = ip_mat_create (a->h, a->w, a->k, 0.f);
        for (i = 0; i < a -> k; i ++){
            for (j = 0; j < a -> h; j ++){
                for (z = 0; z < a -> w; z ++){
                    set_val (out, j, z, i, get_val(a, j, z, i) + c);
                }
            }
        }
        compute_stats(out);
        return out;
    } else {
        printf ("ip_mat_add_scalar non riuscita: matrice passata come parametro vuota");
        exit (1);
    }
}
   
ip_mat * ip_mat_mean(ip_mat * a, ip_mat * b){
    if (a && b && a -> w == b -> w && a -> h == b -> h && a -> k == b -> k){
        unsigned int i, j, z;
        ip_mat * out;
        out = ip_mat_create (a->h, a->w, a->k, 0.f);
        for (i = 0; i < a -> k; i ++){
            for (j = 0; j < a -> h; j ++){
                for (z = 0; z < a -> w; z ++){
                    set_val (out, j, z, i, (get_val(a, j, z, i) + get_val(b, j, z, i))/2);
                }
            }
        }
        compute_stats(out);
        return out;
    } else {
        printf ("ip_mat_mean non riuscita: almeno una delle matrici passate come parametro è vuota o matrici differiscono per almeno una dimensione");
        exit (1);
    }
}

ip_mat * ip_mat_to_gray_scale(ip_mat * in){
    /* assumiamo che ci siano 3 canali */
    if (in){
        unsigned int i, j, z;
        float temp; /* media dei valori di ogni canale */
        ip_mat * out;
        out = ip_mat_create (in -> h, in -> w, in -> k, 0.f);
        for (j = 0; j < in -> h; j ++){
            for (z = 0; z < in -> w; z ++){
                temp = 0;
                for(i = 0; i < 3; i ++){
                    temp += get_val(in, j, z, i);
                }
                temp /= 3;
                for(i = 0; i < 3; i ++){
                    set_val (out, j, z, i, temp);
                }
            }
        }
        compute_stats(out);
        return out;
    } else {
        printf ("ip_mat_to_gray_scale non riuscita: matrice passata come parametro vuota");
        exit(1);
    }
}

ip_mat * ip_mat_blend(ip_mat * a, ip_mat * b, float alpha){
    if (a && b && alpha >= 0.f && alpha <= 1.f){
        ip_mat * out;
        float blend;
        unsigned int h, w;
        unsigned int i, j, z;
        h = (a->h > b->h) ? b->h : a->h;
        w = (a->w > b->w) ? b->w : a->w;
        out = ip_mat_create (h, w, a -> k, 0.f);
        for (i = 0; i < out -> k; i ++){
            for (j = 0; j < h; j ++){
                for (z = 0; z < w; z ++){
                    blend = (alpha * get_val(a, j, z, i)) + ((1.0f - alpha) * get_val(b, j, z, i));
                    set_val (out, j, z, i, blend);
                }
            }
        }
        compute_stats(out);
        return out;
    } else {
        printf ("ip_mat_blend non riuscita: almeno una delle matrici passate come parametro è vuota o alpha non compreso tra 0 e 1");
        exit (1);
    }
}

ip_mat * ip_mat_brighten(ip_mat * a, float bright){
    ip_mat * out;
    out = ip_mat_add_scalar (a, bright);
    return out;
}

ip_mat * ip_mat_corrupt(ip_mat * a, float amount){
    if (a){
        ip_mat * out;
        unsigned int i, j, z;
        float gauss_noise; /* rumore applicato ad ogni pixel */
        out = ip_mat_create(a -> h, a -> w, a -> k, 0.f);
        for (i = 0; i < a -> k; i ++){
            for (j = 0; j < a -> h; j ++){
                for (z = 0; z < a -> w; z ++){
                    gauss_noise = get_normal_random(0, amount/2);
                    set_val (out, j, z, i, get_val(a, j, z, i) + gauss_noise);
                }
            }
        }
        compute_stats(out);
        return out;
    } else {
        printf ("ip_mat_to_corrupt non riuscita: matrice passata come parametro vuota");
        exit(1);
    }
}

ip_mat * ip_mat_padding(ip_mat * a, unsigned int pad_h, unsigned int pad_w){
    if (a){ /* assumiamo che pad_h e pad_w siano sempre positivi */
        unsigned int i, j, z;
        ip_mat * out;
        out = ip_mat_create(a->h + 2*pad_h, a->w + 2*pad_w, a->k, 0.f);
        for (i = 0; i < out->k; i ++){
            for (j = pad_h; j < out->h - pad_h; j ++){
                for (z = pad_w; z < out->w - pad_w; z ++){
                    set_val (out, j, z, i, get_val(a, j - pad_h, z - pad_w, i));
                }
            }
        } 
        compute_stats(out);
        return out;
    } else {
        printf ("ip_mat_padding non riuscita: matrice passata come parametro vuota");
        exit(1);
    }
}

float convpixel (ip_mat * sub_a, ip_mat * f, int k){
    unsigned int j, z;
    float out = 0.f;
    for (j = 0; j < f->h; j ++){
        for (z = 0; z < f->w; z ++){
            out += get_val(sub_a, j, z, k) * get_val(f, j, z, 0);
        }
    }
    return out;
}

ip_mat * ip_mat_convolve(ip_mat * a, ip_mat * f){
    if (a && f){
        ip_mat *out;
        ip_mat *conv; /* matrice ottenuta dal padding di a */
        ip_mat *sub; /* sottomatrice di conv */
        int vpad, opad;
        unsigned int i, j, z;
        out = ip_mat_create(a->h, a->w, a->k, 0.f);
        vpad = (f->h - 1)/2;
        opad = (f->w - 1)/2;
        conv = ip_mat_padding(a, vpad, opad);
        for (i = 0; i < a->k; i ++){
            for (j = 0; j < conv->h - f->h + 1; j ++){
                for (z = 0; z < conv->w - f->w + 1; z ++){
                    sub = ip_mat_subset(conv, j, j + f->h, z, z + f->w);
                    set_val (out, j, z, i, convpixel(sub, f, i));
                    ip_mat_free(sub);
                }
            }
        }
        compute_stats(out);
        ip_mat_free(conv);
        return out;
    } else {
        printf ("ip_mat_convolve non riuscita: almeno una delle matrici passate come parametro è vuota");
        exit(1);
    }
}

void riempiMat (ip_mat * a, float v []){
    unsigned int i, j, z;
    for (j = 0, i = 0; j < a->h; j ++){
        for (z = 0; z < a->w; z ++, i++){
            set_val (a, j, z, 0, v[i]);
        }
    }
}

ip_mat * create_sharpen_filter(){
    ip_mat * out;
    float valori [9] = {0.f, -1.f, 0.f, -1.f, 5.f, -1.f, 0.f, -1.f, 0.f};
    out = ip_mat_create(3, 3, 1, 0.f);
    riempiMat(out, valori);
    compute_stats(out);
    return out;
}

ip_mat * create_edge_filter(){
    ip_mat * out;
    float valori [9] = {-1.f, -1.f, -1.f, -1.f, 8.f, -1.f, -1.f, -1.f, -1.f};
    out = ip_mat_create(3, 3, 1, 0.f);
    riempiMat(out, valori);
    compute_stats(out);
    return out;
}

ip_mat * create_emboss_filter(){
    ip_mat * out;
    float valori [9] = {-2.f, -1.f, 0.f, -1.f, 1.f, 1.f, 0.f, 1.f, 2.f};
    out = ip_mat_create(3, 3, 1, 0.f);
    riempiMat(out, valori);
    compute_stats(out);
    return out;
}

ip_mat * create_average_filter(unsigned int w, unsigned int h, unsigned int k){
    unsigned int i, j, z;
    float c;
    ip_mat * out;
    out = ip_mat_create(h, w, k, 0.f);
    c = 1.f/(w * h);
    for (j = 0, i = 0; j < h; j ++){
        for (z = 0; z < w; z ++, i++){
            set_val (out, j, z, 0, c);
        }
    }
    compute_stats(out);
    return out;
}

float gaussiana2v (int x, int y, float sigma){
    return ((1.f/(2.f * PI * sigma * sigma)) * exp(-(x*x+y*y)/(2.f*sigma*sigma)));
}

ip_mat * create_gaussian_filter(unsigned int w, unsigned int h, unsigned int k, float sigma){
    int cx, cy, x, y;
    unsigned int z, i, j;
    float *somme, somma;
    ip_mat * gauss_filter;
    somme = (float*) malloc((sizeof(float))*k);
    if(somme==NULL){
        printf("Allocazione non riuscita");
        exit(1);
    }
    gauss_filter = ip_mat_create (h, w, k, 0.f);
    cx = h/2;
    cy = w/2;
    for (z = 0; z < k; z ++){
        somma = 0.f;
        for (i = 0; i < h; i ++){
            for (j = 0; j < w; j ++){
                x = i - cx;
                y = j - cy;
                set_val (gauss_filter, i, j, z, gaussiana2v(x, y, sigma));
                somma += get_val(gauss_filter, i, j, z);
            }
        }
        somme[z]=somma;
    }
    for (z = 0; z < k; z ++){
        for (i = 0; i < h; i ++){
            for (j = 0; j < w; j ++){
                set_val (gauss_filter, i, j, z, get_val(gauss_filter, i, j, z)/somme[z]);
            }
        }
    }
    free(somme);
    compute_stats(gauss_filter);
    return gauss_filter;
}

void rescale(ip_mat * t, float new_max){
    if(t){
        unsigned int i, j, z;
        float val;
        compute_stats(t);
        for (z = 0; z < t->k; z ++){
            for (i = 0; i < t->h; i ++){
                for (j = 0; j < t->w; j ++){
                    val =  (get_val(t, i, j, z) - t->stat[z].min)/(t->stat[z].max - t->stat[z].min);
                    set_val(t, i, j, z, val * new_max);
                }
            }
        }
        compute_stats(t);
    } else {
        printf("rescale non riuscita: matrice passata come parametro vuota");
        exit(1);
    }
}

void clamp(ip_mat * t, float low, float high){
    if (t){
        float v;
        unsigned int i, j, z;
        for (z = 0; z < t->k; z ++){
            for (i = 0; i < t->h; i ++){
                for (j = 0; j < t->w; j ++){
                    v = get_val(t, i, j, z);
                    if (v < low) set_val(t, i, j, z, low);
                    if (v > high) set_val(t, i, j, z, high);
                }
            }
        }
        compute_stats(t);
    } else {
        printf("clamp non riuscita: matrice passata come parametro vuota");
        exit(1);
    }
}


void ip_mat_show(ip_mat * t){
    unsigned int i,l,j;
    printf("Matrix of size %d x %d x %d (hxwxk)\n",t->h,t->w,t->k);
    for (l = 0; l < t->k; l++) {
        printf("Slice %d\n", l);
        for(i=0;i<t->h;i++) {
            for (j = 0; j < t->w; j++) {
                printf("%f ", get_val(t,i,j,l));
            }
            printf("\n");
        }
        printf("\n");
    }
}

void ip_mat_show_stats(ip_mat * t){
    unsigned int k;

    compute_stats(t);

    for(k=0;k<t->k;k++){
        printf("Channel %d:\n", k);
        printf("\t Min: %f\n", t->stat[k].min);
        printf("\t Max: %f\n", t->stat[k].max);
        printf("\t Mean: %f\n", t->stat[k].mean);
    }
}

ip_mat * bitmap_to_ip_mat(Bitmap * img){
    unsigned int i=0,j=0;

    unsigned char R,G,B;

    unsigned int h = img->h;
    unsigned int w = img->w;

    ip_mat * out = ip_mat_create(h, w,3,0);

    for (i = 0; i < h; i++)              /* rows */
    {
        for (j = 0; j < w; j++)          /* columns */
        {
            bm_get_pixel(img, j,i,&R, &G, &B);
            set_val(out,i,j,0,(float) R);
            set_val(out,i,j,1,(float) G);
            set_val(out,i,j,2,(float) B);
        }
    }
    compute_stats(out);
    return out;
}

Bitmap * ip_mat_to_bitmap(ip_mat * t){

    Bitmap *b = bm_create(t->w,t->h);

    unsigned int i, j;
    for (i = 0; i < t->h; i++)              /* rows */
    {
        for (j = 0; j < t->w; j++)          /* columns */
        {
            bm_set_pixel(b, j,i, (unsigned char) get_val(t,i,j,0),
                    (unsigned char) get_val(t,i,j,1),
                    (unsigned char) get_val(t,i,j,2));
        }
    }
    return b;
}

float get_normal_random(float media, float std){

    float y1 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float y2 = ( (float)(rand()) + 1. )/( (float)(RAND_MAX) + 1. );
    float num = cos(2*PI*y2)*sqrt(-2.*log(y1));

    return media + num*std;
}
