// Copyright, Matthieu Guerquin-Kern, 2012

#include <cmath>

#ifndef FAST
//#define FAST // Uncomment this if you want the simpler/more-reliable/slower code
#endif

/* Min and Max */
#ifndef MAX
#define	MAX(A, B) ((A) > (B) ? (A) : (B))
#endif

// Global variables
#ifndef M_PI
#define M_PI 3.1415926535897931159979634685441851615906
#endif

/************************* H ***************************************/
#ifndef FAST
static void H_complex(const double X_real[], const double X_imag[], const double kx[], const double ky[], double *m_real, double *m_imag, const unsigned int m, const unsigned int n, const unsigned int Nb_samples)
{
    int ind_pos,ind_sample,posx,posy;
    double scalprod,c,s;
    for (ind_sample=0;ind_sample < Nb_samples; ind_sample++){
        for (posy=0,ind_pos=0;posy < n; posy++){
            for (posx=0;posx < m; posx++,ind_pos++){
                scalprod = posx*kx[ind_sample]/m+posy*ky[ind_sample]/n;
                c = cos(2*M_PI*scalprod);
                s = sin(2*M_PI*scalprod);
                m_real[ind_sample] += X_real[ind_pos]*c - X_imag[ind_pos]*s;
                m_imag[ind_sample] += X_real[ind_pos]*s + X_imag[ind_pos]*c;
            }
        }
    }
    return;
}
#else // FAST
static void H_complex( const double X_real[], const double X_imag[], const double kx[], const double ky[], double *m_real, double *m_imag, const unsigned int m, const unsigned int n, const unsigned int Nb_samples)
{
    int ind_pos,ind_sample,posx,posy;
    double rx,ry,ix,iy,cy,cyn,sy,c,s;
    double *cx=(double*)malloc(m*sizeof(double));
    double *sx=(double*)malloc(m*sizeof(double));
    for (ind_sample=0;ind_sample < Nb_samples; ind_sample++){
        rx = cos(2*M_PI*kx[ind_sample]/m);
        ry = cos(2*M_PI*ky[ind_sample]/n);
        ix = sin(2*M_PI*kx[ind_sample]/m);
        iy = sin(2*M_PI*ky[ind_sample]/n);
        // posy = 0
        for (posx=0,cx[0]=1,sx[0]=0,ind_pos=0;posx < m-1; posx++,ind_pos++){
            m_real[ind_sample] += X_real[ind_pos]*cx[posx] - X_imag[ind_pos]*sx[posx];//c = cos(2*M_PI*scalprod);
            m_imag[ind_sample] += X_real[ind_pos]*sx[posx] + X_imag[ind_pos]*cx[posx];//s = sin(2*M_PI*scalprod);
            cx[posx+1] = cx[posx]*rx-sx[posx]*ix;
            sx[posx+1] = cx[posx]*ix+sx[posx]*rx;
        }
        m_real[ind_sample] += X_real[ind_pos]*cx[posx+1] - X_imag[ind_pos]*sx[posx+1];
        m_imag[ind_sample] += X_real[ind_pos]*sx[posx+1] + X_imag[ind_pos]*cx[posx+1];
        ind_pos++;
        // posy>0
        for (posy=1,cy=ry,sy=iy;posy < n; posy++){
            for (posx=0;posx < m; posx++,ind_pos++){
                //scalprod = posx*kx[ind_sample]/m + posy*ky[ind_sample]/n;
                c = cx[posx]*cy-sx[posx]*sy;//c = cos(2*M_PI*scalprod);
                s = sx[posx]*cy+cx[posx]*sy;//s = sin(2*M_PI*scalprod);
                m_real[ind_sample] += X_real[ind_pos]*c - X_imag[ind_pos]*s;
                m_imag[ind_sample] += X_real[ind_pos]*s + X_imag[ind_pos]*c;
            }
            cyn = cy*ry-sy*iy;
            sy = cy*iy+sy*ry;
            cy = cyn;
        }
    }
    free(cx);
    free(sx);
    return;
}
#endif //FAST

#ifndef FAST
static void H_real(const double X_real[], const double kx[], const double ky[], double *m_real, double *m_imag, const unsigned int m, const unsigned int n, const unsigned int Nb_samples)
{
    int ind_pos,ind_sample,posx,posy;
    double scalprod,c,s;
    for (ind_sample=0;ind_sample < Nb_samples; ind_sample++){
        for (posy=0,ind_pos=0;posy < n; posy++){
            for (posx=0;posx < m; posx++,ind_pos++){
                scalprod = posx*kx[ind_sample]/m+posy*ky[ind_sample]/n;
                c = cos(2*M_PI*scalprod);
                s = sin(2*M_PI*scalprod);
                m_real[ind_sample] += X_real[ind_pos]*c;
                m_imag[ind_sample] += X_real[ind_pos]*s;
            }
        }
    }
    return;
}
#else // FAST
static void H_real(const double X_real[], const double kx[], const double ky[], double *m_real, double *m_imag, const unsigned int m, const unsigned int n, const unsigned int Nb_samples)
{
    int ind_pos,ind_sample,posx,posy;
    double rx,ry,ix,iy,cy,cyn,sy;
    double *cx=(double*)malloc(m*sizeof(double));
    double *sx=(double*)malloc(m*sizeof(double));
    
    for (ind_sample=0;ind_sample < Nb_samples; ind_sample++){
        rx = cos(2*M_PI*kx[ind_sample]/m);
        ry = cos(2*M_PI*ky[ind_sample]/n);
        ix = sin(2*M_PI*kx[ind_sample]/m);
        iy = sin(2*M_PI*ky[ind_sample]/n);
        // posy = 0 and precomputing cx and cs
        for (posx=0,cx[0]=1,sx[0]=0,ind_pos=0;posx < m-1; posx++,ind_pos++){
            m_real[ind_sample] += X_real[ind_pos]*cx[posx];
            m_imag[ind_sample] += X_real[ind_pos]*sx[posx];
            cx[posx+1] = cx[posx]*rx-sx[posx]*ix;
            sx[posx+1] = cx[posx]*ix+sx[posx]*rx;
        }
        m_real[ind_sample] += X_real[ind_pos]*cx[m-1];
        m_imag[ind_sample] += X_real[ind_pos]*sx[m-1];
        ind_pos++;
        // posy>0
        for (posy=1,cy=ry,sy=iy;posy < n; posy++){
            for (posx=0;posx < m; posx++,ind_pos++){
                //scalprod = posx*kx[ind_sample]/m + posy*ky[ind_sample]/n;
                m_real[ind_sample] += X_real[ind_pos]*(cx[posx]*cy-sx[posx]*sy);//c = cos(2*M_PI*scalprod);
                m_imag[ind_sample] += X_real[ind_pos]*(sx[posx]*cy+cx[posx]*sy);//s = sin(2*M_PI*scalprod);
            }
            cyn = cy*ry-sy*iy;
            sy = cy*iy+sy*ry;
            cy = cyn;
        }
    }
    free(cx);free(sx);
    return;
}
#endif //FAST

/************************* HH ***************************************/
#ifndef FAST
static void HH( const double m_real[], const double m_imag[], const double kx[], const double ky[], double *X_real, double *X_imag, const unsigned int m, const unsigned int n, const unsigned int Nb_samples)
{
    int ind_pos,ind_sample,posx,posy;
    double c,s,scalprod;
    for (ind_sample=0;ind_sample < Nb_samples; ind_sample++){
        for (posy=0,ind_pos=0;posy < n; posy++){
            for (posx=0;posx < m; posx++,ind_pos++){
                scalprod = -posx*kx[ind_sample]/m-posy*ky[ind_sample]/n;
                c = cos(2*M_PI*scalprod);
                s = sin(2*M_PI*scalprod);
                X_real[ind_pos] += m_real[ind_sample]*c - m_imag[ind_sample]*s;
                X_imag[ind_pos] += m_real[ind_sample]*s + m_imag[ind_sample]*c;
            }
        }
    }
    return;
}
#else //FAST
static void HH( const double m_real[], const double m_imag[], const double kx[], const double ky[], double *X_real, double *X_imag, const unsigned int m, const unsigned int n, const unsigned int Nb_samples)
{
    int ind_pos,ind_sample,posx,posy;
    double rx,ry,ix,iy,cy,cyn,sy;
    double *cx=(double*)malloc(m*sizeof(double));
    double *sx=(double*)malloc(m*sizeof(double));
    double c,s,scalprod;
    for (ind_sample=0;ind_sample < Nb_samples; ind_sample++){
        rx = cos(-2*M_PI*kx[ind_sample]/m);
        ry = cos(-2*M_PI*ky[ind_sample]/n);
        ix = sin(-2*M_PI*kx[ind_sample]/m);
        iy = sin(-2*M_PI*ky[ind_sample]/n);
        // posy = 0 and precomputing cx and cs
        for (posx=0,ind_pos=0,cx[0]=1,sx[0]=0;posx < m-1; posx++,ind_pos++){
            X_real[ind_pos] += m_real[ind_sample]*cx[posx] - m_imag[ind_sample]*sx[posx];
            X_imag[ind_pos] += m_real[ind_sample]*sx[posx] + m_imag[ind_sample]*cx[posx];
            cx[posx+1] = cx[posx]*rx-sx[posx]*ix;
            sx[posx+1] = cx[posx]*ix+sx[posx]*rx;
        }
        X_real[ind_pos] += m_real[ind_sample]*cx[m-1] - m_imag[ind_sample]*sx[m-1];
        X_imag[ind_pos] += m_real[ind_sample]*sx[m-1] + m_imag[ind_sample]*cx[m-1];
        ind_pos++;
        // posy>0
        for (posy=1,cy=ry,sy=iy;posy < n; posy++){
            for (posx=0;posx < m; posx++,ind_pos++){
                //scalprod = -posx*kx[ind_sample]/m-posy*ky[ind_sample]/n;
                c = cx[posx]*cy-sx[posx]*sy;//cos(2*M_PI*scalprod);//
                s = sx[posx]*cy+cx[posx]*sy;//sin(2*M_PI*scalprod);//
                X_real[ind_pos] += m_real[ind_sample]*c - m_imag[ind_sample]*s;
                X_imag[ind_pos] += m_real[ind_sample]*s + m_imag[ind_sample]*c;
            }
            cyn = cy*ry-sy*iy;
            sy = cy*iy+sy*ry;
            cy = cyn;
        }
    }
    return;
}
#endif //FAST

/************************* G ***************************************/
#ifndef FAST
static void G( const double kx[], const double ky[], double *G_real, double *G_imag, const double weight[], const unsigned int m, const unsigned int n, const unsigned int Nb_samples)
{
    int ind_pos,ind_pos_opp,ind_sample,posx,posy,endx;
    double scalprod;
    
    for (posy=-n+1,ind_pos_opp=(2*m-1)*(2*n-1)-1,ind_pos=0;posy<1; posy++){
        if (posy == 0){ endx=1;}else{ endx = m;} // we scan half of the space, exploiting Hermitian symmetry to fill the rest
        for (posx=-m+1;posx < endx; posx++,ind_pos++,ind_pos_opp--){
            for (ind_sample=0;ind_sample < Nb_samples; ind_sample++){
                scalprod = posx*kx[ind_sample]/m+posy*ky[ind_sample]/n;
                G_real[ind_pos] += cos(-2*M_PI*scalprod)*weight[ind_sample];
                G_imag[ind_pos] += sin(-2*M_PI*scalprod)*weight[ind_sample];
            }
            // Hermitian symmetry
            G_real[ind_pos_opp] = G_real[ind_pos];
            G_imag[ind_pos_opp] = -G_imag[ind_pos];
        }
    }
    return;
}
#else // FAST
static void G( const double kx[], const double ky[], double *G_real, double *G_imag, const double weight[], const unsigned int m, const unsigned int n, const unsigned int Nb_samples)
{
    int ind_pos,ind_pos_opp,ind_sample,posx,posy,endx;
    double rx,ry,ix,iy,cy,cyn,sy,c,s;
    double *cx=(double*)malloc(m*sizeof(double));
    double *sx=(double*)malloc(m*sizeof(double));
    double scalprod;
    for (ind_sample=0;ind_sample < Nb_samples; ind_sample++){
        rx = cos(-2*M_PI*kx[ind_sample]/m);
        ry = cos(-2*M_PI*ky[ind_sample]/n);
        ix = sin(-2*M_PI*kx[ind_sample]/m);
        iy = sin(-2*M_PI*ky[ind_sample]/n);
        // Case posy==0
        for (posx=0,ind_pos_opp=2*m*n-m-n,ind_pos=2*m*n-m-n,cx[0]=weight[ind_sample],sx[0]=0;posx < m-1; posx++,ind_pos--){
            //scalprod = -posx*kx[ind_sample]/m;
            G_real[ind_pos] += cx[posx];//weight[ind_sample]*cos(-2*M_PI*scalprod);//
            G_imag[ind_pos] += sx[posx];//weight[ind_sample]*sin(-2*M_PI*scalprod);//
            cx[posx+1] = cx[posx]*rx-sx[posx]*ix;
            sx[posx+1] = cx[posx]*ix+sx[posx]*rx;
        }
        G_real[ind_pos] += cx[posx+1];
        G_imag[ind_pos] += cx[posx+1];
        ind_pos--;
        // posy<0
        for (posy=-1,cy=ry,sy=iy;posy>-n; posy--){
            for (posx=m-1;posx > 0; posx--,ind_pos--){
                //scalprod = posx*kx[ind_sample]/m+posy*ky[ind_sample]/n;
                G_real[ind_pos] += (cx[posx]*cy-sx[posx]*sy);//weight[ind_sample]*cos(-2*M_PI*scalprod);//
                G_imag[ind_pos] += (sx[posx]*cy+cx[posx]*sy);//weight[ind_sample]*sin(-2*M_PI*scalprod);//
            }
            for (posx=0;posx < m; posx++,ind_pos--){
                //scalprod = -posx*kx[ind_sample]/m+posy*ky[ind_sample]/n;
                G_real[ind_pos] += (cx[posx]*cy-sx[posx]*sy);//weight[ind_sample]*cos(-2*M_PI*scalprod);//
                G_imag[ind_pos] += (sx[posx]*cy+cx[posx]*sy);//weight[ind_sample]*sin(-2*M_PI*scalprod);//
            }
            cyn = cy*ry-sy*iy;
            sy = cy*iy+sy*ry;
            cy = cyn;
        }
    }
    // Hermitian symmetry
    for (ind_pos=0,ind_pos_opp=(2*m-1)*(2*n-1)-1;ind_pos<2*m*n-m-n+1;ind_pos++,ind_pos_opp--){
        G_real[ind_pos_opp] = G_real[ind_pos];
        G_imag[ind_pos_opp] = -G_imag[ind_pos];
    }
    return;
}
#endif // FAST

double* create_periodically_padded_array(const double *a,const int mrx,const int mry,const int msp){
    // SHOULDN'T BE CALLED IF msp>mrx OR msp>mry
    if (msp>mrx||msp>mry){
        return NULL;
    }
    double *pad=(double*)malloc((mrx+2*msp)*(mry+2*msp)*sizeof(double));
    double *a_p, *pad_p;
    // Copy array in center
    for (int j=0;j<mry;j++){
        a_p = (double*)a + mrx*j;
        pad_p = pad + (mrx+2*msp)*(j+msp) + msp;
        for (int i=0;i<mrx;i++){ // Is memcpy faster?
            *(pad_p+i) = *(a_p+i);
        }
    }
    // Central borders in y
    for (int j=0;j<msp;j++){
        // Top
        a_p = (double*)a + mrx*(mry-msp+j);
        pad_p = pad + (mrx+2*msp)*j + msp;
        for (int i=0;i<mrx;i++){
            *(pad_p+i) = *(a_p+i);
        }
        // Bottom
        a_p = (double*)a + mrx*j;
        pad_p = pad + (mrx+2*msp)*(mry+msp+j) + msp;
        for (int i=0;i<mrx;i++){
            *(pad_p+i) = *(a_p+i);
        }
    }
    // Central borders in x
    for (int j=0;j<mry;j++){
        // Left
        a_p = (double*)a + mrx*j+mrx-msp;
        pad_p = pad + (mrx+2*msp)*(j+msp);
        for (int i=0;i<msp;i++){
            *(pad_p+i) = *(a_p+i);
        }
        // Right
        a_p = (double*)a + mrx*j;
        pad_p = pad + (mrx+2*msp)*(j+msp)+mrx+msp;
        for (int i=0;i<msp;i++){
            *(pad_p+i) = *(a_p+i);
        }
    }
    // Corners
    for (int j=0;j<msp;j++){
        // Top - Left
        a_p = (double*)a + mrx*(mry-msp+j)+mrx-msp;
        pad_p = pad + (mrx+2*msp)*j;
        for (int i=0;i<msp;i++){
            *(pad_p+i) = *(a_p+i);
        }
        // Bottom - Right
        a_p = (double*)a + mrx*j;
        pad_p = pad + (mrx+2*msp)*(mry+msp+j)+mrx+msp;
        for (int i=0;i<msp;i++){
            *(pad_p+i) = *(a_p+i);
        }
        // Bottom - Left
        a_p = (double*)a + mrx*(mry-msp+j);
        pad_p = pad + (mrx+2*msp)*j+mrx+msp;
        for (int i=0;i<msp;i++){
            *(pad_p+i) = *(a_p+i);
        }
        // Top - Right
        a_p = (double*)a + mrx*j+mrx-msp;
        pad_p = pad + (mrx+2*msp)*(mry+msp+j);
        for (int i=0;i<msp;i++){
            *(pad_p+i) = *(a_p+i);
        }
    }
    return pad;
}

/****************** forward_gaussian_gridding ************************/
/*static void forward_gaussian_gridding(
        double *v0r,
        double *v0i,
        const double *yr,
        const double *yi,
        const int msp,
        const int *px,
        const int *py,
        const int mrx,
        const int mry,
        const int nb,
        const double *e2x,
        const double *e2y,
        const double *e3x,
        const double *e3y)
{
    double *v0r_p=v0r,*v0i_p=v0i;
    double *yrP_p,*yrM_p,*yiP_p,*yiM_p;
    int *px_p=(int*)px,*py_p=(int*)py;
    int colindM,colindP;
    int rowindM,rowindP;
    double *e2x_p=(double*)e2x,*e2y_p=(double*)e2y;
    double *e2xPl=(double*)malloc(msp*sizeof(double));
    double *e2xMl=(double*)malloc(msp*sizeof(double));
    double e2xPinc,exP; // in 'xPl', 'xMl', 'yPl', 'yMl',
    double e2xMinc,exM; // x means 'first dimension', y is the second
    double e2yPl,e2yPinc,eyP; // P means 'plus' and M means 'minus'
    double e2yMl,e2yMinc,eyM; // l is the variable varying from -Msp+1 to Msp
    
    double *yPADr = create_periodically_padded_array(yr,mrx,mry,msp);
    double *yPADi = create_periodically_padded_array(yi,mrx,mry,msp);
    
    for(int i=0;i<nb;i++,v0r_p++,v0i_p++,px_p++,py_p++,e2x_p++,e2y_p++){
        e2xPinc = *e2x_p;
        e2xMinc = 1.0/(*e2x_p);
        e2yPinc = *e2y_p;
        e2yMinc = 1.0/(*e2y_p);
        e2yPl = e2yPinc;
        e2yMl = 1.0;
        e2xPl[0] = e2xPinc;
        e2xMl[0] = 1.0;
        for(int lx=1;lx<msp;lx++){
            e2xPl[lx]=e2xPl[lx-1]*e2xPinc;
            e2xMl[lx]=e2xMl[lx-1]*e2xMinc;
        }
        for(int ly=0;ly<msp;ly++,e2yPl*=e2yPinc,e2yMl*=e2yMinc){
            eyP = e2yPl*e3y[ly+1];
            eyM = e2yMl*e3y[ly];
            colindP = *py_p + ly + 1+msp;
            //if (colindP>=mry){colindP -= mry;}
            colindM = *py_p - ly+msp;
            //if (colindM<0){colindM += mry;}
            yrP_p = (double*)yPADr+colindP*(mrx+2*msp);
            yrM_p = (double*)yPADr+colindM*(mrx+2*msp);
            yiP_p = (double*)yPADi+colindP*(mrx+2*msp);
            yiM_p = (double*)yPADi+colindM*(mrx+2*msp);
            for(int lx=0;lx<msp;lx++){
                exP = e2xPl[lx]*e3x[lx+1];
                exM = e2xMl[lx]*e3x[lx];
                rowindP = *px_p + lx + 1+msp;
                //if (rowindP>=mrx){rowindP -= mrx;} // Looks like these periodic boundary conditions
                rowindM = *px_p - lx+msp;              // take half of the computing time
                //if (rowindM<0){rowindM += mrx;}
                *v0r_p += (*(yrP_p+rowindP)*eyP+*(yrM_p+rowindP)*eyM)*exP;
                *v0i_p += (*(yiP_p+rowindP)*eyP+*(yiM_p+rowindP)*eyM)*exP;
                *v0r_p += (*(yrP_p+rowindM)*eyP+*(yrM_p+rowindM)*eyM)*exM;
                *v0i_p += (*(yiP_p+rowindM)*eyP+*(yiM_p+rowindM)*eyM)*exM;
            }
        }
    }
    free(e2xPl);
    free(e2xMl);
    free(yPADr);
    free(yPADi);
    return;
}*/

static void forward_gaussian_gridding(
        double *v0r,
        double *v0i,
        const double *yr,
        const double *yi,
        const int msp,
        const int *px,
        const int *py,
        const int mrx,
        const int mry,
        const int nb,
        const double *e2x,
        const double *e2y,
        const double *e3x,
        const double *e3y)
{
    double *v0r_p=v0r,*v0i_p=v0i;
    double *yrP_p,*yrM_p,*yiP_p,*yiM_p;
    int *px_p=(int*)px,*py_p=(int*)py;
    int colindM,colindP;
    int rowindM,rowindP;
    double *e2x_p=(double*)e2x,*e2y_p=(double*)e2y;
    double *e2xPl=(double*)malloc(msp*sizeof(double));
    double *e2xMl=(double*)malloc(msp*sizeof(double));
    double e2xPinc,exP; // in 'xPl', 'xMl', 'yPl', 'yMl',
    double e2xMinc,exM; // x means 'first dimension', y is the second
    double e2yPl,e2yPinc,eyP; // P means 'plus' and M means 'minus'
    double e2yMl,e2yMinc,eyM; // l is the variable varying from -Msp+1 to Msp
    
    for(int i=0;i<nb;i++,v0r_p++,v0i_p++,px_p++,py_p++,e2x_p++,e2y_p++){
        e2xPinc = *e2x_p;
        e2xMinc = 1.0/(*e2x_p);
        e2yPinc = *e2y_p;
        e2yMinc = 1.0/(*e2y_p);
        e2yPl = e2yPinc;
        e2yMl = 1.0;
        e2xPl[0] = e2xPinc;
        e2xMl[0] = 1.0;
        for(int lx=1;lx<msp;lx++){
            e2xPl[lx]=e2xPl[lx-1]*e2xPinc;
            e2xMl[lx]=e2xMl[lx-1]*e2xMinc;
        }
        for(int ly=0;ly<msp;ly++,e2yPl*=e2yPinc,e2yMl*=e2yMinc){
            eyP = e2yPl*e3y[ly+1];
            eyM = e2yMl*e3y[ly];
            colindP = *py_p + ly + 1;
            if (colindP>=mry){colindP -= mry;}
            colindM = *py_p - ly;
            if (colindM<0){colindM += mry;}
            yrP_p = (double*)yr+colindP*mrx;
            yrM_p = (double*)yr+colindM*mrx;
            yiP_p = (double*)yi+colindP*mrx;
            yiM_p = (double*)yi+colindM*mrx;
            for(int lx=0;lx<msp;lx++){
                exP = e2xPl[lx]*e3x[lx+1];
                exM = e2xMl[lx]*e3x[lx];
                rowindP = *px_p + lx + 1;
                if (rowindP>=mrx){rowindP -= mrx;} // Looks like these periodic boundary conditions
                rowindM = *px_p - lx;              // take half of the computing time
                if (rowindM<0){rowindM += mrx;}
                *v0r_p += (*(yrP_p+rowindP)*eyP+*(yrM_p+rowindP)*eyM)*exP;
                *v0i_p += (*(yiP_p+rowindP)*eyP+*(yiM_p+rowindP)*eyM)*exP;
                *v0r_p += (*(yrP_p+rowindM)*eyP+*(yrM_p+rowindM)*eyM)*exM;
                *v0i_p += (*(yiP_p+rowindM)*eyP+*(yiM_p+rowindM)*eyM)*exM;
            }
        }
    }
    free(e2xPl);
    free(e2xMl);
    return;
}

static void forward_gaussian_gridding_real(
        double *v0r,
        const double *yr,
        const int msp,
        const int *px,
        const int *py,
        const int mrx,
        const int mry,
        const int nb,
        const double *e2x,
        const double *e2y,
        const double *e3x,
        const double *e3y)
{
    double *v0r_p=v0r;
    double *yrP_p,*yrM_p;
    int *px_p=(int*)px,*py_p=(int*)py;
    int colindM,colindP;
    int rowindM,rowindP;
    double *e2x_p=(double*)e2x,*e2y_p=(double*)e2y;
    double *e2xPl=(double*)malloc(msp*sizeof(double));
    double *e2xMl=(double*)malloc(msp*sizeof(double));
    double e2xPinc,exP; // in 'xPl', 'xMl', 'yPl', 'yMl',
    double e2xMinc,exM; // x means 'first dimension', y is the second
    double e2yPl,e2yPinc,eyP; // P means 'plus' and M means 'minus'
    double e2yMl,e2yMinc,eyM; // l is the variable varying from -Msp+1 to Msp
    
    for(int i=0;i<nb;i++,v0r_p++,px_p++,py_p++,e2x_p++,e2y_p++){
        e2xPinc = *e2x_p;
        e2xMinc = 1.0/(*e2x_p);
        e2yPinc = *e2y_p;
        e2yMinc = 1.0/(*e2y_p);
        e2yPl = e2yPinc;
        e2yMl = 1.0;
        e2xPl[0] = e2xPinc;
        e2xMl[0] = 1.0;
        for(int lx=1;lx<msp;lx++){
            e2xPl[lx]=e2xPl[lx-1]*e2xPinc;
            e2xMl[lx]=e2xMl[lx-1]*e2xMinc;
        }
        for(int ly=0;ly<msp;ly++,e2yPl*=e2yPinc,e2yMl*=e2yMinc){
            eyP = e2yPl*e3y[ly+1];
            eyM = e2yMl*e3y[ly];
            colindP = *py_p + ly + 1;
            if (colindP>=mry){colindP -= mry;}
            colindM = *py_p - ly;
            if (colindM<0){colindM += mry;}
            yrP_p = (double*)yr+colindP*mrx;
            yrM_p = (double*)yr+colindM*mrx;
            for(int lx=0;lx<msp;lx++){
                exP = e2xPl[lx]*e3x[lx+1];
                exM = e2xMl[lx]*e3x[lx];
                rowindP = *px_p + lx + 1;
                if (rowindP>=mrx){rowindP -= mrx;} // Looks like these periodic boundary conditions
                rowindM = *px_p - lx;              // take half of the computing time
                if (rowindM<0){rowindM += mrx;}
                *v0r_p += (*(yrP_p+rowindP)*eyP+*(yrM_p+rowindP)*eyM)*exP;
                *v0r_p += (*(yrP_p+rowindM)*eyP+*(yrM_p+rowindM)*eyM)*exM;
            }
        }
    }
    free(e2xPl);
    free(e2xMl);
    return;
}

/****************** backward_gaussian_gridding ************************/
static void backward_gaussian_gridding(
        double *cr,
        double *ci,
        const double *v0r,
        const double *v0i,
        const int msp,
        const int *px,
        const int *py,
        const int mrx,
        const int mry,
        const int nb,
        const double *e2x,
        const double *e2y,
        const double *e3x,
        const double *e3y)
{
    double *crP_p,*crM_p,*ciP_p,*ciM_p;
    double *v0r_p=(double*)v0r,*v0i_p=(double*)v0i;
    int *px_p=(int*)px,*py_p=(int*)py;
    int colindM,colindP;
    int rowindM,rowindP;
    double *e2x_p=(double*)e2x,*e2y_p=(double*)e2y;
    double *e2xPl=(double*)malloc(msp*sizeof(double));
    double *e2xMl=(double*)malloc(msp*sizeof(double));
    double e2xPinc,exP; // in 'xPl', 'xMl', 'yPl', 'yMl',
    double e2xMinc,exM; // x means 'first dimension', y is the second
    double e2yPl,e2yPinc,eyP; // P means 'plus' and M means 'minus'
    double e2yMl,e2yMinc,eyM; // l is the variable varying from -Msp+1 to Msp
    
    for(int i=0;i<nb;i++,v0r_p++,v0i_p++,px_p++,py_p++,e2x_p++,e2y_p++){
        e2xPinc = *e2x_p;
        e2xMinc = 1.0/(*e2x_p);
        e2yPinc = *e2y_p;
        e2yMinc = 1.0/(*e2y_p);
        e2yPl = e2yPinc;
        e2yMl = 1.0;
        e2xPl[0] = e2xPinc;
        e2xMl[0] = 1.0;
        for(int lx=1;lx<msp;lx++){
            e2xPl[lx]=e2xPl[lx-1]*e2xPinc;
            e2xMl[lx]=e2xMl[lx-1]*e2xMinc;
        }
        for(int ly=0;ly<msp;ly++,e2yPl*=e2yPinc,e2yMl*=e2yMinc){
            eyP = e2yPl*e3y[ly+1];
            eyM = e2yMl*e3y[ly];
            colindP = *py_p + ly + 1;
            if (colindP>=mry){colindP -= mry;}
            colindM = *py_p - ly;
            if (colindM<0){colindM += mry;}
            crP_p = (double*)cr+colindP*mrx;
            crM_p = (double*)cr+colindM*mrx;
            ciP_p = (double*)ci+colindP*mrx;
            ciM_p = (double*)ci+colindM*mrx;
            for(int lx=0;lx<msp;lx++){
                exP = e2xPl[lx]*e3x[lx+1];
                exM = e2xMl[lx]*e3x[lx];
                rowindP = *px_p + lx + 1;
                if (rowindP>=mrx){rowindP -= mrx;} // Looks like these periodic boundary conditions
                rowindM = *px_p - lx;              // take half of the computing time
                if (rowindM<0){rowindM += mrx;}
                *(crP_p+rowindP) += (*v0r_p)*eyP*exP;
                *(crP_p+rowindM) += (*v0r_p)*eyP*exM;
                *(crM_p+rowindP) += (*v0r_p)*eyM*exP;
                *(crM_p+rowindM) += (*v0r_p)*eyM*exM;
                *(ciP_p+rowindP) += (*v0i_p)*eyP*exP;
                *(ciP_p+rowindM) += (*v0i_p)*eyP*exM;
                *(ciM_p+rowindP) += (*v0i_p)*eyM*exP;
                *(ciM_p+rowindM) += (*v0i_p)*eyM*exM;
            }
        }
    }
    free(e2xPl);
    free(e2xMl);
    return;
}

static void backward_gaussian_gridding_real(
        double *cr,
        const double *v0r,
        const int msp,
        const int *px,
        const int *py,
        const int mrx,
        const int mry,
        const int nb,
        const double *e2x,
        const double *e2y,
        const double *e3x,
        const double *e3y)
{
    double *crP_p,*crM_p;
    double *v0r_p=(double*)v0r;
    int *px_p=(int*)px,*py_p=(int*)py;
    int colindM,colindP;
    int rowindM,rowindP;
    double *e2x_p=(double*)e2x,*e2y_p=(double*)e2y;
    double *e2xPl=(double*)malloc(msp*sizeof(double));
    double *e2xMl=(double*)malloc(msp*sizeof(double));
    double e2xPinc,exP; // in 'xPl', 'xMl', 'yPl', 'yMl',
    double e2xMinc,exM; // x means 'first dimension', y is the second
    double e2yPl,e2yPinc,eyP; // P means 'plus' and M means 'minus'
    double e2yMl,e2yMinc,eyM; // l is the variable varying from -Msp+1 to Msp
    
    for(int i=0;i<nb;i++,v0r_p++,px_p++,py_p++,e2x_p++,e2y_p++){
        e2xPinc = *e2x_p;
        e2xMinc = 1.0/(*e2x_p);
        e2yPinc = *e2y_p;
        e2yMinc = 1.0/(*e2y_p);
        e2yPl = e2yPinc;
        e2yMl = 1.0;
        e2xPl[0] = e2xPinc;
        e2xMl[0] = 1.0;
        for(int lx=1;lx<msp;lx++){
            e2xPl[lx]=e2xPl[lx-1]*e2xPinc;
            e2xMl[lx]=e2xMl[lx-1]*e2xMinc;
        }
        for(int ly=0;ly<msp;ly++,e2yPl*=e2yPinc,e2yMl*=e2yMinc){
            eyP = e2yPl*e3y[ly+1];
            eyM = e2yMl*e3y[ly];
            colindP = *py_p + ly + 1;
            if (colindP>=mry){colindP -= mry;}
            colindM = *py_p - ly;
            if (colindM<0){colindM += mry;}
            crP_p = (double*)cr+colindP*mrx;
            crM_p = (double*)cr+colindM*mrx;
            for(int lx=0;lx<msp;lx++){
                exP = e2xPl[lx]*e3x[lx+1];
                exM = e2xMl[lx]*e3x[lx];
                rowindP = *px_p + lx + 1;
                if (rowindP>=mrx){rowindP -= mrx;} // Looks like these periodic boundary conditions
                rowindM = *px_p - lx;              // take half of the computing time
                if (rowindM<0){rowindM += mrx;}
                *(crP_p+rowindP) += (*v0r_p)*eyP*exP;
                *(crP_p+rowindM) += (*v0r_p)*eyP*exM;
                *(crM_p+rowindP) += (*v0r_p)*eyM*exP;
                *(crM_p+rowindM) += (*v0r_p)*eyM*exM;
            }
        }
    }
    free(e2xPl);
    free(e2xMl);
    return;
}