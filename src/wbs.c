#include <R.h>
#include <Rmath.h>

#define IDX(i,j,ld) ((((j)-1) * (ld))+((i)-1))
void ipi_arg_max(double *res, int *n, int *ipargmax, double *ipmax){
	/*consider a version with median!*/
	int i;
	double tmp,max_fabs;

	ipargmax[0] = 0;
	max_fabs = 0;

	for(i=0; i < (*n)-1; i++){
		tmp = fabs(res[i]);

		if(tmp>max_fabs){
			ipargmax[0] = i;
			max_fabs = tmp;
		}
	}


		


	ipmax[0] = res[*ipargmax];
	
}


void wbs_ipi(double *x, int *n, double *res, double *iplus, double *iminus, int *ipargmax, double *ipmax) {
     double sumx, factor;
     int i;
     double n_dbl = (double)(*n);
     double one_over_n = 1.0 / n_dbl;
     double n_squared = n_dbl * n_dbl;
     double i_dbl;
     double iplusone_inv;

     
     sumx = 0;
     for (i=1; i < *n; i++) {
         sumx += x[i];
     }
     
     iminus[0] = 1.0/sqrt(n_squared - n_dbl) * sumx;
     iplus[0] = sqrt(1.0 - one_over_n) * x[0];
     res[0] = iplus[0] - iminus[0];
     
     for (i=1; i < (*n) - 1; i++) {
	 i_dbl = (double)i;
	 iplusone_inv =  1.0/(i_dbl+1.0);

         factor = sqrt((n_dbl-i_dbl-1.0) * i_dbl * iplusone_inv / ( n_dbl-i_dbl) );
         iplus[i] = iplus[i-1] * factor + x[i] * sqrt(iplusone_inv - one_over_n);
         iminus[i] = iminus[i-1] / factor - x[i] / sqrt(n_squared * iplusone_inv - n_dbl);
         res[i] = iplus[i] - iminus[i];
     }
     
     ipi_arg_max(res,n,ipargmax,ipmax);

}



void bs_rec(double *x, int *n, int s, int e, int *len, double *res,double *iplus, double *iminus, double *ipres,int *ipargmax, double *ipmax, double minth){

	*len = e-s+1;
	int cptcand;

	if(*len>1){
			
		
		wbs_ipi(&x[s-1], len, ipres, iplus, iminus,ipargmax,ipmax);

		cptcand= *ipargmax +s;
		
		res[IDX(cptcand,1,*n-1)] = (double)s;
		res[IDX(cptcand,2,*n-1)] = (double)e;
		res[IDX(cptcand,3,*n-1)] = (double)cptcand;
		res[IDX(cptcand,4,*n-1)] = (double)(*ipmax);

		if(minth>fabs(*ipmax) || minth < 0){
			minth = fabs(*ipmax);
		}

		res[IDX(cptcand,5,*n-1)] = minth;

		bs_rec(x, n, s, cptcand, len, res,iplus, iminus,ipres,ipargmax, ipmax, minth);
		bs_rec(x, n, cptcand+1, e, len, res,iplus, iminus,ipres,ipargmax, ipmax, minth);

	}

}

void bs_rec_wrapper(double *x, int *n, double *res){
	
	double *iplus= Calloc(*n-1,double);
	double *iminus= Calloc(*n-1,double);
	double *ipres= Calloc(*n-1,double);
	double *ipmax = Calloc(1,double);
	int *ipargmax = Calloc(1,int);
	int *len = Calloc(1,int);
	/*negative value of minth serves as infinity*/
	bs_rec(x, n, 1, *n, len, res,iplus, iminus,ipres,ipargmax, ipmax, -1.0);

	Free(iplus);
	Free(iminus);
	Free(ipres);
	Free(ipmax);
	Free(ipargmax);
	Free(len);
}

void wbs_int_rec(double *x, int *n, int s, int e, int *len, double *res, double *iplus, double *iminus, double *ipres, int *ipargmax, double *ipmax, double *wbsres, int *index, int indexn, int *M, double minth){
	*len = e-s+1;

	if(*len>1){
		if(indexn > 0){
			int cptcand;
			wbs_ipi(&x[s-1], len, ipres, iplus, iminus,ipargmax,ipmax);

			if(fabs(*ipmax)< (wbsres[IDX(index[0],5,*M)])){
				cptcand= wbsres[IDX(index[0],3,*M)];
				res[IDX(cptcand,1,*n-1)] = (double)s;
				res[IDX(cptcand,2,*n-1)] = (double)e;
				res[IDX(cptcand,3,*n-1)] = (double)cptcand;
				res[IDX(cptcand,4,*n-1)] = (double)(wbsres[IDX(index[0],4,*M)]);
				
				if(minth>wbsres[IDX(index[0],5,*M)] || minth < 0){
					minth = wbsres[IDX(index[0],5,*M)];
				}

			}else{
				cptcand= *ipargmax +s;
				res[IDX(cptcand,1,*n-1)] = (double)s;
				res[IDX(cptcand,2,*n-1)] = (double)e;
				res[IDX(cptcand,3,*n-1)] = (double)cptcand;
				res[IDX(cptcand,4,*n-1)] = (double)(*ipmax);

				if(minth>fabs(*ipmax) || minth < 0){
					minth = fabs(*ipmax);
				}
			}

			res[IDX(cptcand,5,*n-1)] = minth;

			/*left*/
			int *indexl= Calloc(indexn,int);
			int *indexr= Calloc(indexn,int);
			int indexnl = 0;
			int indexnr = 0;
			int i;

			for(i=1; i<=indexn; i++){
				if((wbsres[IDX(index[i-1],1,*M)] >= s) & (wbsres[IDX(index[i-1],2,*M)]  <= cptcand) ){
					indexl[indexnl] = index[i-1];
					indexnl ++;
				}else if((wbsres[IDX(index[i-1],1,*M)] >= (cptcand+1)) & (wbsres[IDX(index[i-1],2,*M)]  <= e) ){
					indexr[indexnr] = index[i-1];
					indexnr ++;
				}
			}

			if(indexnl){
				indexl = Realloc(indexl,indexnl,int);
	
				wbs_int_rec(x, n, s, cptcand, len, res, iplus, iminus, ipres, ipargmax, ipmax, wbsres, indexl, indexnl, M,minth);
				Free(indexl);
			}else{
				Free(indexl);

				bs_rec(x, n, s, cptcand, len, res,iplus, iminus,ipres,ipargmax, ipmax,minth);
			}

			if(indexnr){
				indexr = Realloc(indexr,indexnr,int);
		
				wbs_int_rec(x, n, cptcand +1, e, len, res, iplus, iminus, ipres, ipargmax, ipmax, wbsres, indexr, indexnr, M,minth);
				Free(indexr);
			}else{
				Free(indexr);
		
				bs_rec(x, n, cptcand +1, e, len, res,iplus, iminus,ipres,ipargmax, ipmax,minth);
			}
			
		}else{
	
			bs_rec(x, n, s, e, len, res,iplus, iminus,ipres,ipargmax, ipmax,minth);
		}
	}
}
void wbs_int_rec_wrapper(double *x, int *n, double *res,  int *intervals, int *M){
	double *iplus= Calloc(*n-1,double);
	double *iminus= Calloc(*n-1,double);
	double *ipres= Calloc(*n-1,double);
	double *wbsres= Calloc(*M * 5,double);
	double *ipmax = Calloc(1,double);
	int *index = Calloc(*M,int);
	int *ipargmax = Calloc(1,int);
	int *len = Calloc(1,int);
	int i,s,e,cptcand;
	
	/* find cpt candidates on given intervals*/

	for(i=1; i<=*M;i++){
		s = intervals[IDX(i,1,*M)];
		e = intervals[IDX(i,2,*M)];
		len[0] = e-s+1;

		wbs_ipi(&x[s-1], len, ipres, iplus, iminus,ipargmax,ipmax);
		cptcand= *ipargmax +s;

		wbsres[IDX(i,1,*M)] = (double)s;
		wbsres[IDX(i,2,*M)] = (double)e;
		wbsres[IDX(i,3,*M)] = (double)cptcand;
		wbsres[IDX(i,4,*M)] = (double)(*ipmax);
		wbsres[IDX(i,5,*M)] = fabs(*ipmax);
		index[i-1] = i;
	}
	/* sort elementd from the one with the largest abs(cusum)*/
	double *tmp= Calloc(*M,double);

	memcpy(tmp,&wbsres[IDX(1,5,*M)],*M * sizeof(double));

	revsort(tmp,index,*M);

	Free(tmp);
	


	/*standard wbs part*/
	
	wbs_int_rec(x, n, 1, *n, len, res, iplus, iminus, ipres, ipargmax, ipmax, wbsres, index, *M, M,-1.0);

	Free(iplus);
	Free(iminus);
	Free(ipres);
	Free(ipmax);
	Free(index);
	Free(ipargmax);
	Free(len);
	Free(wbsres);
}



void wbs(double *x, int *n, double *res, int *intervals, int *M){
	double *iplus= Calloc(*n-1,double);
	double *iminus= Calloc(*n-1,double);
	double *ipres= Calloc(*n-1,double);
	double *ipmax = Calloc(1,double);
	int *ipargmax = Calloc(1,int);
	int *len = Calloc(1,int);
	int i,s,e,cptcand;

	for(i=1; i<=*M;i++){
		s = intervals[IDX(i,1,*M)];
		e = intervals[IDX(i,2,*M)];
		len[0] = e-s+1;

		wbs_ipi(&x[s-1], len, ipres, iplus, iminus,ipargmax,ipmax);
		cptcand= *ipargmax +s;

		res[IDX(i,1,*M)] = (double)s;
		res[IDX(i,2,*M)] = (double)e;
		res[IDX(i,3,*M)] = (double)cptcand;
		res[IDX(i,4,*M)] = (double)(*ipmax);
		res[IDX(i,5,*M)] = fabs((*ipmax));
	}

	Free(iplus);
	Free(iminus);
	Free(ipres);
	Free(ipmax);
	Free(ipargmax);
	Free(len);
}

