#include "wbs.h"

void ipi_arg_max(double *res, int n, int *ipargmax, double *ipmax){

	int i,k;
	int max_count=0;
	double tmp,max_fabs;
	 

	ipargmax[0] = 0;
	max_fabs = -1;

	for(i=0; i < n-1; i++){
		tmp = fabs(res[i]);

		if(tmp>max_fabs){
			ipargmax[0] = i;
			max_fabs = tmp;
			max_count = 1;
		}else if(tmp==max_fabs) max_count++;
	}
	
	/*if there are multiple points maximizing cusums, we take the median*/

	if(max_count>1){
		max_count = max_count/2 + (int)(max_count%2);
		
		k=0;
		i=0;

		while((i< (n-1))  && (k<max_count)){
			i++;
			if(fabs(res[i])==max_fabs) k++;
		}

		ipargmax[0] = i;
	}

		


	ipmax[0] = res[*ipargmax];
	
}

void wbs_ipi(double *x, int n, double *res, double *iplus, double *iminus, int *ipargmax, double *ipmax) {
  
     double sumx, factor;
     int i;
     double n_dbl = (double)(n);
     double one_over_n = 1.0 / n_dbl;
     double n_squared = n_dbl * n_dbl;
     double i_dbl;
     double iplusone_inv;
     


     
     sumx = 0;
     for (i=1; i < n; i++) {
         sumx += x[i];
     }
     
     iminus[0] = 1.0/sqrt(n_squared - n_dbl) * sumx;
     iplus[0] = sqrt(1.0 - one_over_n) * x[0];
     res[0] = iplus[0] - iminus[0];
     
     for (i=1; i < n - 1; i++) {
	       i_dbl = (double)i;
	       iplusone_inv =  1.0/(i_dbl+1.0);

         factor = sqrt((n_dbl-i_dbl-1.0) * i_dbl * iplusone_inv / ( n_dbl-i_dbl) );
         iplus[i] = iplus[i-1] * factor + x[i] * sqrt(iplusone_inv - one_over_n);
         iminus[i] = iminus[i-1] / factor - x[i] / sqrt(n_squared * iplusone_inv - n_dbl);
         res[i] = iplus[i] - iminus[i];
     }
     
     ipi_arg_max(res,n, ipargmax, ipmax);

}



void bs_rec(double *x, int n, int s, int e, double *res, double *iplus, double *iminus, double *ipres,  double minth, int scale){

	int len = e-s+1;
	int cptcand;
  int ipargmax;
  double ipmax;

	if(len>1){
			
		
		wbs_ipi(&x[s-1], len, ipres, iplus, iminus,&ipargmax,&ipmax);

		cptcand= ipargmax +s;
		
		res[IDX(cptcand,1,n-1)] = (double)s;
		res[IDX(cptcand,2,n-1)] = (double)e;
		res[IDX(cptcand,3,n-1)] = (double)cptcand;
		res[IDX(cptcand,4,n-1)] = (double)(ipmax);

		if(minth>fabs(ipmax) || minth < 0){
			minth = fabs(ipmax);
		}

		res[IDX(cptcand,5, n-1)] = (double)minth;
		res[IDX(cptcand,6, n-1)] = (double)scale;

		bs_rec(x, n, s, cptcand, res,iplus, iminus,ipres, minth, scale+1);
		bs_rec(x, n, cptcand+1, e,  res,iplus, iminus,ipres, minth, scale+1);

	}

}

void bs_rec_wrapper(double *x, int *n, double *res){
	
	double *iplus= Calloc(*n-1,double);
	double *iminus= Calloc(*n-1,double);
	double *ipres= Calloc(*n-1,double);
	
	/*negative value of minth serves as infinity*/
	bs_rec(x, *n, 1, *n,  res,iplus, iminus,ipres, -1.0, 1);

	Free(iplus);
	Free(iminus);
	Free(ipres);
  
}

void wbs_int_rec(double *x, int n, int s, int e, double *res, double *iplus, double *iminus, double *ipres, double *wbsres, int *index, int indexn, int M, double minth, int scale){
  
	int len = e-s+1;
  int ipargmax;
  double ipmax;
  
	if(len>1){
		if(indexn > 0){
			int cptcand;
		
      wbs_ipi(&x[s-1], len, ipres, iplus, iminus,&ipargmax,&ipmax);

			if(fabs(ipmax) < (wbsres[IDX(index[0],5,M)])){
        
				cptcand= wbsres[IDX(index[0],3,M)];
				
        res[IDX(cptcand,1,n-1)] = (double)s;
				res[IDX(cptcand,2,n-1)] = (double)e;
				res[IDX(cptcand,3,n-1)] = (double)cptcand;
				res[IDX(cptcand,4,n-1)] = (double)(wbsres[IDX(index[0],4,M)]);
				
				if(minth>wbsres[IDX(index[0],5,M)] || minth < 0){
					minth = wbsres[IDX(index[0],5,M)];
				}
				
			}else{
        
				cptcand= ipargmax +s;
				res[IDX(cptcand,1,n-1)] = (double)s;
				res[IDX(cptcand,2,n-1)] = (double)e;
				res[IDX(cptcand,3,n-1)] = (double)cptcand;
				res[IDX(cptcand,4,n-1)] = (double)(ipmax);

				if(minth>fabs(ipmax) || minth < 0){
					minth = fabs(ipmax);
				}
			}

			res[IDX(cptcand,5,n-1)] = minth;
			res[IDX(cptcand,6,n-1)] = (double)scale;

			/*left*/
			int *indexl= Calloc(indexn,int);
			int *indexr= Calloc(indexn,int);
			int indexnl = 0;
			int indexnr = 0;
			int i;

			for(i=1; i<=indexn; i++){
				if((wbsres[IDX(index[i-1],1,M)] >= s) & (wbsres[IDX(index[i-1],2,M)]  <= cptcand) ){
					indexl[indexnl] = index[i-1];
					indexnl ++;
				}else if((wbsres[IDX(index[i-1],1,M)] >= (cptcand+1)) & (wbsres[IDX(index[i-1],2,M)]  <= e) ){
					indexr[indexnr] = index[i-1];
					indexnr ++;
				}
			}

			if(indexnl){
				indexl = Realloc(indexl,indexnl,int);
        wbs_int_rec(x, n, s, cptcand, res, iplus, iminus, ipres,  wbsres, indexl, indexnl, M,minth, scale+1);
				Free(indexl);
        
			}else{
				Free(indexl);

				bs_rec(x, n, s, cptcand,  res,iplus, iminus,ipres,minth, scale+1);
			}

			if(indexnr){
				indexr = Realloc(indexr,indexnr,int);
		
				wbs_int_rec(x, n, cptcand +1, e, res, iplus, iminus, ipres, wbsres, indexr, indexnr, M,minth,scale+1);
				Free(indexr);
        
			}else{
				
        Free(indexr);		
				bs_rec(x, n, cptcand +1, e, res, iplus, iminus,ipres,minth,scale+1);
        
			}
			
		}else{
	
			bs_rec(x, n, s, e, res,iplus, iminus,ipres,minth, scale);
		}
	}
}
void wbs_int_rec_wrapper(double *x, int *n, double *res,  int *intervals, int *M){
	double *iplus= Calloc(*n-1,double);
	double *iminus= Calloc(*n-1,double);
	double *ipres= Calloc(*n-1,double);
	double *wbsres= Calloc(*M * 5,double);
	int *index = Calloc(*M,int);
	int i,s,e,cptcand;
	int ipargmax;
  double ipmax;
	/* find cpt candidates on given intervals*/

	for(i=1; i<=*M;i++){
		s = intervals[IDX(i,1,*M)];
		e = intervals[IDX(i,2,*M)];

		wbs_ipi(&x[s-1], e-s+1, ipres, iplus, iminus,&ipargmax,&ipmax);
		
    cptcand= ipargmax +s;

		wbsres[IDX(i,1,*M)] = (double)s;
		wbsres[IDX(i,2,*M)] = (double)e;
		wbsres[IDX(i,3,*M)] = (double)cptcand;
		wbsres[IDX(i,4,*M)] = (double)(ipmax);
		wbsres[IDX(i,5,*M)] = fabs(ipmax);
		index[i-1] = i;
	}
	/* sort elementd from the one with the largest abs(cusum)*/
  
	double *tmp= Calloc(*M,double);

	memcpy(tmp,&wbsres[IDX(1,5,*M)],*M * sizeof(double));

	revsort(tmp,index,*M);

	Free(tmp);
	


	/*standard wbs part*/
	
	wbs_int_rec(x, *n, 1, *n, res, iplus, iminus, ipres, wbsres, index, *M, *M,-1.0,1);

	Free(iplus);
	Free(iminus);
	Free(ipres);
  Free(index);	
	Free(wbsres);
  
}


void wbs_rec(double *x, int n, int s, int e,  double *res, double *wbsres, int *index, int indexn, int M, int scale){
	
  int len = e-s+1;

	if(len>1){
	
    if(indexn > 0){
		
      int cptcand;
	
			cptcand= (int) wbsres[IDX(index[0],3,M)];
	
			res[IDX(cptcand,1,n-1)] = (double)(wbsres[IDX(index[0],1,M)]);
			res[IDX(cptcand,2,n-1)] = (double)(wbsres[IDX(index[0],2,M)]);
			res[IDX(cptcand,3,n-1)] = (double)cptcand;
			res[IDX(cptcand,4,n-1)] = (double)(wbsres[IDX(index[0],4,M)]);	
			res[IDX(cptcand,5,n-1)] = wbsres[IDX(index[0],5,M)];
			res[IDX(cptcand,6,n-1)] = (double)scale;

			/*left*/
			int *indexl = Calloc(indexn,int);
			int *indexr = Calloc(indexn,int);
			int indexnl = 0;
			int indexnr = 0;
			int i;

			for(i=1; i<=indexn; i++){
				if((wbsres[IDX(index[i-1],1,M)] >= s) & (wbsres[IDX(index[i-1],2,M)]  <= cptcand) ){
					indexl[indexnl] = index[i-1];
					indexnl ++;
				}else if((wbsres[IDX(index[i-1],1,M)] >= (cptcand+1)) & (wbsres[IDX(index[i-1],2,M)]  <= e) ){
					indexr[indexnr] = index[i-1];
					indexnr ++;
				}
			}

			if(indexnl){
				indexl = Realloc(indexl,indexnl,int);
				wbs_rec(x, n, s, cptcand, res, wbsres, indexl, indexnl, M, scale+1);
				Free(indexl);
			}
      
			if(indexnr){
				indexr = Realloc(indexr,indexnr,int);
				wbs_rec(x, n, cptcand +1, e, res,  wbsres, indexr, indexnr, M, scale+1);
				Free(indexr);
			}
			
		}
	}
}

void wbs_rec_wrapper(double *x, int *n, double *res,  int *intervals, int *M){
  
	double *iplus= Calloc(*n-1,double);
	double *iminus= Calloc(*n-1,double);
	double *ipres= Calloc(*n-1,double);
	double *wbsres= Calloc(*M * 5,double);
	double ipmax;
	int *index = Calloc(*M,int);
	int ipargmax, i,s,e,cptcand;
	
	/* find cpt candidates on given intervals*/

	for(i=1; i<=*M;i++){
    
		s = intervals[IDX(i,1,*M)];
		e = intervals[IDX(i,2,*M)];

		wbs_ipi(&x[s-1], e-s+1, ipres, iplus, iminus, &ipargmax,& ipmax);
		cptcand= ipargmax +s;

		wbsres[IDX(i,1,*M)] = (double)s;
		wbsres[IDX(i,2,*M)] = (double)e;
		wbsres[IDX(i,3,*M)] = (double)cptcand;
		wbsres[IDX(i,4,*M)] = (double)(ipmax);
		wbsres[IDX(i,5,*M)] = fabs(ipmax);
		index[i-1] = i;
    
	}
	/* sort elementd from the one with the largest abs(cusum)*/
	double *tmp= Calloc(*M,double);

	memcpy(tmp,&wbsres[IDX(1,5,*M)],*M * sizeof(double));

	revsort(tmp,index,*M);

	Free(tmp);
	


	/*standard wbs part*/
	
	wbs_rec(x, *n, 1, *n, res, wbsres, index, *M, *M, 1);

	Free(iplus);
	Free(iminus);
	Free(ipres);	
	Free(index);
	Free(wbsres);
  
}
