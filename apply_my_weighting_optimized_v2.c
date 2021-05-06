#include <cstdlib>
#include <cmath>
#include <ctype.h>   // for toupper()  */
#include <string.h>  /* for strings */
#include <time.h>
#include <ctime>
#include <string>
#include <iostream>  //  C++ stream IO
#include <sstream>   // string streams
#include <fstream>
#include <iomanip>   //  to format the output
#include <vector>
#include "EMData.h"
//#include <math.h>
#include "List.h"
#include "util.h"
#include "XYData.h"
#include <cstdio>
#include <stdlib.h>
#include <cstdlib>
#include <unistd.h>
using namespace std;

float CTF_AST (int x1, int y1, int ny, float ds, float dfu, float dfv, float dfdiff, float dfang ,float lambda, float cs, float ampconst, int mode){
	float v,ss,ag,gamma,df_ast;
	if (mode==0){
		ss=hypot((float)x1,(float)y1-ny/2)*ds;
		ag=atan2(float(y1-ny/2),float(x1));
		df_ast=0.5*(dfu+dfv+2*dfdiff*cos(2*(dfang*PI/180-ag)));
		gamma=-2*PI*(cs*2.5e6*lambda*lambda*lambda*ss*ss*ss*ss+df_ast*5000.0*lambda*ss*ss);
		v=(sqrt(1.0-ampconst*ampconst)*sin(gamma)+ampconst*cos(gamma))>0?1.0:-1.0;		//do flipphase
	}
	else if (mode==2){
		ss=hypot((float)x1,(float)y1-ny/2)*ds;
		ag=atan2(float(y1-ny/2),float(x1));
		df_ast=0.5*(dfu+dfv+2*dfdiff*cos(2*(dfang*PI/180-ag)));
		gamma=-2*PI*(cs*2.5e6*lambda*lambda*lambda*ss*ss*ss*ss+df_ast*5000.0*lambda*ss*ss);
		v=fabs(sqrt(1.0-ampconst*ampconst)*sin(gamma)+ampconst*cos(gamma));		//	return abs ctf value
	}
	else{
		ss=hypot((float)x1,(float)y1-ny/2)*ds;
		ag=atan2(float(y1-ny/2),float(x1));
		df_ast=0.5*(dfu+dfv+2*dfdiff*cos(2*(dfang*PI/180-ag)));
		gamma=-2*PI*(cs*2.5e6*lambda*lambda*lambda*ss*ss*ss*ss+df_ast*5000.0*lambda*ss*ss);
		v=(sqrt(1.0-ampconst*ampconst)*sin(gamma)+ampconst*cos(gamma));		//	return ctf value
	}
	return v;
}

int main(int argc, char *argv[]){
	char inlst[512],outlst[512],snr[512];
	string option;
	float apix,edge_half_width=4,outer_Radius=0,defocus=0,dfdiff=0,dfang=90;
	float energy=300,cs=2.7,ampconst=0.07,ds,lambda,gamma,ss,v,dfu,dfv,ag;
	float *parm,cc_score=0,lowres=30,highres=8,thres,phi_step;
	float a,b,b2,bfactor,bfactor2,bfactor3,k=0,kk=3,Bf;
	int count=0,j,x1,y1,x,y,n,first,last,ny_2d,nx2,ny2,nz,nx,ny,yy,d_m;
	int sizex,sizey,tmp_nx;
	float theta0,phi0,delta_theta,delta_phi,r,d;
	float x_norm,y_norm,z_norm,xx0,yy0,zz0;
	int current=0;
	const char *filespec;
	if(argc==1){
		printf("This program apply optimized weights to raw images.\n");
		printf("--input            input 2d raw images lstfile with ctf information\n");
		printf("--angpix            input pixel size in angstroms\n");
		printf("--weight           optimized weight with SSNR parmeters.\n");
		printf("--kk               overlapping density parameter, default is 3 \n");
		printf("--first            the first image to process.\n");
		printf("--last             the last image to process.\n");
		printf("--energy           accerlerating voltage in kV.\n");
		printf("--cs               spherical aberration in um.\n");
		printf("--Highres          high resolution cut \n");
		printf("--Lowres           low resolution cut\n");
		printf("--radius         mask radius in pixel\n");
		printf("--edge_half_width        soft edge width \n");
		printf("--output           output file rootname\n");
		return 0;
	}
	for(int i=1;i<argc;i++){
		option=argv[i];
		if(option.compare("--input")==0){
			i++;
			strcpy(inlst, argv[i]);
			continue;
		}
		option=argv[i];
		if(option.compare("--output")==0){
			i++;
			strcpy(outlst, argv[i]);
			continue;
		}	
		if(option.compare("--weight")==0){
			i++;
			strcpy(snr, argv[i]);
			continue;
		}
		if(option.compare("--kk")==0){
			i++;
			kk=atof(argv[i]);
			continue;
		}
		if(option.compare("--angpix")==0){
			i++;
			apix=atof(argv[i]);
			continue;
		}
		if(option.compare("--first")==0){
			i++;
			first=atoi(argv[i]);
			continue;
		}
		if(option.compare("--last")==0){
			i++;
			last=atoi(argv[i]);
			continue;
		}
		if(option.compare("--energy")==0){
			i++;
			energy=atof(argv[i]);
			continue;
		}
		if(option.compare("--cs")==0){
			i++;
			cs=atof(argv[i]);
			continue;
		}
		if(option.compare("--Highres")==0){
			i++;
			highres=atof(argv[i]);
			continue;
		}
		if(option.compare("--Lowres")==0){
			i++;
			lowres=atof(argv[i]);
			continue;
		}
		if(option.compare("--radius")==0){
			i++;
			outer_Radius=atof(argv[i]);
			continue;
		}
		if(option.compare("--edge_half_width")==0){
			i++;
			edge_half_width=atof(argv[i]);
			continue;
		}
		printf("Undefined option: %s .Abort.\n",argv[i]);
	}
	std::ifstream wiener_parm(snr);
	float m[5000]={0},rc[5000]={0},signal[5000]={0};
	string str1="_input.lst";
	string strlst=outlst+str1;
	string str2;
	FILE *g=fopen(strlst.c_str(),"wb");
    for(n=first;n<last;n++){
		while(!wiener_parm.eof()){
			if(current==n){
				getline(wiener_parm, str2,'\n');
				sscanf(str2.c_str(),"%f\t%f\t%f\t%f\t%f\t%f",&a,&b,&b2,&bfactor,&bfactor2,&bfactor3);//cout<<a<<"  "<<b<<"  "<<b2<<"  "<<bfactor<<endl;
				break;
			}
			current += 1;
		}
	    char s[MAXPATHLEN+10],t[MAXPATHLEN+1],u[256],*c;
		filespec=inlst;
		FILE *in = fopen(filespec,"rb");
		if (!in) { printf("'%s' not found\n",filespec); return (-3); }
		int i,nn,ret=0;
		u[0] = 0;
		fgets(s,MAXPATHLEN+10,in);
		if (strncmp(s,"#LSX",4)==0) {
			int step=0;
			fgets(s,MAXPATHLEN+10,in);
			fgets(s,MAXPATHLEN+10,in);
			sscanf(s+1," %d",&step);
			if (step==0) {fclose(in); return -1; }
			int tl=ftell(in);
			if (n>0) fseek(in,step*n,SEEK_CUR);
			if (!fgets(s,MAXPATHLEN+10,in)) { fclose(in); return -2; }
			if ((ftell(in)-tl)%step !=0) {
				printf("Error in #LSX file detected. Please rerun lstfast.py\n");
				return(-3);
			}
			u[0]='\0';
			sscanf(s," %d %s %[^\n]",&nn,t,u);
			fclose(in);
		}
		else {
			if (n<0) n=0;
			t[0]=0;
			for (i=0; i<=n; i++) {
				u[0]='\0';
				if (!fgets(s,MAXPATHLEN+10,in)) break;
				if (s[0]=='#') i--;
				else sscanf(s," %d %s %[^\n]",&nn,t,u);
			}
			fclose(in);
			if (i<=n) return (-1);
		}
		vector<string> pairs;
		split(u,pairs,0,0);
		EMData *i2d_emdata=new EMData();
		i2d_emdata->readImage(t,nn);
		cout<<t<<"\t"<<nn<<endl;
		//i2d_emdata->writeImage("debug.hdf",-1);
		nx2=i2d_emdata->xSize();
		ny2=i2d_emdata->ySize();
		cout<<ny2<<endl;
		for(int i=0; i<pairs.size(); i++){
            if(strncmp(pairs[i].c_str(),"defocus=",8)==0){
			  	if(sscanf(pairs[i].substr(8).c_str(),"%f",&defocus)==1){
					defocus=defocus*(-1);
				}
		    }
			else if(strncmp(pairs[i].c_str(),"dfdiff=",7)==0){
				if(sscanf(pairs[i].substr(7).c_str(),"%f",&dfdiff)==1){
				}
			}
         /*   else if(strncmp(pairs[i].c_str(),"euler=",6)==0){
				if(sscanf(pairs[i].substr(6).c_str(),"%f,%f,%f",&euler1,&euler2,&euler3)==3){
				}
			}*/
			else if(strncmp(pairs[i].c_str(),"dfang=",6)==0){
				if(sscanf(pairs[i].substr(6).c_str(),"%f",&dfang)==1){
				}
			}
			else {
			}
	    }
		ampconst=0.07;
		dfu=defocus+dfdiff; //defocus is minus, so abs(dfu) < abs(dfv)
		dfv=defocus-dfdiff;
		lambda=12.2639/sqrt(energy*1000.0+0.97845*energy*energy);
		ds=1/(apix*ny2);
		//i2d_emdata->writeImage("debug.hdf",0);
		EMData *i2d_emdata_FFT=new EMData();
	    i2d_emdata_FFT=i2d_emdata->doFFT();
	    i2d_emdata_FFT->ap2ri();
	    EMData *filter=new EMData();
	    filter->setSize(ny2+2,ny2,1);
	    filter->setComplex(1);
	    filter->setRI(1);
	    filter->one();
        float *filter_data=filter->getData();
	    for (y1=0; y1<ny2; y1++) {
		    for (x1=0; x1<(ny2+2)/2; x1++) {
				r=sqrt(SQR(x1)+SQR(y1-ny2/2));
				ss=r*ds;
		        v=CTF_AST(x1,y1,ny2,ds,dfu,dfv,dfdiff,dfang,lambda,cs,ampconst,0);
		        filter_data[x1*2+y1*(ny2+2)]*=v;
		        filter_data[x1*2+y1*(ny2+2)+1]=filter_data[x1*2+y1*(ny2+2)];
		    }
	    }       //has astimatism right now
	    filter->doneData();
	    filter->update();
	    i2d_emdata_FFT->mult(filter);
	    delete filter;
		EMData *i2d_emdata_FFT_ift=new EMData();
		i2d_emdata_FFT_ift=i2d_emdata_FFT->doIFT();
		i2d_emdata_FFT->gimmeFFT();
		delete i2d_emdata_FFT;
		//i2d_emdata_FFT->toCorner();
		////////////////////apply mask
		float *i2d_emdata_FFT_ift_data=i2d_emdata_FFT_ift->getData();
		for (int jj=0; jj<ny2; jj++){
			for (int ii=0; ii<ny2; ii++){
				r=sqrt(SQR(ii-ny2/2)+SQR(jj-ny2/2));
				if( r>(outer_Radius+2*edge_half_width)){
			        i2d_emdata_FFT_ift_data[ii+jj*ny2]=0;
			    }
			    else if (r<outer_Radius){
				     i2d_emdata_FFT_ift_data[ii+jj*ny2]=i2d_emdata_FFT_ift_data[ii+jj*ny2];
				}
			    else {
		            d=d=0.5*cos(PI*(r-outer_Radius)/(2*edge_half_width))+0.5;
			        i2d_emdata_FFT_ift_data[ii+jj*ny2]=i2d_emdata_FFT_ift_data[ii+jj*ny2]*d;
			    }
		    }
		}
		i2d_emdata_FFT_ift->doneData();
		i2d_emdata_FFT_ift->update();
		////////////////////////////////////////////////////////////
		/*float *rd = (float *)malloc(ny2 / 2 * sizeof(float));
	    float *ra = (float *)malloc(ny2 / 2 * sizeof(float));
		memset(rd, 0, ny2 / 2 * sizeof(float));
	    memset(ra, 0, ny2 / 2 * sizeof(float));*/
		EMData *i2d_emdata_FFT_cpy=new EMData();
		i2d_emdata_FFT_cpy=i2d_emdata_FFT_ift->doFFT();
		i2d_emdata_FFT_cpy->ri2ap();
		i2d_emdata_FFT_ift->gimmeFFT();
		delete i2d_emdata_FFT_ift;
		int fpp_xsize=i2d_emdata_FFT_cpy->xSize();
		int fpp_ysize=i2d_emdata_FFT_cpy->ySize();
		float ra1[5000]={0},rb1[5000]={0};
		float ra2[5000]={0},rb2[5000]={0};
		float Ncurve[5000]={0};
		float *i2d_emdata_FFT_cpy_data=i2d_emdata_FFT_cpy->getData();
		int xt1=ny2/2,xt2=ny2/2;
		for (int jj=0, y2=0;jj<fpp_ysize;jj++){
			for(int ii=0;ii<fpp_xsize;ii+=2,y2+=2){
				r = sqrt(SQR(ii / 2) + SQR(jj - fpp_ysize / 2));
				x = floor(r + 0.5) - 1;
				ss=r*ds;
				signal[x]=(exp(bfactor*ss*ss+bfactor2*ss+bfactor3));
				Ncurve[x] = exp(a*ss*ss+b*ss+b2);
				if (x<ny2/2 && x >= 0) {
					ra1[x]+=SQR(i2d_emdata_FFT_cpy_data[y2]);
					rb1[x]+=1;
				}
			}
		}
        for (int jj=0, y2=0;jj<fpp_ysize;jj++){
			for(int ii=0;ii<fpp_xsize;ii+=2,y2+=2){
				r = sqrt(SQR(ii / 2) + SQR(jj - fpp_ysize / 2));
				x = floor(r + 0.5) - 1;
				ss=r*ds;
				if (x<ny2/2 && x >= 0) {
					v=CTF_AST(ii/2,jj,ny2,ds,dfu,dfv,dfdiff,dfang,lambda,cs,ampconst,2);
					//if(abs(v)<0.01 && x<xt1 && x>0) xt1=x;
				   // if(abs(v)>0.98 && x>(xt1) && x<xt2) xt2=x;
					i2d_emdata_FFT_cpy_data[y2]=i2d_emdata_FFT_cpy_data[y2]*sqrt((signal[x]*v*v+Ncurve[x])/signal[x])/sqrt(ra1[x]/rb1[x]);
                    //if(x<xt1 ) i2d_emdata_FFT_cpy_data[y2]=i2d_emdata_FFT_cpy_data[y2]*sqrt(ra1[xt1]/rb1[xt1])/sqrt(ra1[x]/rb1[x]);
				    if(x>(ny2*apix/6)) i2d_emdata_FFT_cpy_data[y2]=i2d_emdata_FFT_cpy_data[y2]*exp(-100*ss*ss);
					//rd[x] += SQR(i2d_emdata_FFT_cpy_data[y2])*100;
				}
		    }
	    } 
		///////////////////////////////////////////do lowpass and highpass below
		for (int jj=0, y2=0;jj<ny2;jj++){
			for(int ii=0;ii<ny2+2;ii+=2,y2+=2){
				r = sqrt(SQR(ii / 2) + SQR(jj - ny2 / 2));
				x = floor(r + 0.5) - 1;
				if (x<ny2*apix/highres && x >= ny2*apix/lowres) {
			        i2d_emdata_FFT_cpy_data[y2]=i2d_emdata_FFT_cpy_data[y2];
				}
				else if(x>=ny2*apix/highres && x<ny2*apix/highres+8){
					i2d_emdata_FFT_cpy_data[y2]=i2d_emdata_FFT_cpy_data[y2]*(0.5*cos(PI*(x-ny2*apix/highres)/(2*8))+0.5);
				}
				else if(x>=(ny2*apix/lowres-8) && x<ny2*apix/lowres && x>=0){
					i2d_emdata_FFT_cpy_data[y2]=i2d_emdata_FFT_cpy_data[y2]*(0.5*cos(PI*(ny2*apix/lowres-x)/(2*8))+0.5);
				}
				else
					i2d_emdata_FFT_cpy_data[y2]=0;
			}
		}
		i2d_emdata_FFT_cpy->doneData();
		i2d_emdata_FFT_cpy->update();		
		float infile_mean2,total_n,infile_sigma2;
		float *i2d_emdata_FFT_data2=i2d_emdata_FFT_cpy->getData();
		
		
		//////////////////////////////////////////////////////////////////apply weighting function
		
		int x_temp;
		for(int j=0,y2=0;j<ny2;j++){
			for(int i=0;i<ny2+2;i+=2,y2+=2){
				r = sqrt(SQR(i/2) + SQR(j-ny2/2));
				x=floor(r+0.5)-1;
				ss=r*ds;
				if(x<ny2/2 && x>=0){
					v=CTF_AST(i/2,j,ny2,ds,dfu,dfv,dfdiff,dfang,lambda,cs,ampconst,2);
					signal[x]=(exp(bfactor*ss*ss+bfactor2*ss+bfactor3))/(kk+1);
					Ncurve[x] = exp(a*ss*ss+b*ss+b2);
					//euler_w[x]=1.68118*ss;
					i2d_emdata_FFT_data2[y2]=i2d_emdata_FFT_data2[y2]*sqrt(1/(Ncurve[x]/signal[x]+kk*v*v));
					//rd[x] += (v/(Ncurve[x]/signal[x]+kk*v*v));
					//ra[x] += 1.0;
				}
			}
		}
		////////////////////////////////////////////////////////////////////////////////do normalization below
		for (int jj=0, y2=0;jj<ny2;jj++){
			for(int ii=0;ii<ny2+2;ii+=2,y2+=2){
				r = sqrt(SQR(ii / 2) + SQR(jj - ny2 / 2));
				x = floor(r + 0.5) - 1;
			    ss=r*ds;
				if (x<ny2/2 && x >= 0) {
					infile_sigma2+=SQR(i2d_emdata_FFT_data2[y2]);
					total_n+=1;
				}
			}
		}
		infile_mean2=infile_sigma2/SQR(total_n);
		for (int jj=0, y2=0;jj<ny2;jj++){
			for(int ii=0;ii<ny2+2;ii+=2,y2+=2){
				i2d_emdata_FFT_data2[y2]=i2d_emdata_FFT_data2[y2]/sqrt(infile_mean2);
			}
		}
		i2d_emdata_FFT_cpy->doneData();
		i2d_emdata_FFT_cpy->update();
		////////////////////////////////////////////////////////////////////////////////////////
		/*FILE *fp=fopen(snr_parm.c_str(),"wb");
		for(int i=0;i<ny2/2;i++){
		    fprintf(fp,"%f\t%f\n",i/(ny2*apix),(rd[i]/ra[i]));
		}
		fclose(fp);*/
		char cc[2];
		sprintf(cc,"%d",n);
		string str0=".hdf";
		string str3=".";
		string strout=outlst+str3+cc+str0;
		i2d_emdata_FFT_cpy->ap2ri();
		EMData *i2d_emdata_out=i2d_emdata_FFT_cpy->doIFT();
		i2d_emdata_out->writeImage(strout.c_str(),-1);
		fprintf(g,"%d\t%s\tdefocus=%f\tdfdiff=%f\tdfang=%f\n",0,strout.c_str(),(-1)*defocus,dfdiff,dfang);
		/////////////////////////////////////////
		i2d_emdata_FFT_cpy->gimmeFFT();
		delete i2d_emdata_FFT_cpy;
		delete i2d_emdata_out;
		i2d_emdata->gimmeFFT();
		delete i2d_emdata;
		/*free(rd);
		free(ra);
		rd=NULL;
		ra=NULL;*/
	}
	fclose(g);
	return 0;
}
	