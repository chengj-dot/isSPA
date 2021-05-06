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

inline float asub(float x,float y) {
float r;

r=x-y;
if (r>PI) r=r-2*PI;
if (r<(-1)*PI) r=2*PI+r;

return r;
}

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
float *Calculate_cc (EMData *map1, EMData *map2, float rot, float ds, float d_m){
	float score=0,phy=0,r=0,phy1=0,delta_phy=0,sd=0,pres=0,pres_down=0,pres_up=0,ss,ag,v,peak=0;
	int x, total_n=0,nx,ny,nz,Ny;
	Ny=d_m;
	float re,im,cx=0,cy=0;
	float centerx=0,centery=0;
	EMData *cf,*cf_ift,*f1_fft;
	cf=new EMData();cf_ift=new EMData();
	f1_fft=new EMData();
	f1_fft=map1->copy();
	f1_fft->ap2ri();
	cf=map2->copy();
//	cf=cf->copy(0,0);
	cf->ap2ri();
	nx=cf->xSize();
	ny=cf->ySize();
	nz=cf->zSize();
	EMData *cf_ft=new EMData();
	cf_ft=cf->copy();
	float *rdata2=f1_fft->getDataRO();
	float *rdata1=cf_ft->getData();
	for (int i=0; i<nx*ny*nz; i+=2) {
		re=(rdata1[i]*rdata2[i]+rdata1[i+1]*rdata2[i+1]);
		im=(rdata1[i+1]*rdata2[i]-rdata1[i]*rdata2[i+1]);
		rdata1[i]=re;
		rdata1[i+1]=im;
	}
	f1_fft->doneData();		
	cf_ft->doneData();           // ok done with the original data
	cf_ft->toCorner();
	EMData *cf_ft_ift=new EMData();
	cf_ft_ift=cf_ft->doIFT();
	cf_ft->gimmeFFT();
	delete cf_ft;
	float *odata=cf_ft_ift->getData();
	for(int j=Ny/4;j<ny-Ny/4;j++){
		for(int i=Ny/4;i<ny-Ny/4;i++){
			if(odata[i+j*ny] > peak) { peak=odata[i+j*ny]; cx=i; cy=j;}
		}
	}
	//peak=peak/(ny*ny)/(Ny*Ny);
	/////////////////////////////////output CCG values below
	centerx=(cx-ny/2)*cos(rot*PI/180)+(cy-ny/2)*sin(rot*PI/180)+ny/2;
	centery=(cy-ny/2)*cos(rot*PI/180)-(cx-ny/2)*sin(rot*PI/180)+ny/2;cout<<"centerx="<<centerx<<"  centery="<<centery<<"  peak="<<peak<<endl;
	////////////////////////////
	float ra=0,rb=0,rc=0,var=0;
	for(int j=0;j<ny;j++){
		for(int i=0;i<ny;i++){
			if(i==cx and j==cy){
				continue;
			}
			ra += odata[i+j*ny];
			rb += SQR(odata[i+j*ny]);
			rc += 1.0;
		}
	}
	var=rb/rc-SQR(ra/rc);
	sd=sqrt(var);
	score=peak/sd;
	delete f1_fft;
	delete cf_ift;
	cf->gimmeFFT();
	delete cf;
	//if(centerx<Ny/3 || centery<Ny/3 || centerx>(ny-Ny/3) || centery>(ny-Ny/3)) score=0;
	cf_ft_ift->doneData();
	float cc[3]={0},*p;
	cc[0]=centerx;
	cc[1]=centery;
	cc[2]=score;
	p=cc;
	delete cf_ft_ift;
	return p;
}
int main(int argc, char *argv[]){
	char inlst[512],outlst[512],temp2d[512],eulerf[512],snr[512],fscw[512];
	string option;
	float apix,edge_half_width=4,euler1=0,euler2=0,euler3=0,centerx=0,centery=0,x=0,y=0,defocus=0,dfdiff=0,dfang=90;
	float energy=300,cs=2.7,ampconst=0.07,ds,lambda,gamma,ss,v,dfu,dfv,ag;
	float *parm,cc_score=0,lowres=30,highres=8,thres,phi_step;
	float a,b,b2,bfactor,bfactor2,bfactor3,k=0,kk=3,Bf;
	int count=0,j,x1,y1,n,first,last,ny_2d,nx2,ny2,nz,nx,ny,yy,d_m;
	int sizex,sizey,tmp_nx;
	float theta0,phi0,delta_theta,delta_phi,r,d;
	float x_norm,y_norm,z_norm,xx0,yy0,zz0;
	int current=0;
	const char *filespec;
	if(argc==1){
		printf("This program detect targets with orientations and tanslations.\n");
		printf("--input            input filtered images lstfile with ctf information\n");
		printf("--angpix            input pixel size in angstroms\n");
		printf("--template         input 2D projections templates in .hdf format\n");
		printf("--eulerfile        euler file with euler values\n");
		printf("--phistep         inplane rotation sampling\n");
		printf("--weight           optimized weight with SSNR parmeters.\n");
		printf("--kk               overlapping density parameter, default is 3.\n");
		printf("--first            the first image to process.\n");
		printf("--last             the last image to process.\n");
		printf("--energy           accerlerating voltage in kV.\n");
		printf("--cs               spherical aberration in um.\n");
		printf("--Highres          high resolution cut \n");
		printf("--Lowres           low resolution cut\n");
		printf("--diameter         target diameter in pixel\n");
		printf("--threshold        cc threshold value, only output LOCs beyond this value\n");
		printf("--output           output lstfile filaname\n");
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
        if(option.compare("--template")==0){
			i++;
			strcpy(temp2d, argv[i]);
			continue;
		}		
		if(option.compare("--eulerfile")==0){
			i++;
			strcpy(eulerf, argv[i]);
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
		if(option.compare("--phistep")==0){
			i++;
			phi_step=atof(argv[i]);
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
		if(option.compare("--diameter")==0){
			i++;
			d_m=atof(argv[i]);
			continue;
		}
		if(option.compare("--threshold")==0){
			i++;
			thres=atof(argv[i]);
			continue;
		}
		printf("Undefined option: %s .Abort.\n",argv[i]);
	}
	std::ifstream eulerfile(eulerf);
    std::string euler_value[20000];
    int n1=0;
	while(!eulerfile.eof()){	
	    getline(eulerfile,euler_value[n1],'\n');
	    n1++;
	}
	n1=n1-1;
	std::ifstream snr_parm(snr);
	string str2,str0;
	while(!snr_parm.eof()){
	    getline(snr_parm, str2,'\n');
		sscanf(str2.c_str(),"%f\t%f\t%f\t%f\t%f\t%f",&a,&b,&b2,&bfactor,&bfactor2,&bfactor3);//cout<<a<<"  "<<b<<"  "<<b2<<"  "<<bfactor<<endl;
	}
	FILE *fp=fopen(outlst,"wb");
	for(n=first;n<last;n++){
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
		//i2d_emdata->writeImage("debug.hdf",-1);
		nx2=i2d_emdata->xSize();
		ny2=i2d_emdata->ySize();
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
			else if(strncmp(pairs[i].c_str(),"dfang=",6)==0){
				if(sscanf(pairs[i].substr(6).c_str(),"%f",&dfang)==1){
				}
			}
			else {
			}
	    }
		dfu=defocus+dfdiff; //defocus is minus, so abs(dfu) < abs(dfv)
		dfv=defocus-dfdiff;
		lambda=12.2639/sqrt(energy*1000.0+0.97845*energy*energy);
		ds=1/(apix*ny2);
		for(int J=0;J<n1;J++){
        	sscanf(euler_value[J].c_str(),"%f %f",&euler1,&euler2);
			EMData *infile_temp=new EMData();
	        infile_temp->readImage(temp2d,J);
			ny=infile_temp->ySize();
			///////////////////////////////////////////////padding the template
			float ra=0,rb=0;
			EMData *infile_temp_clip=new EMData();
			infile_temp_clip=infile_temp->clip((ny-ny2)/2,(ny-ny2)/2,ny2,ny2);
			infile_temp_clip->edgeNormalize(1);
			/////////////////////////////////////////////////////////apply whitening filter
			float ra2[5000]={0},rb2[5000]={0};
			float Ncurve[5000]={0},signal[5000]={0};
			EMData *infile_FFT=new EMData();
			infile_FFT=infile_temp_clip->doFFT();
			infile_temp_clip->gimmeFFT();
			delete infile_temp_clip;
			infile_FFT->ri2ap();
			EMData *infile_cpy=new EMData();
		    infile_cpy=infile_FFT->copy();
		    float *infile_cpy_data=infile_cpy->getData();
		    float *infile_data=infile_FFT->getDataRO();
			for(int jj=0, y2=0;jj<ny2;jj++){
			    for(int ii=0;ii<ny2+2;ii+=2,y2+=2){
				    r = sqrt(SQR(ii / 2) + SQR(jj - ny2 / 2));
				    int xx = (floor(r + 0.5) - 1);
				    if (xx<ny2 / 2 && xx >= 0) {
					    ra2[xx]+=SQR(infile_data[y2]);
					    rb2[xx] += 1.0;
				    }
			    }
			}
            for (int jj=0, y2=0;jj<ny2;jj++){
			    for(int ii=0;ii<ny2+2;ii+=2,y2+=2){
				     r = sqrt(SQR(ii / 2) + SQR(jj - ny2 / 2));
				     int xx = floor(r + 0.5) - 1;
					 ss = r*ds;
				     if (x<ny2 / 2 && xx >= 0) {
					     float fb_infile=ra2[xx]/rb2[xx];
					     v=CTF_AST(ii/2,jj,ny2,ds,dfu,dfv,dfdiff,dfang,lambda,cs,ampconst,2);
		                 infile_cpy_data[y2]=infile_cpy_data[y2]/sqrt(fb_infile);
				     }
		        }
	        } 
		    infile_FFT->doneData();
		    infile_FFT->update();
			delete infile_FFT;
			infile_cpy->doneData();
		    infile_cpy->update();
			///////////////////////////////////////////////////////////////////////////////apply mask
			infile_cpy->ap2ri();
			EMData *infile_FFT_ift=new EMData();
			infile_FFT_ift=infile_cpy->doIFT();
			infile_cpy->gimmeFFT();
			delete infile_cpy;
			float *infile_filtered_data=infile_FFT_ift->getData();
			for (int jj=0; jj<ny2; jj++){
				for (int ii=0; ii<ny2; ii++){
					r=sqrt(SQR(ii-ny2/2)+SQR(jj-ny2/2));
					if( r>(d_m/2+2*edge_half_width)){
			             infile_filtered_data[ii+jj*ny2]=0;
			        }
			        else if (r<d_m/2){
				         infile_filtered_data[ii+jj*ny2]=infile_filtered_data[ii+jj*ny2];
				            }
			        else {
		                 d=0.5*cos(PI*(r-d_m/2)/(2*edge_half_width))+0.5;
			             infile_filtered_data[ii+jj*ny2]=infile_filtered_data[ii+jj*ny2]*d;
			        }
				}
			}
			infile_FFT_ift->doneData();
			infile_FFT_ift->update();//infile_FFT_ift->writeImage("debug-4.mrc",-1);
			///////////////////////////////////////////do lowpass below
			EMData *infile_filtered_norm=new EMData();
			infile_filtered_norm=infile_FFT_ift->doFFT();
			infile_FFT_ift->gimmeFFT();
			infile_filtered_norm->ri2ap();
			float *infile_filtered_norm_data=infile_filtered_norm->getData();
			for (int jj=0, y2=0;jj<ny2;jj++){
			    for(int ii=0;ii<ny2+2;ii+=2,y2+=2){
				     r = sqrt(SQR(ii / 2) + SQR(jj - ny2 / 2));
				     x = floor(r + 0.5) - 1;
				     if (x<ny2*apix/highres && x >= 0) {
						 infile_filtered_norm_data[y2]=infile_filtered_norm_data[y2];
					 }
					else if(x>=ny2*apix/highres && x<ny2*apix/highres+8){
						 infile_filtered_norm_data[y2]=infile_filtered_norm_data[y2]*(0.5*cos(PI*(x-ny2*apix/highres)/(2*8))+0.5);
					}
					else if(x>=(ny2*apix/lowres-8) && x<ny2*apix/lowres && x>=0){
						 infile_filtered_norm_data[y2]=infile_filtered_norm_data[y2]*(0.5*cos(PI*(ny2*apix/lowres-x)/(2*8))+0.5);
					}
					else
						 infile_filtered_norm_data[y2]=0;
				}
			}
			float infile_mean2=0,infile_sigma2=0,total_n=0;			
			//////////////////////////////////////////////////////////apply my weighting function
			float weight_norm=0,count_w;
			for(int jj=0,y2=0;jj<ny2;jj++){
			    for(int ii=0;ii<ny2+2;ii+=2,y2+=2){
				    r = sqrt(SQR(ii/2) + SQR(jj-ny2/2));
				    int xx=floor(r+0.5)-1;
				    ss=r*ds;
				    if(xx<ny2/2 && xx>=0){
					    v=CTF_AST(ii/2,jj,ny2,ds,dfu,dfv,dfdiff,dfang,lambda,cs,ampconst,2);
					    signal[xx]=(exp(bfactor*ss*ss+bfactor2*ss+bfactor3))/(kk+1);
					    Ncurve[xx]=exp(a*ss*ss+b*ss+b2)/signal[xx];
						//euler_w[x]=1.68118*ss;
					    infile_filtered_norm_data[y2]=infile_filtered_norm_data[y2]*v*sqrt(1/(Ncurve[xx]+kk*v*v));
				    }
			    }
		    }
			//////////////////////////////////do normalization
			for (int jj=0, y2=0;jj<ny2;jj++){
			    for(int ii=0;ii<ny2+2;ii+=2,y2+=2){
				     r = sqrt(SQR(ii / 2) + SQR(jj - ny2 / 2));
				     int xx = floor(r + 0.5) - 1;
					 ss=r*ds;
					 
				     if (xx<ny2/2 && xx >= 0) {
						 infile_sigma2+=SQR(infile_filtered_norm_data[y2]);
						 total_n+=1;
					 }
				}
			}
			infile_mean2=infile_sigma2/SQR(total_n);
			for (int jj=0, y2=0;jj<ny2;jj++){
			    for(int ii=0;ii<ny2+2;ii+=2,y2+=2){
					infile_filtered_norm_data[y2]=infile_filtered_norm_data[y2]/sqrt(infile_mean2);
				}
			}
			infile_filtered_norm->doneData();
			infile_filtered_norm->update();
			///////////////////////////////////////////////////////////////////////////////////////////
			for(euler3=0.0;euler3<360.0;euler3+=phi_step){	
			    if(1){
					EMData *i2d_emdata_cpy=new EMData();
				    i2d_emdata_cpy=i2d_emdata->copy();
					i2d_emdata_cpy->setRAlign(euler3*PI/180);
		            i2d_emdata_cpy->rotateAndTranslate();
		            i2d_emdata_cpy->doneData();
		            i2d_emdata_cpy->update();
					EMData *i2d_filtered_norm=new EMData();
					i2d_filtered_norm=i2d_emdata_cpy->doFFT();
					i2d_emdata_cpy->gimmeFFT();
					i2d_filtered_norm->ri2ap();
					///////////////////////////////////////////////////
					float *cc_parm;
					cc_parm=Calculate_cc(infile_filtered_norm,i2d_filtered_norm,euler3,ds,d_m);
					centerx=*cc_parm;
					centery=*(cc_parm+1);
					cc_score=*(cc_parm+2);
					//cout<<J<<"\t"<<centerx<<"\t"<<centery<<"\t"<<euler1<<"\t"<<euler2<<"\t"<<euler3<<"\t"<<cc_score<<endl;
					if((cc_score)>thres){
						fprintf(fp,"%d\t%s\tdefocus=%f\tdfdiff=%f\tdfang=%f\teuler=%f,%f,%f\tcenter=%f,%f\tscore=%f\n",nn,t,(-1)*defocus,
						dfdiff,dfang,euler1,euler2,euler3,centerx,centery,cc_score);
					}		   
	                delete i2d_emdata_cpy;
					delete i2d_filtered_norm;
				}
			}
			delete infile_filtered_norm;
			delete infile_FFT_ift;
			delete infile_temp;
		}
		delete i2d_emdata;
	}
	fclose(fp);
	return 0;
}	
