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

void display_usage(){
	cout<<"This program calculate cc of two images."<<endl;
	cout<<"-F <input FSC file>"<<endl;
	cout<<"-a -v <apix, accleration voltage>"<<endl;
	cout<<"-H -L  <FSC high threshold and low threshold>"<<endl;
	cout<<"-r <3D reference map>"<<endl;
	cout<<"-i <Input 2d lst>"<<endl;
	cout<<"-f -l <first to last in input 2d lst>"<<endl;
	cout<<"-R -u <mask radius in pixel, -R = outer, -u = soft edge width>"<<endl;
	cout<<"-o -c <output lst, debug parm , 1 for cc, 0 for phase residual>"<<endl;
}
static const char *optString = "F:a:v:H:L:r:i:f:l:R:u:o:c:h?";

inline float asub(float x, float y){
	float r;
	r=x-y;
	if(r>PI) r=r-2*PI;
	if(r<(-1)*PI) r=r+2*PI;
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

float calculate_cc(EMData *map1, EMData *map2, int ny, float dfu, float dfv, float dfdiff, float dfang, float fsc_y[], float highres, float lowres, float apix, float cs, 
float lambda, float ampconst, int debug){
	EMData *map1_fft=new EMData();
	EMData *map2_fft=new EMData();
	map1_fft=map1->doFFT();
	map2_fft=map2->doFFT();
	map1->gimmeFFT();
	map2->gimmeFFT();
	//delete map1;
	//delete map2;
	float ds=1/(apix*ny);
/*	map1_fft->ri2ap();
	map2_fft->ri2ap();
	float *d1=map1_fft->getData();
	float *d2=map2_fft->getData();
	float sum1[2000],count[2000],sum2[2000];
	for(int j=0, y=0; j<ny; j++){
		for(int i=0; i<(ny+2); i+=2, y+=2){
			float r=sqrt(SQR(i/2)+SQR(j-ny/2));
			int x=floor(r+.5)-1;
			if (x<ny/2 && x>=0) {
				sum1[x]+=SQR(d1[y]);
				sum2[x]+=SQR(d2[y]);
				count[x]+=1.0;
			}
		}
	}
	for(int j=0, y=0; j<ny; j++){
		for(int i=0; i<(ny+2); i+=2, y+=2){
			float r=sqrt(SQR(i/2)+SQR(j-ny/2));
			int x=floor(r+.5)-1;
			if (x<ny/2 && x>=0) {
				float fa=sum1[x]/count[x];
				float fb=sum2[x]/count[x];
				d1[y]=d1[y]/sqrt(fa);
				d2[y]=d2[y]/sqrt(fb);
			}
		}
	}
	map1_fft->doneData();
	map2_fft->doneData();
	map1_fft->update();
	map2_fft->update();*/
	if(debug==1){
	    map1_fft->ap2ri();
	    map2_fft->ap2ri();
	}
	if(debug==0){
		map1_fft->ri2ap();
		map2_fft->ri2ap();
	}
	float *data1=map1_fft->getDataRO();
	float *data2=map2_fft->getDataRO();
	float score=0,norm1=0,norm2=0,norm3=0;
	/*float *rd = (float *)malloc(ny / 2 * sizeof(float));
	float *ra = (float *)malloc(ny / 2 * sizeof(float));
	float *rb = (float *)malloc(ny / 2 * sizeof(float));
	float *rc = (float *)malloc(ny / 2 * sizeof(float));
	memset(rd, 0, ny / 2 * sizeof(float));
	memset(ra, 0, ny / 2 * sizeof(float));
	memset(rb, 0, ny / 2 * sizeof(float));
	memset(rc, 0, ny / 2 * sizeof(float));*/
	for(int j=0, y=0; j<ny; j++){
		for(int i=0; i<(ny+2); i+=2, y+=2){
			float r=sqrt(SQR(i/2)+SQR(j-ny/2));
			int x=floor(r+.5)-1;
			if (x<(ny*apix/highres) && x>(ny*apix/lowres)) { 
			    float v=CTF_AST(i/2,j,ny,ds,dfu,dfv,dfdiff,dfang,lambda,cs,ampconst,2);
				if(debug==1){
				    score+=(data1[y]*data2[y]+data1[y+1]*data2[y+1])*fsc_y[x]*v;
				    norm1+=data1[y]*data1[y]+data1[y+1]*data1[y+1]; 
				    norm2+=data2[y]*data2[y]+data2[y+1]*data2[y+1];
				    norm3+=SQR(fsc_y[x]*v);
				}
				if(debug==0){
					float delta=asub(data1[y+1],data2[y+1]);
					score+=fsc_y[x]*v*cos(delta);
					norm1+=(fsc_y[x]*v);
				}
			}
		}
	}
	map1_fft->doneData();
	map2_fft->doneData();
	delete map1_fft;
	delete map2_fft;
	if(debug==1) score=score/sqrt(norm1*norm2*norm3);
	if(debug==0) score=score/(norm1);
	/*free(ra);
	free(rb);
	free(rc);
	free(rd);
	ra=NULL;
	rb=NULL;
	rc=NULL;
	rd=NULL;*/
	return score;
}

int main(int argc, char *argv[]){
	string input3d,inlstname,fscfile,outname;
	float highres,lowres,apix,radius,softedge;
	float cs=2.7,lambda,dfdiff,defocus,dfu,dfv,dfang,ds,euler1,euler2,euler3,centerx,centery,ampconst=0.1,energy;
	int nx_3d,ny_3d,nz_3d,nx,ny,first,last;
	int n=0;
	int debug=0;
	float score=0;
	const char *filespec;
	string tmp;
    int opt;
	opt = getopt(argc,argv,optString);
	while(opt!=-1) {
		switch(opt) {
			case 'F':
				fscfile=optarg;
				break;
			case 'a':
				tmp=optarg;
				apix=atof(tmp.c_str());
				break;
			case 'v':
			    tmp=optarg;
				energy=atof(tmp.c_str());
				break;
			case 'H':
				tmp=optarg;
				highres=atof(tmp.c_str());
				break;
			case 'L':
				tmp=optarg;
				lowres=atof(tmp.c_str());
				break;
			case 'r':
				input3d=optarg;
				break;
			case 'i':
				inlstname=optarg;
				break;
			case 'f':
				tmp=optarg;
				first=atoi(tmp.c_str());
				break;
			case 'l':
				tmp=optarg;
				last=atoi(tmp.c_str());
				break;
			case 'R':
				tmp=optarg;
				radius=atof(tmp.c_str());
				break;
			case 'u':	
				tmp=optarg;
				softedge=atof(tmp.c_str());
				break;
			case 'o':
				outname=optarg;
				break;
			case 'c':
				tmp=optarg;
				debug=atoi(tmp.c_str());
				break;
			case '?':
                display_usage();
				return 0;
                break;
			case 'h':
				display_usage();
				return 0;
                break;
		}
		opt = getopt( argc, argv, optString );
	}
	std::ifstream fsc3d(fscfile.c_str());
	float fsc_x[1000],fsc_y[1000];
	int count_fsc=0;
    string str0;
	while(1){
		getline(fsc3d, str0,'\n');
		if(fsc3d.eof()){
			break;
			str0="";
		}
		sscanf(str0.c_str(),"%f\t%f",&fsc_x[count_fsc],&fsc_y[count_fsc]);
		count_fsc++;
		str0="";
	}
	ofstream fpout;
	fpout.open(outname.c_str());
	EMData *em3d=new EMData();
	em3d->readImage(input3d.c_str(),-1);
	nx_3d=em3d->xSize();
	ny_3d=em3d->ySize();
	nz_3d=em3d->zSize();
	for(n=first;n<last;n++){
		//  Read lst and return dfdiff, dfang, euler angle, center
		filespec=inlstname.c_str();
		char s[MAXPATHLEN+10],t[MAXPATHLEN+1],u[256],*c;
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
		for(int i=0; i<pairs.size(); i++){
			if(strncmp(pairs[i].c_str(),"euler=",6)==0) {
				if(sscanf(pairs[i].substr(6).c_str(),"%f,%f,%f",&euler1,&euler2,&euler3)==3){
				}
			}
			else if(strncmp(pairs[i].c_str(),"center=",7)==0) {
				if(sscanf(pairs[i].substr(7).c_str(),"%f,%f",&centerx,&centery)==2){
				}
			}
			else if(strncmp(pairs[i].c_str(),"defocus=",8)==0) {
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
		EMData *em2d_image=new EMData();
		em2d_image->readImage(t,nn);
		nx=em2d_image->xSize();
		ny=em2d_image->ySize();
		if(nx != nx_3d || ny != ny_3d) {
			cout<<"2d image is not the same size as 3d map!!"<<endl;
			break;
		}
		dfu=defocus+dfdiff;
		dfv=defocus-dfdiff;
		lambda=12.2639/sqrt(energy*1000+0.97845*energy*energy);
		ds=1/(apix*ny);
		/////// translate 2d image to the center
		EMData *em2d_cpy=new EMData();
		em2d_cpy=em2d_image->copy();
		float *data1=em2d_image->getDataRO();
		float *do_i2d_translate=em2d_cpy->getData();
		for (int j=0,y=-ny/2.0; j<ny; j++,y+=1.0) {
		    for (int i=0,x=-nx/2.0; i<nx; i++,x+=1.0) {
			    float x2=x+centerx;
			    float y2=y+centery;
			    if(x2<0)
				    x2=nx+x2;
			    if(x2>nx)
				    x2=x2-nx;
			    if(y2<0)
				    y2=ny+y2;
			    if(y2>ny)
				    y2=y2-ny;
//				if (x2<0||x2>nx2-1.0||y2<0||y2>ny2-1.0) {do_i2d_translate[i+j*nx]=0;continue;}
			    int ii=floor(x2);
			    int jj=floor(y2);
			    int k0=ii+jj*nx;
			    int k1=k0+1;
			    int k2=k0+nx+1;
			    int k3=k0+nx;
			    if (ii==nx-1) { k1--; k2--; }
			    if (jj==ny-1) { k2-=nx; k3-=nx; }
			    float t0=(x2-(float)ii);
			    float u0=(y2-(float)jj);
			    float tt=1-t0;
			    float uu=1-u0;
			    float p0=data1[k0]*tt*uu;
			    float p1=data1[k1]*t0*uu;
			    float p3=data1[k3]*tt*u0;
			    float p2=data1[k2]*t0*u0;
			    do_i2d_translate[i+j*nx]=p0+p1+p2+p3;
		    }
	    }
		em2d_cpy->doneData();
		em2d_cpy->update();
		em2d_image->doneData();
		delete em2d_image;
		//////////// project 3d map to projection in a special orientation
		EMData *em3d_cpy=new EMData();
		em3d_cpy=em3d->copy();
		em3d_cpy->setTAlign(0,0,0);
	    em3d_cpy->setRAlign(euler1*PI/180,euler2*PI/180,euler3*PI/180);
	    em3d_cpy->rotateAndTranslate();
		float *em3d_data=em3d_cpy->getDataRO();
		EMData *proj=new EMData();
		proj->setSize(nx,ny,1);
		float *proj_data=proj->getData();
		for(int i=0; i<nx*ny; i++) proj_data[i]=0;
		for(int k=0; k<nz_3d; k++){
			for(int j=0; j<ny_3d; j++){
				for(int i=0; i<nx_3d; i++){
					proj_data[i+j*nx_3d]=proj_data[i+j*nx_3d]+em3d_data[i+j*nx_3d+k*nx_3d*ny_3d];
				}
			}
		}
		proj->doneData();
		proj->update();
        em3d_cpy->doneData();
		delete em3d_cpy;
		//////////////////apply ctf flipphase to 2d images
		EMData *filter=new EMData();
		filter->setSize(nx+2,ny,1);
		filter->setComplex(1);
		filter->setRI(1);
		filter->one();
		float *filter_data=filter->getData();
		for(int j=0; j<ny; j++){
			for(int i=0; i<(nx+2)/2; i++){
				float ctf=CTF_AST(i,j,ny,ds,dfu,dfv,dfdiff,dfang,lambda,cs,ampconst,0);
				filter_data[i*2+j*(nx+2)]*=ctf;
				filter_data[i*2+j*(nx+2)+1]=filter_data[i*2+j*(nx+2)];
			}
		}
		filter->doneData();
		filter->update();
		EMData *em2d_fft=new EMData();
		em2d_fft=em2d_cpy->doFFT();
		em2d_fft->ap2ri();
		em2d_cpy->gimmeFFT();
		delete em2d_cpy;
		em2d_fft->mult(filter);
		delete filter;
		EMData *em2d_flip=new EMData();
		em2d_flip=em2d_fft->doIFT();
		em2d_fft->gimmeFFT();
		delete em2d_fft;
		///////////////// apply circular mask in real space to 2d images
		float *em2d_flip_data=em2d_flip->getData();
		for(int j=0; j<ny; j++){
			for(int i=0; i<nx; i++){
				float r=sqrt(SQR(i-nx/2)+SQR(j-ny/2));
				if(r>(radius+2*softedge)){
					em2d_flip_data[i+j*nx]=0;
				}
				else if(r<radius){
					em2d_flip_data[i+j*nx]=em2d_flip_data[i+j*nx];
				}
				else{
					float d=0.5*cos(PI*(r-radius)/(2*softedge))+0.5;
					em2d_flip_data[i+j*nx]=em2d_flip_data[i+j*nx]*d;
				}
			}
		}
		em2d_flip->doneData();
		em2d_flip->update();
		////////////////////////calcuate weighted cc
		score=calculate_cc(em2d_flip, proj, ny, dfu, dfv, dfdiff, dfang, fsc_y, highres, lowres, apix, cs, lambda, ampconst, debug);
		fpout<<nn<<"\t"<<t<<"\t"<<"defocus="<<(-1)*defocus<<"\t"<<"dfdiff="<<dfdiff<<"\t"<<"dfang="<<dfang<<"\t"<<"euler="<<euler1<<","<<euler2<<","<<euler3<<"\t"<<"center="<<centerx<<","<<centery<<
			"\t"<<"score="<<score<<endl;
		delete proj;
		delete em2d_flip;
	}
	fpout.close();
	delete em3d;
	return 0;
}
