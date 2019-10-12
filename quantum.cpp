#include <fftw3.h>
#include <iostream>
#include <cmath>
#include <math.h>
#include <fstream>
#include <cstdlib>  
#include <graphics.h>
#include <complex.h>

using namespace std;

int main(void){
	
	//int gd = DETECT,gm;
	int gd = 6;
	int gm = 3;
    initgraph(&gd,&gm,NULL);
    //initgraph(IBM8514,IBM8514,NULL);
    
	double m = 1.0;
	double h = 1.0;
	double dx = 12.0;//szerokosc
	double t = 1000.;
	//double Nt = 1000;
	double dt = 10.0;
	int N = 800;
	double p = 0;
	double dx2 = 1.0;//krok
	fftw_complex psi0[N];
	fftw_complex psi1[N];
	fftw_complex psi2[N];
	double w0[N];
	double w1[N];
	double w2[N];
	double w3[N];
	double w[N];
	double Z2;
	double xp = 200;
	fftw_plan plan;
	fftw_plan plan2;
	float pm[N];
	double Ak[N];
	
	//setcolor(BLACK);
	outtextxy(120,580,"p=0, V=0, klakta nr:90");
	//setbkcolor(WHITE);
	
	for(int i=0 ; i<N ; i++){
		psi0[i][0] = pow((2 * M_PI * dx * dx),-0.25) * exp(-1 * (pow(i*dx2-xp,2)) / (4 * dx * dx)) * cos(p*dx2*i/h); //rzeczywiste
		psi0[i][1] = pow((2 * M_PI * dx * dx),-0.25) * exp(-1 * (pow(i*dx2-xp,2)) / (4 * dx * dx)) * sin(p*dx2*i/h); //urojone
	}
		
	for(int i=0 ; i<N ; i++){ 
		w0[i] = 0;
	}
	for(int i=0 ; i<N ; i++){ 
		w1[i] = 0.5*(6 * pow((double)i/N,5) - 15 * pow((double)i/N,4) + 10 * pow((double)i/N,3));
		Z2 = w1[i];
	}
	for(int i=0 ; i<N ; i++){ 
		w2[i] = Z2;
	}
	for(int i=0 ; i<N ; i++){
		w3[i] = -0.5*(6 * pow((double)i/N,5) - 15 * pow((double)i/N,4) + 10 * pow((double)i/N,3)) + w2[i];
	}
	
	
		for(int i=0 ; i<N ; i++){
			if(i<=N*0.5)
				pm[i] = (2.*M_PI * i)/(N * dx2);
			else
				pm[i] = 2.*M_PI*(i - N)/(N * dx2);//12.18 | k=pm
			}
	
	for(int a = 0 ; a*dt<t ; a++){
		
		for(int i=0 ; i<N ; i++){ 
			//w[i]=w0[i];
		if(((double)i<=0.5*(double)N) || (double)i>=0.8*(double)N){
			//putpixel(i,w0[i],14);
			w[i]=w0[i];
		}
		if((double)i>=0.5*(double)N && (double)i<=0.6*(double)N){
			//putpixel(i,w1[10*(i-(int)(0.5*(double)N))],14);
			w[i]=w1[10*(i-(int)(0.5*(double)N))];
		}
		if((double)i>=0.6*(double)N && (double)i<=0.7*(double)N){
			//putpixel(i,w2[10*(i-(int)(0.6*(double)N))],14);
			w[i]=w2[10*(i-(int)(0.6*(double)N))];	
		}
		if((double)i>=0.7*(double)N && (double)i<=0.8*(double)N){
			//putpixel(i,w3[10*(i-(int)(0.7*(double)N))],14);
			w[i]=w3[10*(i-(int)(0.7*(double)N))];
		}
		putpixel(i,-150*w[i]+560,BLUE);
		
		}
		
		for(int i = 0 ; i<N ; i++)
		{
			psi1[i][0] = psi0[i][0]*cos((w[i]-w0[i])/2*dt/h) + psi0[i][1]*sin((w[i]-w0[i])/2*dt/h);
			psi1[i][1] = -psi0[i][0]*sin((w[i]-w0[i])/2*dt/h) + psi0[i][1]*cos((w[i]-w0[i])/2*dt/h);
		}	
		
		plan = fftw_plan_dft_1d(N, psi1, psi2, FFTW_FORWARD, FFTW_ESTIMATE);
		fftw_execute(plan);
		fftw_destroy_plan(plan);
		fftw_cleanup();	
	
		for(int i=0 ; i<N ; i++){
			psi1[i][0] = psi2[i][0]*cos(h*dt/(2*m)*pm[i]*pm[i]) + psi2[i][1]*sin(h*dt/(2*m)*pm[i]*pm[i]); //liczby z 12.18
			psi1[i][1] = -psi2[i][0]*sin(h*dt/(2*m)*pm[i]*pm[i]) + psi2[i][1]*cos(h*dt/(2*m)*pm[i]*pm[i]);
		}
		
		plan2 = fftw_plan_dft_1d(N, psi1, psi2, FFTW_BACKWARD, FFTW_ESTIMATE);
		fftw_execute(plan2);
		fftw_destroy_plan(plan2);
		fftw_cleanup();
		
		for(int i = 0 ; i<N ; i++){
			putpixel(i,-3000*(pow(psi0[i][0],2)+pow(psi0[i][1],2))+560,BLACK);
		}
		
		for(int i = 0 ; i<N ; i++){
			psi2[i][0] = psi2[i][0]/N;
			psi2[i][1] = psi2[i][1]/N;
			psi0[i][0] = psi2[i][0]*cos((w[i]-w0[i])/2*dt/h) + psi2[i][1]*sin((w[i]-w0[i])/2*dt/h);
			psi0[i][1] = -psi2[i][0]*sin((w[i]-w0[i])/2*dt/h) + psi2[i][1]*cos((w[i]-w0[i])/2*dt/h);
			putpixel(i,-3000*(pow(psi0[i][0],2)+pow(psi0[i][1],2))+560,GREEN);
			Ak[i] = pow((2*M_PI*(dx*dx+h*h*a*dt*a*dt/(4*m*m*dx*dx))),-0.5)*exp(-(i*dx2-xp)*(i*dx2-xp)/(2*(dx*dx+h*h*a*dt*a*dt/(4*m*m*dx*dx))));
			//putpixel(i,-2000*Ak[i]+590,RED);
			//cout<<-2000*(pow(psi0[i][0],2)+pow(psi0[i][1],2))+470<<endl;
		
		}
		cout<<a<<endl;
		if(a==90) break;
	}
	
	
	getch();
	closegraph();
	
	return 0;

}
