



//working sound analysis program windows 8.1 version

    //source:
    //https://www.google.ch/patents/US6957241
    //http://www.google.ch/patents/US20020083107
    //https://www.beechwood.eu/fft-implementation-r2-dit-r4-dif-r8-dif/
    //http://www.cmlab.csie.ntu.edu.tw/cml/dsp/training/coding/transform/fft.html
    //http://dsp.stackexchange.com/questions/3481/radix-4-fft-implementation

    //https://community.arm.com/graphics/b/blog/posts/speeding-up-fast-fourier-transform-mixed-radix-on-mobile-arm-mali-gpu-by-means-of-opencl---part-1
    //book: "Cyfrowe przetwarzanie sygnalow" - Tomasz Zielinski it has quick version of radix-2 because it calculates less sin() and cos()

    //Rabiner L.R., Gold B. Theory and application of digital signal processing p 378 mixed radix


//http://mixedradixfastfouriertransformifft.blogspot.com/

//author marcin matysek (r)ewertyn.PL 2marcin56@gmail.com
//open-source

//algorytm mixed radix FFT  radix-7 * radix-3 * radix-2 * radix-4 * radix-11 N = points witch shift fi v 1.0


#include <iostream>
#include "conio.h"
#include <stdlib.h>
#include <math.h>
#include <cmath>
#include <time.h>
#include <complex>
#include <fstream>

using namespace std;

//complex number method:

void fun_fourier_transform_FFT_radix_2_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp6,std::complex<double> z[]);
void fun_fourier_transform_FFT_radix_3_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp6,std::complex<double> z[]);
void fun_fourier_transform_FFT_radix_4_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp6,std::complex<double> z[]);
void fun_fourier_transform_FFT_radix_5_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp6,std::complex<double> z[]);
void fun_fourier_transform_FFT_radix_7_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp99,std::complex<double> z[]);
void fun_fourier_transform_FFT_radix_11_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp99,std::complex<double> z[]);
void fun_fourier_transform_FFT_radix_97_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp99,std::complex<double> z[]);

void fun_fourier_transform_FFT_radix_2_stg_first(std::complex<double> tab[],int stg_first,double fi2,std::complex<double> z[]);
void fun_fourier_transform_FFT_radix_3_stg_first(std::complex<double> tab[],int stg_first,double fi2,std::complex<double> z[]);
void fun_fourier_transform_FFT_radix_4_stg_first(std::complex<double> tab[],int stg_first,double fi2,std::complex<double> z[]);
void fun_fourier_transform_FFT_radix_5_stg_first(std::complex<double> tab[],int stg_first,double fi2,std::complex<double> z_rx5[]);
void fun_fourier_transform_FFT_radix_7_stg_first(std::complex<double> tab[],int stg_first,double fi2,std::complex<double> z[]);
void fun_fourier_transform_FFT_radix_11_stg_first(std::complex<double> tab[],int stg_first,double fi2,std::complex<double> z[]);
void fun_fourier_transform_FFT_radix_97_stg_first(std::complex<double> tab[],int stg_first,double fi2,std::complex<double> z[]);

void fun_inverse_fourier_transform_FFT_radix_2_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp6,std::complex<double> z[]);
void fun_inverse_fourier_transform_FFT_radix_3_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp6,std::complex<double> z[]);
void fun_inverse_fourier_transform_FFT_radix_4_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp6,std::complex<double> z[]);
void fun_inverse_fourier_transform_FFT_radix_5_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp6,std::complex<double> z[]);
void fun_inverse_fourier_transform_FFT_radix_7_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp99,std::complex<double> z[]);
void fun_inverse_fourier_transform_FFT_radix_11_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp99,std::complex<double> z[]);
void fun_inverse_fourier_transform_FFT_radix_97_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp99,std::complex<double> z[]);

void fun_inverse_fourier_transform_FFT_radix_2_stg_first(std::complex<double> tab[],int stg_first,double fi2,std::complex<double> z[]);
void fun_inverse_fourier_transform_FFT_radix_3_stg_first(std::complex<double> tab[],int stg_first,double fi2,std::complex<double> z[]);
void fun_inverse_fourier_transform_FFT_radix_4_stg_first(std::complex<double> tab[],int stg_first,double fi2,std::complex<double> z[]);
void fun_inverse_fourier_transform_FFT_radix_5_stg_first(std::complex<double> tab[],int stg_first,double fi2,std::complex<double> z_rx5[]);
void fun_inverse_fourier_transform_FFT_radix_7_stg_first(std::complex<double> tab[],int stg_first,double fi2,std::complex<double> z[]);
void fun_inverse_fourier_transform_FFT_radix_11_stg_first(std::complex<double> tab[],int stg_first,double fi2,std::complex<double> z[]);
void fun_inverse_fourier_transform_FFT_radix_97_stg_first(std::complex<double> tab[],int stg_first,double fi2,std::complex<double> z[]);





void fun_inverse_table_FFT(int M,std::complex<double> tab[]);
int radix_base(int N,int stg[]);

void fun_fourier_transform_FFT_mixed_radix(int N,std::complex<double> tab[]);
void fun_inverse_fourier_transform_FFT_mixed_radix(int N,std::complex<double> tab[]);
void fun_fourier_transform_DFT(int N,std::complex<double> tab[]);





////////////////////////////////////////////////////
/////////////////////
double fi=0.00;//for first stage FFT and for DFT
/////////////////////
///////////////////////////////////////////////////


int tb[600][3]={};





static double diffclock(clock_t clock2,clock_t clock1)
{
    double diffticks=clock1-clock2;
    double diffms=(diffticks)/(CLOCKS_PER_SEC/1000);
    return diffms;
}






//==============
int main()
{
    int N;


srand( time( NULL ) );

ifstream file1("11.txt");
ofstream file3("1.txt");
ofstream file4("2.txt");
int a1=0,b1=0;

int aaa=0;
for(int i=0;i<100;i++)
{

    file1>>aaa;
    file1>>aaa;
    file1>>tb[i][0];
    file1>>tb[i][1];
    file1>>tb[i][2];
  //  cout<<aaa<<" "<<tb[i][0]<<" "<<tb[i][1]<<" "<<tb[i][2]<<endl;system("pause");
}





////////////////////
///////////////////
    N=97*2;
    //if N==period of signal in table tab[] then resolution = 1 Hz
///////////////////
///////////////////





    double c1=0,c2=0;
    std::complex<double> *tab2=new std::complex<double>[N];
    std::complex<double> *tab3=new std::complex<double>[N];
    std::complex<double> *tab4=new std::complex<double>[N];
    std::complex<double> tmp1=0;

    const double pi=3.141592653589793238462;

    for(int i=0;i<N;i++)
    {
          c1=rand()%1255/4.0;
          c2=rand()%13591/9.0;
          //  c1=255;
          //  c2=255;
          //c1=i;
          //c2=N+i;
          //c1=4*sin(1*2*pi*i/float(N));
          //c2=5435*sin(21*2*pi*i/float(N));
          //c1=c1+33.19*sin(7*2*pi*i/float(N));
          //c2=0;
        tab2[i].real(c1);
        tab2[i].imag(c2);
        tab3[i].real(tab2[i].real());
        tab3[i].imag(tab2[i].imag());
        tab4[i].real(tab2[i].real());
        tab4[i].imag(tab2[i].imag());
    }


    double time2;
    double zmienna=0;

    cout<<"signal 1="<<endl;
    cout<<"tab2[j].real():"<<endl;
    for(int j=0;j<N;j++)
    {
    cout.precision(4);
    cout<<round(tab2[j].real()*1000)/1000<<"  ";
    }
    cout<<endl;
    cout<<endl;
    cout<<"tab2[j].imag():"<<endl;
    for(int j=0;j<N;j++)
    {
    cout.precision(4);
    cout<<round(tab2[j].imag()*1000)/1000<<"  ";
    }
    cout<<endl;
    cout<<endl;

    cout<<"signal 2="<<endl;
    cout<<"tab3[j].real():"<<endl;
    for(int j=0;j<N;j++)
    {
    cout.precision(4);
    cout<<round(tab3[j].real()*1000)/1000<<"  ";
    }
    cout<<endl;
    cout<<endl;
    cout<<"tab3[j].imag():"<<endl;
    for(int j=0;j<N;j++)
    {
    cout.precision(4);
    cout<<round(tab3[j].imag()*1000)/1000<<"  ";
    }
    cout<<endl;
    cout<<endl;
    clock_t start = clock();

      cout<<"calculating first: DFT"<<endl;
    start = clock();
    fun_fourier_transform_DFT(N,tab3);
    time2=diffclock( start, clock() );
    cout<<"time: [s] "<<time2/1000<<endl;
    system("pause");
    //////////////////////////////////////////////////////////
     cout<<"frequency Hz from DFT"<<endl;
    cout<<"tab3[j].real():"<<endl;
    for(int j=0;j<N;j++)
    {
    cout.precision(4);
    cout<<round(tab3[j].real()*1000)/1000<<"  ";
    }
    cout<<endl;
    cout<<endl;//system("pause");
    cout<<"tab3[j].imag():"<<endl;
    for(int j=0;j<N;j++)
    {
    cout.precision(4);
    cout<<round(tab3[j].imag()*1000)/1000<<"  ";
    }
    cout<<endl;
    cout<<endl;//system("pause");
    system("pause");
int mx=0;
int fg=0;

    for(int i=0;i<N;i++)
    {
        tab2[i].real(tab4[i].real());
        tab2[i].imag(tab4[i].imag());
    }

    cout<<"calculating: FFT"<<endl;
    start = clock();
    fun_fourier_transform_FFT_mixed_radix(N,tab2);
    time2=diffclock( start, clock() );
    cout<<"time: [ms] "<<time2<<endl;
    cout<<"inverse table"<<endl;
    start = clock();
    fun_inverse_table_FFT(N,tab2);
    time2=diffclock( start, clock() );
    cout<<"time: [ms] "<<time2<<endl;
    system("pause");
    ///////////////////////////////////////////////////////////

/////////////////////////////////////////////////////
    time2=diffclock( start, clock() );


    int flag=0,fl=0;
       for(int j=0;j<N;j++)
        {
            for(int i=0;i<N;i++)
        {
          //if((fabs(tab3[i].real()-tab2[j].real())<=0.03)&&fabs(tab3[0].real()-tab2[0].real())<0.03)
          if((fabs(tab3[i].real()-tab2[j].real())<=0.03)&&(fabs(tab3[i].imag()-tab2[j].imag())<=0.03))
          {

              flag++;
/*
            for(int k=0;k<N;k++)
              {
                  if(fabs(tab3[k].real()-tab2[j].real())<=0.03)
                  {fl++;}

              }
              if(fl==1){
                flag++;}
                fl=0;
.*/
            cout.precision(4);
            //cout<<"DFT["<<i<<"]==FFT["<<j<<"]= "<<round(tab2[j].real()*1000)/1000<<"  ";//system("pause");
          }
            else {
               // cout<<-1<<" . ";//system("pause");
            }
        }
        }
cout<<"if DFT[i]==FFT[j]"<<endl;
    //cout<<"flag="<<flag<<endl;
    if(flag>=N/2){
cout<<endl<<"flag= "<<flag<<endl;
fg++;


/*
    cout<<"frequency Hz FFT radix-4 radix-2"<<endl;
    for(int j=0;j<N;j++)
    {
    cout.precision(4);
    cout<<round(tab2[j].real()*1000)/1000<<"  ";
    }
    cout<<endl;
    cout<<endl;
    for(int j=0;j<N;j++)
    {
    cout.precision(4);
    cout<<round(tab2[j].imag()*1000)/1000<<"  ";
    }
    cout<<endl;
    cout<<endl;system("pause");

    cout<<"frequency Hz DFT"<<endl;
    for(int j=0;j<N;j++)
    {
    cout.precision(4);
    cout<<round(tab3[j].real()*1000)/1000<<"  ";
    }
    cout<<endl;
    cout<<endl;system("pause");
    for(int j=0;j<N;j++)
    {
    cout.precision(4);
    cout<<round(tab3[j].imag()*1000)/1000<<"  ";
    }
    cout<<endl;
    cout<<endl;
    cout<<"if radix-4 == DFT tab2[j].real(): "<<endl;system("pause");
*/
       for(int j=0;j<N;j++)
        {
            if((j-1)%10000==0){system("pause")  ;}

            for(int i=0;i<N;i++)
        {
          //if((fabs(tab3[i].real()-tab2[j].real())<=0.03)&&fabs(tab3[0].real()-tab2[0].real())<0.03)
         if((fabs(tab3[i].real()-tab2[j].real())<=0.03)&&(fabs(tab3[i].imag()-tab2[j].imag())<=0.03))
          {
              if(i==30){   cout<<" ...   DFT["<<i<<"]==FFT["<<j<<"]= "<<round(tab2[j].real()*1000)/1000<<"  "<<endl;system("pause");}
            cout.precision(4);
            file3<<"DFT["<<i<<"]==FFT["<<j<<"]= "<<round(tab2[j].real()*10000)/10000<<"  "<<endl;;//system("pause");
           cout<<"DFT["<<i<<"]==FFT["<<j<<"]= "<<round(tab2[j].real()*10000)/10000<<"  "<<endl;;//system("pause");
          }
            else {
               // cout<<-1<<" . ";//system("pause");
            }
        }
        }
    system("pause");
    file3<<endl;

               for(int j=0;j<N;j++)
        {
            for(int i=0;i<N;i++)
        {
          if((fabs(tab3[i].real()-tab2[j].real())<=0.03)&&fabs(tab3[0].real()-tab2[0].real())<0.03)
          {
            cout.precision(4);
            file3<<" "<<i;//system("pause");
          }
            else {
               // cout<<-1<<" . ";//system("pause");
            }
        }
        }
        /////////////////////////////////////////
        tmp1=0;
        //cout<<"tmp1"<<tmp1<<endl;
        for(int i=1;i<N;i++)
        {
            tmp1.real(tmp1.real()+round(tab2[i].real()*100000000)/100000000);
            tmp1.imag(tmp1.imag()+round(tab2[i].imag()*100000000)/100000000);
        }
         cout<<"tab[0] "<<tab2[0]<<endl;
         cout<<"Sum "<<tmp1<<endl;
        /////////////////////////////////////////
    file3<<endl<<endl;

    }


    if(flag>mx){mx=flag;}
    flag=0;

if(a1==10000000){cout<<b1<<endl;a1=0;b1++;}
a1++;

  cout<<endl<<"max="<<mx<<endl;    file3.close();file1.close();file4.close();
    system("pause");
    cout<<endl;

    cout<<"calculating: inverse FFT"<<endl;
    start = clock();
    fun_inverse_fourier_transform_FFT_mixed_radix(N,tab2);
    time2=diffclock( start, clock() );
    cout<<"time: [ms] "<<time2<<endl;
    cout<<"inverse table"<<endl;
    start = clock();
    fun_inverse_table_FFT(N,tab2);
    time2=diffclock( start, clock() );
    cout<<"time: [ms] "<<time2<<endl;
    system("pause");

    cout<<endl;
/*
    cout<<"inverse FFT "<<endl;
    for(int j=0;j<N;j++)
    {
    cout.precision(4);
    cout<<round(tab2[j].real()*1000)/1000<<"  ";
    }
    cout<<endl;
    cout<<endl;
    for(int j=0;j<N;j++)
    {
    cout.precision(4);
    cout<<round(tab2[j].imag()*1000)/1000<<"  ";
    }
    cout<<endl;
    cout<<endl;
*/
    flag=0;
           for(int j=0;j<N;j++)
        {
            if((j-1)%10000==0){system("pause")  ;}

            for(int i=0;i<N;i++)
        {
          //if((fabs(tab3[i].real()-tab2[j].real())<=0.03)&&fabs(tab3[0].real()-tab2[0].real())<0.03)
         if((fabs(tab4[i].real()-tab2[j].real())<=0.03)&&(fabs(tab4[i].imag()-tab2[j].imag())<=0.03))
          {
              flag++;
              if(i==30){   cout<<" ...   DFT["<<i<<"]==FFT["<<j<<"]= "<<round(tab2[j].real()*1000)/1000<<"  "<<endl;system("pause");}
            cout.precision(4);
            file3<<"DFT["<<i<<"]==FFT["<<j<<"]= "<<round(tab2[j].real()*10000)/10000<<"  "<<endl;;//system("pause");
           cout<<"DFT["<<i<<"]==FFT["<<j<<"]= "<<round(tab2[j].real()*10000)/10000<<"  "<<endl;;//system("pause");
          }
            else {
               // cout<<-1<<" . ";//system("pause");
            }
        }
        }
    cout<<endl<<"flag= "<<flag<<endl;
    system("pause");




    delete [] tab2;
    delete [] tab3;
    delete [] tab4;
    system("pause");
    return 0;
}

void fun_inverse_bits_radix_4(int N,std::complex<double> tab[])
{
//code by Sidney Burrus
//http://dsp.stackexchange.com/questions/3481/radix-4-fft-implementation
// Radix-4 bit-reverse
}
///////////////
void fun_fourier_transform_DFT(int N,std::complex<double> tab[])
{
    const double pi=3.141592653589793238462;
    std::complex<double> *tab2 = new std::complex<double>[N];    // tab2[]==N
    std::complex<double>  w[1]={{1,1}};


double zmienna1=2*pi/(float)N;
double fi2=fi;

for (int i=0;i<N;i++)
{
    for(int j=0;j<N;j++)
    {
          //complex number method:
          w[0].real(cos(i*j*zmienna1+fi2));
          w[0].imag(-sin(i*j*zmienna1+fi2));
          tab2[i]=tab2[i]+tab[j]*w[0];

    }
}

//new:
    for(int j=0;j<N;j++)
    {
      tab[j].real(tab2[j].real()*2/N);
     tab[j].imag(tab2[j].imag()*2/N);
    }

    delete tab2;
}
//////////////////

void fun_fourier_transform_FFT_mixed_radix(int N,std::complex<double> tab[])
{

//author marcin matysek (r)ewertyn.PL 2marcin56@gmail.com
//open-source
    const double pi=3.141592653589793238462;
    std::complex<double> *tab2 = new std::complex<double>[N];    // tab2[]==N

    double tmp2;
    double tmp3;
    double tmp6;
    double tmp7;
    double tmp8;
    double tmp9;
    double tmp10;
    double tmp11;
    double tmp97;

    double fi2=fi; //shift only for stage nr 1

    int i=0;
    int rx5=5,rx4=4,rx3=3,rx2=2,rx7=7,rx11=11,rx97=97;
    int stg[100]={};
    int nb_stages=0;
    int nb1,nb2,nb3,nb4,nb5_stg_previous,stg_first;

    tmp6=2*pi/(N/1);

    tmp2=2*pi/(2/1);
    tmp3=2*pi/(3/1);
    tmp10=2*pi/(4/1);
    tmp8=2*pi/(5/1);
    tmp7=2*pi/(7/1);
    tmp11=2*pi/(11/1);
    tmp97=2*pi/(97/1);



    std::complex<double>  z_rx2[2]={{1,0}};
    std::complex<double>  z_rx3[3]={{1,0}};
    std::complex<double>  z_rx4[2]={{1,0}};
    std::complex<double>  z_rx5[5]={{1,0}};
    std::complex<double>  z_rx7[7]={{1,0}};
    std::complex<double>  z_rx11[11]={{1,0}};
    std::complex<double>  z_rx97[10000]={{1,0}};

//radix 2 fundament
          z_rx2[0].real(cos(0*tmp2));
          z_rx2[0].imag(-sin(0*tmp2));
          z_rx2[1].real(cos(1*tmp2));
          z_rx2[1].imag(-sin(1*tmp2));
//radix 3 fundament
          z_rx3[0].real(cos(0*tmp3));
          z_rx3[0].imag(-sin(0*tmp3));
          z_rx3[1].real(cos(1*tmp3));
          z_rx3[1].imag(-sin(1*tmp3));
          z_rx3[2].real(cos(2*tmp3));
          z_rx3[2].imag(-sin(2*tmp3));
//radix 4 fundament
          z_rx4[0].real(cos(0*tmp10));
          z_rx4[0].imag(-sin(0*tmp10));
          z_rx4[1].real(cos(1*tmp10));
          z_rx4[1].imag(-sin(1*tmp10));
//radix 5 fundament
          z_rx5[0].real(cos(0*tmp8));
          z_rx5[0].imag(-sin(0*tmp8));
          z_rx5[1].real(cos(1*tmp8));
          z_rx5[1].imag(-sin(1*tmp8));
          z_rx5[2].real(cos(2*tmp8));
          z_rx5[2].imag(-sin(2*tmp8));
          z_rx5[3].real(cos(3*tmp8));
          z_rx5[3].imag(-sin(3*tmp8));
          z_rx5[4].real(cos(4*tmp8));
          z_rx5[4].imag(-sin(4*tmp8));
//radix 7 fundament
          z_rx7[0].real(cos(0*tmp7));
          z_rx7[0].imag(-sin(0*tmp7));
          z_rx7[1].real(cos(1*tmp7));
          z_rx7[1].imag(-sin(1*tmp7));
          z_rx7[2].real(cos(2*tmp7));
          z_rx7[2].imag(-sin(2*tmp7));
          z_rx7[3].real(cos(3*tmp7));
          z_rx7[3].imag(-sin(3*tmp7));
          z_rx7[4].real(cos(4*tmp7));
          z_rx7[4].imag(-sin(4*tmp7));
          z_rx7[5].real(cos(5*tmp7));
          z_rx7[5].imag(-sin(5*tmp7));
          z_rx7[6].real(cos(6*tmp7));
          z_rx7[6].imag(-sin(6*tmp7));
//radix 11 fundament
          z_rx11[0].real(cos(0*tmp11));
          z_rx11[0].imag(-sin(0*tmp11));
          z_rx11[1].real(cos(1*tmp11));
          z_rx11[1].imag(-sin(1*tmp11));
          z_rx11[2].real(cos(2*tmp11));
          z_rx11[2].imag(-sin(2*tmp11));
          z_rx11[3].real(cos(3*tmp11));
          z_rx11[3].imag(-sin(3*tmp11));
          z_rx11[4].real(cos(4*tmp11));
          z_rx11[4].imag(-sin(4*tmp11));
          z_rx11[5].real(cos(5*tmp11));
          z_rx11[5].imag(-sin(5*tmp11));
          z_rx11[6].real(cos(6*tmp11));
          z_rx11[6].imag(-sin(6*tmp11));
          z_rx11[7].real(cos(7*tmp11));
          z_rx11[7].imag(-sin(7*tmp11));
          z_rx11[8].real(cos(8*tmp11));
          z_rx11[8].imag(-sin(8*tmp11));
          z_rx11[9].real(cos(9*tmp11));
          z_rx11[9].imag(-sin(9*tmp11));
          z_rx11[10].real(cos(10*tmp11));
          z_rx11[10].imag(-sin(10*tmp11));



        for (int i2=0;i2<97;i2++)
        {
            for (int j2=0;j2<97;j2++)
            {
                z_rx97[i2*j2].real(cos(i2*j2*tmp97));
                z_rx97[i2*j2].imag(-sin(i2*j2*tmp97));
            }
        }

/*
for(int j=0;j<102;j++)
{
    for(int i=0;i<102;i++)
    {
      if(((fabs(round(z_rx11[j].imag()*1000)/1000-round(z_rx11[i].imag()*1000)/1000)<0.001)
         &&(fabs(round(z_rx11[j].real()*1000)/1000-round(z_rx11[i].real()*1000)/1000)<0.001))
         ||((fabs(round(z_rx11[j].imag()*1000)/1000+round(z_rx11[i].imag()*1000)/1000)<0.001)
         &&(fabs(round(z_rx11[j].real()*1000)/1000+round(z_rx11[i].real()*1000)/1000)<0.001)))
         {
             cout<<j<<" "<<round(z_rx11[j].real()*1000)/1000<<" "<<round(z_rx11[j].imag()*1000)/1000<<"   ";
             cout<<i<<" "<<round(z_rx11[i].real()*1000)/1000<<" "<<round(z_rx11[i].imag()*1000)/1000<<endl;
         }
             //cout<<endl;
    }
    //system("pause");
}
*/
	nb_stages=radix_base(N,stg);

        if(nb_stages>=1){cout<<"N= "<<N<<endl;}
        for(int m=1;m<=nb_stages;m++)
        {
        cout <<"stage:"<<m<<" = radix-"<<stg[m] << endl;
        }
        cout << endl;



    if(nb_stages>=1)
    {
        stg_first=N/stg[1];
        if(stg[1]==2)
        {
            fun_fourier_transform_FFT_radix_2_stg_first(tab,stg_first,fi2,z_rx2);
        }
        else if(stg[1]==3)
        {
            fun_fourier_transform_FFT_radix_3_stg_first(tab,stg_first,fi2,z_rx3);
        }
        else if(stg[1]==4)
        {
            fun_fourier_transform_FFT_radix_4_stg_first(tab,stg_first,fi2,z_rx4);
        }
        else if(stg[1]==5)
        {
            fun_fourier_transform_FFT_radix_5_stg_first(tab,stg_first,fi2,z_rx5);
        }
        else if(stg[1]==7)
        {
            fun_fourier_transform_FFT_radix_7_stg_first(tab,stg_first,fi2,z_rx7);
        }
        else if(stg[1]==11)
        {
            fun_fourier_transform_FFT_radix_11_stg_first(tab,stg_first,fi2,z_rx11);
        }
        else if(stg[1]==97)
        {
            fun_fourier_transform_FFT_radix_97_stg_first(tab,stg_first,fi2,z_rx97);
        }
        else{}
        nb1=N;
        nb4=1;
        for(int i=0;i<nb_stages-1;i++)
        {
            nb1=nb1/stg[0+i];
            nb2=nb1/stg[1+i];
            nb3=nb2/stg[2+i];
            nb4=nb4*stg[0+i];
            nb5_stg_previous=stg[1+i];

            if(stg[i+2]==2)
            {
                fun_fourier_transform_FFT_radix_2_stg_rest(tab,nb1,nb2,nb3,nb4,nb5_stg_previous,tmp6,z_rx2);
            }
            else if(stg[i+2]==3)
            {
                fun_fourier_transform_FFT_radix_3_stg_rest(tab,nb1,nb2,nb3,nb4,nb5_stg_previous,tmp6,z_rx3);
            }
            else if(stg[i+2]==4)
            {
                fun_fourier_transform_FFT_radix_4_stg_rest(tab,nb1,nb2,nb3,nb4,nb5_stg_previous,tmp6,z_rx4);
            }
            else if(stg[i+2]==5)
            {
                fun_fourier_transform_FFT_radix_5_stg_rest(tab,nb1,nb2,nb3,nb4,nb5_stg_previous,tmp6,z_rx5);
            }
            else if(stg[i+2]==7)
            {
                fun_fourier_transform_FFT_radix_7_stg_rest(tab,nb1,nb2,nb3,nb4,nb5_stg_previous,tmp6,z_rx7);
            }
            else if(stg[i+2]==11)
            {
                fun_fourier_transform_FFT_radix_11_stg_rest(tab,nb1,nb2,nb3,nb4,nb5_stg_previous,tmp6,z_rx11);
            }
            else if(stg[i+2]==97)
            {
                fun_fourier_transform_FFT_radix_97_stg_rest(tab,nb1,nb2,nb3,nb4,nb5_stg_previous,tmp6,z_rx97);
            }
            else{}

        }
    }

//new:
    for(int j=0;j<N;j++)
    {
     tab[j].real(tab[j].real()*2/N);
     tab[j].imag(tab[j].imag()*2/N);
    }
    delete [] tab2;
}
///////////////////////////////////////////////////


  void fun_fourier_transform_FFT_radix_2_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp6,std::complex<double> z[])
    {
        std::complex<double> tmp1,tmp2,tmp40,tmp50;
        std::complex<double>  w[2]={{1,0}};
        double tmp100=0.0;
        double tmp200=0.0;
        double tmp300=0.0;
        int nb_tmp1=0.0;
        int nb_tmp2=0.0;
        int nb_tmp3=0.0;

        for(int b=0;b<nb5_stg_previous;b=b+1)
        {
            tmp300=nb4*b*tmp6;
            tmp100=nb3*tmp300;
            for(int j=0;j<nb3;j=j+1)
            {
                tmp200=j*tmp300;
                w[0].real( cos(0*tmp100+tmp200));
                w[0].imag(-sin(0*tmp100+tmp200));
                w[1].real( cos(1*tmp100+tmp200));
                //w[1].imag(-sin(nb4*b*(1*nb3+j)*tmp5));
                w[1].imag(-sin(1*tmp100+tmp200));

                nb_tmp1=b*nb2+j;
                for(int i=0;i<nb4;i=i+1)
                {
                    nb_tmp2=i*nb1;
                    nb_tmp3=nb_tmp1+nb_tmp2;
                    tmp1=w[0]*tab[nb_tmp3+0*nb3];
                    tmp2=w[1]*tab[nb_tmp3+1*nb3];

                    tmp40=z[0]*(tmp1+tmp2);

                    tab[nb_tmp3+0*nb3]=tmp40;
                    tab[nb_tmp3+1*nb3]=z[0]*tmp1+z[1]*tmp2;
                }
            }
        }
    }
///////////////////////////////////////////////////////////////


    void fun_fourier_transform_FFT_radix_3_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp6,std::complex<double> z[])
    {
        std::complex<double> tmp1,tmp2,tmp3,tmp4,tmp5;
        std::complex<double>  w[3]={{1,0}};
        double tmp100=0.0;
        double tmp200=0.0;
        double tmp300=0.0;
        int nb_tmp1=0.0;
        int nb_tmp2=0.0;
        int nb_tmp3=0.0;

        for(int b=0;b<nb5_stg_previous;b=b+1)
        {
            tmp300=nb4*b*tmp6;
            tmp100=nb3*tmp300;
            for(int j=0;j<nb3;j=j+1)
            {
                tmp200=j*tmp300;
                w[0].real(cos(0*tmp100+tmp200));
                w[0].imag(-sin(0*tmp100+tmp200));
                w[1].real(cos(1*tmp100+tmp200));
                w[1].imag(-sin(1*tmp100+tmp200));
                //w[2].real(cos(nb4*b*(2*nb3+j)*tmp6));
                w[2].real(cos(2*tmp100+tmp200));
                w[2].imag(-sin(2*tmp100+tmp200));

                nb_tmp1=b*nb2+j;
                for(int i=0;i<nb4;i=i+1)
                {
                    nb_tmp2=i*nb1;
                    nb_tmp3=nb_tmp1+nb_tmp2;
                  tmp1=w[0]*tab[nb_tmp3+0*nb3];
                  tmp2=w[1]*tab[nb_tmp3+1*nb3];
                  tmp3=w[2]*tab[nb_tmp3+2*nb3];

                  tmp4=z[0]*tmp1;
                  tmp5=z[0]*(tmp2+tmp3);

                 //radix-3
                  tab[nb_tmp3+0*nb3]  = tmp4+tmp5;
                  tab[nb_tmp3+1*nb3]  = tmp4+z[1]*tmp2+z[2]*tmp3;
                  tab[nb_tmp3+2*nb3]  = tmp4+z[2]*tmp2+z[1]*tmp3;
                }
            }
        }
    }
///////////////////////////////////////////////////////////////


    void fun_fourier_transform_FFT_radix_4_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp6,std::complex<double> z[])
    {
        std::complex<double> tmp1,tmp2,tmp3,tmp4;
        std::complex<double> tmp101,tmp102,tmp103,tmp104;
        std::complex<double> tmp111,tmp112,tmp113,tmp114;
        std::complex<double>  w[4]={{1,0}};
        double tmp100=0.0;
        double tmp200=0.0;
        double tmp300=0.0;
        int nb_tmp1=0.0;
        int nb_tmp2=0.0;
        int nb_tmp3=0.0;

        for(int b=0;b<nb5_stg_previous;b=b+1)
        {
            tmp300=nb4*b*tmp6;
            tmp100=nb3*tmp300;
            for(int j=0;j<nb3;j=j+1)
            {
              tmp200=j*tmp300;
              w[0].real(cos(0*tmp100+tmp200));
              w[0].imag(-sin(0*tmp100+tmp200));
              w[1].real(cos(1*tmp100+tmp200));
              w[1].imag(-sin(1*tmp100+tmp200));
              w[2].real(cos(2*tmp100+tmp200));
              w[2].imag(-sin(2*tmp100+tmp200));
              //w[3].real(cos(nb4*b*(3*nb3+j)*tmp5));
              w[3].real(cos(3*tmp100+tmp200));
              w[3].imag(-sin(3*tmp100+tmp200));

                nb_tmp1=b*nb2+j;
                for(int i=0;i<nb4;i=i+1)
                {
                    nb_tmp2=i*nb1;
                    nb_tmp3=nb_tmp1+nb_tmp2;
                    tmp1=w[0]*tab[nb_tmp3+0*nb3];
                    tmp2=w[1]*tab[nb_tmp3+1*nb3];
                    tmp3=w[2]*tab[nb_tmp3+2*nb3];
                    tmp4=w[3]*tab[nb_tmp3+3*nb3];


                    tmp101=tmp2-tmp4;
                    tmp102=tmp2+tmp4;
                    tmp103=tmp1-tmp3;
                    tmp104=tmp1+tmp3;

                    tmp111=z[0]*(tmp104+tmp102);
                    tmp112=z[0]*tmp103;
                    tmp113=z[0]*(tmp104-tmp102);
                    tmp114=z[0]*tmp103;

                    //radix-4
                    tab[nb_tmp3+0*nb3]   =tmp111;
                    tab[nb_tmp3+1*nb3]   =tmp112+z[1]*tmp101;
                    tab[nb_tmp3+2*nb3]   =tmp113;
                    tab[nb_tmp3+3*nb3]   =tmp114-z[1]*tmp101;
                }
            }
        }
    }
/////////////////////////////////////////////////////////////////

   void fun_fourier_transform_FFT_radix_5_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp6,std::complex<double> z[])
    {
        std::complex<double> tmp1,tmp2,tmp3,tmp4,tmp5,tmp10,tmp20;
        std::complex<double>  w[5]={{1,0}};
        double tmp100=0.0;
        double tmp200=0.0;
        double tmp300=0.0;
        int nb_tmp1=0.0;
        int nb_tmp2=0.0;
        int nb_tmp3=0.0;

        for(int b=0;b<nb5_stg_previous;b=b+1)
        {
            tmp300=nb4*b*tmp6;
            tmp100=nb3*tmp300;
            for(int j=0;j<nb3;j=j+1)
            {
                tmp200=j*tmp300;
                w[0].real(cos(0*tmp100+tmp200));
                w[0].imag(-sin(0*tmp100+tmp200));
                w[1].real(cos(1*tmp100+tmp200));
                w[1].imag(-sin(1*tmp100+tmp200));
                //w[2].real(cos(nb4*b*(2*nb3+j)*tmp6));
                w[2].real(cos(2*tmp100+tmp200));
                w[2].imag(-sin(2*tmp100+tmp200));
                w[3].real(cos(3*tmp100+tmp200));
                w[3].imag(-sin(3*tmp100+tmp200));
                w[4].real(cos(4*tmp100+tmp200));
                w[4].imag(-sin(4*tmp100+tmp200));

                nb_tmp1=b*nb2+j;
                for(int i=0;i<nb4;i=i+1)
                {
                    nb_tmp2=i*nb1;
                    nb_tmp3=nb_tmp1+nb_tmp2;
                  tmp1=w[0]*tab[nb_tmp3+0*nb3];
                  tmp2=w[1]*tab[nb_tmp3+1*nb3];
                  tmp3=w[2]*tab[nb_tmp3+2*nb3];
                  tmp4=w[3]*tab[nb_tmp3+3*nb3];
                  tmp5=w[4]*tab[nb_tmp3+4*nb3];

                  tmp10=z[0]*tmp1;
                  tmp20=z[0]*(tmp1+tmp2+tmp3+tmp4+tmp5);

                 //radix-5
                  tab[nb_tmp3+0*nb3] =tmp20;
                  tab[nb_tmp3+1*nb3] =tmp10+z[1]*tmp2+z[2]*tmp3+z[3]*tmp4+z[4]*tmp5;
                  tab[nb_tmp3+2*nb3] =tmp10+z[2]*tmp2+z[4]*tmp3+z[1]*tmp4+z[3]*tmp5;
                  tab[nb_tmp3+3*nb3] =tmp10+z[3]*tmp2+z[1]*tmp3+z[4]*tmp4+z[2]*tmp5;
                  tab[nb_tmp3+4*nb3] =tmp10+z[4]*tmp2+z[3]*tmp3+z[2]*tmp4+z[1]*tmp5;
                }
            }
        }
    }
///////////////////////////////////////////////////////////////////////

   void fun_fourier_transform_FFT_radix_7_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp99,std::complex<double> z[])
    {
        std::complex<double> tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp10,tmp20;
        std::complex<double>  w[7]={{1,0}};
        double tmp100=0.0;
        double tmp200=0.0;
        double tmp300=0.0;
        int nb_tmp1=0.0;
        int nb_tmp2=0.0;
        int nb_tmp3=0.0;

        for(int b=0;b<nb5_stg_previous;b=b+1)
        {
            tmp300=nb4*b*tmp99;
            tmp100=nb3*tmp300;
            for(int j=0;j<nb3;j=j+1)
            {
                tmp200=j*tmp300;
                w[0].real(cos(0*tmp100+tmp200));
                w[0].imag(-sin(0*tmp100+tmp200));
                w[1].real(cos(1*tmp100+tmp200));
                w[1].imag(-sin(1*tmp100+tmp200));
                //w[2].real(cos(nb4*b*(2*nb3+j)*tmp99));
                w[2].real(cos(2*tmp100+tmp200));
                w[2].imag(-sin(2*tmp100+tmp200));
                w[3].real(cos(3*tmp100+tmp200));
                w[3].imag(-sin(3*tmp100+tmp200));
                w[4].real(cos(4*tmp100+tmp200));
                w[4].imag(-sin(4*tmp100+tmp200));
                w[5].real(cos(5*tmp100+tmp200));
                w[5].imag(-sin(5*tmp100+tmp200));
                w[6].real(cos(6*tmp100+tmp200));
                w[6].imag(-sin(6*tmp100+tmp200));

                nb_tmp1=b*nb2+j;
                for(int i=0;i<nb4;i=i+1)
                {
                    nb_tmp2=i*nb1;
                    nb_tmp3=nb_tmp1+nb_tmp2;
                  tmp1=w[0]*tab[nb_tmp3+0*nb3];
                  tmp2=w[1]*tab[nb_tmp3+1*nb3];
                  tmp3=w[2]*tab[nb_tmp3+2*nb3];
                  tmp4=w[3]*tab[nb_tmp3+3*nb3];
                  tmp5=w[4]*tab[nb_tmp3+4*nb3];
                  tmp6=w[5]*tab[nb_tmp3+5*nb3];
                  tmp7=w[6]*tab[nb_tmp3+6*nb3];

                  tmp10=z[0]*tmp1;
                  tmp20=z[0]*(tmp1+tmp2+tmp3+tmp4+tmp5+tmp6+tmp7);

                 //radix-7
                  tab[nb_tmp3+0*nb3] =tmp20;
                  tab[nb_tmp3+1*nb3] =tmp10+z[1]*tmp2+z[2]*tmp3+z[3]*tmp4+z[4]*tmp5+z[5]*tmp6+z[6]*tmp7;
                  tab[nb_tmp3+2*nb3] =tmp10+z[2]*tmp2+z[4]*tmp3+z[6]*tmp4+z[1]*tmp5+z[3]*tmp6+z[5]*tmp7;
                  tab[nb_tmp3+3*nb3] =tmp10+z[3]*tmp2+z[6]*tmp3+z[2]*tmp4+z[5]*tmp5+z[1]*tmp6+z[4]*tmp7;
                  tab[nb_tmp3+4*nb3] =tmp10+z[4]*tmp2+z[1]*tmp3+z[5]*tmp4+z[2]*tmp5+z[6]*tmp6+z[3]*tmp7;
                  tab[nb_tmp3+5*nb3] =tmp10+z[5]*tmp2+z[3]*tmp3+z[1]*tmp4+z[6]*tmp5+z[4]*tmp6+z[2]*tmp7;
                  tab[nb_tmp3+6*nb3] =tmp10+z[6]*tmp2+z[5]*tmp3+z[4]*tmp4+z[3]*tmp5+z[2]*tmp6+z[1]*tmp7;
                }
            }
        }
    }
///////////////////////////////////////////////////////////////////////

   void fun_fourier_transform_FFT_radix_11_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp99,std::complex<double> z[])
    {
        std::complex<double> tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10,tmp11,tmp30,tmp20;
        std::complex<double>  w[11]={{1,0}};
        double tmp100=0.0;
        double tmp200=0.0;
        double tmp300=0.0;
        int nb_tmp1=0.0;
        int nb_tmp2=0.0;
        int nb_tmp3=0.0;

        for(int b=0;b<nb5_stg_previous;b=b+1)
        {
            tmp300=nb4*b*tmp99;
            tmp100=nb3*tmp300;
            for(int j=0;j<nb3;j=j+1)
            {
                tmp200=j*tmp300;
                w[0].real(cos(0*tmp100+tmp200));
                w[0].imag(-sin(0*tmp100+tmp200));
                w[1].real(cos(1*tmp100+tmp200));
                w[1].imag(-sin(1*tmp100+tmp200));
                //w[2].real(cos(nb4*b*(2*nb3+j)*tmp99));
                w[2].real(cos(2*tmp100+tmp200));
                w[2].imag(-sin(2*tmp100+tmp200));
                w[3].real(cos(3*tmp100+tmp200));
                w[3].imag(-sin(3*tmp100+tmp200));
                w[4].real(cos(4*tmp100+tmp200));
                w[4].imag(-sin(4*tmp100+tmp200));
                w[5].real(cos(5*tmp100+tmp200));
                w[5].imag(-sin(5*tmp100+tmp200));
                w[6].real(cos(6*tmp100+tmp200));
                w[6].imag(-sin(6*tmp100+tmp200));
                w[7].real(cos(7*tmp100+tmp200));
                w[7].imag(-sin(7*tmp100+tmp200));
                w[8].real(cos(8*tmp100+tmp200));
                w[8].imag(-sin(8*tmp100+tmp200));
                w[9].real(cos(9*tmp100+tmp200));
                w[9].imag(-sin(9*tmp100+tmp200));
                w[10].real(cos(10*tmp100+tmp200));
                w[10].imag(-sin(10*tmp100+tmp200));

                nb_tmp1=b*nb2+j;
                for(int i=0;i<nb4;i=i+1)
                {
                    nb_tmp2=i*nb1;
                    nb_tmp3=nb_tmp1+nb_tmp2;
                  tmp1=w[0]*tab[nb_tmp3+0*nb3];
                  tmp2=w[1]*tab[nb_tmp3+1*nb3];
                  tmp3=w[2]*tab[nb_tmp3+2*nb3];
                  tmp4=w[3]*tab[nb_tmp3+3*nb3];
                  tmp5=w[4]*tab[nb_tmp3+4*nb3];
                  tmp6=w[5]*tab[nb_tmp3+5*nb3];
                  tmp7=w[6]*tab[nb_tmp3+6*nb3];
                  tmp8=w[7]*tab[nb_tmp3+7*nb3];
                  tmp9=w[8]*tab[nb_tmp3+8*nb3];
                  tmp10=w[9]*tab[nb_tmp3+9*nb3];
                  tmp11=w[10]*tab[nb_tmp3+10*nb3];

                  tmp30=z[0]*tmp1;
                  tmp20=z[0]*(tmp1+tmp2+tmp3+tmp4+tmp5+tmp6+tmp7+tmp8+tmp9+tmp10+tmp11);

                 //radix-11
                  tab[nb_tmp3+0*nb3] =tmp20;
                  tab[nb_tmp3+1*nb3] =tmp30+z[1]*tmp2+z[2]*tmp3+z[3]*tmp4+z[4]*tmp5+z[5]*tmp6+z[6]*tmp7+z[7]*tmp8+z[8]*tmp9+z[9]*tmp10+z[10]*tmp11;
                  tab[nb_tmp3+2*nb3] =tmp30+z[2]*tmp2+z[4]*tmp3+z[6]*tmp4+z[8]*tmp5+z[10]*tmp6+z[1]*tmp7+z[3]*tmp8+z[5]*tmp9+z[7]*tmp10+z[9]*tmp11;
                  tab[nb_tmp3+3*nb3] =tmp30+z[3]*tmp2+z[6]*tmp3+z[9]*tmp4+z[1]*tmp5+z[4]*tmp6+z[7]*tmp7+z[10]*tmp8+z[2]*tmp9+z[5]*tmp10+z[8]*tmp11;
                  tab[nb_tmp3+4*nb3] =tmp30+z[4]*tmp2+z[8]*tmp3+z[1]*tmp4+z[5]*tmp5+z[9]*tmp6+z[2]*tmp7+z[6]*tmp8+z[10]*tmp9+z[3]*tmp10+z[7]*tmp11;
                  tab[nb_tmp3+5*nb3] =tmp30+z[5]*tmp2+z[10]*tmp3+z[4]*tmp4+z[9]*tmp5+z[3]*tmp6+z[8]*tmp7+z[2]*tmp8+z[7]*tmp9+z[1]*tmp10+z[6]*tmp11;
                  tab[nb_tmp3+6*nb3] =tmp30+z[6]*tmp2+z[1]*tmp3+z[7]*tmp4+z[2]*tmp5+z[8]*tmp6+z[3]*tmp7+z[9]*tmp8+z[4]*tmp9+z[10]*tmp10+z[5]*tmp11;
                  tab[nb_tmp3+7*nb3] =tmp30+z[7]*tmp2+z[3]*tmp3+z[10]*tmp4+z[6]*tmp5+z[2]*tmp6+z[9]*tmp7+z[5]*tmp8+z[1]*tmp9+z[8]*tmp10+z[4]*tmp11;
                  tab[nb_tmp3+8*nb3] =tmp30+z[8]*tmp2+z[5]*tmp3+z[2]*tmp4+z[10]*tmp5+z[7]*tmp6+z[4]*tmp7+z[1]*tmp8+z[9]*tmp9+z[6]*tmp10+z[3]*tmp11;
                  tab[nb_tmp3+9*nb3] =tmp30+z[9]*tmp2+z[7]*tmp3+z[5]*tmp4+z[3]*tmp5+z[1]*tmp6+z[10]*tmp7+z[8]*tmp8+z[6]*tmp9+z[4]*tmp10+z[2]*tmp11;
                  tab[nb_tmp3+10*nb3] =tmp30+z[10]*tmp2+z[9]*tmp3+z[8]*tmp4+z[7]*tmp5+z[6]*tmp6+z[5]*tmp7+z[4]*tmp8+z[3]*tmp9+z[2]*tmp10+z[1]*tmp11;
                }

            }
        }
    }
///////////////////////////////////////////////////////////////////////

void fun_fourier_transform_FFT_radix_97_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp99,std::complex<double> z[])
    {
        std::complex<double> tmp1[97],tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10,tmp11,tmp30,tmp20;
        std::complex<double>  w[97]={{1,0}};
        double tmp100=0.0;
        double tmp200=0.0;
        double tmp300=0.0;
        int nb_tmp1=0.0;
        int nb_tmp2=0.0;
        int nb_tmp3=0.0;

        for(int b=0;b<nb5_stg_previous;b=b+1)
        {
            tmp300=nb4*b*tmp99;
            tmp100=nb3*tmp300;
            for(int j=0;j<nb3;j=j+1)
            {
                tmp200=j*tmp300;
                for(int j9=0;j9<97;j9++)
                {
                     w[j9].real(cos(j9*tmp100+tmp200));
                     w[j9].imag(-sin(j9*tmp100+tmp200));
                }

                nb_tmp1=b*nb2+j;
                for(int i=0;i<nb4;i=i+1)
                {
                        nb_tmp2=i*nb1;
                        nb_tmp3=nb_tmp1+nb_tmp2;
                for(int j9=0;j9<97;j9++)
                {
                       tmp1[j9]=w[j9]*tab[nb_tmp3+j9*nb3];
                }


                         //radix-97
                   for(int i9=0;i9<97;i9++)
                  {
                    for(int j9=0;j9<97;j9++)
                    {
                        tab[nb_tmp3+i9*nb3].real(0);
                        tab[nb_tmp3+i9*nb3].imag(0);
                    }
                  }


                     tab[nb_tmp3+0*nb3]=tab[nb_tmp3+0*nb3]+z[0]*tmp1[0]
                    +z[0]*tmp1[1]+z[0]*tmp1[2]+z[0]*tmp1[3]+z[0]*tmp1[4]+z[0]*tmp1[5]+z[0]*tmp1[6]+z[0]*tmp1[7]+z[0]*tmp1[8]+z[0]*tmp1[9]+z[0]*tmp1[10]
                    +z[0]*tmp1[11]+z[0]*tmp1[12]+z[0]*tmp1[13]+z[0]*tmp1[14]+z[0]*tmp1[15]+z[0]*tmp1[16]+z[0]*tmp1[17]+z[0]*tmp1[18]+z[0]*tmp1[19]+z[0]*tmp1[20]
                    +z[0]*tmp1[21]+z[0]*tmp1[22]+z[0]*tmp1[23]+z[0]*tmp1[24]+z[0]*tmp1[25]+z[0]*tmp1[26]+z[0]*tmp1[27]+z[0]*tmp1[28]+z[0]*tmp1[29]+z[0]*tmp1[30]
                    +z[0]*tmp1[31]+z[0]*tmp1[32]+z[0]*tmp1[33]+z[0]*tmp1[34]+z[0]*tmp1[35]+z[0]*tmp1[36]+z[0]*tmp1[37]+z[0]*tmp1[38]+z[0]*tmp1[39]+z[0]*tmp1[40]
                    +z[0]*tmp1[41]+z[0]*tmp1[42]+z[0]*tmp1[43]+z[0]*tmp1[44]+z[0]*tmp1[45]+z[0]*tmp1[46]+z[0]*tmp1[47]+z[0]*tmp1[48]+z[0]*tmp1[49]+z[0]*tmp1[50]
                    +z[0]*tmp1[51]+z[0]*tmp1[52]+z[0]*tmp1[53]+z[0]*tmp1[54]+z[0]*tmp1[55]+z[0]*tmp1[56]+z[0]*tmp1[57]+z[0]*tmp1[58]+z[0]*tmp1[59]+z[0]*tmp1[60]
                    +z[0]*tmp1[61]+z[0]*tmp1[62]+z[0]*tmp1[63]+z[0]*tmp1[64]+z[0]*tmp1[65]+z[0]*tmp1[66]+z[0]*tmp1[67]+z[0]*tmp1[68]+z[0]*tmp1[69]+z[0]*tmp1[70]
                    +z[0]*tmp1[71]+z[0]*tmp1[72]+z[0]*tmp1[73]+z[0]*tmp1[74]+z[0]*tmp1[75]+z[0]*tmp1[76]+z[0]*tmp1[77]+z[0]*tmp1[78]+z[0]*tmp1[79]+z[0]*tmp1[80]
                    +z[0]*tmp1[81]+z[0]*tmp1[82]+z[0]*tmp1[83]+z[0]*tmp1[84]+z[0]*tmp1[85]+z[0]*tmp1[86]+z[0]*tmp1[87]+z[0]*tmp1[88]+z[0]*tmp1[89]+z[0]*tmp1[90]
                    +z[0]*tmp1[91]+z[0]*tmp1[92]+z[0]*tmp1[93]+z[0]*tmp1[94]+z[0]*tmp1[95]+z[0]*tmp1[96];
                    tab[nb_tmp3+1*nb3]=tab[nb_tmp3+1*nb3]+z[0]*tmp1[0]
                    +z[1]*tmp1[1]+z[2]*tmp1[2]+z[3]*tmp1[3]+z[4]*tmp1[4]+z[5]*tmp1[5]+z[6]*tmp1[6]+z[7]*tmp1[7]+z[8]*tmp1[8]+z[9]*tmp1[9]+z[10]*tmp1[10]
                    +z[11]*tmp1[11]+z[12]*tmp1[12]+z[13]*tmp1[13]+z[14]*tmp1[14]+z[15]*tmp1[15]+z[16]*tmp1[16]+z[17]*tmp1[17]+z[18]*tmp1[18]+z[19]*tmp1[19]+z[20]*tmp1[20]
                    +z[21]*tmp1[21]+z[22]*tmp1[22]+z[23]*tmp1[23]+z[24]*tmp1[24]+z[25]*tmp1[25]+z[26]*tmp1[26]+z[27]*tmp1[27]+z[28]*tmp1[28]+z[29]*tmp1[29]+z[30]*tmp1[30]
                    +z[31]*tmp1[31]+z[32]*tmp1[32]+z[33]*tmp1[33]+z[34]*tmp1[34]+z[35]*tmp1[35]+z[36]*tmp1[36]+z[37]*tmp1[37]+z[38]*tmp1[38]+z[39]*tmp1[39]+z[40]*tmp1[40]
                    +z[41]*tmp1[41]+z[42]*tmp1[42]+z[43]*tmp1[43]+z[44]*tmp1[44]+z[45]*tmp1[45]+z[46]*tmp1[46]+z[47]*tmp1[47]+z[48]*tmp1[48]+z[49]*tmp1[49]+z[50]*tmp1[50]
                    +z[51]*tmp1[51]+z[52]*tmp1[52]+z[53]*tmp1[53]+z[54]*tmp1[54]+z[55]*tmp1[55]+z[56]*tmp1[56]+z[57]*tmp1[57]+z[58]*tmp1[58]+z[59]*tmp1[59]+z[60]*tmp1[60]
                    +z[61]*tmp1[61]+z[62]*tmp1[62]+z[63]*tmp1[63]+z[64]*tmp1[64]+z[65]*tmp1[65]+z[66]*tmp1[66]+z[67]*tmp1[67]+z[68]*tmp1[68]+z[69]*tmp1[69]+z[70]*tmp1[70]
                    +z[71]*tmp1[71]+z[72]*tmp1[72]+z[73]*tmp1[73]+z[74]*tmp1[74]+z[75]*tmp1[75]+z[76]*tmp1[76]+z[77]*tmp1[77]+z[78]*tmp1[78]+z[79]*tmp1[79]+z[80]*tmp1[80]
                    +z[81]*tmp1[81]+z[82]*tmp1[82]+z[83]*tmp1[83]+z[84]*tmp1[84]+z[85]*tmp1[85]+z[86]*tmp1[86]+z[87]*tmp1[87]+z[88]*tmp1[88]+z[89]*tmp1[89]+z[90]*tmp1[90]
                    +z[91]*tmp1[91]+z[92]*tmp1[92]+z[93]*tmp1[93]+z[94]*tmp1[94]+z[95]*tmp1[95]+z[96]*tmp1[96];
                    tab[nb_tmp3+2*nb3]=tab[nb_tmp3+2*nb3]+z[0]*tmp1[0]
                    +z[2]*tmp1[1]+z[4]*tmp1[2]+z[6]*tmp1[3]+z[8]*tmp1[4]+z[10]*tmp1[5]+z[12]*tmp1[6]+z[14]*tmp1[7]+z[16]*tmp1[8]+z[18]*tmp1[9]+z[20]*tmp1[10]
                    +z[22]*tmp1[11]+z[24]*tmp1[12]+z[26]*tmp1[13]+z[28]*tmp1[14]+z[30]*tmp1[15]+z[32]*tmp1[16]+z[34]*tmp1[17]+z[36]*tmp1[18]+z[38]*tmp1[19]+z[40]*tmp1[20]
                    +z[42]*tmp1[21]+z[44]*tmp1[22]+z[46]*tmp1[23]+z[48]*tmp1[24]+z[50]*tmp1[25]+z[52]*tmp1[26]+z[54]*tmp1[27]+z[56]*tmp1[28]+z[58]*tmp1[29]+z[60]*tmp1[30]
                    +z[62]*tmp1[31]+z[64]*tmp1[32]+z[66]*tmp1[33]+z[68]*tmp1[34]+z[70]*tmp1[35]+z[72]*tmp1[36]+z[74]*tmp1[37]+z[76]*tmp1[38]+z[78]*tmp1[39]+z[80]*tmp1[40]
                    +z[82]*tmp1[41]+z[84]*tmp1[42]+z[86]*tmp1[43]+z[88]*tmp1[44]+z[90]*tmp1[45]+z[92]*tmp1[46]+z[94]*tmp1[47]+z[96]*tmp1[48]+z[1]*tmp1[49]+z[3]*tmp1[50]
                    +z[5]*tmp1[51]+z[7]*tmp1[52]+z[9]*tmp1[53]+z[11]*tmp1[54]+z[13]*tmp1[55]+z[15]*tmp1[56]+z[17]*tmp1[57]+z[19]*tmp1[58]+z[21]*tmp1[59]+z[23]*tmp1[60]
                    +z[25]*tmp1[61]+z[27]*tmp1[62]+z[29]*tmp1[63]+z[31]*tmp1[64]+z[33]*tmp1[65]+z[35]*tmp1[66]+z[37]*tmp1[67]+z[39]*tmp1[68]+z[41]*tmp1[69]+z[43]*tmp1[70]
                    +z[45]*tmp1[71]+z[47]*tmp1[72]+z[49]*tmp1[73]+z[51]*tmp1[74]+z[53]*tmp1[75]+z[55]*tmp1[76]+z[57]*tmp1[77]+z[59]*tmp1[78]+z[61]*tmp1[79]+z[63]*tmp1[80]
                    +z[65]*tmp1[81]+z[67]*tmp1[82]+z[69]*tmp1[83]+z[71]*tmp1[84]+z[73]*tmp1[85]+z[75]*tmp1[86]+z[77]*tmp1[87]+z[79]*tmp1[88]+z[81]*tmp1[89]+z[83]*tmp1[90]
                    +z[85]*tmp1[91]+z[87]*tmp1[92]+z[89]*tmp1[93]+z[91]*tmp1[94]+z[93]*tmp1[95]+z[95]*tmp1[96];
                    tab[nb_tmp3+3*nb3]=tab[nb_tmp3+3*nb3]+z[0]*tmp1[0]
                    +z[3]*tmp1[1]+z[6]*tmp1[2]+z[9]*tmp1[3]+z[12]*tmp1[4]+z[15]*tmp1[5]+z[18]*tmp1[6]+z[21]*tmp1[7]+z[24]*tmp1[8]+z[27]*tmp1[9]+z[30]*tmp1[10]
                    +z[33]*tmp1[11]+z[36]*tmp1[12]+z[39]*tmp1[13]+z[42]*tmp1[14]+z[45]*tmp1[15]+z[48]*tmp1[16]+z[51]*tmp1[17]+z[54]*tmp1[18]+z[57]*tmp1[19]+z[60]*tmp1[20]
                    +z[63]*tmp1[21]+z[66]*tmp1[22]+z[69]*tmp1[23]+z[72]*tmp1[24]+z[75]*tmp1[25]+z[78]*tmp1[26]+z[81]*tmp1[27]+z[84]*tmp1[28]+z[87]*tmp1[29]+z[90]*tmp1[30]
                    +z[93]*tmp1[31]+z[96]*tmp1[32]+z[2]*tmp1[33]+z[5]*tmp1[34]+z[8]*tmp1[35]+z[11]*tmp1[36]+z[14]*tmp1[37]+z[17]*tmp1[38]+z[20]*tmp1[39]+z[23]*tmp1[40]
                    +z[26]*tmp1[41]+z[29]*tmp1[42]+z[32]*tmp1[43]+z[35]*tmp1[44]+z[38]*tmp1[45]+z[41]*tmp1[46]+z[44]*tmp1[47]+z[47]*tmp1[48]+z[50]*tmp1[49]+z[53]*tmp1[50]
                    +z[56]*tmp1[51]+z[59]*tmp1[52]+z[62]*tmp1[53]+z[65]*tmp1[54]+z[68]*tmp1[55]+z[71]*tmp1[56]+z[74]*tmp1[57]+z[77]*tmp1[58]+z[80]*tmp1[59]+z[83]*tmp1[60]
                    +z[86]*tmp1[61]+z[89]*tmp1[62]+z[92]*tmp1[63]+z[95]*tmp1[64]+z[1]*tmp1[65]+z[4]*tmp1[66]+z[7]*tmp1[67]+z[10]*tmp1[68]+z[13]*tmp1[69]+z[16]*tmp1[70]
                    +z[19]*tmp1[71]+z[22]*tmp1[72]+z[25]*tmp1[73]+z[28]*tmp1[74]+z[31]*tmp1[75]+z[34]*tmp1[76]+z[37]*tmp1[77]+z[40]*tmp1[78]+z[43]*tmp1[79]+z[46]*tmp1[80]
                    +z[49]*tmp1[81]+z[52]*tmp1[82]+z[55]*tmp1[83]+z[58]*tmp1[84]+z[61]*tmp1[85]+z[64]*tmp1[86]+z[67]*tmp1[87]+z[70]*tmp1[88]+z[73]*tmp1[89]+z[76]*tmp1[90]
                    +z[79]*tmp1[91]+z[82]*tmp1[92]+z[85]*tmp1[93]+z[88]*tmp1[94]+z[91]*tmp1[95]+z[94]*tmp1[96];
                    tab[nb_tmp3+4*nb3]=tab[nb_tmp3+4*nb3]+z[0]*tmp1[0]
                    +z[4]*tmp1[1]+z[8]*tmp1[2]+z[12]*tmp1[3]+z[16]*tmp1[4]+z[20]*tmp1[5]+z[24]*tmp1[6]+z[28]*tmp1[7]+z[32]*tmp1[8]+z[36]*tmp1[9]+z[40]*tmp1[10]
                    +z[44]*tmp1[11]+z[48]*tmp1[12]+z[52]*tmp1[13]+z[56]*tmp1[14]+z[60]*tmp1[15]+z[64]*tmp1[16]+z[68]*tmp1[17]+z[72]*tmp1[18]+z[76]*tmp1[19]+z[80]*tmp1[20]
                    +z[84]*tmp1[21]+z[88]*tmp1[22]+z[92]*tmp1[23]+z[96]*tmp1[24]+z[3]*tmp1[25]+z[7]*tmp1[26]+z[11]*tmp1[27]+z[15]*tmp1[28]+z[19]*tmp1[29]+z[23]*tmp1[30]
                    +z[27]*tmp1[31]+z[31]*tmp1[32]+z[35]*tmp1[33]+z[39]*tmp1[34]+z[43]*tmp1[35]+z[47]*tmp1[36]+z[51]*tmp1[37]+z[55]*tmp1[38]+z[59]*tmp1[39]+z[63]*tmp1[40]
                    +z[67]*tmp1[41]+z[71]*tmp1[42]+z[75]*tmp1[43]+z[79]*tmp1[44]+z[83]*tmp1[45]+z[87]*tmp1[46]+z[91]*tmp1[47]+z[95]*tmp1[48]+z[2]*tmp1[49]+z[6]*tmp1[50]
                    +z[10]*tmp1[51]+z[14]*tmp1[52]+z[18]*tmp1[53]+z[22]*tmp1[54]+z[26]*tmp1[55]+z[30]*tmp1[56]+z[34]*tmp1[57]+z[38]*tmp1[58]+z[42]*tmp1[59]+z[46]*tmp1[60]
                    +z[50]*tmp1[61]+z[54]*tmp1[62]+z[58]*tmp1[63]+z[62]*tmp1[64]+z[66]*tmp1[65]+z[70]*tmp1[66]+z[74]*tmp1[67]+z[78]*tmp1[68]+z[82]*tmp1[69]+z[86]*tmp1[70]
                    +z[90]*tmp1[71]+z[94]*tmp1[72]+z[1]*tmp1[73]+z[5]*tmp1[74]+z[9]*tmp1[75]+z[13]*tmp1[76]+z[17]*tmp1[77]+z[21]*tmp1[78]+z[25]*tmp1[79]+z[29]*tmp1[80]
                    +z[33]*tmp1[81]+z[37]*tmp1[82]+z[41]*tmp1[83]+z[45]*tmp1[84]+z[49]*tmp1[85]+z[53]*tmp1[86]+z[57]*tmp1[87]+z[61]*tmp1[88]+z[65]*tmp1[89]+z[69]*tmp1[90]
                    +z[73]*tmp1[91]+z[77]*tmp1[92]+z[81]*tmp1[93]+z[85]*tmp1[94]+z[89]*tmp1[95]+z[93]*tmp1[96];
                    tab[nb_tmp3+5*nb3]=tab[nb_tmp3+5*nb3]+z[0]*tmp1[0]
                    +z[5]*tmp1[1]+z[10]*tmp1[2]+z[15]*tmp1[3]+z[20]*tmp1[4]+z[25]*tmp1[5]+z[30]*tmp1[6]+z[35]*tmp1[7]+z[40]*tmp1[8]+z[45]*tmp1[9]+z[50]*tmp1[10]
                    +z[55]*tmp1[11]+z[60]*tmp1[12]+z[65]*tmp1[13]+z[70]*tmp1[14]+z[75]*tmp1[15]+z[80]*tmp1[16]+z[85]*tmp1[17]+z[90]*tmp1[18]+z[95]*tmp1[19]+z[3]*tmp1[20]
                    +z[8]*tmp1[21]+z[13]*tmp1[22]+z[18]*tmp1[23]+z[23]*tmp1[24]+z[28]*tmp1[25]+z[33]*tmp1[26]+z[38]*tmp1[27]+z[43]*tmp1[28]+z[48]*tmp1[29]+z[53]*tmp1[30]
                    +z[58]*tmp1[31]+z[63]*tmp1[32]+z[68]*tmp1[33]+z[73]*tmp1[34]+z[78]*tmp1[35]+z[83]*tmp1[36]+z[88]*tmp1[37]+z[93]*tmp1[38]+z[1]*tmp1[39]+z[6]*tmp1[40]
                    +z[11]*tmp1[41]+z[16]*tmp1[42]+z[21]*tmp1[43]+z[26]*tmp1[44]+z[31]*tmp1[45]+z[36]*tmp1[46]+z[41]*tmp1[47]+z[46]*tmp1[48]+z[51]*tmp1[49]+z[56]*tmp1[50]
                    +z[61]*tmp1[51]+z[66]*tmp1[52]+z[71]*tmp1[53]+z[76]*tmp1[54]+z[81]*tmp1[55]+z[86]*tmp1[56]+z[91]*tmp1[57]+z[96]*tmp1[58]+z[4]*tmp1[59]+z[9]*tmp1[60]
                    +z[14]*tmp1[61]+z[19]*tmp1[62]+z[24]*tmp1[63]+z[29]*tmp1[64]+z[34]*tmp1[65]+z[39]*tmp1[66]+z[44]*tmp1[67]+z[49]*tmp1[68]+z[54]*tmp1[69]+z[59]*tmp1[70]
                    +z[64]*tmp1[71]+z[69]*tmp1[72]+z[74]*tmp1[73]+z[79]*tmp1[74]+z[84]*tmp1[75]+z[89]*tmp1[76]+z[94]*tmp1[77]+z[2]*tmp1[78]+z[7]*tmp1[79]+z[12]*tmp1[80]
                    +z[17]*tmp1[81]+z[22]*tmp1[82]+z[27]*tmp1[83]+z[32]*tmp1[84]+z[37]*tmp1[85]+z[42]*tmp1[86]+z[47]*tmp1[87]+z[52]*tmp1[88]+z[57]*tmp1[89]+z[62]*tmp1[90]
                    +z[67]*tmp1[91]+z[72]*tmp1[92]+z[77]*tmp1[93]+z[82]*tmp1[94]+z[87]*tmp1[95]+z[92]*tmp1[96];
                    tab[nb_tmp3+6*nb3]=tab[nb_tmp3+6*nb3]+z[0]*tmp1[0]
                    +z[6]*tmp1[1]+z[12]*tmp1[2]+z[18]*tmp1[3]+z[24]*tmp1[4]+z[30]*tmp1[5]+z[36]*tmp1[6]+z[42]*tmp1[7]+z[48]*tmp1[8]+z[54]*tmp1[9]+z[60]*tmp1[10]
                    +z[66]*tmp1[11]+z[72]*tmp1[12]+z[78]*tmp1[13]+z[84]*tmp1[14]+z[90]*tmp1[15]+z[96]*tmp1[16]+z[5]*tmp1[17]+z[11]*tmp1[18]+z[17]*tmp1[19]+z[23]*tmp1[20]
                    +z[29]*tmp1[21]+z[35]*tmp1[22]+z[41]*tmp1[23]+z[47]*tmp1[24]+z[53]*tmp1[25]+z[59]*tmp1[26]+z[65]*tmp1[27]+z[71]*tmp1[28]+z[77]*tmp1[29]+z[83]*tmp1[30]
                    +z[89]*tmp1[31]+z[95]*tmp1[32]+z[4]*tmp1[33]+z[10]*tmp1[34]+z[16]*tmp1[35]+z[22]*tmp1[36]+z[28]*tmp1[37]+z[34]*tmp1[38]+z[40]*tmp1[39]+z[46]*tmp1[40]
                    +z[52]*tmp1[41]+z[58]*tmp1[42]+z[64]*tmp1[43]+z[70]*tmp1[44]+z[76]*tmp1[45]+z[82]*tmp1[46]+z[88]*tmp1[47]+z[94]*tmp1[48]+z[3]*tmp1[49]+z[9]*tmp1[50]
                    +z[15]*tmp1[51]+z[21]*tmp1[52]+z[27]*tmp1[53]+z[33]*tmp1[54]+z[39]*tmp1[55]+z[45]*tmp1[56]+z[51]*tmp1[57]+z[57]*tmp1[58]+z[63]*tmp1[59]+z[69]*tmp1[60]
                    +z[75]*tmp1[61]+z[81]*tmp1[62]+z[87]*tmp1[63]+z[93]*tmp1[64]+z[2]*tmp1[65]+z[8]*tmp1[66]+z[14]*tmp1[67]+z[20]*tmp1[68]+z[26]*tmp1[69]+z[32]*tmp1[70]
                    +z[38]*tmp1[71]+z[44]*tmp1[72]+z[50]*tmp1[73]+z[56]*tmp1[74]+z[62]*tmp1[75]+z[68]*tmp1[76]+z[74]*tmp1[77]+z[80]*tmp1[78]+z[86]*tmp1[79]+z[92]*tmp1[80]
                    +z[1]*tmp1[81]+z[7]*tmp1[82]+z[13]*tmp1[83]+z[19]*tmp1[84]+z[25]*tmp1[85]+z[31]*tmp1[86]+z[37]*tmp1[87]+z[43]*tmp1[88]+z[49]*tmp1[89]+z[55]*tmp1[90]
                    +z[61]*tmp1[91]+z[67]*tmp1[92]+z[73]*tmp1[93]+z[79]*tmp1[94]+z[85]*tmp1[95]+z[91]*tmp1[96];
                    tab[nb_tmp3+7*nb3]=tab[nb_tmp3+7*nb3]+z[0]*tmp1[0]
                    +z[7]*tmp1[1]+z[14]*tmp1[2]+z[21]*tmp1[3]+z[28]*tmp1[4]+z[35]*tmp1[5]+z[42]*tmp1[6]+z[49]*tmp1[7]+z[56]*tmp1[8]+z[63]*tmp1[9]+z[70]*tmp1[10]
                    +z[77]*tmp1[11]+z[84]*tmp1[12]+z[91]*tmp1[13]+z[1]*tmp1[14]+z[8]*tmp1[15]+z[15]*tmp1[16]+z[22]*tmp1[17]+z[29]*tmp1[18]+z[36]*tmp1[19]+z[43]*tmp1[20]
                    +z[50]*tmp1[21]+z[57]*tmp1[22]+z[64]*tmp1[23]+z[71]*tmp1[24]+z[78]*tmp1[25]+z[85]*tmp1[26]+z[92]*tmp1[27]+z[2]*tmp1[28]+z[9]*tmp1[29]+z[16]*tmp1[30]
                    +z[23]*tmp1[31]+z[30]*tmp1[32]+z[37]*tmp1[33]+z[44]*tmp1[34]+z[51]*tmp1[35]+z[58]*tmp1[36]+z[65]*tmp1[37]+z[72]*tmp1[38]+z[79]*tmp1[39]+z[86]*tmp1[40]
                    +z[93]*tmp1[41]+z[3]*tmp1[42]+z[10]*tmp1[43]+z[17]*tmp1[44]+z[24]*tmp1[45]+z[31]*tmp1[46]+z[38]*tmp1[47]+z[45]*tmp1[48]+z[52]*tmp1[49]+z[59]*tmp1[50]
                    +z[66]*tmp1[51]+z[73]*tmp1[52]+z[80]*tmp1[53]+z[87]*tmp1[54]+z[94]*tmp1[55]+z[4]*tmp1[56]+z[11]*tmp1[57]+z[18]*tmp1[58]+z[25]*tmp1[59]+z[32]*tmp1[60]
                    +z[39]*tmp1[61]+z[46]*tmp1[62]+z[53]*tmp1[63]+z[60]*tmp1[64]+z[67]*tmp1[65]+z[74]*tmp1[66]+z[81]*tmp1[67]+z[88]*tmp1[68]+z[95]*tmp1[69]+z[5]*tmp1[70]
                    +z[12]*tmp1[71]+z[19]*tmp1[72]+z[26]*tmp1[73]+z[33]*tmp1[74]+z[40]*tmp1[75]+z[47]*tmp1[76]+z[54]*tmp1[77]+z[61]*tmp1[78]+z[68]*tmp1[79]+z[75]*tmp1[80]
                    +z[82]*tmp1[81]+z[89]*tmp1[82]+z[96]*tmp1[83]+z[6]*tmp1[84]+z[13]*tmp1[85]+z[20]*tmp1[86]+z[27]*tmp1[87]+z[34]*tmp1[88]+z[41]*tmp1[89]+z[48]*tmp1[90]
                    +z[55]*tmp1[91]+z[62]*tmp1[92]+z[69]*tmp1[93]+z[76]*tmp1[94]+z[83]*tmp1[95]+z[90]*tmp1[96];
                    tab[nb_tmp3+8*nb3]=tab[nb_tmp3+8*nb3]+z[0]*tmp1[0]
                    +z[8]*tmp1[1]+z[16]*tmp1[2]+z[24]*tmp1[3]+z[32]*tmp1[4]+z[40]*tmp1[5]+z[48]*tmp1[6]+z[56]*tmp1[7]+z[64]*tmp1[8]+z[72]*tmp1[9]+z[80]*tmp1[10]
                    +z[88]*tmp1[11]+z[96]*tmp1[12]+z[7]*tmp1[13]+z[15]*tmp1[14]+z[23]*tmp1[15]+z[31]*tmp1[16]+z[39]*tmp1[17]+z[47]*tmp1[18]+z[55]*tmp1[19]+z[63]*tmp1[20]
                    +z[71]*tmp1[21]+z[79]*tmp1[22]+z[87]*tmp1[23]+z[95]*tmp1[24]+z[6]*tmp1[25]+z[14]*tmp1[26]+z[22]*tmp1[27]+z[30]*tmp1[28]+z[38]*tmp1[29]+z[46]*tmp1[30]
                    +z[54]*tmp1[31]+z[62]*tmp1[32]+z[70]*tmp1[33]+z[78]*tmp1[34]+z[86]*tmp1[35]+z[94]*tmp1[36]+z[5]*tmp1[37]+z[13]*tmp1[38]+z[21]*tmp1[39]+z[29]*tmp1[40]
                    +z[37]*tmp1[41]+z[45]*tmp1[42]+z[53]*tmp1[43]+z[61]*tmp1[44]+z[69]*tmp1[45]+z[77]*tmp1[46]+z[85]*tmp1[47]+z[93]*tmp1[48]+z[4]*tmp1[49]+z[12]*tmp1[50]
                    +z[20]*tmp1[51]+z[28]*tmp1[52]+z[36]*tmp1[53]+z[44]*tmp1[54]+z[52]*tmp1[55]+z[60]*tmp1[56]+z[68]*tmp1[57]+z[76]*tmp1[58]+z[84]*tmp1[59]+z[92]*tmp1[60]
                    +z[3]*tmp1[61]+z[11]*tmp1[62]+z[19]*tmp1[63]+z[27]*tmp1[64]+z[35]*tmp1[65]+z[43]*tmp1[66]+z[51]*tmp1[67]+z[59]*tmp1[68]+z[67]*tmp1[69]+z[75]*tmp1[70]
                    +z[83]*tmp1[71]+z[91]*tmp1[72]+z[2]*tmp1[73]+z[10]*tmp1[74]+z[18]*tmp1[75]+z[26]*tmp1[76]+z[34]*tmp1[77]+z[42]*tmp1[78]+z[50]*tmp1[79]+z[58]*tmp1[80]
                    +z[66]*tmp1[81]+z[74]*tmp1[82]+z[82]*tmp1[83]+z[90]*tmp1[84]+z[1]*tmp1[85]+z[9]*tmp1[86]+z[17]*tmp1[87]+z[25]*tmp1[88]+z[33]*tmp1[89]+z[41]*tmp1[90]
                    +z[49]*tmp1[91]+z[57]*tmp1[92]+z[65]*tmp1[93]+z[73]*tmp1[94]+z[81]*tmp1[95]+z[89]*tmp1[96];
                    tab[nb_tmp3+9*nb3]=tab[nb_tmp3+9*nb3]+z[0]*tmp1[0]
                    +z[9]*tmp1[1]+z[18]*tmp1[2]+z[27]*tmp1[3]+z[36]*tmp1[4]+z[45]*tmp1[5]+z[54]*tmp1[6]+z[63]*tmp1[7]+z[72]*tmp1[8]+z[81]*tmp1[9]+z[90]*tmp1[10]
                    +z[2]*tmp1[11]+z[11]*tmp1[12]+z[20]*tmp1[13]+z[29]*tmp1[14]+z[38]*tmp1[15]+z[47]*tmp1[16]+z[56]*tmp1[17]+z[65]*tmp1[18]+z[74]*tmp1[19]+z[83]*tmp1[20]
                    +z[92]*tmp1[21]+z[4]*tmp1[22]+z[13]*tmp1[23]+z[22]*tmp1[24]+z[31]*tmp1[25]+z[40]*tmp1[26]+z[49]*tmp1[27]+z[58]*tmp1[28]+z[67]*tmp1[29]+z[76]*tmp1[30]
                    +z[85]*tmp1[31]+z[94]*tmp1[32]+z[6]*tmp1[33]+z[15]*tmp1[34]+z[24]*tmp1[35]+z[33]*tmp1[36]+z[42]*tmp1[37]+z[51]*tmp1[38]+z[60]*tmp1[39]+z[69]*tmp1[40]
                    +z[78]*tmp1[41]+z[87]*tmp1[42]+z[96]*tmp1[43]+z[8]*tmp1[44]+z[17]*tmp1[45]+z[26]*tmp1[46]+z[35]*tmp1[47]+z[44]*tmp1[48]+z[53]*tmp1[49]+z[62]*tmp1[50]
                    +z[71]*tmp1[51]+z[80]*tmp1[52]+z[89]*tmp1[53]+z[1]*tmp1[54]+z[10]*tmp1[55]+z[19]*tmp1[56]+z[28]*tmp1[57]+z[37]*tmp1[58]+z[46]*tmp1[59]+z[55]*tmp1[60]
                    +z[64]*tmp1[61]+z[73]*tmp1[62]+z[82]*tmp1[63]+z[91]*tmp1[64]+z[3]*tmp1[65]+z[12]*tmp1[66]+z[21]*tmp1[67]+z[30]*tmp1[68]+z[39]*tmp1[69]+z[48]*tmp1[70]
                    +z[57]*tmp1[71]+z[66]*tmp1[72]+z[75]*tmp1[73]+z[84]*tmp1[74]+z[93]*tmp1[75]+z[5]*tmp1[76]+z[14]*tmp1[77]+z[23]*tmp1[78]+z[32]*tmp1[79]+z[41]*tmp1[80]
                    +z[50]*tmp1[81]+z[59]*tmp1[82]+z[68]*tmp1[83]+z[77]*tmp1[84]+z[86]*tmp1[85]+z[95]*tmp1[86]+z[7]*tmp1[87]+z[16]*tmp1[88]+z[25]*tmp1[89]+z[34]*tmp1[90]
                    +z[43]*tmp1[91]+z[52]*tmp1[92]+z[61]*tmp1[93]+z[70]*tmp1[94]+z[79]*tmp1[95]+z[88]*tmp1[96];
                    tab[nb_tmp3+10*nb3]=tab[nb_tmp3+10*nb3]+z[0]*tmp1[0]
                    +z[10]*tmp1[1]+z[20]*tmp1[2]+z[30]*tmp1[3]+z[40]*tmp1[4]+z[50]*tmp1[5]+z[60]*tmp1[6]+z[70]*tmp1[7]+z[80]*tmp1[8]+z[90]*tmp1[9]+z[3]*tmp1[10]
                    +z[13]*tmp1[11]+z[23]*tmp1[12]+z[33]*tmp1[13]+z[43]*tmp1[14]+z[53]*tmp1[15]+z[63]*tmp1[16]+z[73]*tmp1[17]+z[83]*tmp1[18]+z[93]*tmp1[19]+z[6]*tmp1[20]
                    +z[16]*tmp1[21]+z[26]*tmp1[22]+z[36]*tmp1[23]+z[46]*tmp1[24]+z[56]*tmp1[25]+z[66]*tmp1[26]+z[76]*tmp1[27]+z[86]*tmp1[28]+z[96]*tmp1[29]+z[9]*tmp1[30]
                    +z[19]*tmp1[31]+z[29]*tmp1[32]+z[39]*tmp1[33]+z[49]*tmp1[34]+z[59]*tmp1[35]+z[69]*tmp1[36]+z[79]*tmp1[37]+z[89]*tmp1[38]+z[2]*tmp1[39]+z[12]*tmp1[40]
                    +z[22]*tmp1[41]+z[32]*tmp1[42]+z[42]*tmp1[43]+z[52]*tmp1[44]+z[62]*tmp1[45]+z[72]*tmp1[46]+z[82]*tmp1[47]+z[92]*tmp1[48]+z[5]*tmp1[49]+z[15]*tmp1[50]
                    +z[25]*tmp1[51]+z[35]*tmp1[52]+z[45]*tmp1[53]+z[55]*tmp1[54]+z[65]*tmp1[55]+z[75]*tmp1[56]+z[85]*tmp1[57]+z[95]*tmp1[58]+z[8]*tmp1[59]+z[18]*tmp1[60]
                    +z[28]*tmp1[61]+z[38]*tmp1[62]+z[48]*tmp1[63]+z[58]*tmp1[64]+z[68]*tmp1[65]+z[78]*tmp1[66]+z[88]*tmp1[67]+z[1]*tmp1[68]+z[11]*tmp1[69]+z[21]*tmp1[70]
                    +z[31]*tmp1[71]+z[41]*tmp1[72]+z[51]*tmp1[73]+z[61]*tmp1[74]+z[71]*tmp1[75]+z[81]*tmp1[76]+z[91]*tmp1[77]+z[4]*tmp1[78]+z[14]*tmp1[79]+z[24]*tmp1[80]
                    +z[34]*tmp1[81]+z[44]*tmp1[82]+z[54]*tmp1[83]+z[64]*tmp1[84]+z[74]*tmp1[85]+z[84]*tmp1[86]+z[94]*tmp1[87]+z[7]*tmp1[88]+z[17]*tmp1[89]+z[27]*tmp1[90]
                    +z[37]*tmp1[91]+z[47]*tmp1[92]+z[57]*tmp1[93]+z[67]*tmp1[94]+z[77]*tmp1[95]+z[87]*tmp1[96];
                    tab[nb_tmp3+11*nb3]=tab[nb_tmp3+11*nb3]+z[0]*tmp1[0]
                    +z[11]*tmp1[1]+z[22]*tmp1[2]+z[33]*tmp1[3]+z[44]*tmp1[4]+z[55]*tmp1[5]+z[66]*tmp1[6]+z[77]*tmp1[7]+z[88]*tmp1[8]+z[2]*tmp1[9]+z[13]*tmp1[10]
                    +z[24]*tmp1[11]+z[35]*tmp1[12]+z[46]*tmp1[13]+z[57]*tmp1[14]+z[68]*tmp1[15]+z[79]*tmp1[16]+z[90]*tmp1[17]+z[4]*tmp1[18]+z[15]*tmp1[19]+z[26]*tmp1[20]
                    +z[37]*tmp1[21]+z[48]*tmp1[22]+z[59]*tmp1[23]+z[70]*tmp1[24]+z[81]*tmp1[25]+z[92]*tmp1[26]+z[6]*tmp1[27]+z[17]*tmp1[28]+z[28]*tmp1[29]+z[39]*tmp1[30]
                    +z[50]*tmp1[31]+z[61]*tmp1[32]+z[72]*tmp1[33]+z[83]*tmp1[34]+z[94]*tmp1[35]+z[8]*tmp1[36]+z[19]*tmp1[37]+z[30]*tmp1[38]+z[41]*tmp1[39]+z[52]*tmp1[40]
                    +z[63]*tmp1[41]+z[74]*tmp1[42]+z[85]*tmp1[43]+z[96]*tmp1[44]+z[10]*tmp1[45]+z[21]*tmp1[46]+z[32]*tmp1[47]+z[43]*tmp1[48]+z[54]*tmp1[49]+z[65]*tmp1[50]
                    +z[76]*tmp1[51]+z[87]*tmp1[52]+z[1]*tmp1[53]+z[12]*tmp1[54]+z[23]*tmp1[55]+z[34]*tmp1[56]+z[45]*tmp1[57]+z[56]*tmp1[58]+z[67]*tmp1[59]+z[78]*tmp1[60]
                    +z[89]*tmp1[61]+z[3]*tmp1[62]+z[14]*tmp1[63]+z[25]*tmp1[64]+z[36]*tmp1[65]+z[47]*tmp1[66]+z[58]*tmp1[67]+z[69]*tmp1[68]+z[80]*tmp1[69]+z[91]*tmp1[70]
                    +z[5]*tmp1[71]+z[16]*tmp1[72]+z[27]*tmp1[73]+z[38]*tmp1[74]+z[49]*tmp1[75]+z[60]*tmp1[76]+z[71]*tmp1[77]+z[82]*tmp1[78]+z[93]*tmp1[79]+z[7]*tmp1[80]
                    +z[18]*tmp1[81]+z[29]*tmp1[82]+z[40]*tmp1[83]+z[51]*tmp1[84]+z[62]*tmp1[85]+z[73]*tmp1[86]+z[84]*tmp1[87]+z[95]*tmp1[88]+z[9]*tmp1[89]+z[20]*tmp1[90]
                    +z[31]*tmp1[91]+z[42]*tmp1[92]+z[53]*tmp1[93]+z[64]*tmp1[94]+z[75]*tmp1[95]+z[86]*tmp1[96];
                    tab[nb_tmp3+12*nb3]=tab[nb_tmp3+12*nb3]+z[0]*tmp1[0]
                    +z[12]*tmp1[1]+z[24]*tmp1[2]+z[36]*tmp1[3]+z[48]*tmp1[4]+z[60]*tmp1[5]+z[72]*tmp1[6]+z[84]*tmp1[7]+z[96]*tmp1[8]+z[11]*tmp1[9]+z[23]*tmp1[10]
                    +z[35]*tmp1[11]+z[47]*tmp1[12]+z[59]*tmp1[13]+z[71]*tmp1[14]+z[83]*tmp1[15]+z[95]*tmp1[16]+z[10]*tmp1[17]+z[22]*tmp1[18]+z[34]*tmp1[19]+z[46]*tmp1[20]
                    +z[58]*tmp1[21]+z[70]*tmp1[22]+z[82]*tmp1[23]+z[94]*tmp1[24]+z[9]*tmp1[25]+z[21]*tmp1[26]+z[33]*tmp1[27]+z[45]*tmp1[28]+z[57]*tmp1[29]+z[69]*tmp1[30]
                    +z[81]*tmp1[31]+z[93]*tmp1[32]+z[8]*tmp1[33]+z[20]*tmp1[34]+z[32]*tmp1[35]+z[44]*tmp1[36]+z[56]*tmp1[37]+z[68]*tmp1[38]+z[80]*tmp1[39]+z[92]*tmp1[40]
                    +z[7]*tmp1[41]+z[19]*tmp1[42]+z[31]*tmp1[43]+z[43]*tmp1[44]+z[55]*tmp1[45]+z[67]*tmp1[46]+z[79]*tmp1[47]+z[91]*tmp1[48]+z[6]*tmp1[49]+z[18]*tmp1[50]
                    +z[30]*tmp1[51]+z[42]*tmp1[52]+z[54]*tmp1[53]+z[66]*tmp1[54]+z[78]*tmp1[55]+z[90]*tmp1[56]+z[5]*tmp1[57]+z[17]*tmp1[58]+z[29]*tmp1[59]+z[41]*tmp1[60]
                    +z[53]*tmp1[61]+z[65]*tmp1[62]+z[77]*tmp1[63]+z[89]*tmp1[64]+z[4]*tmp1[65]+z[16]*tmp1[66]+z[28]*tmp1[67]+z[40]*tmp1[68]+z[52]*tmp1[69]+z[64]*tmp1[70]
                    +z[76]*tmp1[71]+z[88]*tmp1[72]+z[3]*tmp1[73]+z[15]*tmp1[74]+z[27]*tmp1[75]+z[39]*tmp1[76]+z[51]*tmp1[77]+z[63]*tmp1[78]+z[75]*tmp1[79]+z[87]*tmp1[80]
                    +z[2]*tmp1[81]+z[14]*tmp1[82]+z[26]*tmp1[83]+z[38]*tmp1[84]+z[50]*tmp1[85]+z[62]*tmp1[86]+z[74]*tmp1[87]+z[86]*tmp1[88]+z[1]*tmp1[89]+z[13]*tmp1[90]
                    +z[25]*tmp1[91]+z[37]*tmp1[92]+z[49]*tmp1[93]+z[61]*tmp1[94]+z[73]*tmp1[95]+z[85]*tmp1[96];
                    tab[nb_tmp3+13*nb3]=tab[nb_tmp3+13*nb3]+z[0]*tmp1[0]
                    +z[13]*tmp1[1]+z[26]*tmp1[2]+z[39]*tmp1[3]+z[52]*tmp1[4]+z[65]*tmp1[5]+z[78]*tmp1[6]+z[91]*tmp1[7]+z[7]*tmp1[8]+z[20]*tmp1[9]+z[33]*tmp1[10]
                    +z[46]*tmp1[11]+z[59]*tmp1[12]+z[72]*tmp1[13]+z[85]*tmp1[14]+z[1]*tmp1[15]+z[14]*tmp1[16]+z[27]*tmp1[17]+z[40]*tmp1[18]+z[53]*tmp1[19]+z[66]*tmp1[20]
                    +z[79]*tmp1[21]+z[92]*tmp1[22]+z[8]*tmp1[23]+z[21]*tmp1[24]+z[34]*tmp1[25]+z[47]*tmp1[26]+z[60]*tmp1[27]+z[73]*tmp1[28]+z[86]*tmp1[29]+z[2]*tmp1[30]
                    +z[15]*tmp1[31]+z[28]*tmp1[32]+z[41]*tmp1[33]+z[54]*tmp1[34]+z[67]*tmp1[35]+z[80]*tmp1[36]+z[93]*tmp1[37]+z[9]*tmp1[38]+z[22]*tmp1[39]+z[35]*tmp1[40]
                    +z[48]*tmp1[41]+z[61]*tmp1[42]+z[74]*tmp1[43]+z[87]*tmp1[44]+z[3]*tmp1[45]+z[16]*tmp1[46]+z[29]*tmp1[47]+z[42]*tmp1[48]+z[55]*tmp1[49]+z[68]*tmp1[50]
                    +z[81]*tmp1[51]+z[94]*tmp1[52]+z[10]*tmp1[53]+z[23]*tmp1[54]+z[36]*tmp1[55]+z[49]*tmp1[56]+z[62]*tmp1[57]+z[75]*tmp1[58]+z[88]*tmp1[59]+z[4]*tmp1[60]
                    +z[17]*tmp1[61]+z[30]*tmp1[62]+z[43]*tmp1[63]+z[56]*tmp1[64]+z[69]*tmp1[65]+z[82]*tmp1[66]+z[95]*tmp1[67]+z[11]*tmp1[68]+z[24]*tmp1[69]+z[37]*tmp1[70]
                    +z[50]*tmp1[71]+z[63]*tmp1[72]+z[76]*tmp1[73]+z[89]*tmp1[74]+z[5]*tmp1[75]+z[18]*tmp1[76]+z[31]*tmp1[77]+z[44]*tmp1[78]+z[57]*tmp1[79]+z[70]*tmp1[80]
                    +z[83]*tmp1[81]+z[96]*tmp1[82]+z[12]*tmp1[83]+z[25]*tmp1[84]+z[38]*tmp1[85]+z[51]*tmp1[86]+z[64]*tmp1[87]+z[77]*tmp1[88]+z[90]*tmp1[89]+z[6]*tmp1[90]
                    +z[19]*tmp1[91]+z[32]*tmp1[92]+z[45]*tmp1[93]+z[58]*tmp1[94]+z[71]*tmp1[95]+z[84]*tmp1[96];
                    tab[nb_tmp3+14*nb3]=tab[nb_tmp3+14*nb3]+z[0]*tmp1[0]
                    +z[14]*tmp1[1]+z[28]*tmp1[2]+z[42]*tmp1[3]+z[56]*tmp1[4]+z[70]*tmp1[5]+z[84]*tmp1[6]+z[1]*tmp1[7]+z[15]*tmp1[8]+z[29]*tmp1[9]+z[43]*tmp1[10]
                    +z[57]*tmp1[11]+z[71]*tmp1[12]+z[85]*tmp1[13]+z[2]*tmp1[14]+z[16]*tmp1[15]+z[30]*tmp1[16]+z[44]*tmp1[17]+z[58]*tmp1[18]+z[72]*tmp1[19]+z[86]*tmp1[20]
                    +z[3]*tmp1[21]+z[17]*tmp1[22]+z[31]*tmp1[23]+z[45]*tmp1[24]+z[59]*tmp1[25]+z[73]*tmp1[26]+z[87]*tmp1[27]+z[4]*tmp1[28]+z[18]*tmp1[29]+z[32]*tmp1[30]
                    +z[46]*tmp1[31]+z[60]*tmp1[32]+z[74]*tmp1[33]+z[88]*tmp1[34]+z[5]*tmp1[35]+z[19]*tmp1[36]+z[33]*tmp1[37]+z[47]*tmp1[38]+z[61]*tmp1[39]+z[75]*tmp1[40]
                    +z[89]*tmp1[41]+z[6]*tmp1[42]+z[20]*tmp1[43]+z[34]*tmp1[44]+z[48]*tmp1[45]+z[62]*tmp1[46]+z[76]*tmp1[47]+z[90]*tmp1[48]+z[7]*tmp1[49]+z[21]*tmp1[50]
                    +z[35]*tmp1[51]+z[49]*tmp1[52]+z[63]*tmp1[53]+z[77]*tmp1[54]+z[91]*tmp1[55]+z[8]*tmp1[56]+z[22]*tmp1[57]+z[36]*tmp1[58]+z[50]*tmp1[59]+z[64]*tmp1[60]
                    +z[78]*tmp1[61]+z[92]*tmp1[62]+z[9]*tmp1[63]+z[23]*tmp1[64]+z[37]*tmp1[65]+z[51]*tmp1[66]+z[65]*tmp1[67]+z[79]*tmp1[68]+z[93]*tmp1[69]+z[10]*tmp1[70]
                    +z[24]*tmp1[71]+z[38]*tmp1[72]+z[52]*tmp1[73]+z[66]*tmp1[74]+z[80]*tmp1[75]+z[94]*tmp1[76]+z[11]*tmp1[77]+z[25]*tmp1[78]+z[39]*tmp1[79]+z[53]*tmp1[80]
                    +z[67]*tmp1[81]+z[81]*tmp1[82]+z[95]*tmp1[83]+z[12]*tmp1[84]+z[26]*tmp1[85]+z[40]*tmp1[86]+z[54]*tmp1[87]+z[68]*tmp1[88]+z[82]*tmp1[89]+z[96]*tmp1[90]
                    +z[13]*tmp1[91]+z[27]*tmp1[92]+z[41]*tmp1[93]+z[55]*tmp1[94]+z[69]*tmp1[95]+z[83]*tmp1[96];
                    tab[nb_tmp3+15*nb3]=tab[nb_tmp3+15*nb3]+z[0]*tmp1[0]
                    +z[15]*tmp1[1]+z[30]*tmp1[2]+z[45]*tmp1[3]+z[60]*tmp1[4]+z[75]*tmp1[5]+z[90]*tmp1[6]+z[8]*tmp1[7]+z[23]*tmp1[8]+z[38]*tmp1[9]+z[53]*tmp1[10]
                    +z[68]*tmp1[11]+z[83]*tmp1[12]+z[1]*tmp1[13]+z[16]*tmp1[14]+z[31]*tmp1[15]+z[46]*tmp1[16]+z[61]*tmp1[17]+z[76]*tmp1[18]+z[91]*tmp1[19]+z[9]*tmp1[20]
                    +z[24]*tmp1[21]+z[39]*tmp1[22]+z[54]*tmp1[23]+z[69]*tmp1[24]+z[84]*tmp1[25]+z[2]*tmp1[26]+z[17]*tmp1[27]+z[32]*tmp1[28]+z[47]*tmp1[29]+z[62]*tmp1[30]
                    +z[77]*tmp1[31]+z[92]*tmp1[32]+z[10]*tmp1[33]+z[25]*tmp1[34]+z[40]*tmp1[35]+z[55]*tmp1[36]+z[70]*tmp1[37]+z[85]*tmp1[38]+z[3]*tmp1[39]+z[18]*tmp1[40]
                    +z[33]*tmp1[41]+z[48]*tmp1[42]+z[63]*tmp1[43]+z[78]*tmp1[44]+z[93]*tmp1[45]+z[11]*tmp1[46]+z[26]*tmp1[47]+z[41]*tmp1[48]+z[56]*tmp1[49]+z[71]*tmp1[50]
                    +z[86]*tmp1[51]+z[4]*tmp1[52]+z[19]*tmp1[53]+z[34]*tmp1[54]+z[49]*tmp1[55]+z[64]*tmp1[56]+z[79]*tmp1[57]+z[94]*tmp1[58]+z[12]*tmp1[59]+z[27]*tmp1[60]
                    +z[42]*tmp1[61]+z[57]*tmp1[62]+z[72]*tmp1[63]+z[87]*tmp1[64]+z[5]*tmp1[65]+z[20]*tmp1[66]+z[35]*tmp1[67]+z[50]*tmp1[68]+z[65]*tmp1[69]+z[80]*tmp1[70]
                    +z[95]*tmp1[71]+z[13]*tmp1[72]+z[28]*tmp1[73]+z[43]*tmp1[74]+z[58]*tmp1[75]+z[73]*tmp1[76]+z[88]*tmp1[77]+z[6]*tmp1[78]+z[21]*tmp1[79]+z[36]*tmp1[80]
                    +z[51]*tmp1[81]+z[66]*tmp1[82]+z[81]*tmp1[83]+z[96]*tmp1[84]+z[14]*tmp1[85]+z[29]*tmp1[86]+z[44]*tmp1[87]+z[59]*tmp1[88]+z[74]*tmp1[89]+z[89]*tmp1[90]
                    +z[7]*tmp1[91]+z[22]*tmp1[92]+z[37]*tmp1[93]+z[52]*tmp1[94]+z[67]*tmp1[95]+z[82]*tmp1[96];
                    tab[nb_tmp3+16*nb3]=tab[nb_tmp3+16*nb3]+z[0]*tmp1[0]
                    +z[16]*tmp1[1]+z[32]*tmp1[2]+z[48]*tmp1[3]+z[64]*tmp1[4]+z[80]*tmp1[5]+z[96]*tmp1[6]+z[15]*tmp1[7]+z[31]*tmp1[8]+z[47]*tmp1[9]+z[63]*tmp1[10]
                    +z[79]*tmp1[11]+z[95]*tmp1[12]+z[14]*tmp1[13]+z[30]*tmp1[14]+z[46]*tmp1[15]+z[62]*tmp1[16]+z[78]*tmp1[17]+z[94]*tmp1[18]+z[13]*tmp1[19]+z[29]*tmp1[20]
                    +z[45]*tmp1[21]+z[61]*tmp1[22]+z[77]*tmp1[23]+z[93]*tmp1[24]+z[12]*tmp1[25]+z[28]*tmp1[26]+z[44]*tmp1[27]+z[60]*tmp1[28]+z[76]*tmp1[29]+z[92]*tmp1[30]
                    +z[11]*tmp1[31]+z[27]*tmp1[32]+z[43]*tmp1[33]+z[59]*tmp1[34]+z[75]*tmp1[35]+z[91]*tmp1[36]+z[10]*tmp1[37]+z[26]*tmp1[38]+z[42]*tmp1[39]+z[58]*tmp1[40]
                    +z[74]*tmp1[41]+z[90]*tmp1[42]+z[9]*tmp1[43]+z[25]*tmp1[44]+z[41]*tmp1[45]+z[57]*tmp1[46]+z[73]*tmp1[47]+z[89]*tmp1[48]+z[8]*tmp1[49]+z[24]*tmp1[50]
                    +z[40]*tmp1[51]+z[56]*tmp1[52]+z[72]*tmp1[53]+z[88]*tmp1[54]+z[7]*tmp1[55]+z[23]*tmp1[56]+z[39]*tmp1[57]+z[55]*tmp1[58]+z[71]*tmp1[59]+z[87]*tmp1[60]
                    +z[6]*tmp1[61]+z[22]*tmp1[62]+z[38]*tmp1[63]+z[54]*tmp1[64]+z[70]*tmp1[65]+z[86]*tmp1[66]+z[5]*tmp1[67]+z[21]*tmp1[68]+z[37]*tmp1[69]+z[53]*tmp1[70]
                    +z[69]*tmp1[71]+z[85]*tmp1[72]+z[4]*tmp1[73]+z[20]*tmp1[74]+z[36]*tmp1[75]+z[52]*tmp1[76]+z[68]*tmp1[77]+z[84]*tmp1[78]+z[3]*tmp1[79]+z[19]*tmp1[80]
                    +z[35]*tmp1[81]+z[51]*tmp1[82]+z[67]*tmp1[83]+z[83]*tmp1[84]+z[2]*tmp1[85]+z[18]*tmp1[86]+z[34]*tmp1[87]+z[50]*tmp1[88]+z[66]*tmp1[89]+z[82]*tmp1[90]
                    +z[1]*tmp1[91]+z[17]*tmp1[92]+z[33]*tmp1[93]+z[49]*tmp1[94]+z[65]*tmp1[95]+z[81]*tmp1[96];
                    tab[nb_tmp3+17*nb3]=tab[nb_tmp3+17*nb3]+z[0]*tmp1[0]
                    +z[17]*tmp1[1]+z[34]*tmp1[2]+z[51]*tmp1[3]+z[68]*tmp1[4]+z[85]*tmp1[5]+z[5]*tmp1[6]+z[22]*tmp1[7]+z[39]*tmp1[8]+z[56]*tmp1[9]+z[73]*tmp1[10]
                    +z[90]*tmp1[11]+z[10]*tmp1[12]+z[27]*tmp1[13]+z[44]*tmp1[14]+z[61]*tmp1[15]+z[78]*tmp1[16]+z[95]*tmp1[17]+z[15]*tmp1[18]+z[32]*tmp1[19]+z[49]*tmp1[20]
                    +z[66]*tmp1[21]+z[83]*tmp1[22]+z[3]*tmp1[23]+z[20]*tmp1[24]+z[37]*tmp1[25]+z[54]*tmp1[26]+z[71]*tmp1[27]+z[88]*tmp1[28]+z[8]*tmp1[29]+z[25]*tmp1[30]
                    +z[42]*tmp1[31]+z[59]*tmp1[32]+z[76]*tmp1[33]+z[93]*tmp1[34]+z[13]*tmp1[35]+z[30]*tmp1[36]+z[47]*tmp1[37]+z[64]*tmp1[38]+z[81]*tmp1[39]+z[1]*tmp1[40]
                    +z[18]*tmp1[41]+z[35]*tmp1[42]+z[52]*tmp1[43]+z[69]*tmp1[44]+z[86]*tmp1[45]+z[6]*tmp1[46]+z[23]*tmp1[47]+z[40]*tmp1[48]+z[57]*tmp1[49]+z[74]*tmp1[50]
                    +z[91]*tmp1[51]+z[11]*tmp1[52]+z[28]*tmp1[53]+z[45]*tmp1[54]+z[62]*tmp1[55]+z[79]*tmp1[56]+z[96]*tmp1[57]+z[16]*tmp1[58]+z[33]*tmp1[59]+z[50]*tmp1[60]
                    +z[67]*tmp1[61]+z[84]*tmp1[62]+z[4]*tmp1[63]+z[21]*tmp1[64]+z[38]*tmp1[65]+z[55]*tmp1[66]+z[72]*tmp1[67]+z[89]*tmp1[68]+z[9]*tmp1[69]+z[26]*tmp1[70]
                    +z[43]*tmp1[71]+z[60]*tmp1[72]+z[77]*tmp1[73]+z[94]*tmp1[74]+z[14]*tmp1[75]+z[31]*tmp1[76]+z[48]*tmp1[77]+z[65]*tmp1[78]+z[82]*tmp1[79]+z[2]*tmp1[80]
                    +z[19]*tmp1[81]+z[36]*tmp1[82]+z[53]*tmp1[83]+z[70]*tmp1[84]+z[87]*tmp1[85]+z[7]*tmp1[86]+z[24]*tmp1[87]+z[41]*tmp1[88]+z[58]*tmp1[89]+z[75]*tmp1[90]
                    +z[92]*tmp1[91]+z[12]*tmp1[92]+z[29]*tmp1[93]+z[46]*tmp1[94]+z[63]*tmp1[95]+z[80]*tmp1[96];
                    tab[nb_tmp3+18*nb3]=tab[nb_tmp3+18*nb3]+z[0]*tmp1[0]
                    +z[18]*tmp1[1]+z[36]*tmp1[2]+z[54]*tmp1[3]+z[72]*tmp1[4]+z[90]*tmp1[5]+z[11]*tmp1[6]+z[29]*tmp1[7]+z[47]*tmp1[8]+z[65]*tmp1[9]+z[83]*tmp1[10]
                    +z[4]*tmp1[11]+z[22]*tmp1[12]+z[40]*tmp1[13]+z[58]*tmp1[14]+z[76]*tmp1[15]+z[94]*tmp1[16]+z[15]*tmp1[17]+z[33]*tmp1[18]+z[51]*tmp1[19]+z[69]*tmp1[20]
                    +z[87]*tmp1[21]+z[8]*tmp1[22]+z[26]*tmp1[23]+z[44]*tmp1[24]+z[62]*tmp1[25]+z[80]*tmp1[26]+z[1]*tmp1[27]+z[19]*tmp1[28]+z[37]*tmp1[29]+z[55]*tmp1[30]
                    +z[73]*tmp1[31]+z[91]*tmp1[32]+z[12]*tmp1[33]+z[30]*tmp1[34]+z[48]*tmp1[35]+z[66]*tmp1[36]+z[84]*tmp1[37]+z[5]*tmp1[38]+z[23]*tmp1[39]+z[41]*tmp1[40]
                    +z[59]*tmp1[41]+z[77]*tmp1[42]+z[95]*tmp1[43]+z[16]*tmp1[44]+z[34]*tmp1[45]+z[52]*tmp1[46]+z[70]*tmp1[47]+z[88]*tmp1[48]+z[9]*tmp1[49]+z[27]*tmp1[50]
                    +z[45]*tmp1[51]+z[63]*tmp1[52]+z[81]*tmp1[53]+z[2]*tmp1[54]+z[20]*tmp1[55]+z[38]*tmp1[56]+z[56]*tmp1[57]+z[74]*tmp1[58]+z[92]*tmp1[59]+z[13]*tmp1[60]
                    +z[31]*tmp1[61]+z[49]*tmp1[62]+z[67]*tmp1[63]+z[85]*tmp1[64]+z[6]*tmp1[65]+z[24]*tmp1[66]+z[42]*tmp1[67]+z[60]*tmp1[68]+z[78]*tmp1[69]+z[96]*tmp1[70]
                    +z[17]*tmp1[71]+z[35]*tmp1[72]+z[53]*tmp1[73]+z[71]*tmp1[74]+z[89]*tmp1[75]+z[10]*tmp1[76]+z[28]*tmp1[77]+z[46]*tmp1[78]+z[64]*tmp1[79]+z[82]*tmp1[80]
                    +z[3]*tmp1[81]+z[21]*tmp1[82]+z[39]*tmp1[83]+z[57]*tmp1[84]+z[75]*tmp1[85]+z[93]*tmp1[86]+z[14]*tmp1[87]+z[32]*tmp1[88]+z[50]*tmp1[89]+z[68]*tmp1[90]
                    +z[86]*tmp1[91]+z[7]*tmp1[92]+z[25]*tmp1[93]+z[43]*tmp1[94]+z[61]*tmp1[95]+z[79]*tmp1[96];
                    tab[nb_tmp3+19*nb3]=tab[nb_tmp3+19*nb3]+z[0]*tmp1[0]
                    +z[19]*tmp1[1]+z[38]*tmp1[2]+z[57]*tmp1[3]+z[76]*tmp1[4]+z[95]*tmp1[5]+z[17]*tmp1[6]+z[36]*tmp1[7]+z[55]*tmp1[8]+z[74]*tmp1[9]+z[93]*tmp1[10]
                    +z[15]*tmp1[11]+z[34]*tmp1[12]+z[53]*tmp1[13]+z[72]*tmp1[14]+z[91]*tmp1[15]+z[13]*tmp1[16]+z[32]*tmp1[17]+z[51]*tmp1[18]+z[70]*tmp1[19]+z[89]*tmp1[20]
                    +z[11]*tmp1[21]+z[30]*tmp1[22]+z[49]*tmp1[23]+z[68]*tmp1[24]+z[87]*tmp1[25]+z[9]*tmp1[26]+z[28]*tmp1[27]+z[47]*tmp1[28]+z[66]*tmp1[29]+z[85]*tmp1[30]
                    +z[7]*tmp1[31]+z[26]*tmp1[32]+z[45]*tmp1[33]+z[64]*tmp1[34]+z[83]*tmp1[35]+z[5]*tmp1[36]+z[24]*tmp1[37]+z[43]*tmp1[38]+z[62]*tmp1[39]+z[81]*tmp1[40]
                    +z[3]*tmp1[41]+z[22]*tmp1[42]+z[41]*tmp1[43]+z[60]*tmp1[44]+z[79]*tmp1[45]+z[1]*tmp1[46]+z[20]*tmp1[47]+z[39]*tmp1[48]+z[58]*tmp1[49]+z[77]*tmp1[50]
                    +z[96]*tmp1[51]+z[18]*tmp1[52]+z[37]*tmp1[53]+z[56]*tmp1[54]+z[75]*tmp1[55]+z[94]*tmp1[56]+z[16]*tmp1[57]+z[35]*tmp1[58]+z[54]*tmp1[59]+z[73]*tmp1[60]
                    +z[92]*tmp1[61]+z[14]*tmp1[62]+z[33]*tmp1[63]+z[52]*tmp1[64]+z[71]*tmp1[65]+z[90]*tmp1[66]+z[12]*tmp1[67]+z[31]*tmp1[68]+z[50]*tmp1[69]+z[69]*tmp1[70]
                    +z[88]*tmp1[71]+z[10]*tmp1[72]+z[29]*tmp1[73]+z[48]*tmp1[74]+z[67]*tmp1[75]+z[86]*tmp1[76]+z[8]*tmp1[77]+z[27]*tmp1[78]+z[46]*tmp1[79]+z[65]*tmp1[80]
                    +z[84]*tmp1[81]+z[6]*tmp1[82]+z[25]*tmp1[83]+z[44]*tmp1[84]+z[63]*tmp1[85]+z[82]*tmp1[86]+z[4]*tmp1[87]+z[23]*tmp1[88]+z[42]*tmp1[89]+z[61]*tmp1[90]
                    +z[80]*tmp1[91]+z[2]*tmp1[92]+z[21]*tmp1[93]+z[40]*tmp1[94]+z[59]*tmp1[95]+z[78]*tmp1[96];
                    tab[nb_tmp3+20*nb3]=tab[nb_tmp3+20*nb3]+z[0]*tmp1[0]
                    +z[20]*tmp1[1]+z[40]*tmp1[2]+z[60]*tmp1[3]+z[80]*tmp1[4]+z[3]*tmp1[5]+z[23]*tmp1[6]+z[43]*tmp1[7]+z[63]*tmp1[8]+z[83]*tmp1[9]+z[6]*tmp1[10]
                    +z[26]*tmp1[11]+z[46]*tmp1[12]+z[66]*tmp1[13]+z[86]*tmp1[14]+z[9]*tmp1[15]+z[29]*tmp1[16]+z[49]*tmp1[17]+z[69]*tmp1[18]+z[89]*tmp1[19]+z[12]*tmp1[20]
                    +z[32]*tmp1[21]+z[52]*tmp1[22]+z[72]*tmp1[23]+z[92]*tmp1[24]+z[15]*tmp1[25]+z[35]*tmp1[26]+z[55]*tmp1[27]+z[75]*tmp1[28]+z[95]*tmp1[29]+z[18]*tmp1[30]
                    +z[38]*tmp1[31]+z[58]*tmp1[32]+z[78]*tmp1[33]+z[1]*tmp1[34]+z[21]*tmp1[35]+z[41]*tmp1[36]+z[61]*tmp1[37]+z[81]*tmp1[38]+z[4]*tmp1[39]+z[24]*tmp1[40]
                    +z[44]*tmp1[41]+z[64]*tmp1[42]+z[84]*tmp1[43]+z[7]*tmp1[44]+z[27]*tmp1[45]+z[47]*tmp1[46]+z[67]*tmp1[47]+z[87]*tmp1[48]+z[10]*tmp1[49]+z[30]*tmp1[50]
                    +z[50]*tmp1[51]+z[70]*tmp1[52]+z[90]*tmp1[53]+z[13]*tmp1[54]+z[33]*tmp1[55]+z[53]*tmp1[56]+z[73]*tmp1[57]+z[93]*tmp1[58]+z[16]*tmp1[59]+z[36]*tmp1[60]
                    +z[56]*tmp1[61]+z[76]*tmp1[62]+z[96]*tmp1[63]+z[19]*tmp1[64]+z[39]*tmp1[65]+z[59]*tmp1[66]+z[79]*tmp1[67]+z[2]*tmp1[68]+z[22]*tmp1[69]+z[42]*tmp1[70]
                    +z[62]*tmp1[71]+z[82]*tmp1[72]+z[5]*tmp1[73]+z[25]*tmp1[74]+z[45]*tmp1[75]+z[65]*tmp1[76]+z[85]*tmp1[77]+z[8]*tmp1[78]+z[28]*tmp1[79]+z[48]*tmp1[80]
                    +z[68]*tmp1[81]+z[88]*tmp1[82]+z[11]*tmp1[83]+z[31]*tmp1[84]+z[51]*tmp1[85]+z[71]*tmp1[86]+z[91]*tmp1[87]+z[14]*tmp1[88]+z[34]*tmp1[89]+z[54]*tmp1[90]
                    +z[74]*tmp1[91]+z[94]*tmp1[92]+z[17]*tmp1[93]+z[37]*tmp1[94]+z[57]*tmp1[95]+z[77]*tmp1[96];
                    tab[nb_tmp3+21*nb3]=tab[nb_tmp3+21*nb3]+z[0]*tmp1[0]
                    +z[21]*tmp1[1]+z[42]*tmp1[2]+z[63]*tmp1[3]+z[84]*tmp1[4]+z[8]*tmp1[5]+z[29]*tmp1[6]+z[50]*tmp1[7]+z[71]*tmp1[8]+z[92]*tmp1[9]+z[16]*tmp1[10]
                    +z[37]*tmp1[11]+z[58]*tmp1[12]+z[79]*tmp1[13]+z[3]*tmp1[14]+z[24]*tmp1[15]+z[45]*tmp1[16]+z[66]*tmp1[17]+z[87]*tmp1[18]+z[11]*tmp1[19]+z[32]*tmp1[20]
                    +z[53]*tmp1[21]+z[74]*tmp1[22]+z[95]*tmp1[23]+z[19]*tmp1[24]+z[40]*tmp1[25]+z[61]*tmp1[26]+z[82]*tmp1[27]+z[6]*tmp1[28]+z[27]*tmp1[29]+z[48]*tmp1[30]
                    +z[69]*tmp1[31]+z[90]*tmp1[32]+z[14]*tmp1[33]+z[35]*tmp1[34]+z[56]*tmp1[35]+z[77]*tmp1[36]+z[1]*tmp1[37]+z[22]*tmp1[38]+z[43]*tmp1[39]+z[64]*tmp1[40]
                    +z[85]*tmp1[41]+z[9]*tmp1[42]+z[30]*tmp1[43]+z[51]*tmp1[44]+z[72]*tmp1[45]+z[93]*tmp1[46]+z[17]*tmp1[47]+z[38]*tmp1[48]+z[59]*tmp1[49]+z[80]*tmp1[50]
                    +z[4]*tmp1[51]+z[25]*tmp1[52]+z[46]*tmp1[53]+z[67]*tmp1[54]+z[88]*tmp1[55]+z[12]*tmp1[56]+z[33]*tmp1[57]+z[54]*tmp1[58]+z[75]*tmp1[59]+z[96]*tmp1[60]
                    +z[20]*tmp1[61]+z[41]*tmp1[62]+z[62]*tmp1[63]+z[83]*tmp1[64]+z[7]*tmp1[65]+z[28]*tmp1[66]+z[49]*tmp1[67]+z[70]*tmp1[68]+z[91]*tmp1[69]+z[15]*tmp1[70]
                    +z[36]*tmp1[71]+z[57]*tmp1[72]+z[78]*tmp1[73]+z[2]*tmp1[74]+z[23]*tmp1[75]+z[44]*tmp1[76]+z[65]*tmp1[77]+z[86]*tmp1[78]+z[10]*tmp1[79]+z[31]*tmp1[80]
                    +z[52]*tmp1[81]+z[73]*tmp1[82]+z[94]*tmp1[83]+z[18]*tmp1[84]+z[39]*tmp1[85]+z[60]*tmp1[86]+z[81]*tmp1[87]+z[5]*tmp1[88]+z[26]*tmp1[89]+z[47]*tmp1[90]
                    +z[68]*tmp1[91]+z[89]*tmp1[92]+z[13]*tmp1[93]+z[34]*tmp1[94]+z[55]*tmp1[95]+z[76]*tmp1[96];
                    tab[nb_tmp3+22*nb3]=tab[nb_tmp3+22*nb3]+z[0]*tmp1[0]
                    +z[22]*tmp1[1]+z[44]*tmp1[2]+z[66]*tmp1[3]+z[88]*tmp1[4]+z[13]*tmp1[5]+z[35]*tmp1[6]+z[57]*tmp1[7]+z[79]*tmp1[8]+z[4]*tmp1[9]+z[26]*tmp1[10]
                    +z[48]*tmp1[11]+z[70]*tmp1[12]+z[92]*tmp1[13]+z[17]*tmp1[14]+z[39]*tmp1[15]+z[61]*tmp1[16]+z[83]*tmp1[17]+z[8]*tmp1[18]+z[30]*tmp1[19]+z[52]*tmp1[20]
                    +z[74]*tmp1[21]+z[96]*tmp1[22]+z[21]*tmp1[23]+z[43]*tmp1[24]+z[65]*tmp1[25]+z[87]*tmp1[26]+z[12]*tmp1[27]+z[34]*tmp1[28]+z[56]*tmp1[29]+z[78]*tmp1[30]
                    +z[3]*tmp1[31]+z[25]*tmp1[32]+z[47]*tmp1[33]+z[69]*tmp1[34]+z[91]*tmp1[35]+z[16]*tmp1[36]+z[38]*tmp1[37]+z[60]*tmp1[38]+z[82]*tmp1[39]+z[7]*tmp1[40]
                    +z[29]*tmp1[41]+z[51]*tmp1[42]+z[73]*tmp1[43]+z[95]*tmp1[44]+z[20]*tmp1[45]+z[42]*tmp1[46]+z[64]*tmp1[47]+z[86]*tmp1[48]+z[11]*tmp1[49]+z[33]*tmp1[50]
                    +z[55]*tmp1[51]+z[77]*tmp1[52]+z[2]*tmp1[53]+z[24]*tmp1[54]+z[46]*tmp1[55]+z[68]*tmp1[56]+z[90]*tmp1[57]+z[15]*tmp1[58]+z[37]*tmp1[59]+z[59]*tmp1[60]
                    +z[81]*tmp1[61]+z[6]*tmp1[62]+z[28]*tmp1[63]+z[50]*tmp1[64]+z[72]*tmp1[65]+z[94]*tmp1[66]+z[19]*tmp1[67]+z[41]*tmp1[68]+z[63]*tmp1[69]+z[85]*tmp1[70]
                    +z[10]*tmp1[71]+z[32]*tmp1[72]+z[54]*tmp1[73]+z[76]*tmp1[74]+z[1]*tmp1[75]+z[23]*tmp1[76]+z[45]*tmp1[77]+z[67]*tmp1[78]+z[89]*tmp1[79]+z[14]*tmp1[80]
                    +z[36]*tmp1[81]+z[58]*tmp1[82]+z[80]*tmp1[83]+z[5]*tmp1[84]+z[27]*tmp1[85]+z[49]*tmp1[86]+z[71]*tmp1[87]+z[93]*tmp1[88]+z[18]*tmp1[89]+z[40]*tmp1[90]
                    +z[62]*tmp1[91]+z[84]*tmp1[92]+z[9]*tmp1[93]+z[31]*tmp1[94]+z[53]*tmp1[95]+z[75]*tmp1[96];
                    tab[nb_tmp3+23*nb3]=tab[nb_tmp3+23*nb3]+z[0]*tmp1[0]
                    +z[23]*tmp1[1]+z[46]*tmp1[2]+z[69]*tmp1[3]+z[92]*tmp1[4]+z[18]*tmp1[5]+z[41]*tmp1[6]+z[64]*tmp1[7]+z[87]*tmp1[8]+z[13]*tmp1[9]+z[36]*tmp1[10]
                    +z[59]*tmp1[11]+z[82]*tmp1[12]+z[8]*tmp1[13]+z[31]*tmp1[14]+z[54]*tmp1[15]+z[77]*tmp1[16]+z[3]*tmp1[17]+z[26]*tmp1[18]+z[49]*tmp1[19]+z[72]*tmp1[20]
                    +z[95]*tmp1[21]+z[21]*tmp1[22]+z[44]*tmp1[23]+z[67]*tmp1[24]+z[90]*tmp1[25]+z[16]*tmp1[26]+z[39]*tmp1[27]+z[62]*tmp1[28]+z[85]*tmp1[29]+z[11]*tmp1[30]
                    +z[34]*tmp1[31]+z[57]*tmp1[32]+z[80]*tmp1[33]+z[6]*tmp1[34]+z[29]*tmp1[35]+z[52]*tmp1[36]+z[75]*tmp1[37]+z[1]*tmp1[38]+z[24]*tmp1[39]+z[47]*tmp1[40]
                    +z[70]*tmp1[41]+z[93]*tmp1[42]+z[19]*tmp1[43]+z[42]*tmp1[44]+z[65]*tmp1[45]+z[88]*tmp1[46]+z[14]*tmp1[47]+z[37]*tmp1[48]+z[60]*tmp1[49]+z[83]*tmp1[50]
                    +z[9]*tmp1[51]+z[32]*tmp1[52]+z[55]*tmp1[53]+z[78]*tmp1[54]+z[4]*tmp1[55]+z[27]*tmp1[56]+z[50]*tmp1[57]+z[73]*tmp1[58]+z[96]*tmp1[59]+z[22]*tmp1[60]
                    +z[45]*tmp1[61]+z[68]*tmp1[62]+z[91]*tmp1[63]+z[17]*tmp1[64]+z[40]*tmp1[65]+z[63]*tmp1[66]+z[86]*tmp1[67]+z[12]*tmp1[68]+z[35]*tmp1[69]+z[58]*tmp1[70]
                    +z[81]*tmp1[71]+z[7]*tmp1[72]+z[30]*tmp1[73]+z[53]*tmp1[74]+z[76]*tmp1[75]+z[2]*tmp1[76]+z[25]*tmp1[77]+z[48]*tmp1[78]+z[71]*tmp1[79]+z[94]*tmp1[80]
                    +z[20]*tmp1[81]+z[43]*tmp1[82]+z[66]*tmp1[83]+z[89]*tmp1[84]+z[15]*tmp1[85]+z[38]*tmp1[86]+z[61]*tmp1[87]+z[84]*tmp1[88]+z[10]*tmp1[89]+z[33]*tmp1[90]
                    +z[56]*tmp1[91]+z[79]*tmp1[92]+z[5]*tmp1[93]+z[28]*tmp1[94]+z[51]*tmp1[95]+z[74]*tmp1[96];
                    tab[nb_tmp3+24*nb3]=tab[nb_tmp3+24*nb3]+z[0]*tmp1[0]
                    +z[24]*tmp1[1]+z[48]*tmp1[2]+z[72]*tmp1[3]+z[96]*tmp1[4]+z[23]*tmp1[5]+z[47]*tmp1[6]+z[71]*tmp1[7]+z[95]*tmp1[8]+z[22]*tmp1[9]+z[46]*tmp1[10]
                    +z[70]*tmp1[11]+z[94]*tmp1[12]+z[21]*tmp1[13]+z[45]*tmp1[14]+z[69]*tmp1[15]+z[93]*tmp1[16]+z[20]*tmp1[17]+z[44]*tmp1[18]+z[68]*tmp1[19]+z[92]*tmp1[20]
                    +z[19]*tmp1[21]+z[43]*tmp1[22]+z[67]*tmp1[23]+z[91]*tmp1[24]+z[18]*tmp1[25]+z[42]*tmp1[26]+z[66]*tmp1[27]+z[90]*tmp1[28]+z[17]*tmp1[29]+z[41]*tmp1[30]
                    +z[65]*tmp1[31]+z[89]*tmp1[32]+z[16]*tmp1[33]+z[40]*tmp1[34]+z[64]*tmp1[35]+z[88]*tmp1[36]+z[15]*tmp1[37]+z[39]*tmp1[38]+z[63]*tmp1[39]+z[87]*tmp1[40]
                    +z[14]*tmp1[41]+z[38]*tmp1[42]+z[62]*tmp1[43]+z[86]*tmp1[44]+z[13]*tmp1[45]+z[37]*tmp1[46]+z[61]*tmp1[47]+z[85]*tmp1[48]+z[12]*tmp1[49]+z[36]*tmp1[50]
                    +z[60]*tmp1[51]+z[84]*tmp1[52]+z[11]*tmp1[53]+z[35]*tmp1[54]+z[59]*tmp1[55]+z[83]*tmp1[56]+z[10]*tmp1[57]+z[34]*tmp1[58]+z[58]*tmp1[59]+z[82]*tmp1[60]
                    +z[9]*tmp1[61]+z[33]*tmp1[62]+z[57]*tmp1[63]+z[81]*tmp1[64]+z[8]*tmp1[65]+z[32]*tmp1[66]+z[56]*tmp1[67]+z[80]*tmp1[68]+z[7]*tmp1[69]+z[31]*tmp1[70]
                    +z[55]*tmp1[71]+z[79]*tmp1[72]+z[6]*tmp1[73]+z[30]*tmp1[74]+z[54]*tmp1[75]+z[78]*tmp1[76]+z[5]*tmp1[77]+z[29]*tmp1[78]+z[53]*tmp1[79]+z[77]*tmp1[80]
                    +z[4]*tmp1[81]+z[28]*tmp1[82]+z[52]*tmp1[83]+z[76]*tmp1[84]+z[3]*tmp1[85]+z[27]*tmp1[86]+z[51]*tmp1[87]+z[75]*tmp1[88]+z[2]*tmp1[89]+z[26]*tmp1[90]
                    +z[50]*tmp1[91]+z[74]*tmp1[92]+z[1]*tmp1[93]+z[25]*tmp1[94]+z[49]*tmp1[95]+z[73]*tmp1[96];
                    tab[nb_tmp3+25*nb3]=tab[nb_tmp3+25*nb3]+z[0]*tmp1[0]
                    +z[25]*tmp1[1]+z[50]*tmp1[2]+z[75]*tmp1[3]+z[3]*tmp1[4]+z[28]*tmp1[5]+z[53]*tmp1[6]+z[78]*tmp1[7]+z[6]*tmp1[8]+z[31]*tmp1[9]+z[56]*tmp1[10]
                    +z[81]*tmp1[11]+z[9]*tmp1[12]+z[34]*tmp1[13]+z[59]*tmp1[14]+z[84]*tmp1[15]+z[12]*tmp1[16]+z[37]*tmp1[17]+z[62]*tmp1[18]+z[87]*tmp1[19]+z[15]*tmp1[20]
                    +z[40]*tmp1[21]+z[65]*tmp1[22]+z[90]*tmp1[23]+z[18]*tmp1[24]+z[43]*tmp1[25]+z[68]*tmp1[26]+z[93]*tmp1[27]+z[21]*tmp1[28]+z[46]*tmp1[29]+z[71]*tmp1[30]
                    +z[96]*tmp1[31]+z[24]*tmp1[32]+z[49]*tmp1[33]+z[74]*tmp1[34]+z[2]*tmp1[35]+z[27]*tmp1[36]+z[52]*tmp1[37]+z[77]*tmp1[38]+z[5]*tmp1[39]+z[30]*tmp1[40]
                    +z[55]*tmp1[41]+z[80]*tmp1[42]+z[8]*tmp1[43]+z[33]*tmp1[44]+z[58]*tmp1[45]+z[83]*tmp1[46]+z[11]*tmp1[47]+z[36]*tmp1[48]+z[61]*tmp1[49]+z[86]*tmp1[50]
                    +z[14]*tmp1[51]+z[39]*tmp1[52]+z[64]*tmp1[53]+z[89]*tmp1[54]+z[17]*tmp1[55]+z[42]*tmp1[56]+z[67]*tmp1[57]+z[92]*tmp1[58]+z[20]*tmp1[59]+z[45]*tmp1[60]
                    +z[70]*tmp1[61]+z[95]*tmp1[62]+z[23]*tmp1[63]+z[48]*tmp1[64]+z[73]*tmp1[65]+z[1]*tmp1[66]+z[26]*tmp1[67]+z[51]*tmp1[68]+z[76]*tmp1[69]+z[4]*tmp1[70]
                    +z[29]*tmp1[71]+z[54]*tmp1[72]+z[79]*tmp1[73]+z[7]*tmp1[74]+z[32]*tmp1[75]+z[57]*tmp1[76]+z[82]*tmp1[77]+z[10]*tmp1[78]+z[35]*tmp1[79]+z[60]*tmp1[80]
                    +z[85]*tmp1[81]+z[13]*tmp1[82]+z[38]*tmp1[83]+z[63]*tmp1[84]+z[88]*tmp1[85]+z[16]*tmp1[86]+z[41]*tmp1[87]+z[66]*tmp1[88]+z[91]*tmp1[89]+z[19]*tmp1[90]
                    +z[44]*tmp1[91]+z[69]*tmp1[92]+z[94]*tmp1[93]+z[22]*tmp1[94]+z[47]*tmp1[95]+z[72]*tmp1[96];
                    tab[nb_tmp3+26*nb3]=tab[nb_tmp3+26*nb3]+z[0]*tmp1[0]
                    +z[26]*tmp1[1]+z[52]*tmp1[2]+z[78]*tmp1[3]+z[7]*tmp1[4]+z[33]*tmp1[5]+z[59]*tmp1[6]+z[85]*tmp1[7]+z[14]*tmp1[8]+z[40]*tmp1[9]+z[66]*tmp1[10]
                    +z[92]*tmp1[11]+z[21]*tmp1[12]+z[47]*tmp1[13]+z[73]*tmp1[14]+z[2]*tmp1[15]+z[28]*tmp1[16]+z[54]*tmp1[17]+z[80]*tmp1[18]+z[9]*tmp1[19]+z[35]*tmp1[20]
                    +z[61]*tmp1[21]+z[87]*tmp1[22]+z[16]*tmp1[23]+z[42]*tmp1[24]+z[68]*tmp1[25]+z[94]*tmp1[26]+z[23]*tmp1[27]+z[49]*tmp1[28]+z[75]*tmp1[29]+z[4]*tmp1[30]
                    +z[30]*tmp1[31]+z[56]*tmp1[32]+z[82]*tmp1[33]+z[11]*tmp1[34]+z[37]*tmp1[35]+z[63]*tmp1[36]+z[89]*tmp1[37]+z[18]*tmp1[38]+z[44]*tmp1[39]+z[70]*tmp1[40]
                    +z[96]*tmp1[41]+z[25]*tmp1[42]+z[51]*tmp1[43]+z[77]*tmp1[44]+z[6]*tmp1[45]+z[32]*tmp1[46]+z[58]*tmp1[47]+z[84]*tmp1[48]+z[13]*tmp1[49]+z[39]*tmp1[50]
                    +z[65]*tmp1[51]+z[91]*tmp1[52]+z[20]*tmp1[53]+z[46]*tmp1[54]+z[72]*tmp1[55]+z[1]*tmp1[56]+z[27]*tmp1[57]+z[53]*tmp1[58]+z[79]*tmp1[59]+z[8]*tmp1[60]
                    +z[34]*tmp1[61]+z[60]*tmp1[62]+z[86]*tmp1[63]+z[15]*tmp1[64]+z[41]*tmp1[65]+z[67]*tmp1[66]+z[93]*tmp1[67]+z[22]*tmp1[68]+z[48]*tmp1[69]+z[74]*tmp1[70]
                    +z[3]*tmp1[71]+z[29]*tmp1[72]+z[55]*tmp1[73]+z[81]*tmp1[74]+z[10]*tmp1[75]+z[36]*tmp1[76]+z[62]*tmp1[77]+z[88]*tmp1[78]+z[17]*tmp1[79]+z[43]*tmp1[80]
                    +z[69]*tmp1[81]+z[95]*tmp1[82]+z[24]*tmp1[83]+z[50]*tmp1[84]+z[76]*tmp1[85]+z[5]*tmp1[86]+z[31]*tmp1[87]+z[57]*tmp1[88]+z[83]*tmp1[89]+z[12]*tmp1[90]
                    +z[38]*tmp1[91]+z[64]*tmp1[92]+z[90]*tmp1[93]+z[19]*tmp1[94]+z[45]*tmp1[95]+z[71]*tmp1[96];
                    tab[nb_tmp3+27*nb3]=tab[nb_tmp3+27*nb3]+z[0]*tmp1[0]
                    +z[27]*tmp1[1]+z[54]*tmp1[2]+z[81]*tmp1[3]+z[11]*tmp1[4]+z[38]*tmp1[5]+z[65]*tmp1[6]+z[92]*tmp1[7]+z[22]*tmp1[8]+z[49]*tmp1[9]+z[76]*tmp1[10]
                    +z[6]*tmp1[11]+z[33]*tmp1[12]+z[60]*tmp1[13]+z[87]*tmp1[14]+z[17]*tmp1[15]+z[44]*tmp1[16]+z[71]*tmp1[17]+z[1]*tmp1[18]+z[28]*tmp1[19]+z[55]*tmp1[20]
                    +z[82]*tmp1[21]+z[12]*tmp1[22]+z[39]*tmp1[23]+z[66]*tmp1[24]+z[93]*tmp1[25]+z[23]*tmp1[26]+z[50]*tmp1[27]+z[77]*tmp1[28]+z[7]*tmp1[29]+z[34]*tmp1[30]
                    +z[61]*tmp1[31]+z[88]*tmp1[32]+z[18]*tmp1[33]+z[45]*tmp1[34]+z[72]*tmp1[35]+z[2]*tmp1[36]+z[29]*tmp1[37]+z[56]*tmp1[38]+z[83]*tmp1[39]+z[13]*tmp1[40]
                    +z[40]*tmp1[41]+z[67]*tmp1[42]+z[94]*tmp1[43]+z[24]*tmp1[44]+z[51]*tmp1[45]+z[78]*tmp1[46]+z[8]*tmp1[47]+z[35]*tmp1[48]+z[62]*tmp1[49]+z[89]*tmp1[50]
                    +z[19]*tmp1[51]+z[46]*tmp1[52]+z[73]*tmp1[53]+z[3]*tmp1[54]+z[30]*tmp1[55]+z[57]*tmp1[56]+z[84]*tmp1[57]+z[14]*tmp1[58]+z[41]*tmp1[59]+z[68]*tmp1[60]
                    +z[95]*tmp1[61]+z[25]*tmp1[62]+z[52]*tmp1[63]+z[79]*tmp1[64]+z[9]*tmp1[65]+z[36]*tmp1[66]+z[63]*tmp1[67]+z[90]*tmp1[68]+z[20]*tmp1[69]+z[47]*tmp1[70]
                    +z[74]*tmp1[71]+z[4]*tmp1[72]+z[31]*tmp1[73]+z[58]*tmp1[74]+z[85]*tmp1[75]+z[15]*tmp1[76]+z[42]*tmp1[77]+z[69]*tmp1[78]+z[96]*tmp1[79]+z[26]*tmp1[80]
                    +z[53]*tmp1[81]+z[80]*tmp1[82]+z[10]*tmp1[83]+z[37]*tmp1[84]+z[64]*tmp1[85]+z[91]*tmp1[86]+z[21]*tmp1[87]+z[48]*tmp1[88]+z[75]*tmp1[89]+z[5]*tmp1[90]
                    +z[32]*tmp1[91]+z[59]*tmp1[92]+z[86]*tmp1[93]+z[16]*tmp1[94]+z[43]*tmp1[95]+z[70]*tmp1[96];
                    tab[nb_tmp3+28*nb3]=tab[nb_tmp3+28*nb3]+z[0]*tmp1[0]
                    +z[28]*tmp1[1]+z[56]*tmp1[2]+z[84]*tmp1[3]+z[15]*tmp1[4]+z[43]*tmp1[5]+z[71]*tmp1[6]+z[2]*tmp1[7]+z[30]*tmp1[8]+z[58]*tmp1[9]+z[86]*tmp1[10]
                    +z[17]*tmp1[11]+z[45]*tmp1[12]+z[73]*tmp1[13]+z[4]*tmp1[14]+z[32]*tmp1[15]+z[60]*tmp1[16]+z[88]*tmp1[17]+z[19]*tmp1[18]+z[47]*tmp1[19]+z[75]*tmp1[20]
                    +z[6]*tmp1[21]+z[34]*tmp1[22]+z[62]*tmp1[23]+z[90]*tmp1[24]+z[21]*tmp1[25]+z[49]*tmp1[26]+z[77]*tmp1[27]+z[8]*tmp1[28]+z[36]*tmp1[29]+z[64]*tmp1[30]
                    +z[92]*tmp1[31]+z[23]*tmp1[32]+z[51]*tmp1[33]+z[79]*tmp1[34]+z[10]*tmp1[35]+z[38]*tmp1[36]+z[66]*tmp1[37]+z[94]*tmp1[38]+z[25]*tmp1[39]+z[53]*tmp1[40]
                    +z[81]*tmp1[41]+z[12]*tmp1[42]+z[40]*tmp1[43]+z[68]*tmp1[44]+z[96]*tmp1[45]+z[27]*tmp1[46]+z[55]*tmp1[47]+z[83]*tmp1[48]+z[14]*tmp1[49]+z[42]*tmp1[50]
                    +z[70]*tmp1[51]+z[1]*tmp1[52]+z[29]*tmp1[53]+z[57]*tmp1[54]+z[85]*tmp1[55]+z[16]*tmp1[56]+z[44]*tmp1[57]+z[72]*tmp1[58]+z[3]*tmp1[59]+z[31]*tmp1[60]
                    +z[59]*tmp1[61]+z[87]*tmp1[62]+z[18]*tmp1[63]+z[46]*tmp1[64]+z[74]*tmp1[65]+z[5]*tmp1[66]+z[33]*tmp1[67]+z[61]*tmp1[68]+z[89]*tmp1[69]+z[20]*tmp1[70]
                    +z[48]*tmp1[71]+z[76]*tmp1[72]+z[7]*tmp1[73]+z[35]*tmp1[74]+z[63]*tmp1[75]+z[91]*tmp1[76]+z[22]*tmp1[77]+z[50]*tmp1[78]+z[78]*tmp1[79]+z[9]*tmp1[80]
                    +z[37]*tmp1[81]+z[65]*tmp1[82]+z[93]*tmp1[83]+z[24]*tmp1[84]+z[52]*tmp1[85]+z[80]*tmp1[86]+z[11]*tmp1[87]+z[39]*tmp1[88]+z[67]*tmp1[89]+z[95]*tmp1[90]
                    +z[26]*tmp1[91]+z[54]*tmp1[92]+z[82]*tmp1[93]+z[13]*tmp1[94]+z[41]*tmp1[95]+z[69]*tmp1[96];
                    tab[nb_tmp3+29*nb3]=tab[nb_tmp3+29*nb3]+z[0]*tmp1[0]
                    +z[29]*tmp1[1]+z[58]*tmp1[2]+z[87]*tmp1[3]+z[19]*tmp1[4]+z[48]*tmp1[5]+z[77]*tmp1[6]+z[9]*tmp1[7]+z[38]*tmp1[8]+z[67]*tmp1[9]+z[96]*tmp1[10]
                    +z[28]*tmp1[11]+z[57]*tmp1[12]+z[86]*tmp1[13]+z[18]*tmp1[14]+z[47]*tmp1[15]+z[76]*tmp1[16]+z[8]*tmp1[17]+z[37]*tmp1[18]+z[66]*tmp1[19]+z[95]*tmp1[20]
                    +z[27]*tmp1[21]+z[56]*tmp1[22]+z[85]*tmp1[23]+z[17]*tmp1[24]+z[46]*tmp1[25]+z[75]*tmp1[26]+z[7]*tmp1[27]+z[36]*tmp1[28]+z[65]*tmp1[29]+z[94]*tmp1[30]
                    +z[26]*tmp1[31]+z[55]*tmp1[32]+z[84]*tmp1[33]+z[16]*tmp1[34]+z[45]*tmp1[35]+z[74]*tmp1[36]+z[6]*tmp1[37]+z[35]*tmp1[38]+z[64]*tmp1[39]+z[93]*tmp1[40]
                    +z[25]*tmp1[41]+z[54]*tmp1[42]+z[83]*tmp1[43]+z[15]*tmp1[44]+z[44]*tmp1[45]+z[73]*tmp1[46]+z[5]*tmp1[47]+z[34]*tmp1[48]+z[63]*tmp1[49]+z[92]*tmp1[50]
                    +z[24]*tmp1[51]+z[53]*tmp1[52]+z[82]*tmp1[53]+z[14]*tmp1[54]+z[43]*tmp1[55]+z[72]*tmp1[56]+z[4]*tmp1[57]+z[33]*tmp1[58]+z[62]*tmp1[59]+z[91]*tmp1[60]
                    +z[23]*tmp1[61]+z[52]*tmp1[62]+z[81]*tmp1[63]+z[13]*tmp1[64]+z[42]*tmp1[65]+z[71]*tmp1[66]+z[3]*tmp1[67]+z[32]*tmp1[68]+z[61]*tmp1[69]+z[90]*tmp1[70]
                    +z[22]*tmp1[71]+z[51]*tmp1[72]+z[80]*tmp1[73]+z[12]*tmp1[74]+z[41]*tmp1[75]+z[70]*tmp1[76]+z[2]*tmp1[77]+z[31]*tmp1[78]+z[60]*tmp1[79]+z[89]*tmp1[80]
                    +z[21]*tmp1[81]+z[50]*tmp1[82]+z[79]*tmp1[83]+z[11]*tmp1[84]+z[40]*tmp1[85]+z[69]*tmp1[86]+z[1]*tmp1[87]+z[30]*tmp1[88]+z[59]*tmp1[89]+z[88]*tmp1[90]
                    +z[20]*tmp1[91]+z[49]*tmp1[92]+z[78]*tmp1[93]+z[10]*tmp1[94]+z[39]*tmp1[95]+z[68]*tmp1[96];
                    tab[nb_tmp3+30*nb3]=tab[nb_tmp3+30*nb3]+z[0]*tmp1[0]
                    +z[30]*tmp1[1]+z[60]*tmp1[2]+z[90]*tmp1[3]+z[23]*tmp1[4]+z[53]*tmp1[5]+z[83]*tmp1[6]+z[16]*tmp1[7]+z[46]*tmp1[8]+z[76]*tmp1[9]+z[9]*tmp1[10]
                    +z[39]*tmp1[11]+z[69]*tmp1[12]+z[2]*tmp1[13]+z[32]*tmp1[14]+z[62]*tmp1[15]+z[92]*tmp1[16]+z[25]*tmp1[17]+z[55]*tmp1[18]+z[85]*tmp1[19]+z[18]*tmp1[20]
                    +z[48]*tmp1[21]+z[78]*tmp1[22]+z[11]*tmp1[23]+z[41]*tmp1[24]+z[71]*tmp1[25]+z[4]*tmp1[26]+z[34]*tmp1[27]+z[64]*tmp1[28]+z[94]*tmp1[29]+z[27]*tmp1[30]
                    +z[57]*tmp1[31]+z[87]*tmp1[32]+z[20]*tmp1[33]+z[50]*tmp1[34]+z[80]*tmp1[35]+z[13]*tmp1[36]+z[43]*tmp1[37]+z[73]*tmp1[38]+z[6]*tmp1[39]+z[36]*tmp1[40]
                    +z[66]*tmp1[41]+z[96]*tmp1[42]+z[29]*tmp1[43]+z[59]*tmp1[44]+z[89]*tmp1[45]+z[22]*tmp1[46]+z[52]*tmp1[47]+z[82]*tmp1[48]+z[15]*tmp1[49]+z[45]*tmp1[50]
                    +z[75]*tmp1[51]+z[8]*tmp1[52]+z[38]*tmp1[53]+z[68]*tmp1[54]+z[1]*tmp1[55]+z[31]*tmp1[56]+z[61]*tmp1[57]+z[91]*tmp1[58]+z[24]*tmp1[59]+z[54]*tmp1[60]
                    +z[84]*tmp1[61]+z[17]*tmp1[62]+z[47]*tmp1[63]+z[77]*tmp1[64]+z[10]*tmp1[65]+z[40]*tmp1[66]+z[70]*tmp1[67]+z[3]*tmp1[68]+z[33]*tmp1[69]+z[63]*tmp1[70]
                    +z[93]*tmp1[71]+z[26]*tmp1[72]+z[56]*tmp1[73]+z[86]*tmp1[74]+z[19]*tmp1[75]+z[49]*tmp1[76]+z[79]*tmp1[77]+z[12]*tmp1[78]+z[42]*tmp1[79]+z[72]*tmp1[80]
                    +z[5]*tmp1[81]+z[35]*tmp1[82]+z[65]*tmp1[83]+z[95]*tmp1[84]+z[28]*tmp1[85]+z[58]*tmp1[86]+z[88]*tmp1[87]+z[21]*tmp1[88]+z[51]*tmp1[89]+z[81]*tmp1[90]
                    +z[14]*tmp1[91]+z[44]*tmp1[92]+z[74]*tmp1[93]+z[7]*tmp1[94]+z[37]*tmp1[95]+z[67]*tmp1[96];
                    tab[nb_tmp3+31*nb3]=tab[nb_tmp3+31*nb3]+z[0]*tmp1[0]
                    +z[31]*tmp1[1]+z[62]*tmp1[2]+z[93]*tmp1[3]+z[27]*tmp1[4]+z[58]*tmp1[5]+z[89]*tmp1[6]+z[23]*tmp1[7]+z[54]*tmp1[8]+z[85]*tmp1[9]+z[19]*tmp1[10]
                    +z[50]*tmp1[11]+z[81]*tmp1[12]+z[15]*tmp1[13]+z[46]*tmp1[14]+z[77]*tmp1[15]+z[11]*tmp1[16]+z[42]*tmp1[17]+z[73]*tmp1[18]+z[7]*tmp1[19]+z[38]*tmp1[20]
                    +z[69]*tmp1[21]+z[3]*tmp1[22]+z[34]*tmp1[23]+z[65]*tmp1[24]+z[96]*tmp1[25]+z[30]*tmp1[26]+z[61]*tmp1[27]+z[92]*tmp1[28]+z[26]*tmp1[29]+z[57]*tmp1[30]
                    +z[88]*tmp1[31]+z[22]*tmp1[32]+z[53]*tmp1[33]+z[84]*tmp1[34]+z[18]*tmp1[35]+z[49]*tmp1[36]+z[80]*tmp1[37]+z[14]*tmp1[38]+z[45]*tmp1[39]+z[76]*tmp1[40]
                    +z[10]*tmp1[41]+z[41]*tmp1[42]+z[72]*tmp1[43]+z[6]*tmp1[44]+z[37]*tmp1[45]+z[68]*tmp1[46]+z[2]*tmp1[47]+z[33]*tmp1[48]+z[64]*tmp1[49]+z[95]*tmp1[50]
                    +z[29]*tmp1[51]+z[60]*tmp1[52]+z[91]*tmp1[53]+z[25]*tmp1[54]+z[56]*tmp1[55]+z[87]*tmp1[56]+z[21]*tmp1[57]+z[52]*tmp1[58]+z[83]*tmp1[59]+z[17]*tmp1[60]
                    +z[48]*tmp1[61]+z[79]*tmp1[62]+z[13]*tmp1[63]+z[44]*tmp1[64]+z[75]*tmp1[65]+z[9]*tmp1[66]+z[40]*tmp1[67]+z[71]*tmp1[68]+z[5]*tmp1[69]+z[36]*tmp1[70]
                    +z[67]*tmp1[71]+z[1]*tmp1[72]+z[32]*tmp1[73]+z[63]*tmp1[74]+z[94]*tmp1[75]+z[28]*tmp1[76]+z[59]*tmp1[77]+z[90]*tmp1[78]+z[24]*tmp1[79]+z[55]*tmp1[80]
                    +z[86]*tmp1[81]+z[20]*tmp1[82]+z[51]*tmp1[83]+z[82]*tmp1[84]+z[16]*tmp1[85]+z[47]*tmp1[86]+z[78]*tmp1[87]+z[12]*tmp1[88]+z[43]*tmp1[89]+z[74]*tmp1[90]
                    +z[8]*tmp1[91]+z[39]*tmp1[92]+z[70]*tmp1[93]+z[4]*tmp1[94]+z[35]*tmp1[95]+z[66]*tmp1[96];
                    tab[nb_tmp3+32*nb3]=tab[nb_tmp3+32*nb3]+z[0]*tmp1[0]
                    +z[32]*tmp1[1]+z[64]*tmp1[2]+z[96]*tmp1[3]+z[31]*tmp1[4]+z[63]*tmp1[5]+z[95]*tmp1[6]+z[30]*tmp1[7]+z[62]*tmp1[8]+z[94]*tmp1[9]+z[29]*tmp1[10]
                    +z[61]*tmp1[11]+z[93]*tmp1[12]+z[28]*tmp1[13]+z[60]*tmp1[14]+z[92]*tmp1[15]+z[27]*tmp1[16]+z[59]*tmp1[17]+z[91]*tmp1[18]+z[26]*tmp1[19]+z[58]*tmp1[20]
                    +z[90]*tmp1[21]+z[25]*tmp1[22]+z[57]*tmp1[23]+z[89]*tmp1[24]+z[24]*tmp1[25]+z[56]*tmp1[26]+z[88]*tmp1[27]+z[23]*tmp1[28]+z[55]*tmp1[29]+z[87]*tmp1[30]
                    +z[22]*tmp1[31]+z[54]*tmp1[32]+z[86]*tmp1[33]+z[21]*tmp1[34]+z[53]*tmp1[35]+z[85]*tmp1[36]+z[20]*tmp1[37]+z[52]*tmp1[38]+z[84]*tmp1[39]+z[19]*tmp1[40]
                    +z[51]*tmp1[41]+z[83]*tmp1[42]+z[18]*tmp1[43]+z[50]*tmp1[44]+z[82]*tmp1[45]+z[17]*tmp1[46]+z[49]*tmp1[47]+z[81]*tmp1[48]+z[16]*tmp1[49]+z[48]*tmp1[50]
                    +z[80]*tmp1[51]+z[15]*tmp1[52]+z[47]*tmp1[53]+z[79]*tmp1[54]+z[14]*tmp1[55]+z[46]*tmp1[56]+z[78]*tmp1[57]+z[13]*tmp1[58]+z[45]*tmp1[59]+z[77]*tmp1[60]
                    +z[12]*tmp1[61]+z[44]*tmp1[62]+z[76]*tmp1[63]+z[11]*tmp1[64]+z[43]*tmp1[65]+z[75]*tmp1[66]+z[10]*tmp1[67]+z[42]*tmp1[68]+z[74]*tmp1[69]+z[9]*tmp1[70]
                    +z[41]*tmp1[71]+z[73]*tmp1[72]+z[8]*tmp1[73]+z[40]*tmp1[74]+z[72]*tmp1[75]+z[7]*tmp1[76]+z[39]*tmp1[77]+z[71]*tmp1[78]+z[6]*tmp1[79]+z[38]*tmp1[80]
                    +z[70]*tmp1[81]+z[5]*tmp1[82]+z[37]*tmp1[83]+z[69]*tmp1[84]+z[4]*tmp1[85]+z[36]*tmp1[86]+z[68]*tmp1[87]+z[3]*tmp1[88]+z[35]*tmp1[89]+z[67]*tmp1[90]
                    +z[2]*tmp1[91]+z[34]*tmp1[92]+z[66]*tmp1[93]+z[1]*tmp1[94]+z[33]*tmp1[95]+z[65]*tmp1[96];
                    tab[nb_tmp3+33*nb3]=tab[nb_tmp3+33*nb3]+z[0]*tmp1[0]
                    +z[33]*tmp1[1]+z[66]*tmp1[2]+z[2]*tmp1[3]+z[35]*tmp1[4]+z[68]*tmp1[5]+z[4]*tmp1[6]+z[37]*tmp1[7]+z[70]*tmp1[8]+z[6]*tmp1[9]+z[39]*tmp1[10]
                    +z[72]*tmp1[11]+z[8]*tmp1[12]+z[41]*tmp1[13]+z[74]*tmp1[14]+z[10]*tmp1[15]+z[43]*tmp1[16]+z[76]*tmp1[17]+z[12]*tmp1[18]+z[45]*tmp1[19]+z[78]*tmp1[20]
                    +z[14]*tmp1[21]+z[47]*tmp1[22]+z[80]*tmp1[23]+z[16]*tmp1[24]+z[49]*tmp1[25]+z[82]*tmp1[26]+z[18]*tmp1[27]+z[51]*tmp1[28]+z[84]*tmp1[29]+z[20]*tmp1[30]
                    +z[53]*tmp1[31]+z[86]*tmp1[32]+z[22]*tmp1[33]+z[55]*tmp1[34]+z[88]*tmp1[35]+z[24]*tmp1[36]+z[57]*tmp1[37]+z[90]*tmp1[38]+z[26]*tmp1[39]+z[59]*tmp1[40]
                    +z[92]*tmp1[41]+z[28]*tmp1[42]+z[61]*tmp1[43]+z[94]*tmp1[44]+z[30]*tmp1[45]+z[63]*tmp1[46]+z[96]*tmp1[47]+z[32]*tmp1[48]+z[65]*tmp1[49]+z[1]*tmp1[50]
                    +z[34]*tmp1[51]+z[67]*tmp1[52]+z[3]*tmp1[53]+z[36]*tmp1[54]+z[69]*tmp1[55]+z[5]*tmp1[56]+z[38]*tmp1[57]+z[71]*tmp1[58]+z[7]*tmp1[59]+z[40]*tmp1[60]
                    +z[73]*tmp1[61]+z[9]*tmp1[62]+z[42]*tmp1[63]+z[75]*tmp1[64]+z[11]*tmp1[65]+z[44]*tmp1[66]+z[77]*tmp1[67]+z[13]*tmp1[68]+z[46]*tmp1[69]+z[79]*tmp1[70]
                    +z[15]*tmp1[71]+z[48]*tmp1[72]+z[81]*tmp1[73]+z[17]*tmp1[74]+z[50]*tmp1[75]+z[83]*tmp1[76]+z[19]*tmp1[77]+z[52]*tmp1[78]+z[85]*tmp1[79]+z[21]*tmp1[80]
                    +z[54]*tmp1[81]+z[87]*tmp1[82]+z[23]*tmp1[83]+z[56]*tmp1[84]+z[89]*tmp1[85]+z[25]*tmp1[86]+z[58]*tmp1[87]+z[91]*tmp1[88]+z[27]*tmp1[89]+z[60]*tmp1[90]
                    +z[93]*tmp1[91]+z[29]*tmp1[92]+z[62]*tmp1[93]+z[95]*tmp1[94]+z[31]*tmp1[95]+z[64]*tmp1[96];
                    tab[nb_tmp3+34*nb3]=tab[nb_tmp3+34*nb3]+z[0]*tmp1[0]
                    +z[34]*tmp1[1]+z[68]*tmp1[2]+z[5]*tmp1[3]+z[39]*tmp1[4]+z[73]*tmp1[5]+z[10]*tmp1[6]+z[44]*tmp1[7]+z[78]*tmp1[8]+z[15]*tmp1[9]+z[49]*tmp1[10]
                    +z[83]*tmp1[11]+z[20]*tmp1[12]+z[54]*tmp1[13]+z[88]*tmp1[14]+z[25]*tmp1[15]+z[59]*tmp1[16]+z[93]*tmp1[17]+z[30]*tmp1[18]+z[64]*tmp1[19]+z[1]*tmp1[20]
                    +z[35]*tmp1[21]+z[69]*tmp1[22]+z[6]*tmp1[23]+z[40]*tmp1[24]+z[74]*tmp1[25]+z[11]*tmp1[26]+z[45]*tmp1[27]+z[79]*tmp1[28]+z[16]*tmp1[29]+z[50]*tmp1[30]
                    +z[84]*tmp1[31]+z[21]*tmp1[32]+z[55]*tmp1[33]+z[89]*tmp1[34]+z[26]*tmp1[35]+z[60]*tmp1[36]+z[94]*tmp1[37]+z[31]*tmp1[38]+z[65]*tmp1[39]+z[2]*tmp1[40]
                    +z[36]*tmp1[41]+z[70]*tmp1[42]+z[7]*tmp1[43]+z[41]*tmp1[44]+z[75]*tmp1[45]+z[12]*tmp1[46]+z[46]*tmp1[47]+z[80]*tmp1[48]+z[17]*tmp1[49]+z[51]*tmp1[50]
                    +z[85]*tmp1[51]+z[22]*tmp1[52]+z[56]*tmp1[53]+z[90]*tmp1[54]+z[27]*tmp1[55]+z[61]*tmp1[56]+z[95]*tmp1[57]+z[32]*tmp1[58]+z[66]*tmp1[59]+z[3]*tmp1[60]
                    +z[37]*tmp1[61]+z[71]*tmp1[62]+z[8]*tmp1[63]+z[42]*tmp1[64]+z[76]*tmp1[65]+z[13]*tmp1[66]+z[47]*tmp1[67]+z[81]*tmp1[68]+z[18]*tmp1[69]+z[52]*tmp1[70]
                    +z[86]*tmp1[71]+z[23]*tmp1[72]+z[57]*tmp1[73]+z[91]*tmp1[74]+z[28]*tmp1[75]+z[62]*tmp1[76]+z[96]*tmp1[77]+z[33]*tmp1[78]+z[67]*tmp1[79]+z[4]*tmp1[80]
                    +z[38]*tmp1[81]+z[72]*tmp1[82]+z[9]*tmp1[83]+z[43]*tmp1[84]+z[77]*tmp1[85]+z[14]*tmp1[86]+z[48]*tmp1[87]+z[82]*tmp1[88]+z[19]*tmp1[89]+z[53]*tmp1[90]
                    +z[87]*tmp1[91]+z[24]*tmp1[92]+z[58]*tmp1[93]+z[92]*tmp1[94]+z[29]*tmp1[95]+z[63]*tmp1[96];
                    tab[nb_tmp3+35*nb3]=tab[nb_tmp3+35*nb3]+z[0]*tmp1[0]
                    +z[35]*tmp1[1]+z[70]*tmp1[2]+z[8]*tmp1[3]+z[43]*tmp1[4]+z[78]*tmp1[5]+z[16]*tmp1[6]+z[51]*tmp1[7]+z[86]*tmp1[8]+z[24]*tmp1[9]+z[59]*tmp1[10]
                    +z[94]*tmp1[11]+z[32]*tmp1[12]+z[67]*tmp1[13]+z[5]*tmp1[14]+z[40]*tmp1[15]+z[75]*tmp1[16]+z[13]*tmp1[17]+z[48]*tmp1[18]+z[83]*tmp1[19]+z[21]*tmp1[20]
                    +z[56]*tmp1[21]+z[91]*tmp1[22]+z[29]*tmp1[23]+z[64]*tmp1[24]+z[2]*tmp1[25]+z[37]*tmp1[26]+z[72]*tmp1[27]+z[10]*tmp1[28]+z[45]*tmp1[29]+z[80]*tmp1[30]
                    +z[18]*tmp1[31]+z[53]*tmp1[32]+z[88]*tmp1[33]+z[26]*tmp1[34]+z[61]*tmp1[35]+z[96]*tmp1[36]+z[34]*tmp1[37]+z[69]*tmp1[38]+z[7]*tmp1[39]+z[42]*tmp1[40]
                    +z[77]*tmp1[41]+z[15]*tmp1[42]+z[50]*tmp1[43]+z[85]*tmp1[44]+z[23]*tmp1[45]+z[58]*tmp1[46]+z[93]*tmp1[47]+z[31]*tmp1[48]+z[66]*tmp1[49]+z[4]*tmp1[50]
                    +z[39]*tmp1[51]+z[74]*tmp1[52]+z[12]*tmp1[53]+z[47]*tmp1[54]+z[82]*tmp1[55]+z[20]*tmp1[56]+z[55]*tmp1[57]+z[90]*tmp1[58]+z[28]*tmp1[59]+z[63]*tmp1[60]
                    +z[1]*tmp1[61]+z[36]*tmp1[62]+z[71]*tmp1[63]+z[9]*tmp1[64]+z[44]*tmp1[65]+z[79]*tmp1[66]+z[17]*tmp1[67]+z[52]*tmp1[68]+z[87]*tmp1[69]+z[25]*tmp1[70]
                    +z[60]*tmp1[71]+z[95]*tmp1[72]+z[33]*tmp1[73]+z[68]*tmp1[74]+z[6]*tmp1[75]+z[41]*tmp1[76]+z[76]*tmp1[77]+z[14]*tmp1[78]+z[49]*tmp1[79]+z[84]*tmp1[80]
                    +z[22]*tmp1[81]+z[57]*tmp1[82]+z[92]*tmp1[83]+z[30]*tmp1[84]+z[65]*tmp1[85]+z[3]*tmp1[86]+z[38]*tmp1[87]+z[73]*tmp1[88]+z[11]*tmp1[89]+z[46]*tmp1[90]
                    +z[81]*tmp1[91]+z[19]*tmp1[92]+z[54]*tmp1[93]+z[89]*tmp1[94]+z[27]*tmp1[95]+z[62]*tmp1[96];
                    tab[nb_tmp3+36*nb3]=tab[nb_tmp3+36*nb3]+z[0]*tmp1[0]
                    +z[36]*tmp1[1]+z[72]*tmp1[2]+z[11]*tmp1[3]+z[47]*tmp1[4]+z[83]*tmp1[5]+z[22]*tmp1[6]+z[58]*tmp1[7]+z[94]*tmp1[8]+z[33]*tmp1[9]+z[69]*tmp1[10]
                    +z[8]*tmp1[11]+z[44]*tmp1[12]+z[80]*tmp1[13]+z[19]*tmp1[14]+z[55]*tmp1[15]+z[91]*tmp1[16]+z[30]*tmp1[17]+z[66]*tmp1[18]+z[5]*tmp1[19]+z[41]*tmp1[20]
                    +z[77]*tmp1[21]+z[16]*tmp1[22]+z[52]*tmp1[23]+z[88]*tmp1[24]+z[27]*tmp1[25]+z[63]*tmp1[26]+z[2]*tmp1[27]+z[38]*tmp1[28]+z[74]*tmp1[29]+z[13]*tmp1[30]
                    +z[49]*tmp1[31]+z[85]*tmp1[32]+z[24]*tmp1[33]+z[60]*tmp1[34]+z[96]*tmp1[35]+z[35]*tmp1[36]+z[71]*tmp1[37]+z[10]*tmp1[38]+z[46]*tmp1[39]+z[82]*tmp1[40]
                    +z[21]*tmp1[41]+z[57]*tmp1[42]+z[93]*tmp1[43]+z[32]*tmp1[44]+z[68]*tmp1[45]+z[7]*tmp1[46]+z[43]*tmp1[47]+z[79]*tmp1[48]+z[18]*tmp1[49]+z[54]*tmp1[50]
                    +z[90]*tmp1[51]+z[29]*tmp1[52]+z[65]*tmp1[53]+z[4]*tmp1[54]+z[40]*tmp1[55]+z[76]*tmp1[56]+z[15]*tmp1[57]+z[51]*tmp1[58]+z[87]*tmp1[59]+z[26]*tmp1[60]
                    +z[62]*tmp1[61]+z[1]*tmp1[62]+z[37]*tmp1[63]+z[73]*tmp1[64]+z[12]*tmp1[65]+z[48]*tmp1[66]+z[84]*tmp1[67]+z[23]*tmp1[68]+z[59]*tmp1[69]+z[95]*tmp1[70]
                    +z[34]*tmp1[71]+z[70]*tmp1[72]+z[9]*tmp1[73]+z[45]*tmp1[74]+z[81]*tmp1[75]+z[20]*tmp1[76]+z[56]*tmp1[77]+z[92]*tmp1[78]+z[31]*tmp1[79]+z[67]*tmp1[80]
                    +z[6]*tmp1[81]+z[42]*tmp1[82]+z[78]*tmp1[83]+z[17]*tmp1[84]+z[53]*tmp1[85]+z[89]*tmp1[86]+z[28]*tmp1[87]+z[64]*tmp1[88]+z[3]*tmp1[89]+z[39]*tmp1[90]
                    +z[75]*tmp1[91]+z[14]*tmp1[92]+z[50]*tmp1[93]+z[86]*tmp1[94]+z[25]*tmp1[95]+z[61]*tmp1[96];
                    tab[nb_tmp3+37*nb3]=tab[nb_tmp3+37*nb3]+z[0]*tmp1[0]
                    +z[37]*tmp1[1]+z[74]*tmp1[2]+z[14]*tmp1[3]+z[51]*tmp1[4]+z[88]*tmp1[5]+z[28]*tmp1[6]+z[65]*tmp1[7]+z[5]*tmp1[8]+z[42]*tmp1[9]+z[79]*tmp1[10]
                    +z[19]*tmp1[11]+z[56]*tmp1[12]+z[93]*tmp1[13]+z[33]*tmp1[14]+z[70]*tmp1[15]+z[10]*tmp1[16]+z[47]*tmp1[17]+z[84]*tmp1[18]+z[24]*tmp1[19]+z[61]*tmp1[20]
                    +z[1]*tmp1[21]+z[38]*tmp1[22]+z[75]*tmp1[23]+z[15]*tmp1[24]+z[52]*tmp1[25]+z[89]*tmp1[26]+z[29]*tmp1[27]+z[66]*tmp1[28]+z[6]*tmp1[29]+z[43]*tmp1[30]
                    +z[80]*tmp1[31]+z[20]*tmp1[32]+z[57]*tmp1[33]+z[94]*tmp1[34]+z[34]*tmp1[35]+z[71]*tmp1[36]+z[11]*tmp1[37]+z[48]*tmp1[38]+z[85]*tmp1[39]+z[25]*tmp1[40]
                    +z[62]*tmp1[41]+z[2]*tmp1[42]+z[39]*tmp1[43]+z[76]*tmp1[44]+z[16]*tmp1[45]+z[53]*tmp1[46]+z[90]*tmp1[47]+z[30]*tmp1[48]+z[67]*tmp1[49]+z[7]*tmp1[50]
                    +z[44]*tmp1[51]+z[81]*tmp1[52]+z[21]*tmp1[53]+z[58]*tmp1[54]+z[95]*tmp1[55]+z[35]*tmp1[56]+z[72]*tmp1[57]+z[12]*tmp1[58]+z[49]*tmp1[59]+z[86]*tmp1[60]
                    +z[26]*tmp1[61]+z[63]*tmp1[62]+z[3]*tmp1[63]+z[40]*tmp1[64]+z[77]*tmp1[65]+z[17]*tmp1[66]+z[54]*tmp1[67]+z[91]*tmp1[68]+z[31]*tmp1[69]+z[68]*tmp1[70]
                    +z[8]*tmp1[71]+z[45]*tmp1[72]+z[82]*tmp1[73]+z[22]*tmp1[74]+z[59]*tmp1[75]+z[96]*tmp1[76]+z[36]*tmp1[77]+z[73]*tmp1[78]+z[13]*tmp1[79]+z[50]*tmp1[80]
                    +z[87]*tmp1[81]+z[27]*tmp1[82]+z[64]*tmp1[83]+z[4]*tmp1[84]+z[41]*tmp1[85]+z[78]*tmp1[86]+z[18]*tmp1[87]+z[55]*tmp1[88]+z[92]*tmp1[89]+z[32]*tmp1[90]
                    +z[69]*tmp1[91]+z[9]*tmp1[92]+z[46]*tmp1[93]+z[83]*tmp1[94]+z[23]*tmp1[95]+z[60]*tmp1[96];
                    tab[nb_tmp3+38*nb3]=tab[nb_tmp3+38*nb3]+z[0]*tmp1[0]
                    +z[38]*tmp1[1]+z[76]*tmp1[2]+z[17]*tmp1[3]+z[55]*tmp1[4]+z[93]*tmp1[5]+z[34]*tmp1[6]+z[72]*tmp1[7]+z[13]*tmp1[8]+z[51]*tmp1[9]+z[89]*tmp1[10]
                    +z[30]*tmp1[11]+z[68]*tmp1[12]+z[9]*tmp1[13]+z[47]*tmp1[14]+z[85]*tmp1[15]+z[26]*tmp1[16]+z[64]*tmp1[17]+z[5]*tmp1[18]+z[43]*tmp1[19]+z[81]*tmp1[20]
                    +z[22]*tmp1[21]+z[60]*tmp1[22]+z[1]*tmp1[23]+z[39]*tmp1[24]+z[77]*tmp1[25]+z[18]*tmp1[26]+z[56]*tmp1[27]+z[94]*tmp1[28]+z[35]*tmp1[29]+z[73]*tmp1[30]
                    +z[14]*tmp1[31]+z[52]*tmp1[32]+z[90]*tmp1[33]+z[31]*tmp1[34]+z[69]*tmp1[35]+z[10]*tmp1[36]+z[48]*tmp1[37]+z[86]*tmp1[38]+z[27]*tmp1[39]+z[65]*tmp1[40]
                    +z[6]*tmp1[41]+z[44]*tmp1[42]+z[82]*tmp1[43]+z[23]*tmp1[44]+z[61]*tmp1[45]+z[2]*tmp1[46]+z[40]*tmp1[47]+z[78]*tmp1[48]+z[19]*tmp1[49]+z[57]*tmp1[50]
                    +z[95]*tmp1[51]+z[36]*tmp1[52]+z[74]*tmp1[53]+z[15]*tmp1[54]+z[53]*tmp1[55]+z[91]*tmp1[56]+z[32]*tmp1[57]+z[70]*tmp1[58]+z[11]*tmp1[59]+z[49]*tmp1[60]
                    +z[87]*tmp1[61]+z[28]*tmp1[62]+z[66]*tmp1[63]+z[7]*tmp1[64]+z[45]*tmp1[65]+z[83]*tmp1[66]+z[24]*tmp1[67]+z[62]*tmp1[68]+z[3]*tmp1[69]+z[41]*tmp1[70]
                    +z[79]*tmp1[71]+z[20]*tmp1[72]+z[58]*tmp1[73]+z[96]*tmp1[74]+z[37]*tmp1[75]+z[75]*tmp1[76]+z[16]*tmp1[77]+z[54]*tmp1[78]+z[92]*tmp1[79]+z[33]*tmp1[80]
                    +z[71]*tmp1[81]+z[12]*tmp1[82]+z[50]*tmp1[83]+z[88]*tmp1[84]+z[29]*tmp1[85]+z[67]*tmp1[86]+z[8]*tmp1[87]+z[46]*tmp1[88]+z[84]*tmp1[89]+z[25]*tmp1[90]
                    +z[63]*tmp1[91]+z[4]*tmp1[92]+z[42]*tmp1[93]+z[80]*tmp1[94]+z[21]*tmp1[95]+z[59]*tmp1[96];
                    tab[nb_tmp3+39*nb3]=tab[nb_tmp3+39*nb3]+z[0]*tmp1[0]
                    +z[39]*tmp1[1]+z[78]*tmp1[2]+z[20]*tmp1[3]+z[59]*tmp1[4]+z[1]*tmp1[5]+z[40]*tmp1[6]+z[79]*tmp1[7]+z[21]*tmp1[8]+z[60]*tmp1[9]+z[2]*tmp1[10]
                    +z[41]*tmp1[11]+z[80]*tmp1[12]+z[22]*tmp1[13]+z[61]*tmp1[14]+z[3]*tmp1[15]+z[42]*tmp1[16]+z[81]*tmp1[17]+z[23]*tmp1[18]+z[62]*tmp1[19]+z[4]*tmp1[20]
                    +z[43]*tmp1[21]+z[82]*tmp1[22]+z[24]*tmp1[23]+z[63]*tmp1[24]+z[5]*tmp1[25]+z[44]*tmp1[26]+z[83]*tmp1[27]+z[25]*tmp1[28]+z[64]*tmp1[29]+z[6]*tmp1[30]
                    +z[45]*tmp1[31]+z[84]*tmp1[32]+z[26]*tmp1[33]+z[65]*tmp1[34]+z[7]*tmp1[35]+z[46]*tmp1[36]+z[85]*tmp1[37]+z[27]*tmp1[38]+z[66]*tmp1[39]+z[8]*tmp1[40]
                    +z[47]*tmp1[41]+z[86]*tmp1[42]+z[28]*tmp1[43]+z[67]*tmp1[44]+z[9]*tmp1[45]+z[48]*tmp1[46]+z[87]*tmp1[47]+z[29]*tmp1[48]+z[68]*tmp1[49]+z[10]*tmp1[50]
                    +z[49]*tmp1[51]+z[88]*tmp1[52]+z[30]*tmp1[53]+z[69]*tmp1[54]+z[11]*tmp1[55]+z[50]*tmp1[56]+z[89]*tmp1[57]+z[31]*tmp1[58]+z[70]*tmp1[59]+z[12]*tmp1[60]
                    +z[51]*tmp1[61]+z[90]*tmp1[62]+z[32]*tmp1[63]+z[71]*tmp1[64]+z[13]*tmp1[65]+z[52]*tmp1[66]+z[91]*tmp1[67]+z[33]*tmp1[68]+z[72]*tmp1[69]+z[14]*tmp1[70]
                    +z[53]*tmp1[71]+z[92]*tmp1[72]+z[34]*tmp1[73]+z[73]*tmp1[74]+z[15]*tmp1[75]+z[54]*tmp1[76]+z[93]*tmp1[77]+z[35]*tmp1[78]+z[74]*tmp1[79]+z[16]*tmp1[80]
                    +z[55]*tmp1[81]+z[94]*tmp1[82]+z[36]*tmp1[83]+z[75]*tmp1[84]+z[17]*tmp1[85]+z[56]*tmp1[86]+z[95]*tmp1[87]+z[37]*tmp1[88]+z[76]*tmp1[89]+z[18]*tmp1[90]
                    +z[57]*tmp1[91]+z[96]*tmp1[92]+z[38]*tmp1[93]+z[77]*tmp1[94]+z[19]*tmp1[95]+z[58]*tmp1[96];
                    tab[nb_tmp3+40*nb3]=tab[nb_tmp3+40*nb3]+z[0]*tmp1[0]
                    +z[40]*tmp1[1]+z[80]*tmp1[2]+z[23]*tmp1[3]+z[63]*tmp1[4]+z[6]*tmp1[5]+z[46]*tmp1[6]+z[86]*tmp1[7]+z[29]*tmp1[8]+z[69]*tmp1[9]+z[12]*tmp1[10]
                    +z[52]*tmp1[11]+z[92]*tmp1[12]+z[35]*tmp1[13]+z[75]*tmp1[14]+z[18]*tmp1[15]+z[58]*tmp1[16]+z[1]*tmp1[17]+z[41]*tmp1[18]+z[81]*tmp1[19]+z[24]*tmp1[20]
                    +z[64]*tmp1[21]+z[7]*tmp1[22]+z[47]*tmp1[23]+z[87]*tmp1[24]+z[30]*tmp1[25]+z[70]*tmp1[26]+z[13]*tmp1[27]+z[53]*tmp1[28]+z[93]*tmp1[29]+z[36]*tmp1[30]
                    +z[76]*tmp1[31]+z[19]*tmp1[32]+z[59]*tmp1[33]+z[2]*tmp1[34]+z[42]*tmp1[35]+z[82]*tmp1[36]+z[25]*tmp1[37]+z[65]*tmp1[38]+z[8]*tmp1[39]+z[48]*tmp1[40]
                    +z[88]*tmp1[41]+z[31]*tmp1[42]+z[71]*tmp1[43]+z[14]*tmp1[44]+z[54]*tmp1[45]+z[94]*tmp1[46]+z[37]*tmp1[47]+z[77]*tmp1[48]+z[20]*tmp1[49]+z[60]*tmp1[50]
                    +z[3]*tmp1[51]+z[43]*tmp1[52]+z[83]*tmp1[53]+z[26]*tmp1[54]+z[66]*tmp1[55]+z[9]*tmp1[56]+z[49]*tmp1[57]+z[89]*tmp1[58]+z[32]*tmp1[59]+z[72]*tmp1[60]
                    +z[15]*tmp1[61]+z[55]*tmp1[62]+z[95]*tmp1[63]+z[38]*tmp1[64]+z[78]*tmp1[65]+z[21]*tmp1[66]+z[61]*tmp1[67]+z[4]*tmp1[68]+z[44]*tmp1[69]+z[84]*tmp1[70]
                    +z[27]*tmp1[71]+z[67]*tmp1[72]+z[10]*tmp1[73]+z[50]*tmp1[74]+z[90]*tmp1[75]+z[33]*tmp1[76]+z[73]*tmp1[77]+z[16]*tmp1[78]+z[56]*tmp1[79]+z[96]*tmp1[80]
                    +z[39]*tmp1[81]+z[79]*tmp1[82]+z[22]*tmp1[83]+z[62]*tmp1[84]+z[5]*tmp1[85]+z[45]*tmp1[86]+z[85]*tmp1[87]+z[28]*tmp1[88]+z[68]*tmp1[89]+z[11]*tmp1[90]
                    +z[51]*tmp1[91]+z[91]*tmp1[92]+z[34]*tmp1[93]+z[74]*tmp1[94]+z[17]*tmp1[95]+z[57]*tmp1[96];
                    tab[nb_tmp3+41*nb3]=tab[nb_tmp3+41*nb3]+z[0]*tmp1[0]
                    +z[41]*tmp1[1]+z[82]*tmp1[2]+z[26]*tmp1[3]+z[67]*tmp1[4]+z[11]*tmp1[5]+z[52]*tmp1[6]+z[93]*tmp1[7]+z[37]*tmp1[8]+z[78]*tmp1[9]+z[22]*tmp1[10]
                    +z[63]*tmp1[11]+z[7]*tmp1[12]+z[48]*tmp1[13]+z[89]*tmp1[14]+z[33]*tmp1[15]+z[74]*tmp1[16]+z[18]*tmp1[17]+z[59]*tmp1[18]+z[3]*tmp1[19]+z[44]*tmp1[20]
                    +z[85]*tmp1[21]+z[29]*tmp1[22]+z[70]*tmp1[23]+z[14]*tmp1[24]+z[55]*tmp1[25]+z[96]*tmp1[26]+z[40]*tmp1[27]+z[81]*tmp1[28]+z[25]*tmp1[29]+z[66]*tmp1[30]
                    +z[10]*tmp1[31]+z[51]*tmp1[32]+z[92]*tmp1[33]+z[36]*tmp1[34]+z[77]*tmp1[35]+z[21]*tmp1[36]+z[62]*tmp1[37]+z[6]*tmp1[38]+z[47]*tmp1[39]+z[88]*tmp1[40]
                    +z[32]*tmp1[41]+z[73]*tmp1[42]+z[17]*tmp1[43]+z[58]*tmp1[44]+z[2]*tmp1[45]+z[43]*tmp1[46]+z[84]*tmp1[47]+z[28]*tmp1[48]+z[69]*tmp1[49]+z[13]*tmp1[50]
                    +z[54]*tmp1[51]+z[95]*tmp1[52]+z[39]*tmp1[53]+z[80]*tmp1[54]+z[24]*tmp1[55]+z[65]*tmp1[56]+z[9]*tmp1[57]+z[50]*tmp1[58]+z[91]*tmp1[59]+z[35]*tmp1[60]
                    +z[76]*tmp1[61]+z[20]*tmp1[62]+z[61]*tmp1[63]+z[5]*tmp1[64]+z[46]*tmp1[65]+z[87]*tmp1[66]+z[31]*tmp1[67]+z[72]*tmp1[68]+z[16]*tmp1[69]+z[57]*tmp1[70]
                    +z[1]*tmp1[71]+z[42]*tmp1[72]+z[83]*tmp1[73]+z[27]*tmp1[74]+z[68]*tmp1[75]+z[12]*tmp1[76]+z[53]*tmp1[77]+z[94]*tmp1[78]+z[38]*tmp1[79]+z[79]*tmp1[80]
                    +z[23]*tmp1[81]+z[64]*tmp1[82]+z[8]*tmp1[83]+z[49]*tmp1[84]+z[90]*tmp1[85]+z[34]*tmp1[86]+z[75]*tmp1[87]+z[19]*tmp1[88]+z[60]*tmp1[89]+z[4]*tmp1[90]
                    +z[45]*tmp1[91]+z[86]*tmp1[92]+z[30]*tmp1[93]+z[71]*tmp1[94]+z[15]*tmp1[95]+z[56]*tmp1[96];
                    tab[nb_tmp3+42*nb3]=tab[nb_tmp3+42*nb3]+z[0]*tmp1[0]
                    +z[42]*tmp1[1]+z[84]*tmp1[2]+z[29]*tmp1[3]+z[71]*tmp1[4]+z[16]*tmp1[5]+z[58]*tmp1[6]+z[3]*tmp1[7]+z[45]*tmp1[8]+z[87]*tmp1[9]+z[32]*tmp1[10]
                    +z[74]*tmp1[11]+z[19]*tmp1[12]+z[61]*tmp1[13]+z[6]*tmp1[14]+z[48]*tmp1[15]+z[90]*tmp1[16]+z[35]*tmp1[17]+z[77]*tmp1[18]+z[22]*tmp1[19]+z[64]*tmp1[20]
                    +z[9]*tmp1[21]+z[51]*tmp1[22]+z[93]*tmp1[23]+z[38]*tmp1[24]+z[80]*tmp1[25]+z[25]*tmp1[26]+z[67]*tmp1[27]+z[12]*tmp1[28]+z[54]*tmp1[29]+z[96]*tmp1[30]
                    +z[41]*tmp1[31]+z[83]*tmp1[32]+z[28]*tmp1[33]+z[70]*tmp1[34]+z[15]*tmp1[35]+z[57]*tmp1[36]+z[2]*tmp1[37]+z[44]*tmp1[38]+z[86]*tmp1[39]+z[31]*tmp1[40]
                    +z[73]*tmp1[41]+z[18]*tmp1[42]+z[60]*tmp1[43]+z[5]*tmp1[44]+z[47]*tmp1[45]+z[89]*tmp1[46]+z[34]*tmp1[47]+z[76]*tmp1[48]+z[21]*tmp1[49]+z[63]*tmp1[50]
                    +z[8]*tmp1[51]+z[50]*tmp1[52]+z[92]*tmp1[53]+z[37]*tmp1[54]+z[79]*tmp1[55]+z[24]*tmp1[56]+z[66]*tmp1[57]+z[11]*tmp1[58]+z[53]*tmp1[59]+z[95]*tmp1[60]
                    +z[40]*tmp1[61]+z[82]*tmp1[62]+z[27]*tmp1[63]+z[69]*tmp1[64]+z[14]*tmp1[65]+z[56]*tmp1[66]+z[1]*tmp1[67]+z[43]*tmp1[68]+z[85]*tmp1[69]+z[30]*tmp1[70]
                    +z[72]*tmp1[71]+z[17]*tmp1[72]+z[59]*tmp1[73]+z[4]*tmp1[74]+z[46]*tmp1[75]+z[88]*tmp1[76]+z[33]*tmp1[77]+z[75]*tmp1[78]+z[20]*tmp1[79]+z[62]*tmp1[80]
                    +z[7]*tmp1[81]+z[49]*tmp1[82]+z[91]*tmp1[83]+z[36]*tmp1[84]+z[78]*tmp1[85]+z[23]*tmp1[86]+z[65]*tmp1[87]+z[10]*tmp1[88]+z[52]*tmp1[89]+z[94]*tmp1[90]
                    +z[39]*tmp1[91]+z[81]*tmp1[92]+z[26]*tmp1[93]+z[68]*tmp1[94]+z[13]*tmp1[95]+z[55]*tmp1[96];
                    tab[nb_tmp3+43*nb3]=tab[nb_tmp3+43*nb3]+z[0]*tmp1[0]
                    +z[43]*tmp1[1]+z[86]*tmp1[2]+z[32]*tmp1[3]+z[75]*tmp1[4]+z[21]*tmp1[5]+z[64]*tmp1[6]+z[10]*tmp1[7]+z[53]*tmp1[8]+z[96]*tmp1[9]+z[42]*tmp1[10]
                    +z[85]*tmp1[11]+z[31]*tmp1[12]+z[74]*tmp1[13]+z[20]*tmp1[14]+z[63]*tmp1[15]+z[9]*tmp1[16]+z[52]*tmp1[17]+z[95]*tmp1[18]+z[41]*tmp1[19]+z[84]*tmp1[20]
                    +z[30]*tmp1[21]+z[73]*tmp1[22]+z[19]*tmp1[23]+z[62]*tmp1[24]+z[8]*tmp1[25]+z[51]*tmp1[26]+z[94]*tmp1[27]+z[40]*tmp1[28]+z[83]*tmp1[29]+z[29]*tmp1[30]
                    +z[72]*tmp1[31]+z[18]*tmp1[32]+z[61]*tmp1[33]+z[7]*tmp1[34]+z[50]*tmp1[35]+z[93]*tmp1[36]+z[39]*tmp1[37]+z[82]*tmp1[38]+z[28]*tmp1[39]+z[71]*tmp1[40]
                    +z[17]*tmp1[41]+z[60]*tmp1[42]+z[6]*tmp1[43]+z[49]*tmp1[44]+z[92]*tmp1[45]+z[38]*tmp1[46]+z[81]*tmp1[47]+z[27]*tmp1[48]+z[70]*tmp1[49]+z[16]*tmp1[50]
                    +z[59]*tmp1[51]+z[5]*tmp1[52]+z[48]*tmp1[53]+z[91]*tmp1[54]+z[37]*tmp1[55]+z[80]*tmp1[56]+z[26]*tmp1[57]+z[69]*tmp1[58]+z[15]*tmp1[59]+z[58]*tmp1[60]
                    +z[4]*tmp1[61]+z[47]*tmp1[62]+z[90]*tmp1[63]+z[36]*tmp1[64]+z[79]*tmp1[65]+z[25]*tmp1[66]+z[68]*tmp1[67]+z[14]*tmp1[68]+z[57]*tmp1[69]+z[3]*tmp1[70]
                    +z[46]*tmp1[71]+z[89]*tmp1[72]+z[35]*tmp1[73]+z[78]*tmp1[74]+z[24]*tmp1[75]+z[67]*tmp1[76]+z[13]*tmp1[77]+z[56]*tmp1[78]+z[2]*tmp1[79]+z[45]*tmp1[80]
                    +z[88]*tmp1[81]+z[34]*tmp1[82]+z[77]*tmp1[83]+z[23]*tmp1[84]+z[66]*tmp1[85]+z[12]*tmp1[86]+z[55]*tmp1[87]+z[1]*tmp1[88]+z[44]*tmp1[89]+z[87]*tmp1[90]
                    +z[33]*tmp1[91]+z[76]*tmp1[92]+z[22]*tmp1[93]+z[65]*tmp1[94]+z[11]*tmp1[95]+z[54]*tmp1[96];
                    tab[nb_tmp3+44*nb3]=tab[nb_tmp3+44*nb3]+z[0]*tmp1[0]
                    +z[44]*tmp1[1]+z[88]*tmp1[2]+z[35]*tmp1[3]+z[79]*tmp1[4]+z[26]*tmp1[5]+z[70]*tmp1[6]+z[17]*tmp1[7]+z[61]*tmp1[8]+z[8]*tmp1[9]+z[52]*tmp1[10]
                    +z[96]*tmp1[11]+z[43]*tmp1[12]+z[87]*tmp1[13]+z[34]*tmp1[14]+z[78]*tmp1[15]+z[25]*tmp1[16]+z[69]*tmp1[17]+z[16]*tmp1[18]+z[60]*tmp1[19]+z[7]*tmp1[20]
                    +z[51]*tmp1[21]+z[95]*tmp1[22]+z[42]*tmp1[23]+z[86]*tmp1[24]+z[33]*tmp1[25]+z[77]*tmp1[26]+z[24]*tmp1[27]+z[68]*tmp1[28]+z[15]*tmp1[29]+z[59]*tmp1[30]
                    +z[6]*tmp1[31]+z[50]*tmp1[32]+z[94]*tmp1[33]+z[41]*tmp1[34]+z[85]*tmp1[35]+z[32]*tmp1[36]+z[76]*tmp1[37]+z[23]*tmp1[38]+z[67]*tmp1[39]+z[14]*tmp1[40]
                    +z[58]*tmp1[41]+z[5]*tmp1[42]+z[49]*tmp1[43]+z[93]*tmp1[44]+z[40]*tmp1[45]+z[84]*tmp1[46]+z[31]*tmp1[47]+z[75]*tmp1[48]+z[22]*tmp1[49]+z[66]*tmp1[50]
                    +z[13]*tmp1[51]+z[57]*tmp1[52]+z[4]*tmp1[53]+z[48]*tmp1[54]+z[92]*tmp1[55]+z[39]*tmp1[56]+z[83]*tmp1[57]+z[30]*tmp1[58]+z[74]*tmp1[59]+z[21]*tmp1[60]
                    +z[65]*tmp1[61]+z[12]*tmp1[62]+z[56]*tmp1[63]+z[3]*tmp1[64]+z[47]*tmp1[65]+z[91]*tmp1[66]+z[38]*tmp1[67]+z[82]*tmp1[68]+z[29]*tmp1[69]+z[73]*tmp1[70]
                    +z[20]*tmp1[71]+z[64]*tmp1[72]+z[11]*tmp1[73]+z[55]*tmp1[74]+z[2]*tmp1[75]+z[46]*tmp1[76]+z[90]*tmp1[77]+z[37]*tmp1[78]+z[81]*tmp1[79]+z[28]*tmp1[80]
                    +z[72]*tmp1[81]+z[19]*tmp1[82]+z[63]*tmp1[83]+z[10]*tmp1[84]+z[54]*tmp1[85]+z[1]*tmp1[86]+z[45]*tmp1[87]+z[89]*tmp1[88]+z[36]*tmp1[89]+z[80]*tmp1[90]
                    +z[27]*tmp1[91]+z[71]*tmp1[92]+z[18]*tmp1[93]+z[62]*tmp1[94]+z[9]*tmp1[95]+z[53]*tmp1[96];
                    tab[nb_tmp3+45*nb3]=tab[nb_tmp3+45*nb3]+z[0]*tmp1[0]
                    +z[45]*tmp1[1]+z[90]*tmp1[2]+z[38]*tmp1[3]+z[83]*tmp1[4]+z[31]*tmp1[5]+z[76]*tmp1[6]+z[24]*tmp1[7]+z[69]*tmp1[8]+z[17]*tmp1[9]+z[62]*tmp1[10]
                    +z[10]*tmp1[11]+z[55]*tmp1[12]+z[3]*tmp1[13]+z[48]*tmp1[14]+z[93]*tmp1[15]+z[41]*tmp1[16]+z[86]*tmp1[17]+z[34]*tmp1[18]+z[79]*tmp1[19]+z[27]*tmp1[20]
                    +z[72]*tmp1[21]+z[20]*tmp1[22]+z[65]*tmp1[23]+z[13]*tmp1[24]+z[58]*tmp1[25]+z[6]*tmp1[26]+z[51]*tmp1[27]+z[96]*tmp1[28]+z[44]*tmp1[29]+z[89]*tmp1[30]
                    +z[37]*tmp1[31]+z[82]*tmp1[32]+z[30]*tmp1[33]+z[75]*tmp1[34]+z[23]*tmp1[35]+z[68]*tmp1[36]+z[16]*tmp1[37]+z[61]*tmp1[38]+z[9]*tmp1[39]+z[54]*tmp1[40]
                    +z[2]*tmp1[41]+z[47]*tmp1[42]+z[92]*tmp1[43]+z[40]*tmp1[44]+z[85]*tmp1[45]+z[33]*tmp1[46]+z[78]*tmp1[47]+z[26]*tmp1[48]+z[71]*tmp1[49]+z[19]*tmp1[50]
                    +z[64]*tmp1[51]+z[12]*tmp1[52]+z[57]*tmp1[53]+z[5]*tmp1[54]+z[50]*tmp1[55]+z[95]*tmp1[56]+z[43]*tmp1[57]+z[88]*tmp1[58]+z[36]*tmp1[59]+z[81]*tmp1[60]
                    +z[29]*tmp1[61]+z[74]*tmp1[62]+z[22]*tmp1[63]+z[67]*tmp1[64]+z[15]*tmp1[65]+z[60]*tmp1[66]+z[8]*tmp1[67]+z[53]*tmp1[68]+z[1]*tmp1[69]+z[46]*tmp1[70]
                    +z[91]*tmp1[71]+z[39]*tmp1[72]+z[84]*tmp1[73]+z[32]*tmp1[74]+z[77]*tmp1[75]+z[25]*tmp1[76]+z[70]*tmp1[77]+z[18]*tmp1[78]+z[63]*tmp1[79]+z[11]*tmp1[80]
                    +z[56]*tmp1[81]+z[4]*tmp1[82]+z[49]*tmp1[83]+z[94]*tmp1[84]+z[42]*tmp1[85]+z[87]*tmp1[86]+z[35]*tmp1[87]+z[80]*tmp1[88]+z[28]*tmp1[89]+z[73]*tmp1[90]
                    +z[21]*tmp1[91]+z[66]*tmp1[92]+z[14]*tmp1[93]+z[59]*tmp1[94]+z[7]*tmp1[95]+z[52]*tmp1[96];
                    tab[nb_tmp3+46*nb3]=tab[nb_tmp3+46*nb3]+z[0]*tmp1[0]
                    +z[46]*tmp1[1]+z[92]*tmp1[2]+z[41]*tmp1[3]+z[87]*tmp1[4]+z[36]*tmp1[5]+z[82]*tmp1[6]+z[31]*tmp1[7]+z[77]*tmp1[8]+z[26]*tmp1[9]+z[72]*tmp1[10]
                    +z[21]*tmp1[11]+z[67]*tmp1[12]+z[16]*tmp1[13]+z[62]*tmp1[14]+z[11]*tmp1[15]+z[57]*tmp1[16]+z[6]*tmp1[17]+z[52]*tmp1[18]+z[1]*tmp1[19]+z[47]*tmp1[20]
                    +z[93]*tmp1[21]+z[42]*tmp1[22]+z[88]*tmp1[23]+z[37]*tmp1[24]+z[83]*tmp1[25]+z[32]*tmp1[26]+z[78]*tmp1[27]+z[27]*tmp1[28]+z[73]*tmp1[29]+z[22]*tmp1[30]
                    +z[68]*tmp1[31]+z[17]*tmp1[32]+z[63]*tmp1[33]+z[12]*tmp1[34]+z[58]*tmp1[35]+z[7]*tmp1[36]+z[53]*tmp1[37]+z[2]*tmp1[38]+z[48]*tmp1[39]+z[94]*tmp1[40]
                    +z[43]*tmp1[41]+z[89]*tmp1[42]+z[38]*tmp1[43]+z[84]*tmp1[44]+z[33]*tmp1[45]+z[79]*tmp1[46]+z[28]*tmp1[47]+z[74]*tmp1[48]+z[23]*tmp1[49]+z[69]*tmp1[50]
                    +z[18]*tmp1[51]+z[64]*tmp1[52]+z[13]*tmp1[53]+z[59]*tmp1[54]+z[8]*tmp1[55]+z[54]*tmp1[56]+z[3]*tmp1[57]+z[49]*tmp1[58]+z[95]*tmp1[59]+z[44]*tmp1[60]
                    +z[90]*tmp1[61]+z[39]*tmp1[62]+z[85]*tmp1[63]+z[34]*tmp1[64]+z[80]*tmp1[65]+z[29]*tmp1[66]+z[75]*tmp1[67]+z[24]*tmp1[68]+z[70]*tmp1[69]+z[19]*tmp1[70]
                    +z[65]*tmp1[71]+z[14]*tmp1[72]+z[60]*tmp1[73]+z[9]*tmp1[74]+z[55]*tmp1[75]+z[4]*tmp1[76]+z[50]*tmp1[77]+z[96]*tmp1[78]+z[45]*tmp1[79]+z[91]*tmp1[80]
                    +z[40]*tmp1[81]+z[86]*tmp1[82]+z[35]*tmp1[83]+z[81]*tmp1[84]+z[30]*tmp1[85]+z[76]*tmp1[86]+z[25]*tmp1[87]+z[71]*tmp1[88]+z[20]*tmp1[89]+z[66]*tmp1[90]
                    +z[15]*tmp1[91]+z[61]*tmp1[92]+z[10]*tmp1[93]+z[56]*tmp1[94]+z[5]*tmp1[95]+z[51]*tmp1[96];
                    tab[nb_tmp3+47*nb3]=tab[nb_tmp3+47*nb3]+z[0]*tmp1[0]
                    +z[47]*tmp1[1]+z[94]*tmp1[2]+z[44]*tmp1[3]+z[91]*tmp1[4]+z[41]*tmp1[5]+z[88]*tmp1[6]+z[38]*tmp1[7]+z[85]*tmp1[8]+z[35]*tmp1[9]+z[82]*tmp1[10]
                    +z[32]*tmp1[11]+z[79]*tmp1[12]+z[29]*tmp1[13]+z[76]*tmp1[14]+z[26]*tmp1[15]+z[73]*tmp1[16]+z[23]*tmp1[17]+z[70]*tmp1[18]+z[20]*tmp1[19]+z[67]*tmp1[20]
                    +z[17]*tmp1[21]+z[64]*tmp1[22]+z[14]*tmp1[23]+z[61]*tmp1[24]+z[11]*tmp1[25]+z[58]*tmp1[26]+z[8]*tmp1[27]+z[55]*tmp1[28]+z[5]*tmp1[29]+z[52]*tmp1[30]
                    +z[2]*tmp1[31]+z[49]*tmp1[32]+z[96]*tmp1[33]+z[46]*tmp1[34]+z[93]*tmp1[35]+z[43]*tmp1[36]+z[90]*tmp1[37]+z[40]*tmp1[38]+z[87]*tmp1[39]+z[37]*tmp1[40]
                    +z[84]*tmp1[41]+z[34]*tmp1[42]+z[81]*tmp1[43]+z[31]*tmp1[44]+z[78]*tmp1[45]+z[28]*tmp1[46]+z[75]*tmp1[47]+z[25]*tmp1[48]+z[72]*tmp1[49]+z[22]*tmp1[50]
                    +z[69]*tmp1[51]+z[19]*tmp1[52]+z[66]*tmp1[53]+z[16]*tmp1[54]+z[63]*tmp1[55]+z[13]*tmp1[56]+z[60]*tmp1[57]+z[10]*tmp1[58]+z[57]*tmp1[59]+z[7]*tmp1[60]
                    +z[54]*tmp1[61]+z[4]*tmp1[62]+z[51]*tmp1[63]+z[1]*tmp1[64]+z[48]*tmp1[65]+z[95]*tmp1[66]+z[45]*tmp1[67]+z[92]*tmp1[68]+z[42]*tmp1[69]+z[89]*tmp1[70]
                    +z[39]*tmp1[71]+z[86]*tmp1[72]+z[36]*tmp1[73]+z[83]*tmp1[74]+z[33]*tmp1[75]+z[80]*tmp1[76]+z[30]*tmp1[77]+z[77]*tmp1[78]+z[27]*tmp1[79]+z[74]*tmp1[80]
                    +z[24]*tmp1[81]+z[71]*tmp1[82]+z[21]*tmp1[83]+z[68]*tmp1[84]+z[18]*tmp1[85]+z[65]*tmp1[86]+z[15]*tmp1[87]+z[62]*tmp1[88]+z[12]*tmp1[89]+z[59]*tmp1[90]
                    +z[9]*tmp1[91]+z[56]*tmp1[92]+z[6]*tmp1[93]+z[53]*tmp1[94]+z[3]*tmp1[95]+z[50]*tmp1[96];
                    tab[nb_tmp3+48*nb3]=tab[nb_tmp3+48*nb3]+z[0]*tmp1[0]
                    +z[48]*tmp1[1]+z[96]*tmp1[2]+z[47]*tmp1[3]+z[95]*tmp1[4]+z[46]*tmp1[5]+z[94]*tmp1[6]+z[45]*tmp1[7]+z[93]*tmp1[8]+z[44]*tmp1[9]+z[92]*tmp1[10]
                    +z[43]*tmp1[11]+z[91]*tmp1[12]+z[42]*tmp1[13]+z[90]*tmp1[14]+z[41]*tmp1[15]+z[89]*tmp1[16]+z[40]*tmp1[17]+z[88]*tmp1[18]+z[39]*tmp1[19]+z[87]*tmp1[20]
                    +z[38]*tmp1[21]+z[86]*tmp1[22]+z[37]*tmp1[23]+z[85]*tmp1[24]+z[36]*tmp1[25]+z[84]*tmp1[26]+z[35]*tmp1[27]+z[83]*tmp1[28]+z[34]*tmp1[29]+z[82]*tmp1[30]
                    +z[33]*tmp1[31]+z[81]*tmp1[32]+z[32]*tmp1[33]+z[80]*tmp1[34]+z[31]*tmp1[35]+z[79]*tmp1[36]+z[30]*tmp1[37]+z[78]*tmp1[38]+z[29]*tmp1[39]+z[77]*tmp1[40]
                    +z[28]*tmp1[41]+z[76]*tmp1[42]+z[27]*tmp1[43]+z[75]*tmp1[44]+z[26]*tmp1[45]+z[74]*tmp1[46]+z[25]*tmp1[47]+z[73]*tmp1[48]+z[24]*tmp1[49]+z[72]*tmp1[50]
                    +z[23]*tmp1[51]+z[71]*tmp1[52]+z[22]*tmp1[53]+z[70]*tmp1[54]+z[21]*tmp1[55]+z[69]*tmp1[56]+z[20]*tmp1[57]+z[68]*tmp1[58]+z[19]*tmp1[59]+z[67]*tmp1[60]
                    +z[18]*tmp1[61]+z[66]*tmp1[62]+z[17]*tmp1[63]+z[65]*tmp1[64]+z[16]*tmp1[65]+z[64]*tmp1[66]+z[15]*tmp1[67]+z[63]*tmp1[68]+z[14]*tmp1[69]+z[62]*tmp1[70]
                    +z[13]*tmp1[71]+z[61]*tmp1[72]+z[12]*tmp1[73]+z[60]*tmp1[74]+z[11]*tmp1[75]+z[59]*tmp1[76]+z[10]*tmp1[77]+z[58]*tmp1[78]+z[9]*tmp1[79]+z[57]*tmp1[80]
                    +z[8]*tmp1[81]+z[56]*tmp1[82]+z[7]*tmp1[83]+z[55]*tmp1[84]+z[6]*tmp1[85]+z[54]*tmp1[86]+z[5]*tmp1[87]+z[53]*tmp1[88]+z[4]*tmp1[89]+z[52]*tmp1[90]
                    +z[3]*tmp1[91]+z[51]*tmp1[92]+z[2]*tmp1[93]+z[50]*tmp1[94]+z[1]*tmp1[95]+z[49]*tmp1[96];
                    tab[nb_tmp3+49*nb3]=tab[nb_tmp3+49*nb3]+z[0]*tmp1[0]
                    +z[49]*tmp1[1]+z[1]*tmp1[2]+z[50]*tmp1[3]+z[2]*tmp1[4]+z[51]*tmp1[5]+z[3]*tmp1[6]+z[52]*tmp1[7]+z[4]*tmp1[8]+z[53]*tmp1[9]+z[5]*tmp1[10]
                    +z[54]*tmp1[11]+z[6]*tmp1[12]+z[55]*tmp1[13]+z[7]*tmp1[14]+z[56]*tmp1[15]+z[8]*tmp1[16]+z[57]*tmp1[17]+z[9]*tmp1[18]+z[58]*tmp1[19]+z[10]*tmp1[20]
                    +z[59]*tmp1[21]+z[11]*tmp1[22]+z[60]*tmp1[23]+z[12]*tmp1[24]+z[61]*tmp1[25]+z[13]*tmp1[26]+z[62]*tmp1[27]+z[14]*tmp1[28]+z[63]*tmp1[29]+z[15]*tmp1[30]
                    +z[64]*tmp1[31]+z[16]*tmp1[32]+z[65]*tmp1[33]+z[17]*tmp1[34]+z[66]*tmp1[35]+z[18]*tmp1[36]+z[67]*tmp1[37]+z[19]*tmp1[38]+z[68]*tmp1[39]+z[20]*tmp1[40]
                    +z[69]*tmp1[41]+z[21]*tmp1[42]+z[70]*tmp1[43]+z[22]*tmp1[44]+z[71]*tmp1[45]+z[23]*tmp1[46]+z[72]*tmp1[47]+z[24]*tmp1[48]+z[73]*tmp1[49]+z[25]*tmp1[50]
                    +z[74]*tmp1[51]+z[26]*tmp1[52]+z[75]*tmp1[53]+z[27]*tmp1[54]+z[76]*tmp1[55]+z[28]*tmp1[56]+z[77]*tmp1[57]+z[29]*tmp1[58]+z[78]*tmp1[59]+z[30]*tmp1[60]
                    +z[79]*tmp1[61]+z[31]*tmp1[62]+z[80]*tmp1[63]+z[32]*tmp1[64]+z[81]*tmp1[65]+z[33]*tmp1[66]+z[82]*tmp1[67]+z[34]*tmp1[68]+z[83]*tmp1[69]+z[35]*tmp1[70]
                    +z[84]*tmp1[71]+z[36]*tmp1[72]+z[85]*tmp1[73]+z[37]*tmp1[74]+z[86]*tmp1[75]+z[38]*tmp1[76]+z[87]*tmp1[77]+z[39]*tmp1[78]+z[88]*tmp1[79]+z[40]*tmp1[80]
                    +z[89]*tmp1[81]+z[41]*tmp1[82]+z[90]*tmp1[83]+z[42]*tmp1[84]+z[91]*tmp1[85]+z[43]*tmp1[86]+z[92]*tmp1[87]+z[44]*tmp1[88]+z[93]*tmp1[89]+z[45]*tmp1[90]
                    +z[94]*tmp1[91]+z[46]*tmp1[92]+z[95]*tmp1[93]+z[47]*tmp1[94]+z[96]*tmp1[95]+z[48]*tmp1[96];
                    tab[nb_tmp3+50*nb3]=tab[nb_tmp3+50*nb3]+z[0]*tmp1[0]
                    +z[50]*tmp1[1]+z[3]*tmp1[2]+z[53]*tmp1[3]+z[6]*tmp1[4]+z[56]*tmp1[5]+z[9]*tmp1[6]+z[59]*tmp1[7]+z[12]*tmp1[8]+z[62]*tmp1[9]+z[15]*tmp1[10]
                    +z[65]*tmp1[11]+z[18]*tmp1[12]+z[68]*tmp1[13]+z[21]*tmp1[14]+z[71]*tmp1[15]+z[24]*tmp1[16]+z[74]*tmp1[17]+z[27]*tmp1[18]+z[77]*tmp1[19]+z[30]*tmp1[20]
                    +z[80]*tmp1[21]+z[33]*tmp1[22]+z[83]*tmp1[23]+z[36]*tmp1[24]+z[86]*tmp1[25]+z[39]*tmp1[26]+z[89]*tmp1[27]+z[42]*tmp1[28]+z[92]*tmp1[29]+z[45]*tmp1[30]
                    +z[95]*tmp1[31]+z[48]*tmp1[32]+z[1]*tmp1[33]+z[51]*tmp1[34]+z[4]*tmp1[35]+z[54]*tmp1[36]+z[7]*tmp1[37]+z[57]*tmp1[38]+z[10]*tmp1[39]+z[60]*tmp1[40]
                    +z[13]*tmp1[41]+z[63]*tmp1[42]+z[16]*tmp1[43]+z[66]*tmp1[44]+z[19]*tmp1[45]+z[69]*tmp1[46]+z[22]*tmp1[47]+z[72]*tmp1[48]+z[25]*tmp1[49]+z[75]*tmp1[50]
                    +z[28]*tmp1[51]+z[78]*tmp1[52]+z[31]*tmp1[53]+z[81]*tmp1[54]+z[34]*tmp1[55]+z[84]*tmp1[56]+z[37]*tmp1[57]+z[87]*tmp1[58]+z[40]*tmp1[59]+z[90]*tmp1[60]
                    +z[43]*tmp1[61]+z[93]*tmp1[62]+z[46]*tmp1[63]+z[96]*tmp1[64]+z[49]*tmp1[65]+z[2]*tmp1[66]+z[52]*tmp1[67]+z[5]*tmp1[68]+z[55]*tmp1[69]+z[8]*tmp1[70]
                    +z[58]*tmp1[71]+z[11]*tmp1[72]+z[61]*tmp1[73]+z[14]*tmp1[74]+z[64]*tmp1[75]+z[17]*tmp1[76]+z[67]*tmp1[77]+z[20]*tmp1[78]+z[70]*tmp1[79]+z[23]*tmp1[80]
                    +z[73]*tmp1[81]+z[26]*tmp1[82]+z[76]*tmp1[83]+z[29]*tmp1[84]+z[79]*tmp1[85]+z[32]*tmp1[86]+z[82]*tmp1[87]+z[35]*tmp1[88]+z[85]*tmp1[89]+z[38]*tmp1[90]
                    +z[88]*tmp1[91]+z[41]*tmp1[92]+z[91]*tmp1[93]+z[44]*tmp1[94]+z[94]*tmp1[95]+z[47]*tmp1[96];
                    tab[nb_tmp3+51*nb3]=tab[nb_tmp3+51*nb3]+z[0]*tmp1[0]
                    +z[51]*tmp1[1]+z[5]*tmp1[2]+z[56]*tmp1[3]+z[10]*tmp1[4]+z[61]*tmp1[5]+z[15]*tmp1[6]+z[66]*tmp1[7]+z[20]*tmp1[8]+z[71]*tmp1[9]+z[25]*tmp1[10]
                    +z[76]*tmp1[11]+z[30]*tmp1[12]+z[81]*tmp1[13]+z[35]*tmp1[14]+z[86]*tmp1[15]+z[40]*tmp1[16]+z[91]*tmp1[17]+z[45]*tmp1[18]+z[96]*tmp1[19]+z[50]*tmp1[20]
                    +z[4]*tmp1[21]+z[55]*tmp1[22]+z[9]*tmp1[23]+z[60]*tmp1[24]+z[14]*tmp1[25]+z[65]*tmp1[26]+z[19]*tmp1[27]+z[70]*tmp1[28]+z[24]*tmp1[29]+z[75]*tmp1[30]
                    +z[29]*tmp1[31]+z[80]*tmp1[32]+z[34]*tmp1[33]+z[85]*tmp1[34]+z[39]*tmp1[35]+z[90]*tmp1[36]+z[44]*tmp1[37]+z[95]*tmp1[38]+z[49]*tmp1[39]+z[3]*tmp1[40]
                    +z[54]*tmp1[41]+z[8]*tmp1[42]+z[59]*tmp1[43]+z[13]*tmp1[44]+z[64]*tmp1[45]+z[18]*tmp1[46]+z[69]*tmp1[47]+z[23]*tmp1[48]+z[74]*tmp1[49]+z[28]*tmp1[50]
                    +z[79]*tmp1[51]+z[33]*tmp1[52]+z[84]*tmp1[53]+z[38]*tmp1[54]+z[89]*tmp1[55]+z[43]*tmp1[56]+z[94]*tmp1[57]+z[48]*tmp1[58]+z[2]*tmp1[59]+z[53]*tmp1[60]
                    +z[7]*tmp1[61]+z[58]*tmp1[62]+z[12]*tmp1[63]+z[63]*tmp1[64]+z[17]*tmp1[65]+z[68]*tmp1[66]+z[22]*tmp1[67]+z[73]*tmp1[68]+z[27]*tmp1[69]+z[78]*tmp1[70]
                    +z[32]*tmp1[71]+z[83]*tmp1[72]+z[37]*tmp1[73]+z[88]*tmp1[74]+z[42]*tmp1[75]+z[93]*tmp1[76]+z[47]*tmp1[77]+z[1]*tmp1[78]+z[52]*tmp1[79]+z[6]*tmp1[80]
                    +z[57]*tmp1[81]+z[11]*tmp1[82]+z[62]*tmp1[83]+z[16]*tmp1[84]+z[67]*tmp1[85]+z[21]*tmp1[86]+z[72]*tmp1[87]+z[26]*tmp1[88]+z[77]*tmp1[89]+z[31]*tmp1[90]
                    +z[82]*tmp1[91]+z[36]*tmp1[92]+z[87]*tmp1[93]+z[41]*tmp1[94]+z[92]*tmp1[95]+z[46]*tmp1[96];
                    tab[nb_tmp3+52*nb3]=tab[nb_tmp3+52*nb3]+z[0]*tmp1[0]
                    +z[52]*tmp1[1]+z[7]*tmp1[2]+z[59]*tmp1[3]+z[14]*tmp1[4]+z[66]*tmp1[5]+z[21]*tmp1[6]+z[73]*tmp1[7]+z[28]*tmp1[8]+z[80]*tmp1[9]+z[35]*tmp1[10]
                    +z[87]*tmp1[11]+z[42]*tmp1[12]+z[94]*tmp1[13]+z[49]*tmp1[14]+z[4]*tmp1[15]+z[56]*tmp1[16]+z[11]*tmp1[17]+z[63]*tmp1[18]+z[18]*tmp1[19]+z[70]*tmp1[20]
                    +z[25]*tmp1[21]+z[77]*tmp1[22]+z[32]*tmp1[23]+z[84]*tmp1[24]+z[39]*tmp1[25]+z[91]*tmp1[26]+z[46]*tmp1[27]+z[1]*tmp1[28]+z[53]*tmp1[29]+z[8]*tmp1[30]
                    +z[60]*tmp1[31]+z[15]*tmp1[32]+z[67]*tmp1[33]+z[22]*tmp1[34]+z[74]*tmp1[35]+z[29]*tmp1[36]+z[81]*tmp1[37]+z[36]*tmp1[38]+z[88]*tmp1[39]+z[43]*tmp1[40]
                    +z[95]*tmp1[41]+z[50]*tmp1[42]+z[5]*tmp1[43]+z[57]*tmp1[44]+z[12]*tmp1[45]+z[64]*tmp1[46]+z[19]*tmp1[47]+z[71]*tmp1[48]+z[26]*tmp1[49]+z[78]*tmp1[50]
                    +z[33]*tmp1[51]+z[85]*tmp1[52]+z[40]*tmp1[53]+z[92]*tmp1[54]+z[47]*tmp1[55]+z[2]*tmp1[56]+z[54]*tmp1[57]+z[9]*tmp1[58]+z[61]*tmp1[59]+z[16]*tmp1[60]
                    +z[68]*tmp1[61]+z[23]*tmp1[62]+z[75]*tmp1[63]+z[30]*tmp1[64]+z[82]*tmp1[65]+z[37]*tmp1[66]+z[89]*tmp1[67]+z[44]*tmp1[68]+z[96]*tmp1[69]+z[51]*tmp1[70]
                    +z[6]*tmp1[71]+z[58]*tmp1[72]+z[13]*tmp1[73]+z[65]*tmp1[74]+z[20]*tmp1[75]+z[72]*tmp1[76]+z[27]*tmp1[77]+z[79]*tmp1[78]+z[34]*tmp1[79]+z[86]*tmp1[80]
                    +z[41]*tmp1[81]+z[93]*tmp1[82]+z[48]*tmp1[83]+z[3]*tmp1[84]+z[55]*tmp1[85]+z[10]*tmp1[86]+z[62]*tmp1[87]+z[17]*tmp1[88]+z[69]*tmp1[89]+z[24]*tmp1[90]
                    +z[76]*tmp1[91]+z[31]*tmp1[92]+z[83]*tmp1[93]+z[38]*tmp1[94]+z[90]*tmp1[95]+z[45]*tmp1[96];
                    tab[nb_tmp3+53*nb3]=tab[nb_tmp3+53*nb3]+z[0]*tmp1[0]
                    +z[53]*tmp1[1]+z[9]*tmp1[2]+z[62]*tmp1[3]+z[18]*tmp1[4]+z[71]*tmp1[5]+z[27]*tmp1[6]+z[80]*tmp1[7]+z[36]*tmp1[8]+z[89]*tmp1[9]+z[45]*tmp1[10]
                    +z[1]*tmp1[11]+z[54]*tmp1[12]+z[10]*tmp1[13]+z[63]*tmp1[14]+z[19]*tmp1[15]+z[72]*tmp1[16]+z[28]*tmp1[17]+z[81]*tmp1[18]+z[37]*tmp1[19]+z[90]*tmp1[20]
                    +z[46]*tmp1[21]+z[2]*tmp1[22]+z[55]*tmp1[23]+z[11]*tmp1[24]+z[64]*tmp1[25]+z[20]*tmp1[26]+z[73]*tmp1[27]+z[29]*tmp1[28]+z[82]*tmp1[29]+z[38]*tmp1[30]
                    +z[91]*tmp1[31]+z[47]*tmp1[32]+z[3]*tmp1[33]+z[56]*tmp1[34]+z[12]*tmp1[35]+z[65]*tmp1[36]+z[21]*tmp1[37]+z[74]*tmp1[38]+z[30]*tmp1[39]+z[83]*tmp1[40]
                    +z[39]*tmp1[41]+z[92]*tmp1[42]+z[48]*tmp1[43]+z[4]*tmp1[44]+z[57]*tmp1[45]+z[13]*tmp1[46]+z[66]*tmp1[47]+z[22]*tmp1[48]+z[75]*tmp1[49]+z[31]*tmp1[50]
                    +z[84]*tmp1[51]+z[40]*tmp1[52]+z[93]*tmp1[53]+z[49]*tmp1[54]+z[5]*tmp1[55]+z[58]*tmp1[56]+z[14]*tmp1[57]+z[67]*tmp1[58]+z[23]*tmp1[59]+z[76]*tmp1[60]
                    +z[32]*tmp1[61]+z[85]*tmp1[62]+z[41]*tmp1[63]+z[94]*tmp1[64]+z[50]*tmp1[65]+z[6]*tmp1[66]+z[59]*tmp1[67]+z[15]*tmp1[68]+z[68]*tmp1[69]+z[24]*tmp1[70]
                    +z[77]*tmp1[71]+z[33]*tmp1[72]+z[86]*tmp1[73]+z[42]*tmp1[74]+z[95]*tmp1[75]+z[51]*tmp1[76]+z[7]*tmp1[77]+z[60]*tmp1[78]+z[16]*tmp1[79]+z[69]*tmp1[80]
                    +z[25]*tmp1[81]+z[78]*tmp1[82]+z[34]*tmp1[83]+z[87]*tmp1[84]+z[43]*tmp1[85]+z[96]*tmp1[86]+z[52]*tmp1[87]+z[8]*tmp1[88]+z[61]*tmp1[89]+z[17]*tmp1[90]
                    +z[70]*tmp1[91]+z[26]*tmp1[92]+z[79]*tmp1[93]+z[35]*tmp1[94]+z[88]*tmp1[95]+z[44]*tmp1[96];
                    tab[nb_tmp3+54*nb3]=tab[nb_tmp3+54*nb3]+z[0]*tmp1[0]
                    +z[54]*tmp1[1]+z[11]*tmp1[2]+z[65]*tmp1[3]+z[22]*tmp1[4]+z[76]*tmp1[5]+z[33]*tmp1[6]+z[87]*tmp1[7]+z[44]*tmp1[8]+z[1]*tmp1[9]+z[55]*tmp1[10]
                    +z[12]*tmp1[11]+z[66]*tmp1[12]+z[23]*tmp1[13]+z[77]*tmp1[14]+z[34]*tmp1[15]+z[88]*tmp1[16]+z[45]*tmp1[17]+z[2]*tmp1[18]+z[56]*tmp1[19]+z[13]*tmp1[20]
                    +z[67]*tmp1[21]+z[24]*tmp1[22]+z[78]*tmp1[23]+z[35]*tmp1[24]+z[89]*tmp1[25]+z[46]*tmp1[26]+z[3]*tmp1[27]+z[57]*tmp1[28]+z[14]*tmp1[29]+z[68]*tmp1[30]
                    +z[25]*tmp1[31]+z[79]*tmp1[32]+z[36]*tmp1[33]+z[90]*tmp1[34]+z[47]*tmp1[35]+z[4]*tmp1[36]+z[58]*tmp1[37]+z[15]*tmp1[38]+z[69]*tmp1[39]+z[26]*tmp1[40]
                    +z[80]*tmp1[41]+z[37]*tmp1[42]+z[91]*tmp1[43]+z[48]*tmp1[44]+z[5]*tmp1[45]+z[59]*tmp1[46]+z[16]*tmp1[47]+z[70]*tmp1[48]+z[27]*tmp1[49]+z[81]*tmp1[50]
                    +z[38]*tmp1[51]+z[92]*tmp1[52]+z[49]*tmp1[53]+z[6]*tmp1[54]+z[60]*tmp1[55]+z[17]*tmp1[56]+z[71]*tmp1[57]+z[28]*tmp1[58]+z[82]*tmp1[59]+z[39]*tmp1[60]
                    +z[93]*tmp1[61]+z[50]*tmp1[62]+z[7]*tmp1[63]+z[61]*tmp1[64]+z[18]*tmp1[65]+z[72]*tmp1[66]+z[29]*tmp1[67]+z[83]*tmp1[68]+z[40]*tmp1[69]+z[94]*tmp1[70]
                    +z[51]*tmp1[71]+z[8]*tmp1[72]+z[62]*tmp1[73]+z[19]*tmp1[74]+z[73]*tmp1[75]+z[30]*tmp1[76]+z[84]*tmp1[77]+z[41]*tmp1[78]+z[95]*tmp1[79]+z[52]*tmp1[80]
                    +z[9]*tmp1[81]+z[63]*tmp1[82]+z[20]*tmp1[83]+z[74]*tmp1[84]+z[31]*tmp1[85]+z[85]*tmp1[86]+z[42]*tmp1[87]+z[96]*tmp1[88]+z[53]*tmp1[89]+z[10]*tmp1[90]
                    +z[64]*tmp1[91]+z[21]*tmp1[92]+z[75]*tmp1[93]+z[32]*tmp1[94]+z[86]*tmp1[95]+z[43]*tmp1[96];
                    tab[nb_tmp3+55*nb3]=tab[nb_tmp3+55*nb3]+z[0]*tmp1[0]
                    +z[55]*tmp1[1]+z[13]*tmp1[2]+z[68]*tmp1[3]+z[26]*tmp1[4]+z[81]*tmp1[5]+z[39]*tmp1[6]+z[94]*tmp1[7]+z[52]*tmp1[8]+z[10]*tmp1[9]+z[65]*tmp1[10]
                    +z[23]*tmp1[11]+z[78]*tmp1[12]+z[36]*tmp1[13]+z[91]*tmp1[14]+z[49]*tmp1[15]+z[7]*tmp1[16]+z[62]*tmp1[17]+z[20]*tmp1[18]+z[75]*tmp1[19]+z[33]*tmp1[20]
                    +z[88]*tmp1[21]+z[46]*tmp1[22]+z[4]*tmp1[23]+z[59]*tmp1[24]+z[17]*tmp1[25]+z[72]*tmp1[26]+z[30]*tmp1[27]+z[85]*tmp1[28]+z[43]*tmp1[29]+z[1]*tmp1[30]
                    +z[56]*tmp1[31]+z[14]*tmp1[32]+z[69]*tmp1[33]+z[27]*tmp1[34]+z[82]*tmp1[35]+z[40]*tmp1[36]+z[95]*tmp1[37]+z[53]*tmp1[38]+z[11]*tmp1[39]+z[66]*tmp1[40]
                    +z[24]*tmp1[41]+z[79]*tmp1[42]+z[37]*tmp1[43]+z[92]*tmp1[44]+z[50]*tmp1[45]+z[8]*tmp1[46]+z[63]*tmp1[47]+z[21]*tmp1[48]+z[76]*tmp1[49]+z[34]*tmp1[50]
                    +z[89]*tmp1[51]+z[47]*tmp1[52]+z[5]*tmp1[53]+z[60]*tmp1[54]+z[18]*tmp1[55]+z[73]*tmp1[56]+z[31]*tmp1[57]+z[86]*tmp1[58]+z[44]*tmp1[59]+z[2]*tmp1[60]
                    +z[57]*tmp1[61]+z[15]*tmp1[62]+z[70]*tmp1[63]+z[28]*tmp1[64]+z[83]*tmp1[65]+z[41]*tmp1[66]+z[96]*tmp1[67]+z[54]*tmp1[68]+z[12]*tmp1[69]+z[67]*tmp1[70]
                    +z[25]*tmp1[71]+z[80]*tmp1[72]+z[38]*tmp1[73]+z[93]*tmp1[74]+z[51]*tmp1[75]+z[9]*tmp1[76]+z[64]*tmp1[77]+z[22]*tmp1[78]+z[77]*tmp1[79]+z[35]*tmp1[80]
                    +z[90]*tmp1[81]+z[48]*tmp1[82]+z[6]*tmp1[83]+z[61]*tmp1[84]+z[19]*tmp1[85]+z[74]*tmp1[86]+z[32]*tmp1[87]+z[87]*tmp1[88]+z[45]*tmp1[89]+z[3]*tmp1[90]
                    +z[58]*tmp1[91]+z[16]*tmp1[92]+z[71]*tmp1[93]+z[29]*tmp1[94]+z[84]*tmp1[95]+z[42]*tmp1[96];
                    tab[nb_tmp3+56*nb3]=tab[nb_tmp3+56*nb3]+z[0]*tmp1[0]
                    +z[56]*tmp1[1]+z[15]*tmp1[2]+z[71]*tmp1[3]+z[30]*tmp1[4]+z[86]*tmp1[5]+z[45]*tmp1[6]+z[4]*tmp1[7]+z[60]*tmp1[8]+z[19]*tmp1[9]+z[75]*tmp1[10]
                    +z[34]*tmp1[11]+z[90]*tmp1[12]+z[49]*tmp1[13]+z[8]*tmp1[14]+z[64]*tmp1[15]+z[23]*tmp1[16]+z[79]*tmp1[17]+z[38]*tmp1[18]+z[94]*tmp1[19]+z[53]*tmp1[20]
                    +z[12]*tmp1[21]+z[68]*tmp1[22]+z[27]*tmp1[23]+z[83]*tmp1[24]+z[42]*tmp1[25]+z[1]*tmp1[26]+z[57]*tmp1[27]+z[16]*tmp1[28]+z[72]*tmp1[29]+z[31]*tmp1[30]
                    +z[87]*tmp1[31]+z[46]*tmp1[32]+z[5]*tmp1[33]+z[61]*tmp1[34]+z[20]*tmp1[35]+z[76]*tmp1[36]+z[35]*tmp1[37]+z[91]*tmp1[38]+z[50]*tmp1[39]+z[9]*tmp1[40]
                    +z[65]*tmp1[41]+z[24]*tmp1[42]+z[80]*tmp1[43]+z[39]*tmp1[44]+z[95]*tmp1[45]+z[54]*tmp1[46]+z[13]*tmp1[47]+z[69]*tmp1[48]+z[28]*tmp1[49]+z[84]*tmp1[50]
                    +z[43]*tmp1[51]+z[2]*tmp1[52]+z[58]*tmp1[53]+z[17]*tmp1[54]+z[73]*tmp1[55]+z[32]*tmp1[56]+z[88]*tmp1[57]+z[47]*tmp1[58]+z[6]*tmp1[59]+z[62]*tmp1[60]
                    +z[21]*tmp1[61]+z[77]*tmp1[62]+z[36]*tmp1[63]+z[92]*tmp1[64]+z[51]*tmp1[65]+z[10]*tmp1[66]+z[66]*tmp1[67]+z[25]*tmp1[68]+z[81]*tmp1[69]+z[40]*tmp1[70]
                    +z[96]*tmp1[71]+z[55]*tmp1[72]+z[14]*tmp1[73]+z[70]*tmp1[74]+z[29]*tmp1[75]+z[85]*tmp1[76]+z[44]*tmp1[77]+z[3]*tmp1[78]+z[59]*tmp1[79]+z[18]*tmp1[80]
                    +z[74]*tmp1[81]+z[33]*tmp1[82]+z[89]*tmp1[83]+z[48]*tmp1[84]+z[7]*tmp1[85]+z[63]*tmp1[86]+z[22]*tmp1[87]+z[78]*tmp1[88]+z[37]*tmp1[89]+z[93]*tmp1[90]
                    +z[52]*tmp1[91]+z[11]*tmp1[92]+z[67]*tmp1[93]+z[26]*tmp1[94]+z[82]*tmp1[95]+z[41]*tmp1[96];
                    tab[nb_tmp3+57*nb3]=tab[nb_tmp3+57*nb3]+z[0]*tmp1[0]
                    +z[57]*tmp1[1]+z[17]*tmp1[2]+z[74]*tmp1[3]+z[34]*tmp1[4]+z[91]*tmp1[5]+z[51]*tmp1[6]+z[11]*tmp1[7]+z[68]*tmp1[8]+z[28]*tmp1[9]+z[85]*tmp1[10]
                    +z[45]*tmp1[11]+z[5]*tmp1[12]+z[62]*tmp1[13]+z[22]*tmp1[14]+z[79]*tmp1[15]+z[39]*tmp1[16]+z[96]*tmp1[17]+z[56]*tmp1[18]+z[16]*tmp1[19]+z[73]*tmp1[20]
                    +z[33]*tmp1[21]+z[90]*tmp1[22]+z[50]*tmp1[23]+z[10]*tmp1[24]+z[67]*tmp1[25]+z[27]*tmp1[26]+z[84]*tmp1[27]+z[44]*tmp1[28]+z[4]*tmp1[29]+z[61]*tmp1[30]
                    +z[21]*tmp1[31]+z[78]*tmp1[32]+z[38]*tmp1[33]+z[95]*tmp1[34]+z[55]*tmp1[35]+z[15]*tmp1[36]+z[72]*tmp1[37]+z[32]*tmp1[38]+z[89]*tmp1[39]+z[49]*tmp1[40]
                    +z[9]*tmp1[41]+z[66]*tmp1[42]+z[26]*tmp1[43]+z[83]*tmp1[44]+z[43]*tmp1[45]+z[3]*tmp1[46]+z[60]*tmp1[47]+z[20]*tmp1[48]+z[77]*tmp1[49]+z[37]*tmp1[50]
                    +z[94]*tmp1[51]+z[54]*tmp1[52]+z[14]*tmp1[53]+z[71]*tmp1[54]+z[31]*tmp1[55]+z[88]*tmp1[56]+z[48]*tmp1[57]+z[8]*tmp1[58]+z[65]*tmp1[59]+z[25]*tmp1[60]
                    +z[82]*tmp1[61]+z[42]*tmp1[62]+z[2]*tmp1[63]+z[59]*tmp1[64]+z[19]*tmp1[65]+z[76]*tmp1[66]+z[36]*tmp1[67]+z[93]*tmp1[68]+z[53]*tmp1[69]+z[13]*tmp1[70]
                    +z[70]*tmp1[71]+z[30]*tmp1[72]+z[87]*tmp1[73]+z[47]*tmp1[74]+z[7]*tmp1[75]+z[64]*tmp1[76]+z[24]*tmp1[77]+z[81]*tmp1[78]+z[41]*tmp1[79]+z[1]*tmp1[80]
                    +z[58]*tmp1[81]+z[18]*tmp1[82]+z[75]*tmp1[83]+z[35]*tmp1[84]+z[92]*tmp1[85]+z[52]*tmp1[86]+z[12]*tmp1[87]+z[69]*tmp1[88]+z[29]*tmp1[89]+z[86]*tmp1[90]
                    +z[46]*tmp1[91]+z[6]*tmp1[92]+z[63]*tmp1[93]+z[23]*tmp1[94]+z[80]*tmp1[95]+z[40]*tmp1[96];
                    tab[nb_tmp3+58*nb3]=tab[nb_tmp3+58*nb3]+z[0]*tmp1[0]
                    +z[58]*tmp1[1]+z[19]*tmp1[2]+z[77]*tmp1[3]+z[38]*tmp1[4]+z[96]*tmp1[5]+z[57]*tmp1[6]+z[18]*tmp1[7]+z[76]*tmp1[8]+z[37]*tmp1[9]+z[95]*tmp1[10]
                    +z[56]*tmp1[11]+z[17]*tmp1[12]+z[75]*tmp1[13]+z[36]*tmp1[14]+z[94]*tmp1[15]+z[55]*tmp1[16]+z[16]*tmp1[17]+z[74]*tmp1[18]+z[35]*tmp1[19]+z[93]*tmp1[20]
                    +z[54]*tmp1[21]+z[15]*tmp1[22]+z[73]*tmp1[23]+z[34]*tmp1[24]+z[92]*tmp1[25]+z[53]*tmp1[26]+z[14]*tmp1[27]+z[72]*tmp1[28]+z[33]*tmp1[29]+z[91]*tmp1[30]
                    +z[52]*tmp1[31]+z[13]*tmp1[32]+z[71]*tmp1[33]+z[32]*tmp1[34]+z[90]*tmp1[35]+z[51]*tmp1[36]+z[12]*tmp1[37]+z[70]*tmp1[38]+z[31]*tmp1[39]+z[89]*tmp1[40]
                    +z[50]*tmp1[41]+z[11]*tmp1[42]+z[69]*tmp1[43]+z[30]*tmp1[44]+z[88]*tmp1[45]+z[49]*tmp1[46]+z[10]*tmp1[47]+z[68]*tmp1[48]+z[29]*tmp1[49]+z[87]*tmp1[50]
                    +z[48]*tmp1[51]+z[9]*tmp1[52]+z[67]*tmp1[53]+z[28]*tmp1[54]+z[86]*tmp1[55]+z[47]*tmp1[56]+z[8]*tmp1[57]+z[66]*tmp1[58]+z[27]*tmp1[59]+z[85]*tmp1[60]
                    +z[46]*tmp1[61]+z[7]*tmp1[62]+z[65]*tmp1[63]+z[26]*tmp1[64]+z[84]*tmp1[65]+z[45]*tmp1[66]+z[6]*tmp1[67]+z[64]*tmp1[68]+z[25]*tmp1[69]+z[83]*tmp1[70]
                    +z[44]*tmp1[71]+z[5]*tmp1[72]+z[63]*tmp1[73]+z[24]*tmp1[74]+z[82]*tmp1[75]+z[43]*tmp1[76]+z[4]*tmp1[77]+z[62]*tmp1[78]+z[23]*tmp1[79]+z[81]*tmp1[80]
                    +z[42]*tmp1[81]+z[3]*tmp1[82]+z[61]*tmp1[83]+z[22]*tmp1[84]+z[80]*tmp1[85]+z[41]*tmp1[86]+z[2]*tmp1[87]+z[60]*tmp1[88]+z[21]*tmp1[89]+z[79]*tmp1[90]
                    +z[40]*tmp1[91]+z[1]*tmp1[92]+z[59]*tmp1[93]+z[20]*tmp1[94]+z[78]*tmp1[95]+z[39]*tmp1[96];
                    tab[nb_tmp3+59*nb3]=tab[nb_tmp3+59*nb3]+z[0]*tmp1[0]
                    +z[59]*tmp1[1]+z[21]*tmp1[2]+z[80]*tmp1[3]+z[42]*tmp1[4]+z[4]*tmp1[5]+z[63]*tmp1[6]+z[25]*tmp1[7]+z[84]*tmp1[8]+z[46]*tmp1[9]+z[8]*tmp1[10]
                    +z[67]*tmp1[11]+z[29]*tmp1[12]+z[88]*tmp1[13]+z[50]*tmp1[14]+z[12]*tmp1[15]+z[71]*tmp1[16]+z[33]*tmp1[17]+z[92]*tmp1[18]+z[54]*tmp1[19]+z[16]*tmp1[20]
                    +z[75]*tmp1[21]+z[37]*tmp1[22]+z[96]*tmp1[23]+z[58]*tmp1[24]+z[20]*tmp1[25]+z[79]*tmp1[26]+z[41]*tmp1[27]+z[3]*tmp1[28]+z[62]*tmp1[29]+z[24]*tmp1[30]
                    +z[83]*tmp1[31]+z[45]*tmp1[32]+z[7]*tmp1[33]+z[66]*tmp1[34]+z[28]*tmp1[35]+z[87]*tmp1[36]+z[49]*tmp1[37]+z[11]*tmp1[38]+z[70]*tmp1[39]+z[32]*tmp1[40]
                    +z[91]*tmp1[41]+z[53]*tmp1[42]+z[15]*tmp1[43]+z[74]*tmp1[44]+z[36]*tmp1[45]+z[95]*tmp1[46]+z[57]*tmp1[47]+z[19]*tmp1[48]+z[78]*tmp1[49]+z[40]*tmp1[50]
                    +z[2]*tmp1[51]+z[61]*tmp1[52]+z[23]*tmp1[53]+z[82]*tmp1[54]+z[44]*tmp1[55]+z[6]*tmp1[56]+z[65]*tmp1[57]+z[27]*tmp1[58]+z[86]*tmp1[59]+z[48]*tmp1[60]
                    +z[10]*tmp1[61]+z[69]*tmp1[62]+z[31]*tmp1[63]+z[90]*tmp1[64]+z[52]*tmp1[65]+z[14]*tmp1[66]+z[73]*tmp1[67]+z[35]*tmp1[68]+z[94]*tmp1[69]+z[56]*tmp1[70]
                    +z[18]*tmp1[71]+z[77]*tmp1[72]+z[39]*tmp1[73]+z[1]*tmp1[74]+z[60]*tmp1[75]+z[22]*tmp1[76]+z[81]*tmp1[77]+z[43]*tmp1[78]+z[5]*tmp1[79]+z[64]*tmp1[80]
                    +z[26]*tmp1[81]+z[85]*tmp1[82]+z[47]*tmp1[83]+z[9]*tmp1[84]+z[68]*tmp1[85]+z[30]*tmp1[86]+z[89]*tmp1[87]+z[51]*tmp1[88]+z[13]*tmp1[89]+z[72]*tmp1[90]
                    +z[34]*tmp1[91]+z[93]*tmp1[92]+z[55]*tmp1[93]+z[17]*tmp1[94]+z[76]*tmp1[95]+z[38]*tmp1[96];
                    tab[nb_tmp3+60*nb3]=tab[nb_tmp3+60*nb3]+z[0]*tmp1[0]
                    +z[60]*tmp1[1]+z[23]*tmp1[2]+z[83]*tmp1[3]+z[46]*tmp1[4]+z[9]*tmp1[5]+z[69]*tmp1[6]+z[32]*tmp1[7]+z[92]*tmp1[8]+z[55]*tmp1[9]+z[18]*tmp1[10]
                    +z[78]*tmp1[11]+z[41]*tmp1[12]+z[4]*tmp1[13]+z[64]*tmp1[14]+z[27]*tmp1[15]+z[87]*tmp1[16]+z[50]*tmp1[17]+z[13]*tmp1[18]+z[73]*tmp1[19]+z[36]*tmp1[20]
                    +z[96]*tmp1[21]+z[59]*tmp1[22]+z[22]*tmp1[23]+z[82]*tmp1[24]+z[45]*tmp1[25]+z[8]*tmp1[26]+z[68]*tmp1[27]+z[31]*tmp1[28]+z[91]*tmp1[29]+z[54]*tmp1[30]
                    +z[17]*tmp1[31]+z[77]*tmp1[32]+z[40]*tmp1[33]+z[3]*tmp1[34]+z[63]*tmp1[35]+z[26]*tmp1[36]+z[86]*tmp1[37]+z[49]*tmp1[38]+z[12]*tmp1[39]+z[72]*tmp1[40]
                    +z[35]*tmp1[41]+z[95]*tmp1[42]+z[58]*tmp1[43]+z[21]*tmp1[44]+z[81]*tmp1[45]+z[44]*tmp1[46]+z[7]*tmp1[47]+z[67]*tmp1[48]+z[30]*tmp1[49]+z[90]*tmp1[50]
                    +z[53]*tmp1[51]+z[16]*tmp1[52]+z[76]*tmp1[53]+z[39]*tmp1[54]+z[2]*tmp1[55]+z[62]*tmp1[56]+z[25]*tmp1[57]+z[85]*tmp1[58]+z[48]*tmp1[59]+z[11]*tmp1[60]
                    +z[71]*tmp1[61]+z[34]*tmp1[62]+z[94]*tmp1[63]+z[57]*tmp1[64]+z[20]*tmp1[65]+z[80]*tmp1[66]+z[43]*tmp1[67]+z[6]*tmp1[68]+z[66]*tmp1[69]+z[29]*tmp1[70]
                    +z[89]*tmp1[71]+z[52]*tmp1[72]+z[15]*tmp1[73]+z[75]*tmp1[74]+z[38]*tmp1[75]+z[1]*tmp1[76]+z[61]*tmp1[77]+z[24]*tmp1[78]+z[84]*tmp1[79]+z[47]*tmp1[80]
                    +z[10]*tmp1[81]+z[70]*tmp1[82]+z[33]*tmp1[83]+z[93]*tmp1[84]+z[56]*tmp1[85]+z[19]*tmp1[86]+z[79]*tmp1[87]+z[42]*tmp1[88]+z[5]*tmp1[89]+z[65]*tmp1[90]
                    +z[28]*tmp1[91]+z[88]*tmp1[92]+z[51]*tmp1[93]+z[14]*tmp1[94]+z[74]*tmp1[95]+z[37]*tmp1[96];
                    tab[nb_tmp3+61*nb3]=tab[nb_tmp3+61*nb3]+z[0]*tmp1[0]
                    +z[61]*tmp1[1]+z[25]*tmp1[2]+z[86]*tmp1[3]+z[50]*tmp1[4]+z[14]*tmp1[5]+z[75]*tmp1[6]+z[39]*tmp1[7]+z[3]*tmp1[8]+z[64]*tmp1[9]+z[28]*tmp1[10]
                    +z[89]*tmp1[11]+z[53]*tmp1[12]+z[17]*tmp1[13]+z[78]*tmp1[14]+z[42]*tmp1[15]+z[6]*tmp1[16]+z[67]*tmp1[17]+z[31]*tmp1[18]+z[92]*tmp1[19]+z[56]*tmp1[20]
                    +z[20]*tmp1[21]+z[81]*tmp1[22]+z[45]*tmp1[23]+z[9]*tmp1[24]+z[70]*tmp1[25]+z[34]*tmp1[26]+z[95]*tmp1[27]+z[59]*tmp1[28]+z[23]*tmp1[29]+z[84]*tmp1[30]
                    +z[48]*tmp1[31]+z[12]*tmp1[32]+z[73]*tmp1[33]+z[37]*tmp1[34]+z[1]*tmp1[35]+z[62]*tmp1[36]+z[26]*tmp1[37]+z[87]*tmp1[38]+z[51]*tmp1[39]+z[15]*tmp1[40]
                    +z[76]*tmp1[41]+z[40]*tmp1[42]+z[4]*tmp1[43]+z[65]*tmp1[44]+z[29]*tmp1[45]+z[90]*tmp1[46]+z[54]*tmp1[47]+z[18]*tmp1[48]+z[79]*tmp1[49]+z[43]*tmp1[50]
                    +z[7]*tmp1[51]+z[68]*tmp1[52]+z[32]*tmp1[53]+z[93]*tmp1[54]+z[57]*tmp1[55]+z[21]*tmp1[56]+z[82]*tmp1[57]+z[46]*tmp1[58]+z[10]*tmp1[59]+z[71]*tmp1[60]
                    +z[35]*tmp1[61]+z[96]*tmp1[62]+z[60]*tmp1[63]+z[24]*tmp1[64]+z[85]*tmp1[65]+z[49]*tmp1[66]+z[13]*tmp1[67]+z[74]*tmp1[68]+z[38]*tmp1[69]+z[2]*tmp1[70]
                    +z[63]*tmp1[71]+z[27]*tmp1[72]+z[88]*tmp1[73]+z[52]*tmp1[74]+z[16]*tmp1[75]+z[77]*tmp1[76]+z[41]*tmp1[77]+z[5]*tmp1[78]+z[66]*tmp1[79]+z[30]*tmp1[80]
                    +z[91]*tmp1[81]+z[55]*tmp1[82]+z[19]*tmp1[83]+z[80]*tmp1[84]+z[44]*tmp1[85]+z[8]*tmp1[86]+z[69]*tmp1[87]+z[33]*tmp1[88]+z[94]*tmp1[89]+z[58]*tmp1[90]
                    +z[22]*tmp1[91]+z[83]*tmp1[92]+z[47]*tmp1[93]+z[11]*tmp1[94]+z[72]*tmp1[95]+z[36]*tmp1[96];
                    tab[nb_tmp3+62*nb3]=tab[nb_tmp3+62*nb3]+z[0]*tmp1[0]
                    +z[62]*tmp1[1]+z[27]*tmp1[2]+z[89]*tmp1[3]+z[54]*tmp1[4]+z[19]*tmp1[5]+z[81]*tmp1[6]+z[46]*tmp1[7]+z[11]*tmp1[8]+z[73]*tmp1[9]+z[38]*tmp1[10]
                    +z[3]*tmp1[11]+z[65]*tmp1[12]+z[30]*tmp1[13]+z[92]*tmp1[14]+z[57]*tmp1[15]+z[22]*tmp1[16]+z[84]*tmp1[17]+z[49]*tmp1[18]+z[14]*tmp1[19]+z[76]*tmp1[20]
                    +z[41]*tmp1[21]+z[6]*tmp1[22]+z[68]*tmp1[23]+z[33]*tmp1[24]+z[95]*tmp1[25]+z[60]*tmp1[26]+z[25]*tmp1[27]+z[87]*tmp1[28]+z[52]*tmp1[29]+z[17]*tmp1[30]
                    +z[79]*tmp1[31]+z[44]*tmp1[32]+z[9]*tmp1[33]+z[71]*tmp1[34]+z[36]*tmp1[35]+z[1]*tmp1[36]+z[63]*tmp1[37]+z[28]*tmp1[38]+z[90]*tmp1[39]+z[55]*tmp1[40]
                    +z[20]*tmp1[41]+z[82]*tmp1[42]+z[47]*tmp1[43]+z[12]*tmp1[44]+z[74]*tmp1[45]+z[39]*tmp1[46]+z[4]*tmp1[47]+z[66]*tmp1[48]+z[31]*tmp1[49]+z[93]*tmp1[50]
                    +z[58]*tmp1[51]+z[23]*tmp1[52]+z[85]*tmp1[53]+z[50]*tmp1[54]+z[15]*tmp1[55]+z[77]*tmp1[56]+z[42]*tmp1[57]+z[7]*tmp1[58]+z[69]*tmp1[59]+z[34]*tmp1[60]
                    +z[96]*tmp1[61]+z[61]*tmp1[62]+z[26]*tmp1[63]+z[88]*tmp1[64]+z[53]*tmp1[65]+z[18]*tmp1[66]+z[80]*tmp1[67]+z[45]*tmp1[68]+z[10]*tmp1[69]+z[72]*tmp1[70]
                    +z[37]*tmp1[71]+z[2]*tmp1[72]+z[64]*tmp1[73]+z[29]*tmp1[74]+z[91]*tmp1[75]+z[56]*tmp1[76]+z[21]*tmp1[77]+z[83]*tmp1[78]+z[48]*tmp1[79]+z[13]*tmp1[80]
                    +z[75]*tmp1[81]+z[40]*tmp1[82]+z[5]*tmp1[83]+z[67]*tmp1[84]+z[32]*tmp1[85]+z[94]*tmp1[86]+z[59]*tmp1[87]+z[24]*tmp1[88]+z[86]*tmp1[89]+z[51]*tmp1[90]
                    +z[16]*tmp1[91]+z[78]*tmp1[92]+z[43]*tmp1[93]+z[8]*tmp1[94]+z[70]*tmp1[95]+z[35]*tmp1[96];
                    tab[nb_tmp3+63*nb3]=tab[nb_tmp3+63*nb3]+z[0]*tmp1[0]
                    +z[63]*tmp1[1]+z[29]*tmp1[2]+z[92]*tmp1[3]+z[58]*tmp1[4]+z[24]*tmp1[5]+z[87]*tmp1[6]+z[53]*tmp1[7]+z[19]*tmp1[8]+z[82]*tmp1[9]+z[48]*tmp1[10]
                    +z[14]*tmp1[11]+z[77]*tmp1[12]+z[43]*tmp1[13]+z[9]*tmp1[14]+z[72]*tmp1[15]+z[38]*tmp1[16]+z[4]*tmp1[17]+z[67]*tmp1[18]+z[33]*tmp1[19]+z[96]*tmp1[20]
                    +z[62]*tmp1[21]+z[28]*tmp1[22]+z[91]*tmp1[23]+z[57]*tmp1[24]+z[23]*tmp1[25]+z[86]*tmp1[26]+z[52]*tmp1[27]+z[18]*tmp1[28]+z[81]*tmp1[29]+z[47]*tmp1[30]
                    +z[13]*tmp1[31]+z[76]*tmp1[32]+z[42]*tmp1[33]+z[8]*tmp1[34]+z[71]*tmp1[35]+z[37]*tmp1[36]+z[3]*tmp1[37]+z[66]*tmp1[38]+z[32]*tmp1[39]+z[95]*tmp1[40]
                    +z[61]*tmp1[41]+z[27]*tmp1[42]+z[90]*tmp1[43]+z[56]*tmp1[44]+z[22]*tmp1[45]+z[85]*tmp1[46]+z[51]*tmp1[47]+z[17]*tmp1[48]+z[80]*tmp1[49]+z[46]*tmp1[50]
                    +z[12]*tmp1[51]+z[75]*tmp1[52]+z[41]*tmp1[53]+z[7]*tmp1[54]+z[70]*tmp1[55]+z[36]*tmp1[56]+z[2]*tmp1[57]+z[65]*tmp1[58]+z[31]*tmp1[59]+z[94]*tmp1[60]
                    +z[60]*tmp1[61]+z[26]*tmp1[62]+z[89]*tmp1[63]+z[55]*tmp1[64]+z[21]*tmp1[65]+z[84]*tmp1[66]+z[50]*tmp1[67]+z[16]*tmp1[68]+z[79]*tmp1[69]+z[45]*tmp1[70]
                    +z[11]*tmp1[71]+z[74]*tmp1[72]+z[40]*tmp1[73]+z[6]*tmp1[74]+z[69]*tmp1[75]+z[35]*tmp1[76]+z[1]*tmp1[77]+z[64]*tmp1[78]+z[30]*tmp1[79]+z[93]*tmp1[80]
                    +z[59]*tmp1[81]+z[25]*tmp1[82]+z[88]*tmp1[83]+z[54]*tmp1[84]+z[20]*tmp1[85]+z[83]*tmp1[86]+z[49]*tmp1[87]+z[15]*tmp1[88]+z[78]*tmp1[89]+z[44]*tmp1[90]
                    +z[10]*tmp1[91]+z[73]*tmp1[92]+z[39]*tmp1[93]+z[5]*tmp1[94]+z[68]*tmp1[95]+z[34]*tmp1[96];
                    tab[nb_tmp3+64*nb3]=tab[nb_tmp3+64*nb3]+z[0]*tmp1[0]
                    +z[64]*tmp1[1]+z[31]*tmp1[2]+z[95]*tmp1[3]+z[62]*tmp1[4]+z[29]*tmp1[5]+z[93]*tmp1[6]+z[60]*tmp1[7]+z[27]*tmp1[8]+z[91]*tmp1[9]+z[58]*tmp1[10]
                    +z[25]*tmp1[11]+z[89]*tmp1[12]+z[56]*tmp1[13]+z[23]*tmp1[14]+z[87]*tmp1[15]+z[54]*tmp1[16]+z[21]*tmp1[17]+z[85]*tmp1[18]+z[52]*tmp1[19]+z[19]*tmp1[20]
                    +z[83]*tmp1[21]+z[50]*tmp1[22]+z[17]*tmp1[23]+z[81]*tmp1[24]+z[48]*tmp1[25]+z[15]*tmp1[26]+z[79]*tmp1[27]+z[46]*tmp1[28]+z[13]*tmp1[29]+z[77]*tmp1[30]
                    +z[44]*tmp1[31]+z[11]*tmp1[32]+z[75]*tmp1[33]+z[42]*tmp1[34]+z[9]*tmp1[35]+z[73]*tmp1[36]+z[40]*tmp1[37]+z[7]*tmp1[38]+z[71]*tmp1[39]+z[38]*tmp1[40]
                    +z[5]*tmp1[41]+z[69]*tmp1[42]+z[36]*tmp1[43]+z[3]*tmp1[44]+z[67]*tmp1[45]+z[34]*tmp1[46]+z[1]*tmp1[47]+z[65]*tmp1[48]+z[32]*tmp1[49]+z[96]*tmp1[50]
                    +z[63]*tmp1[51]+z[30]*tmp1[52]+z[94]*tmp1[53]+z[61]*tmp1[54]+z[28]*tmp1[55]+z[92]*tmp1[56]+z[59]*tmp1[57]+z[26]*tmp1[58]+z[90]*tmp1[59]+z[57]*tmp1[60]
                    +z[24]*tmp1[61]+z[88]*tmp1[62]+z[55]*tmp1[63]+z[22]*tmp1[64]+z[86]*tmp1[65]+z[53]*tmp1[66]+z[20]*tmp1[67]+z[84]*tmp1[68]+z[51]*tmp1[69]+z[18]*tmp1[70]
                    +z[82]*tmp1[71]+z[49]*tmp1[72]+z[16]*tmp1[73]+z[80]*tmp1[74]+z[47]*tmp1[75]+z[14]*tmp1[76]+z[78]*tmp1[77]+z[45]*tmp1[78]+z[12]*tmp1[79]+z[76]*tmp1[80]
                    +z[43]*tmp1[81]+z[10]*tmp1[82]+z[74]*tmp1[83]+z[41]*tmp1[84]+z[8]*tmp1[85]+z[72]*tmp1[86]+z[39]*tmp1[87]+z[6]*tmp1[88]+z[70]*tmp1[89]+z[37]*tmp1[90]
                    +z[4]*tmp1[91]+z[68]*tmp1[92]+z[35]*tmp1[93]+z[2]*tmp1[94]+z[66]*tmp1[95]+z[33]*tmp1[96];
                    tab[nb_tmp3+65*nb3]=tab[nb_tmp3+65*nb3]+z[0]*tmp1[0]
                    +z[65]*tmp1[1]+z[33]*tmp1[2]+z[1]*tmp1[3]+z[66]*tmp1[4]+z[34]*tmp1[5]+z[2]*tmp1[6]+z[67]*tmp1[7]+z[35]*tmp1[8]+z[3]*tmp1[9]+z[68]*tmp1[10]
                    +z[36]*tmp1[11]+z[4]*tmp1[12]+z[69]*tmp1[13]+z[37]*tmp1[14]+z[5]*tmp1[15]+z[70]*tmp1[16]+z[38]*tmp1[17]+z[6]*tmp1[18]+z[71]*tmp1[19]+z[39]*tmp1[20]
                    +z[7]*tmp1[21]+z[72]*tmp1[22]+z[40]*tmp1[23]+z[8]*tmp1[24]+z[73]*tmp1[25]+z[41]*tmp1[26]+z[9]*tmp1[27]+z[74]*tmp1[28]+z[42]*tmp1[29]+z[10]*tmp1[30]
                    +z[75]*tmp1[31]+z[43]*tmp1[32]+z[11]*tmp1[33]+z[76]*tmp1[34]+z[44]*tmp1[35]+z[12]*tmp1[36]+z[77]*tmp1[37]+z[45]*tmp1[38]+z[13]*tmp1[39]+z[78]*tmp1[40]
                    +z[46]*tmp1[41]+z[14]*tmp1[42]+z[79]*tmp1[43]+z[47]*tmp1[44]+z[15]*tmp1[45]+z[80]*tmp1[46]+z[48]*tmp1[47]+z[16]*tmp1[48]+z[81]*tmp1[49]+z[49]*tmp1[50]
                    +z[17]*tmp1[51]+z[82]*tmp1[52]+z[50]*tmp1[53]+z[18]*tmp1[54]+z[83]*tmp1[55]+z[51]*tmp1[56]+z[19]*tmp1[57]+z[84]*tmp1[58]+z[52]*tmp1[59]+z[20]*tmp1[60]
                    +z[85]*tmp1[61]+z[53]*tmp1[62]+z[21]*tmp1[63]+z[86]*tmp1[64]+z[54]*tmp1[65]+z[22]*tmp1[66]+z[87]*tmp1[67]+z[55]*tmp1[68]+z[23]*tmp1[69]+z[88]*tmp1[70]
                    +z[56]*tmp1[71]+z[24]*tmp1[72]+z[89]*tmp1[73]+z[57]*tmp1[74]+z[25]*tmp1[75]+z[90]*tmp1[76]+z[58]*tmp1[77]+z[26]*tmp1[78]+z[91]*tmp1[79]+z[59]*tmp1[80]
                    +z[27]*tmp1[81]+z[92]*tmp1[82]+z[60]*tmp1[83]+z[28]*tmp1[84]+z[93]*tmp1[85]+z[61]*tmp1[86]+z[29]*tmp1[87]+z[94]*tmp1[88]+z[62]*tmp1[89]+z[30]*tmp1[90]
                    +z[95]*tmp1[91]+z[63]*tmp1[92]+z[31]*tmp1[93]+z[96]*tmp1[94]+z[64]*tmp1[95]+z[32]*tmp1[96];
                    tab[nb_tmp3+66*nb3]=tab[nb_tmp3+66*nb3]+z[0]*tmp1[0]
                    +z[66]*tmp1[1]+z[35]*tmp1[2]+z[4]*tmp1[3]+z[70]*tmp1[4]+z[39]*tmp1[5]+z[8]*tmp1[6]+z[74]*tmp1[7]+z[43]*tmp1[8]+z[12]*tmp1[9]+z[78]*tmp1[10]
                    +z[47]*tmp1[11]+z[16]*tmp1[12]+z[82]*tmp1[13]+z[51]*tmp1[14]+z[20]*tmp1[15]+z[86]*tmp1[16]+z[55]*tmp1[17]+z[24]*tmp1[18]+z[90]*tmp1[19]+z[59]*tmp1[20]
                    +z[28]*tmp1[21]+z[94]*tmp1[22]+z[63]*tmp1[23]+z[32]*tmp1[24]+z[1]*tmp1[25]+z[67]*tmp1[26]+z[36]*tmp1[27]+z[5]*tmp1[28]+z[71]*tmp1[29]+z[40]*tmp1[30]
                    +z[9]*tmp1[31]+z[75]*tmp1[32]+z[44]*tmp1[33]+z[13]*tmp1[34]+z[79]*tmp1[35]+z[48]*tmp1[36]+z[17]*tmp1[37]+z[83]*tmp1[38]+z[52]*tmp1[39]+z[21]*tmp1[40]
                    +z[87]*tmp1[41]+z[56]*tmp1[42]+z[25]*tmp1[43]+z[91]*tmp1[44]+z[60]*tmp1[45]+z[29]*tmp1[46]+z[95]*tmp1[47]+z[64]*tmp1[48]+z[33]*tmp1[49]+z[2]*tmp1[50]
                    +z[68]*tmp1[51]+z[37]*tmp1[52]+z[6]*tmp1[53]+z[72]*tmp1[54]+z[41]*tmp1[55]+z[10]*tmp1[56]+z[76]*tmp1[57]+z[45]*tmp1[58]+z[14]*tmp1[59]+z[80]*tmp1[60]
                    +z[49]*tmp1[61]+z[18]*tmp1[62]+z[84]*tmp1[63]+z[53]*tmp1[64]+z[22]*tmp1[65]+z[88]*tmp1[66]+z[57]*tmp1[67]+z[26]*tmp1[68]+z[92]*tmp1[69]+z[61]*tmp1[70]
                    +z[30]*tmp1[71]+z[96]*tmp1[72]+z[65]*tmp1[73]+z[34]*tmp1[74]+z[3]*tmp1[75]+z[69]*tmp1[76]+z[38]*tmp1[77]+z[7]*tmp1[78]+z[73]*tmp1[79]+z[42]*tmp1[80]
                    +z[11]*tmp1[81]+z[77]*tmp1[82]+z[46]*tmp1[83]+z[15]*tmp1[84]+z[81]*tmp1[85]+z[50]*tmp1[86]+z[19]*tmp1[87]+z[85]*tmp1[88]+z[54]*tmp1[89]+z[23]*tmp1[90]
                    +z[89]*tmp1[91]+z[58]*tmp1[92]+z[27]*tmp1[93]+z[93]*tmp1[94]+z[62]*tmp1[95]+z[31]*tmp1[96];
                    tab[nb_tmp3+67*nb3]=tab[nb_tmp3+67*nb3]+z[0]*tmp1[0]
                    +z[67]*tmp1[1]+z[37]*tmp1[2]+z[7]*tmp1[3]+z[74]*tmp1[4]+z[44]*tmp1[5]+z[14]*tmp1[6]+z[81]*tmp1[7]+z[51]*tmp1[8]+z[21]*tmp1[9]+z[88]*tmp1[10]
                    +z[58]*tmp1[11]+z[28]*tmp1[12]+z[95]*tmp1[13]+z[65]*tmp1[14]+z[35]*tmp1[15]+z[5]*tmp1[16]+z[72]*tmp1[17]+z[42]*tmp1[18]+z[12]*tmp1[19]+z[79]*tmp1[20]
                    +z[49]*tmp1[21]+z[19]*tmp1[22]+z[86]*tmp1[23]+z[56]*tmp1[24]+z[26]*tmp1[25]+z[93]*tmp1[26]+z[63]*tmp1[27]+z[33]*tmp1[28]+z[3]*tmp1[29]+z[70]*tmp1[30]
                    +z[40]*tmp1[31]+z[10]*tmp1[32]+z[77]*tmp1[33]+z[47]*tmp1[34]+z[17]*tmp1[35]+z[84]*tmp1[36]+z[54]*tmp1[37]+z[24]*tmp1[38]+z[91]*tmp1[39]+z[61]*tmp1[40]
                    +z[31]*tmp1[41]+z[1]*tmp1[42]+z[68]*tmp1[43]+z[38]*tmp1[44]+z[8]*tmp1[45]+z[75]*tmp1[46]+z[45]*tmp1[47]+z[15]*tmp1[48]+z[82]*tmp1[49]+z[52]*tmp1[50]
                    +z[22]*tmp1[51]+z[89]*tmp1[52]+z[59]*tmp1[53]+z[29]*tmp1[54]+z[96]*tmp1[55]+z[66]*tmp1[56]+z[36]*tmp1[57]+z[6]*tmp1[58]+z[73]*tmp1[59]+z[43]*tmp1[60]
                    +z[13]*tmp1[61]+z[80]*tmp1[62]+z[50]*tmp1[63]+z[20]*tmp1[64]+z[87]*tmp1[65]+z[57]*tmp1[66]+z[27]*tmp1[67]+z[94]*tmp1[68]+z[64]*tmp1[69]+z[34]*tmp1[70]
                    +z[4]*tmp1[71]+z[71]*tmp1[72]+z[41]*tmp1[73]+z[11]*tmp1[74]+z[78]*tmp1[75]+z[48]*tmp1[76]+z[18]*tmp1[77]+z[85]*tmp1[78]+z[55]*tmp1[79]+z[25]*tmp1[80]
                    +z[92]*tmp1[81]+z[62]*tmp1[82]+z[32]*tmp1[83]+z[2]*tmp1[84]+z[69]*tmp1[85]+z[39]*tmp1[86]+z[9]*tmp1[87]+z[76]*tmp1[88]+z[46]*tmp1[89]+z[16]*tmp1[90]
                    +z[83]*tmp1[91]+z[53]*tmp1[92]+z[23]*tmp1[93]+z[90]*tmp1[94]+z[60]*tmp1[95]+z[30]*tmp1[96];
                    tab[nb_tmp3+68*nb3]=tab[nb_tmp3+68*nb3]+z[0]*tmp1[0]
                    +z[68]*tmp1[1]+z[39]*tmp1[2]+z[10]*tmp1[3]+z[78]*tmp1[4]+z[49]*tmp1[5]+z[20]*tmp1[6]+z[88]*tmp1[7]+z[59]*tmp1[8]+z[30]*tmp1[9]+z[1]*tmp1[10]
                    +z[69]*tmp1[11]+z[40]*tmp1[12]+z[11]*tmp1[13]+z[79]*tmp1[14]+z[50]*tmp1[15]+z[21]*tmp1[16]+z[89]*tmp1[17]+z[60]*tmp1[18]+z[31]*tmp1[19]+z[2]*tmp1[20]
                    +z[70]*tmp1[21]+z[41]*tmp1[22]+z[12]*tmp1[23]+z[80]*tmp1[24]+z[51]*tmp1[25]+z[22]*tmp1[26]+z[90]*tmp1[27]+z[61]*tmp1[28]+z[32]*tmp1[29]+z[3]*tmp1[30]
                    +z[71]*tmp1[31]+z[42]*tmp1[32]+z[13]*tmp1[33]+z[81]*tmp1[34]+z[52]*tmp1[35]+z[23]*tmp1[36]+z[91]*tmp1[37]+z[62]*tmp1[38]+z[33]*tmp1[39]+z[4]*tmp1[40]
                    +z[72]*tmp1[41]+z[43]*tmp1[42]+z[14]*tmp1[43]+z[82]*tmp1[44]+z[53]*tmp1[45]+z[24]*tmp1[46]+z[92]*tmp1[47]+z[63]*tmp1[48]+z[34]*tmp1[49]+z[5]*tmp1[50]
                    +z[73]*tmp1[51]+z[44]*tmp1[52]+z[15]*tmp1[53]+z[83]*tmp1[54]+z[54]*tmp1[55]+z[25]*tmp1[56]+z[93]*tmp1[57]+z[64]*tmp1[58]+z[35]*tmp1[59]+z[6]*tmp1[60]
                    +z[74]*tmp1[61]+z[45]*tmp1[62]+z[16]*tmp1[63]+z[84]*tmp1[64]+z[55]*tmp1[65]+z[26]*tmp1[66]+z[94]*tmp1[67]+z[65]*tmp1[68]+z[36]*tmp1[69]+z[7]*tmp1[70]
                    +z[75]*tmp1[71]+z[46]*tmp1[72]+z[17]*tmp1[73]+z[85]*tmp1[74]+z[56]*tmp1[75]+z[27]*tmp1[76]+z[95]*tmp1[77]+z[66]*tmp1[78]+z[37]*tmp1[79]+z[8]*tmp1[80]
                    +z[76]*tmp1[81]+z[47]*tmp1[82]+z[18]*tmp1[83]+z[86]*tmp1[84]+z[57]*tmp1[85]+z[28]*tmp1[86]+z[96]*tmp1[87]+z[67]*tmp1[88]+z[38]*tmp1[89]+z[9]*tmp1[90]
                    +z[77]*tmp1[91]+z[48]*tmp1[92]+z[19]*tmp1[93]+z[87]*tmp1[94]+z[58]*tmp1[95]+z[29]*tmp1[96];
                    tab[nb_tmp3+69*nb3]=tab[nb_tmp3+69*nb3]+z[0]*tmp1[0]
                    +z[69]*tmp1[1]+z[41]*tmp1[2]+z[13]*tmp1[3]+z[82]*tmp1[4]+z[54]*tmp1[5]+z[26]*tmp1[6]+z[95]*tmp1[7]+z[67]*tmp1[8]+z[39]*tmp1[9]+z[11]*tmp1[10]
                    +z[80]*tmp1[11]+z[52]*tmp1[12]+z[24]*tmp1[13]+z[93]*tmp1[14]+z[65]*tmp1[15]+z[37]*tmp1[16]+z[9]*tmp1[17]+z[78]*tmp1[18]+z[50]*tmp1[19]+z[22]*tmp1[20]
                    +z[91]*tmp1[21]+z[63]*tmp1[22]+z[35]*tmp1[23]+z[7]*tmp1[24]+z[76]*tmp1[25]+z[48]*tmp1[26]+z[20]*tmp1[27]+z[89]*tmp1[28]+z[61]*tmp1[29]+z[33]*tmp1[30]
                    +z[5]*tmp1[31]+z[74]*tmp1[32]+z[46]*tmp1[33]+z[18]*tmp1[34]+z[87]*tmp1[35]+z[59]*tmp1[36]+z[31]*tmp1[37]+z[3]*tmp1[38]+z[72]*tmp1[39]+z[44]*tmp1[40]
                    +z[16]*tmp1[41]+z[85]*tmp1[42]+z[57]*tmp1[43]+z[29]*tmp1[44]+z[1]*tmp1[45]+z[70]*tmp1[46]+z[42]*tmp1[47]+z[14]*tmp1[48]+z[83]*tmp1[49]+z[55]*tmp1[50]
                    +z[27]*tmp1[51]+z[96]*tmp1[52]+z[68]*tmp1[53]+z[40]*tmp1[54]+z[12]*tmp1[55]+z[81]*tmp1[56]+z[53]*tmp1[57]+z[25]*tmp1[58]+z[94]*tmp1[59]+z[66]*tmp1[60]
                    +z[38]*tmp1[61]+z[10]*tmp1[62]+z[79]*tmp1[63]+z[51]*tmp1[64]+z[23]*tmp1[65]+z[92]*tmp1[66]+z[64]*tmp1[67]+z[36]*tmp1[68]+z[8]*tmp1[69]+z[77]*tmp1[70]
                    +z[49]*tmp1[71]+z[21]*tmp1[72]+z[90]*tmp1[73]+z[62]*tmp1[74]+z[34]*tmp1[75]+z[6]*tmp1[76]+z[75]*tmp1[77]+z[47]*tmp1[78]+z[19]*tmp1[79]+z[88]*tmp1[80]
                    +z[60]*tmp1[81]+z[32]*tmp1[82]+z[4]*tmp1[83]+z[73]*tmp1[84]+z[45]*tmp1[85]+z[17]*tmp1[86]+z[86]*tmp1[87]+z[58]*tmp1[88]+z[30]*tmp1[89]+z[2]*tmp1[90]
                    +z[71]*tmp1[91]+z[43]*tmp1[92]+z[15]*tmp1[93]+z[84]*tmp1[94]+z[56]*tmp1[95]+z[28]*tmp1[96];
                    tab[nb_tmp3+70*nb3]=tab[nb_tmp3+70*nb3]+z[0]*tmp1[0]
                    +z[70]*tmp1[1]+z[43]*tmp1[2]+z[16]*tmp1[3]+z[86]*tmp1[4]+z[59]*tmp1[5]+z[32]*tmp1[6]+z[5]*tmp1[7]+z[75]*tmp1[8]+z[48]*tmp1[9]+z[21]*tmp1[10]
                    +z[91]*tmp1[11]+z[64]*tmp1[12]+z[37]*tmp1[13]+z[10]*tmp1[14]+z[80]*tmp1[15]+z[53]*tmp1[16]+z[26]*tmp1[17]+z[96]*tmp1[18]+z[69]*tmp1[19]+z[42]*tmp1[20]
                    +z[15]*tmp1[21]+z[85]*tmp1[22]+z[58]*tmp1[23]+z[31]*tmp1[24]+z[4]*tmp1[25]+z[74]*tmp1[26]+z[47]*tmp1[27]+z[20]*tmp1[28]+z[90]*tmp1[29]+z[63]*tmp1[30]
                    +z[36]*tmp1[31]+z[9]*tmp1[32]+z[79]*tmp1[33]+z[52]*tmp1[34]+z[25]*tmp1[35]+z[95]*tmp1[36]+z[68]*tmp1[37]+z[41]*tmp1[38]+z[14]*tmp1[39]+z[84]*tmp1[40]
                    +z[57]*tmp1[41]+z[30]*tmp1[42]+z[3]*tmp1[43]+z[73]*tmp1[44]+z[46]*tmp1[45]+z[19]*tmp1[46]+z[89]*tmp1[47]+z[62]*tmp1[48]+z[35]*tmp1[49]+z[8]*tmp1[50]
                    +z[78]*tmp1[51]+z[51]*tmp1[52]+z[24]*tmp1[53]+z[94]*tmp1[54]+z[67]*tmp1[55]+z[40]*tmp1[56]+z[13]*tmp1[57]+z[83]*tmp1[58]+z[56]*tmp1[59]+z[29]*tmp1[60]
                    +z[2]*tmp1[61]+z[72]*tmp1[62]+z[45]*tmp1[63]+z[18]*tmp1[64]+z[88]*tmp1[65]+z[61]*tmp1[66]+z[34]*tmp1[67]+z[7]*tmp1[68]+z[77]*tmp1[69]+z[50]*tmp1[70]
                    +z[23]*tmp1[71]+z[93]*tmp1[72]+z[66]*tmp1[73]+z[39]*tmp1[74]+z[12]*tmp1[75]+z[82]*tmp1[76]+z[55]*tmp1[77]+z[28]*tmp1[78]+z[1]*tmp1[79]+z[71]*tmp1[80]
                    +z[44]*tmp1[81]+z[17]*tmp1[82]+z[87]*tmp1[83]+z[60]*tmp1[84]+z[33]*tmp1[85]+z[6]*tmp1[86]+z[76]*tmp1[87]+z[49]*tmp1[88]+z[22]*tmp1[89]+z[92]*tmp1[90]
                    +z[65]*tmp1[91]+z[38]*tmp1[92]+z[11]*tmp1[93]+z[81]*tmp1[94]+z[54]*tmp1[95]+z[27]*tmp1[96];
                    tab[nb_tmp3+71*nb3]=tab[nb_tmp3+71*nb3]+z[0]*tmp1[0]
                    +z[71]*tmp1[1]+z[45]*tmp1[2]+z[19]*tmp1[3]+z[90]*tmp1[4]+z[64]*tmp1[5]+z[38]*tmp1[6]+z[12]*tmp1[7]+z[83]*tmp1[8]+z[57]*tmp1[9]+z[31]*tmp1[10]
                    +z[5]*tmp1[11]+z[76]*tmp1[12]+z[50]*tmp1[13]+z[24]*tmp1[14]+z[95]*tmp1[15]+z[69]*tmp1[16]+z[43]*tmp1[17]+z[17]*tmp1[18]+z[88]*tmp1[19]+z[62]*tmp1[20]
                    +z[36]*tmp1[21]+z[10]*tmp1[22]+z[81]*tmp1[23]+z[55]*tmp1[24]+z[29]*tmp1[25]+z[3]*tmp1[26]+z[74]*tmp1[27]+z[48]*tmp1[28]+z[22]*tmp1[29]+z[93]*tmp1[30]
                    +z[67]*tmp1[31]+z[41]*tmp1[32]+z[15]*tmp1[33]+z[86]*tmp1[34]+z[60]*tmp1[35]+z[34]*tmp1[36]+z[8]*tmp1[37]+z[79]*tmp1[38]+z[53]*tmp1[39]+z[27]*tmp1[40]
                    +z[1]*tmp1[41]+z[72]*tmp1[42]+z[46]*tmp1[43]+z[20]*tmp1[44]+z[91]*tmp1[45]+z[65]*tmp1[46]+z[39]*tmp1[47]+z[13]*tmp1[48]+z[84]*tmp1[49]+z[58]*tmp1[50]
                    +z[32]*tmp1[51]+z[6]*tmp1[52]+z[77]*tmp1[53]+z[51]*tmp1[54]+z[25]*tmp1[55]+z[96]*tmp1[56]+z[70]*tmp1[57]+z[44]*tmp1[58]+z[18]*tmp1[59]+z[89]*tmp1[60]
                    +z[63]*tmp1[61]+z[37]*tmp1[62]+z[11]*tmp1[63]+z[82]*tmp1[64]+z[56]*tmp1[65]+z[30]*tmp1[66]+z[4]*tmp1[67]+z[75]*tmp1[68]+z[49]*tmp1[69]+z[23]*tmp1[70]
                    +z[94]*tmp1[71]+z[68]*tmp1[72]+z[42]*tmp1[73]+z[16]*tmp1[74]+z[87]*tmp1[75]+z[61]*tmp1[76]+z[35]*tmp1[77]+z[9]*tmp1[78]+z[80]*tmp1[79]+z[54]*tmp1[80]
                    +z[28]*tmp1[81]+z[2]*tmp1[82]+z[73]*tmp1[83]+z[47]*tmp1[84]+z[21]*tmp1[85]+z[92]*tmp1[86]+z[66]*tmp1[87]+z[40]*tmp1[88]+z[14]*tmp1[89]+z[85]*tmp1[90]
                    +z[59]*tmp1[91]+z[33]*tmp1[92]+z[7]*tmp1[93]+z[78]*tmp1[94]+z[52]*tmp1[95]+z[26]*tmp1[96];
                    tab[nb_tmp3+72*nb3]=tab[nb_tmp3+72*nb3]+z[0]*tmp1[0]
                    +z[72]*tmp1[1]+z[47]*tmp1[2]+z[22]*tmp1[3]+z[94]*tmp1[4]+z[69]*tmp1[5]+z[44]*tmp1[6]+z[19]*tmp1[7]+z[91]*tmp1[8]+z[66]*tmp1[9]+z[41]*tmp1[10]
                    +z[16]*tmp1[11]+z[88]*tmp1[12]+z[63]*tmp1[13]+z[38]*tmp1[14]+z[13]*tmp1[15]+z[85]*tmp1[16]+z[60]*tmp1[17]+z[35]*tmp1[18]+z[10]*tmp1[19]+z[82]*tmp1[20]
                    +z[57]*tmp1[21]+z[32]*tmp1[22]+z[7]*tmp1[23]+z[79]*tmp1[24]+z[54]*tmp1[25]+z[29]*tmp1[26]+z[4]*tmp1[27]+z[76]*tmp1[28]+z[51]*tmp1[29]+z[26]*tmp1[30]
                    +z[1]*tmp1[31]+z[73]*tmp1[32]+z[48]*tmp1[33]+z[23]*tmp1[34]+z[95]*tmp1[35]+z[70]*tmp1[36]+z[45]*tmp1[37]+z[20]*tmp1[38]+z[92]*tmp1[39]+z[67]*tmp1[40]
                    +z[42]*tmp1[41]+z[17]*tmp1[42]+z[89]*tmp1[43]+z[64]*tmp1[44]+z[39]*tmp1[45]+z[14]*tmp1[46]+z[86]*tmp1[47]+z[61]*tmp1[48]+z[36]*tmp1[49]+z[11]*tmp1[50]
                    +z[83]*tmp1[51]+z[58]*tmp1[52]+z[33]*tmp1[53]+z[8]*tmp1[54]+z[80]*tmp1[55]+z[55]*tmp1[56]+z[30]*tmp1[57]+z[5]*tmp1[58]+z[77]*tmp1[59]+z[52]*tmp1[60]
                    +z[27]*tmp1[61]+z[2]*tmp1[62]+z[74]*tmp1[63]+z[49]*tmp1[64]+z[24]*tmp1[65]+z[96]*tmp1[66]+z[71]*tmp1[67]+z[46]*tmp1[68]+z[21]*tmp1[69]+z[93]*tmp1[70]
                    +z[68]*tmp1[71]+z[43]*tmp1[72]+z[18]*tmp1[73]+z[90]*tmp1[74]+z[65]*tmp1[75]+z[40]*tmp1[76]+z[15]*tmp1[77]+z[87]*tmp1[78]+z[62]*tmp1[79]+z[37]*tmp1[80]
                    +z[12]*tmp1[81]+z[84]*tmp1[82]+z[59]*tmp1[83]+z[34]*tmp1[84]+z[9]*tmp1[85]+z[81]*tmp1[86]+z[56]*tmp1[87]+z[31]*tmp1[88]+z[6]*tmp1[89]+z[78]*tmp1[90]
                    +z[53]*tmp1[91]+z[28]*tmp1[92]+z[3]*tmp1[93]+z[75]*tmp1[94]+z[50]*tmp1[95]+z[25]*tmp1[96];
                    tab[nb_tmp3+73*nb3]=tab[nb_tmp3+73*nb3]+z[0]*tmp1[0]
                    +z[73]*tmp1[1]+z[49]*tmp1[2]+z[25]*tmp1[3]+z[1]*tmp1[4]+z[74]*tmp1[5]+z[50]*tmp1[6]+z[26]*tmp1[7]+z[2]*tmp1[8]+z[75]*tmp1[9]+z[51]*tmp1[10]
                    +z[27]*tmp1[11]+z[3]*tmp1[12]+z[76]*tmp1[13]+z[52]*tmp1[14]+z[28]*tmp1[15]+z[4]*tmp1[16]+z[77]*tmp1[17]+z[53]*tmp1[18]+z[29]*tmp1[19]+z[5]*tmp1[20]
                    +z[78]*tmp1[21]+z[54]*tmp1[22]+z[30]*tmp1[23]+z[6]*tmp1[24]+z[79]*tmp1[25]+z[55]*tmp1[26]+z[31]*tmp1[27]+z[7]*tmp1[28]+z[80]*tmp1[29]+z[56]*tmp1[30]
                    +z[32]*tmp1[31]+z[8]*tmp1[32]+z[81]*tmp1[33]+z[57]*tmp1[34]+z[33]*tmp1[35]+z[9]*tmp1[36]+z[82]*tmp1[37]+z[58]*tmp1[38]+z[34]*tmp1[39]+z[10]*tmp1[40]
                    +z[83]*tmp1[41]+z[59]*tmp1[42]+z[35]*tmp1[43]+z[11]*tmp1[44]+z[84]*tmp1[45]+z[60]*tmp1[46]+z[36]*tmp1[47]+z[12]*tmp1[48]+z[85]*tmp1[49]+z[61]*tmp1[50]
                    +z[37]*tmp1[51]+z[13]*tmp1[52]+z[86]*tmp1[53]+z[62]*tmp1[54]+z[38]*tmp1[55]+z[14]*tmp1[56]+z[87]*tmp1[57]+z[63]*tmp1[58]+z[39]*tmp1[59]+z[15]*tmp1[60]
                    +z[88]*tmp1[61]+z[64]*tmp1[62]+z[40]*tmp1[63]+z[16]*tmp1[64]+z[89]*tmp1[65]+z[65]*tmp1[66]+z[41]*tmp1[67]+z[17]*tmp1[68]+z[90]*tmp1[69]+z[66]*tmp1[70]
                    +z[42]*tmp1[71]+z[18]*tmp1[72]+z[91]*tmp1[73]+z[67]*tmp1[74]+z[43]*tmp1[75]+z[19]*tmp1[76]+z[92]*tmp1[77]+z[68]*tmp1[78]+z[44]*tmp1[79]+z[20]*tmp1[80]
                    +z[93]*tmp1[81]+z[69]*tmp1[82]+z[45]*tmp1[83]+z[21]*tmp1[84]+z[94]*tmp1[85]+z[70]*tmp1[86]+z[46]*tmp1[87]+z[22]*tmp1[88]+z[95]*tmp1[89]+z[71]*tmp1[90]
                    +z[47]*tmp1[91]+z[23]*tmp1[92]+z[96]*tmp1[93]+z[72]*tmp1[94]+z[48]*tmp1[95]+z[24]*tmp1[96];
                    tab[nb_tmp3+74*nb3]=tab[nb_tmp3+74*nb3]+z[0]*tmp1[0]
                    +z[74]*tmp1[1]+z[51]*tmp1[2]+z[28]*tmp1[3]+z[5]*tmp1[4]+z[79]*tmp1[5]+z[56]*tmp1[6]+z[33]*tmp1[7]+z[10]*tmp1[8]+z[84]*tmp1[9]+z[61]*tmp1[10]
                    +z[38]*tmp1[11]+z[15]*tmp1[12]+z[89]*tmp1[13]+z[66]*tmp1[14]+z[43]*tmp1[15]+z[20]*tmp1[16]+z[94]*tmp1[17]+z[71]*tmp1[18]+z[48]*tmp1[19]+z[25]*tmp1[20]
                    +z[2]*tmp1[21]+z[76]*tmp1[22]+z[53]*tmp1[23]+z[30]*tmp1[24]+z[7]*tmp1[25]+z[81]*tmp1[26]+z[58]*tmp1[27]+z[35]*tmp1[28]+z[12]*tmp1[29]+z[86]*tmp1[30]
                    +z[63]*tmp1[31]+z[40]*tmp1[32]+z[17]*tmp1[33]+z[91]*tmp1[34]+z[68]*tmp1[35]+z[45]*tmp1[36]+z[22]*tmp1[37]+z[96]*tmp1[38]+z[73]*tmp1[39]+z[50]*tmp1[40]
                    +z[27]*tmp1[41]+z[4]*tmp1[42]+z[78]*tmp1[43]+z[55]*tmp1[44]+z[32]*tmp1[45]+z[9]*tmp1[46]+z[83]*tmp1[47]+z[60]*tmp1[48]+z[37]*tmp1[49]+z[14]*tmp1[50]
                    +z[88]*tmp1[51]+z[65]*tmp1[52]+z[42]*tmp1[53]+z[19]*tmp1[54]+z[93]*tmp1[55]+z[70]*tmp1[56]+z[47]*tmp1[57]+z[24]*tmp1[58]+z[1]*tmp1[59]+z[75]*tmp1[60]
                    +z[52]*tmp1[61]+z[29]*tmp1[62]+z[6]*tmp1[63]+z[80]*tmp1[64]+z[57]*tmp1[65]+z[34]*tmp1[66]+z[11]*tmp1[67]+z[85]*tmp1[68]+z[62]*tmp1[69]+z[39]*tmp1[70]
                    +z[16]*tmp1[71]+z[90]*tmp1[72]+z[67]*tmp1[73]+z[44]*tmp1[74]+z[21]*tmp1[75]+z[95]*tmp1[76]+z[72]*tmp1[77]+z[49]*tmp1[78]+z[26]*tmp1[79]+z[3]*tmp1[80]
                    +z[77]*tmp1[81]+z[54]*tmp1[82]+z[31]*tmp1[83]+z[8]*tmp1[84]+z[82]*tmp1[85]+z[59]*tmp1[86]+z[36]*tmp1[87]+z[13]*tmp1[88]+z[87]*tmp1[89]+z[64]*tmp1[90]
                    +z[41]*tmp1[91]+z[18]*tmp1[92]+z[92]*tmp1[93]+z[69]*tmp1[94]+z[46]*tmp1[95]+z[23]*tmp1[96];
                    tab[nb_tmp3+75*nb3]=tab[nb_tmp3+75*nb3]+z[0]*tmp1[0]
                    +z[75]*tmp1[1]+z[53]*tmp1[2]+z[31]*tmp1[3]+z[9]*tmp1[4]+z[84]*tmp1[5]+z[62]*tmp1[6]+z[40]*tmp1[7]+z[18]*tmp1[8]+z[93]*tmp1[9]+z[71]*tmp1[10]
                    +z[49]*tmp1[11]+z[27]*tmp1[12]+z[5]*tmp1[13]+z[80]*tmp1[14]+z[58]*tmp1[15]+z[36]*tmp1[16]+z[14]*tmp1[17]+z[89]*tmp1[18]+z[67]*tmp1[19]+z[45]*tmp1[20]
                    +z[23]*tmp1[21]+z[1]*tmp1[22]+z[76]*tmp1[23]+z[54]*tmp1[24]+z[32]*tmp1[25]+z[10]*tmp1[26]+z[85]*tmp1[27]+z[63]*tmp1[28]+z[41]*tmp1[29]+z[19]*tmp1[30]
                    +z[94]*tmp1[31]+z[72]*tmp1[32]+z[50]*tmp1[33]+z[28]*tmp1[34]+z[6]*tmp1[35]+z[81]*tmp1[36]+z[59]*tmp1[37]+z[37]*tmp1[38]+z[15]*tmp1[39]+z[90]*tmp1[40]
                    +z[68]*tmp1[41]+z[46]*tmp1[42]+z[24]*tmp1[43]+z[2]*tmp1[44]+z[77]*tmp1[45]+z[55]*tmp1[46]+z[33]*tmp1[47]+z[11]*tmp1[48]+z[86]*tmp1[49]+z[64]*tmp1[50]
                    +z[42]*tmp1[51]+z[20]*tmp1[52]+z[95]*tmp1[53]+z[73]*tmp1[54]+z[51]*tmp1[55]+z[29]*tmp1[56]+z[7]*tmp1[57]+z[82]*tmp1[58]+z[60]*tmp1[59]+z[38]*tmp1[60]
                    +z[16]*tmp1[61]+z[91]*tmp1[62]+z[69]*tmp1[63]+z[47]*tmp1[64]+z[25]*tmp1[65]+z[3]*tmp1[66]+z[78]*tmp1[67]+z[56]*tmp1[68]+z[34]*tmp1[69]+z[12]*tmp1[70]
                    +z[87]*tmp1[71]+z[65]*tmp1[72]+z[43]*tmp1[73]+z[21]*tmp1[74]+z[96]*tmp1[75]+z[74]*tmp1[76]+z[52]*tmp1[77]+z[30]*tmp1[78]+z[8]*tmp1[79]+z[83]*tmp1[80]
                    +z[61]*tmp1[81]+z[39]*tmp1[82]+z[17]*tmp1[83]+z[92]*tmp1[84]+z[70]*tmp1[85]+z[48]*tmp1[86]+z[26]*tmp1[87]+z[4]*tmp1[88]+z[79]*tmp1[89]+z[57]*tmp1[90]
                    +z[35]*tmp1[91]+z[13]*tmp1[92]+z[88]*tmp1[93]+z[66]*tmp1[94]+z[44]*tmp1[95]+z[22]*tmp1[96];
                    tab[nb_tmp3+76*nb3]=tab[nb_tmp3+76*nb3]+z[0]*tmp1[0]
                    +z[76]*tmp1[1]+z[55]*tmp1[2]+z[34]*tmp1[3]+z[13]*tmp1[4]+z[89]*tmp1[5]+z[68]*tmp1[6]+z[47]*tmp1[7]+z[26]*tmp1[8]+z[5]*tmp1[9]+z[81]*tmp1[10]
                    +z[60]*tmp1[11]+z[39]*tmp1[12]+z[18]*tmp1[13]+z[94]*tmp1[14]+z[73]*tmp1[15]+z[52]*tmp1[16]+z[31]*tmp1[17]+z[10]*tmp1[18]+z[86]*tmp1[19]+z[65]*tmp1[20]
                    +z[44]*tmp1[21]+z[23]*tmp1[22]+z[2]*tmp1[23]+z[78]*tmp1[24]+z[57]*tmp1[25]+z[36]*tmp1[26]+z[15]*tmp1[27]+z[91]*tmp1[28]+z[70]*tmp1[29]+z[49]*tmp1[30]
                    +z[28]*tmp1[31]+z[7]*tmp1[32]+z[83]*tmp1[33]+z[62]*tmp1[34]+z[41]*tmp1[35]+z[20]*tmp1[36]+z[96]*tmp1[37]+z[75]*tmp1[38]+z[54]*tmp1[39]+z[33]*tmp1[40]
                    +z[12]*tmp1[41]+z[88]*tmp1[42]+z[67]*tmp1[43]+z[46]*tmp1[44]+z[25]*tmp1[45]+z[4]*tmp1[46]+z[80]*tmp1[47]+z[59]*tmp1[48]+z[38]*tmp1[49]+z[17]*tmp1[50]
                    +z[93]*tmp1[51]+z[72]*tmp1[52]+z[51]*tmp1[53]+z[30]*tmp1[54]+z[9]*tmp1[55]+z[85]*tmp1[56]+z[64]*tmp1[57]+z[43]*tmp1[58]+z[22]*tmp1[59]+z[1]*tmp1[60]
                    +z[77]*tmp1[61]+z[56]*tmp1[62]+z[35]*tmp1[63]+z[14]*tmp1[64]+z[90]*tmp1[65]+z[69]*tmp1[66]+z[48]*tmp1[67]+z[27]*tmp1[68]+z[6]*tmp1[69]+z[82]*tmp1[70]
                    +z[61]*tmp1[71]+z[40]*tmp1[72]+z[19]*tmp1[73]+z[95]*tmp1[74]+z[74]*tmp1[75]+z[53]*tmp1[76]+z[32]*tmp1[77]+z[11]*tmp1[78]+z[87]*tmp1[79]+z[66]*tmp1[80]
                    +z[45]*tmp1[81]+z[24]*tmp1[82]+z[3]*tmp1[83]+z[79]*tmp1[84]+z[58]*tmp1[85]+z[37]*tmp1[86]+z[16]*tmp1[87]+z[92]*tmp1[88]+z[71]*tmp1[89]+z[50]*tmp1[90]
                    +z[29]*tmp1[91]+z[8]*tmp1[92]+z[84]*tmp1[93]+z[63]*tmp1[94]+z[42]*tmp1[95]+z[21]*tmp1[96];
                    tab[nb_tmp3+77*nb3]=tab[nb_tmp3+77*nb3]+z[0]*tmp1[0]
                    +z[77]*tmp1[1]+z[57]*tmp1[2]+z[37]*tmp1[3]+z[17]*tmp1[4]+z[94]*tmp1[5]+z[74]*tmp1[6]+z[54]*tmp1[7]+z[34]*tmp1[8]+z[14]*tmp1[9]+z[91]*tmp1[10]
                    +z[71]*tmp1[11]+z[51]*tmp1[12]+z[31]*tmp1[13]+z[11]*tmp1[14]+z[88]*tmp1[15]+z[68]*tmp1[16]+z[48]*tmp1[17]+z[28]*tmp1[18]+z[8]*tmp1[19]+z[85]*tmp1[20]
                    +z[65]*tmp1[21]+z[45]*tmp1[22]+z[25]*tmp1[23]+z[5]*tmp1[24]+z[82]*tmp1[25]+z[62]*tmp1[26]+z[42]*tmp1[27]+z[22]*tmp1[28]+z[2]*tmp1[29]+z[79]*tmp1[30]
                    +z[59]*tmp1[31]+z[39]*tmp1[32]+z[19]*tmp1[33]+z[96]*tmp1[34]+z[76]*tmp1[35]+z[56]*tmp1[36]+z[36]*tmp1[37]+z[16]*tmp1[38]+z[93]*tmp1[39]+z[73]*tmp1[40]
                    +z[53]*tmp1[41]+z[33]*tmp1[42]+z[13]*tmp1[43]+z[90]*tmp1[44]+z[70]*tmp1[45]+z[50]*tmp1[46]+z[30]*tmp1[47]+z[10]*tmp1[48]+z[87]*tmp1[49]+z[67]*tmp1[50]
                    +z[47]*tmp1[51]+z[27]*tmp1[52]+z[7]*tmp1[53]+z[84]*tmp1[54]+z[64]*tmp1[55]+z[44]*tmp1[56]+z[24]*tmp1[57]+z[4]*tmp1[58]+z[81]*tmp1[59]+z[61]*tmp1[60]
                    +z[41]*tmp1[61]+z[21]*tmp1[62]+z[1]*tmp1[63]+z[78]*tmp1[64]+z[58]*tmp1[65]+z[38]*tmp1[66]+z[18]*tmp1[67]+z[95]*tmp1[68]+z[75]*tmp1[69]+z[55]*tmp1[70]
                    +z[35]*tmp1[71]+z[15]*tmp1[72]+z[92]*tmp1[73]+z[72]*tmp1[74]+z[52]*tmp1[75]+z[32]*tmp1[76]+z[12]*tmp1[77]+z[89]*tmp1[78]+z[69]*tmp1[79]+z[49]*tmp1[80]
                    +z[29]*tmp1[81]+z[9]*tmp1[82]+z[86]*tmp1[83]+z[66]*tmp1[84]+z[46]*tmp1[85]+z[26]*tmp1[86]+z[6]*tmp1[87]+z[83]*tmp1[88]+z[63]*tmp1[89]+z[43]*tmp1[90]
                    +z[23]*tmp1[91]+z[3]*tmp1[92]+z[80]*tmp1[93]+z[60]*tmp1[94]+z[40]*tmp1[95]+z[20]*tmp1[96];
                    tab[nb_tmp3+78*nb3]=tab[nb_tmp3+78*nb3]+z[0]*tmp1[0]
                    +z[78]*tmp1[1]+z[59]*tmp1[2]+z[40]*tmp1[3]+z[21]*tmp1[4]+z[2]*tmp1[5]+z[80]*tmp1[6]+z[61]*tmp1[7]+z[42]*tmp1[8]+z[23]*tmp1[9]+z[4]*tmp1[10]
                    +z[82]*tmp1[11]+z[63]*tmp1[12]+z[44]*tmp1[13]+z[25]*tmp1[14]+z[6]*tmp1[15]+z[84]*tmp1[16]+z[65]*tmp1[17]+z[46]*tmp1[18]+z[27]*tmp1[19]+z[8]*tmp1[20]
                    +z[86]*tmp1[21]+z[67]*tmp1[22]+z[48]*tmp1[23]+z[29]*tmp1[24]+z[10]*tmp1[25]+z[88]*tmp1[26]+z[69]*tmp1[27]+z[50]*tmp1[28]+z[31]*tmp1[29]+z[12]*tmp1[30]
                    +z[90]*tmp1[31]+z[71]*tmp1[32]+z[52]*tmp1[33]+z[33]*tmp1[34]+z[14]*tmp1[35]+z[92]*tmp1[36]+z[73]*tmp1[37]+z[54]*tmp1[38]+z[35]*tmp1[39]+z[16]*tmp1[40]
                    +z[94]*tmp1[41]+z[75]*tmp1[42]+z[56]*tmp1[43]+z[37]*tmp1[44]+z[18]*tmp1[45]+z[96]*tmp1[46]+z[77]*tmp1[47]+z[58]*tmp1[48]+z[39]*tmp1[49]+z[20]*tmp1[50]
                    +z[1]*tmp1[51]+z[79]*tmp1[52]+z[60]*tmp1[53]+z[41]*tmp1[54]+z[22]*tmp1[55]+z[3]*tmp1[56]+z[81]*tmp1[57]+z[62]*tmp1[58]+z[43]*tmp1[59]+z[24]*tmp1[60]
                    +z[5]*tmp1[61]+z[83]*tmp1[62]+z[64]*tmp1[63]+z[45]*tmp1[64]+z[26]*tmp1[65]+z[7]*tmp1[66]+z[85]*tmp1[67]+z[66]*tmp1[68]+z[47]*tmp1[69]+z[28]*tmp1[70]
                    +z[9]*tmp1[71]+z[87]*tmp1[72]+z[68]*tmp1[73]+z[49]*tmp1[74]+z[30]*tmp1[75]+z[11]*tmp1[76]+z[89]*tmp1[77]+z[70]*tmp1[78]+z[51]*tmp1[79]+z[32]*tmp1[80]
                    +z[13]*tmp1[81]+z[91]*tmp1[82]+z[72]*tmp1[83]+z[53]*tmp1[84]+z[34]*tmp1[85]+z[15]*tmp1[86]+z[93]*tmp1[87]+z[74]*tmp1[88]+z[55]*tmp1[89]+z[36]*tmp1[90]
                    +z[17]*tmp1[91]+z[95]*tmp1[92]+z[76]*tmp1[93]+z[57]*tmp1[94]+z[38]*tmp1[95]+z[19]*tmp1[96];
                    tab[nb_tmp3+79*nb3]=tab[nb_tmp3+79*nb3]+z[0]*tmp1[0]
                    +z[79]*tmp1[1]+z[61]*tmp1[2]+z[43]*tmp1[3]+z[25]*tmp1[4]+z[7]*tmp1[5]+z[86]*tmp1[6]+z[68]*tmp1[7]+z[50]*tmp1[8]+z[32]*tmp1[9]+z[14]*tmp1[10]
                    +z[93]*tmp1[11]+z[75]*tmp1[12]+z[57]*tmp1[13]+z[39]*tmp1[14]+z[21]*tmp1[15]+z[3]*tmp1[16]+z[82]*tmp1[17]+z[64]*tmp1[18]+z[46]*tmp1[19]+z[28]*tmp1[20]
                    +z[10]*tmp1[21]+z[89]*tmp1[22]+z[71]*tmp1[23]+z[53]*tmp1[24]+z[35]*tmp1[25]+z[17]*tmp1[26]+z[96]*tmp1[27]+z[78]*tmp1[28]+z[60]*tmp1[29]+z[42]*tmp1[30]
                    +z[24]*tmp1[31]+z[6]*tmp1[32]+z[85]*tmp1[33]+z[67]*tmp1[34]+z[49]*tmp1[35]+z[31]*tmp1[36]+z[13]*tmp1[37]+z[92]*tmp1[38]+z[74]*tmp1[39]+z[56]*tmp1[40]
                    +z[38]*tmp1[41]+z[20]*tmp1[42]+z[2]*tmp1[43]+z[81]*tmp1[44]+z[63]*tmp1[45]+z[45]*tmp1[46]+z[27]*tmp1[47]+z[9]*tmp1[48]+z[88]*tmp1[49]+z[70]*tmp1[50]
                    +z[52]*tmp1[51]+z[34]*tmp1[52]+z[16]*tmp1[53]+z[95]*tmp1[54]+z[77]*tmp1[55]+z[59]*tmp1[56]+z[41]*tmp1[57]+z[23]*tmp1[58]+z[5]*tmp1[59]+z[84]*tmp1[60]
                    +z[66]*tmp1[61]+z[48]*tmp1[62]+z[30]*tmp1[63]+z[12]*tmp1[64]+z[91]*tmp1[65]+z[73]*tmp1[66]+z[55]*tmp1[67]+z[37]*tmp1[68]+z[19]*tmp1[69]+z[1]*tmp1[70]
                    +z[80]*tmp1[71]+z[62]*tmp1[72]+z[44]*tmp1[73]+z[26]*tmp1[74]+z[8]*tmp1[75]+z[87]*tmp1[76]+z[69]*tmp1[77]+z[51]*tmp1[78]+z[33]*tmp1[79]+z[15]*tmp1[80]
                    +z[94]*tmp1[81]+z[76]*tmp1[82]+z[58]*tmp1[83]+z[40]*tmp1[84]+z[22]*tmp1[85]+z[4]*tmp1[86]+z[83]*tmp1[87]+z[65]*tmp1[88]+z[47]*tmp1[89]+z[29]*tmp1[90]
                    +z[11]*tmp1[91]+z[90]*tmp1[92]+z[72]*tmp1[93]+z[54]*tmp1[94]+z[36]*tmp1[95]+z[18]*tmp1[96];
                    tab[nb_tmp3+80*nb3]=tab[nb_tmp3+80*nb3]+z[0]*tmp1[0]
                    +z[80]*tmp1[1]+z[63]*tmp1[2]+z[46]*tmp1[3]+z[29]*tmp1[4]+z[12]*tmp1[5]+z[92]*tmp1[6]+z[75]*tmp1[7]+z[58]*tmp1[8]+z[41]*tmp1[9]+z[24]*tmp1[10]
                    +z[7]*tmp1[11]+z[87]*tmp1[12]+z[70]*tmp1[13]+z[53]*tmp1[14]+z[36]*tmp1[15]+z[19]*tmp1[16]+z[2]*tmp1[17]+z[82]*tmp1[18]+z[65]*tmp1[19]+z[48]*tmp1[20]
                    +z[31]*tmp1[21]+z[14]*tmp1[22]+z[94]*tmp1[23]+z[77]*tmp1[24]+z[60]*tmp1[25]+z[43]*tmp1[26]+z[26]*tmp1[27]+z[9]*tmp1[28]+z[89]*tmp1[29]+z[72]*tmp1[30]
                    +z[55]*tmp1[31]+z[38]*tmp1[32]+z[21]*tmp1[33]+z[4]*tmp1[34]+z[84]*tmp1[35]+z[67]*tmp1[36]+z[50]*tmp1[37]+z[33]*tmp1[38]+z[16]*tmp1[39]+z[96]*tmp1[40]
                    +z[79]*tmp1[41]+z[62]*tmp1[42]+z[45]*tmp1[43]+z[28]*tmp1[44]+z[11]*tmp1[45]+z[91]*tmp1[46]+z[74]*tmp1[47]+z[57]*tmp1[48]+z[40]*tmp1[49]+z[23]*tmp1[50]
                    +z[6]*tmp1[51]+z[86]*tmp1[52]+z[69]*tmp1[53]+z[52]*tmp1[54]+z[35]*tmp1[55]+z[18]*tmp1[56]+z[1]*tmp1[57]+z[81]*tmp1[58]+z[64]*tmp1[59]+z[47]*tmp1[60]
                    +z[30]*tmp1[61]+z[13]*tmp1[62]+z[93]*tmp1[63]+z[76]*tmp1[64]+z[59]*tmp1[65]+z[42]*tmp1[66]+z[25]*tmp1[67]+z[8]*tmp1[68]+z[88]*tmp1[69]+z[71]*tmp1[70]
                    +z[54]*tmp1[71]+z[37]*tmp1[72]+z[20]*tmp1[73]+z[3]*tmp1[74]+z[83]*tmp1[75]+z[66]*tmp1[76]+z[49]*tmp1[77]+z[32]*tmp1[78]+z[15]*tmp1[79]+z[95]*tmp1[80]
                    +z[78]*tmp1[81]+z[61]*tmp1[82]+z[44]*tmp1[83]+z[27]*tmp1[84]+z[10]*tmp1[85]+z[90]*tmp1[86]+z[73]*tmp1[87]+z[56]*tmp1[88]+z[39]*tmp1[89]+z[22]*tmp1[90]
                    +z[5]*tmp1[91]+z[85]*tmp1[92]+z[68]*tmp1[93]+z[51]*tmp1[94]+z[34]*tmp1[95]+z[17]*tmp1[96];
                    tab[nb_tmp3+81*nb3]=tab[nb_tmp3+81*nb3]+z[0]*tmp1[0]
                    +z[81]*tmp1[1]+z[65]*tmp1[2]+z[49]*tmp1[3]+z[33]*tmp1[4]+z[17]*tmp1[5]+z[1]*tmp1[6]+z[82]*tmp1[7]+z[66]*tmp1[8]+z[50]*tmp1[9]+z[34]*tmp1[10]
                    +z[18]*tmp1[11]+z[2]*tmp1[12]+z[83]*tmp1[13]+z[67]*tmp1[14]+z[51]*tmp1[15]+z[35]*tmp1[16]+z[19]*tmp1[17]+z[3]*tmp1[18]+z[84]*tmp1[19]+z[68]*tmp1[20]
                    +z[52]*tmp1[21]+z[36]*tmp1[22]+z[20]*tmp1[23]+z[4]*tmp1[24]+z[85]*tmp1[25]+z[69]*tmp1[26]+z[53]*tmp1[27]+z[37]*tmp1[28]+z[21]*tmp1[29]+z[5]*tmp1[30]
                    +z[86]*tmp1[31]+z[70]*tmp1[32]+z[54]*tmp1[33]+z[38]*tmp1[34]+z[22]*tmp1[35]+z[6]*tmp1[36]+z[87]*tmp1[37]+z[71]*tmp1[38]+z[55]*tmp1[39]+z[39]*tmp1[40]
                    +z[23]*tmp1[41]+z[7]*tmp1[42]+z[88]*tmp1[43]+z[72]*tmp1[44]+z[56]*tmp1[45]+z[40]*tmp1[46]+z[24]*tmp1[47]+z[8]*tmp1[48]+z[89]*tmp1[49]+z[73]*tmp1[50]
                    +z[57]*tmp1[51]+z[41]*tmp1[52]+z[25]*tmp1[53]+z[9]*tmp1[54]+z[90]*tmp1[55]+z[74]*tmp1[56]+z[58]*tmp1[57]+z[42]*tmp1[58]+z[26]*tmp1[59]+z[10]*tmp1[60]
                    +z[91]*tmp1[61]+z[75]*tmp1[62]+z[59]*tmp1[63]+z[43]*tmp1[64]+z[27]*tmp1[65]+z[11]*tmp1[66]+z[92]*tmp1[67]+z[76]*tmp1[68]+z[60]*tmp1[69]+z[44]*tmp1[70]
                    +z[28]*tmp1[71]+z[12]*tmp1[72]+z[93]*tmp1[73]+z[77]*tmp1[74]+z[61]*tmp1[75]+z[45]*tmp1[76]+z[29]*tmp1[77]+z[13]*tmp1[78]+z[94]*tmp1[79]+z[78]*tmp1[80]
                    +z[62]*tmp1[81]+z[46]*tmp1[82]+z[30]*tmp1[83]+z[14]*tmp1[84]+z[95]*tmp1[85]+z[79]*tmp1[86]+z[63]*tmp1[87]+z[47]*tmp1[88]+z[31]*tmp1[89]+z[15]*tmp1[90]
                    +z[96]*tmp1[91]+z[80]*tmp1[92]+z[64]*tmp1[93]+z[48]*tmp1[94]+z[32]*tmp1[95]+z[16]*tmp1[96];
                    tab[nb_tmp3+82*nb3]=tab[nb_tmp3+82*nb3]+z[0]*tmp1[0]
                    +z[82]*tmp1[1]+z[67]*tmp1[2]+z[52]*tmp1[3]+z[37]*tmp1[4]+z[22]*tmp1[5]+z[7]*tmp1[6]+z[89]*tmp1[7]+z[74]*tmp1[8]+z[59]*tmp1[9]+z[44]*tmp1[10]
                    +z[29]*tmp1[11]+z[14]*tmp1[12]+z[96]*tmp1[13]+z[81]*tmp1[14]+z[66]*tmp1[15]+z[51]*tmp1[16]+z[36]*tmp1[17]+z[21]*tmp1[18]+z[6]*tmp1[19]+z[88]*tmp1[20]
                    +z[73]*tmp1[21]+z[58]*tmp1[22]+z[43]*tmp1[23]+z[28]*tmp1[24]+z[13]*tmp1[25]+z[95]*tmp1[26]+z[80]*tmp1[27]+z[65]*tmp1[28]+z[50]*tmp1[29]+z[35]*tmp1[30]
                    +z[20]*tmp1[31]+z[5]*tmp1[32]+z[87]*tmp1[33]+z[72]*tmp1[34]+z[57]*tmp1[35]+z[42]*tmp1[36]+z[27]*tmp1[37]+z[12]*tmp1[38]+z[94]*tmp1[39]+z[79]*tmp1[40]
                    +z[64]*tmp1[41]+z[49]*tmp1[42]+z[34]*tmp1[43]+z[19]*tmp1[44]+z[4]*tmp1[45]+z[86]*tmp1[46]+z[71]*tmp1[47]+z[56]*tmp1[48]+z[41]*tmp1[49]+z[26]*tmp1[50]
                    +z[11]*tmp1[51]+z[93]*tmp1[52]+z[78]*tmp1[53]+z[63]*tmp1[54]+z[48]*tmp1[55]+z[33]*tmp1[56]+z[18]*tmp1[57]+z[3]*tmp1[58]+z[85]*tmp1[59]+z[70]*tmp1[60]
                    +z[55]*tmp1[61]+z[40]*tmp1[62]+z[25]*tmp1[63]+z[10]*tmp1[64]+z[92]*tmp1[65]+z[77]*tmp1[66]+z[62]*tmp1[67]+z[47]*tmp1[68]+z[32]*tmp1[69]+z[17]*tmp1[70]
                    +z[2]*tmp1[71]+z[84]*tmp1[72]+z[69]*tmp1[73]+z[54]*tmp1[74]+z[39]*tmp1[75]+z[24]*tmp1[76]+z[9]*tmp1[77]+z[91]*tmp1[78]+z[76]*tmp1[79]+z[61]*tmp1[80]
                    +z[46]*tmp1[81]+z[31]*tmp1[82]+z[16]*tmp1[83]+z[1]*tmp1[84]+z[83]*tmp1[85]+z[68]*tmp1[86]+z[53]*tmp1[87]+z[38]*tmp1[88]+z[23]*tmp1[89]+z[8]*tmp1[90]
                    +z[90]*tmp1[91]+z[75]*tmp1[92]+z[60]*tmp1[93]+z[45]*tmp1[94]+z[30]*tmp1[95]+z[15]*tmp1[96];
                    tab[nb_tmp3+83*nb3]=tab[nb_tmp3+83*nb3]+z[0]*tmp1[0]
                    +z[83]*tmp1[1]+z[69]*tmp1[2]+z[55]*tmp1[3]+z[41]*tmp1[4]+z[27]*tmp1[5]+z[13]*tmp1[6]+z[96]*tmp1[7]+z[82]*tmp1[8]+z[68]*tmp1[9]+z[54]*tmp1[10]
                    +z[40]*tmp1[11]+z[26]*tmp1[12]+z[12]*tmp1[13]+z[95]*tmp1[14]+z[81]*tmp1[15]+z[67]*tmp1[16]+z[53]*tmp1[17]+z[39]*tmp1[18]+z[25]*tmp1[19]+z[11]*tmp1[20]
                    +z[94]*tmp1[21]+z[80]*tmp1[22]+z[66]*tmp1[23]+z[52]*tmp1[24]+z[38]*tmp1[25]+z[24]*tmp1[26]+z[10]*tmp1[27]+z[93]*tmp1[28]+z[79]*tmp1[29]+z[65]*tmp1[30]
                    +z[51]*tmp1[31]+z[37]*tmp1[32]+z[23]*tmp1[33]+z[9]*tmp1[34]+z[92]*tmp1[35]+z[78]*tmp1[36]+z[64]*tmp1[37]+z[50]*tmp1[38]+z[36]*tmp1[39]+z[22]*tmp1[40]
                    +z[8]*tmp1[41]+z[91]*tmp1[42]+z[77]*tmp1[43]+z[63]*tmp1[44]+z[49]*tmp1[45]+z[35]*tmp1[46]+z[21]*tmp1[47]+z[7]*tmp1[48]+z[90]*tmp1[49]+z[76]*tmp1[50]
                    +z[62]*tmp1[51]+z[48]*tmp1[52]+z[34]*tmp1[53]+z[20]*tmp1[54]+z[6]*tmp1[55]+z[89]*tmp1[56]+z[75]*tmp1[57]+z[61]*tmp1[58]+z[47]*tmp1[59]+z[33]*tmp1[60]
                    +z[19]*tmp1[61]+z[5]*tmp1[62]+z[88]*tmp1[63]+z[74]*tmp1[64]+z[60]*tmp1[65]+z[46]*tmp1[66]+z[32]*tmp1[67]+z[18]*tmp1[68]+z[4]*tmp1[69]+z[87]*tmp1[70]
                    +z[73]*tmp1[71]+z[59]*tmp1[72]+z[45]*tmp1[73]+z[31]*tmp1[74]+z[17]*tmp1[75]+z[3]*tmp1[76]+z[86]*tmp1[77]+z[72]*tmp1[78]+z[58]*tmp1[79]+z[44]*tmp1[80]
                    +z[30]*tmp1[81]+z[16]*tmp1[82]+z[2]*tmp1[83]+z[85]*tmp1[84]+z[71]*tmp1[85]+z[57]*tmp1[86]+z[43]*tmp1[87]+z[29]*tmp1[88]+z[15]*tmp1[89]+z[1]*tmp1[90]
                    +z[84]*tmp1[91]+z[70]*tmp1[92]+z[56]*tmp1[93]+z[42]*tmp1[94]+z[28]*tmp1[95]+z[14]*tmp1[96];
                    tab[nb_tmp3+84*nb3]=tab[nb_tmp3+84*nb3]+z[0]*tmp1[0]
                    +z[84]*tmp1[1]+z[71]*tmp1[2]+z[58]*tmp1[3]+z[45]*tmp1[4]+z[32]*tmp1[5]+z[19]*tmp1[6]+z[6]*tmp1[7]+z[90]*tmp1[8]+z[77]*tmp1[9]+z[64]*tmp1[10]
                    +z[51]*tmp1[11]+z[38]*tmp1[12]+z[25]*tmp1[13]+z[12]*tmp1[14]+z[96]*tmp1[15]+z[83]*tmp1[16]+z[70]*tmp1[17]+z[57]*tmp1[18]+z[44]*tmp1[19]+z[31]*tmp1[20]
                    +z[18]*tmp1[21]+z[5]*tmp1[22]+z[89]*tmp1[23]+z[76]*tmp1[24]+z[63]*tmp1[25]+z[50]*tmp1[26]+z[37]*tmp1[27]+z[24]*tmp1[28]+z[11]*tmp1[29]+z[95]*tmp1[30]
                    +z[82]*tmp1[31]+z[69]*tmp1[32]+z[56]*tmp1[33]+z[43]*tmp1[34]+z[30]*tmp1[35]+z[17]*tmp1[36]+z[4]*tmp1[37]+z[88]*tmp1[38]+z[75]*tmp1[39]+z[62]*tmp1[40]
                    +z[49]*tmp1[41]+z[36]*tmp1[42]+z[23]*tmp1[43]+z[10]*tmp1[44]+z[94]*tmp1[45]+z[81]*tmp1[46]+z[68]*tmp1[47]+z[55]*tmp1[48]+z[42]*tmp1[49]+z[29]*tmp1[50]
                    +z[16]*tmp1[51]+z[3]*tmp1[52]+z[87]*tmp1[53]+z[74]*tmp1[54]+z[61]*tmp1[55]+z[48]*tmp1[56]+z[35]*tmp1[57]+z[22]*tmp1[58]+z[9]*tmp1[59]+z[93]*tmp1[60]
                    +z[80]*tmp1[61]+z[67]*tmp1[62]+z[54]*tmp1[63]+z[41]*tmp1[64]+z[28]*tmp1[65]+z[15]*tmp1[66]+z[2]*tmp1[67]+z[86]*tmp1[68]+z[73]*tmp1[69]+z[60]*tmp1[70]
                    +z[47]*tmp1[71]+z[34]*tmp1[72]+z[21]*tmp1[73]+z[8]*tmp1[74]+z[92]*tmp1[75]+z[79]*tmp1[76]+z[66]*tmp1[77]+z[53]*tmp1[78]+z[40]*tmp1[79]+z[27]*tmp1[80]
                    +z[14]*tmp1[81]+z[1]*tmp1[82]+z[85]*tmp1[83]+z[72]*tmp1[84]+z[59]*tmp1[85]+z[46]*tmp1[86]+z[33]*tmp1[87]+z[20]*tmp1[88]+z[7]*tmp1[89]+z[91]*tmp1[90]
                    +z[78]*tmp1[91]+z[65]*tmp1[92]+z[52]*tmp1[93]+z[39]*tmp1[94]+z[26]*tmp1[95]+z[13]*tmp1[96];
                    tab[nb_tmp3+85*nb3]=tab[nb_tmp3+85*nb3]+z[0]*tmp1[0]
                    +z[85]*tmp1[1]+z[73]*tmp1[2]+z[61]*tmp1[3]+z[49]*tmp1[4]+z[37]*tmp1[5]+z[25]*tmp1[6]+z[13]*tmp1[7]+z[1]*tmp1[8]+z[86]*tmp1[9]+z[74]*tmp1[10]
                    +z[62]*tmp1[11]+z[50]*tmp1[12]+z[38]*tmp1[13]+z[26]*tmp1[14]+z[14]*tmp1[15]+z[2]*tmp1[16]+z[87]*tmp1[17]+z[75]*tmp1[18]+z[63]*tmp1[19]+z[51]*tmp1[20]
                    +z[39]*tmp1[21]+z[27]*tmp1[22]+z[15]*tmp1[23]+z[3]*tmp1[24]+z[88]*tmp1[25]+z[76]*tmp1[26]+z[64]*tmp1[27]+z[52]*tmp1[28]+z[40]*tmp1[29]+z[28]*tmp1[30]
                    +z[16]*tmp1[31]+z[4]*tmp1[32]+z[89]*tmp1[33]+z[77]*tmp1[34]+z[65]*tmp1[35]+z[53]*tmp1[36]+z[41]*tmp1[37]+z[29]*tmp1[38]+z[17]*tmp1[39]+z[5]*tmp1[40]
                    +z[90]*tmp1[41]+z[78]*tmp1[42]+z[66]*tmp1[43]+z[54]*tmp1[44]+z[42]*tmp1[45]+z[30]*tmp1[46]+z[18]*tmp1[47]+z[6]*tmp1[48]+z[91]*tmp1[49]+z[79]*tmp1[50]
                    +z[67]*tmp1[51]+z[55]*tmp1[52]+z[43]*tmp1[53]+z[31]*tmp1[54]+z[19]*tmp1[55]+z[7]*tmp1[56]+z[92]*tmp1[57]+z[80]*tmp1[58]+z[68]*tmp1[59]+z[56]*tmp1[60]
                    +z[44]*tmp1[61]+z[32]*tmp1[62]+z[20]*tmp1[63]+z[8]*tmp1[64]+z[93]*tmp1[65]+z[81]*tmp1[66]+z[69]*tmp1[67]+z[57]*tmp1[68]+z[45]*tmp1[69]+z[33]*tmp1[70]
                    +z[21]*tmp1[71]+z[9]*tmp1[72]+z[94]*tmp1[73]+z[82]*tmp1[74]+z[70]*tmp1[75]+z[58]*tmp1[76]+z[46]*tmp1[77]+z[34]*tmp1[78]+z[22]*tmp1[79]+z[10]*tmp1[80]
                    +z[95]*tmp1[81]+z[83]*tmp1[82]+z[71]*tmp1[83]+z[59]*tmp1[84]+z[47]*tmp1[85]+z[35]*tmp1[86]+z[23]*tmp1[87]+z[11]*tmp1[88]+z[96]*tmp1[89]+z[84]*tmp1[90]
                    +z[72]*tmp1[91]+z[60]*tmp1[92]+z[48]*tmp1[93]+z[36]*tmp1[94]+z[24]*tmp1[95]+z[12]*tmp1[96];
                    tab[nb_tmp3+86*nb3]=tab[nb_tmp3+86*nb3]+z[0]*tmp1[0]
                    +z[86]*tmp1[1]+z[75]*tmp1[2]+z[64]*tmp1[3]+z[53]*tmp1[4]+z[42]*tmp1[5]+z[31]*tmp1[6]+z[20]*tmp1[7]+z[9]*tmp1[8]+z[95]*tmp1[9]+z[84]*tmp1[10]
                    +z[73]*tmp1[11]+z[62]*tmp1[12]+z[51]*tmp1[13]+z[40]*tmp1[14]+z[29]*tmp1[15]+z[18]*tmp1[16]+z[7]*tmp1[17]+z[93]*tmp1[18]+z[82]*tmp1[19]+z[71]*tmp1[20]
                    +z[60]*tmp1[21]+z[49]*tmp1[22]+z[38]*tmp1[23]+z[27]*tmp1[24]+z[16]*tmp1[25]+z[5]*tmp1[26]+z[91]*tmp1[27]+z[80]*tmp1[28]+z[69]*tmp1[29]+z[58]*tmp1[30]
                    +z[47]*tmp1[31]+z[36]*tmp1[32]+z[25]*tmp1[33]+z[14]*tmp1[34]+z[3]*tmp1[35]+z[89]*tmp1[36]+z[78]*tmp1[37]+z[67]*tmp1[38]+z[56]*tmp1[39]+z[45]*tmp1[40]
                    +z[34]*tmp1[41]+z[23]*tmp1[42]+z[12]*tmp1[43]+z[1]*tmp1[44]+z[87]*tmp1[45]+z[76]*tmp1[46]+z[65]*tmp1[47]+z[54]*tmp1[48]+z[43]*tmp1[49]+z[32]*tmp1[50]
                    +z[21]*tmp1[51]+z[10]*tmp1[52]+z[96]*tmp1[53]+z[85]*tmp1[54]+z[74]*tmp1[55]+z[63]*tmp1[56]+z[52]*tmp1[57]+z[41]*tmp1[58]+z[30]*tmp1[59]+z[19]*tmp1[60]
                    +z[8]*tmp1[61]+z[94]*tmp1[62]+z[83]*tmp1[63]+z[72]*tmp1[64]+z[61]*tmp1[65]+z[50]*tmp1[66]+z[39]*tmp1[67]+z[28]*tmp1[68]+z[17]*tmp1[69]+z[6]*tmp1[70]
                    +z[92]*tmp1[71]+z[81]*tmp1[72]+z[70]*tmp1[73]+z[59]*tmp1[74]+z[48]*tmp1[75]+z[37]*tmp1[76]+z[26]*tmp1[77]+z[15]*tmp1[78]+z[4]*tmp1[79]+z[90]*tmp1[80]
                    +z[79]*tmp1[81]+z[68]*tmp1[82]+z[57]*tmp1[83]+z[46]*tmp1[84]+z[35]*tmp1[85]+z[24]*tmp1[86]+z[13]*tmp1[87]+z[2]*tmp1[88]+z[88]*tmp1[89]+z[77]*tmp1[90]
                    +z[66]*tmp1[91]+z[55]*tmp1[92]+z[44]*tmp1[93]+z[33]*tmp1[94]+z[22]*tmp1[95]+z[11]*tmp1[96];
                    tab[nb_tmp3+87*nb3]=tab[nb_tmp3+87*nb3]+z[0]*tmp1[0]
                    +z[87]*tmp1[1]+z[77]*tmp1[2]+z[67]*tmp1[3]+z[57]*tmp1[4]+z[47]*tmp1[5]+z[37]*tmp1[6]+z[27]*tmp1[7]+z[17]*tmp1[8]+z[7]*tmp1[9]+z[94]*tmp1[10]
                    +z[84]*tmp1[11]+z[74]*tmp1[12]+z[64]*tmp1[13]+z[54]*tmp1[14]+z[44]*tmp1[15]+z[34]*tmp1[16]+z[24]*tmp1[17]+z[14]*tmp1[18]+z[4]*tmp1[19]+z[91]*tmp1[20]
                    +z[81]*tmp1[21]+z[71]*tmp1[22]+z[61]*tmp1[23]+z[51]*tmp1[24]+z[41]*tmp1[25]+z[31]*tmp1[26]+z[21]*tmp1[27]+z[11]*tmp1[28]+z[1]*tmp1[29]+z[88]*tmp1[30]
                    +z[78]*tmp1[31]+z[68]*tmp1[32]+z[58]*tmp1[33]+z[48]*tmp1[34]+z[38]*tmp1[35]+z[28]*tmp1[36]+z[18]*tmp1[37]+z[8]*tmp1[38]+z[95]*tmp1[39]+z[85]*tmp1[40]
                    +z[75]*tmp1[41]+z[65]*tmp1[42]+z[55]*tmp1[43]+z[45]*tmp1[44]+z[35]*tmp1[45]+z[25]*tmp1[46]+z[15]*tmp1[47]+z[5]*tmp1[48]+z[92]*tmp1[49]+z[82]*tmp1[50]
                    +z[72]*tmp1[51]+z[62]*tmp1[52]+z[52]*tmp1[53]+z[42]*tmp1[54]+z[32]*tmp1[55]+z[22]*tmp1[56]+z[12]*tmp1[57]+z[2]*tmp1[58]+z[89]*tmp1[59]+z[79]*tmp1[60]
                    +z[69]*tmp1[61]+z[59]*tmp1[62]+z[49]*tmp1[63]+z[39]*tmp1[64]+z[29]*tmp1[65]+z[19]*tmp1[66]+z[9]*tmp1[67]+z[96]*tmp1[68]+z[86]*tmp1[69]+z[76]*tmp1[70]
                    +z[66]*tmp1[71]+z[56]*tmp1[72]+z[46]*tmp1[73]+z[36]*tmp1[74]+z[26]*tmp1[75]+z[16]*tmp1[76]+z[6]*tmp1[77]+z[93]*tmp1[78]+z[83]*tmp1[79]+z[73]*tmp1[80]
                    +z[63]*tmp1[81]+z[53]*tmp1[82]+z[43]*tmp1[83]+z[33]*tmp1[84]+z[23]*tmp1[85]+z[13]*tmp1[86]+z[3]*tmp1[87]+z[90]*tmp1[88]+z[80]*tmp1[89]+z[70]*tmp1[90]
                    +z[60]*tmp1[91]+z[50]*tmp1[92]+z[40]*tmp1[93]+z[30]*tmp1[94]+z[20]*tmp1[95]+z[10]*tmp1[96];
                    tab[nb_tmp3+88*nb3]=tab[nb_tmp3+88*nb3]+z[0]*tmp1[0]
                    +z[88]*tmp1[1]+z[79]*tmp1[2]+z[70]*tmp1[3]+z[61]*tmp1[4]+z[52]*tmp1[5]+z[43]*tmp1[6]+z[34]*tmp1[7]+z[25]*tmp1[8]+z[16]*tmp1[9]+z[7]*tmp1[10]
                    +z[95]*tmp1[11]+z[86]*tmp1[12]+z[77]*tmp1[13]+z[68]*tmp1[14]+z[59]*tmp1[15]+z[50]*tmp1[16]+z[41]*tmp1[17]+z[32]*tmp1[18]+z[23]*tmp1[19]+z[14]*tmp1[20]
                    +z[5]*tmp1[21]+z[93]*tmp1[22]+z[84]*tmp1[23]+z[75]*tmp1[24]+z[66]*tmp1[25]+z[57]*tmp1[26]+z[48]*tmp1[27]+z[39]*tmp1[28]+z[30]*tmp1[29]+z[21]*tmp1[30]
                    +z[12]*tmp1[31]+z[3]*tmp1[32]+z[91]*tmp1[33]+z[82]*tmp1[34]+z[73]*tmp1[35]+z[64]*tmp1[36]+z[55]*tmp1[37]+z[46]*tmp1[38]+z[37]*tmp1[39]+z[28]*tmp1[40]
                    +z[19]*tmp1[41]+z[10]*tmp1[42]+z[1]*tmp1[43]+z[89]*tmp1[44]+z[80]*tmp1[45]+z[71]*tmp1[46]+z[62]*tmp1[47]+z[53]*tmp1[48]+z[44]*tmp1[49]+z[35]*tmp1[50]
                    +z[26]*tmp1[51]+z[17]*tmp1[52]+z[8]*tmp1[53]+z[96]*tmp1[54]+z[87]*tmp1[55]+z[78]*tmp1[56]+z[69]*tmp1[57]+z[60]*tmp1[58]+z[51]*tmp1[59]+z[42]*tmp1[60]
                    +z[33]*tmp1[61]+z[24]*tmp1[62]+z[15]*tmp1[63]+z[6]*tmp1[64]+z[94]*tmp1[65]+z[85]*tmp1[66]+z[76]*tmp1[67]+z[67]*tmp1[68]+z[58]*tmp1[69]+z[49]*tmp1[70]
                    +z[40]*tmp1[71]+z[31]*tmp1[72]+z[22]*tmp1[73]+z[13]*tmp1[74]+z[4]*tmp1[75]+z[92]*tmp1[76]+z[83]*tmp1[77]+z[74]*tmp1[78]+z[65]*tmp1[79]+z[56]*tmp1[80]
                    +z[47]*tmp1[81]+z[38]*tmp1[82]+z[29]*tmp1[83]+z[20]*tmp1[84]+z[11]*tmp1[85]+z[2]*tmp1[86]+z[90]*tmp1[87]+z[81]*tmp1[88]+z[72]*tmp1[89]+z[63]*tmp1[90]
                    +z[54]*tmp1[91]+z[45]*tmp1[92]+z[36]*tmp1[93]+z[27]*tmp1[94]+z[18]*tmp1[95]+z[9]*tmp1[96];
                    tab[nb_tmp3+89*nb3]=tab[nb_tmp3+89*nb3]+z[0]*tmp1[0]
                    +z[89]*tmp1[1]+z[81]*tmp1[2]+z[73]*tmp1[3]+z[65]*tmp1[4]+z[57]*tmp1[5]+z[49]*tmp1[6]+z[41]*tmp1[7]+z[33]*tmp1[8]+z[25]*tmp1[9]+z[17]*tmp1[10]
                    +z[9]*tmp1[11]+z[1]*tmp1[12]+z[90]*tmp1[13]+z[82]*tmp1[14]+z[74]*tmp1[15]+z[66]*tmp1[16]+z[58]*tmp1[17]+z[50]*tmp1[18]+z[42]*tmp1[19]+z[34]*tmp1[20]
                    +z[26]*tmp1[21]+z[18]*tmp1[22]+z[10]*tmp1[23]+z[2]*tmp1[24]+z[91]*tmp1[25]+z[83]*tmp1[26]+z[75]*tmp1[27]+z[67]*tmp1[28]+z[59]*tmp1[29]+z[51]*tmp1[30]
                    +z[43]*tmp1[31]+z[35]*tmp1[32]+z[27]*tmp1[33]+z[19]*tmp1[34]+z[11]*tmp1[35]+z[3]*tmp1[36]+z[92]*tmp1[37]+z[84]*tmp1[38]+z[76]*tmp1[39]+z[68]*tmp1[40]
                    +z[60]*tmp1[41]+z[52]*tmp1[42]+z[44]*tmp1[43]+z[36]*tmp1[44]+z[28]*tmp1[45]+z[20]*tmp1[46]+z[12]*tmp1[47]+z[4]*tmp1[48]+z[93]*tmp1[49]+z[85]*tmp1[50]
                    +z[77]*tmp1[51]+z[69]*tmp1[52]+z[61]*tmp1[53]+z[53]*tmp1[54]+z[45]*tmp1[55]+z[37]*tmp1[56]+z[29]*tmp1[57]+z[21]*tmp1[58]+z[13]*tmp1[59]+z[5]*tmp1[60]
                    +z[94]*tmp1[61]+z[86]*tmp1[62]+z[78]*tmp1[63]+z[70]*tmp1[64]+z[62]*tmp1[65]+z[54]*tmp1[66]+z[46]*tmp1[67]+z[38]*tmp1[68]+z[30]*tmp1[69]+z[22]*tmp1[70]
                    +z[14]*tmp1[71]+z[6]*tmp1[72]+z[95]*tmp1[73]+z[87]*tmp1[74]+z[79]*tmp1[75]+z[71]*tmp1[76]+z[63]*tmp1[77]+z[55]*tmp1[78]+z[47]*tmp1[79]+z[39]*tmp1[80]
                    +z[31]*tmp1[81]+z[23]*tmp1[82]+z[15]*tmp1[83]+z[7]*tmp1[84]+z[96]*tmp1[85]+z[88]*tmp1[86]+z[80]*tmp1[87]+z[72]*tmp1[88]+z[64]*tmp1[89]+z[56]*tmp1[90]
                    +z[48]*tmp1[91]+z[40]*tmp1[92]+z[32]*tmp1[93]+z[24]*tmp1[94]+z[16]*tmp1[95]+z[8]*tmp1[96];
                    tab[nb_tmp3+90*nb3]=tab[nb_tmp3+90*nb3]+z[0]*tmp1[0]
                    +z[90]*tmp1[1]+z[83]*tmp1[2]+z[76]*tmp1[3]+z[69]*tmp1[4]+z[62]*tmp1[5]+z[55]*tmp1[6]+z[48]*tmp1[7]+z[41]*tmp1[8]+z[34]*tmp1[9]+z[27]*tmp1[10]
                    +z[20]*tmp1[11]+z[13]*tmp1[12]+z[6]*tmp1[13]+z[96]*tmp1[14]+z[89]*tmp1[15]+z[82]*tmp1[16]+z[75]*tmp1[17]+z[68]*tmp1[18]+z[61]*tmp1[19]+z[54]*tmp1[20]
                    +z[47]*tmp1[21]+z[40]*tmp1[22]+z[33]*tmp1[23]+z[26]*tmp1[24]+z[19]*tmp1[25]+z[12]*tmp1[26]+z[5]*tmp1[27]+z[95]*tmp1[28]+z[88]*tmp1[29]+z[81]*tmp1[30]
                    +z[74]*tmp1[31]+z[67]*tmp1[32]+z[60]*tmp1[33]+z[53]*tmp1[34]+z[46]*tmp1[35]+z[39]*tmp1[36]+z[32]*tmp1[37]+z[25]*tmp1[38]+z[18]*tmp1[39]+z[11]*tmp1[40]
                    +z[4]*tmp1[41]+z[94]*tmp1[42]+z[87]*tmp1[43]+z[80]*tmp1[44]+z[73]*tmp1[45]+z[66]*tmp1[46]+z[59]*tmp1[47]+z[52]*tmp1[48]+z[45]*tmp1[49]+z[38]*tmp1[50]
                    +z[31]*tmp1[51]+z[24]*tmp1[52]+z[17]*tmp1[53]+z[10]*tmp1[54]+z[3]*tmp1[55]+z[93]*tmp1[56]+z[86]*tmp1[57]+z[79]*tmp1[58]+z[72]*tmp1[59]+z[65]*tmp1[60]
                    +z[58]*tmp1[61]+z[51]*tmp1[62]+z[44]*tmp1[63]+z[37]*tmp1[64]+z[30]*tmp1[65]+z[23]*tmp1[66]+z[16]*tmp1[67]+z[9]*tmp1[68]+z[2]*tmp1[69]+z[92]*tmp1[70]
                    +z[85]*tmp1[71]+z[78]*tmp1[72]+z[71]*tmp1[73]+z[64]*tmp1[74]+z[57]*tmp1[75]+z[50]*tmp1[76]+z[43]*tmp1[77]+z[36]*tmp1[78]+z[29]*tmp1[79]+z[22]*tmp1[80]
                    +z[15]*tmp1[81]+z[8]*tmp1[82]+z[1]*tmp1[83]+z[91]*tmp1[84]+z[84]*tmp1[85]+z[77]*tmp1[86]+z[70]*tmp1[87]+z[63]*tmp1[88]+z[56]*tmp1[89]+z[49]*tmp1[90]
                    +z[42]*tmp1[91]+z[35]*tmp1[92]+z[28]*tmp1[93]+z[21]*tmp1[94]+z[14]*tmp1[95]+z[7]*tmp1[96];
                    tab[nb_tmp3+91*nb3]=tab[nb_tmp3+91*nb3]+z[0]*tmp1[0]
                    +z[91]*tmp1[1]+z[85]*tmp1[2]+z[79]*tmp1[3]+z[73]*tmp1[4]+z[67]*tmp1[5]+z[61]*tmp1[6]+z[55]*tmp1[7]+z[49]*tmp1[8]+z[43]*tmp1[9]+z[37]*tmp1[10]
                    +z[31]*tmp1[11]+z[25]*tmp1[12]+z[19]*tmp1[13]+z[13]*tmp1[14]+z[7]*tmp1[15]+z[1]*tmp1[16]+z[92]*tmp1[17]+z[86]*tmp1[18]+z[80]*tmp1[19]+z[74]*tmp1[20]
                    +z[68]*tmp1[21]+z[62]*tmp1[22]+z[56]*tmp1[23]+z[50]*tmp1[24]+z[44]*tmp1[25]+z[38]*tmp1[26]+z[32]*tmp1[27]+z[26]*tmp1[28]+z[20]*tmp1[29]+z[14]*tmp1[30]
                    +z[8]*tmp1[31]+z[2]*tmp1[32]+z[93]*tmp1[33]+z[87]*tmp1[34]+z[81]*tmp1[35]+z[75]*tmp1[36]+z[69]*tmp1[37]+z[63]*tmp1[38]+z[57]*tmp1[39]+z[51]*tmp1[40]
                    +z[45]*tmp1[41]+z[39]*tmp1[42]+z[33]*tmp1[43]+z[27]*tmp1[44]+z[21]*tmp1[45]+z[15]*tmp1[46]+z[9]*tmp1[47]+z[3]*tmp1[48]+z[94]*tmp1[49]+z[88]*tmp1[50]
                    +z[82]*tmp1[51]+z[76]*tmp1[52]+z[70]*tmp1[53]+z[64]*tmp1[54]+z[58]*tmp1[55]+z[52]*tmp1[56]+z[46]*tmp1[57]+z[40]*tmp1[58]+z[34]*tmp1[59]+z[28]*tmp1[60]
                    +z[22]*tmp1[61]+z[16]*tmp1[62]+z[10]*tmp1[63]+z[4]*tmp1[64]+z[95]*tmp1[65]+z[89]*tmp1[66]+z[83]*tmp1[67]+z[77]*tmp1[68]+z[71]*tmp1[69]+z[65]*tmp1[70]
                    +z[59]*tmp1[71]+z[53]*tmp1[72]+z[47]*tmp1[73]+z[41]*tmp1[74]+z[35]*tmp1[75]+z[29]*tmp1[76]+z[23]*tmp1[77]+z[17]*tmp1[78]+z[11]*tmp1[79]+z[5]*tmp1[80]
                    +z[96]*tmp1[81]+z[90]*tmp1[82]+z[84]*tmp1[83]+z[78]*tmp1[84]+z[72]*tmp1[85]+z[66]*tmp1[86]+z[60]*tmp1[87]+z[54]*tmp1[88]+z[48]*tmp1[89]+z[42]*tmp1[90]
                    +z[36]*tmp1[91]+z[30]*tmp1[92]+z[24]*tmp1[93]+z[18]*tmp1[94]+z[12]*tmp1[95]+z[6]*tmp1[96];
                    tab[nb_tmp3+92*nb3]=tab[nb_tmp3+92*nb3]+z[0]*tmp1[0]
                    +z[92]*tmp1[1]+z[87]*tmp1[2]+z[82]*tmp1[3]+z[77]*tmp1[4]+z[72]*tmp1[5]+z[67]*tmp1[6]+z[62]*tmp1[7]+z[57]*tmp1[8]+z[52]*tmp1[9]+z[47]*tmp1[10]
                    +z[42]*tmp1[11]+z[37]*tmp1[12]+z[32]*tmp1[13]+z[27]*tmp1[14]+z[22]*tmp1[15]+z[17]*tmp1[16]+z[12]*tmp1[17]+z[7]*tmp1[18]+z[2]*tmp1[19]+z[94]*tmp1[20]
                    +z[89]*tmp1[21]+z[84]*tmp1[22]+z[79]*tmp1[23]+z[74]*tmp1[24]+z[69]*tmp1[25]+z[64]*tmp1[26]+z[59]*tmp1[27]+z[54]*tmp1[28]+z[49]*tmp1[29]+z[44]*tmp1[30]
                    +z[39]*tmp1[31]+z[34]*tmp1[32]+z[29]*tmp1[33]+z[24]*tmp1[34]+z[19]*tmp1[35]+z[14]*tmp1[36]+z[9]*tmp1[37]+z[4]*tmp1[38]+z[96]*tmp1[39]+z[91]*tmp1[40]
                    +z[86]*tmp1[41]+z[81]*tmp1[42]+z[76]*tmp1[43]+z[71]*tmp1[44]+z[66]*tmp1[45]+z[61]*tmp1[46]+z[56]*tmp1[47]+z[51]*tmp1[48]+z[46]*tmp1[49]+z[41]*tmp1[50]
                    +z[36]*tmp1[51]+z[31]*tmp1[52]+z[26]*tmp1[53]+z[21]*tmp1[54]+z[16]*tmp1[55]+z[11]*tmp1[56]+z[6]*tmp1[57]+z[1]*tmp1[58]+z[93]*tmp1[59]+z[88]*tmp1[60]
                    +z[83]*tmp1[61]+z[78]*tmp1[62]+z[73]*tmp1[63]+z[68]*tmp1[64]+z[63]*tmp1[65]+z[58]*tmp1[66]+z[53]*tmp1[67]+z[48]*tmp1[68]+z[43]*tmp1[69]+z[38]*tmp1[70]
                    +z[33]*tmp1[71]+z[28]*tmp1[72]+z[23]*tmp1[73]+z[18]*tmp1[74]+z[13]*tmp1[75]+z[8]*tmp1[76]+z[3]*tmp1[77]+z[95]*tmp1[78]+z[90]*tmp1[79]+z[85]*tmp1[80]
                    +z[80]*tmp1[81]+z[75]*tmp1[82]+z[70]*tmp1[83]+z[65]*tmp1[84]+z[60]*tmp1[85]+z[55]*tmp1[86]+z[50]*tmp1[87]+z[45]*tmp1[88]+z[40]*tmp1[89]+z[35]*tmp1[90]
                    +z[30]*tmp1[91]+z[25]*tmp1[92]+z[20]*tmp1[93]+z[15]*tmp1[94]+z[10]*tmp1[95]+z[5]*tmp1[96];
                    tab[nb_tmp3+93*nb3]=tab[nb_tmp3+93*nb3]+z[0]*tmp1[0]
                    +z[93]*tmp1[1]+z[89]*tmp1[2]+z[85]*tmp1[3]+z[81]*tmp1[4]+z[77]*tmp1[5]+z[73]*tmp1[6]+z[69]*tmp1[7]+z[65]*tmp1[8]+z[61]*tmp1[9]+z[57]*tmp1[10]
                    +z[53]*tmp1[11]+z[49]*tmp1[12]+z[45]*tmp1[13]+z[41]*tmp1[14]+z[37]*tmp1[15]+z[33]*tmp1[16]+z[29]*tmp1[17]+z[25]*tmp1[18]+z[21]*tmp1[19]+z[17]*tmp1[20]
                    +z[13]*tmp1[21]+z[9]*tmp1[22]+z[5]*tmp1[23]+z[1]*tmp1[24]+z[94]*tmp1[25]+z[90]*tmp1[26]+z[86]*tmp1[27]+z[82]*tmp1[28]+z[78]*tmp1[29]+z[74]*tmp1[30]
                    +z[70]*tmp1[31]+z[66]*tmp1[32]+z[62]*tmp1[33]+z[58]*tmp1[34]+z[54]*tmp1[35]+z[50]*tmp1[36]+z[46]*tmp1[37]+z[42]*tmp1[38]+z[38]*tmp1[39]+z[34]*tmp1[40]
                    +z[30]*tmp1[41]+z[26]*tmp1[42]+z[22]*tmp1[43]+z[18]*tmp1[44]+z[14]*tmp1[45]+z[10]*tmp1[46]+z[6]*tmp1[47]+z[2]*tmp1[48]+z[95]*tmp1[49]+z[91]*tmp1[50]
                    +z[87]*tmp1[51]+z[83]*tmp1[52]+z[79]*tmp1[53]+z[75]*tmp1[54]+z[71]*tmp1[55]+z[67]*tmp1[56]+z[63]*tmp1[57]+z[59]*tmp1[58]+z[55]*tmp1[59]+z[51]*tmp1[60]
                    +z[47]*tmp1[61]+z[43]*tmp1[62]+z[39]*tmp1[63]+z[35]*tmp1[64]+z[31]*tmp1[65]+z[27]*tmp1[66]+z[23]*tmp1[67]+z[19]*tmp1[68]+z[15]*tmp1[69]+z[11]*tmp1[70]
                    +z[7]*tmp1[71]+z[3]*tmp1[72]+z[96]*tmp1[73]+z[92]*tmp1[74]+z[88]*tmp1[75]+z[84]*tmp1[76]+z[80]*tmp1[77]+z[76]*tmp1[78]+z[72]*tmp1[79]+z[68]*tmp1[80]
                    +z[64]*tmp1[81]+z[60]*tmp1[82]+z[56]*tmp1[83]+z[52]*tmp1[84]+z[48]*tmp1[85]+z[44]*tmp1[86]+z[40]*tmp1[87]+z[36]*tmp1[88]+z[32]*tmp1[89]+z[28]*tmp1[90]
                    +z[24]*tmp1[91]+z[20]*tmp1[92]+z[16]*tmp1[93]+z[12]*tmp1[94]+z[8]*tmp1[95]+z[4]*tmp1[96];
                    tab[nb_tmp3+94*nb3]=tab[nb_tmp3+94*nb3]+z[0]*tmp1[0]
                    +z[94]*tmp1[1]+z[91]*tmp1[2]+z[88]*tmp1[3]+z[85]*tmp1[4]+z[82]*tmp1[5]+z[79]*tmp1[6]+z[76]*tmp1[7]+z[73]*tmp1[8]+z[70]*tmp1[9]+z[67]*tmp1[10]
                    +z[64]*tmp1[11]+z[61]*tmp1[12]+z[58]*tmp1[13]+z[55]*tmp1[14]+z[52]*tmp1[15]+z[49]*tmp1[16]+z[46]*tmp1[17]+z[43]*tmp1[18]+z[40]*tmp1[19]+z[37]*tmp1[20]
                    +z[34]*tmp1[21]+z[31]*tmp1[22]+z[28]*tmp1[23]+z[25]*tmp1[24]+z[22]*tmp1[25]+z[19]*tmp1[26]+z[16]*tmp1[27]+z[13]*tmp1[28]+z[10]*tmp1[29]+z[7]*tmp1[30]
                    +z[4]*tmp1[31]+z[1]*tmp1[32]+z[95]*tmp1[33]+z[92]*tmp1[34]+z[89]*tmp1[35]+z[86]*tmp1[36]+z[83]*tmp1[37]+z[80]*tmp1[38]+z[77]*tmp1[39]+z[74]*tmp1[40]
                    +z[71]*tmp1[41]+z[68]*tmp1[42]+z[65]*tmp1[43]+z[62]*tmp1[44]+z[59]*tmp1[45]+z[56]*tmp1[46]+z[53]*tmp1[47]+z[50]*tmp1[48]+z[47]*tmp1[49]+z[44]*tmp1[50]
                    +z[41]*tmp1[51]+z[38]*tmp1[52]+z[35]*tmp1[53]+z[32]*tmp1[54]+z[29]*tmp1[55]+z[26]*tmp1[56]+z[23]*tmp1[57]+z[20]*tmp1[58]+z[17]*tmp1[59]+z[14]*tmp1[60]
                    +z[11]*tmp1[61]+z[8]*tmp1[62]+z[5]*tmp1[63]+z[2]*tmp1[64]+z[96]*tmp1[65]+z[93]*tmp1[66]+z[90]*tmp1[67]+z[87]*tmp1[68]+z[84]*tmp1[69]+z[81]*tmp1[70]
                    +z[78]*tmp1[71]+z[75]*tmp1[72]+z[72]*tmp1[73]+z[69]*tmp1[74]+z[66]*tmp1[75]+z[63]*tmp1[76]+z[60]*tmp1[77]+z[57]*tmp1[78]+z[54]*tmp1[79]+z[51]*tmp1[80]
                    +z[48]*tmp1[81]+z[45]*tmp1[82]+z[42]*tmp1[83]+z[39]*tmp1[84]+z[36]*tmp1[85]+z[33]*tmp1[86]+z[30]*tmp1[87]+z[27]*tmp1[88]+z[24]*tmp1[89]+z[21]*tmp1[90]
                    +z[18]*tmp1[91]+z[15]*tmp1[92]+z[12]*tmp1[93]+z[9]*tmp1[94]+z[6]*tmp1[95]+z[3]*tmp1[96];
                    tab[nb_tmp3+95*nb3]=tab[nb_tmp3+95*nb3]+z[0]*tmp1[0]
                    +z[95]*tmp1[1]+z[93]*tmp1[2]+z[91]*tmp1[3]+z[89]*tmp1[4]+z[87]*tmp1[5]+z[85]*tmp1[6]+z[83]*tmp1[7]+z[81]*tmp1[8]+z[79]*tmp1[9]+z[77]*tmp1[10]
                    +z[75]*tmp1[11]+z[73]*tmp1[12]+z[71]*tmp1[13]+z[69]*tmp1[14]+z[67]*tmp1[15]+z[65]*tmp1[16]+z[63]*tmp1[17]+z[61]*tmp1[18]+z[59]*tmp1[19]+z[57]*tmp1[20]
                    +z[55]*tmp1[21]+z[53]*tmp1[22]+z[51]*tmp1[23]+z[49]*tmp1[24]+z[47]*tmp1[25]+z[45]*tmp1[26]+z[43]*tmp1[27]+z[41]*tmp1[28]+z[39]*tmp1[29]+z[37]*tmp1[30]
                    +z[35]*tmp1[31]+z[33]*tmp1[32]+z[31]*tmp1[33]+z[29]*tmp1[34]+z[27]*tmp1[35]+z[25]*tmp1[36]+z[23]*tmp1[37]+z[21]*tmp1[38]+z[19]*tmp1[39]+z[17]*tmp1[40]
                    +z[15]*tmp1[41]+z[13]*tmp1[42]+z[11]*tmp1[43]+z[9]*tmp1[44]+z[7]*tmp1[45]+z[5]*tmp1[46]+z[3]*tmp1[47]+z[1]*tmp1[48]+z[96]*tmp1[49]+z[94]*tmp1[50]
                    +z[92]*tmp1[51]+z[90]*tmp1[52]+z[88]*tmp1[53]+z[86]*tmp1[54]+z[84]*tmp1[55]+z[82]*tmp1[56]+z[80]*tmp1[57]+z[78]*tmp1[58]+z[76]*tmp1[59]+z[74]*tmp1[60]
                    +z[72]*tmp1[61]+z[70]*tmp1[62]+z[68]*tmp1[63]+z[66]*tmp1[64]+z[64]*tmp1[65]+z[62]*tmp1[66]+z[60]*tmp1[67]+z[58]*tmp1[68]+z[56]*tmp1[69]+z[54]*tmp1[70]
                    +z[52]*tmp1[71]+z[50]*tmp1[72]+z[48]*tmp1[73]+z[46]*tmp1[74]+z[44]*tmp1[75]+z[42]*tmp1[76]+z[40]*tmp1[77]+z[38]*tmp1[78]+z[36]*tmp1[79]+z[34]*tmp1[80]
                    +z[32]*tmp1[81]+z[30]*tmp1[82]+z[28]*tmp1[83]+z[26]*tmp1[84]+z[24]*tmp1[85]+z[22]*tmp1[86]+z[20]*tmp1[87]+z[18]*tmp1[88]+z[16]*tmp1[89]+z[14]*tmp1[90]
                    +z[12]*tmp1[91]+z[10]*tmp1[92]+z[8]*tmp1[93]+z[6]*tmp1[94]+z[4]*tmp1[95]+z[2]*tmp1[96];
                    tab[nb_tmp3+96*nb3]=tab[nb_tmp3+96*nb3]+z[0]*tmp1[0]
                    +z[96]*tmp1[1]+z[95]*tmp1[2]+z[94]*tmp1[3]+z[93]*tmp1[4]+z[92]*tmp1[5]+z[91]*tmp1[6]+z[90]*tmp1[7]+z[89]*tmp1[8]+z[88]*tmp1[9]+z[87]*tmp1[10]
                    +z[86]*tmp1[11]+z[85]*tmp1[12]+z[84]*tmp1[13]+z[83]*tmp1[14]+z[82]*tmp1[15]+z[81]*tmp1[16]+z[80]*tmp1[17]+z[79]*tmp1[18]+z[78]*tmp1[19]+z[77]*tmp1[20]
                    +z[76]*tmp1[21]+z[75]*tmp1[22]+z[74]*tmp1[23]+z[73]*tmp1[24]+z[72]*tmp1[25]+z[71]*tmp1[26]+z[70]*tmp1[27]+z[69]*tmp1[28]+z[68]*tmp1[29]+z[67]*tmp1[30]
                    +z[66]*tmp1[31]+z[65]*tmp1[32]+z[64]*tmp1[33]+z[63]*tmp1[34]+z[62]*tmp1[35]+z[61]*tmp1[36]+z[60]*tmp1[37]+z[59]*tmp1[38]+z[58]*tmp1[39]+z[57]*tmp1[40]
                    +z[56]*tmp1[41]+z[55]*tmp1[42]+z[54]*tmp1[43]+z[53]*tmp1[44]+z[52]*tmp1[45]+z[51]*tmp1[46]+z[50]*tmp1[47]+z[49]*tmp1[48]+z[48]*tmp1[49]+z[47]*tmp1[50]
                    +z[46]*tmp1[51]+z[45]*tmp1[52]+z[44]*tmp1[53]+z[43]*tmp1[54]+z[42]*tmp1[55]+z[41]*tmp1[56]+z[40]*tmp1[57]+z[39]*tmp1[58]+z[38]*tmp1[59]+z[37]*tmp1[60]
                    +z[36]*tmp1[61]+z[35]*tmp1[62]+z[34]*tmp1[63]+z[33]*tmp1[64]+z[32]*tmp1[65]+z[31]*tmp1[66]+z[30]*tmp1[67]+z[29]*tmp1[68]+z[28]*tmp1[69]+z[27]*tmp1[70]
                    +z[26]*tmp1[71]+z[25]*tmp1[72]+z[24]*tmp1[73]+z[23]*tmp1[74]+z[22]*tmp1[75]+z[21]*tmp1[76]+z[20]*tmp1[77]+z[19]*tmp1[78]+z[18]*tmp1[79]+z[17]*tmp1[80]
                    +z[16]*tmp1[81]+z[15]*tmp1[82]+z[14]*tmp1[83]+z[13]*tmp1[84]+z[12]*tmp1[85]+z[11]*tmp1[86]+z[10]*tmp1[87]+z[9]*tmp1[88]+z[8]*tmp1[89]+z[7]*tmp1[90]
                    +z[6]*tmp1[91]+z[5]*tmp1[92]+z[4]*tmp1[93]+z[3]*tmp1[94]+z[2]*tmp1[95]+z[1]*tmp1[96];



                }

            }
        }
    }
///////////////////////////////////////////////////////////////////////


    void fun_inverse_table_FFT(int M,std::complex<double> tab[])
    {
        int rx5=5,rx4=4,rx3=3,rx2=2,rx7=7,rx11=11;
        int stg[100]={};
        int *tab8 = new int[M];
        int *tab9 = new int[M];
        std::complex<double> *tab11 = new std::complex<double>[M];
        int nb_stages=0;
        int nb1=0;
        int nb2=0;
        int nb3=0;

        nb_stages=6;

	nb_stages=radix_base(M,stg);

        for(int i=0;i<M;i++)
        {
            //tab9[i]=tab2[i];
            //tab8[i]=tab2[i];
            tab9[i]=i;
            tab8[i]=i;
        }

        nb3=1;
        for(int op=nb_stages;op>=2;op--)
        {
            nb1=stg[op];
            nb3=nb3*stg[op];

            if(op==nb_stages)
            {
                nb2=stg[0];
            }
            else
            {
                nb2=nb2*stg[op+1];
            }

               for(int i=0,n=0,p=0;i<M;i=i+M/nb3,n++)
            {
                if(n>=nb1)
                {
                    n=0,p=p+M/nb2;
                }
                for(int j=0,k=0;j<M/nb3;j++,k=k+nb1)
                {
                    if(op%2==0)
                    {
                        tab8[i+j]=tab9[k+n+p];
                    }
                    else
                    {
                        tab9[i+j]=tab8[k+n+p];
                    }
                }
            }
        }

        for(int i=0;i<M;i++)
        {
          tab11[i]=tab[tab8[i]];

        }
        for(int i=0;i<M;i++)
        {
          tab[i]=tab11[i];
        }

        delete [] tab8;
        delete [] tab9;
        delete [] tab11;
    }
///////////////////////////////////////

     int radix_base(int N,int stg[])
        {
        int k=0;
        double M=(double)N;
        double epsilon1;
        stg[0]=1;
        //*flg2=0;
        //cout<<"M= "<<M<<endl;
        epsilon1=0.000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000000001;

        for(int j=0;j<200;j++)
        {
            if(fmod(M,97)<=epsilon1)
            {
                k++;
                M=M/97.0;
                stg[k]=97;
            }
            else if(fmod(M,11)<=epsilon1)
            {
                k++;
                M=M/11.0;
                stg[k]=11;
            }
            else if(fmod(M,7)<=epsilon1)
            {
                k++;
                M=M/7.0;
                stg[k]=7;
            }
            else if(fmod(M,5)<=epsilon1)
            {
                k++;
                M=M/5.0;
                stg[k]=5;
            }
            else if(fmod(M,4)<=epsilon1)
            {
                k++;
                M=M/4.0;
                stg[k]=4;
            }
            else if(fmod(M,3)<=epsilon1)
            {
                k++;
                M=M/3.0;
                stg[k]=3;
            }
            else if(fmod(M,2)<=epsilon1)
            {
                k++;
                M=M/2.0;
                stg[k]=2;
            }
            else if(M>=1.0-epsilon1&&M<=1.0+epsilon1)
            {
                //*flag=*flag+1;
                 //cout<<"*flag= "<<*flag<<" N= "<<N<<endl;
                break;
            }
            else
            {

                cout<<endl<< "Unsupported signal: N= "<<N<< endl;
                if(k>0)
                {
                    for(int m=1;m<=k;m++)
                    {
                     cout <<"stage:"<<m<<" = radix-"<<stg[m] << endl;
                    }
                }
                cout <<"stage:"<<k+1<<" = radix-??" << endl;

                k=0;
                break;
            }
            //*flg2=*flg2+1;
        }
       return k;
    }
/////////////////////////////////////////////////////////////


	void fun_fourier_transform_FFT_radix_2_stg_first(std::complex<double> tab[],int stg_first,double fi2,std::complex<double> z[])
    {
        std::complex<double> tmp1,tmp2,tmp40,tmp50;
        std::complex<double>  w[2]={{1,0}};

        w[0].real(cos(0+fi2));
        w[0].imag(-sin(0+fi2));
        w[1].real(cos(0+fi2));
        w[1].imag(-sin(0+fi2));

        for(int i=0;i<stg_first;i=i+1)
        {
            tmp1=w[0]*tab[i+0*stg_first];
            tmp2=w[1]*tab[i+1*stg_first];

            tmp40=z[0]*(tmp1+tmp2);

            tab[i+0*stg_first]=tmp40;
            tab[i+1*stg_first]=z[0]*tmp1+z[1]*tmp2;
        }
    }
///////////////////////////////////////////////////////////////


   void fun_fourier_transform_FFT_radix_3_stg_first(std::complex<double> tab[],int stg_first,double fi2,std::complex<double> z[])
    {
        std::complex<double> tmp1,tmp2,tmp3,tmp4,tmp5;
        std::complex<double>  w[3]={{1,0}};

        w[0].real(cos(0+fi2));
        w[0].imag(-sin(0+fi2));
        w[1].real(cos(0+fi2));
        w[1].imag(-sin(0+fi2));
        w[2].real(cos(0+fi2));
        w[2].imag(-sin(0+fi2));

        for(int i=0;i<stg_first;i=i+1)
        {
          tmp1=w[0]*tab[i+0*stg_first];
          tmp2=w[1]*tab[i+1*stg_first];
          tmp3=w[2]*tab[i+2*stg_first];

          tmp4=z[0]*tmp1;
          tmp5=z[0]*(tmp2+tmp3);

         //radix-3
          tab[i+0*stg_first]  = tmp4+tmp5;
          tab[i+1*stg_first]  = tmp4+z[1]*tmp2+z[2]*tmp3;
          tab[i+2*stg_first]  = tmp4+z[2]*tmp2+z[1]*tmp3;
        }
    }
///////////////////////////////////////////////////////////////


   void fun_fourier_transform_FFT_radix_4_stg_first(std::complex<double> tab[],int stg_first,double fi2,std::complex<double> z[])
    {
        std::complex<double> tmp1,tmp2,tmp3,tmp4;
        std::complex<double> tmp101,tmp102,tmp103,tmp104;
        std::complex<double> tmp111,tmp112,tmp113,tmp114;
        std::complex<double>  w[4]={{1,0}};

        w[0].real(cos(0+fi2));
        w[0].imag(-sin(0+fi2));
        w[1].real(cos(0+fi2));
        w[1].imag(-sin(0+fi2));
        w[2].real(cos(0+fi2));
        w[2].imag(-sin(0+fi2));
        w[3].real(cos(0+fi2));
        w[3].imag(-sin(0+fi2));
        for(int i=0;i<stg_first;i=i+1)
        {
            tmp1=w[0]*tab[i+0*stg_first];
            tmp2=w[1]*tab[i+1*stg_first];
            tmp3=w[2]*tab[i+2*stg_first];
            tmp4=w[3]*tab[i+3*stg_first];


            tmp101=tmp2-tmp4;
            tmp102=tmp2+tmp4;
            tmp103=tmp1-tmp3;
            tmp104=tmp1+tmp3;

            tmp111=z[0]*(tmp104+tmp102);
            tmp112=z[0]*tmp103;
            tmp113=z[0]*(tmp104-tmp102);
            tmp114=z[0]*tmp103;

            //radix-4
            tab[i+0*stg_first]   =tmp111;
            tab[i+1*stg_first]   =tmp112+z[1]*tmp101;
            tab[i+2*stg_first]   =tmp113;
            tab[i+3*stg_first]   =tmp114-z[1]*tmp101;
        }
    }
/////////////////////////////////////////////////////////////////


void fun_fourier_transform_FFT_radix_5_stg_first(std::complex<double> tab[],int stg_first,double fi2,std::complex<double> z_rx5[])
  {
        std::complex<double> tmp1,tmp2,tmp3,tmp4,tmp5,tmp33,tmp34;
        std::complex<double>  w[5]={{1,0}};

        w[0].real(cos(0+fi2));
        w[0].imag(-sin(0+fi2));
        w[1].real(cos(0+fi2));
        w[1].imag(-sin(0+fi2));
        w[2].real(cos(0+fi2));
        w[2].imag(-sin(0+fi2));
        w[3].real(cos(0+fi2));
        w[3].imag(-sin(0+fi2));
        w[4].real(cos(0+fi2));
        w[4].imag(-sin(0+fi2));

        for(int i=0;i<stg_first;i++)
        {
          tmp1=w[0]*tab[i+0*stg_first];
          tmp2=w[1]*tab[i+1*stg_first];
          tmp3=w[2]*tab[i+2*stg_first];
          tmp4=w[4]*tab[i+3*stg_first];
          tmp5=w[4]*tab[i+4*stg_first];

          tmp33=z_rx5[0]*tmp1;
          tmp34=z_rx5[0]*(tmp1+tmp2+tmp3+tmp4+tmp5);
         //radix-5
          tab[i+0*stg_first] =tmp34;
          tab[i+1*stg_first] =tmp33+z_rx5[1]*tmp2+z_rx5[2]*tmp3+z_rx5[3]*tmp4+z_rx5[4]*tmp5;
          tab[i+2*stg_first] =tmp33+z_rx5[2]*tmp2+z_rx5[4]*tmp3+z_rx5[1]*tmp4+z_rx5[3]*tmp5;
          tab[i+3*stg_first] =tmp33+z_rx5[3]*tmp2+z_rx5[1]*tmp3+z_rx5[4]*tmp4+z_rx5[2]*tmp5;
          tab[i+4*stg_first] =tmp33+z_rx5[4]*tmp2+z_rx5[3]*tmp3+z_rx5[2]*tmp4+z_rx5[1]*tmp5;
        }
  }
/////////////////////////////////////////


   void fun_fourier_transform_FFT_radix_7_stg_first(std::complex<double> tab[],int stg_first,double fi2,std::complex<double> z[])
  {
        std::complex<double> tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp10,tmp20;
        std::complex<double>  w[7]={{1,0}};


        w[0].real(cos(0+fi2));
        w[0].imag(-sin(0+fi2));
        w[1].real(cos(0+fi2));
        w[1].imag(-sin(0+fi2));
        w[2].real(cos(0+fi2));
        w[2].imag(-sin(0+fi2));
        w[3].real(cos(0+fi2));
        w[3].imag(-sin(0+fi2));
        w[4].real(cos(0+fi2));
        w[4].imag(-sin(0+fi2));
        w[5].real(cos(0+fi2));
        w[5].imag(-sin(0+fi2));
        w[6].real(cos(0+fi2));
        w[6].imag(-sin(0+fi2));

        for(int i=0;i<stg_first;i=i+1)
        {
          tmp1=w[0]*tab[i+0*stg_first];
          tmp2=w[1]*tab[i+1*stg_first];
          tmp3=w[2]*tab[i+2*stg_first];
          tmp4=w[3]*tab[i+3*stg_first];
          tmp5=w[4]*tab[i+4*stg_first];
          tmp6=w[5]*tab[i+5*stg_first];
          tmp7=w[6]*tab[i+6*stg_first];

          tmp10=z[0]*tmp1;
          tmp20=z[0]*(tmp1+tmp2+tmp3+tmp4+tmp5+tmp6+tmp7);

         //radix-7
          tab[i+0*stg_first] =tmp20;
          tab[i+1*stg_first] =tmp10+z[1]*tmp2+z[2]*tmp3+z[3]*tmp4+z[4]*tmp5+z[5]*tmp6+z[6]*tmp7;
          tab[i+2*stg_first] =tmp10+z[2]*tmp2+z[4]*tmp3+z[6]*tmp4+z[1]*tmp5+z[3]*tmp6+z[5]*tmp7;
          tab[i+3*stg_first] =tmp10+z[3]*tmp2+z[6]*tmp3+z[2]*tmp4+z[5]*tmp5+z[1]*tmp6+z[4]*tmp7;
          tab[i+4*stg_first] =tmp10+z[4]*tmp2+z[1]*tmp3+z[5]*tmp4+z[2]*tmp5+z[6]*tmp6+z[3]*tmp7;
          tab[i+5*stg_first] =tmp10+z[5]*tmp2+z[3]*tmp3+z[1]*tmp4+z[6]*tmp5+z[4]*tmp6+z[2]*tmp7;
          tab[i+6*stg_first] =tmp10+z[6]*tmp2+z[5]*tmp3+z[4]*tmp4+z[3]*tmp5+z[2]*tmp6+z[1]*tmp7;
        }

    }
///////////////////////////////////////////////////////////////////////


 void fun_fourier_transform_FFT_radix_11_stg_first(std::complex<double> tab[],int stg_first,double fi2,std::complex<double> z[])
  {
        std::complex<double> tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10,tmp11,tmp33,tmp34;
        std::complex<double>  w[11]={{1,0}};

        w[0].real(cos(0+fi2));
        w[0].imag(-sin(0+fi2));
        w[1].real(cos(0+fi2));
        w[1].imag(-sin(0+fi2));
        w[2].real(cos(0+fi2));
        w[2].imag(-sin(0+fi2));
        w[3].real(cos(0+fi2));
        w[3].imag(-sin(0+fi2));
        w[4].real(cos(0+fi2));
        w[4].imag(-sin(0+fi2));
        w[5].real(cos(0+fi2));
        w[5].imag(-sin(0+fi2));
        w[6].real(cos(0+fi2));
        w[6].imag(-sin(0+fi2));
        w[7].real(cos(0+fi2));
        w[7].imag(-sin(0+fi2));
        w[8].real(cos(0+fi2));
        w[8].imag(-sin(0+fi2));
        w[9].real(cos(0+fi2));
        w[9].imag(-sin(0+fi2));
        w[10].real(cos(0+fi2));
        w[10].imag(-sin(0+fi2));


        for(int i=0;i<stg_first;i=i+1)
        {
          tmp1=w[0]*tab[i+0*stg_first];
          tmp2=w[1]*tab[i+1*stg_first];
          tmp3=w[2]*tab[i+2*stg_first];
          tmp4=w[3]*tab[i+3*stg_first];
          tmp5=w[4]*tab[i+4*stg_first];
          tmp6=w[5]*tab[i+5*stg_first];
          tmp7=w[6]*tab[i+6*stg_first];
          tmp8=w[7]*tab[i+7*stg_first];
          tmp9=w[8]*tab[i+8*stg_first];
          tmp10=w[9]*tab[i+9*stg_first];
          tmp11=w[10]*tab[i+10*stg_first];

          tmp33=z[0]*tmp1;
          tmp34=z[0]*(tmp1+tmp2+tmp3+tmp4+tmp5+tmp6+tmp7+tmp8+tmp9+tmp10+tmp11);

         //radix-11
          tab[i+0*stg_first] =tmp34;
          tab[i+1*stg_first] =tmp33+z[1]*tmp2+z[2]*tmp3+z[3]*tmp4+z[4]*tmp5+z[5]*tmp6+z[6]*tmp7+z[7]*tmp8+z[8]*tmp9+z[9]*tmp10+z[10]*tmp11;
          tab[i+2*stg_first] =tmp33+z[2]*tmp2+z[4]*tmp3+z[6]*tmp4+z[8]*tmp5+z[10]*tmp6+z[1]*tmp7+z[3]*tmp8+z[5]*tmp9+z[7]*tmp10+z[9]*tmp11;
          tab[i+3*stg_first] =tmp33+z[3]*tmp2+z[6]*tmp3+z[9]*tmp4+z[1]*tmp5+z[4]*tmp6+z[7]*tmp7+z[10]*tmp8+z[2]*tmp9+z[5]*tmp10+z[8]*tmp11;
          tab[i+4*stg_first] =tmp33+z[4]*tmp2+z[8]*tmp3+z[1]*tmp4+z[5]*tmp5+z[9]*tmp6+z[2]*tmp7+z[6]*tmp8+z[10]*tmp9+z[3]*tmp10+z[7]*tmp11;
          tab[i+5*stg_first] =tmp33+z[5]*tmp2+z[10]*tmp3+z[4]*tmp4+z[9]*tmp5+z[3]*tmp6+z[8]*tmp7+z[2]*tmp8+z[7]*tmp9+z[1]*tmp10+z[6]*tmp11;
          tab[i+6*stg_first] =tmp33+z[6]*tmp2+z[1]*tmp3+z[7]*tmp4+z[2]*tmp5+z[8]*tmp6+z[3]*tmp7+z[9]*tmp8+z[4]*tmp9+z[10]*tmp10+z[5]*tmp11;
          tab[i+7*stg_first] =tmp33+z[7]*tmp2+z[3]*tmp3+z[10]*tmp4+z[6]*tmp5+z[2]*tmp6+z[9]*tmp7+z[5]*tmp8+z[1]*tmp9+z[8]*tmp10+z[4]*tmp11;
          tab[i+8*stg_first] =tmp33+z[8]*tmp2+z[5]*tmp3+z[2]*tmp4+z[10]*tmp5+z[7]*tmp6+z[4]*tmp7+z[1]*tmp8+z[9]*tmp9+z[6]*tmp10+z[3]*tmp11;
          tab[i+9*stg_first] =tmp33+z[9]*tmp2+z[7]*tmp3+z[5]*tmp4+z[3]*tmp5+z[1]*tmp6+z[10]*tmp7+z[8]*tmp8+z[6]*tmp9+z[4]*tmp10+z[2]*tmp11;
          tab[i+10*stg_first] =tmp33+z[10]*tmp2+z[9]*tmp3+z[8]*tmp4+z[7]*tmp5+z[6]*tmp6+z[5]*tmp7+z[4]*tmp8+z[3]*tmp9+z[2]*tmp10+z[1]*tmp11;
        }
    }
///////////////////////////////////////////////////

 void fun_fourier_transform_FFT_radix_97_stg_first(std::complex<double> tab[],int stg_first,double fi2,std::complex<double> z[])
  {
        std::complex<double> tmp1[97];
        std::complex<double>  w[97]={{1,0}};


	for(int j9=0;j9<97;j9++)
	{
        w[j9].real(cos(0+fi2));
        w[j9].imag(-sin(0+fi2));
	}

        for(int i=0;i<stg_first;i=i+1)
        {
          for(int j9=0;j9<97;j9++)
          {
                tmp1[j9]=w[j9]*tab[i+j9*stg_first];
          }

          //radix-97
          for(int i9=0;i9<97;i9++)
          {
            for(int j9=0;j9<97;j9++)
            {
                tab[i+i9*stg_first].real(0);
                tab[i+i9*stg_first].imag(0);
            }
          }

            tab[i+0*stg_first]=tab[i+0*stg_first]+z[0]*tmp1[0]
            +z[0]*tmp1[1]+z[0]*tmp1[2]+z[0]*tmp1[3]+z[0]*tmp1[4]+z[0]*tmp1[5]+z[0]*tmp1[6]+z[0]*tmp1[7]+z[0]*tmp1[8]+z[0]*tmp1[9]+z[0]*tmp1[10]
            +z[0]*tmp1[11]+z[0]*tmp1[12]+z[0]*tmp1[13]+z[0]*tmp1[14]+z[0]*tmp1[15]+z[0]*tmp1[16]+z[0]*tmp1[17]+z[0]*tmp1[18]+z[0]*tmp1[19]+z[0]*tmp1[20]
            +z[0]*tmp1[21]+z[0]*tmp1[22]+z[0]*tmp1[23]+z[0]*tmp1[24]+z[0]*tmp1[25]+z[0]*tmp1[26]+z[0]*tmp1[27]+z[0]*tmp1[28]+z[0]*tmp1[29]+z[0]*tmp1[30]
            +z[0]*tmp1[31]+z[0]*tmp1[32]+z[0]*tmp1[33]+z[0]*tmp1[34]+z[0]*tmp1[35]+z[0]*tmp1[36]+z[0]*tmp1[37]+z[0]*tmp1[38]+z[0]*tmp1[39]+z[0]*tmp1[40]
            +z[0]*tmp1[41]+z[0]*tmp1[42]+z[0]*tmp1[43]+z[0]*tmp1[44]+z[0]*tmp1[45]+z[0]*tmp1[46]+z[0]*tmp1[47]+z[0]*tmp1[48]+z[0]*tmp1[49]+z[0]*tmp1[50]
            +z[0]*tmp1[51]+z[0]*tmp1[52]+z[0]*tmp1[53]+z[0]*tmp1[54]+z[0]*tmp1[55]+z[0]*tmp1[56]+z[0]*tmp1[57]+z[0]*tmp1[58]+z[0]*tmp1[59]+z[0]*tmp1[60]
            +z[0]*tmp1[61]+z[0]*tmp1[62]+z[0]*tmp1[63]+z[0]*tmp1[64]+z[0]*tmp1[65]+z[0]*tmp1[66]+z[0]*tmp1[67]+z[0]*tmp1[68]+z[0]*tmp1[69]+z[0]*tmp1[70]
            +z[0]*tmp1[71]+z[0]*tmp1[72]+z[0]*tmp1[73]+z[0]*tmp1[74]+z[0]*tmp1[75]+z[0]*tmp1[76]+z[0]*tmp1[77]+z[0]*tmp1[78]+z[0]*tmp1[79]+z[0]*tmp1[80]
            +z[0]*tmp1[81]+z[0]*tmp1[82]+z[0]*tmp1[83]+z[0]*tmp1[84]+z[0]*tmp1[85]+z[0]*tmp1[86]+z[0]*tmp1[87]+z[0]*tmp1[88]+z[0]*tmp1[89]+z[0]*tmp1[90]
            +z[0]*tmp1[91]+z[0]*tmp1[92]+z[0]*tmp1[93]+z[0]*tmp1[94]+z[0]*tmp1[95]+z[0]*tmp1[96];
            tab[i+1*stg_first]=tab[i+1*stg_first]+z[0]*tmp1[0]
            +z[1]*tmp1[1]+z[2]*tmp1[2]+z[3]*tmp1[3]+z[4]*tmp1[4]+z[5]*tmp1[5]+z[6]*tmp1[6]+z[7]*tmp1[7]+z[8]*tmp1[8]+z[9]*tmp1[9]+z[10]*tmp1[10]
            +z[11]*tmp1[11]+z[12]*tmp1[12]+z[13]*tmp1[13]+z[14]*tmp1[14]+z[15]*tmp1[15]+z[16]*tmp1[16]+z[17]*tmp1[17]+z[18]*tmp1[18]+z[19]*tmp1[19]+z[20]*tmp1[20]
            +z[21]*tmp1[21]+z[22]*tmp1[22]+z[23]*tmp1[23]+z[24]*tmp1[24]+z[25]*tmp1[25]+z[26]*tmp1[26]+z[27]*tmp1[27]+z[28]*tmp1[28]+z[29]*tmp1[29]+z[30]*tmp1[30]
            +z[31]*tmp1[31]+z[32]*tmp1[32]+z[33]*tmp1[33]+z[34]*tmp1[34]+z[35]*tmp1[35]+z[36]*tmp1[36]+z[37]*tmp1[37]+z[38]*tmp1[38]+z[39]*tmp1[39]+z[40]*tmp1[40]
            +z[41]*tmp1[41]+z[42]*tmp1[42]+z[43]*tmp1[43]+z[44]*tmp1[44]+z[45]*tmp1[45]+z[46]*tmp1[46]+z[47]*tmp1[47]+z[48]*tmp1[48]+z[49]*tmp1[49]+z[50]*tmp1[50]
            +z[51]*tmp1[51]+z[52]*tmp1[52]+z[53]*tmp1[53]+z[54]*tmp1[54]+z[55]*tmp1[55]+z[56]*tmp1[56]+z[57]*tmp1[57]+z[58]*tmp1[58]+z[59]*tmp1[59]+z[60]*tmp1[60]
            +z[61]*tmp1[61]+z[62]*tmp1[62]+z[63]*tmp1[63]+z[64]*tmp1[64]+z[65]*tmp1[65]+z[66]*tmp1[66]+z[67]*tmp1[67]+z[68]*tmp1[68]+z[69]*tmp1[69]+z[70]*tmp1[70]
            +z[71]*tmp1[71]+z[72]*tmp1[72]+z[73]*tmp1[73]+z[74]*tmp1[74]+z[75]*tmp1[75]+z[76]*tmp1[76]+z[77]*tmp1[77]+z[78]*tmp1[78]+z[79]*tmp1[79]+z[80]*tmp1[80]
            +z[81]*tmp1[81]+z[82]*tmp1[82]+z[83]*tmp1[83]+z[84]*tmp1[84]+z[85]*tmp1[85]+z[86]*tmp1[86]+z[87]*tmp1[87]+z[88]*tmp1[88]+z[89]*tmp1[89]+z[90]*tmp1[90]
            +z[91]*tmp1[91]+z[92]*tmp1[92]+z[93]*tmp1[93]+z[94]*tmp1[94]+z[95]*tmp1[95]+z[96]*tmp1[96];
            tab[i+2*stg_first]=tab[i+2*stg_first]+z[0]*tmp1[0]
            +z[2]*tmp1[1]+z[4]*tmp1[2]+z[6]*tmp1[3]+z[8]*tmp1[4]+z[10]*tmp1[5]+z[12]*tmp1[6]+z[14]*tmp1[7]+z[16]*tmp1[8]+z[18]*tmp1[9]+z[20]*tmp1[10]
            +z[22]*tmp1[11]+z[24]*tmp1[12]+z[26]*tmp1[13]+z[28]*tmp1[14]+z[30]*tmp1[15]+z[32]*tmp1[16]+z[34]*tmp1[17]+z[36]*tmp1[18]+z[38]*tmp1[19]+z[40]*tmp1[20]
            +z[42]*tmp1[21]+z[44]*tmp1[22]+z[46]*tmp1[23]+z[48]*tmp1[24]+z[50]*tmp1[25]+z[52]*tmp1[26]+z[54]*tmp1[27]+z[56]*tmp1[28]+z[58]*tmp1[29]+z[60]*tmp1[30]
            +z[62]*tmp1[31]+z[64]*tmp1[32]+z[66]*tmp1[33]+z[68]*tmp1[34]+z[70]*tmp1[35]+z[72]*tmp1[36]+z[74]*tmp1[37]+z[76]*tmp1[38]+z[78]*tmp1[39]+z[80]*tmp1[40]
            +z[82]*tmp1[41]+z[84]*tmp1[42]+z[86]*tmp1[43]+z[88]*tmp1[44]+z[90]*tmp1[45]+z[92]*tmp1[46]+z[94]*tmp1[47]+z[96]*tmp1[48]+z[1]*tmp1[49]+z[3]*tmp1[50]
            +z[5]*tmp1[51]+z[7]*tmp1[52]+z[9]*tmp1[53]+z[11]*tmp1[54]+z[13]*tmp1[55]+z[15]*tmp1[56]+z[17]*tmp1[57]+z[19]*tmp1[58]+z[21]*tmp1[59]+z[23]*tmp1[60]
            +z[25]*tmp1[61]+z[27]*tmp1[62]+z[29]*tmp1[63]+z[31]*tmp1[64]+z[33]*tmp1[65]+z[35]*tmp1[66]+z[37]*tmp1[67]+z[39]*tmp1[68]+z[41]*tmp1[69]+z[43]*tmp1[70]
            +z[45]*tmp1[71]+z[47]*tmp1[72]+z[49]*tmp1[73]+z[51]*tmp1[74]+z[53]*tmp1[75]+z[55]*tmp1[76]+z[57]*tmp1[77]+z[59]*tmp1[78]+z[61]*tmp1[79]+z[63]*tmp1[80]
            +z[65]*tmp1[81]+z[67]*tmp1[82]+z[69]*tmp1[83]+z[71]*tmp1[84]+z[73]*tmp1[85]+z[75]*tmp1[86]+z[77]*tmp1[87]+z[79]*tmp1[88]+z[81]*tmp1[89]+z[83]*tmp1[90]
            +z[85]*tmp1[91]+z[87]*tmp1[92]+z[89]*tmp1[93]+z[91]*tmp1[94]+z[93]*tmp1[95]+z[95]*tmp1[96];
            tab[i+3*stg_first]=tab[i+3*stg_first]+z[0]*tmp1[0]
            +z[3]*tmp1[1]+z[6]*tmp1[2]+z[9]*tmp1[3]+z[12]*tmp1[4]+z[15]*tmp1[5]+z[18]*tmp1[6]+z[21]*tmp1[7]+z[24]*tmp1[8]+z[27]*tmp1[9]+z[30]*tmp1[10]
            +z[33]*tmp1[11]+z[36]*tmp1[12]+z[39]*tmp1[13]+z[42]*tmp1[14]+z[45]*tmp1[15]+z[48]*tmp1[16]+z[51]*tmp1[17]+z[54]*tmp1[18]+z[57]*tmp1[19]+z[60]*tmp1[20]
            +z[63]*tmp1[21]+z[66]*tmp1[22]+z[69]*tmp1[23]+z[72]*tmp1[24]+z[75]*tmp1[25]+z[78]*tmp1[26]+z[81]*tmp1[27]+z[84]*tmp1[28]+z[87]*tmp1[29]+z[90]*tmp1[30]
            +z[93]*tmp1[31]+z[96]*tmp1[32]+z[2]*tmp1[33]+z[5]*tmp1[34]+z[8]*tmp1[35]+z[11]*tmp1[36]+z[14]*tmp1[37]+z[17]*tmp1[38]+z[20]*tmp1[39]+z[23]*tmp1[40]
            +z[26]*tmp1[41]+z[29]*tmp1[42]+z[32]*tmp1[43]+z[35]*tmp1[44]+z[38]*tmp1[45]+z[41]*tmp1[46]+z[44]*tmp1[47]+z[47]*tmp1[48]+z[50]*tmp1[49]+z[53]*tmp1[50]
            +z[56]*tmp1[51]+z[59]*tmp1[52]+z[62]*tmp1[53]+z[65]*tmp1[54]+z[68]*tmp1[55]+z[71]*tmp1[56]+z[74]*tmp1[57]+z[77]*tmp1[58]+z[80]*tmp1[59]+z[83]*tmp1[60]
            +z[86]*tmp1[61]+z[89]*tmp1[62]+z[92]*tmp1[63]+z[95]*tmp1[64]+z[1]*tmp1[65]+z[4]*tmp1[66]+z[7]*tmp1[67]+z[10]*tmp1[68]+z[13]*tmp1[69]+z[16]*tmp1[70]
            +z[19]*tmp1[71]+z[22]*tmp1[72]+z[25]*tmp1[73]+z[28]*tmp1[74]+z[31]*tmp1[75]+z[34]*tmp1[76]+z[37]*tmp1[77]+z[40]*tmp1[78]+z[43]*tmp1[79]+z[46]*tmp1[80]
            +z[49]*tmp1[81]+z[52]*tmp1[82]+z[55]*tmp1[83]+z[58]*tmp1[84]+z[61]*tmp1[85]+z[64]*tmp1[86]+z[67]*tmp1[87]+z[70]*tmp1[88]+z[73]*tmp1[89]+z[76]*tmp1[90]
            +z[79]*tmp1[91]+z[82]*tmp1[92]+z[85]*tmp1[93]+z[88]*tmp1[94]+z[91]*tmp1[95]+z[94]*tmp1[96];
            tab[i+4*stg_first]=tab[i+4*stg_first]+z[0]*tmp1[0]
            +z[4]*tmp1[1]+z[8]*tmp1[2]+z[12]*tmp1[3]+z[16]*tmp1[4]+z[20]*tmp1[5]+z[24]*tmp1[6]+z[28]*tmp1[7]+z[32]*tmp1[8]+z[36]*tmp1[9]+z[40]*tmp1[10]
            +z[44]*tmp1[11]+z[48]*tmp1[12]+z[52]*tmp1[13]+z[56]*tmp1[14]+z[60]*tmp1[15]+z[64]*tmp1[16]+z[68]*tmp1[17]+z[72]*tmp1[18]+z[76]*tmp1[19]+z[80]*tmp1[20]
            +z[84]*tmp1[21]+z[88]*tmp1[22]+z[92]*tmp1[23]+z[96]*tmp1[24]+z[3]*tmp1[25]+z[7]*tmp1[26]+z[11]*tmp1[27]+z[15]*tmp1[28]+z[19]*tmp1[29]+z[23]*tmp1[30]
            +z[27]*tmp1[31]+z[31]*tmp1[32]+z[35]*tmp1[33]+z[39]*tmp1[34]+z[43]*tmp1[35]+z[47]*tmp1[36]+z[51]*tmp1[37]+z[55]*tmp1[38]+z[59]*tmp1[39]+z[63]*tmp1[40]
            +z[67]*tmp1[41]+z[71]*tmp1[42]+z[75]*tmp1[43]+z[79]*tmp1[44]+z[83]*tmp1[45]+z[87]*tmp1[46]+z[91]*tmp1[47]+z[95]*tmp1[48]+z[2]*tmp1[49]+z[6]*tmp1[50]
            +z[10]*tmp1[51]+z[14]*tmp1[52]+z[18]*tmp1[53]+z[22]*tmp1[54]+z[26]*tmp1[55]+z[30]*tmp1[56]+z[34]*tmp1[57]+z[38]*tmp1[58]+z[42]*tmp1[59]+z[46]*tmp1[60]
            +z[50]*tmp1[61]+z[54]*tmp1[62]+z[58]*tmp1[63]+z[62]*tmp1[64]+z[66]*tmp1[65]+z[70]*tmp1[66]+z[74]*tmp1[67]+z[78]*tmp1[68]+z[82]*tmp1[69]+z[86]*tmp1[70]
            +z[90]*tmp1[71]+z[94]*tmp1[72]+z[1]*tmp1[73]+z[5]*tmp1[74]+z[9]*tmp1[75]+z[13]*tmp1[76]+z[17]*tmp1[77]+z[21]*tmp1[78]+z[25]*tmp1[79]+z[29]*tmp1[80]
            +z[33]*tmp1[81]+z[37]*tmp1[82]+z[41]*tmp1[83]+z[45]*tmp1[84]+z[49]*tmp1[85]+z[53]*tmp1[86]+z[57]*tmp1[87]+z[61]*tmp1[88]+z[65]*tmp1[89]+z[69]*tmp1[90]
            +z[73]*tmp1[91]+z[77]*tmp1[92]+z[81]*tmp1[93]+z[85]*tmp1[94]+z[89]*tmp1[95]+z[93]*tmp1[96];
            tab[i+5*stg_first]=tab[i+5*stg_first]+z[0]*tmp1[0]
            +z[5]*tmp1[1]+z[10]*tmp1[2]+z[15]*tmp1[3]+z[20]*tmp1[4]+z[25]*tmp1[5]+z[30]*tmp1[6]+z[35]*tmp1[7]+z[40]*tmp1[8]+z[45]*tmp1[9]+z[50]*tmp1[10]
            +z[55]*tmp1[11]+z[60]*tmp1[12]+z[65]*tmp1[13]+z[70]*tmp1[14]+z[75]*tmp1[15]+z[80]*tmp1[16]+z[85]*tmp1[17]+z[90]*tmp1[18]+z[95]*tmp1[19]+z[3]*tmp1[20]
            +z[8]*tmp1[21]+z[13]*tmp1[22]+z[18]*tmp1[23]+z[23]*tmp1[24]+z[28]*tmp1[25]+z[33]*tmp1[26]+z[38]*tmp1[27]+z[43]*tmp1[28]+z[48]*tmp1[29]+z[53]*tmp1[30]
            +z[58]*tmp1[31]+z[63]*tmp1[32]+z[68]*tmp1[33]+z[73]*tmp1[34]+z[78]*tmp1[35]+z[83]*tmp1[36]+z[88]*tmp1[37]+z[93]*tmp1[38]+z[1]*tmp1[39]+z[6]*tmp1[40]
            +z[11]*tmp1[41]+z[16]*tmp1[42]+z[21]*tmp1[43]+z[26]*tmp1[44]+z[31]*tmp1[45]+z[36]*tmp1[46]+z[41]*tmp1[47]+z[46]*tmp1[48]+z[51]*tmp1[49]+z[56]*tmp1[50]
            +z[61]*tmp1[51]+z[66]*tmp1[52]+z[71]*tmp1[53]+z[76]*tmp1[54]+z[81]*tmp1[55]+z[86]*tmp1[56]+z[91]*tmp1[57]+z[96]*tmp1[58]+z[4]*tmp1[59]+z[9]*tmp1[60]
            +z[14]*tmp1[61]+z[19]*tmp1[62]+z[24]*tmp1[63]+z[29]*tmp1[64]+z[34]*tmp1[65]+z[39]*tmp1[66]+z[44]*tmp1[67]+z[49]*tmp1[68]+z[54]*tmp1[69]+z[59]*tmp1[70]
            +z[64]*tmp1[71]+z[69]*tmp1[72]+z[74]*tmp1[73]+z[79]*tmp1[74]+z[84]*tmp1[75]+z[89]*tmp1[76]+z[94]*tmp1[77]+z[2]*tmp1[78]+z[7]*tmp1[79]+z[12]*tmp1[80]
            +z[17]*tmp1[81]+z[22]*tmp1[82]+z[27]*tmp1[83]+z[32]*tmp1[84]+z[37]*tmp1[85]+z[42]*tmp1[86]+z[47]*tmp1[87]+z[52]*tmp1[88]+z[57]*tmp1[89]+z[62]*tmp1[90]
            +z[67]*tmp1[91]+z[72]*tmp1[92]+z[77]*tmp1[93]+z[82]*tmp1[94]+z[87]*tmp1[95]+z[92]*tmp1[96];
            tab[i+6*stg_first]=tab[i+6*stg_first]+z[0]*tmp1[0]
            +z[6]*tmp1[1]+z[12]*tmp1[2]+z[18]*tmp1[3]+z[24]*tmp1[4]+z[30]*tmp1[5]+z[36]*tmp1[6]+z[42]*tmp1[7]+z[48]*tmp1[8]+z[54]*tmp1[9]+z[60]*tmp1[10]
            +z[66]*tmp1[11]+z[72]*tmp1[12]+z[78]*tmp1[13]+z[84]*tmp1[14]+z[90]*tmp1[15]+z[96]*tmp1[16]+z[5]*tmp1[17]+z[11]*tmp1[18]+z[17]*tmp1[19]+z[23]*tmp1[20]
            +z[29]*tmp1[21]+z[35]*tmp1[22]+z[41]*tmp1[23]+z[47]*tmp1[24]+z[53]*tmp1[25]+z[59]*tmp1[26]+z[65]*tmp1[27]+z[71]*tmp1[28]+z[77]*tmp1[29]+z[83]*tmp1[30]
            +z[89]*tmp1[31]+z[95]*tmp1[32]+z[4]*tmp1[33]+z[10]*tmp1[34]+z[16]*tmp1[35]+z[22]*tmp1[36]+z[28]*tmp1[37]+z[34]*tmp1[38]+z[40]*tmp1[39]+z[46]*tmp1[40]
            +z[52]*tmp1[41]+z[58]*tmp1[42]+z[64]*tmp1[43]+z[70]*tmp1[44]+z[76]*tmp1[45]+z[82]*tmp1[46]+z[88]*tmp1[47]+z[94]*tmp1[48]+z[3]*tmp1[49]+z[9]*tmp1[50]
            +z[15]*tmp1[51]+z[21]*tmp1[52]+z[27]*tmp1[53]+z[33]*tmp1[54]+z[39]*tmp1[55]+z[45]*tmp1[56]+z[51]*tmp1[57]+z[57]*tmp1[58]+z[63]*tmp1[59]+z[69]*tmp1[60]
            +z[75]*tmp1[61]+z[81]*tmp1[62]+z[87]*tmp1[63]+z[93]*tmp1[64]+z[2]*tmp1[65]+z[8]*tmp1[66]+z[14]*tmp1[67]+z[20]*tmp1[68]+z[26]*tmp1[69]+z[32]*tmp1[70]
            +z[38]*tmp1[71]+z[44]*tmp1[72]+z[50]*tmp1[73]+z[56]*tmp1[74]+z[62]*tmp1[75]+z[68]*tmp1[76]+z[74]*tmp1[77]+z[80]*tmp1[78]+z[86]*tmp1[79]+z[92]*tmp1[80]
            +z[1]*tmp1[81]+z[7]*tmp1[82]+z[13]*tmp1[83]+z[19]*tmp1[84]+z[25]*tmp1[85]+z[31]*tmp1[86]+z[37]*tmp1[87]+z[43]*tmp1[88]+z[49]*tmp1[89]+z[55]*tmp1[90]
            +z[61]*tmp1[91]+z[67]*tmp1[92]+z[73]*tmp1[93]+z[79]*tmp1[94]+z[85]*tmp1[95]+z[91]*tmp1[96];
            tab[i+7*stg_first]=tab[i+7*stg_first]+z[0]*tmp1[0]
            +z[7]*tmp1[1]+z[14]*tmp1[2]+z[21]*tmp1[3]+z[28]*tmp1[4]+z[35]*tmp1[5]+z[42]*tmp1[6]+z[49]*tmp1[7]+z[56]*tmp1[8]+z[63]*tmp1[9]+z[70]*tmp1[10]
            +z[77]*tmp1[11]+z[84]*tmp1[12]+z[91]*tmp1[13]+z[1]*tmp1[14]+z[8]*tmp1[15]+z[15]*tmp1[16]+z[22]*tmp1[17]+z[29]*tmp1[18]+z[36]*tmp1[19]+z[43]*tmp1[20]
            +z[50]*tmp1[21]+z[57]*tmp1[22]+z[64]*tmp1[23]+z[71]*tmp1[24]+z[78]*tmp1[25]+z[85]*tmp1[26]+z[92]*tmp1[27]+z[2]*tmp1[28]+z[9]*tmp1[29]+z[16]*tmp1[30]
            +z[23]*tmp1[31]+z[30]*tmp1[32]+z[37]*tmp1[33]+z[44]*tmp1[34]+z[51]*tmp1[35]+z[58]*tmp1[36]+z[65]*tmp1[37]+z[72]*tmp1[38]+z[79]*tmp1[39]+z[86]*tmp1[40]
            +z[93]*tmp1[41]+z[3]*tmp1[42]+z[10]*tmp1[43]+z[17]*tmp1[44]+z[24]*tmp1[45]+z[31]*tmp1[46]+z[38]*tmp1[47]+z[45]*tmp1[48]+z[52]*tmp1[49]+z[59]*tmp1[50]
            +z[66]*tmp1[51]+z[73]*tmp1[52]+z[80]*tmp1[53]+z[87]*tmp1[54]+z[94]*tmp1[55]+z[4]*tmp1[56]+z[11]*tmp1[57]+z[18]*tmp1[58]+z[25]*tmp1[59]+z[32]*tmp1[60]
            +z[39]*tmp1[61]+z[46]*tmp1[62]+z[53]*tmp1[63]+z[60]*tmp1[64]+z[67]*tmp1[65]+z[74]*tmp1[66]+z[81]*tmp1[67]+z[88]*tmp1[68]+z[95]*tmp1[69]+z[5]*tmp1[70]
            +z[12]*tmp1[71]+z[19]*tmp1[72]+z[26]*tmp1[73]+z[33]*tmp1[74]+z[40]*tmp1[75]+z[47]*tmp1[76]+z[54]*tmp1[77]+z[61]*tmp1[78]+z[68]*tmp1[79]+z[75]*tmp1[80]
            +z[82]*tmp1[81]+z[89]*tmp1[82]+z[96]*tmp1[83]+z[6]*tmp1[84]+z[13]*tmp1[85]+z[20]*tmp1[86]+z[27]*tmp1[87]+z[34]*tmp1[88]+z[41]*tmp1[89]+z[48]*tmp1[90]
            +z[55]*tmp1[91]+z[62]*tmp1[92]+z[69]*tmp1[93]+z[76]*tmp1[94]+z[83]*tmp1[95]+z[90]*tmp1[96];
            tab[i+8*stg_first]=tab[i+8*stg_first]+z[0]*tmp1[0]
            +z[8]*tmp1[1]+z[16]*tmp1[2]+z[24]*tmp1[3]+z[32]*tmp1[4]+z[40]*tmp1[5]+z[48]*tmp1[6]+z[56]*tmp1[7]+z[64]*tmp1[8]+z[72]*tmp1[9]+z[80]*tmp1[10]
            +z[88]*tmp1[11]+z[96]*tmp1[12]+z[7]*tmp1[13]+z[15]*tmp1[14]+z[23]*tmp1[15]+z[31]*tmp1[16]+z[39]*tmp1[17]+z[47]*tmp1[18]+z[55]*tmp1[19]+z[63]*tmp1[20]
            +z[71]*tmp1[21]+z[79]*tmp1[22]+z[87]*tmp1[23]+z[95]*tmp1[24]+z[6]*tmp1[25]+z[14]*tmp1[26]+z[22]*tmp1[27]+z[30]*tmp1[28]+z[38]*tmp1[29]+z[46]*tmp1[30]
            +z[54]*tmp1[31]+z[62]*tmp1[32]+z[70]*tmp1[33]+z[78]*tmp1[34]+z[86]*tmp1[35]+z[94]*tmp1[36]+z[5]*tmp1[37]+z[13]*tmp1[38]+z[21]*tmp1[39]+z[29]*tmp1[40]
            +z[37]*tmp1[41]+z[45]*tmp1[42]+z[53]*tmp1[43]+z[61]*tmp1[44]+z[69]*tmp1[45]+z[77]*tmp1[46]+z[85]*tmp1[47]+z[93]*tmp1[48]+z[4]*tmp1[49]+z[12]*tmp1[50]
            +z[20]*tmp1[51]+z[28]*tmp1[52]+z[36]*tmp1[53]+z[44]*tmp1[54]+z[52]*tmp1[55]+z[60]*tmp1[56]+z[68]*tmp1[57]+z[76]*tmp1[58]+z[84]*tmp1[59]+z[92]*tmp1[60]
            +z[3]*tmp1[61]+z[11]*tmp1[62]+z[19]*tmp1[63]+z[27]*tmp1[64]+z[35]*tmp1[65]+z[43]*tmp1[66]+z[51]*tmp1[67]+z[59]*tmp1[68]+z[67]*tmp1[69]+z[75]*tmp1[70]
            +z[83]*tmp1[71]+z[91]*tmp1[72]+z[2]*tmp1[73]+z[10]*tmp1[74]+z[18]*tmp1[75]+z[26]*tmp1[76]+z[34]*tmp1[77]+z[42]*tmp1[78]+z[50]*tmp1[79]+z[58]*tmp1[80]
            +z[66]*tmp1[81]+z[74]*tmp1[82]+z[82]*tmp1[83]+z[90]*tmp1[84]+z[1]*tmp1[85]+z[9]*tmp1[86]+z[17]*tmp1[87]+z[25]*tmp1[88]+z[33]*tmp1[89]+z[41]*tmp1[90]
            +z[49]*tmp1[91]+z[57]*tmp1[92]+z[65]*tmp1[93]+z[73]*tmp1[94]+z[81]*tmp1[95]+z[89]*tmp1[96];
            tab[i+9*stg_first]=tab[i+9*stg_first]+z[0]*tmp1[0]
            +z[9]*tmp1[1]+z[18]*tmp1[2]+z[27]*tmp1[3]+z[36]*tmp1[4]+z[45]*tmp1[5]+z[54]*tmp1[6]+z[63]*tmp1[7]+z[72]*tmp1[8]+z[81]*tmp1[9]+z[90]*tmp1[10]
            +z[2]*tmp1[11]+z[11]*tmp1[12]+z[20]*tmp1[13]+z[29]*tmp1[14]+z[38]*tmp1[15]+z[47]*tmp1[16]+z[56]*tmp1[17]+z[65]*tmp1[18]+z[74]*tmp1[19]+z[83]*tmp1[20]
            +z[92]*tmp1[21]+z[4]*tmp1[22]+z[13]*tmp1[23]+z[22]*tmp1[24]+z[31]*tmp1[25]+z[40]*tmp1[26]+z[49]*tmp1[27]+z[58]*tmp1[28]+z[67]*tmp1[29]+z[76]*tmp1[30]
            +z[85]*tmp1[31]+z[94]*tmp1[32]+z[6]*tmp1[33]+z[15]*tmp1[34]+z[24]*tmp1[35]+z[33]*tmp1[36]+z[42]*tmp1[37]+z[51]*tmp1[38]+z[60]*tmp1[39]+z[69]*tmp1[40]
            +z[78]*tmp1[41]+z[87]*tmp1[42]+z[96]*tmp1[43]+z[8]*tmp1[44]+z[17]*tmp1[45]+z[26]*tmp1[46]+z[35]*tmp1[47]+z[44]*tmp1[48]+z[53]*tmp1[49]+z[62]*tmp1[50]
            +z[71]*tmp1[51]+z[80]*tmp1[52]+z[89]*tmp1[53]+z[1]*tmp1[54]+z[10]*tmp1[55]+z[19]*tmp1[56]+z[28]*tmp1[57]+z[37]*tmp1[58]+z[46]*tmp1[59]+z[55]*tmp1[60]
            +z[64]*tmp1[61]+z[73]*tmp1[62]+z[82]*tmp1[63]+z[91]*tmp1[64]+z[3]*tmp1[65]+z[12]*tmp1[66]+z[21]*tmp1[67]+z[30]*tmp1[68]+z[39]*tmp1[69]+z[48]*tmp1[70]
            +z[57]*tmp1[71]+z[66]*tmp1[72]+z[75]*tmp1[73]+z[84]*tmp1[74]+z[93]*tmp1[75]+z[5]*tmp1[76]+z[14]*tmp1[77]+z[23]*tmp1[78]+z[32]*tmp1[79]+z[41]*tmp1[80]
            +z[50]*tmp1[81]+z[59]*tmp1[82]+z[68]*tmp1[83]+z[77]*tmp1[84]+z[86]*tmp1[85]+z[95]*tmp1[86]+z[7]*tmp1[87]+z[16]*tmp1[88]+z[25]*tmp1[89]+z[34]*tmp1[90]
            +z[43]*tmp1[91]+z[52]*tmp1[92]+z[61]*tmp1[93]+z[70]*tmp1[94]+z[79]*tmp1[95]+z[88]*tmp1[96];
            tab[i+10*stg_first]=tab[i+10*stg_first]+z[0]*tmp1[0]
            +z[10]*tmp1[1]+z[20]*tmp1[2]+z[30]*tmp1[3]+z[40]*tmp1[4]+z[50]*tmp1[5]+z[60]*tmp1[6]+z[70]*tmp1[7]+z[80]*tmp1[8]+z[90]*tmp1[9]+z[3]*tmp1[10]
            +z[13]*tmp1[11]+z[23]*tmp1[12]+z[33]*tmp1[13]+z[43]*tmp1[14]+z[53]*tmp1[15]+z[63]*tmp1[16]+z[73]*tmp1[17]+z[83]*tmp1[18]+z[93]*tmp1[19]+z[6]*tmp1[20]
            +z[16]*tmp1[21]+z[26]*tmp1[22]+z[36]*tmp1[23]+z[46]*tmp1[24]+z[56]*tmp1[25]+z[66]*tmp1[26]+z[76]*tmp1[27]+z[86]*tmp1[28]+z[96]*tmp1[29]+z[9]*tmp1[30]
            +z[19]*tmp1[31]+z[29]*tmp1[32]+z[39]*tmp1[33]+z[49]*tmp1[34]+z[59]*tmp1[35]+z[69]*tmp1[36]+z[79]*tmp1[37]+z[89]*tmp1[38]+z[2]*tmp1[39]+z[12]*tmp1[40]
            +z[22]*tmp1[41]+z[32]*tmp1[42]+z[42]*tmp1[43]+z[52]*tmp1[44]+z[62]*tmp1[45]+z[72]*tmp1[46]+z[82]*tmp1[47]+z[92]*tmp1[48]+z[5]*tmp1[49]+z[15]*tmp1[50]
            +z[25]*tmp1[51]+z[35]*tmp1[52]+z[45]*tmp1[53]+z[55]*tmp1[54]+z[65]*tmp1[55]+z[75]*tmp1[56]+z[85]*tmp1[57]+z[95]*tmp1[58]+z[8]*tmp1[59]+z[18]*tmp1[60]
            +z[28]*tmp1[61]+z[38]*tmp1[62]+z[48]*tmp1[63]+z[58]*tmp1[64]+z[68]*tmp1[65]+z[78]*tmp1[66]+z[88]*tmp1[67]+z[1]*tmp1[68]+z[11]*tmp1[69]+z[21]*tmp1[70]
            +z[31]*tmp1[71]+z[41]*tmp1[72]+z[51]*tmp1[73]+z[61]*tmp1[74]+z[71]*tmp1[75]+z[81]*tmp1[76]+z[91]*tmp1[77]+z[4]*tmp1[78]+z[14]*tmp1[79]+z[24]*tmp1[80]
            +z[34]*tmp1[81]+z[44]*tmp1[82]+z[54]*tmp1[83]+z[64]*tmp1[84]+z[74]*tmp1[85]+z[84]*tmp1[86]+z[94]*tmp1[87]+z[7]*tmp1[88]+z[17]*tmp1[89]+z[27]*tmp1[90]
            +z[37]*tmp1[91]+z[47]*tmp1[92]+z[57]*tmp1[93]+z[67]*tmp1[94]+z[77]*tmp1[95]+z[87]*tmp1[96];
            tab[i+11*stg_first]=tab[i+11*stg_first]+z[0]*tmp1[0]
            +z[11]*tmp1[1]+z[22]*tmp1[2]+z[33]*tmp1[3]+z[44]*tmp1[4]+z[55]*tmp1[5]+z[66]*tmp1[6]+z[77]*tmp1[7]+z[88]*tmp1[8]+z[2]*tmp1[9]+z[13]*tmp1[10]
            +z[24]*tmp1[11]+z[35]*tmp1[12]+z[46]*tmp1[13]+z[57]*tmp1[14]+z[68]*tmp1[15]+z[79]*tmp1[16]+z[90]*tmp1[17]+z[4]*tmp1[18]+z[15]*tmp1[19]+z[26]*tmp1[20]
            +z[37]*tmp1[21]+z[48]*tmp1[22]+z[59]*tmp1[23]+z[70]*tmp1[24]+z[81]*tmp1[25]+z[92]*tmp1[26]+z[6]*tmp1[27]+z[17]*tmp1[28]+z[28]*tmp1[29]+z[39]*tmp1[30]
            +z[50]*tmp1[31]+z[61]*tmp1[32]+z[72]*tmp1[33]+z[83]*tmp1[34]+z[94]*tmp1[35]+z[8]*tmp1[36]+z[19]*tmp1[37]+z[30]*tmp1[38]+z[41]*tmp1[39]+z[52]*tmp1[40]
            +z[63]*tmp1[41]+z[74]*tmp1[42]+z[85]*tmp1[43]+z[96]*tmp1[44]+z[10]*tmp1[45]+z[21]*tmp1[46]+z[32]*tmp1[47]+z[43]*tmp1[48]+z[54]*tmp1[49]+z[65]*tmp1[50]
            +z[76]*tmp1[51]+z[87]*tmp1[52]+z[1]*tmp1[53]+z[12]*tmp1[54]+z[23]*tmp1[55]+z[34]*tmp1[56]+z[45]*tmp1[57]+z[56]*tmp1[58]+z[67]*tmp1[59]+z[78]*tmp1[60]
            +z[89]*tmp1[61]+z[3]*tmp1[62]+z[14]*tmp1[63]+z[25]*tmp1[64]+z[36]*tmp1[65]+z[47]*tmp1[66]+z[58]*tmp1[67]+z[69]*tmp1[68]+z[80]*tmp1[69]+z[91]*tmp1[70]
            +z[5]*tmp1[71]+z[16]*tmp1[72]+z[27]*tmp1[73]+z[38]*tmp1[74]+z[49]*tmp1[75]+z[60]*tmp1[76]+z[71]*tmp1[77]+z[82]*tmp1[78]+z[93]*tmp1[79]+z[7]*tmp1[80]
            +z[18]*tmp1[81]+z[29]*tmp1[82]+z[40]*tmp1[83]+z[51]*tmp1[84]+z[62]*tmp1[85]+z[73]*tmp1[86]+z[84]*tmp1[87]+z[95]*tmp1[88]+z[9]*tmp1[89]+z[20]*tmp1[90]
            +z[31]*tmp1[91]+z[42]*tmp1[92]+z[53]*tmp1[93]+z[64]*tmp1[94]+z[75]*tmp1[95]+z[86]*tmp1[96];
            tab[i+12*stg_first]=tab[i+12*stg_first]+z[0]*tmp1[0]
            +z[12]*tmp1[1]+z[24]*tmp1[2]+z[36]*tmp1[3]+z[48]*tmp1[4]+z[60]*tmp1[5]+z[72]*tmp1[6]+z[84]*tmp1[7]+z[96]*tmp1[8]+z[11]*tmp1[9]+z[23]*tmp1[10]
            +z[35]*tmp1[11]+z[47]*tmp1[12]+z[59]*tmp1[13]+z[71]*tmp1[14]+z[83]*tmp1[15]+z[95]*tmp1[16]+z[10]*tmp1[17]+z[22]*tmp1[18]+z[34]*tmp1[19]+z[46]*tmp1[20]
            +z[58]*tmp1[21]+z[70]*tmp1[22]+z[82]*tmp1[23]+z[94]*tmp1[24]+z[9]*tmp1[25]+z[21]*tmp1[26]+z[33]*tmp1[27]+z[45]*tmp1[28]+z[57]*tmp1[29]+z[69]*tmp1[30]
            +z[81]*tmp1[31]+z[93]*tmp1[32]+z[8]*tmp1[33]+z[20]*tmp1[34]+z[32]*tmp1[35]+z[44]*tmp1[36]+z[56]*tmp1[37]+z[68]*tmp1[38]+z[80]*tmp1[39]+z[92]*tmp1[40]
            +z[7]*tmp1[41]+z[19]*tmp1[42]+z[31]*tmp1[43]+z[43]*tmp1[44]+z[55]*tmp1[45]+z[67]*tmp1[46]+z[79]*tmp1[47]+z[91]*tmp1[48]+z[6]*tmp1[49]+z[18]*tmp1[50]
            +z[30]*tmp1[51]+z[42]*tmp1[52]+z[54]*tmp1[53]+z[66]*tmp1[54]+z[78]*tmp1[55]+z[90]*tmp1[56]+z[5]*tmp1[57]+z[17]*tmp1[58]+z[29]*tmp1[59]+z[41]*tmp1[60]
            +z[53]*tmp1[61]+z[65]*tmp1[62]+z[77]*tmp1[63]+z[89]*tmp1[64]+z[4]*tmp1[65]+z[16]*tmp1[66]+z[28]*tmp1[67]+z[40]*tmp1[68]+z[52]*tmp1[69]+z[64]*tmp1[70]
            +z[76]*tmp1[71]+z[88]*tmp1[72]+z[3]*tmp1[73]+z[15]*tmp1[74]+z[27]*tmp1[75]+z[39]*tmp1[76]+z[51]*tmp1[77]+z[63]*tmp1[78]+z[75]*tmp1[79]+z[87]*tmp1[80]
            +z[2]*tmp1[81]+z[14]*tmp1[82]+z[26]*tmp1[83]+z[38]*tmp1[84]+z[50]*tmp1[85]+z[62]*tmp1[86]+z[74]*tmp1[87]+z[86]*tmp1[88]+z[1]*tmp1[89]+z[13]*tmp1[90]
            +z[25]*tmp1[91]+z[37]*tmp1[92]+z[49]*tmp1[93]+z[61]*tmp1[94]+z[73]*tmp1[95]+z[85]*tmp1[96];
            tab[i+13*stg_first]=tab[i+13*stg_first]+z[0]*tmp1[0]
            +z[13]*tmp1[1]+z[26]*tmp1[2]+z[39]*tmp1[3]+z[52]*tmp1[4]+z[65]*tmp1[5]+z[78]*tmp1[6]+z[91]*tmp1[7]+z[7]*tmp1[8]+z[20]*tmp1[9]+z[33]*tmp1[10]
            +z[46]*tmp1[11]+z[59]*tmp1[12]+z[72]*tmp1[13]+z[85]*tmp1[14]+z[1]*tmp1[15]+z[14]*tmp1[16]+z[27]*tmp1[17]+z[40]*tmp1[18]+z[53]*tmp1[19]+z[66]*tmp1[20]
            +z[79]*tmp1[21]+z[92]*tmp1[22]+z[8]*tmp1[23]+z[21]*tmp1[24]+z[34]*tmp1[25]+z[47]*tmp1[26]+z[60]*tmp1[27]+z[73]*tmp1[28]+z[86]*tmp1[29]+z[2]*tmp1[30]
            +z[15]*tmp1[31]+z[28]*tmp1[32]+z[41]*tmp1[33]+z[54]*tmp1[34]+z[67]*tmp1[35]+z[80]*tmp1[36]+z[93]*tmp1[37]+z[9]*tmp1[38]+z[22]*tmp1[39]+z[35]*tmp1[40]
            +z[48]*tmp1[41]+z[61]*tmp1[42]+z[74]*tmp1[43]+z[87]*tmp1[44]+z[3]*tmp1[45]+z[16]*tmp1[46]+z[29]*tmp1[47]+z[42]*tmp1[48]+z[55]*tmp1[49]+z[68]*tmp1[50]
            +z[81]*tmp1[51]+z[94]*tmp1[52]+z[10]*tmp1[53]+z[23]*tmp1[54]+z[36]*tmp1[55]+z[49]*tmp1[56]+z[62]*tmp1[57]+z[75]*tmp1[58]+z[88]*tmp1[59]+z[4]*tmp1[60]
            +z[17]*tmp1[61]+z[30]*tmp1[62]+z[43]*tmp1[63]+z[56]*tmp1[64]+z[69]*tmp1[65]+z[82]*tmp1[66]+z[95]*tmp1[67]+z[11]*tmp1[68]+z[24]*tmp1[69]+z[37]*tmp1[70]
            +z[50]*tmp1[71]+z[63]*tmp1[72]+z[76]*tmp1[73]+z[89]*tmp1[74]+z[5]*tmp1[75]+z[18]*tmp1[76]+z[31]*tmp1[77]+z[44]*tmp1[78]+z[57]*tmp1[79]+z[70]*tmp1[80]
            +z[83]*tmp1[81]+z[96]*tmp1[82]+z[12]*tmp1[83]+z[25]*tmp1[84]+z[38]*tmp1[85]+z[51]*tmp1[86]+z[64]*tmp1[87]+z[77]*tmp1[88]+z[90]*tmp1[89]+z[6]*tmp1[90]
            +z[19]*tmp1[91]+z[32]*tmp1[92]+z[45]*tmp1[93]+z[58]*tmp1[94]+z[71]*tmp1[95]+z[84]*tmp1[96];
            tab[i+14*stg_first]=tab[i+14*stg_first]+z[0]*tmp1[0]
            +z[14]*tmp1[1]+z[28]*tmp1[2]+z[42]*tmp1[3]+z[56]*tmp1[4]+z[70]*tmp1[5]+z[84]*tmp1[6]+z[1]*tmp1[7]+z[15]*tmp1[8]+z[29]*tmp1[9]+z[43]*tmp1[10]
            +z[57]*tmp1[11]+z[71]*tmp1[12]+z[85]*tmp1[13]+z[2]*tmp1[14]+z[16]*tmp1[15]+z[30]*tmp1[16]+z[44]*tmp1[17]+z[58]*tmp1[18]+z[72]*tmp1[19]+z[86]*tmp1[20]
            +z[3]*tmp1[21]+z[17]*tmp1[22]+z[31]*tmp1[23]+z[45]*tmp1[24]+z[59]*tmp1[25]+z[73]*tmp1[26]+z[87]*tmp1[27]+z[4]*tmp1[28]+z[18]*tmp1[29]+z[32]*tmp1[30]
            +z[46]*tmp1[31]+z[60]*tmp1[32]+z[74]*tmp1[33]+z[88]*tmp1[34]+z[5]*tmp1[35]+z[19]*tmp1[36]+z[33]*tmp1[37]+z[47]*tmp1[38]+z[61]*tmp1[39]+z[75]*tmp1[40]
            +z[89]*tmp1[41]+z[6]*tmp1[42]+z[20]*tmp1[43]+z[34]*tmp1[44]+z[48]*tmp1[45]+z[62]*tmp1[46]+z[76]*tmp1[47]+z[90]*tmp1[48]+z[7]*tmp1[49]+z[21]*tmp1[50]
            +z[35]*tmp1[51]+z[49]*tmp1[52]+z[63]*tmp1[53]+z[77]*tmp1[54]+z[91]*tmp1[55]+z[8]*tmp1[56]+z[22]*tmp1[57]+z[36]*tmp1[58]+z[50]*tmp1[59]+z[64]*tmp1[60]
            +z[78]*tmp1[61]+z[92]*tmp1[62]+z[9]*tmp1[63]+z[23]*tmp1[64]+z[37]*tmp1[65]+z[51]*tmp1[66]+z[65]*tmp1[67]+z[79]*tmp1[68]+z[93]*tmp1[69]+z[10]*tmp1[70]
            +z[24]*tmp1[71]+z[38]*tmp1[72]+z[52]*tmp1[73]+z[66]*tmp1[74]+z[80]*tmp1[75]+z[94]*tmp1[76]+z[11]*tmp1[77]+z[25]*tmp1[78]+z[39]*tmp1[79]+z[53]*tmp1[80]
            +z[67]*tmp1[81]+z[81]*tmp1[82]+z[95]*tmp1[83]+z[12]*tmp1[84]+z[26]*tmp1[85]+z[40]*tmp1[86]+z[54]*tmp1[87]+z[68]*tmp1[88]+z[82]*tmp1[89]+z[96]*tmp1[90]
            +z[13]*tmp1[91]+z[27]*tmp1[92]+z[41]*tmp1[93]+z[55]*tmp1[94]+z[69]*tmp1[95]+z[83]*tmp1[96];
            tab[i+15*stg_first]=tab[i+15*stg_first]+z[0]*tmp1[0]
            +z[15]*tmp1[1]+z[30]*tmp1[2]+z[45]*tmp1[3]+z[60]*tmp1[4]+z[75]*tmp1[5]+z[90]*tmp1[6]+z[8]*tmp1[7]+z[23]*tmp1[8]+z[38]*tmp1[9]+z[53]*tmp1[10]
            +z[68]*tmp1[11]+z[83]*tmp1[12]+z[1]*tmp1[13]+z[16]*tmp1[14]+z[31]*tmp1[15]+z[46]*tmp1[16]+z[61]*tmp1[17]+z[76]*tmp1[18]+z[91]*tmp1[19]+z[9]*tmp1[20]
            +z[24]*tmp1[21]+z[39]*tmp1[22]+z[54]*tmp1[23]+z[69]*tmp1[24]+z[84]*tmp1[25]+z[2]*tmp1[26]+z[17]*tmp1[27]+z[32]*tmp1[28]+z[47]*tmp1[29]+z[62]*tmp1[30]
            +z[77]*tmp1[31]+z[92]*tmp1[32]+z[10]*tmp1[33]+z[25]*tmp1[34]+z[40]*tmp1[35]+z[55]*tmp1[36]+z[70]*tmp1[37]+z[85]*tmp1[38]+z[3]*tmp1[39]+z[18]*tmp1[40]
            +z[33]*tmp1[41]+z[48]*tmp1[42]+z[63]*tmp1[43]+z[78]*tmp1[44]+z[93]*tmp1[45]+z[11]*tmp1[46]+z[26]*tmp1[47]+z[41]*tmp1[48]+z[56]*tmp1[49]+z[71]*tmp1[50]
            +z[86]*tmp1[51]+z[4]*tmp1[52]+z[19]*tmp1[53]+z[34]*tmp1[54]+z[49]*tmp1[55]+z[64]*tmp1[56]+z[79]*tmp1[57]+z[94]*tmp1[58]+z[12]*tmp1[59]+z[27]*tmp1[60]
            +z[42]*tmp1[61]+z[57]*tmp1[62]+z[72]*tmp1[63]+z[87]*tmp1[64]+z[5]*tmp1[65]+z[20]*tmp1[66]+z[35]*tmp1[67]+z[50]*tmp1[68]+z[65]*tmp1[69]+z[80]*tmp1[70]
            +z[95]*tmp1[71]+z[13]*tmp1[72]+z[28]*tmp1[73]+z[43]*tmp1[74]+z[58]*tmp1[75]+z[73]*tmp1[76]+z[88]*tmp1[77]+z[6]*tmp1[78]+z[21]*tmp1[79]+z[36]*tmp1[80]
            +z[51]*tmp1[81]+z[66]*tmp1[82]+z[81]*tmp1[83]+z[96]*tmp1[84]+z[14]*tmp1[85]+z[29]*tmp1[86]+z[44]*tmp1[87]+z[59]*tmp1[88]+z[74]*tmp1[89]+z[89]*tmp1[90]
            +z[7]*tmp1[91]+z[22]*tmp1[92]+z[37]*tmp1[93]+z[52]*tmp1[94]+z[67]*tmp1[95]+z[82]*tmp1[96];
            tab[i+16*stg_first]=tab[i+16*stg_first]+z[0]*tmp1[0]
            +z[16]*tmp1[1]+z[32]*tmp1[2]+z[48]*tmp1[3]+z[64]*tmp1[4]+z[80]*tmp1[5]+z[96]*tmp1[6]+z[15]*tmp1[7]+z[31]*tmp1[8]+z[47]*tmp1[9]+z[63]*tmp1[10]
            +z[79]*tmp1[11]+z[95]*tmp1[12]+z[14]*tmp1[13]+z[30]*tmp1[14]+z[46]*tmp1[15]+z[62]*tmp1[16]+z[78]*tmp1[17]+z[94]*tmp1[18]+z[13]*tmp1[19]+z[29]*tmp1[20]
            +z[45]*tmp1[21]+z[61]*tmp1[22]+z[77]*tmp1[23]+z[93]*tmp1[24]+z[12]*tmp1[25]+z[28]*tmp1[26]+z[44]*tmp1[27]+z[60]*tmp1[28]+z[76]*tmp1[29]+z[92]*tmp1[30]
            +z[11]*tmp1[31]+z[27]*tmp1[32]+z[43]*tmp1[33]+z[59]*tmp1[34]+z[75]*tmp1[35]+z[91]*tmp1[36]+z[10]*tmp1[37]+z[26]*tmp1[38]+z[42]*tmp1[39]+z[58]*tmp1[40]
            +z[74]*tmp1[41]+z[90]*tmp1[42]+z[9]*tmp1[43]+z[25]*tmp1[44]+z[41]*tmp1[45]+z[57]*tmp1[46]+z[73]*tmp1[47]+z[89]*tmp1[48]+z[8]*tmp1[49]+z[24]*tmp1[50]
            +z[40]*tmp1[51]+z[56]*tmp1[52]+z[72]*tmp1[53]+z[88]*tmp1[54]+z[7]*tmp1[55]+z[23]*tmp1[56]+z[39]*tmp1[57]+z[55]*tmp1[58]+z[71]*tmp1[59]+z[87]*tmp1[60]
            +z[6]*tmp1[61]+z[22]*tmp1[62]+z[38]*tmp1[63]+z[54]*tmp1[64]+z[70]*tmp1[65]+z[86]*tmp1[66]+z[5]*tmp1[67]+z[21]*tmp1[68]+z[37]*tmp1[69]+z[53]*tmp1[70]
            +z[69]*tmp1[71]+z[85]*tmp1[72]+z[4]*tmp1[73]+z[20]*tmp1[74]+z[36]*tmp1[75]+z[52]*tmp1[76]+z[68]*tmp1[77]+z[84]*tmp1[78]+z[3]*tmp1[79]+z[19]*tmp1[80]
            +z[35]*tmp1[81]+z[51]*tmp1[82]+z[67]*tmp1[83]+z[83]*tmp1[84]+z[2]*tmp1[85]+z[18]*tmp1[86]+z[34]*tmp1[87]+z[50]*tmp1[88]+z[66]*tmp1[89]+z[82]*tmp1[90]
            +z[1]*tmp1[91]+z[17]*tmp1[92]+z[33]*tmp1[93]+z[49]*tmp1[94]+z[65]*tmp1[95]+z[81]*tmp1[96];
            tab[i+17*stg_first]=tab[i+17*stg_first]+z[0]*tmp1[0]
            +z[17]*tmp1[1]+z[34]*tmp1[2]+z[51]*tmp1[3]+z[68]*tmp1[4]+z[85]*tmp1[5]+z[5]*tmp1[6]+z[22]*tmp1[7]+z[39]*tmp1[8]+z[56]*tmp1[9]+z[73]*tmp1[10]
            +z[90]*tmp1[11]+z[10]*tmp1[12]+z[27]*tmp1[13]+z[44]*tmp1[14]+z[61]*tmp1[15]+z[78]*tmp1[16]+z[95]*tmp1[17]+z[15]*tmp1[18]+z[32]*tmp1[19]+z[49]*tmp1[20]
            +z[66]*tmp1[21]+z[83]*tmp1[22]+z[3]*tmp1[23]+z[20]*tmp1[24]+z[37]*tmp1[25]+z[54]*tmp1[26]+z[71]*tmp1[27]+z[88]*tmp1[28]+z[8]*tmp1[29]+z[25]*tmp1[30]
            +z[42]*tmp1[31]+z[59]*tmp1[32]+z[76]*tmp1[33]+z[93]*tmp1[34]+z[13]*tmp1[35]+z[30]*tmp1[36]+z[47]*tmp1[37]+z[64]*tmp1[38]+z[81]*tmp1[39]+z[1]*tmp1[40]
            +z[18]*tmp1[41]+z[35]*tmp1[42]+z[52]*tmp1[43]+z[69]*tmp1[44]+z[86]*tmp1[45]+z[6]*tmp1[46]+z[23]*tmp1[47]+z[40]*tmp1[48]+z[57]*tmp1[49]+z[74]*tmp1[50]
            +z[91]*tmp1[51]+z[11]*tmp1[52]+z[28]*tmp1[53]+z[45]*tmp1[54]+z[62]*tmp1[55]+z[79]*tmp1[56]+z[96]*tmp1[57]+z[16]*tmp1[58]+z[33]*tmp1[59]+z[50]*tmp1[60]
            +z[67]*tmp1[61]+z[84]*tmp1[62]+z[4]*tmp1[63]+z[21]*tmp1[64]+z[38]*tmp1[65]+z[55]*tmp1[66]+z[72]*tmp1[67]+z[89]*tmp1[68]+z[9]*tmp1[69]+z[26]*tmp1[70]
            +z[43]*tmp1[71]+z[60]*tmp1[72]+z[77]*tmp1[73]+z[94]*tmp1[74]+z[14]*tmp1[75]+z[31]*tmp1[76]+z[48]*tmp1[77]+z[65]*tmp1[78]+z[82]*tmp1[79]+z[2]*tmp1[80]
            +z[19]*tmp1[81]+z[36]*tmp1[82]+z[53]*tmp1[83]+z[70]*tmp1[84]+z[87]*tmp1[85]+z[7]*tmp1[86]+z[24]*tmp1[87]+z[41]*tmp1[88]+z[58]*tmp1[89]+z[75]*tmp1[90]
            +z[92]*tmp1[91]+z[12]*tmp1[92]+z[29]*tmp1[93]+z[46]*tmp1[94]+z[63]*tmp1[95]+z[80]*tmp1[96];
            tab[i+18*stg_first]=tab[i+18*stg_first]+z[0]*tmp1[0]
            +z[18]*tmp1[1]+z[36]*tmp1[2]+z[54]*tmp1[3]+z[72]*tmp1[4]+z[90]*tmp1[5]+z[11]*tmp1[6]+z[29]*tmp1[7]+z[47]*tmp1[8]+z[65]*tmp1[9]+z[83]*tmp1[10]
            +z[4]*tmp1[11]+z[22]*tmp1[12]+z[40]*tmp1[13]+z[58]*tmp1[14]+z[76]*tmp1[15]+z[94]*tmp1[16]+z[15]*tmp1[17]+z[33]*tmp1[18]+z[51]*tmp1[19]+z[69]*tmp1[20]
            +z[87]*tmp1[21]+z[8]*tmp1[22]+z[26]*tmp1[23]+z[44]*tmp1[24]+z[62]*tmp1[25]+z[80]*tmp1[26]+z[1]*tmp1[27]+z[19]*tmp1[28]+z[37]*tmp1[29]+z[55]*tmp1[30]
            +z[73]*tmp1[31]+z[91]*tmp1[32]+z[12]*tmp1[33]+z[30]*tmp1[34]+z[48]*tmp1[35]+z[66]*tmp1[36]+z[84]*tmp1[37]+z[5]*tmp1[38]+z[23]*tmp1[39]+z[41]*tmp1[40]
            +z[59]*tmp1[41]+z[77]*tmp1[42]+z[95]*tmp1[43]+z[16]*tmp1[44]+z[34]*tmp1[45]+z[52]*tmp1[46]+z[70]*tmp1[47]+z[88]*tmp1[48]+z[9]*tmp1[49]+z[27]*tmp1[50]
            +z[45]*tmp1[51]+z[63]*tmp1[52]+z[81]*tmp1[53]+z[2]*tmp1[54]+z[20]*tmp1[55]+z[38]*tmp1[56]+z[56]*tmp1[57]+z[74]*tmp1[58]+z[92]*tmp1[59]+z[13]*tmp1[60]
            +z[31]*tmp1[61]+z[49]*tmp1[62]+z[67]*tmp1[63]+z[85]*tmp1[64]+z[6]*tmp1[65]+z[24]*tmp1[66]+z[42]*tmp1[67]+z[60]*tmp1[68]+z[78]*tmp1[69]+z[96]*tmp1[70]
            +z[17]*tmp1[71]+z[35]*tmp1[72]+z[53]*tmp1[73]+z[71]*tmp1[74]+z[89]*tmp1[75]+z[10]*tmp1[76]+z[28]*tmp1[77]+z[46]*tmp1[78]+z[64]*tmp1[79]+z[82]*tmp1[80]
            +z[3]*tmp1[81]+z[21]*tmp1[82]+z[39]*tmp1[83]+z[57]*tmp1[84]+z[75]*tmp1[85]+z[93]*tmp1[86]+z[14]*tmp1[87]+z[32]*tmp1[88]+z[50]*tmp1[89]+z[68]*tmp1[90]
            +z[86]*tmp1[91]+z[7]*tmp1[92]+z[25]*tmp1[93]+z[43]*tmp1[94]+z[61]*tmp1[95]+z[79]*tmp1[96];
            tab[i+19*stg_first]=tab[i+19*stg_first]+z[0]*tmp1[0]
            +z[19]*tmp1[1]+z[38]*tmp1[2]+z[57]*tmp1[3]+z[76]*tmp1[4]+z[95]*tmp1[5]+z[17]*tmp1[6]+z[36]*tmp1[7]+z[55]*tmp1[8]+z[74]*tmp1[9]+z[93]*tmp1[10]
            +z[15]*tmp1[11]+z[34]*tmp1[12]+z[53]*tmp1[13]+z[72]*tmp1[14]+z[91]*tmp1[15]+z[13]*tmp1[16]+z[32]*tmp1[17]+z[51]*tmp1[18]+z[70]*tmp1[19]+z[89]*tmp1[20]
            +z[11]*tmp1[21]+z[30]*tmp1[22]+z[49]*tmp1[23]+z[68]*tmp1[24]+z[87]*tmp1[25]+z[9]*tmp1[26]+z[28]*tmp1[27]+z[47]*tmp1[28]+z[66]*tmp1[29]+z[85]*tmp1[30]
            +z[7]*tmp1[31]+z[26]*tmp1[32]+z[45]*tmp1[33]+z[64]*tmp1[34]+z[83]*tmp1[35]+z[5]*tmp1[36]+z[24]*tmp1[37]+z[43]*tmp1[38]+z[62]*tmp1[39]+z[81]*tmp1[40]
            +z[3]*tmp1[41]+z[22]*tmp1[42]+z[41]*tmp1[43]+z[60]*tmp1[44]+z[79]*tmp1[45]+z[1]*tmp1[46]+z[20]*tmp1[47]+z[39]*tmp1[48]+z[58]*tmp1[49]+z[77]*tmp1[50]
            +z[96]*tmp1[51]+z[18]*tmp1[52]+z[37]*tmp1[53]+z[56]*tmp1[54]+z[75]*tmp1[55]+z[94]*tmp1[56]+z[16]*tmp1[57]+z[35]*tmp1[58]+z[54]*tmp1[59]+z[73]*tmp1[60]
            +z[92]*tmp1[61]+z[14]*tmp1[62]+z[33]*tmp1[63]+z[52]*tmp1[64]+z[71]*tmp1[65]+z[90]*tmp1[66]+z[12]*tmp1[67]+z[31]*tmp1[68]+z[50]*tmp1[69]+z[69]*tmp1[70]
            +z[88]*tmp1[71]+z[10]*tmp1[72]+z[29]*tmp1[73]+z[48]*tmp1[74]+z[67]*tmp1[75]+z[86]*tmp1[76]+z[8]*tmp1[77]+z[27]*tmp1[78]+z[46]*tmp1[79]+z[65]*tmp1[80]
            +z[84]*tmp1[81]+z[6]*tmp1[82]+z[25]*tmp1[83]+z[44]*tmp1[84]+z[63]*tmp1[85]+z[82]*tmp1[86]+z[4]*tmp1[87]+z[23]*tmp1[88]+z[42]*tmp1[89]+z[61]*tmp1[90]
            +z[80]*tmp1[91]+z[2]*tmp1[92]+z[21]*tmp1[93]+z[40]*tmp1[94]+z[59]*tmp1[95]+z[78]*tmp1[96];
            tab[i+20*stg_first]=tab[i+20*stg_first]+z[0]*tmp1[0]
            +z[20]*tmp1[1]+z[40]*tmp1[2]+z[60]*tmp1[3]+z[80]*tmp1[4]+z[3]*tmp1[5]+z[23]*tmp1[6]+z[43]*tmp1[7]+z[63]*tmp1[8]+z[83]*tmp1[9]+z[6]*tmp1[10]
            +z[26]*tmp1[11]+z[46]*tmp1[12]+z[66]*tmp1[13]+z[86]*tmp1[14]+z[9]*tmp1[15]+z[29]*tmp1[16]+z[49]*tmp1[17]+z[69]*tmp1[18]+z[89]*tmp1[19]+z[12]*tmp1[20]
            +z[32]*tmp1[21]+z[52]*tmp1[22]+z[72]*tmp1[23]+z[92]*tmp1[24]+z[15]*tmp1[25]+z[35]*tmp1[26]+z[55]*tmp1[27]+z[75]*tmp1[28]+z[95]*tmp1[29]+z[18]*tmp1[30]
            +z[38]*tmp1[31]+z[58]*tmp1[32]+z[78]*tmp1[33]+z[1]*tmp1[34]+z[21]*tmp1[35]+z[41]*tmp1[36]+z[61]*tmp1[37]+z[81]*tmp1[38]+z[4]*tmp1[39]+z[24]*tmp1[40]
            +z[44]*tmp1[41]+z[64]*tmp1[42]+z[84]*tmp1[43]+z[7]*tmp1[44]+z[27]*tmp1[45]+z[47]*tmp1[46]+z[67]*tmp1[47]+z[87]*tmp1[48]+z[10]*tmp1[49]+z[30]*tmp1[50]
            +z[50]*tmp1[51]+z[70]*tmp1[52]+z[90]*tmp1[53]+z[13]*tmp1[54]+z[33]*tmp1[55]+z[53]*tmp1[56]+z[73]*tmp1[57]+z[93]*tmp1[58]+z[16]*tmp1[59]+z[36]*tmp1[60]
            +z[56]*tmp1[61]+z[76]*tmp1[62]+z[96]*tmp1[63]+z[19]*tmp1[64]+z[39]*tmp1[65]+z[59]*tmp1[66]+z[79]*tmp1[67]+z[2]*tmp1[68]+z[22]*tmp1[69]+z[42]*tmp1[70]
            +z[62]*tmp1[71]+z[82]*tmp1[72]+z[5]*tmp1[73]+z[25]*tmp1[74]+z[45]*tmp1[75]+z[65]*tmp1[76]+z[85]*tmp1[77]+z[8]*tmp1[78]+z[28]*tmp1[79]+z[48]*tmp1[80]
            +z[68]*tmp1[81]+z[88]*tmp1[82]+z[11]*tmp1[83]+z[31]*tmp1[84]+z[51]*tmp1[85]+z[71]*tmp1[86]+z[91]*tmp1[87]+z[14]*tmp1[88]+z[34]*tmp1[89]+z[54]*tmp1[90]
            +z[74]*tmp1[91]+z[94]*tmp1[92]+z[17]*tmp1[93]+z[37]*tmp1[94]+z[57]*tmp1[95]+z[77]*tmp1[96];
            tab[i+21*stg_first]=tab[i+21*stg_first]+z[0]*tmp1[0]
            +z[21]*tmp1[1]+z[42]*tmp1[2]+z[63]*tmp1[3]+z[84]*tmp1[4]+z[8]*tmp1[5]+z[29]*tmp1[6]+z[50]*tmp1[7]+z[71]*tmp1[8]+z[92]*tmp1[9]+z[16]*tmp1[10]
            +z[37]*tmp1[11]+z[58]*tmp1[12]+z[79]*tmp1[13]+z[3]*tmp1[14]+z[24]*tmp1[15]+z[45]*tmp1[16]+z[66]*tmp1[17]+z[87]*tmp1[18]+z[11]*tmp1[19]+z[32]*tmp1[20]
            +z[53]*tmp1[21]+z[74]*tmp1[22]+z[95]*tmp1[23]+z[19]*tmp1[24]+z[40]*tmp1[25]+z[61]*tmp1[26]+z[82]*tmp1[27]+z[6]*tmp1[28]+z[27]*tmp1[29]+z[48]*tmp1[30]
            +z[69]*tmp1[31]+z[90]*tmp1[32]+z[14]*tmp1[33]+z[35]*tmp1[34]+z[56]*tmp1[35]+z[77]*tmp1[36]+z[1]*tmp1[37]+z[22]*tmp1[38]+z[43]*tmp1[39]+z[64]*tmp1[40]
            +z[85]*tmp1[41]+z[9]*tmp1[42]+z[30]*tmp1[43]+z[51]*tmp1[44]+z[72]*tmp1[45]+z[93]*tmp1[46]+z[17]*tmp1[47]+z[38]*tmp1[48]+z[59]*tmp1[49]+z[80]*tmp1[50]
            +z[4]*tmp1[51]+z[25]*tmp1[52]+z[46]*tmp1[53]+z[67]*tmp1[54]+z[88]*tmp1[55]+z[12]*tmp1[56]+z[33]*tmp1[57]+z[54]*tmp1[58]+z[75]*tmp1[59]+z[96]*tmp1[60]
            +z[20]*tmp1[61]+z[41]*tmp1[62]+z[62]*tmp1[63]+z[83]*tmp1[64]+z[7]*tmp1[65]+z[28]*tmp1[66]+z[49]*tmp1[67]+z[70]*tmp1[68]+z[91]*tmp1[69]+z[15]*tmp1[70]
            +z[36]*tmp1[71]+z[57]*tmp1[72]+z[78]*tmp1[73]+z[2]*tmp1[74]+z[23]*tmp1[75]+z[44]*tmp1[76]+z[65]*tmp1[77]+z[86]*tmp1[78]+z[10]*tmp1[79]+z[31]*tmp1[80]
            +z[52]*tmp1[81]+z[73]*tmp1[82]+z[94]*tmp1[83]+z[18]*tmp1[84]+z[39]*tmp1[85]+z[60]*tmp1[86]+z[81]*tmp1[87]+z[5]*tmp1[88]+z[26]*tmp1[89]+z[47]*tmp1[90]
            +z[68]*tmp1[91]+z[89]*tmp1[92]+z[13]*tmp1[93]+z[34]*tmp1[94]+z[55]*tmp1[95]+z[76]*tmp1[96];
            tab[i+22*stg_first]=tab[i+22*stg_first]+z[0]*tmp1[0]
            +z[22]*tmp1[1]+z[44]*tmp1[2]+z[66]*tmp1[3]+z[88]*tmp1[4]+z[13]*tmp1[5]+z[35]*tmp1[6]+z[57]*tmp1[7]+z[79]*tmp1[8]+z[4]*tmp1[9]+z[26]*tmp1[10]
            +z[48]*tmp1[11]+z[70]*tmp1[12]+z[92]*tmp1[13]+z[17]*tmp1[14]+z[39]*tmp1[15]+z[61]*tmp1[16]+z[83]*tmp1[17]+z[8]*tmp1[18]+z[30]*tmp1[19]+z[52]*tmp1[20]
            +z[74]*tmp1[21]+z[96]*tmp1[22]+z[21]*tmp1[23]+z[43]*tmp1[24]+z[65]*tmp1[25]+z[87]*tmp1[26]+z[12]*tmp1[27]+z[34]*tmp1[28]+z[56]*tmp1[29]+z[78]*tmp1[30]
            +z[3]*tmp1[31]+z[25]*tmp1[32]+z[47]*tmp1[33]+z[69]*tmp1[34]+z[91]*tmp1[35]+z[16]*tmp1[36]+z[38]*tmp1[37]+z[60]*tmp1[38]+z[82]*tmp1[39]+z[7]*tmp1[40]
            +z[29]*tmp1[41]+z[51]*tmp1[42]+z[73]*tmp1[43]+z[95]*tmp1[44]+z[20]*tmp1[45]+z[42]*tmp1[46]+z[64]*tmp1[47]+z[86]*tmp1[48]+z[11]*tmp1[49]+z[33]*tmp1[50]
            +z[55]*tmp1[51]+z[77]*tmp1[52]+z[2]*tmp1[53]+z[24]*tmp1[54]+z[46]*tmp1[55]+z[68]*tmp1[56]+z[90]*tmp1[57]+z[15]*tmp1[58]+z[37]*tmp1[59]+z[59]*tmp1[60]
            +z[81]*tmp1[61]+z[6]*tmp1[62]+z[28]*tmp1[63]+z[50]*tmp1[64]+z[72]*tmp1[65]+z[94]*tmp1[66]+z[19]*tmp1[67]+z[41]*tmp1[68]+z[63]*tmp1[69]+z[85]*tmp1[70]
            +z[10]*tmp1[71]+z[32]*tmp1[72]+z[54]*tmp1[73]+z[76]*tmp1[74]+z[1]*tmp1[75]+z[23]*tmp1[76]+z[45]*tmp1[77]+z[67]*tmp1[78]+z[89]*tmp1[79]+z[14]*tmp1[80]
            +z[36]*tmp1[81]+z[58]*tmp1[82]+z[80]*tmp1[83]+z[5]*tmp1[84]+z[27]*tmp1[85]+z[49]*tmp1[86]+z[71]*tmp1[87]+z[93]*tmp1[88]+z[18]*tmp1[89]+z[40]*tmp1[90]
            +z[62]*tmp1[91]+z[84]*tmp1[92]+z[9]*tmp1[93]+z[31]*tmp1[94]+z[53]*tmp1[95]+z[75]*tmp1[96];
            tab[i+23*stg_first]=tab[i+23*stg_first]+z[0]*tmp1[0]
            +z[23]*tmp1[1]+z[46]*tmp1[2]+z[69]*tmp1[3]+z[92]*tmp1[4]+z[18]*tmp1[5]+z[41]*tmp1[6]+z[64]*tmp1[7]+z[87]*tmp1[8]+z[13]*tmp1[9]+z[36]*tmp1[10]
            +z[59]*tmp1[11]+z[82]*tmp1[12]+z[8]*tmp1[13]+z[31]*tmp1[14]+z[54]*tmp1[15]+z[77]*tmp1[16]+z[3]*tmp1[17]+z[26]*tmp1[18]+z[49]*tmp1[19]+z[72]*tmp1[20]
            +z[95]*tmp1[21]+z[21]*tmp1[22]+z[44]*tmp1[23]+z[67]*tmp1[24]+z[90]*tmp1[25]+z[16]*tmp1[26]+z[39]*tmp1[27]+z[62]*tmp1[28]+z[85]*tmp1[29]+z[11]*tmp1[30]
            +z[34]*tmp1[31]+z[57]*tmp1[32]+z[80]*tmp1[33]+z[6]*tmp1[34]+z[29]*tmp1[35]+z[52]*tmp1[36]+z[75]*tmp1[37]+z[1]*tmp1[38]+z[24]*tmp1[39]+z[47]*tmp1[40]
            +z[70]*tmp1[41]+z[93]*tmp1[42]+z[19]*tmp1[43]+z[42]*tmp1[44]+z[65]*tmp1[45]+z[88]*tmp1[46]+z[14]*tmp1[47]+z[37]*tmp1[48]+z[60]*tmp1[49]+z[83]*tmp1[50]
            +z[9]*tmp1[51]+z[32]*tmp1[52]+z[55]*tmp1[53]+z[78]*tmp1[54]+z[4]*tmp1[55]+z[27]*tmp1[56]+z[50]*tmp1[57]+z[73]*tmp1[58]+z[96]*tmp1[59]+z[22]*tmp1[60]
            +z[45]*tmp1[61]+z[68]*tmp1[62]+z[91]*tmp1[63]+z[17]*tmp1[64]+z[40]*tmp1[65]+z[63]*tmp1[66]+z[86]*tmp1[67]+z[12]*tmp1[68]+z[35]*tmp1[69]+z[58]*tmp1[70]
            +z[81]*tmp1[71]+z[7]*tmp1[72]+z[30]*tmp1[73]+z[53]*tmp1[74]+z[76]*tmp1[75]+z[2]*tmp1[76]+z[25]*tmp1[77]+z[48]*tmp1[78]+z[71]*tmp1[79]+z[94]*tmp1[80]
            +z[20]*tmp1[81]+z[43]*tmp1[82]+z[66]*tmp1[83]+z[89]*tmp1[84]+z[15]*tmp1[85]+z[38]*tmp1[86]+z[61]*tmp1[87]+z[84]*tmp1[88]+z[10]*tmp1[89]+z[33]*tmp1[90]
            +z[56]*tmp1[91]+z[79]*tmp1[92]+z[5]*tmp1[93]+z[28]*tmp1[94]+z[51]*tmp1[95]+z[74]*tmp1[96];
            tab[i+24*stg_first]=tab[i+24*stg_first]+z[0]*tmp1[0]
            +z[24]*tmp1[1]+z[48]*tmp1[2]+z[72]*tmp1[3]+z[96]*tmp1[4]+z[23]*tmp1[5]+z[47]*tmp1[6]+z[71]*tmp1[7]+z[95]*tmp1[8]+z[22]*tmp1[9]+z[46]*tmp1[10]
            +z[70]*tmp1[11]+z[94]*tmp1[12]+z[21]*tmp1[13]+z[45]*tmp1[14]+z[69]*tmp1[15]+z[93]*tmp1[16]+z[20]*tmp1[17]+z[44]*tmp1[18]+z[68]*tmp1[19]+z[92]*tmp1[20]
            +z[19]*tmp1[21]+z[43]*tmp1[22]+z[67]*tmp1[23]+z[91]*tmp1[24]+z[18]*tmp1[25]+z[42]*tmp1[26]+z[66]*tmp1[27]+z[90]*tmp1[28]+z[17]*tmp1[29]+z[41]*tmp1[30]
            +z[65]*tmp1[31]+z[89]*tmp1[32]+z[16]*tmp1[33]+z[40]*tmp1[34]+z[64]*tmp1[35]+z[88]*tmp1[36]+z[15]*tmp1[37]+z[39]*tmp1[38]+z[63]*tmp1[39]+z[87]*tmp1[40]
            +z[14]*tmp1[41]+z[38]*tmp1[42]+z[62]*tmp1[43]+z[86]*tmp1[44]+z[13]*tmp1[45]+z[37]*tmp1[46]+z[61]*tmp1[47]+z[85]*tmp1[48]+z[12]*tmp1[49]+z[36]*tmp1[50]
            +z[60]*tmp1[51]+z[84]*tmp1[52]+z[11]*tmp1[53]+z[35]*tmp1[54]+z[59]*tmp1[55]+z[83]*tmp1[56]+z[10]*tmp1[57]+z[34]*tmp1[58]+z[58]*tmp1[59]+z[82]*tmp1[60]
            +z[9]*tmp1[61]+z[33]*tmp1[62]+z[57]*tmp1[63]+z[81]*tmp1[64]+z[8]*tmp1[65]+z[32]*tmp1[66]+z[56]*tmp1[67]+z[80]*tmp1[68]+z[7]*tmp1[69]+z[31]*tmp1[70]
            +z[55]*tmp1[71]+z[79]*tmp1[72]+z[6]*tmp1[73]+z[30]*tmp1[74]+z[54]*tmp1[75]+z[78]*tmp1[76]+z[5]*tmp1[77]+z[29]*tmp1[78]+z[53]*tmp1[79]+z[77]*tmp1[80]
            +z[4]*tmp1[81]+z[28]*tmp1[82]+z[52]*tmp1[83]+z[76]*tmp1[84]+z[3]*tmp1[85]+z[27]*tmp1[86]+z[51]*tmp1[87]+z[75]*tmp1[88]+z[2]*tmp1[89]+z[26]*tmp1[90]
            +z[50]*tmp1[91]+z[74]*tmp1[92]+z[1]*tmp1[93]+z[25]*tmp1[94]+z[49]*tmp1[95]+z[73]*tmp1[96];
            tab[i+25*stg_first]=tab[i+25*stg_first]+z[0]*tmp1[0]
            +z[25]*tmp1[1]+z[50]*tmp1[2]+z[75]*tmp1[3]+z[3]*tmp1[4]+z[28]*tmp1[5]+z[53]*tmp1[6]+z[78]*tmp1[7]+z[6]*tmp1[8]+z[31]*tmp1[9]+z[56]*tmp1[10]
            +z[81]*tmp1[11]+z[9]*tmp1[12]+z[34]*tmp1[13]+z[59]*tmp1[14]+z[84]*tmp1[15]+z[12]*tmp1[16]+z[37]*tmp1[17]+z[62]*tmp1[18]+z[87]*tmp1[19]+z[15]*tmp1[20]
            +z[40]*tmp1[21]+z[65]*tmp1[22]+z[90]*tmp1[23]+z[18]*tmp1[24]+z[43]*tmp1[25]+z[68]*tmp1[26]+z[93]*tmp1[27]+z[21]*tmp1[28]+z[46]*tmp1[29]+z[71]*tmp1[30]
            +z[96]*tmp1[31]+z[24]*tmp1[32]+z[49]*tmp1[33]+z[74]*tmp1[34]+z[2]*tmp1[35]+z[27]*tmp1[36]+z[52]*tmp1[37]+z[77]*tmp1[38]+z[5]*tmp1[39]+z[30]*tmp1[40]
            +z[55]*tmp1[41]+z[80]*tmp1[42]+z[8]*tmp1[43]+z[33]*tmp1[44]+z[58]*tmp1[45]+z[83]*tmp1[46]+z[11]*tmp1[47]+z[36]*tmp1[48]+z[61]*tmp1[49]+z[86]*tmp1[50]
            +z[14]*tmp1[51]+z[39]*tmp1[52]+z[64]*tmp1[53]+z[89]*tmp1[54]+z[17]*tmp1[55]+z[42]*tmp1[56]+z[67]*tmp1[57]+z[92]*tmp1[58]+z[20]*tmp1[59]+z[45]*tmp1[60]
            +z[70]*tmp1[61]+z[95]*tmp1[62]+z[23]*tmp1[63]+z[48]*tmp1[64]+z[73]*tmp1[65]+z[1]*tmp1[66]+z[26]*tmp1[67]+z[51]*tmp1[68]+z[76]*tmp1[69]+z[4]*tmp1[70]
            +z[29]*tmp1[71]+z[54]*tmp1[72]+z[79]*tmp1[73]+z[7]*tmp1[74]+z[32]*tmp1[75]+z[57]*tmp1[76]+z[82]*tmp1[77]+z[10]*tmp1[78]+z[35]*tmp1[79]+z[60]*tmp1[80]
            +z[85]*tmp1[81]+z[13]*tmp1[82]+z[38]*tmp1[83]+z[63]*tmp1[84]+z[88]*tmp1[85]+z[16]*tmp1[86]+z[41]*tmp1[87]+z[66]*tmp1[88]+z[91]*tmp1[89]+z[19]*tmp1[90]
            +z[44]*tmp1[91]+z[69]*tmp1[92]+z[94]*tmp1[93]+z[22]*tmp1[94]+z[47]*tmp1[95]+z[72]*tmp1[96];
            tab[i+26*stg_first]=tab[i+26*stg_first]+z[0]*tmp1[0]
            +z[26]*tmp1[1]+z[52]*tmp1[2]+z[78]*tmp1[3]+z[7]*tmp1[4]+z[33]*tmp1[5]+z[59]*tmp1[6]+z[85]*tmp1[7]+z[14]*tmp1[8]+z[40]*tmp1[9]+z[66]*tmp1[10]
            +z[92]*tmp1[11]+z[21]*tmp1[12]+z[47]*tmp1[13]+z[73]*tmp1[14]+z[2]*tmp1[15]+z[28]*tmp1[16]+z[54]*tmp1[17]+z[80]*tmp1[18]+z[9]*tmp1[19]+z[35]*tmp1[20]
            +z[61]*tmp1[21]+z[87]*tmp1[22]+z[16]*tmp1[23]+z[42]*tmp1[24]+z[68]*tmp1[25]+z[94]*tmp1[26]+z[23]*tmp1[27]+z[49]*tmp1[28]+z[75]*tmp1[29]+z[4]*tmp1[30]
            +z[30]*tmp1[31]+z[56]*tmp1[32]+z[82]*tmp1[33]+z[11]*tmp1[34]+z[37]*tmp1[35]+z[63]*tmp1[36]+z[89]*tmp1[37]+z[18]*tmp1[38]+z[44]*tmp1[39]+z[70]*tmp1[40]
            +z[96]*tmp1[41]+z[25]*tmp1[42]+z[51]*tmp1[43]+z[77]*tmp1[44]+z[6]*tmp1[45]+z[32]*tmp1[46]+z[58]*tmp1[47]+z[84]*tmp1[48]+z[13]*tmp1[49]+z[39]*tmp1[50]
            +z[65]*tmp1[51]+z[91]*tmp1[52]+z[20]*tmp1[53]+z[46]*tmp1[54]+z[72]*tmp1[55]+z[1]*tmp1[56]+z[27]*tmp1[57]+z[53]*tmp1[58]+z[79]*tmp1[59]+z[8]*tmp1[60]
            +z[34]*tmp1[61]+z[60]*tmp1[62]+z[86]*tmp1[63]+z[15]*tmp1[64]+z[41]*tmp1[65]+z[67]*tmp1[66]+z[93]*tmp1[67]+z[22]*tmp1[68]+z[48]*tmp1[69]+z[74]*tmp1[70]
            +z[3]*tmp1[71]+z[29]*tmp1[72]+z[55]*tmp1[73]+z[81]*tmp1[74]+z[10]*tmp1[75]+z[36]*tmp1[76]+z[62]*tmp1[77]+z[88]*tmp1[78]+z[17]*tmp1[79]+z[43]*tmp1[80]
            +z[69]*tmp1[81]+z[95]*tmp1[82]+z[24]*tmp1[83]+z[50]*tmp1[84]+z[76]*tmp1[85]+z[5]*tmp1[86]+z[31]*tmp1[87]+z[57]*tmp1[88]+z[83]*tmp1[89]+z[12]*tmp1[90]
            +z[38]*tmp1[91]+z[64]*tmp1[92]+z[90]*tmp1[93]+z[19]*tmp1[94]+z[45]*tmp1[95]+z[71]*tmp1[96];
            tab[i+27*stg_first]=tab[i+27*stg_first]+z[0]*tmp1[0]
            +z[27]*tmp1[1]+z[54]*tmp1[2]+z[81]*tmp1[3]+z[11]*tmp1[4]+z[38]*tmp1[5]+z[65]*tmp1[6]+z[92]*tmp1[7]+z[22]*tmp1[8]+z[49]*tmp1[9]+z[76]*tmp1[10]
            +z[6]*tmp1[11]+z[33]*tmp1[12]+z[60]*tmp1[13]+z[87]*tmp1[14]+z[17]*tmp1[15]+z[44]*tmp1[16]+z[71]*tmp1[17]+z[1]*tmp1[18]+z[28]*tmp1[19]+z[55]*tmp1[20]
            +z[82]*tmp1[21]+z[12]*tmp1[22]+z[39]*tmp1[23]+z[66]*tmp1[24]+z[93]*tmp1[25]+z[23]*tmp1[26]+z[50]*tmp1[27]+z[77]*tmp1[28]+z[7]*tmp1[29]+z[34]*tmp1[30]
            +z[61]*tmp1[31]+z[88]*tmp1[32]+z[18]*tmp1[33]+z[45]*tmp1[34]+z[72]*tmp1[35]+z[2]*tmp1[36]+z[29]*tmp1[37]+z[56]*tmp1[38]+z[83]*tmp1[39]+z[13]*tmp1[40]
            +z[40]*tmp1[41]+z[67]*tmp1[42]+z[94]*tmp1[43]+z[24]*tmp1[44]+z[51]*tmp1[45]+z[78]*tmp1[46]+z[8]*tmp1[47]+z[35]*tmp1[48]+z[62]*tmp1[49]+z[89]*tmp1[50]
            +z[19]*tmp1[51]+z[46]*tmp1[52]+z[73]*tmp1[53]+z[3]*tmp1[54]+z[30]*tmp1[55]+z[57]*tmp1[56]+z[84]*tmp1[57]+z[14]*tmp1[58]+z[41]*tmp1[59]+z[68]*tmp1[60]
            +z[95]*tmp1[61]+z[25]*tmp1[62]+z[52]*tmp1[63]+z[79]*tmp1[64]+z[9]*tmp1[65]+z[36]*tmp1[66]+z[63]*tmp1[67]+z[90]*tmp1[68]+z[20]*tmp1[69]+z[47]*tmp1[70]
            +z[74]*tmp1[71]+z[4]*tmp1[72]+z[31]*tmp1[73]+z[58]*tmp1[74]+z[85]*tmp1[75]+z[15]*tmp1[76]+z[42]*tmp1[77]+z[69]*tmp1[78]+z[96]*tmp1[79]+z[26]*tmp1[80]
            +z[53]*tmp1[81]+z[80]*tmp1[82]+z[10]*tmp1[83]+z[37]*tmp1[84]+z[64]*tmp1[85]+z[91]*tmp1[86]+z[21]*tmp1[87]+z[48]*tmp1[88]+z[75]*tmp1[89]+z[5]*tmp1[90]
            +z[32]*tmp1[91]+z[59]*tmp1[92]+z[86]*tmp1[93]+z[16]*tmp1[94]+z[43]*tmp1[95]+z[70]*tmp1[96];
            tab[i+28*stg_first]=tab[i+28*stg_first]+z[0]*tmp1[0]
            +z[28]*tmp1[1]+z[56]*tmp1[2]+z[84]*tmp1[3]+z[15]*tmp1[4]+z[43]*tmp1[5]+z[71]*tmp1[6]+z[2]*tmp1[7]+z[30]*tmp1[8]+z[58]*tmp1[9]+z[86]*tmp1[10]
            +z[17]*tmp1[11]+z[45]*tmp1[12]+z[73]*tmp1[13]+z[4]*tmp1[14]+z[32]*tmp1[15]+z[60]*tmp1[16]+z[88]*tmp1[17]+z[19]*tmp1[18]+z[47]*tmp1[19]+z[75]*tmp1[20]
            +z[6]*tmp1[21]+z[34]*tmp1[22]+z[62]*tmp1[23]+z[90]*tmp1[24]+z[21]*tmp1[25]+z[49]*tmp1[26]+z[77]*tmp1[27]+z[8]*tmp1[28]+z[36]*tmp1[29]+z[64]*tmp1[30]
            +z[92]*tmp1[31]+z[23]*tmp1[32]+z[51]*tmp1[33]+z[79]*tmp1[34]+z[10]*tmp1[35]+z[38]*tmp1[36]+z[66]*tmp1[37]+z[94]*tmp1[38]+z[25]*tmp1[39]+z[53]*tmp1[40]
            +z[81]*tmp1[41]+z[12]*tmp1[42]+z[40]*tmp1[43]+z[68]*tmp1[44]+z[96]*tmp1[45]+z[27]*tmp1[46]+z[55]*tmp1[47]+z[83]*tmp1[48]+z[14]*tmp1[49]+z[42]*tmp1[50]
            +z[70]*tmp1[51]+z[1]*tmp1[52]+z[29]*tmp1[53]+z[57]*tmp1[54]+z[85]*tmp1[55]+z[16]*tmp1[56]+z[44]*tmp1[57]+z[72]*tmp1[58]+z[3]*tmp1[59]+z[31]*tmp1[60]
            +z[59]*tmp1[61]+z[87]*tmp1[62]+z[18]*tmp1[63]+z[46]*tmp1[64]+z[74]*tmp1[65]+z[5]*tmp1[66]+z[33]*tmp1[67]+z[61]*tmp1[68]+z[89]*tmp1[69]+z[20]*tmp1[70]
            +z[48]*tmp1[71]+z[76]*tmp1[72]+z[7]*tmp1[73]+z[35]*tmp1[74]+z[63]*tmp1[75]+z[91]*tmp1[76]+z[22]*tmp1[77]+z[50]*tmp1[78]+z[78]*tmp1[79]+z[9]*tmp1[80]
            +z[37]*tmp1[81]+z[65]*tmp1[82]+z[93]*tmp1[83]+z[24]*tmp1[84]+z[52]*tmp1[85]+z[80]*tmp1[86]+z[11]*tmp1[87]+z[39]*tmp1[88]+z[67]*tmp1[89]+z[95]*tmp1[90]
            +z[26]*tmp1[91]+z[54]*tmp1[92]+z[82]*tmp1[93]+z[13]*tmp1[94]+z[41]*tmp1[95]+z[69]*tmp1[96];
            tab[i+29*stg_first]=tab[i+29*stg_first]+z[0]*tmp1[0]
            +z[29]*tmp1[1]+z[58]*tmp1[2]+z[87]*tmp1[3]+z[19]*tmp1[4]+z[48]*tmp1[5]+z[77]*tmp1[6]+z[9]*tmp1[7]+z[38]*tmp1[8]+z[67]*tmp1[9]+z[96]*tmp1[10]
            +z[28]*tmp1[11]+z[57]*tmp1[12]+z[86]*tmp1[13]+z[18]*tmp1[14]+z[47]*tmp1[15]+z[76]*tmp1[16]+z[8]*tmp1[17]+z[37]*tmp1[18]+z[66]*tmp1[19]+z[95]*tmp1[20]
            +z[27]*tmp1[21]+z[56]*tmp1[22]+z[85]*tmp1[23]+z[17]*tmp1[24]+z[46]*tmp1[25]+z[75]*tmp1[26]+z[7]*tmp1[27]+z[36]*tmp1[28]+z[65]*tmp1[29]+z[94]*tmp1[30]
            +z[26]*tmp1[31]+z[55]*tmp1[32]+z[84]*tmp1[33]+z[16]*tmp1[34]+z[45]*tmp1[35]+z[74]*tmp1[36]+z[6]*tmp1[37]+z[35]*tmp1[38]+z[64]*tmp1[39]+z[93]*tmp1[40]
            +z[25]*tmp1[41]+z[54]*tmp1[42]+z[83]*tmp1[43]+z[15]*tmp1[44]+z[44]*tmp1[45]+z[73]*tmp1[46]+z[5]*tmp1[47]+z[34]*tmp1[48]+z[63]*tmp1[49]+z[92]*tmp1[50]
            +z[24]*tmp1[51]+z[53]*tmp1[52]+z[82]*tmp1[53]+z[14]*tmp1[54]+z[43]*tmp1[55]+z[72]*tmp1[56]+z[4]*tmp1[57]+z[33]*tmp1[58]+z[62]*tmp1[59]+z[91]*tmp1[60]
            +z[23]*tmp1[61]+z[52]*tmp1[62]+z[81]*tmp1[63]+z[13]*tmp1[64]+z[42]*tmp1[65]+z[71]*tmp1[66]+z[3]*tmp1[67]+z[32]*tmp1[68]+z[61]*tmp1[69]+z[90]*tmp1[70]
            +z[22]*tmp1[71]+z[51]*tmp1[72]+z[80]*tmp1[73]+z[12]*tmp1[74]+z[41]*tmp1[75]+z[70]*tmp1[76]+z[2]*tmp1[77]+z[31]*tmp1[78]+z[60]*tmp1[79]+z[89]*tmp1[80]
            +z[21]*tmp1[81]+z[50]*tmp1[82]+z[79]*tmp1[83]+z[11]*tmp1[84]+z[40]*tmp1[85]+z[69]*tmp1[86]+z[1]*tmp1[87]+z[30]*tmp1[88]+z[59]*tmp1[89]+z[88]*tmp1[90]
            +z[20]*tmp1[91]+z[49]*tmp1[92]+z[78]*tmp1[93]+z[10]*tmp1[94]+z[39]*tmp1[95]+z[68]*tmp1[96];
            tab[i+30*stg_first]=tab[i+30*stg_first]+z[0]*tmp1[0]
            +z[30]*tmp1[1]+z[60]*tmp1[2]+z[90]*tmp1[3]+z[23]*tmp1[4]+z[53]*tmp1[5]+z[83]*tmp1[6]+z[16]*tmp1[7]+z[46]*tmp1[8]+z[76]*tmp1[9]+z[9]*tmp1[10]
            +z[39]*tmp1[11]+z[69]*tmp1[12]+z[2]*tmp1[13]+z[32]*tmp1[14]+z[62]*tmp1[15]+z[92]*tmp1[16]+z[25]*tmp1[17]+z[55]*tmp1[18]+z[85]*tmp1[19]+z[18]*tmp1[20]
            +z[48]*tmp1[21]+z[78]*tmp1[22]+z[11]*tmp1[23]+z[41]*tmp1[24]+z[71]*tmp1[25]+z[4]*tmp1[26]+z[34]*tmp1[27]+z[64]*tmp1[28]+z[94]*tmp1[29]+z[27]*tmp1[30]
            +z[57]*tmp1[31]+z[87]*tmp1[32]+z[20]*tmp1[33]+z[50]*tmp1[34]+z[80]*tmp1[35]+z[13]*tmp1[36]+z[43]*tmp1[37]+z[73]*tmp1[38]+z[6]*tmp1[39]+z[36]*tmp1[40]
            +z[66]*tmp1[41]+z[96]*tmp1[42]+z[29]*tmp1[43]+z[59]*tmp1[44]+z[89]*tmp1[45]+z[22]*tmp1[46]+z[52]*tmp1[47]+z[82]*tmp1[48]+z[15]*tmp1[49]+z[45]*tmp1[50]
            +z[75]*tmp1[51]+z[8]*tmp1[52]+z[38]*tmp1[53]+z[68]*tmp1[54]+z[1]*tmp1[55]+z[31]*tmp1[56]+z[61]*tmp1[57]+z[91]*tmp1[58]+z[24]*tmp1[59]+z[54]*tmp1[60]
            +z[84]*tmp1[61]+z[17]*tmp1[62]+z[47]*tmp1[63]+z[77]*tmp1[64]+z[10]*tmp1[65]+z[40]*tmp1[66]+z[70]*tmp1[67]+z[3]*tmp1[68]+z[33]*tmp1[69]+z[63]*tmp1[70]
            +z[93]*tmp1[71]+z[26]*tmp1[72]+z[56]*tmp1[73]+z[86]*tmp1[74]+z[19]*tmp1[75]+z[49]*tmp1[76]+z[79]*tmp1[77]+z[12]*tmp1[78]+z[42]*tmp1[79]+z[72]*tmp1[80]
            +z[5]*tmp1[81]+z[35]*tmp1[82]+z[65]*tmp1[83]+z[95]*tmp1[84]+z[28]*tmp1[85]+z[58]*tmp1[86]+z[88]*tmp1[87]+z[21]*tmp1[88]+z[51]*tmp1[89]+z[81]*tmp1[90]
            +z[14]*tmp1[91]+z[44]*tmp1[92]+z[74]*tmp1[93]+z[7]*tmp1[94]+z[37]*tmp1[95]+z[67]*tmp1[96];
            tab[i+31*stg_first]=tab[i+31*stg_first]+z[0]*tmp1[0]
            +z[31]*tmp1[1]+z[62]*tmp1[2]+z[93]*tmp1[3]+z[27]*tmp1[4]+z[58]*tmp1[5]+z[89]*tmp1[6]+z[23]*tmp1[7]+z[54]*tmp1[8]+z[85]*tmp1[9]+z[19]*tmp1[10]
            +z[50]*tmp1[11]+z[81]*tmp1[12]+z[15]*tmp1[13]+z[46]*tmp1[14]+z[77]*tmp1[15]+z[11]*tmp1[16]+z[42]*tmp1[17]+z[73]*tmp1[18]+z[7]*tmp1[19]+z[38]*tmp1[20]
            +z[69]*tmp1[21]+z[3]*tmp1[22]+z[34]*tmp1[23]+z[65]*tmp1[24]+z[96]*tmp1[25]+z[30]*tmp1[26]+z[61]*tmp1[27]+z[92]*tmp1[28]+z[26]*tmp1[29]+z[57]*tmp1[30]
            +z[88]*tmp1[31]+z[22]*tmp1[32]+z[53]*tmp1[33]+z[84]*tmp1[34]+z[18]*tmp1[35]+z[49]*tmp1[36]+z[80]*tmp1[37]+z[14]*tmp1[38]+z[45]*tmp1[39]+z[76]*tmp1[40]
            +z[10]*tmp1[41]+z[41]*tmp1[42]+z[72]*tmp1[43]+z[6]*tmp1[44]+z[37]*tmp1[45]+z[68]*tmp1[46]+z[2]*tmp1[47]+z[33]*tmp1[48]+z[64]*tmp1[49]+z[95]*tmp1[50]
            +z[29]*tmp1[51]+z[60]*tmp1[52]+z[91]*tmp1[53]+z[25]*tmp1[54]+z[56]*tmp1[55]+z[87]*tmp1[56]+z[21]*tmp1[57]+z[52]*tmp1[58]+z[83]*tmp1[59]+z[17]*tmp1[60]
            +z[48]*tmp1[61]+z[79]*tmp1[62]+z[13]*tmp1[63]+z[44]*tmp1[64]+z[75]*tmp1[65]+z[9]*tmp1[66]+z[40]*tmp1[67]+z[71]*tmp1[68]+z[5]*tmp1[69]+z[36]*tmp1[70]
            +z[67]*tmp1[71]+z[1]*tmp1[72]+z[32]*tmp1[73]+z[63]*tmp1[74]+z[94]*tmp1[75]+z[28]*tmp1[76]+z[59]*tmp1[77]+z[90]*tmp1[78]+z[24]*tmp1[79]+z[55]*tmp1[80]
            +z[86]*tmp1[81]+z[20]*tmp1[82]+z[51]*tmp1[83]+z[82]*tmp1[84]+z[16]*tmp1[85]+z[47]*tmp1[86]+z[78]*tmp1[87]+z[12]*tmp1[88]+z[43]*tmp1[89]+z[74]*tmp1[90]
            +z[8]*tmp1[91]+z[39]*tmp1[92]+z[70]*tmp1[93]+z[4]*tmp1[94]+z[35]*tmp1[95]+z[66]*tmp1[96];
            tab[i+32*stg_first]=tab[i+32*stg_first]+z[0]*tmp1[0]
            +z[32]*tmp1[1]+z[64]*tmp1[2]+z[96]*tmp1[3]+z[31]*tmp1[4]+z[63]*tmp1[5]+z[95]*tmp1[6]+z[30]*tmp1[7]+z[62]*tmp1[8]+z[94]*tmp1[9]+z[29]*tmp1[10]
            +z[61]*tmp1[11]+z[93]*tmp1[12]+z[28]*tmp1[13]+z[60]*tmp1[14]+z[92]*tmp1[15]+z[27]*tmp1[16]+z[59]*tmp1[17]+z[91]*tmp1[18]+z[26]*tmp1[19]+z[58]*tmp1[20]
            +z[90]*tmp1[21]+z[25]*tmp1[22]+z[57]*tmp1[23]+z[89]*tmp1[24]+z[24]*tmp1[25]+z[56]*tmp1[26]+z[88]*tmp1[27]+z[23]*tmp1[28]+z[55]*tmp1[29]+z[87]*tmp1[30]
            +z[22]*tmp1[31]+z[54]*tmp1[32]+z[86]*tmp1[33]+z[21]*tmp1[34]+z[53]*tmp1[35]+z[85]*tmp1[36]+z[20]*tmp1[37]+z[52]*tmp1[38]+z[84]*tmp1[39]+z[19]*tmp1[40]
            +z[51]*tmp1[41]+z[83]*tmp1[42]+z[18]*tmp1[43]+z[50]*tmp1[44]+z[82]*tmp1[45]+z[17]*tmp1[46]+z[49]*tmp1[47]+z[81]*tmp1[48]+z[16]*tmp1[49]+z[48]*tmp1[50]
            +z[80]*tmp1[51]+z[15]*tmp1[52]+z[47]*tmp1[53]+z[79]*tmp1[54]+z[14]*tmp1[55]+z[46]*tmp1[56]+z[78]*tmp1[57]+z[13]*tmp1[58]+z[45]*tmp1[59]+z[77]*tmp1[60]
            +z[12]*tmp1[61]+z[44]*tmp1[62]+z[76]*tmp1[63]+z[11]*tmp1[64]+z[43]*tmp1[65]+z[75]*tmp1[66]+z[10]*tmp1[67]+z[42]*tmp1[68]+z[74]*tmp1[69]+z[9]*tmp1[70]
            +z[41]*tmp1[71]+z[73]*tmp1[72]+z[8]*tmp1[73]+z[40]*tmp1[74]+z[72]*tmp1[75]+z[7]*tmp1[76]+z[39]*tmp1[77]+z[71]*tmp1[78]+z[6]*tmp1[79]+z[38]*tmp1[80]
            +z[70]*tmp1[81]+z[5]*tmp1[82]+z[37]*tmp1[83]+z[69]*tmp1[84]+z[4]*tmp1[85]+z[36]*tmp1[86]+z[68]*tmp1[87]+z[3]*tmp1[88]+z[35]*tmp1[89]+z[67]*tmp1[90]
            +z[2]*tmp1[91]+z[34]*tmp1[92]+z[66]*tmp1[93]+z[1]*tmp1[94]+z[33]*tmp1[95]+z[65]*tmp1[96];
            tab[i+33*stg_first]=tab[i+33*stg_first]+z[0]*tmp1[0]
            +z[33]*tmp1[1]+z[66]*tmp1[2]+z[2]*tmp1[3]+z[35]*tmp1[4]+z[68]*tmp1[5]+z[4]*tmp1[6]+z[37]*tmp1[7]+z[70]*tmp1[8]+z[6]*tmp1[9]+z[39]*tmp1[10]
            +z[72]*tmp1[11]+z[8]*tmp1[12]+z[41]*tmp1[13]+z[74]*tmp1[14]+z[10]*tmp1[15]+z[43]*tmp1[16]+z[76]*tmp1[17]+z[12]*tmp1[18]+z[45]*tmp1[19]+z[78]*tmp1[20]
            +z[14]*tmp1[21]+z[47]*tmp1[22]+z[80]*tmp1[23]+z[16]*tmp1[24]+z[49]*tmp1[25]+z[82]*tmp1[26]+z[18]*tmp1[27]+z[51]*tmp1[28]+z[84]*tmp1[29]+z[20]*tmp1[30]
            +z[53]*tmp1[31]+z[86]*tmp1[32]+z[22]*tmp1[33]+z[55]*tmp1[34]+z[88]*tmp1[35]+z[24]*tmp1[36]+z[57]*tmp1[37]+z[90]*tmp1[38]+z[26]*tmp1[39]+z[59]*tmp1[40]
            +z[92]*tmp1[41]+z[28]*tmp1[42]+z[61]*tmp1[43]+z[94]*tmp1[44]+z[30]*tmp1[45]+z[63]*tmp1[46]+z[96]*tmp1[47]+z[32]*tmp1[48]+z[65]*tmp1[49]+z[1]*tmp1[50]
            +z[34]*tmp1[51]+z[67]*tmp1[52]+z[3]*tmp1[53]+z[36]*tmp1[54]+z[69]*tmp1[55]+z[5]*tmp1[56]+z[38]*tmp1[57]+z[71]*tmp1[58]+z[7]*tmp1[59]+z[40]*tmp1[60]
            +z[73]*tmp1[61]+z[9]*tmp1[62]+z[42]*tmp1[63]+z[75]*tmp1[64]+z[11]*tmp1[65]+z[44]*tmp1[66]+z[77]*tmp1[67]+z[13]*tmp1[68]+z[46]*tmp1[69]+z[79]*tmp1[70]
            +z[15]*tmp1[71]+z[48]*tmp1[72]+z[81]*tmp1[73]+z[17]*tmp1[74]+z[50]*tmp1[75]+z[83]*tmp1[76]+z[19]*tmp1[77]+z[52]*tmp1[78]+z[85]*tmp1[79]+z[21]*tmp1[80]
            +z[54]*tmp1[81]+z[87]*tmp1[82]+z[23]*tmp1[83]+z[56]*tmp1[84]+z[89]*tmp1[85]+z[25]*tmp1[86]+z[58]*tmp1[87]+z[91]*tmp1[88]+z[27]*tmp1[89]+z[60]*tmp1[90]
            +z[93]*tmp1[91]+z[29]*tmp1[92]+z[62]*tmp1[93]+z[95]*tmp1[94]+z[31]*tmp1[95]+z[64]*tmp1[96];
            tab[i+34*stg_first]=tab[i+34*stg_first]+z[0]*tmp1[0]
            +z[34]*tmp1[1]+z[68]*tmp1[2]+z[5]*tmp1[3]+z[39]*tmp1[4]+z[73]*tmp1[5]+z[10]*tmp1[6]+z[44]*tmp1[7]+z[78]*tmp1[8]+z[15]*tmp1[9]+z[49]*tmp1[10]
            +z[83]*tmp1[11]+z[20]*tmp1[12]+z[54]*tmp1[13]+z[88]*tmp1[14]+z[25]*tmp1[15]+z[59]*tmp1[16]+z[93]*tmp1[17]+z[30]*tmp1[18]+z[64]*tmp1[19]+z[1]*tmp1[20]
            +z[35]*tmp1[21]+z[69]*tmp1[22]+z[6]*tmp1[23]+z[40]*tmp1[24]+z[74]*tmp1[25]+z[11]*tmp1[26]+z[45]*tmp1[27]+z[79]*tmp1[28]+z[16]*tmp1[29]+z[50]*tmp1[30]
            +z[84]*tmp1[31]+z[21]*tmp1[32]+z[55]*tmp1[33]+z[89]*tmp1[34]+z[26]*tmp1[35]+z[60]*tmp1[36]+z[94]*tmp1[37]+z[31]*tmp1[38]+z[65]*tmp1[39]+z[2]*tmp1[40]
            +z[36]*tmp1[41]+z[70]*tmp1[42]+z[7]*tmp1[43]+z[41]*tmp1[44]+z[75]*tmp1[45]+z[12]*tmp1[46]+z[46]*tmp1[47]+z[80]*tmp1[48]+z[17]*tmp1[49]+z[51]*tmp1[50]
            +z[85]*tmp1[51]+z[22]*tmp1[52]+z[56]*tmp1[53]+z[90]*tmp1[54]+z[27]*tmp1[55]+z[61]*tmp1[56]+z[95]*tmp1[57]+z[32]*tmp1[58]+z[66]*tmp1[59]+z[3]*tmp1[60]
            +z[37]*tmp1[61]+z[71]*tmp1[62]+z[8]*tmp1[63]+z[42]*tmp1[64]+z[76]*tmp1[65]+z[13]*tmp1[66]+z[47]*tmp1[67]+z[81]*tmp1[68]+z[18]*tmp1[69]+z[52]*tmp1[70]
            +z[86]*tmp1[71]+z[23]*tmp1[72]+z[57]*tmp1[73]+z[91]*tmp1[74]+z[28]*tmp1[75]+z[62]*tmp1[76]+z[96]*tmp1[77]+z[33]*tmp1[78]+z[67]*tmp1[79]+z[4]*tmp1[80]
            +z[38]*tmp1[81]+z[72]*tmp1[82]+z[9]*tmp1[83]+z[43]*tmp1[84]+z[77]*tmp1[85]+z[14]*tmp1[86]+z[48]*tmp1[87]+z[82]*tmp1[88]+z[19]*tmp1[89]+z[53]*tmp1[90]
            +z[87]*tmp1[91]+z[24]*tmp1[92]+z[58]*tmp1[93]+z[92]*tmp1[94]+z[29]*tmp1[95]+z[63]*tmp1[96];
            tab[i+35*stg_first]=tab[i+35*stg_first]+z[0]*tmp1[0]
            +z[35]*tmp1[1]+z[70]*tmp1[2]+z[8]*tmp1[3]+z[43]*tmp1[4]+z[78]*tmp1[5]+z[16]*tmp1[6]+z[51]*tmp1[7]+z[86]*tmp1[8]+z[24]*tmp1[9]+z[59]*tmp1[10]
            +z[94]*tmp1[11]+z[32]*tmp1[12]+z[67]*tmp1[13]+z[5]*tmp1[14]+z[40]*tmp1[15]+z[75]*tmp1[16]+z[13]*tmp1[17]+z[48]*tmp1[18]+z[83]*tmp1[19]+z[21]*tmp1[20]
            +z[56]*tmp1[21]+z[91]*tmp1[22]+z[29]*tmp1[23]+z[64]*tmp1[24]+z[2]*tmp1[25]+z[37]*tmp1[26]+z[72]*tmp1[27]+z[10]*tmp1[28]+z[45]*tmp1[29]+z[80]*tmp1[30]
            +z[18]*tmp1[31]+z[53]*tmp1[32]+z[88]*tmp1[33]+z[26]*tmp1[34]+z[61]*tmp1[35]+z[96]*tmp1[36]+z[34]*tmp1[37]+z[69]*tmp1[38]+z[7]*tmp1[39]+z[42]*tmp1[40]
            +z[77]*tmp1[41]+z[15]*tmp1[42]+z[50]*tmp1[43]+z[85]*tmp1[44]+z[23]*tmp1[45]+z[58]*tmp1[46]+z[93]*tmp1[47]+z[31]*tmp1[48]+z[66]*tmp1[49]+z[4]*tmp1[50]
            +z[39]*tmp1[51]+z[74]*tmp1[52]+z[12]*tmp1[53]+z[47]*tmp1[54]+z[82]*tmp1[55]+z[20]*tmp1[56]+z[55]*tmp1[57]+z[90]*tmp1[58]+z[28]*tmp1[59]+z[63]*tmp1[60]
            +z[1]*tmp1[61]+z[36]*tmp1[62]+z[71]*tmp1[63]+z[9]*tmp1[64]+z[44]*tmp1[65]+z[79]*tmp1[66]+z[17]*tmp1[67]+z[52]*tmp1[68]+z[87]*tmp1[69]+z[25]*tmp1[70]
            +z[60]*tmp1[71]+z[95]*tmp1[72]+z[33]*tmp1[73]+z[68]*tmp1[74]+z[6]*tmp1[75]+z[41]*tmp1[76]+z[76]*tmp1[77]+z[14]*tmp1[78]+z[49]*tmp1[79]+z[84]*tmp1[80]
            +z[22]*tmp1[81]+z[57]*tmp1[82]+z[92]*tmp1[83]+z[30]*tmp1[84]+z[65]*tmp1[85]+z[3]*tmp1[86]+z[38]*tmp1[87]+z[73]*tmp1[88]+z[11]*tmp1[89]+z[46]*tmp1[90]
            +z[81]*tmp1[91]+z[19]*tmp1[92]+z[54]*tmp1[93]+z[89]*tmp1[94]+z[27]*tmp1[95]+z[62]*tmp1[96];
            tab[i+36*stg_first]=tab[i+36*stg_first]+z[0]*tmp1[0]
            +z[36]*tmp1[1]+z[72]*tmp1[2]+z[11]*tmp1[3]+z[47]*tmp1[4]+z[83]*tmp1[5]+z[22]*tmp1[6]+z[58]*tmp1[7]+z[94]*tmp1[8]+z[33]*tmp1[9]+z[69]*tmp1[10]
            +z[8]*tmp1[11]+z[44]*tmp1[12]+z[80]*tmp1[13]+z[19]*tmp1[14]+z[55]*tmp1[15]+z[91]*tmp1[16]+z[30]*tmp1[17]+z[66]*tmp1[18]+z[5]*tmp1[19]+z[41]*tmp1[20]
            +z[77]*tmp1[21]+z[16]*tmp1[22]+z[52]*tmp1[23]+z[88]*tmp1[24]+z[27]*tmp1[25]+z[63]*tmp1[26]+z[2]*tmp1[27]+z[38]*tmp1[28]+z[74]*tmp1[29]+z[13]*tmp1[30]
            +z[49]*tmp1[31]+z[85]*tmp1[32]+z[24]*tmp1[33]+z[60]*tmp1[34]+z[96]*tmp1[35]+z[35]*tmp1[36]+z[71]*tmp1[37]+z[10]*tmp1[38]+z[46]*tmp1[39]+z[82]*tmp1[40]
            +z[21]*tmp1[41]+z[57]*tmp1[42]+z[93]*tmp1[43]+z[32]*tmp1[44]+z[68]*tmp1[45]+z[7]*tmp1[46]+z[43]*tmp1[47]+z[79]*tmp1[48]+z[18]*tmp1[49]+z[54]*tmp1[50]
            +z[90]*tmp1[51]+z[29]*tmp1[52]+z[65]*tmp1[53]+z[4]*tmp1[54]+z[40]*tmp1[55]+z[76]*tmp1[56]+z[15]*tmp1[57]+z[51]*tmp1[58]+z[87]*tmp1[59]+z[26]*tmp1[60]
            +z[62]*tmp1[61]+z[1]*tmp1[62]+z[37]*tmp1[63]+z[73]*tmp1[64]+z[12]*tmp1[65]+z[48]*tmp1[66]+z[84]*tmp1[67]+z[23]*tmp1[68]+z[59]*tmp1[69]+z[95]*tmp1[70]
            +z[34]*tmp1[71]+z[70]*tmp1[72]+z[9]*tmp1[73]+z[45]*tmp1[74]+z[81]*tmp1[75]+z[20]*tmp1[76]+z[56]*tmp1[77]+z[92]*tmp1[78]+z[31]*tmp1[79]+z[67]*tmp1[80]
            +z[6]*tmp1[81]+z[42]*tmp1[82]+z[78]*tmp1[83]+z[17]*tmp1[84]+z[53]*tmp1[85]+z[89]*tmp1[86]+z[28]*tmp1[87]+z[64]*tmp1[88]+z[3]*tmp1[89]+z[39]*tmp1[90]
            +z[75]*tmp1[91]+z[14]*tmp1[92]+z[50]*tmp1[93]+z[86]*tmp1[94]+z[25]*tmp1[95]+z[61]*tmp1[96];
            tab[i+37*stg_first]=tab[i+37*stg_first]+z[0]*tmp1[0]
            +z[37]*tmp1[1]+z[74]*tmp1[2]+z[14]*tmp1[3]+z[51]*tmp1[4]+z[88]*tmp1[5]+z[28]*tmp1[6]+z[65]*tmp1[7]+z[5]*tmp1[8]+z[42]*tmp1[9]+z[79]*tmp1[10]
            +z[19]*tmp1[11]+z[56]*tmp1[12]+z[93]*tmp1[13]+z[33]*tmp1[14]+z[70]*tmp1[15]+z[10]*tmp1[16]+z[47]*tmp1[17]+z[84]*tmp1[18]+z[24]*tmp1[19]+z[61]*tmp1[20]
            +z[1]*tmp1[21]+z[38]*tmp1[22]+z[75]*tmp1[23]+z[15]*tmp1[24]+z[52]*tmp1[25]+z[89]*tmp1[26]+z[29]*tmp1[27]+z[66]*tmp1[28]+z[6]*tmp1[29]+z[43]*tmp1[30]
            +z[80]*tmp1[31]+z[20]*tmp1[32]+z[57]*tmp1[33]+z[94]*tmp1[34]+z[34]*tmp1[35]+z[71]*tmp1[36]+z[11]*tmp1[37]+z[48]*tmp1[38]+z[85]*tmp1[39]+z[25]*tmp1[40]
            +z[62]*tmp1[41]+z[2]*tmp1[42]+z[39]*tmp1[43]+z[76]*tmp1[44]+z[16]*tmp1[45]+z[53]*tmp1[46]+z[90]*tmp1[47]+z[30]*tmp1[48]+z[67]*tmp1[49]+z[7]*tmp1[50]
            +z[44]*tmp1[51]+z[81]*tmp1[52]+z[21]*tmp1[53]+z[58]*tmp1[54]+z[95]*tmp1[55]+z[35]*tmp1[56]+z[72]*tmp1[57]+z[12]*tmp1[58]+z[49]*tmp1[59]+z[86]*tmp1[60]
            +z[26]*tmp1[61]+z[63]*tmp1[62]+z[3]*tmp1[63]+z[40]*tmp1[64]+z[77]*tmp1[65]+z[17]*tmp1[66]+z[54]*tmp1[67]+z[91]*tmp1[68]+z[31]*tmp1[69]+z[68]*tmp1[70]
            +z[8]*tmp1[71]+z[45]*tmp1[72]+z[82]*tmp1[73]+z[22]*tmp1[74]+z[59]*tmp1[75]+z[96]*tmp1[76]+z[36]*tmp1[77]+z[73]*tmp1[78]+z[13]*tmp1[79]+z[50]*tmp1[80]
            +z[87]*tmp1[81]+z[27]*tmp1[82]+z[64]*tmp1[83]+z[4]*tmp1[84]+z[41]*tmp1[85]+z[78]*tmp1[86]+z[18]*tmp1[87]+z[55]*tmp1[88]+z[92]*tmp1[89]+z[32]*tmp1[90]
            +z[69]*tmp1[91]+z[9]*tmp1[92]+z[46]*tmp1[93]+z[83]*tmp1[94]+z[23]*tmp1[95]+z[60]*tmp1[96];
            tab[i+38*stg_first]=tab[i+38*stg_first]+z[0]*tmp1[0]
            +z[38]*tmp1[1]+z[76]*tmp1[2]+z[17]*tmp1[3]+z[55]*tmp1[4]+z[93]*tmp1[5]+z[34]*tmp1[6]+z[72]*tmp1[7]+z[13]*tmp1[8]+z[51]*tmp1[9]+z[89]*tmp1[10]
            +z[30]*tmp1[11]+z[68]*tmp1[12]+z[9]*tmp1[13]+z[47]*tmp1[14]+z[85]*tmp1[15]+z[26]*tmp1[16]+z[64]*tmp1[17]+z[5]*tmp1[18]+z[43]*tmp1[19]+z[81]*tmp1[20]
            +z[22]*tmp1[21]+z[60]*tmp1[22]+z[1]*tmp1[23]+z[39]*tmp1[24]+z[77]*tmp1[25]+z[18]*tmp1[26]+z[56]*tmp1[27]+z[94]*tmp1[28]+z[35]*tmp1[29]+z[73]*tmp1[30]
            +z[14]*tmp1[31]+z[52]*tmp1[32]+z[90]*tmp1[33]+z[31]*tmp1[34]+z[69]*tmp1[35]+z[10]*tmp1[36]+z[48]*tmp1[37]+z[86]*tmp1[38]+z[27]*tmp1[39]+z[65]*tmp1[40]
            +z[6]*tmp1[41]+z[44]*tmp1[42]+z[82]*tmp1[43]+z[23]*tmp1[44]+z[61]*tmp1[45]+z[2]*tmp1[46]+z[40]*tmp1[47]+z[78]*tmp1[48]+z[19]*tmp1[49]+z[57]*tmp1[50]
            +z[95]*tmp1[51]+z[36]*tmp1[52]+z[74]*tmp1[53]+z[15]*tmp1[54]+z[53]*tmp1[55]+z[91]*tmp1[56]+z[32]*tmp1[57]+z[70]*tmp1[58]+z[11]*tmp1[59]+z[49]*tmp1[60]
            +z[87]*tmp1[61]+z[28]*tmp1[62]+z[66]*tmp1[63]+z[7]*tmp1[64]+z[45]*tmp1[65]+z[83]*tmp1[66]+z[24]*tmp1[67]+z[62]*tmp1[68]+z[3]*tmp1[69]+z[41]*tmp1[70]
            +z[79]*tmp1[71]+z[20]*tmp1[72]+z[58]*tmp1[73]+z[96]*tmp1[74]+z[37]*tmp1[75]+z[75]*tmp1[76]+z[16]*tmp1[77]+z[54]*tmp1[78]+z[92]*tmp1[79]+z[33]*tmp1[80]
            +z[71]*tmp1[81]+z[12]*tmp1[82]+z[50]*tmp1[83]+z[88]*tmp1[84]+z[29]*tmp1[85]+z[67]*tmp1[86]+z[8]*tmp1[87]+z[46]*tmp1[88]+z[84]*tmp1[89]+z[25]*tmp1[90]
            +z[63]*tmp1[91]+z[4]*tmp1[92]+z[42]*tmp1[93]+z[80]*tmp1[94]+z[21]*tmp1[95]+z[59]*tmp1[96];
            tab[i+39*stg_first]=tab[i+39*stg_first]+z[0]*tmp1[0]
            +z[39]*tmp1[1]+z[78]*tmp1[2]+z[20]*tmp1[3]+z[59]*tmp1[4]+z[1]*tmp1[5]+z[40]*tmp1[6]+z[79]*tmp1[7]+z[21]*tmp1[8]+z[60]*tmp1[9]+z[2]*tmp1[10]
            +z[41]*tmp1[11]+z[80]*tmp1[12]+z[22]*tmp1[13]+z[61]*tmp1[14]+z[3]*tmp1[15]+z[42]*tmp1[16]+z[81]*tmp1[17]+z[23]*tmp1[18]+z[62]*tmp1[19]+z[4]*tmp1[20]
            +z[43]*tmp1[21]+z[82]*tmp1[22]+z[24]*tmp1[23]+z[63]*tmp1[24]+z[5]*tmp1[25]+z[44]*tmp1[26]+z[83]*tmp1[27]+z[25]*tmp1[28]+z[64]*tmp1[29]+z[6]*tmp1[30]
            +z[45]*tmp1[31]+z[84]*tmp1[32]+z[26]*tmp1[33]+z[65]*tmp1[34]+z[7]*tmp1[35]+z[46]*tmp1[36]+z[85]*tmp1[37]+z[27]*tmp1[38]+z[66]*tmp1[39]+z[8]*tmp1[40]
            +z[47]*tmp1[41]+z[86]*tmp1[42]+z[28]*tmp1[43]+z[67]*tmp1[44]+z[9]*tmp1[45]+z[48]*tmp1[46]+z[87]*tmp1[47]+z[29]*tmp1[48]+z[68]*tmp1[49]+z[10]*tmp1[50]
            +z[49]*tmp1[51]+z[88]*tmp1[52]+z[30]*tmp1[53]+z[69]*tmp1[54]+z[11]*tmp1[55]+z[50]*tmp1[56]+z[89]*tmp1[57]+z[31]*tmp1[58]+z[70]*tmp1[59]+z[12]*tmp1[60]
            +z[51]*tmp1[61]+z[90]*tmp1[62]+z[32]*tmp1[63]+z[71]*tmp1[64]+z[13]*tmp1[65]+z[52]*tmp1[66]+z[91]*tmp1[67]+z[33]*tmp1[68]+z[72]*tmp1[69]+z[14]*tmp1[70]
            +z[53]*tmp1[71]+z[92]*tmp1[72]+z[34]*tmp1[73]+z[73]*tmp1[74]+z[15]*tmp1[75]+z[54]*tmp1[76]+z[93]*tmp1[77]+z[35]*tmp1[78]+z[74]*tmp1[79]+z[16]*tmp1[80]
            +z[55]*tmp1[81]+z[94]*tmp1[82]+z[36]*tmp1[83]+z[75]*tmp1[84]+z[17]*tmp1[85]+z[56]*tmp1[86]+z[95]*tmp1[87]+z[37]*tmp1[88]+z[76]*tmp1[89]+z[18]*tmp1[90]
            +z[57]*tmp1[91]+z[96]*tmp1[92]+z[38]*tmp1[93]+z[77]*tmp1[94]+z[19]*tmp1[95]+z[58]*tmp1[96];
            tab[i+40*stg_first]=tab[i+40*stg_first]+z[0]*tmp1[0]
            +z[40]*tmp1[1]+z[80]*tmp1[2]+z[23]*tmp1[3]+z[63]*tmp1[4]+z[6]*tmp1[5]+z[46]*tmp1[6]+z[86]*tmp1[7]+z[29]*tmp1[8]+z[69]*tmp1[9]+z[12]*tmp1[10]
            +z[52]*tmp1[11]+z[92]*tmp1[12]+z[35]*tmp1[13]+z[75]*tmp1[14]+z[18]*tmp1[15]+z[58]*tmp1[16]+z[1]*tmp1[17]+z[41]*tmp1[18]+z[81]*tmp1[19]+z[24]*tmp1[20]
            +z[64]*tmp1[21]+z[7]*tmp1[22]+z[47]*tmp1[23]+z[87]*tmp1[24]+z[30]*tmp1[25]+z[70]*tmp1[26]+z[13]*tmp1[27]+z[53]*tmp1[28]+z[93]*tmp1[29]+z[36]*tmp1[30]
            +z[76]*tmp1[31]+z[19]*tmp1[32]+z[59]*tmp1[33]+z[2]*tmp1[34]+z[42]*tmp1[35]+z[82]*tmp1[36]+z[25]*tmp1[37]+z[65]*tmp1[38]+z[8]*tmp1[39]+z[48]*tmp1[40]
            +z[88]*tmp1[41]+z[31]*tmp1[42]+z[71]*tmp1[43]+z[14]*tmp1[44]+z[54]*tmp1[45]+z[94]*tmp1[46]+z[37]*tmp1[47]+z[77]*tmp1[48]+z[20]*tmp1[49]+z[60]*tmp1[50]
            +z[3]*tmp1[51]+z[43]*tmp1[52]+z[83]*tmp1[53]+z[26]*tmp1[54]+z[66]*tmp1[55]+z[9]*tmp1[56]+z[49]*tmp1[57]+z[89]*tmp1[58]+z[32]*tmp1[59]+z[72]*tmp1[60]
            +z[15]*tmp1[61]+z[55]*tmp1[62]+z[95]*tmp1[63]+z[38]*tmp1[64]+z[78]*tmp1[65]+z[21]*tmp1[66]+z[61]*tmp1[67]+z[4]*tmp1[68]+z[44]*tmp1[69]+z[84]*tmp1[70]
            +z[27]*tmp1[71]+z[67]*tmp1[72]+z[10]*tmp1[73]+z[50]*tmp1[74]+z[90]*tmp1[75]+z[33]*tmp1[76]+z[73]*tmp1[77]+z[16]*tmp1[78]+z[56]*tmp1[79]+z[96]*tmp1[80]
            +z[39]*tmp1[81]+z[79]*tmp1[82]+z[22]*tmp1[83]+z[62]*tmp1[84]+z[5]*tmp1[85]+z[45]*tmp1[86]+z[85]*tmp1[87]+z[28]*tmp1[88]+z[68]*tmp1[89]+z[11]*tmp1[90]
            +z[51]*tmp1[91]+z[91]*tmp1[92]+z[34]*tmp1[93]+z[74]*tmp1[94]+z[17]*tmp1[95]+z[57]*tmp1[96];
            tab[i+41*stg_first]=tab[i+41*stg_first]+z[0]*tmp1[0]
            +z[41]*tmp1[1]+z[82]*tmp1[2]+z[26]*tmp1[3]+z[67]*tmp1[4]+z[11]*tmp1[5]+z[52]*tmp1[6]+z[93]*tmp1[7]+z[37]*tmp1[8]+z[78]*tmp1[9]+z[22]*tmp1[10]
            +z[63]*tmp1[11]+z[7]*tmp1[12]+z[48]*tmp1[13]+z[89]*tmp1[14]+z[33]*tmp1[15]+z[74]*tmp1[16]+z[18]*tmp1[17]+z[59]*tmp1[18]+z[3]*tmp1[19]+z[44]*tmp1[20]
            +z[85]*tmp1[21]+z[29]*tmp1[22]+z[70]*tmp1[23]+z[14]*tmp1[24]+z[55]*tmp1[25]+z[96]*tmp1[26]+z[40]*tmp1[27]+z[81]*tmp1[28]+z[25]*tmp1[29]+z[66]*tmp1[30]
            +z[10]*tmp1[31]+z[51]*tmp1[32]+z[92]*tmp1[33]+z[36]*tmp1[34]+z[77]*tmp1[35]+z[21]*tmp1[36]+z[62]*tmp1[37]+z[6]*tmp1[38]+z[47]*tmp1[39]+z[88]*tmp1[40]
            +z[32]*tmp1[41]+z[73]*tmp1[42]+z[17]*tmp1[43]+z[58]*tmp1[44]+z[2]*tmp1[45]+z[43]*tmp1[46]+z[84]*tmp1[47]+z[28]*tmp1[48]+z[69]*tmp1[49]+z[13]*tmp1[50]
            +z[54]*tmp1[51]+z[95]*tmp1[52]+z[39]*tmp1[53]+z[80]*tmp1[54]+z[24]*tmp1[55]+z[65]*tmp1[56]+z[9]*tmp1[57]+z[50]*tmp1[58]+z[91]*tmp1[59]+z[35]*tmp1[60]
            +z[76]*tmp1[61]+z[20]*tmp1[62]+z[61]*tmp1[63]+z[5]*tmp1[64]+z[46]*tmp1[65]+z[87]*tmp1[66]+z[31]*tmp1[67]+z[72]*tmp1[68]+z[16]*tmp1[69]+z[57]*tmp1[70]
            +z[1]*tmp1[71]+z[42]*tmp1[72]+z[83]*tmp1[73]+z[27]*tmp1[74]+z[68]*tmp1[75]+z[12]*tmp1[76]+z[53]*tmp1[77]+z[94]*tmp1[78]+z[38]*tmp1[79]+z[79]*tmp1[80]
            +z[23]*tmp1[81]+z[64]*tmp1[82]+z[8]*tmp1[83]+z[49]*tmp1[84]+z[90]*tmp1[85]+z[34]*tmp1[86]+z[75]*tmp1[87]+z[19]*tmp1[88]+z[60]*tmp1[89]+z[4]*tmp1[90]
            +z[45]*tmp1[91]+z[86]*tmp1[92]+z[30]*tmp1[93]+z[71]*tmp1[94]+z[15]*tmp1[95]+z[56]*tmp1[96];
            tab[i+42*stg_first]=tab[i+42*stg_first]+z[0]*tmp1[0]
            +z[42]*tmp1[1]+z[84]*tmp1[2]+z[29]*tmp1[3]+z[71]*tmp1[4]+z[16]*tmp1[5]+z[58]*tmp1[6]+z[3]*tmp1[7]+z[45]*tmp1[8]+z[87]*tmp1[9]+z[32]*tmp1[10]
            +z[74]*tmp1[11]+z[19]*tmp1[12]+z[61]*tmp1[13]+z[6]*tmp1[14]+z[48]*tmp1[15]+z[90]*tmp1[16]+z[35]*tmp1[17]+z[77]*tmp1[18]+z[22]*tmp1[19]+z[64]*tmp1[20]
            +z[9]*tmp1[21]+z[51]*tmp1[22]+z[93]*tmp1[23]+z[38]*tmp1[24]+z[80]*tmp1[25]+z[25]*tmp1[26]+z[67]*tmp1[27]+z[12]*tmp1[28]+z[54]*tmp1[29]+z[96]*tmp1[30]
            +z[41]*tmp1[31]+z[83]*tmp1[32]+z[28]*tmp1[33]+z[70]*tmp1[34]+z[15]*tmp1[35]+z[57]*tmp1[36]+z[2]*tmp1[37]+z[44]*tmp1[38]+z[86]*tmp1[39]+z[31]*tmp1[40]
            +z[73]*tmp1[41]+z[18]*tmp1[42]+z[60]*tmp1[43]+z[5]*tmp1[44]+z[47]*tmp1[45]+z[89]*tmp1[46]+z[34]*tmp1[47]+z[76]*tmp1[48]+z[21]*tmp1[49]+z[63]*tmp1[50]
            +z[8]*tmp1[51]+z[50]*tmp1[52]+z[92]*tmp1[53]+z[37]*tmp1[54]+z[79]*tmp1[55]+z[24]*tmp1[56]+z[66]*tmp1[57]+z[11]*tmp1[58]+z[53]*tmp1[59]+z[95]*tmp1[60]
            +z[40]*tmp1[61]+z[82]*tmp1[62]+z[27]*tmp1[63]+z[69]*tmp1[64]+z[14]*tmp1[65]+z[56]*tmp1[66]+z[1]*tmp1[67]+z[43]*tmp1[68]+z[85]*tmp1[69]+z[30]*tmp1[70]
            +z[72]*tmp1[71]+z[17]*tmp1[72]+z[59]*tmp1[73]+z[4]*tmp1[74]+z[46]*tmp1[75]+z[88]*tmp1[76]+z[33]*tmp1[77]+z[75]*tmp1[78]+z[20]*tmp1[79]+z[62]*tmp1[80]
            +z[7]*tmp1[81]+z[49]*tmp1[82]+z[91]*tmp1[83]+z[36]*tmp1[84]+z[78]*tmp1[85]+z[23]*tmp1[86]+z[65]*tmp1[87]+z[10]*tmp1[88]+z[52]*tmp1[89]+z[94]*tmp1[90]
            +z[39]*tmp1[91]+z[81]*tmp1[92]+z[26]*tmp1[93]+z[68]*tmp1[94]+z[13]*tmp1[95]+z[55]*tmp1[96];
            tab[i+43*stg_first]=tab[i+43*stg_first]+z[0]*tmp1[0]
            +z[43]*tmp1[1]+z[86]*tmp1[2]+z[32]*tmp1[3]+z[75]*tmp1[4]+z[21]*tmp1[5]+z[64]*tmp1[6]+z[10]*tmp1[7]+z[53]*tmp1[8]+z[96]*tmp1[9]+z[42]*tmp1[10]
            +z[85]*tmp1[11]+z[31]*tmp1[12]+z[74]*tmp1[13]+z[20]*tmp1[14]+z[63]*tmp1[15]+z[9]*tmp1[16]+z[52]*tmp1[17]+z[95]*tmp1[18]+z[41]*tmp1[19]+z[84]*tmp1[20]
            +z[30]*tmp1[21]+z[73]*tmp1[22]+z[19]*tmp1[23]+z[62]*tmp1[24]+z[8]*tmp1[25]+z[51]*tmp1[26]+z[94]*tmp1[27]+z[40]*tmp1[28]+z[83]*tmp1[29]+z[29]*tmp1[30]
            +z[72]*tmp1[31]+z[18]*tmp1[32]+z[61]*tmp1[33]+z[7]*tmp1[34]+z[50]*tmp1[35]+z[93]*tmp1[36]+z[39]*tmp1[37]+z[82]*tmp1[38]+z[28]*tmp1[39]+z[71]*tmp1[40]
            +z[17]*tmp1[41]+z[60]*tmp1[42]+z[6]*tmp1[43]+z[49]*tmp1[44]+z[92]*tmp1[45]+z[38]*tmp1[46]+z[81]*tmp1[47]+z[27]*tmp1[48]+z[70]*tmp1[49]+z[16]*tmp1[50]
            +z[59]*tmp1[51]+z[5]*tmp1[52]+z[48]*tmp1[53]+z[91]*tmp1[54]+z[37]*tmp1[55]+z[80]*tmp1[56]+z[26]*tmp1[57]+z[69]*tmp1[58]+z[15]*tmp1[59]+z[58]*tmp1[60]
            +z[4]*tmp1[61]+z[47]*tmp1[62]+z[90]*tmp1[63]+z[36]*tmp1[64]+z[79]*tmp1[65]+z[25]*tmp1[66]+z[68]*tmp1[67]+z[14]*tmp1[68]+z[57]*tmp1[69]+z[3]*tmp1[70]
            +z[46]*tmp1[71]+z[89]*tmp1[72]+z[35]*tmp1[73]+z[78]*tmp1[74]+z[24]*tmp1[75]+z[67]*tmp1[76]+z[13]*tmp1[77]+z[56]*tmp1[78]+z[2]*tmp1[79]+z[45]*tmp1[80]
            +z[88]*tmp1[81]+z[34]*tmp1[82]+z[77]*tmp1[83]+z[23]*tmp1[84]+z[66]*tmp1[85]+z[12]*tmp1[86]+z[55]*tmp1[87]+z[1]*tmp1[88]+z[44]*tmp1[89]+z[87]*tmp1[90]
            +z[33]*tmp1[91]+z[76]*tmp1[92]+z[22]*tmp1[93]+z[65]*tmp1[94]+z[11]*tmp1[95]+z[54]*tmp1[96];
            tab[i+44*stg_first]=tab[i+44*stg_first]+z[0]*tmp1[0]
            +z[44]*tmp1[1]+z[88]*tmp1[2]+z[35]*tmp1[3]+z[79]*tmp1[4]+z[26]*tmp1[5]+z[70]*tmp1[6]+z[17]*tmp1[7]+z[61]*tmp1[8]+z[8]*tmp1[9]+z[52]*tmp1[10]
            +z[96]*tmp1[11]+z[43]*tmp1[12]+z[87]*tmp1[13]+z[34]*tmp1[14]+z[78]*tmp1[15]+z[25]*tmp1[16]+z[69]*tmp1[17]+z[16]*tmp1[18]+z[60]*tmp1[19]+z[7]*tmp1[20]
            +z[51]*tmp1[21]+z[95]*tmp1[22]+z[42]*tmp1[23]+z[86]*tmp1[24]+z[33]*tmp1[25]+z[77]*tmp1[26]+z[24]*tmp1[27]+z[68]*tmp1[28]+z[15]*tmp1[29]+z[59]*tmp1[30]
            +z[6]*tmp1[31]+z[50]*tmp1[32]+z[94]*tmp1[33]+z[41]*tmp1[34]+z[85]*tmp1[35]+z[32]*tmp1[36]+z[76]*tmp1[37]+z[23]*tmp1[38]+z[67]*tmp1[39]+z[14]*tmp1[40]
            +z[58]*tmp1[41]+z[5]*tmp1[42]+z[49]*tmp1[43]+z[93]*tmp1[44]+z[40]*tmp1[45]+z[84]*tmp1[46]+z[31]*tmp1[47]+z[75]*tmp1[48]+z[22]*tmp1[49]+z[66]*tmp1[50]
            +z[13]*tmp1[51]+z[57]*tmp1[52]+z[4]*tmp1[53]+z[48]*tmp1[54]+z[92]*tmp1[55]+z[39]*tmp1[56]+z[83]*tmp1[57]+z[30]*tmp1[58]+z[74]*tmp1[59]+z[21]*tmp1[60]
            +z[65]*tmp1[61]+z[12]*tmp1[62]+z[56]*tmp1[63]+z[3]*tmp1[64]+z[47]*tmp1[65]+z[91]*tmp1[66]+z[38]*tmp1[67]+z[82]*tmp1[68]+z[29]*tmp1[69]+z[73]*tmp1[70]
            +z[20]*tmp1[71]+z[64]*tmp1[72]+z[11]*tmp1[73]+z[55]*tmp1[74]+z[2]*tmp1[75]+z[46]*tmp1[76]+z[90]*tmp1[77]+z[37]*tmp1[78]+z[81]*tmp1[79]+z[28]*tmp1[80]
            +z[72]*tmp1[81]+z[19]*tmp1[82]+z[63]*tmp1[83]+z[10]*tmp1[84]+z[54]*tmp1[85]+z[1]*tmp1[86]+z[45]*tmp1[87]+z[89]*tmp1[88]+z[36]*tmp1[89]+z[80]*tmp1[90]
            +z[27]*tmp1[91]+z[71]*tmp1[92]+z[18]*tmp1[93]+z[62]*tmp1[94]+z[9]*tmp1[95]+z[53]*tmp1[96];
            tab[i+45*stg_first]=tab[i+45*stg_first]+z[0]*tmp1[0]
            +z[45]*tmp1[1]+z[90]*tmp1[2]+z[38]*tmp1[3]+z[83]*tmp1[4]+z[31]*tmp1[5]+z[76]*tmp1[6]+z[24]*tmp1[7]+z[69]*tmp1[8]+z[17]*tmp1[9]+z[62]*tmp1[10]
            +z[10]*tmp1[11]+z[55]*tmp1[12]+z[3]*tmp1[13]+z[48]*tmp1[14]+z[93]*tmp1[15]+z[41]*tmp1[16]+z[86]*tmp1[17]+z[34]*tmp1[18]+z[79]*tmp1[19]+z[27]*tmp1[20]
            +z[72]*tmp1[21]+z[20]*tmp1[22]+z[65]*tmp1[23]+z[13]*tmp1[24]+z[58]*tmp1[25]+z[6]*tmp1[26]+z[51]*tmp1[27]+z[96]*tmp1[28]+z[44]*tmp1[29]+z[89]*tmp1[30]
            +z[37]*tmp1[31]+z[82]*tmp1[32]+z[30]*tmp1[33]+z[75]*tmp1[34]+z[23]*tmp1[35]+z[68]*tmp1[36]+z[16]*tmp1[37]+z[61]*tmp1[38]+z[9]*tmp1[39]+z[54]*tmp1[40]
            +z[2]*tmp1[41]+z[47]*tmp1[42]+z[92]*tmp1[43]+z[40]*tmp1[44]+z[85]*tmp1[45]+z[33]*tmp1[46]+z[78]*tmp1[47]+z[26]*tmp1[48]+z[71]*tmp1[49]+z[19]*tmp1[50]
            +z[64]*tmp1[51]+z[12]*tmp1[52]+z[57]*tmp1[53]+z[5]*tmp1[54]+z[50]*tmp1[55]+z[95]*tmp1[56]+z[43]*tmp1[57]+z[88]*tmp1[58]+z[36]*tmp1[59]+z[81]*tmp1[60]
            +z[29]*tmp1[61]+z[74]*tmp1[62]+z[22]*tmp1[63]+z[67]*tmp1[64]+z[15]*tmp1[65]+z[60]*tmp1[66]+z[8]*tmp1[67]+z[53]*tmp1[68]+z[1]*tmp1[69]+z[46]*tmp1[70]
            +z[91]*tmp1[71]+z[39]*tmp1[72]+z[84]*tmp1[73]+z[32]*tmp1[74]+z[77]*tmp1[75]+z[25]*tmp1[76]+z[70]*tmp1[77]+z[18]*tmp1[78]+z[63]*tmp1[79]+z[11]*tmp1[80]
            +z[56]*tmp1[81]+z[4]*tmp1[82]+z[49]*tmp1[83]+z[94]*tmp1[84]+z[42]*tmp1[85]+z[87]*tmp1[86]+z[35]*tmp1[87]+z[80]*tmp1[88]+z[28]*tmp1[89]+z[73]*tmp1[90]
            +z[21]*tmp1[91]+z[66]*tmp1[92]+z[14]*tmp1[93]+z[59]*tmp1[94]+z[7]*tmp1[95]+z[52]*tmp1[96];
            tab[i+46*stg_first]=tab[i+46*stg_first]+z[0]*tmp1[0]
            +z[46]*tmp1[1]+z[92]*tmp1[2]+z[41]*tmp1[3]+z[87]*tmp1[4]+z[36]*tmp1[5]+z[82]*tmp1[6]+z[31]*tmp1[7]+z[77]*tmp1[8]+z[26]*tmp1[9]+z[72]*tmp1[10]
            +z[21]*tmp1[11]+z[67]*tmp1[12]+z[16]*tmp1[13]+z[62]*tmp1[14]+z[11]*tmp1[15]+z[57]*tmp1[16]+z[6]*tmp1[17]+z[52]*tmp1[18]+z[1]*tmp1[19]+z[47]*tmp1[20]
            +z[93]*tmp1[21]+z[42]*tmp1[22]+z[88]*tmp1[23]+z[37]*tmp1[24]+z[83]*tmp1[25]+z[32]*tmp1[26]+z[78]*tmp1[27]+z[27]*tmp1[28]+z[73]*tmp1[29]+z[22]*tmp1[30]
            +z[68]*tmp1[31]+z[17]*tmp1[32]+z[63]*tmp1[33]+z[12]*tmp1[34]+z[58]*tmp1[35]+z[7]*tmp1[36]+z[53]*tmp1[37]+z[2]*tmp1[38]+z[48]*tmp1[39]+z[94]*tmp1[40]
            +z[43]*tmp1[41]+z[89]*tmp1[42]+z[38]*tmp1[43]+z[84]*tmp1[44]+z[33]*tmp1[45]+z[79]*tmp1[46]+z[28]*tmp1[47]+z[74]*tmp1[48]+z[23]*tmp1[49]+z[69]*tmp1[50]
            +z[18]*tmp1[51]+z[64]*tmp1[52]+z[13]*tmp1[53]+z[59]*tmp1[54]+z[8]*tmp1[55]+z[54]*tmp1[56]+z[3]*tmp1[57]+z[49]*tmp1[58]+z[95]*tmp1[59]+z[44]*tmp1[60]
            +z[90]*tmp1[61]+z[39]*tmp1[62]+z[85]*tmp1[63]+z[34]*tmp1[64]+z[80]*tmp1[65]+z[29]*tmp1[66]+z[75]*tmp1[67]+z[24]*tmp1[68]+z[70]*tmp1[69]+z[19]*tmp1[70]
            +z[65]*tmp1[71]+z[14]*tmp1[72]+z[60]*tmp1[73]+z[9]*tmp1[74]+z[55]*tmp1[75]+z[4]*tmp1[76]+z[50]*tmp1[77]+z[96]*tmp1[78]+z[45]*tmp1[79]+z[91]*tmp1[80]
            +z[40]*tmp1[81]+z[86]*tmp1[82]+z[35]*tmp1[83]+z[81]*tmp1[84]+z[30]*tmp1[85]+z[76]*tmp1[86]+z[25]*tmp1[87]+z[71]*tmp1[88]+z[20]*tmp1[89]+z[66]*tmp1[90]
            +z[15]*tmp1[91]+z[61]*tmp1[92]+z[10]*tmp1[93]+z[56]*tmp1[94]+z[5]*tmp1[95]+z[51]*tmp1[96];
            tab[i+47*stg_first]=tab[i+47*stg_first]+z[0]*tmp1[0]
            +z[47]*tmp1[1]+z[94]*tmp1[2]+z[44]*tmp1[3]+z[91]*tmp1[4]+z[41]*tmp1[5]+z[88]*tmp1[6]+z[38]*tmp1[7]+z[85]*tmp1[8]+z[35]*tmp1[9]+z[82]*tmp1[10]
            +z[32]*tmp1[11]+z[79]*tmp1[12]+z[29]*tmp1[13]+z[76]*tmp1[14]+z[26]*tmp1[15]+z[73]*tmp1[16]+z[23]*tmp1[17]+z[70]*tmp1[18]+z[20]*tmp1[19]+z[67]*tmp1[20]
            +z[17]*tmp1[21]+z[64]*tmp1[22]+z[14]*tmp1[23]+z[61]*tmp1[24]+z[11]*tmp1[25]+z[58]*tmp1[26]+z[8]*tmp1[27]+z[55]*tmp1[28]+z[5]*tmp1[29]+z[52]*tmp1[30]
            +z[2]*tmp1[31]+z[49]*tmp1[32]+z[96]*tmp1[33]+z[46]*tmp1[34]+z[93]*tmp1[35]+z[43]*tmp1[36]+z[90]*tmp1[37]+z[40]*tmp1[38]+z[87]*tmp1[39]+z[37]*tmp1[40]
            +z[84]*tmp1[41]+z[34]*tmp1[42]+z[81]*tmp1[43]+z[31]*tmp1[44]+z[78]*tmp1[45]+z[28]*tmp1[46]+z[75]*tmp1[47]+z[25]*tmp1[48]+z[72]*tmp1[49]+z[22]*tmp1[50]
            +z[69]*tmp1[51]+z[19]*tmp1[52]+z[66]*tmp1[53]+z[16]*tmp1[54]+z[63]*tmp1[55]+z[13]*tmp1[56]+z[60]*tmp1[57]+z[10]*tmp1[58]+z[57]*tmp1[59]+z[7]*tmp1[60]
            +z[54]*tmp1[61]+z[4]*tmp1[62]+z[51]*tmp1[63]+z[1]*tmp1[64]+z[48]*tmp1[65]+z[95]*tmp1[66]+z[45]*tmp1[67]+z[92]*tmp1[68]+z[42]*tmp1[69]+z[89]*tmp1[70]
            +z[39]*tmp1[71]+z[86]*tmp1[72]+z[36]*tmp1[73]+z[83]*tmp1[74]+z[33]*tmp1[75]+z[80]*tmp1[76]+z[30]*tmp1[77]+z[77]*tmp1[78]+z[27]*tmp1[79]+z[74]*tmp1[80]
            +z[24]*tmp1[81]+z[71]*tmp1[82]+z[21]*tmp1[83]+z[68]*tmp1[84]+z[18]*tmp1[85]+z[65]*tmp1[86]+z[15]*tmp1[87]+z[62]*tmp1[88]+z[12]*tmp1[89]+z[59]*tmp1[90]
            +z[9]*tmp1[91]+z[56]*tmp1[92]+z[6]*tmp1[93]+z[53]*tmp1[94]+z[3]*tmp1[95]+z[50]*tmp1[96];
            tab[i+48*stg_first]=tab[i+48*stg_first]+z[0]*tmp1[0]
            +z[48]*tmp1[1]+z[96]*tmp1[2]+z[47]*tmp1[3]+z[95]*tmp1[4]+z[46]*tmp1[5]+z[94]*tmp1[6]+z[45]*tmp1[7]+z[93]*tmp1[8]+z[44]*tmp1[9]+z[92]*tmp1[10]
            +z[43]*tmp1[11]+z[91]*tmp1[12]+z[42]*tmp1[13]+z[90]*tmp1[14]+z[41]*tmp1[15]+z[89]*tmp1[16]+z[40]*tmp1[17]+z[88]*tmp1[18]+z[39]*tmp1[19]+z[87]*tmp1[20]
            +z[38]*tmp1[21]+z[86]*tmp1[22]+z[37]*tmp1[23]+z[85]*tmp1[24]+z[36]*tmp1[25]+z[84]*tmp1[26]+z[35]*tmp1[27]+z[83]*tmp1[28]+z[34]*tmp1[29]+z[82]*tmp1[30]
            +z[33]*tmp1[31]+z[81]*tmp1[32]+z[32]*tmp1[33]+z[80]*tmp1[34]+z[31]*tmp1[35]+z[79]*tmp1[36]+z[30]*tmp1[37]+z[78]*tmp1[38]+z[29]*tmp1[39]+z[77]*tmp1[40]
            +z[28]*tmp1[41]+z[76]*tmp1[42]+z[27]*tmp1[43]+z[75]*tmp1[44]+z[26]*tmp1[45]+z[74]*tmp1[46]+z[25]*tmp1[47]+z[73]*tmp1[48]+z[24]*tmp1[49]+z[72]*tmp1[50]
            +z[23]*tmp1[51]+z[71]*tmp1[52]+z[22]*tmp1[53]+z[70]*tmp1[54]+z[21]*tmp1[55]+z[69]*tmp1[56]+z[20]*tmp1[57]+z[68]*tmp1[58]+z[19]*tmp1[59]+z[67]*tmp1[60]
            +z[18]*tmp1[61]+z[66]*tmp1[62]+z[17]*tmp1[63]+z[65]*tmp1[64]+z[16]*tmp1[65]+z[64]*tmp1[66]+z[15]*tmp1[67]+z[63]*tmp1[68]+z[14]*tmp1[69]+z[62]*tmp1[70]
            +z[13]*tmp1[71]+z[61]*tmp1[72]+z[12]*tmp1[73]+z[60]*tmp1[74]+z[11]*tmp1[75]+z[59]*tmp1[76]+z[10]*tmp1[77]+z[58]*tmp1[78]+z[9]*tmp1[79]+z[57]*tmp1[80]
            +z[8]*tmp1[81]+z[56]*tmp1[82]+z[7]*tmp1[83]+z[55]*tmp1[84]+z[6]*tmp1[85]+z[54]*tmp1[86]+z[5]*tmp1[87]+z[53]*tmp1[88]+z[4]*tmp1[89]+z[52]*tmp1[90]
            +z[3]*tmp1[91]+z[51]*tmp1[92]+z[2]*tmp1[93]+z[50]*tmp1[94]+z[1]*tmp1[95]+z[49]*tmp1[96];
            tab[i+49*stg_first]=tab[i+49*stg_first]+z[0]*tmp1[0]
            +z[49]*tmp1[1]+z[1]*tmp1[2]+z[50]*tmp1[3]+z[2]*tmp1[4]+z[51]*tmp1[5]+z[3]*tmp1[6]+z[52]*tmp1[7]+z[4]*tmp1[8]+z[53]*tmp1[9]+z[5]*tmp1[10]
            +z[54]*tmp1[11]+z[6]*tmp1[12]+z[55]*tmp1[13]+z[7]*tmp1[14]+z[56]*tmp1[15]+z[8]*tmp1[16]+z[57]*tmp1[17]+z[9]*tmp1[18]+z[58]*tmp1[19]+z[10]*tmp1[20]
            +z[59]*tmp1[21]+z[11]*tmp1[22]+z[60]*tmp1[23]+z[12]*tmp1[24]+z[61]*tmp1[25]+z[13]*tmp1[26]+z[62]*tmp1[27]+z[14]*tmp1[28]+z[63]*tmp1[29]+z[15]*tmp1[30]
            +z[64]*tmp1[31]+z[16]*tmp1[32]+z[65]*tmp1[33]+z[17]*tmp1[34]+z[66]*tmp1[35]+z[18]*tmp1[36]+z[67]*tmp1[37]+z[19]*tmp1[38]+z[68]*tmp1[39]+z[20]*tmp1[40]
            +z[69]*tmp1[41]+z[21]*tmp1[42]+z[70]*tmp1[43]+z[22]*tmp1[44]+z[71]*tmp1[45]+z[23]*tmp1[46]+z[72]*tmp1[47]+z[24]*tmp1[48]+z[73]*tmp1[49]+z[25]*tmp1[50]
            +z[74]*tmp1[51]+z[26]*tmp1[52]+z[75]*tmp1[53]+z[27]*tmp1[54]+z[76]*tmp1[55]+z[28]*tmp1[56]+z[77]*tmp1[57]+z[29]*tmp1[58]+z[78]*tmp1[59]+z[30]*tmp1[60]
            +z[79]*tmp1[61]+z[31]*tmp1[62]+z[80]*tmp1[63]+z[32]*tmp1[64]+z[81]*tmp1[65]+z[33]*tmp1[66]+z[82]*tmp1[67]+z[34]*tmp1[68]+z[83]*tmp1[69]+z[35]*tmp1[70]
            +z[84]*tmp1[71]+z[36]*tmp1[72]+z[85]*tmp1[73]+z[37]*tmp1[74]+z[86]*tmp1[75]+z[38]*tmp1[76]+z[87]*tmp1[77]+z[39]*tmp1[78]+z[88]*tmp1[79]+z[40]*tmp1[80]
            +z[89]*tmp1[81]+z[41]*tmp1[82]+z[90]*tmp1[83]+z[42]*tmp1[84]+z[91]*tmp1[85]+z[43]*tmp1[86]+z[92]*tmp1[87]+z[44]*tmp1[88]+z[93]*tmp1[89]+z[45]*tmp1[90]
            +z[94]*tmp1[91]+z[46]*tmp1[92]+z[95]*tmp1[93]+z[47]*tmp1[94]+z[96]*tmp1[95]+z[48]*tmp1[96];
            tab[i+50*stg_first]=tab[i+50*stg_first]+z[0]*tmp1[0]
            +z[50]*tmp1[1]+z[3]*tmp1[2]+z[53]*tmp1[3]+z[6]*tmp1[4]+z[56]*tmp1[5]+z[9]*tmp1[6]+z[59]*tmp1[7]+z[12]*tmp1[8]+z[62]*tmp1[9]+z[15]*tmp1[10]
            +z[65]*tmp1[11]+z[18]*tmp1[12]+z[68]*tmp1[13]+z[21]*tmp1[14]+z[71]*tmp1[15]+z[24]*tmp1[16]+z[74]*tmp1[17]+z[27]*tmp1[18]+z[77]*tmp1[19]+z[30]*tmp1[20]
            +z[80]*tmp1[21]+z[33]*tmp1[22]+z[83]*tmp1[23]+z[36]*tmp1[24]+z[86]*tmp1[25]+z[39]*tmp1[26]+z[89]*tmp1[27]+z[42]*tmp1[28]+z[92]*tmp1[29]+z[45]*tmp1[30]
            +z[95]*tmp1[31]+z[48]*tmp1[32]+z[1]*tmp1[33]+z[51]*tmp1[34]+z[4]*tmp1[35]+z[54]*tmp1[36]+z[7]*tmp1[37]+z[57]*tmp1[38]+z[10]*tmp1[39]+z[60]*tmp1[40]
            +z[13]*tmp1[41]+z[63]*tmp1[42]+z[16]*tmp1[43]+z[66]*tmp1[44]+z[19]*tmp1[45]+z[69]*tmp1[46]+z[22]*tmp1[47]+z[72]*tmp1[48]+z[25]*tmp1[49]+z[75]*tmp1[50]
            +z[28]*tmp1[51]+z[78]*tmp1[52]+z[31]*tmp1[53]+z[81]*tmp1[54]+z[34]*tmp1[55]+z[84]*tmp1[56]+z[37]*tmp1[57]+z[87]*tmp1[58]+z[40]*tmp1[59]+z[90]*tmp1[60]
            +z[43]*tmp1[61]+z[93]*tmp1[62]+z[46]*tmp1[63]+z[96]*tmp1[64]+z[49]*tmp1[65]+z[2]*tmp1[66]+z[52]*tmp1[67]+z[5]*tmp1[68]+z[55]*tmp1[69]+z[8]*tmp1[70]
            +z[58]*tmp1[71]+z[11]*tmp1[72]+z[61]*tmp1[73]+z[14]*tmp1[74]+z[64]*tmp1[75]+z[17]*tmp1[76]+z[67]*tmp1[77]+z[20]*tmp1[78]+z[70]*tmp1[79]+z[23]*tmp1[80]
            +z[73]*tmp1[81]+z[26]*tmp1[82]+z[76]*tmp1[83]+z[29]*tmp1[84]+z[79]*tmp1[85]+z[32]*tmp1[86]+z[82]*tmp1[87]+z[35]*tmp1[88]+z[85]*tmp1[89]+z[38]*tmp1[90]
            +z[88]*tmp1[91]+z[41]*tmp1[92]+z[91]*tmp1[93]+z[44]*tmp1[94]+z[94]*tmp1[95]+z[47]*tmp1[96];
            tab[i+51*stg_first]=tab[i+51*stg_first]+z[0]*tmp1[0]
            +z[51]*tmp1[1]+z[5]*tmp1[2]+z[56]*tmp1[3]+z[10]*tmp1[4]+z[61]*tmp1[5]+z[15]*tmp1[6]+z[66]*tmp1[7]+z[20]*tmp1[8]+z[71]*tmp1[9]+z[25]*tmp1[10]
            +z[76]*tmp1[11]+z[30]*tmp1[12]+z[81]*tmp1[13]+z[35]*tmp1[14]+z[86]*tmp1[15]+z[40]*tmp1[16]+z[91]*tmp1[17]+z[45]*tmp1[18]+z[96]*tmp1[19]+z[50]*tmp1[20]
            +z[4]*tmp1[21]+z[55]*tmp1[22]+z[9]*tmp1[23]+z[60]*tmp1[24]+z[14]*tmp1[25]+z[65]*tmp1[26]+z[19]*tmp1[27]+z[70]*tmp1[28]+z[24]*tmp1[29]+z[75]*tmp1[30]
            +z[29]*tmp1[31]+z[80]*tmp1[32]+z[34]*tmp1[33]+z[85]*tmp1[34]+z[39]*tmp1[35]+z[90]*tmp1[36]+z[44]*tmp1[37]+z[95]*tmp1[38]+z[49]*tmp1[39]+z[3]*tmp1[40]
            +z[54]*tmp1[41]+z[8]*tmp1[42]+z[59]*tmp1[43]+z[13]*tmp1[44]+z[64]*tmp1[45]+z[18]*tmp1[46]+z[69]*tmp1[47]+z[23]*tmp1[48]+z[74]*tmp1[49]+z[28]*tmp1[50]
            +z[79]*tmp1[51]+z[33]*tmp1[52]+z[84]*tmp1[53]+z[38]*tmp1[54]+z[89]*tmp1[55]+z[43]*tmp1[56]+z[94]*tmp1[57]+z[48]*tmp1[58]+z[2]*tmp1[59]+z[53]*tmp1[60]
            +z[7]*tmp1[61]+z[58]*tmp1[62]+z[12]*tmp1[63]+z[63]*tmp1[64]+z[17]*tmp1[65]+z[68]*tmp1[66]+z[22]*tmp1[67]+z[73]*tmp1[68]+z[27]*tmp1[69]+z[78]*tmp1[70]
            +z[32]*tmp1[71]+z[83]*tmp1[72]+z[37]*tmp1[73]+z[88]*tmp1[74]+z[42]*tmp1[75]+z[93]*tmp1[76]+z[47]*tmp1[77]+z[1]*tmp1[78]+z[52]*tmp1[79]+z[6]*tmp1[80]
            +z[57]*tmp1[81]+z[11]*tmp1[82]+z[62]*tmp1[83]+z[16]*tmp1[84]+z[67]*tmp1[85]+z[21]*tmp1[86]+z[72]*tmp1[87]+z[26]*tmp1[88]+z[77]*tmp1[89]+z[31]*tmp1[90]
            +z[82]*tmp1[91]+z[36]*tmp1[92]+z[87]*tmp1[93]+z[41]*tmp1[94]+z[92]*tmp1[95]+z[46]*tmp1[96];
            tab[i+52*stg_first]=tab[i+52*stg_first]+z[0]*tmp1[0]
            +z[52]*tmp1[1]+z[7]*tmp1[2]+z[59]*tmp1[3]+z[14]*tmp1[4]+z[66]*tmp1[5]+z[21]*tmp1[6]+z[73]*tmp1[7]+z[28]*tmp1[8]+z[80]*tmp1[9]+z[35]*tmp1[10]
            +z[87]*tmp1[11]+z[42]*tmp1[12]+z[94]*tmp1[13]+z[49]*tmp1[14]+z[4]*tmp1[15]+z[56]*tmp1[16]+z[11]*tmp1[17]+z[63]*tmp1[18]+z[18]*tmp1[19]+z[70]*tmp1[20]
            +z[25]*tmp1[21]+z[77]*tmp1[22]+z[32]*tmp1[23]+z[84]*tmp1[24]+z[39]*tmp1[25]+z[91]*tmp1[26]+z[46]*tmp1[27]+z[1]*tmp1[28]+z[53]*tmp1[29]+z[8]*tmp1[30]
            +z[60]*tmp1[31]+z[15]*tmp1[32]+z[67]*tmp1[33]+z[22]*tmp1[34]+z[74]*tmp1[35]+z[29]*tmp1[36]+z[81]*tmp1[37]+z[36]*tmp1[38]+z[88]*tmp1[39]+z[43]*tmp1[40]
            +z[95]*tmp1[41]+z[50]*tmp1[42]+z[5]*tmp1[43]+z[57]*tmp1[44]+z[12]*tmp1[45]+z[64]*tmp1[46]+z[19]*tmp1[47]+z[71]*tmp1[48]+z[26]*tmp1[49]+z[78]*tmp1[50]
            +z[33]*tmp1[51]+z[85]*tmp1[52]+z[40]*tmp1[53]+z[92]*tmp1[54]+z[47]*tmp1[55]+z[2]*tmp1[56]+z[54]*tmp1[57]+z[9]*tmp1[58]+z[61]*tmp1[59]+z[16]*tmp1[60]
            +z[68]*tmp1[61]+z[23]*tmp1[62]+z[75]*tmp1[63]+z[30]*tmp1[64]+z[82]*tmp1[65]+z[37]*tmp1[66]+z[89]*tmp1[67]+z[44]*tmp1[68]+z[96]*tmp1[69]+z[51]*tmp1[70]
            +z[6]*tmp1[71]+z[58]*tmp1[72]+z[13]*tmp1[73]+z[65]*tmp1[74]+z[20]*tmp1[75]+z[72]*tmp1[76]+z[27]*tmp1[77]+z[79]*tmp1[78]+z[34]*tmp1[79]+z[86]*tmp1[80]
            +z[41]*tmp1[81]+z[93]*tmp1[82]+z[48]*tmp1[83]+z[3]*tmp1[84]+z[55]*tmp1[85]+z[10]*tmp1[86]+z[62]*tmp1[87]+z[17]*tmp1[88]+z[69]*tmp1[89]+z[24]*tmp1[90]
            +z[76]*tmp1[91]+z[31]*tmp1[92]+z[83]*tmp1[93]+z[38]*tmp1[94]+z[90]*tmp1[95]+z[45]*tmp1[96];
            tab[i+53*stg_first]=tab[i+53*stg_first]+z[0]*tmp1[0]
            +z[53]*tmp1[1]+z[9]*tmp1[2]+z[62]*tmp1[3]+z[18]*tmp1[4]+z[71]*tmp1[5]+z[27]*tmp1[6]+z[80]*tmp1[7]+z[36]*tmp1[8]+z[89]*tmp1[9]+z[45]*tmp1[10]
            +z[1]*tmp1[11]+z[54]*tmp1[12]+z[10]*tmp1[13]+z[63]*tmp1[14]+z[19]*tmp1[15]+z[72]*tmp1[16]+z[28]*tmp1[17]+z[81]*tmp1[18]+z[37]*tmp1[19]+z[90]*tmp1[20]
            +z[46]*tmp1[21]+z[2]*tmp1[22]+z[55]*tmp1[23]+z[11]*tmp1[24]+z[64]*tmp1[25]+z[20]*tmp1[26]+z[73]*tmp1[27]+z[29]*tmp1[28]+z[82]*tmp1[29]+z[38]*tmp1[30]
            +z[91]*tmp1[31]+z[47]*tmp1[32]+z[3]*tmp1[33]+z[56]*tmp1[34]+z[12]*tmp1[35]+z[65]*tmp1[36]+z[21]*tmp1[37]+z[74]*tmp1[38]+z[30]*tmp1[39]+z[83]*tmp1[40]
            +z[39]*tmp1[41]+z[92]*tmp1[42]+z[48]*tmp1[43]+z[4]*tmp1[44]+z[57]*tmp1[45]+z[13]*tmp1[46]+z[66]*tmp1[47]+z[22]*tmp1[48]+z[75]*tmp1[49]+z[31]*tmp1[50]
            +z[84]*tmp1[51]+z[40]*tmp1[52]+z[93]*tmp1[53]+z[49]*tmp1[54]+z[5]*tmp1[55]+z[58]*tmp1[56]+z[14]*tmp1[57]+z[67]*tmp1[58]+z[23]*tmp1[59]+z[76]*tmp1[60]
            +z[32]*tmp1[61]+z[85]*tmp1[62]+z[41]*tmp1[63]+z[94]*tmp1[64]+z[50]*tmp1[65]+z[6]*tmp1[66]+z[59]*tmp1[67]+z[15]*tmp1[68]+z[68]*tmp1[69]+z[24]*tmp1[70]
            +z[77]*tmp1[71]+z[33]*tmp1[72]+z[86]*tmp1[73]+z[42]*tmp1[74]+z[95]*tmp1[75]+z[51]*tmp1[76]+z[7]*tmp1[77]+z[60]*tmp1[78]+z[16]*tmp1[79]+z[69]*tmp1[80]
            +z[25]*tmp1[81]+z[78]*tmp1[82]+z[34]*tmp1[83]+z[87]*tmp1[84]+z[43]*tmp1[85]+z[96]*tmp1[86]+z[52]*tmp1[87]+z[8]*tmp1[88]+z[61]*tmp1[89]+z[17]*tmp1[90]
            +z[70]*tmp1[91]+z[26]*tmp1[92]+z[79]*tmp1[93]+z[35]*tmp1[94]+z[88]*tmp1[95]+z[44]*tmp1[96];
            tab[i+54*stg_first]=tab[i+54*stg_first]+z[0]*tmp1[0]
            +z[54]*tmp1[1]+z[11]*tmp1[2]+z[65]*tmp1[3]+z[22]*tmp1[4]+z[76]*tmp1[5]+z[33]*tmp1[6]+z[87]*tmp1[7]+z[44]*tmp1[8]+z[1]*tmp1[9]+z[55]*tmp1[10]
            +z[12]*tmp1[11]+z[66]*tmp1[12]+z[23]*tmp1[13]+z[77]*tmp1[14]+z[34]*tmp1[15]+z[88]*tmp1[16]+z[45]*tmp1[17]+z[2]*tmp1[18]+z[56]*tmp1[19]+z[13]*tmp1[20]
            +z[67]*tmp1[21]+z[24]*tmp1[22]+z[78]*tmp1[23]+z[35]*tmp1[24]+z[89]*tmp1[25]+z[46]*tmp1[26]+z[3]*tmp1[27]+z[57]*tmp1[28]+z[14]*tmp1[29]+z[68]*tmp1[30]
            +z[25]*tmp1[31]+z[79]*tmp1[32]+z[36]*tmp1[33]+z[90]*tmp1[34]+z[47]*tmp1[35]+z[4]*tmp1[36]+z[58]*tmp1[37]+z[15]*tmp1[38]+z[69]*tmp1[39]+z[26]*tmp1[40]
            +z[80]*tmp1[41]+z[37]*tmp1[42]+z[91]*tmp1[43]+z[48]*tmp1[44]+z[5]*tmp1[45]+z[59]*tmp1[46]+z[16]*tmp1[47]+z[70]*tmp1[48]+z[27]*tmp1[49]+z[81]*tmp1[50]
            +z[38]*tmp1[51]+z[92]*tmp1[52]+z[49]*tmp1[53]+z[6]*tmp1[54]+z[60]*tmp1[55]+z[17]*tmp1[56]+z[71]*tmp1[57]+z[28]*tmp1[58]+z[82]*tmp1[59]+z[39]*tmp1[60]
            +z[93]*tmp1[61]+z[50]*tmp1[62]+z[7]*tmp1[63]+z[61]*tmp1[64]+z[18]*tmp1[65]+z[72]*tmp1[66]+z[29]*tmp1[67]+z[83]*tmp1[68]+z[40]*tmp1[69]+z[94]*tmp1[70]
            +z[51]*tmp1[71]+z[8]*tmp1[72]+z[62]*tmp1[73]+z[19]*tmp1[74]+z[73]*tmp1[75]+z[30]*tmp1[76]+z[84]*tmp1[77]+z[41]*tmp1[78]+z[95]*tmp1[79]+z[52]*tmp1[80]
            +z[9]*tmp1[81]+z[63]*tmp1[82]+z[20]*tmp1[83]+z[74]*tmp1[84]+z[31]*tmp1[85]+z[85]*tmp1[86]+z[42]*tmp1[87]+z[96]*tmp1[88]+z[53]*tmp1[89]+z[10]*tmp1[90]
            +z[64]*tmp1[91]+z[21]*tmp1[92]+z[75]*tmp1[93]+z[32]*tmp1[94]+z[86]*tmp1[95]+z[43]*tmp1[96];
            tab[i+55*stg_first]=tab[i+55*stg_first]+z[0]*tmp1[0]
            +z[55]*tmp1[1]+z[13]*tmp1[2]+z[68]*tmp1[3]+z[26]*tmp1[4]+z[81]*tmp1[5]+z[39]*tmp1[6]+z[94]*tmp1[7]+z[52]*tmp1[8]+z[10]*tmp1[9]+z[65]*tmp1[10]
            +z[23]*tmp1[11]+z[78]*tmp1[12]+z[36]*tmp1[13]+z[91]*tmp1[14]+z[49]*tmp1[15]+z[7]*tmp1[16]+z[62]*tmp1[17]+z[20]*tmp1[18]+z[75]*tmp1[19]+z[33]*tmp1[20]
            +z[88]*tmp1[21]+z[46]*tmp1[22]+z[4]*tmp1[23]+z[59]*tmp1[24]+z[17]*tmp1[25]+z[72]*tmp1[26]+z[30]*tmp1[27]+z[85]*tmp1[28]+z[43]*tmp1[29]+z[1]*tmp1[30]
            +z[56]*tmp1[31]+z[14]*tmp1[32]+z[69]*tmp1[33]+z[27]*tmp1[34]+z[82]*tmp1[35]+z[40]*tmp1[36]+z[95]*tmp1[37]+z[53]*tmp1[38]+z[11]*tmp1[39]+z[66]*tmp1[40]
            +z[24]*tmp1[41]+z[79]*tmp1[42]+z[37]*tmp1[43]+z[92]*tmp1[44]+z[50]*tmp1[45]+z[8]*tmp1[46]+z[63]*tmp1[47]+z[21]*tmp1[48]+z[76]*tmp1[49]+z[34]*tmp1[50]
            +z[89]*tmp1[51]+z[47]*tmp1[52]+z[5]*tmp1[53]+z[60]*tmp1[54]+z[18]*tmp1[55]+z[73]*tmp1[56]+z[31]*tmp1[57]+z[86]*tmp1[58]+z[44]*tmp1[59]+z[2]*tmp1[60]
            +z[57]*tmp1[61]+z[15]*tmp1[62]+z[70]*tmp1[63]+z[28]*tmp1[64]+z[83]*tmp1[65]+z[41]*tmp1[66]+z[96]*tmp1[67]+z[54]*tmp1[68]+z[12]*tmp1[69]+z[67]*tmp1[70]
            +z[25]*tmp1[71]+z[80]*tmp1[72]+z[38]*tmp1[73]+z[93]*tmp1[74]+z[51]*tmp1[75]+z[9]*tmp1[76]+z[64]*tmp1[77]+z[22]*tmp1[78]+z[77]*tmp1[79]+z[35]*tmp1[80]
            +z[90]*tmp1[81]+z[48]*tmp1[82]+z[6]*tmp1[83]+z[61]*tmp1[84]+z[19]*tmp1[85]+z[74]*tmp1[86]+z[32]*tmp1[87]+z[87]*tmp1[88]+z[45]*tmp1[89]+z[3]*tmp1[90]
            +z[58]*tmp1[91]+z[16]*tmp1[92]+z[71]*tmp1[93]+z[29]*tmp1[94]+z[84]*tmp1[95]+z[42]*tmp1[96];
            tab[i+56*stg_first]=tab[i+56*stg_first]+z[0]*tmp1[0]
            +z[56]*tmp1[1]+z[15]*tmp1[2]+z[71]*tmp1[3]+z[30]*tmp1[4]+z[86]*tmp1[5]+z[45]*tmp1[6]+z[4]*tmp1[7]+z[60]*tmp1[8]+z[19]*tmp1[9]+z[75]*tmp1[10]
            +z[34]*tmp1[11]+z[90]*tmp1[12]+z[49]*tmp1[13]+z[8]*tmp1[14]+z[64]*tmp1[15]+z[23]*tmp1[16]+z[79]*tmp1[17]+z[38]*tmp1[18]+z[94]*tmp1[19]+z[53]*tmp1[20]
            +z[12]*tmp1[21]+z[68]*tmp1[22]+z[27]*tmp1[23]+z[83]*tmp1[24]+z[42]*tmp1[25]+z[1]*tmp1[26]+z[57]*tmp1[27]+z[16]*tmp1[28]+z[72]*tmp1[29]+z[31]*tmp1[30]
            +z[87]*tmp1[31]+z[46]*tmp1[32]+z[5]*tmp1[33]+z[61]*tmp1[34]+z[20]*tmp1[35]+z[76]*tmp1[36]+z[35]*tmp1[37]+z[91]*tmp1[38]+z[50]*tmp1[39]+z[9]*tmp1[40]
            +z[65]*tmp1[41]+z[24]*tmp1[42]+z[80]*tmp1[43]+z[39]*tmp1[44]+z[95]*tmp1[45]+z[54]*tmp1[46]+z[13]*tmp1[47]+z[69]*tmp1[48]+z[28]*tmp1[49]+z[84]*tmp1[50]
            +z[43]*tmp1[51]+z[2]*tmp1[52]+z[58]*tmp1[53]+z[17]*tmp1[54]+z[73]*tmp1[55]+z[32]*tmp1[56]+z[88]*tmp1[57]+z[47]*tmp1[58]+z[6]*tmp1[59]+z[62]*tmp1[60]
            +z[21]*tmp1[61]+z[77]*tmp1[62]+z[36]*tmp1[63]+z[92]*tmp1[64]+z[51]*tmp1[65]+z[10]*tmp1[66]+z[66]*tmp1[67]+z[25]*tmp1[68]+z[81]*tmp1[69]+z[40]*tmp1[70]
            +z[96]*tmp1[71]+z[55]*tmp1[72]+z[14]*tmp1[73]+z[70]*tmp1[74]+z[29]*tmp1[75]+z[85]*tmp1[76]+z[44]*tmp1[77]+z[3]*tmp1[78]+z[59]*tmp1[79]+z[18]*tmp1[80]
            +z[74]*tmp1[81]+z[33]*tmp1[82]+z[89]*tmp1[83]+z[48]*tmp1[84]+z[7]*tmp1[85]+z[63]*tmp1[86]+z[22]*tmp1[87]+z[78]*tmp1[88]+z[37]*tmp1[89]+z[93]*tmp1[90]
            +z[52]*tmp1[91]+z[11]*tmp1[92]+z[67]*tmp1[93]+z[26]*tmp1[94]+z[82]*tmp1[95]+z[41]*tmp1[96];
            tab[i+57*stg_first]=tab[i+57*stg_first]+z[0]*tmp1[0]
            +z[57]*tmp1[1]+z[17]*tmp1[2]+z[74]*tmp1[3]+z[34]*tmp1[4]+z[91]*tmp1[5]+z[51]*tmp1[6]+z[11]*tmp1[7]+z[68]*tmp1[8]+z[28]*tmp1[9]+z[85]*tmp1[10]
            +z[45]*tmp1[11]+z[5]*tmp1[12]+z[62]*tmp1[13]+z[22]*tmp1[14]+z[79]*tmp1[15]+z[39]*tmp1[16]+z[96]*tmp1[17]+z[56]*tmp1[18]+z[16]*tmp1[19]+z[73]*tmp1[20]
            +z[33]*tmp1[21]+z[90]*tmp1[22]+z[50]*tmp1[23]+z[10]*tmp1[24]+z[67]*tmp1[25]+z[27]*tmp1[26]+z[84]*tmp1[27]+z[44]*tmp1[28]+z[4]*tmp1[29]+z[61]*tmp1[30]
            +z[21]*tmp1[31]+z[78]*tmp1[32]+z[38]*tmp1[33]+z[95]*tmp1[34]+z[55]*tmp1[35]+z[15]*tmp1[36]+z[72]*tmp1[37]+z[32]*tmp1[38]+z[89]*tmp1[39]+z[49]*tmp1[40]
            +z[9]*tmp1[41]+z[66]*tmp1[42]+z[26]*tmp1[43]+z[83]*tmp1[44]+z[43]*tmp1[45]+z[3]*tmp1[46]+z[60]*tmp1[47]+z[20]*tmp1[48]+z[77]*tmp1[49]+z[37]*tmp1[50]
            +z[94]*tmp1[51]+z[54]*tmp1[52]+z[14]*tmp1[53]+z[71]*tmp1[54]+z[31]*tmp1[55]+z[88]*tmp1[56]+z[48]*tmp1[57]+z[8]*tmp1[58]+z[65]*tmp1[59]+z[25]*tmp1[60]
            +z[82]*tmp1[61]+z[42]*tmp1[62]+z[2]*tmp1[63]+z[59]*tmp1[64]+z[19]*tmp1[65]+z[76]*tmp1[66]+z[36]*tmp1[67]+z[93]*tmp1[68]+z[53]*tmp1[69]+z[13]*tmp1[70]
            +z[70]*tmp1[71]+z[30]*tmp1[72]+z[87]*tmp1[73]+z[47]*tmp1[74]+z[7]*tmp1[75]+z[64]*tmp1[76]+z[24]*tmp1[77]+z[81]*tmp1[78]+z[41]*tmp1[79]+z[1]*tmp1[80]
            +z[58]*tmp1[81]+z[18]*tmp1[82]+z[75]*tmp1[83]+z[35]*tmp1[84]+z[92]*tmp1[85]+z[52]*tmp1[86]+z[12]*tmp1[87]+z[69]*tmp1[88]+z[29]*tmp1[89]+z[86]*tmp1[90]
            +z[46]*tmp1[91]+z[6]*tmp1[92]+z[63]*tmp1[93]+z[23]*tmp1[94]+z[80]*tmp1[95]+z[40]*tmp1[96];
            tab[i+58*stg_first]=tab[i+58*stg_first]+z[0]*tmp1[0]
            +z[58]*tmp1[1]+z[19]*tmp1[2]+z[77]*tmp1[3]+z[38]*tmp1[4]+z[96]*tmp1[5]+z[57]*tmp1[6]+z[18]*tmp1[7]+z[76]*tmp1[8]+z[37]*tmp1[9]+z[95]*tmp1[10]
            +z[56]*tmp1[11]+z[17]*tmp1[12]+z[75]*tmp1[13]+z[36]*tmp1[14]+z[94]*tmp1[15]+z[55]*tmp1[16]+z[16]*tmp1[17]+z[74]*tmp1[18]+z[35]*tmp1[19]+z[93]*tmp1[20]
            +z[54]*tmp1[21]+z[15]*tmp1[22]+z[73]*tmp1[23]+z[34]*tmp1[24]+z[92]*tmp1[25]+z[53]*tmp1[26]+z[14]*tmp1[27]+z[72]*tmp1[28]+z[33]*tmp1[29]+z[91]*tmp1[30]
            +z[52]*tmp1[31]+z[13]*tmp1[32]+z[71]*tmp1[33]+z[32]*tmp1[34]+z[90]*tmp1[35]+z[51]*tmp1[36]+z[12]*tmp1[37]+z[70]*tmp1[38]+z[31]*tmp1[39]+z[89]*tmp1[40]
            +z[50]*tmp1[41]+z[11]*tmp1[42]+z[69]*tmp1[43]+z[30]*tmp1[44]+z[88]*tmp1[45]+z[49]*tmp1[46]+z[10]*tmp1[47]+z[68]*tmp1[48]+z[29]*tmp1[49]+z[87]*tmp1[50]
            +z[48]*tmp1[51]+z[9]*tmp1[52]+z[67]*tmp1[53]+z[28]*tmp1[54]+z[86]*tmp1[55]+z[47]*tmp1[56]+z[8]*tmp1[57]+z[66]*tmp1[58]+z[27]*tmp1[59]+z[85]*tmp1[60]
            +z[46]*tmp1[61]+z[7]*tmp1[62]+z[65]*tmp1[63]+z[26]*tmp1[64]+z[84]*tmp1[65]+z[45]*tmp1[66]+z[6]*tmp1[67]+z[64]*tmp1[68]+z[25]*tmp1[69]+z[83]*tmp1[70]
            +z[44]*tmp1[71]+z[5]*tmp1[72]+z[63]*tmp1[73]+z[24]*tmp1[74]+z[82]*tmp1[75]+z[43]*tmp1[76]+z[4]*tmp1[77]+z[62]*tmp1[78]+z[23]*tmp1[79]+z[81]*tmp1[80]
            +z[42]*tmp1[81]+z[3]*tmp1[82]+z[61]*tmp1[83]+z[22]*tmp1[84]+z[80]*tmp1[85]+z[41]*tmp1[86]+z[2]*tmp1[87]+z[60]*tmp1[88]+z[21]*tmp1[89]+z[79]*tmp1[90]
            +z[40]*tmp1[91]+z[1]*tmp1[92]+z[59]*tmp1[93]+z[20]*tmp1[94]+z[78]*tmp1[95]+z[39]*tmp1[96];
            tab[i+59*stg_first]=tab[i+59*stg_first]+z[0]*tmp1[0]
            +z[59]*tmp1[1]+z[21]*tmp1[2]+z[80]*tmp1[3]+z[42]*tmp1[4]+z[4]*tmp1[5]+z[63]*tmp1[6]+z[25]*tmp1[7]+z[84]*tmp1[8]+z[46]*tmp1[9]+z[8]*tmp1[10]
            +z[67]*tmp1[11]+z[29]*tmp1[12]+z[88]*tmp1[13]+z[50]*tmp1[14]+z[12]*tmp1[15]+z[71]*tmp1[16]+z[33]*tmp1[17]+z[92]*tmp1[18]+z[54]*tmp1[19]+z[16]*tmp1[20]
            +z[75]*tmp1[21]+z[37]*tmp1[22]+z[96]*tmp1[23]+z[58]*tmp1[24]+z[20]*tmp1[25]+z[79]*tmp1[26]+z[41]*tmp1[27]+z[3]*tmp1[28]+z[62]*tmp1[29]+z[24]*tmp1[30]
            +z[83]*tmp1[31]+z[45]*tmp1[32]+z[7]*tmp1[33]+z[66]*tmp1[34]+z[28]*tmp1[35]+z[87]*tmp1[36]+z[49]*tmp1[37]+z[11]*tmp1[38]+z[70]*tmp1[39]+z[32]*tmp1[40]
            +z[91]*tmp1[41]+z[53]*tmp1[42]+z[15]*tmp1[43]+z[74]*tmp1[44]+z[36]*tmp1[45]+z[95]*tmp1[46]+z[57]*tmp1[47]+z[19]*tmp1[48]+z[78]*tmp1[49]+z[40]*tmp1[50]
            +z[2]*tmp1[51]+z[61]*tmp1[52]+z[23]*tmp1[53]+z[82]*tmp1[54]+z[44]*tmp1[55]+z[6]*tmp1[56]+z[65]*tmp1[57]+z[27]*tmp1[58]+z[86]*tmp1[59]+z[48]*tmp1[60]
            +z[10]*tmp1[61]+z[69]*tmp1[62]+z[31]*tmp1[63]+z[90]*tmp1[64]+z[52]*tmp1[65]+z[14]*tmp1[66]+z[73]*tmp1[67]+z[35]*tmp1[68]+z[94]*tmp1[69]+z[56]*tmp1[70]
            +z[18]*tmp1[71]+z[77]*tmp1[72]+z[39]*tmp1[73]+z[1]*tmp1[74]+z[60]*tmp1[75]+z[22]*tmp1[76]+z[81]*tmp1[77]+z[43]*tmp1[78]+z[5]*tmp1[79]+z[64]*tmp1[80]
            +z[26]*tmp1[81]+z[85]*tmp1[82]+z[47]*tmp1[83]+z[9]*tmp1[84]+z[68]*tmp1[85]+z[30]*tmp1[86]+z[89]*tmp1[87]+z[51]*tmp1[88]+z[13]*tmp1[89]+z[72]*tmp1[90]
            +z[34]*tmp1[91]+z[93]*tmp1[92]+z[55]*tmp1[93]+z[17]*tmp1[94]+z[76]*tmp1[95]+z[38]*tmp1[96];
            tab[i+60*stg_first]=tab[i+60*stg_first]+z[0]*tmp1[0]
            +z[60]*tmp1[1]+z[23]*tmp1[2]+z[83]*tmp1[3]+z[46]*tmp1[4]+z[9]*tmp1[5]+z[69]*tmp1[6]+z[32]*tmp1[7]+z[92]*tmp1[8]+z[55]*tmp1[9]+z[18]*tmp1[10]
            +z[78]*tmp1[11]+z[41]*tmp1[12]+z[4]*tmp1[13]+z[64]*tmp1[14]+z[27]*tmp1[15]+z[87]*tmp1[16]+z[50]*tmp1[17]+z[13]*tmp1[18]+z[73]*tmp1[19]+z[36]*tmp1[20]
            +z[96]*tmp1[21]+z[59]*tmp1[22]+z[22]*tmp1[23]+z[82]*tmp1[24]+z[45]*tmp1[25]+z[8]*tmp1[26]+z[68]*tmp1[27]+z[31]*tmp1[28]+z[91]*tmp1[29]+z[54]*tmp1[30]
            +z[17]*tmp1[31]+z[77]*tmp1[32]+z[40]*tmp1[33]+z[3]*tmp1[34]+z[63]*tmp1[35]+z[26]*tmp1[36]+z[86]*tmp1[37]+z[49]*tmp1[38]+z[12]*tmp1[39]+z[72]*tmp1[40]
            +z[35]*tmp1[41]+z[95]*tmp1[42]+z[58]*tmp1[43]+z[21]*tmp1[44]+z[81]*tmp1[45]+z[44]*tmp1[46]+z[7]*tmp1[47]+z[67]*tmp1[48]+z[30]*tmp1[49]+z[90]*tmp1[50]
            +z[53]*tmp1[51]+z[16]*tmp1[52]+z[76]*tmp1[53]+z[39]*tmp1[54]+z[2]*tmp1[55]+z[62]*tmp1[56]+z[25]*tmp1[57]+z[85]*tmp1[58]+z[48]*tmp1[59]+z[11]*tmp1[60]
            +z[71]*tmp1[61]+z[34]*tmp1[62]+z[94]*tmp1[63]+z[57]*tmp1[64]+z[20]*tmp1[65]+z[80]*tmp1[66]+z[43]*tmp1[67]+z[6]*tmp1[68]+z[66]*tmp1[69]+z[29]*tmp1[70]
            +z[89]*tmp1[71]+z[52]*tmp1[72]+z[15]*tmp1[73]+z[75]*tmp1[74]+z[38]*tmp1[75]+z[1]*tmp1[76]+z[61]*tmp1[77]+z[24]*tmp1[78]+z[84]*tmp1[79]+z[47]*tmp1[80]
            +z[10]*tmp1[81]+z[70]*tmp1[82]+z[33]*tmp1[83]+z[93]*tmp1[84]+z[56]*tmp1[85]+z[19]*tmp1[86]+z[79]*tmp1[87]+z[42]*tmp1[88]+z[5]*tmp1[89]+z[65]*tmp1[90]
            +z[28]*tmp1[91]+z[88]*tmp1[92]+z[51]*tmp1[93]+z[14]*tmp1[94]+z[74]*tmp1[95]+z[37]*tmp1[96];
            tab[i+61*stg_first]=tab[i+61*stg_first]+z[0]*tmp1[0]
            +z[61]*tmp1[1]+z[25]*tmp1[2]+z[86]*tmp1[3]+z[50]*tmp1[4]+z[14]*tmp1[5]+z[75]*tmp1[6]+z[39]*tmp1[7]+z[3]*tmp1[8]+z[64]*tmp1[9]+z[28]*tmp1[10]
            +z[89]*tmp1[11]+z[53]*tmp1[12]+z[17]*tmp1[13]+z[78]*tmp1[14]+z[42]*tmp1[15]+z[6]*tmp1[16]+z[67]*tmp1[17]+z[31]*tmp1[18]+z[92]*tmp1[19]+z[56]*tmp1[20]
            +z[20]*tmp1[21]+z[81]*tmp1[22]+z[45]*tmp1[23]+z[9]*tmp1[24]+z[70]*tmp1[25]+z[34]*tmp1[26]+z[95]*tmp1[27]+z[59]*tmp1[28]+z[23]*tmp1[29]+z[84]*tmp1[30]
            +z[48]*tmp1[31]+z[12]*tmp1[32]+z[73]*tmp1[33]+z[37]*tmp1[34]+z[1]*tmp1[35]+z[62]*tmp1[36]+z[26]*tmp1[37]+z[87]*tmp1[38]+z[51]*tmp1[39]+z[15]*tmp1[40]
            +z[76]*tmp1[41]+z[40]*tmp1[42]+z[4]*tmp1[43]+z[65]*tmp1[44]+z[29]*tmp1[45]+z[90]*tmp1[46]+z[54]*tmp1[47]+z[18]*tmp1[48]+z[79]*tmp1[49]+z[43]*tmp1[50]
            +z[7]*tmp1[51]+z[68]*tmp1[52]+z[32]*tmp1[53]+z[93]*tmp1[54]+z[57]*tmp1[55]+z[21]*tmp1[56]+z[82]*tmp1[57]+z[46]*tmp1[58]+z[10]*tmp1[59]+z[71]*tmp1[60]
            +z[35]*tmp1[61]+z[96]*tmp1[62]+z[60]*tmp1[63]+z[24]*tmp1[64]+z[85]*tmp1[65]+z[49]*tmp1[66]+z[13]*tmp1[67]+z[74]*tmp1[68]+z[38]*tmp1[69]+z[2]*tmp1[70]
            +z[63]*tmp1[71]+z[27]*tmp1[72]+z[88]*tmp1[73]+z[52]*tmp1[74]+z[16]*tmp1[75]+z[77]*tmp1[76]+z[41]*tmp1[77]+z[5]*tmp1[78]+z[66]*tmp1[79]+z[30]*tmp1[80]
            +z[91]*tmp1[81]+z[55]*tmp1[82]+z[19]*tmp1[83]+z[80]*tmp1[84]+z[44]*tmp1[85]+z[8]*tmp1[86]+z[69]*tmp1[87]+z[33]*tmp1[88]+z[94]*tmp1[89]+z[58]*tmp1[90]
            +z[22]*tmp1[91]+z[83]*tmp1[92]+z[47]*tmp1[93]+z[11]*tmp1[94]+z[72]*tmp1[95]+z[36]*tmp1[96];
            tab[i+62*stg_first]=tab[i+62*stg_first]+z[0]*tmp1[0]
            +z[62]*tmp1[1]+z[27]*tmp1[2]+z[89]*tmp1[3]+z[54]*tmp1[4]+z[19]*tmp1[5]+z[81]*tmp1[6]+z[46]*tmp1[7]+z[11]*tmp1[8]+z[73]*tmp1[9]+z[38]*tmp1[10]
            +z[3]*tmp1[11]+z[65]*tmp1[12]+z[30]*tmp1[13]+z[92]*tmp1[14]+z[57]*tmp1[15]+z[22]*tmp1[16]+z[84]*tmp1[17]+z[49]*tmp1[18]+z[14]*tmp1[19]+z[76]*tmp1[20]
            +z[41]*tmp1[21]+z[6]*tmp1[22]+z[68]*tmp1[23]+z[33]*tmp1[24]+z[95]*tmp1[25]+z[60]*tmp1[26]+z[25]*tmp1[27]+z[87]*tmp1[28]+z[52]*tmp1[29]+z[17]*tmp1[30]
            +z[79]*tmp1[31]+z[44]*tmp1[32]+z[9]*tmp1[33]+z[71]*tmp1[34]+z[36]*tmp1[35]+z[1]*tmp1[36]+z[63]*tmp1[37]+z[28]*tmp1[38]+z[90]*tmp1[39]+z[55]*tmp1[40]
            +z[20]*tmp1[41]+z[82]*tmp1[42]+z[47]*tmp1[43]+z[12]*tmp1[44]+z[74]*tmp1[45]+z[39]*tmp1[46]+z[4]*tmp1[47]+z[66]*tmp1[48]+z[31]*tmp1[49]+z[93]*tmp1[50]
            +z[58]*tmp1[51]+z[23]*tmp1[52]+z[85]*tmp1[53]+z[50]*tmp1[54]+z[15]*tmp1[55]+z[77]*tmp1[56]+z[42]*tmp1[57]+z[7]*tmp1[58]+z[69]*tmp1[59]+z[34]*tmp1[60]
            +z[96]*tmp1[61]+z[61]*tmp1[62]+z[26]*tmp1[63]+z[88]*tmp1[64]+z[53]*tmp1[65]+z[18]*tmp1[66]+z[80]*tmp1[67]+z[45]*tmp1[68]+z[10]*tmp1[69]+z[72]*tmp1[70]
            +z[37]*tmp1[71]+z[2]*tmp1[72]+z[64]*tmp1[73]+z[29]*tmp1[74]+z[91]*tmp1[75]+z[56]*tmp1[76]+z[21]*tmp1[77]+z[83]*tmp1[78]+z[48]*tmp1[79]+z[13]*tmp1[80]
            +z[75]*tmp1[81]+z[40]*tmp1[82]+z[5]*tmp1[83]+z[67]*tmp1[84]+z[32]*tmp1[85]+z[94]*tmp1[86]+z[59]*tmp1[87]+z[24]*tmp1[88]+z[86]*tmp1[89]+z[51]*tmp1[90]
            +z[16]*tmp1[91]+z[78]*tmp1[92]+z[43]*tmp1[93]+z[8]*tmp1[94]+z[70]*tmp1[95]+z[35]*tmp1[96];
            tab[i+63*stg_first]=tab[i+63*stg_first]+z[0]*tmp1[0]
            +z[63]*tmp1[1]+z[29]*tmp1[2]+z[92]*tmp1[3]+z[58]*tmp1[4]+z[24]*tmp1[5]+z[87]*tmp1[6]+z[53]*tmp1[7]+z[19]*tmp1[8]+z[82]*tmp1[9]+z[48]*tmp1[10]
            +z[14]*tmp1[11]+z[77]*tmp1[12]+z[43]*tmp1[13]+z[9]*tmp1[14]+z[72]*tmp1[15]+z[38]*tmp1[16]+z[4]*tmp1[17]+z[67]*tmp1[18]+z[33]*tmp1[19]+z[96]*tmp1[20]
            +z[62]*tmp1[21]+z[28]*tmp1[22]+z[91]*tmp1[23]+z[57]*tmp1[24]+z[23]*tmp1[25]+z[86]*tmp1[26]+z[52]*tmp1[27]+z[18]*tmp1[28]+z[81]*tmp1[29]+z[47]*tmp1[30]
            +z[13]*tmp1[31]+z[76]*tmp1[32]+z[42]*tmp1[33]+z[8]*tmp1[34]+z[71]*tmp1[35]+z[37]*tmp1[36]+z[3]*tmp1[37]+z[66]*tmp1[38]+z[32]*tmp1[39]+z[95]*tmp1[40]
            +z[61]*tmp1[41]+z[27]*tmp1[42]+z[90]*tmp1[43]+z[56]*tmp1[44]+z[22]*tmp1[45]+z[85]*tmp1[46]+z[51]*tmp1[47]+z[17]*tmp1[48]+z[80]*tmp1[49]+z[46]*tmp1[50]
            +z[12]*tmp1[51]+z[75]*tmp1[52]+z[41]*tmp1[53]+z[7]*tmp1[54]+z[70]*tmp1[55]+z[36]*tmp1[56]+z[2]*tmp1[57]+z[65]*tmp1[58]+z[31]*tmp1[59]+z[94]*tmp1[60]
            +z[60]*tmp1[61]+z[26]*tmp1[62]+z[89]*tmp1[63]+z[55]*tmp1[64]+z[21]*tmp1[65]+z[84]*tmp1[66]+z[50]*tmp1[67]+z[16]*tmp1[68]+z[79]*tmp1[69]+z[45]*tmp1[70]
            +z[11]*tmp1[71]+z[74]*tmp1[72]+z[40]*tmp1[73]+z[6]*tmp1[74]+z[69]*tmp1[75]+z[35]*tmp1[76]+z[1]*tmp1[77]+z[64]*tmp1[78]+z[30]*tmp1[79]+z[93]*tmp1[80]
            +z[59]*tmp1[81]+z[25]*tmp1[82]+z[88]*tmp1[83]+z[54]*tmp1[84]+z[20]*tmp1[85]+z[83]*tmp1[86]+z[49]*tmp1[87]+z[15]*tmp1[88]+z[78]*tmp1[89]+z[44]*tmp1[90]
            +z[10]*tmp1[91]+z[73]*tmp1[92]+z[39]*tmp1[93]+z[5]*tmp1[94]+z[68]*tmp1[95]+z[34]*tmp1[96];
            tab[i+64*stg_first]=tab[i+64*stg_first]+z[0]*tmp1[0]
            +z[64]*tmp1[1]+z[31]*tmp1[2]+z[95]*tmp1[3]+z[62]*tmp1[4]+z[29]*tmp1[5]+z[93]*tmp1[6]+z[60]*tmp1[7]+z[27]*tmp1[8]+z[91]*tmp1[9]+z[58]*tmp1[10]
            +z[25]*tmp1[11]+z[89]*tmp1[12]+z[56]*tmp1[13]+z[23]*tmp1[14]+z[87]*tmp1[15]+z[54]*tmp1[16]+z[21]*tmp1[17]+z[85]*tmp1[18]+z[52]*tmp1[19]+z[19]*tmp1[20]
            +z[83]*tmp1[21]+z[50]*tmp1[22]+z[17]*tmp1[23]+z[81]*tmp1[24]+z[48]*tmp1[25]+z[15]*tmp1[26]+z[79]*tmp1[27]+z[46]*tmp1[28]+z[13]*tmp1[29]+z[77]*tmp1[30]
            +z[44]*tmp1[31]+z[11]*tmp1[32]+z[75]*tmp1[33]+z[42]*tmp1[34]+z[9]*tmp1[35]+z[73]*tmp1[36]+z[40]*tmp1[37]+z[7]*tmp1[38]+z[71]*tmp1[39]+z[38]*tmp1[40]
            +z[5]*tmp1[41]+z[69]*tmp1[42]+z[36]*tmp1[43]+z[3]*tmp1[44]+z[67]*tmp1[45]+z[34]*tmp1[46]+z[1]*tmp1[47]+z[65]*tmp1[48]+z[32]*tmp1[49]+z[96]*tmp1[50]
            +z[63]*tmp1[51]+z[30]*tmp1[52]+z[94]*tmp1[53]+z[61]*tmp1[54]+z[28]*tmp1[55]+z[92]*tmp1[56]+z[59]*tmp1[57]+z[26]*tmp1[58]+z[90]*tmp1[59]+z[57]*tmp1[60]
            +z[24]*tmp1[61]+z[88]*tmp1[62]+z[55]*tmp1[63]+z[22]*tmp1[64]+z[86]*tmp1[65]+z[53]*tmp1[66]+z[20]*tmp1[67]+z[84]*tmp1[68]+z[51]*tmp1[69]+z[18]*tmp1[70]
            +z[82]*tmp1[71]+z[49]*tmp1[72]+z[16]*tmp1[73]+z[80]*tmp1[74]+z[47]*tmp1[75]+z[14]*tmp1[76]+z[78]*tmp1[77]+z[45]*tmp1[78]+z[12]*tmp1[79]+z[76]*tmp1[80]
            +z[43]*tmp1[81]+z[10]*tmp1[82]+z[74]*tmp1[83]+z[41]*tmp1[84]+z[8]*tmp1[85]+z[72]*tmp1[86]+z[39]*tmp1[87]+z[6]*tmp1[88]+z[70]*tmp1[89]+z[37]*tmp1[90]
            +z[4]*tmp1[91]+z[68]*tmp1[92]+z[35]*tmp1[93]+z[2]*tmp1[94]+z[66]*tmp1[95]+z[33]*tmp1[96];
            tab[i+65*stg_first]=tab[i+65*stg_first]+z[0]*tmp1[0]
            +z[65]*tmp1[1]+z[33]*tmp1[2]+z[1]*tmp1[3]+z[66]*tmp1[4]+z[34]*tmp1[5]+z[2]*tmp1[6]+z[67]*tmp1[7]+z[35]*tmp1[8]+z[3]*tmp1[9]+z[68]*tmp1[10]
            +z[36]*tmp1[11]+z[4]*tmp1[12]+z[69]*tmp1[13]+z[37]*tmp1[14]+z[5]*tmp1[15]+z[70]*tmp1[16]+z[38]*tmp1[17]+z[6]*tmp1[18]+z[71]*tmp1[19]+z[39]*tmp1[20]
            +z[7]*tmp1[21]+z[72]*tmp1[22]+z[40]*tmp1[23]+z[8]*tmp1[24]+z[73]*tmp1[25]+z[41]*tmp1[26]+z[9]*tmp1[27]+z[74]*tmp1[28]+z[42]*tmp1[29]+z[10]*tmp1[30]
            +z[75]*tmp1[31]+z[43]*tmp1[32]+z[11]*tmp1[33]+z[76]*tmp1[34]+z[44]*tmp1[35]+z[12]*tmp1[36]+z[77]*tmp1[37]+z[45]*tmp1[38]+z[13]*tmp1[39]+z[78]*tmp1[40]
            +z[46]*tmp1[41]+z[14]*tmp1[42]+z[79]*tmp1[43]+z[47]*tmp1[44]+z[15]*tmp1[45]+z[80]*tmp1[46]+z[48]*tmp1[47]+z[16]*tmp1[48]+z[81]*tmp1[49]+z[49]*tmp1[50]
            +z[17]*tmp1[51]+z[82]*tmp1[52]+z[50]*tmp1[53]+z[18]*tmp1[54]+z[83]*tmp1[55]+z[51]*tmp1[56]+z[19]*tmp1[57]+z[84]*tmp1[58]+z[52]*tmp1[59]+z[20]*tmp1[60]
            +z[85]*tmp1[61]+z[53]*tmp1[62]+z[21]*tmp1[63]+z[86]*tmp1[64]+z[54]*tmp1[65]+z[22]*tmp1[66]+z[87]*tmp1[67]+z[55]*tmp1[68]+z[23]*tmp1[69]+z[88]*tmp1[70]
            +z[56]*tmp1[71]+z[24]*tmp1[72]+z[89]*tmp1[73]+z[57]*tmp1[74]+z[25]*tmp1[75]+z[90]*tmp1[76]+z[58]*tmp1[77]+z[26]*tmp1[78]+z[91]*tmp1[79]+z[59]*tmp1[80]
            +z[27]*tmp1[81]+z[92]*tmp1[82]+z[60]*tmp1[83]+z[28]*tmp1[84]+z[93]*tmp1[85]+z[61]*tmp1[86]+z[29]*tmp1[87]+z[94]*tmp1[88]+z[62]*tmp1[89]+z[30]*tmp1[90]
            +z[95]*tmp1[91]+z[63]*tmp1[92]+z[31]*tmp1[93]+z[96]*tmp1[94]+z[64]*tmp1[95]+z[32]*tmp1[96];
            tab[i+66*stg_first]=tab[i+66*stg_first]+z[0]*tmp1[0]
            +z[66]*tmp1[1]+z[35]*tmp1[2]+z[4]*tmp1[3]+z[70]*tmp1[4]+z[39]*tmp1[5]+z[8]*tmp1[6]+z[74]*tmp1[7]+z[43]*tmp1[8]+z[12]*tmp1[9]+z[78]*tmp1[10]
            +z[47]*tmp1[11]+z[16]*tmp1[12]+z[82]*tmp1[13]+z[51]*tmp1[14]+z[20]*tmp1[15]+z[86]*tmp1[16]+z[55]*tmp1[17]+z[24]*tmp1[18]+z[90]*tmp1[19]+z[59]*tmp1[20]
            +z[28]*tmp1[21]+z[94]*tmp1[22]+z[63]*tmp1[23]+z[32]*tmp1[24]+z[1]*tmp1[25]+z[67]*tmp1[26]+z[36]*tmp1[27]+z[5]*tmp1[28]+z[71]*tmp1[29]+z[40]*tmp1[30]
            +z[9]*tmp1[31]+z[75]*tmp1[32]+z[44]*tmp1[33]+z[13]*tmp1[34]+z[79]*tmp1[35]+z[48]*tmp1[36]+z[17]*tmp1[37]+z[83]*tmp1[38]+z[52]*tmp1[39]+z[21]*tmp1[40]
            +z[87]*tmp1[41]+z[56]*tmp1[42]+z[25]*tmp1[43]+z[91]*tmp1[44]+z[60]*tmp1[45]+z[29]*tmp1[46]+z[95]*tmp1[47]+z[64]*tmp1[48]+z[33]*tmp1[49]+z[2]*tmp1[50]
            +z[68]*tmp1[51]+z[37]*tmp1[52]+z[6]*tmp1[53]+z[72]*tmp1[54]+z[41]*tmp1[55]+z[10]*tmp1[56]+z[76]*tmp1[57]+z[45]*tmp1[58]+z[14]*tmp1[59]+z[80]*tmp1[60]
            +z[49]*tmp1[61]+z[18]*tmp1[62]+z[84]*tmp1[63]+z[53]*tmp1[64]+z[22]*tmp1[65]+z[88]*tmp1[66]+z[57]*tmp1[67]+z[26]*tmp1[68]+z[92]*tmp1[69]+z[61]*tmp1[70]
            +z[30]*tmp1[71]+z[96]*tmp1[72]+z[65]*tmp1[73]+z[34]*tmp1[74]+z[3]*tmp1[75]+z[69]*tmp1[76]+z[38]*tmp1[77]+z[7]*tmp1[78]+z[73]*tmp1[79]+z[42]*tmp1[80]
            +z[11]*tmp1[81]+z[77]*tmp1[82]+z[46]*tmp1[83]+z[15]*tmp1[84]+z[81]*tmp1[85]+z[50]*tmp1[86]+z[19]*tmp1[87]+z[85]*tmp1[88]+z[54]*tmp1[89]+z[23]*tmp1[90]
            +z[89]*tmp1[91]+z[58]*tmp1[92]+z[27]*tmp1[93]+z[93]*tmp1[94]+z[62]*tmp1[95]+z[31]*tmp1[96];
            tab[i+67*stg_first]=tab[i+67*stg_first]+z[0]*tmp1[0]
            +z[67]*tmp1[1]+z[37]*tmp1[2]+z[7]*tmp1[3]+z[74]*tmp1[4]+z[44]*tmp1[5]+z[14]*tmp1[6]+z[81]*tmp1[7]+z[51]*tmp1[8]+z[21]*tmp1[9]+z[88]*tmp1[10]
            +z[58]*tmp1[11]+z[28]*tmp1[12]+z[95]*tmp1[13]+z[65]*tmp1[14]+z[35]*tmp1[15]+z[5]*tmp1[16]+z[72]*tmp1[17]+z[42]*tmp1[18]+z[12]*tmp1[19]+z[79]*tmp1[20]
            +z[49]*tmp1[21]+z[19]*tmp1[22]+z[86]*tmp1[23]+z[56]*tmp1[24]+z[26]*tmp1[25]+z[93]*tmp1[26]+z[63]*tmp1[27]+z[33]*tmp1[28]+z[3]*tmp1[29]+z[70]*tmp1[30]
            +z[40]*tmp1[31]+z[10]*tmp1[32]+z[77]*tmp1[33]+z[47]*tmp1[34]+z[17]*tmp1[35]+z[84]*tmp1[36]+z[54]*tmp1[37]+z[24]*tmp1[38]+z[91]*tmp1[39]+z[61]*tmp1[40]
            +z[31]*tmp1[41]+z[1]*tmp1[42]+z[68]*tmp1[43]+z[38]*tmp1[44]+z[8]*tmp1[45]+z[75]*tmp1[46]+z[45]*tmp1[47]+z[15]*tmp1[48]+z[82]*tmp1[49]+z[52]*tmp1[50]
            +z[22]*tmp1[51]+z[89]*tmp1[52]+z[59]*tmp1[53]+z[29]*tmp1[54]+z[96]*tmp1[55]+z[66]*tmp1[56]+z[36]*tmp1[57]+z[6]*tmp1[58]+z[73]*tmp1[59]+z[43]*tmp1[60]
            +z[13]*tmp1[61]+z[80]*tmp1[62]+z[50]*tmp1[63]+z[20]*tmp1[64]+z[87]*tmp1[65]+z[57]*tmp1[66]+z[27]*tmp1[67]+z[94]*tmp1[68]+z[64]*tmp1[69]+z[34]*tmp1[70]
            +z[4]*tmp1[71]+z[71]*tmp1[72]+z[41]*tmp1[73]+z[11]*tmp1[74]+z[78]*tmp1[75]+z[48]*tmp1[76]+z[18]*tmp1[77]+z[85]*tmp1[78]+z[55]*tmp1[79]+z[25]*tmp1[80]
            +z[92]*tmp1[81]+z[62]*tmp1[82]+z[32]*tmp1[83]+z[2]*tmp1[84]+z[69]*tmp1[85]+z[39]*tmp1[86]+z[9]*tmp1[87]+z[76]*tmp1[88]+z[46]*tmp1[89]+z[16]*tmp1[90]
            +z[83]*tmp1[91]+z[53]*tmp1[92]+z[23]*tmp1[93]+z[90]*tmp1[94]+z[60]*tmp1[95]+z[30]*tmp1[96];
            tab[i+68*stg_first]=tab[i+68*stg_first]+z[0]*tmp1[0]
            +z[68]*tmp1[1]+z[39]*tmp1[2]+z[10]*tmp1[3]+z[78]*tmp1[4]+z[49]*tmp1[5]+z[20]*tmp1[6]+z[88]*tmp1[7]+z[59]*tmp1[8]+z[30]*tmp1[9]+z[1]*tmp1[10]
            +z[69]*tmp1[11]+z[40]*tmp1[12]+z[11]*tmp1[13]+z[79]*tmp1[14]+z[50]*tmp1[15]+z[21]*tmp1[16]+z[89]*tmp1[17]+z[60]*tmp1[18]+z[31]*tmp1[19]+z[2]*tmp1[20]
            +z[70]*tmp1[21]+z[41]*tmp1[22]+z[12]*tmp1[23]+z[80]*tmp1[24]+z[51]*tmp1[25]+z[22]*tmp1[26]+z[90]*tmp1[27]+z[61]*tmp1[28]+z[32]*tmp1[29]+z[3]*tmp1[30]
            +z[71]*tmp1[31]+z[42]*tmp1[32]+z[13]*tmp1[33]+z[81]*tmp1[34]+z[52]*tmp1[35]+z[23]*tmp1[36]+z[91]*tmp1[37]+z[62]*tmp1[38]+z[33]*tmp1[39]+z[4]*tmp1[40]
            +z[72]*tmp1[41]+z[43]*tmp1[42]+z[14]*tmp1[43]+z[82]*tmp1[44]+z[53]*tmp1[45]+z[24]*tmp1[46]+z[92]*tmp1[47]+z[63]*tmp1[48]+z[34]*tmp1[49]+z[5]*tmp1[50]
            +z[73]*tmp1[51]+z[44]*tmp1[52]+z[15]*tmp1[53]+z[83]*tmp1[54]+z[54]*tmp1[55]+z[25]*tmp1[56]+z[93]*tmp1[57]+z[64]*tmp1[58]+z[35]*tmp1[59]+z[6]*tmp1[60]
            +z[74]*tmp1[61]+z[45]*tmp1[62]+z[16]*tmp1[63]+z[84]*tmp1[64]+z[55]*tmp1[65]+z[26]*tmp1[66]+z[94]*tmp1[67]+z[65]*tmp1[68]+z[36]*tmp1[69]+z[7]*tmp1[70]
            +z[75]*tmp1[71]+z[46]*tmp1[72]+z[17]*tmp1[73]+z[85]*tmp1[74]+z[56]*tmp1[75]+z[27]*tmp1[76]+z[95]*tmp1[77]+z[66]*tmp1[78]+z[37]*tmp1[79]+z[8]*tmp1[80]
            +z[76]*tmp1[81]+z[47]*tmp1[82]+z[18]*tmp1[83]+z[86]*tmp1[84]+z[57]*tmp1[85]+z[28]*tmp1[86]+z[96]*tmp1[87]+z[67]*tmp1[88]+z[38]*tmp1[89]+z[9]*tmp1[90]
            +z[77]*tmp1[91]+z[48]*tmp1[92]+z[19]*tmp1[93]+z[87]*tmp1[94]+z[58]*tmp1[95]+z[29]*tmp1[96];
            tab[i+69*stg_first]=tab[i+69*stg_first]+z[0]*tmp1[0]
            +z[69]*tmp1[1]+z[41]*tmp1[2]+z[13]*tmp1[3]+z[82]*tmp1[4]+z[54]*tmp1[5]+z[26]*tmp1[6]+z[95]*tmp1[7]+z[67]*tmp1[8]+z[39]*tmp1[9]+z[11]*tmp1[10]
            +z[80]*tmp1[11]+z[52]*tmp1[12]+z[24]*tmp1[13]+z[93]*tmp1[14]+z[65]*tmp1[15]+z[37]*tmp1[16]+z[9]*tmp1[17]+z[78]*tmp1[18]+z[50]*tmp1[19]+z[22]*tmp1[20]
            +z[91]*tmp1[21]+z[63]*tmp1[22]+z[35]*tmp1[23]+z[7]*tmp1[24]+z[76]*tmp1[25]+z[48]*tmp1[26]+z[20]*tmp1[27]+z[89]*tmp1[28]+z[61]*tmp1[29]+z[33]*tmp1[30]
            +z[5]*tmp1[31]+z[74]*tmp1[32]+z[46]*tmp1[33]+z[18]*tmp1[34]+z[87]*tmp1[35]+z[59]*tmp1[36]+z[31]*tmp1[37]+z[3]*tmp1[38]+z[72]*tmp1[39]+z[44]*tmp1[40]
            +z[16]*tmp1[41]+z[85]*tmp1[42]+z[57]*tmp1[43]+z[29]*tmp1[44]+z[1]*tmp1[45]+z[70]*tmp1[46]+z[42]*tmp1[47]+z[14]*tmp1[48]+z[83]*tmp1[49]+z[55]*tmp1[50]
            +z[27]*tmp1[51]+z[96]*tmp1[52]+z[68]*tmp1[53]+z[40]*tmp1[54]+z[12]*tmp1[55]+z[81]*tmp1[56]+z[53]*tmp1[57]+z[25]*tmp1[58]+z[94]*tmp1[59]+z[66]*tmp1[60]
            +z[38]*tmp1[61]+z[10]*tmp1[62]+z[79]*tmp1[63]+z[51]*tmp1[64]+z[23]*tmp1[65]+z[92]*tmp1[66]+z[64]*tmp1[67]+z[36]*tmp1[68]+z[8]*tmp1[69]+z[77]*tmp1[70]
            +z[49]*tmp1[71]+z[21]*tmp1[72]+z[90]*tmp1[73]+z[62]*tmp1[74]+z[34]*tmp1[75]+z[6]*tmp1[76]+z[75]*tmp1[77]+z[47]*tmp1[78]+z[19]*tmp1[79]+z[88]*tmp1[80]
            +z[60]*tmp1[81]+z[32]*tmp1[82]+z[4]*tmp1[83]+z[73]*tmp1[84]+z[45]*tmp1[85]+z[17]*tmp1[86]+z[86]*tmp1[87]+z[58]*tmp1[88]+z[30]*tmp1[89]+z[2]*tmp1[90]
            +z[71]*tmp1[91]+z[43]*tmp1[92]+z[15]*tmp1[93]+z[84]*tmp1[94]+z[56]*tmp1[95]+z[28]*tmp1[96];
            tab[i+70*stg_first]=tab[i+70*stg_first]+z[0]*tmp1[0]
            +z[70]*tmp1[1]+z[43]*tmp1[2]+z[16]*tmp1[3]+z[86]*tmp1[4]+z[59]*tmp1[5]+z[32]*tmp1[6]+z[5]*tmp1[7]+z[75]*tmp1[8]+z[48]*tmp1[9]+z[21]*tmp1[10]
            +z[91]*tmp1[11]+z[64]*tmp1[12]+z[37]*tmp1[13]+z[10]*tmp1[14]+z[80]*tmp1[15]+z[53]*tmp1[16]+z[26]*tmp1[17]+z[96]*tmp1[18]+z[69]*tmp1[19]+z[42]*tmp1[20]
            +z[15]*tmp1[21]+z[85]*tmp1[22]+z[58]*tmp1[23]+z[31]*tmp1[24]+z[4]*tmp1[25]+z[74]*tmp1[26]+z[47]*tmp1[27]+z[20]*tmp1[28]+z[90]*tmp1[29]+z[63]*tmp1[30]
            +z[36]*tmp1[31]+z[9]*tmp1[32]+z[79]*tmp1[33]+z[52]*tmp1[34]+z[25]*tmp1[35]+z[95]*tmp1[36]+z[68]*tmp1[37]+z[41]*tmp1[38]+z[14]*tmp1[39]+z[84]*tmp1[40]
            +z[57]*tmp1[41]+z[30]*tmp1[42]+z[3]*tmp1[43]+z[73]*tmp1[44]+z[46]*tmp1[45]+z[19]*tmp1[46]+z[89]*tmp1[47]+z[62]*tmp1[48]+z[35]*tmp1[49]+z[8]*tmp1[50]
            +z[78]*tmp1[51]+z[51]*tmp1[52]+z[24]*tmp1[53]+z[94]*tmp1[54]+z[67]*tmp1[55]+z[40]*tmp1[56]+z[13]*tmp1[57]+z[83]*tmp1[58]+z[56]*tmp1[59]+z[29]*tmp1[60]
            +z[2]*tmp1[61]+z[72]*tmp1[62]+z[45]*tmp1[63]+z[18]*tmp1[64]+z[88]*tmp1[65]+z[61]*tmp1[66]+z[34]*tmp1[67]+z[7]*tmp1[68]+z[77]*tmp1[69]+z[50]*tmp1[70]
            +z[23]*tmp1[71]+z[93]*tmp1[72]+z[66]*tmp1[73]+z[39]*tmp1[74]+z[12]*tmp1[75]+z[82]*tmp1[76]+z[55]*tmp1[77]+z[28]*tmp1[78]+z[1]*tmp1[79]+z[71]*tmp1[80]
            +z[44]*tmp1[81]+z[17]*tmp1[82]+z[87]*tmp1[83]+z[60]*tmp1[84]+z[33]*tmp1[85]+z[6]*tmp1[86]+z[76]*tmp1[87]+z[49]*tmp1[88]+z[22]*tmp1[89]+z[92]*tmp1[90]
            +z[65]*tmp1[91]+z[38]*tmp1[92]+z[11]*tmp1[93]+z[81]*tmp1[94]+z[54]*tmp1[95]+z[27]*tmp1[96];
            tab[i+71*stg_first]=tab[i+71*stg_first]+z[0]*tmp1[0]
            +z[71]*tmp1[1]+z[45]*tmp1[2]+z[19]*tmp1[3]+z[90]*tmp1[4]+z[64]*tmp1[5]+z[38]*tmp1[6]+z[12]*tmp1[7]+z[83]*tmp1[8]+z[57]*tmp1[9]+z[31]*tmp1[10]
            +z[5]*tmp1[11]+z[76]*tmp1[12]+z[50]*tmp1[13]+z[24]*tmp1[14]+z[95]*tmp1[15]+z[69]*tmp1[16]+z[43]*tmp1[17]+z[17]*tmp1[18]+z[88]*tmp1[19]+z[62]*tmp1[20]
            +z[36]*tmp1[21]+z[10]*tmp1[22]+z[81]*tmp1[23]+z[55]*tmp1[24]+z[29]*tmp1[25]+z[3]*tmp1[26]+z[74]*tmp1[27]+z[48]*tmp1[28]+z[22]*tmp1[29]+z[93]*tmp1[30]
            +z[67]*tmp1[31]+z[41]*tmp1[32]+z[15]*tmp1[33]+z[86]*tmp1[34]+z[60]*tmp1[35]+z[34]*tmp1[36]+z[8]*tmp1[37]+z[79]*tmp1[38]+z[53]*tmp1[39]+z[27]*tmp1[40]
            +z[1]*tmp1[41]+z[72]*tmp1[42]+z[46]*tmp1[43]+z[20]*tmp1[44]+z[91]*tmp1[45]+z[65]*tmp1[46]+z[39]*tmp1[47]+z[13]*tmp1[48]+z[84]*tmp1[49]+z[58]*tmp1[50]
            +z[32]*tmp1[51]+z[6]*tmp1[52]+z[77]*tmp1[53]+z[51]*tmp1[54]+z[25]*tmp1[55]+z[96]*tmp1[56]+z[70]*tmp1[57]+z[44]*tmp1[58]+z[18]*tmp1[59]+z[89]*tmp1[60]
            +z[63]*tmp1[61]+z[37]*tmp1[62]+z[11]*tmp1[63]+z[82]*tmp1[64]+z[56]*tmp1[65]+z[30]*tmp1[66]+z[4]*tmp1[67]+z[75]*tmp1[68]+z[49]*tmp1[69]+z[23]*tmp1[70]
            +z[94]*tmp1[71]+z[68]*tmp1[72]+z[42]*tmp1[73]+z[16]*tmp1[74]+z[87]*tmp1[75]+z[61]*tmp1[76]+z[35]*tmp1[77]+z[9]*tmp1[78]+z[80]*tmp1[79]+z[54]*tmp1[80]
            +z[28]*tmp1[81]+z[2]*tmp1[82]+z[73]*tmp1[83]+z[47]*tmp1[84]+z[21]*tmp1[85]+z[92]*tmp1[86]+z[66]*tmp1[87]+z[40]*tmp1[88]+z[14]*tmp1[89]+z[85]*tmp1[90]
            +z[59]*tmp1[91]+z[33]*tmp1[92]+z[7]*tmp1[93]+z[78]*tmp1[94]+z[52]*tmp1[95]+z[26]*tmp1[96];
            tab[i+72*stg_first]=tab[i+72*stg_first]+z[0]*tmp1[0]
            +z[72]*tmp1[1]+z[47]*tmp1[2]+z[22]*tmp1[3]+z[94]*tmp1[4]+z[69]*tmp1[5]+z[44]*tmp1[6]+z[19]*tmp1[7]+z[91]*tmp1[8]+z[66]*tmp1[9]+z[41]*tmp1[10]
            +z[16]*tmp1[11]+z[88]*tmp1[12]+z[63]*tmp1[13]+z[38]*tmp1[14]+z[13]*tmp1[15]+z[85]*tmp1[16]+z[60]*tmp1[17]+z[35]*tmp1[18]+z[10]*tmp1[19]+z[82]*tmp1[20]
            +z[57]*tmp1[21]+z[32]*tmp1[22]+z[7]*tmp1[23]+z[79]*tmp1[24]+z[54]*tmp1[25]+z[29]*tmp1[26]+z[4]*tmp1[27]+z[76]*tmp1[28]+z[51]*tmp1[29]+z[26]*tmp1[30]
            +z[1]*tmp1[31]+z[73]*tmp1[32]+z[48]*tmp1[33]+z[23]*tmp1[34]+z[95]*tmp1[35]+z[70]*tmp1[36]+z[45]*tmp1[37]+z[20]*tmp1[38]+z[92]*tmp1[39]+z[67]*tmp1[40]
            +z[42]*tmp1[41]+z[17]*tmp1[42]+z[89]*tmp1[43]+z[64]*tmp1[44]+z[39]*tmp1[45]+z[14]*tmp1[46]+z[86]*tmp1[47]+z[61]*tmp1[48]+z[36]*tmp1[49]+z[11]*tmp1[50]
            +z[83]*tmp1[51]+z[58]*tmp1[52]+z[33]*tmp1[53]+z[8]*tmp1[54]+z[80]*tmp1[55]+z[55]*tmp1[56]+z[30]*tmp1[57]+z[5]*tmp1[58]+z[77]*tmp1[59]+z[52]*tmp1[60]
            +z[27]*tmp1[61]+z[2]*tmp1[62]+z[74]*tmp1[63]+z[49]*tmp1[64]+z[24]*tmp1[65]+z[96]*tmp1[66]+z[71]*tmp1[67]+z[46]*tmp1[68]+z[21]*tmp1[69]+z[93]*tmp1[70]
            +z[68]*tmp1[71]+z[43]*tmp1[72]+z[18]*tmp1[73]+z[90]*tmp1[74]+z[65]*tmp1[75]+z[40]*tmp1[76]+z[15]*tmp1[77]+z[87]*tmp1[78]+z[62]*tmp1[79]+z[37]*tmp1[80]
            +z[12]*tmp1[81]+z[84]*tmp1[82]+z[59]*tmp1[83]+z[34]*tmp1[84]+z[9]*tmp1[85]+z[81]*tmp1[86]+z[56]*tmp1[87]+z[31]*tmp1[88]+z[6]*tmp1[89]+z[78]*tmp1[90]
            +z[53]*tmp1[91]+z[28]*tmp1[92]+z[3]*tmp1[93]+z[75]*tmp1[94]+z[50]*tmp1[95]+z[25]*tmp1[96];
            tab[i+73*stg_first]=tab[i+73*stg_first]+z[0]*tmp1[0]
            +z[73]*tmp1[1]+z[49]*tmp1[2]+z[25]*tmp1[3]+z[1]*tmp1[4]+z[74]*tmp1[5]+z[50]*tmp1[6]+z[26]*tmp1[7]+z[2]*tmp1[8]+z[75]*tmp1[9]+z[51]*tmp1[10]
            +z[27]*tmp1[11]+z[3]*tmp1[12]+z[76]*tmp1[13]+z[52]*tmp1[14]+z[28]*tmp1[15]+z[4]*tmp1[16]+z[77]*tmp1[17]+z[53]*tmp1[18]+z[29]*tmp1[19]+z[5]*tmp1[20]
            +z[78]*tmp1[21]+z[54]*tmp1[22]+z[30]*tmp1[23]+z[6]*tmp1[24]+z[79]*tmp1[25]+z[55]*tmp1[26]+z[31]*tmp1[27]+z[7]*tmp1[28]+z[80]*tmp1[29]+z[56]*tmp1[30]
            +z[32]*tmp1[31]+z[8]*tmp1[32]+z[81]*tmp1[33]+z[57]*tmp1[34]+z[33]*tmp1[35]+z[9]*tmp1[36]+z[82]*tmp1[37]+z[58]*tmp1[38]+z[34]*tmp1[39]+z[10]*tmp1[40]
            +z[83]*tmp1[41]+z[59]*tmp1[42]+z[35]*tmp1[43]+z[11]*tmp1[44]+z[84]*tmp1[45]+z[60]*tmp1[46]+z[36]*tmp1[47]+z[12]*tmp1[48]+z[85]*tmp1[49]+z[61]*tmp1[50]
            +z[37]*tmp1[51]+z[13]*tmp1[52]+z[86]*tmp1[53]+z[62]*tmp1[54]+z[38]*tmp1[55]+z[14]*tmp1[56]+z[87]*tmp1[57]+z[63]*tmp1[58]+z[39]*tmp1[59]+z[15]*tmp1[60]
            +z[88]*tmp1[61]+z[64]*tmp1[62]+z[40]*tmp1[63]+z[16]*tmp1[64]+z[89]*tmp1[65]+z[65]*tmp1[66]+z[41]*tmp1[67]+z[17]*tmp1[68]+z[90]*tmp1[69]+z[66]*tmp1[70]
            +z[42]*tmp1[71]+z[18]*tmp1[72]+z[91]*tmp1[73]+z[67]*tmp1[74]+z[43]*tmp1[75]+z[19]*tmp1[76]+z[92]*tmp1[77]+z[68]*tmp1[78]+z[44]*tmp1[79]+z[20]*tmp1[80]
            +z[93]*tmp1[81]+z[69]*tmp1[82]+z[45]*tmp1[83]+z[21]*tmp1[84]+z[94]*tmp1[85]+z[70]*tmp1[86]+z[46]*tmp1[87]+z[22]*tmp1[88]+z[95]*tmp1[89]+z[71]*tmp1[90]
            +z[47]*tmp1[91]+z[23]*tmp1[92]+z[96]*tmp1[93]+z[72]*tmp1[94]+z[48]*tmp1[95]+z[24]*tmp1[96];
            tab[i+74*stg_first]=tab[i+74*stg_first]+z[0]*tmp1[0]
            +z[74]*tmp1[1]+z[51]*tmp1[2]+z[28]*tmp1[3]+z[5]*tmp1[4]+z[79]*tmp1[5]+z[56]*tmp1[6]+z[33]*tmp1[7]+z[10]*tmp1[8]+z[84]*tmp1[9]+z[61]*tmp1[10]
            +z[38]*tmp1[11]+z[15]*tmp1[12]+z[89]*tmp1[13]+z[66]*tmp1[14]+z[43]*tmp1[15]+z[20]*tmp1[16]+z[94]*tmp1[17]+z[71]*tmp1[18]+z[48]*tmp1[19]+z[25]*tmp1[20]
            +z[2]*tmp1[21]+z[76]*tmp1[22]+z[53]*tmp1[23]+z[30]*tmp1[24]+z[7]*tmp1[25]+z[81]*tmp1[26]+z[58]*tmp1[27]+z[35]*tmp1[28]+z[12]*tmp1[29]+z[86]*tmp1[30]
            +z[63]*tmp1[31]+z[40]*tmp1[32]+z[17]*tmp1[33]+z[91]*tmp1[34]+z[68]*tmp1[35]+z[45]*tmp1[36]+z[22]*tmp1[37]+z[96]*tmp1[38]+z[73]*tmp1[39]+z[50]*tmp1[40]
            +z[27]*tmp1[41]+z[4]*tmp1[42]+z[78]*tmp1[43]+z[55]*tmp1[44]+z[32]*tmp1[45]+z[9]*tmp1[46]+z[83]*tmp1[47]+z[60]*tmp1[48]+z[37]*tmp1[49]+z[14]*tmp1[50]
            +z[88]*tmp1[51]+z[65]*tmp1[52]+z[42]*tmp1[53]+z[19]*tmp1[54]+z[93]*tmp1[55]+z[70]*tmp1[56]+z[47]*tmp1[57]+z[24]*tmp1[58]+z[1]*tmp1[59]+z[75]*tmp1[60]
            +z[52]*tmp1[61]+z[29]*tmp1[62]+z[6]*tmp1[63]+z[80]*tmp1[64]+z[57]*tmp1[65]+z[34]*tmp1[66]+z[11]*tmp1[67]+z[85]*tmp1[68]+z[62]*tmp1[69]+z[39]*tmp1[70]
            +z[16]*tmp1[71]+z[90]*tmp1[72]+z[67]*tmp1[73]+z[44]*tmp1[74]+z[21]*tmp1[75]+z[95]*tmp1[76]+z[72]*tmp1[77]+z[49]*tmp1[78]+z[26]*tmp1[79]+z[3]*tmp1[80]
            +z[77]*tmp1[81]+z[54]*tmp1[82]+z[31]*tmp1[83]+z[8]*tmp1[84]+z[82]*tmp1[85]+z[59]*tmp1[86]+z[36]*tmp1[87]+z[13]*tmp1[88]+z[87]*tmp1[89]+z[64]*tmp1[90]
            +z[41]*tmp1[91]+z[18]*tmp1[92]+z[92]*tmp1[93]+z[69]*tmp1[94]+z[46]*tmp1[95]+z[23]*tmp1[96];
            tab[i+75*stg_first]=tab[i+75*stg_first]+z[0]*tmp1[0]
            +z[75]*tmp1[1]+z[53]*tmp1[2]+z[31]*tmp1[3]+z[9]*tmp1[4]+z[84]*tmp1[5]+z[62]*tmp1[6]+z[40]*tmp1[7]+z[18]*tmp1[8]+z[93]*tmp1[9]+z[71]*tmp1[10]
            +z[49]*tmp1[11]+z[27]*tmp1[12]+z[5]*tmp1[13]+z[80]*tmp1[14]+z[58]*tmp1[15]+z[36]*tmp1[16]+z[14]*tmp1[17]+z[89]*tmp1[18]+z[67]*tmp1[19]+z[45]*tmp1[20]
            +z[23]*tmp1[21]+z[1]*tmp1[22]+z[76]*tmp1[23]+z[54]*tmp1[24]+z[32]*tmp1[25]+z[10]*tmp1[26]+z[85]*tmp1[27]+z[63]*tmp1[28]+z[41]*tmp1[29]+z[19]*tmp1[30]
            +z[94]*tmp1[31]+z[72]*tmp1[32]+z[50]*tmp1[33]+z[28]*tmp1[34]+z[6]*tmp1[35]+z[81]*tmp1[36]+z[59]*tmp1[37]+z[37]*tmp1[38]+z[15]*tmp1[39]+z[90]*tmp1[40]
            +z[68]*tmp1[41]+z[46]*tmp1[42]+z[24]*tmp1[43]+z[2]*tmp1[44]+z[77]*tmp1[45]+z[55]*tmp1[46]+z[33]*tmp1[47]+z[11]*tmp1[48]+z[86]*tmp1[49]+z[64]*tmp1[50]
            +z[42]*tmp1[51]+z[20]*tmp1[52]+z[95]*tmp1[53]+z[73]*tmp1[54]+z[51]*tmp1[55]+z[29]*tmp1[56]+z[7]*tmp1[57]+z[82]*tmp1[58]+z[60]*tmp1[59]+z[38]*tmp1[60]
            +z[16]*tmp1[61]+z[91]*tmp1[62]+z[69]*tmp1[63]+z[47]*tmp1[64]+z[25]*tmp1[65]+z[3]*tmp1[66]+z[78]*tmp1[67]+z[56]*tmp1[68]+z[34]*tmp1[69]+z[12]*tmp1[70]
            +z[87]*tmp1[71]+z[65]*tmp1[72]+z[43]*tmp1[73]+z[21]*tmp1[74]+z[96]*tmp1[75]+z[74]*tmp1[76]+z[52]*tmp1[77]+z[30]*tmp1[78]+z[8]*tmp1[79]+z[83]*tmp1[80]
            +z[61]*tmp1[81]+z[39]*tmp1[82]+z[17]*tmp1[83]+z[92]*tmp1[84]+z[70]*tmp1[85]+z[48]*tmp1[86]+z[26]*tmp1[87]+z[4]*tmp1[88]+z[79]*tmp1[89]+z[57]*tmp1[90]
            +z[35]*tmp1[91]+z[13]*tmp1[92]+z[88]*tmp1[93]+z[66]*tmp1[94]+z[44]*tmp1[95]+z[22]*tmp1[96];
            tab[i+76*stg_first]=tab[i+76*stg_first]+z[0]*tmp1[0]
            +z[76]*tmp1[1]+z[55]*tmp1[2]+z[34]*tmp1[3]+z[13]*tmp1[4]+z[89]*tmp1[5]+z[68]*tmp1[6]+z[47]*tmp1[7]+z[26]*tmp1[8]+z[5]*tmp1[9]+z[81]*tmp1[10]
            +z[60]*tmp1[11]+z[39]*tmp1[12]+z[18]*tmp1[13]+z[94]*tmp1[14]+z[73]*tmp1[15]+z[52]*tmp1[16]+z[31]*tmp1[17]+z[10]*tmp1[18]+z[86]*tmp1[19]+z[65]*tmp1[20]
            +z[44]*tmp1[21]+z[23]*tmp1[22]+z[2]*tmp1[23]+z[78]*tmp1[24]+z[57]*tmp1[25]+z[36]*tmp1[26]+z[15]*tmp1[27]+z[91]*tmp1[28]+z[70]*tmp1[29]+z[49]*tmp1[30]
            +z[28]*tmp1[31]+z[7]*tmp1[32]+z[83]*tmp1[33]+z[62]*tmp1[34]+z[41]*tmp1[35]+z[20]*tmp1[36]+z[96]*tmp1[37]+z[75]*tmp1[38]+z[54]*tmp1[39]+z[33]*tmp1[40]
            +z[12]*tmp1[41]+z[88]*tmp1[42]+z[67]*tmp1[43]+z[46]*tmp1[44]+z[25]*tmp1[45]+z[4]*tmp1[46]+z[80]*tmp1[47]+z[59]*tmp1[48]+z[38]*tmp1[49]+z[17]*tmp1[50]
            +z[93]*tmp1[51]+z[72]*tmp1[52]+z[51]*tmp1[53]+z[30]*tmp1[54]+z[9]*tmp1[55]+z[85]*tmp1[56]+z[64]*tmp1[57]+z[43]*tmp1[58]+z[22]*tmp1[59]+z[1]*tmp1[60]
            +z[77]*tmp1[61]+z[56]*tmp1[62]+z[35]*tmp1[63]+z[14]*tmp1[64]+z[90]*tmp1[65]+z[69]*tmp1[66]+z[48]*tmp1[67]+z[27]*tmp1[68]+z[6]*tmp1[69]+z[82]*tmp1[70]
            +z[61]*tmp1[71]+z[40]*tmp1[72]+z[19]*tmp1[73]+z[95]*tmp1[74]+z[74]*tmp1[75]+z[53]*tmp1[76]+z[32]*tmp1[77]+z[11]*tmp1[78]+z[87]*tmp1[79]+z[66]*tmp1[80]
            +z[45]*tmp1[81]+z[24]*tmp1[82]+z[3]*tmp1[83]+z[79]*tmp1[84]+z[58]*tmp1[85]+z[37]*tmp1[86]+z[16]*tmp1[87]+z[92]*tmp1[88]+z[71]*tmp1[89]+z[50]*tmp1[90]
            +z[29]*tmp1[91]+z[8]*tmp1[92]+z[84]*tmp1[93]+z[63]*tmp1[94]+z[42]*tmp1[95]+z[21]*tmp1[96];
            tab[i+77*stg_first]=tab[i+77*stg_first]+z[0]*tmp1[0]
            +z[77]*tmp1[1]+z[57]*tmp1[2]+z[37]*tmp1[3]+z[17]*tmp1[4]+z[94]*tmp1[5]+z[74]*tmp1[6]+z[54]*tmp1[7]+z[34]*tmp1[8]+z[14]*tmp1[9]+z[91]*tmp1[10]
            +z[71]*tmp1[11]+z[51]*tmp1[12]+z[31]*tmp1[13]+z[11]*tmp1[14]+z[88]*tmp1[15]+z[68]*tmp1[16]+z[48]*tmp1[17]+z[28]*tmp1[18]+z[8]*tmp1[19]+z[85]*tmp1[20]
            +z[65]*tmp1[21]+z[45]*tmp1[22]+z[25]*tmp1[23]+z[5]*tmp1[24]+z[82]*tmp1[25]+z[62]*tmp1[26]+z[42]*tmp1[27]+z[22]*tmp1[28]+z[2]*tmp1[29]+z[79]*tmp1[30]
            +z[59]*tmp1[31]+z[39]*tmp1[32]+z[19]*tmp1[33]+z[96]*tmp1[34]+z[76]*tmp1[35]+z[56]*tmp1[36]+z[36]*tmp1[37]+z[16]*tmp1[38]+z[93]*tmp1[39]+z[73]*tmp1[40]
            +z[53]*tmp1[41]+z[33]*tmp1[42]+z[13]*tmp1[43]+z[90]*tmp1[44]+z[70]*tmp1[45]+z[50]*tmp1[46]+z[30]*tmp1[47]+z[10]*tmp1[48]+z[87]*tmp1[49]+z[67]*tmp1[50]
            +z[47]*tmp1[51]+z[27]*tmp1[52]+z[7]*tmp1[53]+z[84]*tmp1[54]+z[64]*tmp1[55]+z[44]*tmp1[56]+z[24]*tmp1[57]+z[4]*tmp1[58]+z[81]*tmp1[59]+z[61]*tmp1[60]
            +z[41]*tmp1[61]+z[21]*tmp1[62]+z[1]*tmp1[63]+z[78]*tmp1[64]+z[58]*tmp1[65]+z[38]*tmp1[66]+z[18]*tmp1[67]+z[95]*tmp1[68]+z[75]*tmp1[69]+z[55]*tmp1[70]
            +z[35]*tmp1[71]+z[15]*tmp1[72]+z[92]*tmp1[73]+z[72]*tmp1[74]+z[52]*tmp1[75]+z[32]*tmp1[76]+z[12]*tmp1[77]+z[89]*tmp1[78]+z[69]*tmp1[79]+z[49]*tmp1[80]
            +z[29]*tmp1[81]+z[9]*tmp1[82]+z[86]*tmp1[83]+z[66]*tmp1[84]+z[46]*tmp1[85]+z[26]*tmp1[86]+z[6]*tmp1[87]+z[83]*tmp1[88]+z[63]*tmp1[89]+z[43]*tmp1[90]
            +z[23]*tmp1[91]+z[3]*tmp1[92]+z[80]*tmp1[93]+z[60]*tmp1[94]+z[40]*tmp1[95]+z[20]*tmp1[96];
            tab[i+78*stg_first]=tab[i+78*stg_first]+z[0]*tmp1[0]
            +z[78]*tmp1[1]+z[59]*tmp1[2]+z[40]*tmp1[3]+z[21]*tmp1[4]+z[2]*tmp1[5]+z[80]*tmp1[6]+z[61]*tmp1[7]+z[42]*tmp1[8]+z[23]*tmp1[9]+z[4]*tmp1[10]
            +z[82]*tmp1[11]+z[63]*tmp1[12]+z[44]*tmp1[13]+z[25]*tmp1[14]+z[6]*tmp1[15]+z[84]*tmp1[16]+z[65]*tmp1[17]+z[46]*tmp1[18]+z[27]*tmp1[19]+z[8]*tmp1[20]
            +z[86]*tmp1[21]+z[67]*tmp1[22]+z[48]*tmp1[23]+z[29]*tmp1[24]+z[10]*tmp1[25]+z[88]*tmp1[26]+z[69]*tmp1[27]+z[50]*tmp1[28]+z[31]*tmp1[29]+z[12]*tmp1[30]
            +z[90]*tmp1[31]+z[71]*tmp1[32]+z[52]*tmp1[33]+z[33]*tmp1[34]+z[14]*tmp1[35]+z[92]*tmp1[36]+z[73]*tmp1[37]+z[54]*tmp1[38]+z[35]*tmp1[39]+z[16]*tmp1[40]
            +z[94]*tmp1[41]+z[75]*tmp1[42]+z[56]*tmp1[43]+z[37]*tmp1[44]+z[18]*tmp1[45]+z[96]*tmp1[46]+z[77]*tmp1[47]+z[58]*tmp1[48]+z[39]*tmp1[49]+z[20]*tmp1[50]
            +z[1]*tmp1[51]+z[79]*tmp1[52]+z[60]*tmp1[53]+z[41]*tmp1[54]+z[22]*tmp1[55]+z[3]*tmp1[56]+z[81]*tmp1[57]+z[62]*tmp1[58]+z[43]*tmp1[59]+z[24]*tmp1[60]
            +z[5]*tmp1[61]+z[83]*tmp1[62]+z[64]*tmp1[63]+z[45]*tmp1[64]+z[26]*tmp1[65]+z[7]*tmp1[66]+z[85]*tmp1[67]+z[66]*tmp1[68]+z[47]*tmp1[69]+z[28]*tmp1[70]
            +z[9]*tmp1[71]+z[87]*tmp1[72]+z[68]*tmp1[73]+z[49]*tmp1[74]+z[30]*tmp1[75]+z[11]*tmp1[76]+z[89]*tmp1[77]+z[70]*tmp1[78]+z[51]*tmp1[79]+z[32]*tmp1[80]
            +z[13]*tmp1[81]+z[91]*tmp1[82]+z[72]*tmp1[83]+z[53]*tmp1[84]+z[34]*tmp1[85]+z[15]*tmp1[86]+z[93]*tmp1[87]+z[74]*tmp1[88]+z[55]*tmp1[89]+z[36]*tmp1[90]
            +z[17]*tmp1[91]+z[95]*tmp1[92]+z[76]*tmp1[93]+z[57]*tmp1[94]+z[38]*tmp1[95]+z[19]*tmp1[96];
            tab[i+79*stg_first]=tab[i+79*stg_first]+z[0]*tmp1[0]
            +z[79]*tmp1[1]+z[61]*tmp1[2]+z[43]*tmp1[3]+z[25]*tmp1[4]+z[7]*tmp1[5]+z[86]*tmp1[6]+z[68]*tmp1[7]+z[50]*tmp1[8]+z[32]*tmp1[9]+z[14]*tmp1[10]
            +z[93]*tmp1[11]+z[75]*tmp1[12]+z[57]*tmp1[13]+z[39]*tmp1[14]+z[21]*tmp1[15]+z[3]*tmp1[16]+z[82]*tmp1[17]+z[64]*tmp1[18]+z[46]*tmp1[19]+z[28]*tmp1[20]
            +z[10]*tmp1[21]+z[89]*tmp1[22]+z[71]*tmp1[23]+z[53]*tmp1[24]+z[35]*tmp1[25]+z[17]*tmp1[26]+z[96]*tmp1[27]+z[78]*tmp1[28]+z[60]*tmp1[29]+z[42]*tmp1[30]
            +z[24]*tmp1[31]+z[6]*tmp1[32]+z[85]*tmp1[33]+z[67]*tmp1[34]+z[49]*tmp1[35]+z[31]*tmp1[36]+z[13]*tmp1[37]+z[92]*tmp1[38]+z[74]*tmp1[39]+z[56]*tmp1[40]
            +z[38]*tmp1[41]+z[20]*tmp1[42]+z[2]*tmp1[43]+z[81]*tmp1[44]+z[63]*tmp1[45]+z[45]*tmp1[46]+z[27]*tmp1[47]+z[9]*tmp1[48]+z[88]*tmp1[49]+z[70]*tmp1[50]
            +z[52]*tmp1[51]+z[34]*tmp1[52]+z[16]*tmp1[53]+z[95]*tmp1[54]+z[77]*tmp1[55]+z[59]*tmp1[56]+z[41]*tmp1[57]+z[23]*tmp1[58]+z[5]*tmp1[59]+z[84]*tmp1[60]
            +z[66]*tmp1[61]+z[48]*tmp1[62]+z[30]*tmp1[63]+z[12]*tmp1[64]+z[91]*tmp1[65]+z[73]*tmp1[66]+z[55]*tmp1[67]+z[37]*tmp1[68]+z[19]*tmp1[69]+z[1]*tmp1[70]
            +z[80]*tmp1[71]+z[62]*tmp1[72]+z[44]*tmp1[73]+z[26]*tmp1[74]+z[8]*tmp1[75]+z[87]*tmp1[76]+z[69]*tmp1[77]+z[51]*tmp1[78]+z[33]*tmp1[79]+z[15]*tmp1[80]
            +z[94]*tmp1[81]+z[76]*tmp1[82]+z[58]*tmp1[83]+z[40]*tmp1[84]+z[22]*tmp1[85]+z[4]*tmp1[86]+z[83]*tmp1[87]+z[65]*tmp1[88]+z[47]*tmp1[89]+z[29]*tmp1[90]
            +z[11]*tmp1[91]+z[90]*tmp1[92]+z[72]*tmp1[93]+z[54]*tmp1[94]+z[36]*tmp1[95]+z[18]*tmp1[96];
            tab[i+80*stg_first]=tab[i+80*stg_first]+z[0]*tmp1[0]
            +z[80]*tmp1[1]+z[63]*tmp1[2]+z[46]*tmp1[3]+z[29]*tmp1[4]+z[12]*tmp1[5]+z[92]*tmp1[6]+z[75]*tmp1[7]+z[58]*tmp1[8]+z[41]*tmp1[9]+z[24]*tmp1[10]
            +z[7]*tmp1[11]+z[87]*tmp1[12]+z[70]*tmp1[13]+z[53]*tmp1[14]+z[36]*tmp1[15]+z[19]*tmp1[16]+z[2]*tmp1[17]+z[82]*tmp1[18]+z[65]*tmp1[19]+z[48]*tmp1[20]
            +z[31]*tmp1[21]+z[14]*tmp1[22]+z[94]*tmp1[23]+z[77]*tmp1[24]+z[60]*tmp1[25]+z[43]*tmp1[26]+z[26]*tmp1[27]+z[9]*tmp1[28]+z[89]*tmp1[29]+z[72]*tmp1[30]
            +z[55]*tmp1[31]+z[38]*tmp1[32]+z[21]*tmp1[33]+z[4]*tmp1[34]+z[84]*tmp1[35]+z[67]*tmp1[36]+z[50]*tmp1[37]+z[33]*tmp1[38]+z[16]*tmp1[39]+z[96]*tmp1[40]
            +z[79]*tmp1[41]+z[62]*tmp1[42]+z[45]*tmp1[43]+z[28]*tmp1[44]+z[11]*tmp1[45]+z[91]*tmp1[46]+z[74]*tmp1[47]+z[57]*tmp1[48]+z[40]*tmp1[49]+z[23]*tmp1[50]
            +z[6]*tmp1[51]+z[86]*tmp1[52]+z[69]*tmp1[53]+z[52]*tmp1[54]+z[35]*tmp1[55]+z[18]*tmp1[56]+z[1]*tmp1[57]+z[81]*tmp1[58]+z[64]*tmp1[59]+z[47]*tmp1[60]
            +z[30]*tmp1[61]+z[13]*tmp1[62]+z[93]*tmp1[63]+z[76]*tmp1[64]+z[59]*tmp1[65]+z[42]*tmp1[66]+z[25]*tmp1[67]+z[8]*tmp1[68]+z[88]*tmp1[69]+z[71]*tmp1[70]
            +z[54]*tmp1[71]+z[37]*tmp1[72]+z[20]*tmp1[73]+z[3]*tmp1[74]+z[83]*tmp1[75]+z[66]*tmp1[76]+z[49]*tmp1[77]+z[32]*tmp1[78]+z[15]*tmp1[79]+z[95]*tmp1[80]
            +z[78]*tmp1[81]+z[61]*tmp1[82]+z[44]*tmp1[83]+z[27]*tmp1[84]+z[10]*tmp1[85]+z[90]*tmp1[86]+z[73]*tmp1[87]+z[56]*tmp1[88]+z[39]*tmp1[89]+z[22]*tmp1[90]
            +z[5]*tmp1[91]+z[85]*tmp1[92]+z[68]*tmp1[93]+z[51]*tmp1[94]+z[34]*tmp1[95]+z[17]*tmp1[96];
            tab[i+81*stg_first]=tab[i+81*stg_first]+z[0]*tmp1[0]
            +z[81]*tmp1[1]+z[65]*tmp1[2]+z[49]*tmp1[3]+z[33]*tmp1[4]+z[17]*tmp1[5]+z[1]*tmp1[6]+z[82]*tmp1[7]+z[66]*tmp1[8]+z[50]*tmp1[9]+z[34]*tmp1[10]
            +z[18]*tmp1[11]+z[2]*tmp1[12]+z[83]*tmp1[13]+z[67]*tmp1[14]+z[51]*tmp1[15]+z[35]*tmp1[16]+z[19]*tmp1[17]+z[3]*tmp1[18]+z[84]*tmp1[19]+z[68]*tmp1[20]
            +z[52]*tmp1[21]+z[36]*tmp1[22]+z[20]*tmp1[23]+z[4]*tmp1[24]+z[85]*tmp1[25]+z[69]*tmp1[26]+z[53]*tmp1[27]+z[37]*tmp1[28]+z[21]*tmp1[29]+z[5]*tmp1[30]
            +z[86]*tmp1[31]+z[70]*tmp1[32]+z[54]*tmp1[33]+z[38]*tmp1[34]+z[22]*tmp1[35]+z[6]*tmp1[36]+z[87]*tmp1[37]+z[71]*tmp1[38]+z[55]*tmp1[39]+z[39]*tmp1[40]
            +z[23]*tmp1[41]+z[7]*tmp1[42]+z[88]*tmp1[43]+z[72]*tmp1[44]+z[56]*tmp1[45]+z[40]*tmp1[46]+z[24]*tmp1[47]+z[8]*tmp1[48]+z[89]*tmp1[49]+z[73]*tmp1[50]
            +z[57]*tmp1[51]+z[41]*tmp1[52]+z[25]*tmp1[53]+z[9]*tmp1[54]+z[90]*tmp1[55]+z[74]*tmp1[56]+z[58]*tmp1[57]+z[42]*tmp1[58]+z[26]*tmp1[59]+z[10]*tmp1[60]
            +z[91]*tmp1[61]+z[75]*tmp1[62]+z[59]*tmp1[63]+z[43]*tmp1[64]+z[27]*tmp1[65]+z[11]*tmp1[66]+z[92]*tmp1[67]+z[76]*tmp1[68]+z[60]*tmp1[69]+z[44]*tmp1[70]
            +z[28]*tmp1[71]+z[12]*tmp1[72]+z[93]*tmp1[73]+z[77]*tmp1[74]+z[61]*tmp1[75]+z[45]*tmp1[76]+z[29]*tmp1[77]+z[13]*tmp1[78]+z[94]*tmp1[79]+z[78]*tmp1[80]
            +z[62]*tmp1[81]+z[46]*tmp1[82]+z[30]*tmp1[83]+z[14]*tmp1[84]+z[95]*tmp1[85]+z[79]*tmp1[86]+z[63]*tmp1[87]+z[47]*tmp1[88]+z[31]*tmp1[89]+z[15]*tmp1[90]
            +z[96]*tmp1[91]+z[80]*tmp1[92]+z[64]*tmp1[93]+z[48]*tmp1[94]+z[32]*tmp1[95]+z[16]*tmp1[96];
            tab[i+82*stg_first]=tab[i+82*stg_first]+z[0]*tmp1[0]
            +z[82]*tmp1[1]+z[67]*tmp1[2]+z[52]*tmp1[3]+z[37]*tmp1[4]+z[22]*tmp1[5]+z[7]*tmp1[6]+z[89]*tmp1[7]+z[74]*tmp1[8]+z[59]*tmp1[9]+z[44]*tmp1[10]
            +z[29]*tmp1[11]+z[14]*tmp1[12]+z[96]*tmp1[13]+z[81]*tmp1[14]+z[66]*tmp1[15]+z[51]*tmp1[16]+z[36]*tmp1[17]+z[21]*tmp1[18]+z[6]*tmp1[19]+z[88]*tmp1[20]
            +z[73]*tmp1[21]+z[58]*tmp1[22]+z[43]*tmp1[23]+z[28]*tmp1[24]+z[13]*tmp1[25]+z[95]*tmp1[26]+z[80]*tmp1[27]+z[65]*tmp1[28]+z[50]*tmp1[29]+z[35]*tmp1[30]
            +z[20]*tmp1[31]+z[5]*tmp1[32]+z[87]*tmp1[33]+z[72]*tmp1[34]+z[57]*tmp1[35]+z[42]*tmp1[36]+z[27]*tmp1[37]+z[12]*tmp1[38]+z[94]*tmp1[39]+z[79]*tmp1[40]
            +z[64]*tmp1[41]+z[49]*tmp1[42]+z[34]*tmp1[43]+z[19]*tmp1[44]+z[4]*tmp1[45]+z[86]*tmp1[46]+z[71]*tmp1[47]+z[56]*tmp1[48]+z[41]*tmp1[49]+z[26]*tmp1[50]
            +z[11]*tmp1[51]+z[93]*tmp1[52]+z[78]*tmp1[53]+z[63]*tmp1[54]+z[48]*tmp1[55]+z[33]*tmp1[56]+z[18]*tmp1[57]+z[3]*tmp1[58]+z[85]*tmp1[59]+z[70]*tmp1[60]
            +z[55]*tmp1[61]+z[40]*tmp1[62]+z[25]*tmp1[63]+z[10]*tmp1[64]+z[92]*tmp1[65]+z[77]*tmp1[66]+z[62]*tmp1[67]+z[47]*tmp1[68]+z[32]*tmp1[69]+z[17]*tmp1[70]
            +z[2]*tmp1[71]+z[84]*tmp1[72]+z[69]*tmp1[73]+z[54]*tmp1[74]+z[39]*tmp1[75]+z[24]*tmp1[76]+z[9]*tmp1[77]+z[91]*tmp1[78]+z[76]*tmp1[79]+z[61]*tmp1[80]
            +z[46]*tmp1[81]+z[31]*tmp1[82]+z[16]*tmp1[83]+z[1]*tmp1[84]+z[83]*tmp1[85]+z[68]*tmp1[86]+z[53]*tmp1[87]+z[38]*tmp1[88]+z[23]*tmp1[89]+z[8]*tmp1[90]
            +z[90]*tmp1[91]+z[75]*tmp1[92]+z[60]*tmp1[93]+z[45]*tmp1[94]+z[30]*tmp1[95]+z[15]*tmp1[96];
            tab[i+83*stg_first]=tab[i+83*stg_first]+z[0]*tmp1[0]
            +z[83]*tmp1[1]+z[69]*tmp1[2]+z[55]*tmp1[3]+z[41]*tmp1[4]+z[27]*tmp1[5]+z[13]*tmp1[6]+z[96]*tmp1[7]+z[82]*tmp1[8]+z[68]*tmp1[9]+z[54]*tmp1[10]
            +z[40]*tmp1[11]+z[26]*tmp1[12]+z[12]*tmp1[13]+z[95]*tmp1[14]+z[81]*tmp1[15]+z[67]*tmp1[16]+z[53]*tmp1[17]+z[39]*tmp1[18]+z[25]*tmp1[19]+z[11]*tmp1[20]
            +z[94]*tmp1[21]+z[80]*tmp1[22]+z[66]*tmp1[23]+z[52]*tmp1[24]+z[38]*tmp1[25]+z[24]*tmp1[26]+z[10]*tmp1[27]+z[93]*tmp1[28]+z[79]*tmp1[29]+z[65]*tmp1[30]
            +z[51]*tmp1[31]+z[37]*tmp1[32]+z[23]*tmp1[33]+z[9]*tmp1[34]+z[92]*tmp1[35]+z[78]*tmp1[36]+z[64]*tmp1[37]+z[50]*tmp1[38]+z[36]*tmp1[39]+z[22]*tmp1[40]
            +z[8]*tmp1[41]+z[91]*tmp1[42]+z[77]*tmp1[43]+z[63]*tmp1[44]+z[49]*tmp1[45]+z[35]*tmp1[46]+z[21]*tmp1[47]+z[7]*tmp1[48]+z[90]*tmp1[49]+z[76]*tmp1[50]
            +z[62]*tmp1[51]+z[48]*tmp1[52]+z[34]*tmp1[53]+z[20]*tmp1[54]+z[6]*tmp1[55]+z[89]*tmp1[56]+z[75]*tmp1[57]+z[61]*tmp1[58]+z[47]*tmp1[59]+z[33]*tmp1[60]
            +z[19]*tmp1[61]+z[5]*tmp1[62]+z[88]*tmp1[63]+z[74]*tmp1[64]+z[60]*tmp1[65]+z[46]*tmp1[66]+z[32]*tmp1[67]+z[18]*tmp1[68]+z[4]*tmp1[69]+z[87]*tmp1[70]
            +z[73]*tmp1[71]+z[59]*tmp1[72]+z[45]*tmp1[73]+z[31]*tmp1[74]+z[17]*tmp1[75]+z[3]*tmp1[76]+z[86]*tmp1[77]+z[72]*tmp1[78]+z[58]*tmp1[79]+z[44]*tmp1[80]
            +z[30]*tmp1[81]+z[16]*tmp1[82]+z[2]*tmp1[83]+z[85]*tmp1[84]+z[71]*tmp1[85]+z[57]*tmp1[86]+z[43]*tmp1[87]+z[29]*tmp1[88]+z[15]*tmp1[89]+z[1]*tmp1[90]
            +z[84]*tmp1[91]+z[70]*tmp1[92]+z[56]*tmp1[93]+z[42]*tmp1[94]+z[28]*tmp1[95]+z[14]*tmp1[96];
            tab[i+84*stg_first]=tab[i+84*stg_first]+z[0]*tmp1[0]
            +z[84]*tmp1[1]+z[71]*tmp1[2]+z[58]*tmp1[3]+z[45]*tmp1[4]+z[32]*tmp1[5]+z[19]*tmp1[6]+z[6]*tmp1[7]+z[90]*tmp1[8]+z[77]*tmp1[9]+z[64]*tmp1[10]
            +z[51]*tmp1[11]+z[38]*tmp1[12]+z[25]*tmp1[13]+z[12]*tmp1[14]+z[96]*tmp1[15]+z[83]*tmp1[16]+z[70]*tmp1[17]+z[57]*tmp1[18]+z[44]*tmp1[19]+z[31]*tmp1[20]
            +z[18]*tmp1[21]+z[5]*tmp1[22]+z[89]*tmp1[23]+z[76]*tmp1[24]+z[63]*tmp1[25]+z[50]*tmp1[26]+z[37]*tmp1[27]+z[24]*tmp1[28]+z[11]*tmp1[29]+z[95]*tmp1[30]
            +z[82]*tmp1[31]+z[69]*tmp1[32]+z[56]*tmp1[33]+z[43]*tmp1[34]+z[30]*tmp1[35]+z[17]*tmp1[36]+z[4]*tmp1[37]+z[88]*tmp1[38]+z[75]*tmp1[39]+z[62]*tmp1[40]
            +z[49]*tmp1[41]+z[36]*tmp1[42]+z[23]*tmp1[43]+z[10]*tmp1[44]+z[94]*tmp1[45]+z[81]*tmp1[46]+z[68]*tmp1[47]+z[55]*tmp1[48]+z[42]*tmp1[49]+z[29]*tmp1[50]
            +z[16]*tmp1[51]+z[3]*tmp1[52]+z[87]*tmp1[53]+z[74]*tmp1[54]+z[61]*tmp1[55]+z[48]*tmp1[56]+z[35]*tmp1[57]+z[22]*tmp1[58]+z[9]*tmp1[59]+z[93]*tmp1[60]
            +z[80]*tmp1[61]+z[67]*tmp1[62]+z[54]*tmp1[63]+z[41]*tmp1[64]+z[28]*tmp1[65]+z[15]*tmp1[66]+z[2]*tmp1[67]+z[86]*tmp1[68]+z[73]*tmp1[69]+z[60]*tmp1[70]
            +z[47]*tmp1[71]+z[34]*tmp1[72]+z[21]*tmp1[73]+z[8]*tmp1[74]+z[92]*tmp1[75]+z[79]*tmp1[76]+z[66]*tmp1[77]+z[53]*tmp1[78]+z[40]*tmp1[79]+z[27]*tmp1[80]
            +z[14]*tmp1[81]+z[1]*tmp1[82]+z[85]*tmp1[83]+z[72]*tmp1[84]+z[59]*tmp1[85]+z[46]*tmp1[86]+z[33]*tmp1[87]+z[20]*tmp1[88]+z[7]*tmp1[89]+z[91]*tmp1[90]
            +z[78]*tmp1[91]+z[65]*tmp1[92]+z[52]*tmp1[93]+z[39]*tmp1[94]+z[26]*tmp1[95]+z[13]*tmp1[96];
            tab[i+85*stg_first]=tab[i+85*stg_first]+z[0]*tmp1[0]
            +z[85]*tmp1[1]+z[73]*tmp1[2]+z[61]*tmp1[3]+z[49]*tmp1[4]+z[37]*tmp1[5]+z[25]*tmp1[6]+z[13]*tmp1[7]+z[1]*tmp1[8]+z[86]*tmp1[9]+z[74]*tmp1[10]
            +z[62]*tmp1[11]+z[50]*tmp1[12]+z[38]*tmp1[13]+z[26]*tmp1[14]+z[14]*tmp1[15]+z[2]*tmp1[16]+z[87]*tmp1[17]+z[75]*tmp1[18]+z[63]*tmp1[19]+z[51]*tmp1[20]
            +z[39]*tmp1[21]+z[27]*tmp1[22]+z[15]*tmp1[23]+z[3]*tmp1[24]+z[88]*tmp1[25]+z[76]*tmp1[26]+z[64]*tmp1[27]+z[52]*tmp1[28]+z[40]*tmp1[29]+z[28]*tmp1[30]
            +z[16]*tmp1[31]+z[4]*tmp1[32]+z[89]*tmp1[33]+z[77]*tmp1[34]+z[65]*tmp1[35]+z[53]*tmp1[36]+z[41]*tmp1[37]+z[29]*tmp1[38]+z[17]*tmp1[39]+z[5]*tmp1[40]
            +z[90]*tmp1[41]+z[78]*tmp1[42]+z[66]*tmp1[43]+z[54]*tmp1[44]+z[42]*tmp1[45]+z[30]*tmp1[46]+z[18]*tmp1[47]+z[6]*tmp1[48]+z[91]*tmp1[49]+z[79]*tmp1[50]
            +z[67]*tmp1[51]+z[55]*tmp1[52]+z[43]*tmp1[53]+z[31]*tmp1[54]+z[19]*tmp1[55]+z[7]*tmp1[56]+z[92]*tmp1[57]+z[80]*tmp1[58]+z[68]*tmp1[59]+z[56]*tmp1[60]
            +z[44]*tmp1[61]+z[32]*tmp1[62]+z[20]*tmp1[63]+z[8]*tmp1[64]+z[93]*tmp1[65]+z[81]*tmp1[66]+z[69]*tmp1[67]+z[57]*tmp1[68]+z[45]*tmp1[69]+z[33]*tmp1[70]
            +z[21]*tmp1[71]+z[9]*tmp1[72]+z[94]*tmp1[73]+z[82]*tmp1[74]+z[70]*tmp1[75]+z[58]*tmp1[76]+z[46]*tmp1[77]+z[34]*tmp1[78]+z[22]*tmp1[79]+z[10]*tmp1[80]
            +z[95]*tmp1[81]+z[83]*tmp1[82]+z[71]*tmp1[83]+z[59]*tmp1[84]+z[47]*tmp1[85]+z[35]*tmp1[86]+z[23]*tmp1[87]+z[11]*tmp1[88]+z[96]*tmp1[89]+z[84]*tmp1[90]
            +z[72]*tmp1[91]+z[60]*tmp1[92]+z[48]*tmp1[93]+z[36]*tmp1[94]+z[24]*tmp1[95]+z[12]*tmp1[96];
            tab[i+86*stg_first]=tab[i+86*stg_first]+z[0]*tmp1[0]
            +z[86]*tmp1[1]+z[75]*tmp1[2]+z[64]*tmp1[3]+z[53]*tmp1[4]+z[42]*tmp1[5]+z[31]*tmp1[6]+z[20]*tmp1[7]+z[9]*tmp1[8]+z[95]*tmp1[9]+z[84]*tmp1[10]
            +z[73]*tmp1[11]+z[62]*tmp1[12]+z[51]*tmp1[13]+z[40]*tmp1[14]+z[29]*tmp1[15]+z[18]*tmp1[16]+z[7]*tmp1[17]+z[93]*tmp1[18]+z[82]*tmp1[19]+z[71]*tmp1[20]
            +z[60]*tmp1[21]+z[49]*tmp1[22]+z[38]*tmp1[23]+z[27]*tmp1[24]+z[16]*tmp1[25]+z[5]*tmp1[26]+z[91]*tmp1[27]+z[80]*tmp1[28]+z[69]*tmp1[29]+z[58]*tmp1[30]
            +z[47]*tmp1[31]+z[36]*tmp1[32]+z[25]*tmp1[33]+z[14]*tmp1[34]+z[3]*tmp1[35]+z[89]*tmp1[36]+z[78]*tmp1[37]+z[67]*tmp1[38]+z[56]*tmp1[39]+z[45]*tmp1[40]
            +z[34]*tmp1[41]+z[23]*tmp1[42]+z[12]*tmp1[43]+z[1]*tmp1[44]+z[87]*tmp1[45]+z[76]*tmp1[46]+z[65]*tmp1[47]+z[54]*tmp1[48]+z[43]*tmp1[49]+z[32]*tmp1[50]
            +z[21]*tmp1[51]+z[10]*tmp1[52]+z[96]*tmp1[53]+z[85]*tmp1[54]+z[74]*tmp1[55]+z[63]*tmp1[56]+z[52]*tmp1[57]+z[41]*tmp1[58]+z[30]*tmp1[59]+z[19]*tmp1[60]
            +z[8]*tmp1[61]+z[94]*tmp1[62]+z[83]*tmp1[63]+z[72]*tmp1[64]+z[61]*tmp1[65]+z[50]*tmp1[66]+z[39]*tmp1[67]+z[28]*tmp1[68]+z[17]*tmp1[69]+z[6]*tmp1[70]
            +z[92]*tmp1[71]+z[81]*tmp1[72]+z[70]*tmp1[73]+z[59]*tmp1[74]+z[48]*tmp1[75]+z[37]*tmp1[76]+z[26]*tmp1[77]+z[15]*tmp1[78]+z[4]*tmp1[79]+z[90]*tmp1[80]
            +z[79]*tmp1[81]+z[68]*tmp1[82]+z[57]*tmp1[83]+z[46]*tmp1[84]+z[35]*tmp1[85]+z[24]*tmp1[86]+z[13]*tmp1[87]+z[2]*tmp1[88]+z[88]*tmp1[89]+z[77]*tmp1[90]
            +z[66]*tmp1[91]+z[55]*tmp1[92]+z[44]*tmp1[93]+z[33]*tmp1[94]+z[22]*tmp1[95]+z[11]*tmp1[96];
            tab[i+87*stg_first]=tab[i+87*stg_first]+z[0]*tmp1[0]
            +z[87]*tmp1[1]+z[77]*tmp1[2]+z[67]*tmp1[3]+z[57]*tmp1[4]+z[47]*tmp1[5]+z[37]*tmp1[6]+z[27]*tmp1[7]+z[17]*tmp1[8]+z[7]*tmp1[9]+z[94]*tmp1[10]
            +z[84]*tmp1[11]+z[74]*tmp1[12]+z[64]*tmp1[13]+z[54]*tmp1[14]+z[44]*tmp1[15]+z[34]*tmp1[16]+z[24]*tmp1[17]+z[14]*tmp1[18]+z[4]*tmp1[19]+z[91]*tmp1[20]
            +z[81]*tmp1[21]+z[71]*tmp1[22]+z[61]*tmp1[23]+z[51]*tmp1[24]+z[41]*tmp1[25]+z[31]*tmp1[26]+z[21]*tmp1[27]+z[11]*tmp1[28]+z[1]*tmp1[29]+z[88]*tmp1[30]
            +z[78]*tmp1[31]+z[68]*tmp1[32]+z[58]*tmp1[33]+z[48]*tmp1[34]+z[38]*tmp1[35]+z[28]*tmp1[36]+z[18]*tmp1[37]+z[8]*tmp1[38]+z[95]*tmp1[39]+z[85]*tmp1[40]
            +z[75]*tmp1[41]+z[65]*tmp1[42]+z[55]*tmp1[43]+z[45]*tmp1[44]+z[35]*tmp1[45]+z[25]*tmp1[46]+z[15]*tmp1[47]+z[5]*tmp1[48]+z[92]*tmp1[49]+z[82]*tmp1[50]
            +z[72]*tmp1[51]+z[62]*tmp1[52]+z[52]*tmp1[53]+z[42]*tmp1[54]+z[32]*tmp1[55]+z[22]*tmp1[56]+z[12]*tmp1[57]+z[2]*tmp1[58]+z[89]*tmp1[59]+z[79]*tmp1[60]
            +z[69]*tmp1[61]+z[59]*tmp1[62]+z[49]*tmp1[63]+z[39]*tmp1[64]+z[29]*tmp1[65]+z[19]*tmp1[66]+z[9]*tmp1[67]+z[96]*tmp1[68]+z[86]*tmp1[69]+z[76]*tmp1[70]
            +z[66]*tmp1[71]+z[56]*tmp1[72]+z[46]*tmp1[73]+z[36]*tmp1[74]+z[26]*tmp1[75]+z[16]*tmp1[76]+z[6]*tmp1[77]+z[93]*tmp1[78]+z[83]*tmp1[79]+z[73]*tmp1[80]
            +z[63]*tmp1[81]+z[53]*tmp1[82]+z[43]*tmp1[83]+z[33]*tmp1[84]+z[23]*tmp1[85]+z[13]*tmp1[86]+z[3]*tmp1[87]+z[90]*tmp1[88]+z[80]*tmp1[89]+z[70]*tmp1[90]
            +z[60]*tmp1[91]+z[50]*tmp1[92]+z[40]*tmp1[93]+z[30]*tmp1[94]+z[20]*tmp1[95]+z[10]*tmp1[96];
            tab[i+88*stg_first]=tab[i+88*stg_first]+z[0]*tmp1[0]
            +z[88]*tmp1[1]+z[79]*tmp1[2]+z[70]*tmp1[3]+z[61]*tmp1[4]+z[52]*tmp1[5]+z[43]*tmp1[6]+z[34]*tmp1[7]+z[25]*tmp1[8]+z[16]*tmp1[9]+z[7]*tmp1[10]
            +z[95]*tmp1[11]+z[86]*tmp1[12]+z[77]*tmp1[13]+z[68]*tmp1[14]+z[59]*tmp1[15]+z[50]*tmp1[16]+z[41]*tmp1[17]+z[32]*tmp1[18]+z[23]*tmp1[19]+z[14]*tmp1[20]
            +z[5]*tmp1[21]+z[93]*tmp1[22]+z[84]*tmp1[23]+z[75]*tmp1[24]+z[66]*tmp1[25]+z[57]*tmp1[26]+z[48]*tmp1[27]+z[39]*tmp1[28]+z[30]*tmp1[29]+z[21]*tmp1[30]
            +z[12]*tmp1[31]+z[3]*tmp1[32]+z[91]*tmp1[33]+z[82]*tmp1[34]+z[73]*tmp1[35]+z[64]*tmp1[36]+z[55]*tmp1[37]+z[46]*tmp1[38]+z[37]*tmp1[39]+z[28]*tmp1[40]
            +z[19]*tmp1[41]+z[10]*tmp1[42]+z[1]*tmp1[43]+z[89]*tmp1[44]+z[80]*tmp1[45]+z[71]*tmp1[46]+z[62]*tmp1[47]+z[53]*tmp1[48]+z[44]*tmp1[49]+z[35]*tmp1[50]
            +z[26]*tmp1[51]+z[17]*tmp1[52]+z[8]*tmp1[53]+z[96]*tmp1[54]+z[87]*tmp1[55]+z[78]*tmp1[56]+z[69]*tmp1[57]+z[60]*tmp1[58]+z[51]*tmp1[59]+z[42]*tmp1[60]
            +z[33]*tmp1[61]+z[24]*tmp1[62]+z[15]*tmp1[63]+z[6]*tmp1[64]+z[94]*tmp1[65]+z[85]*tmp1[66]+z[76]*tmp1[67]+z[67]*tmp1[68]+z[58]*tmp1[69]+z[49]*tmp1[70]
            +z[40]*tmp1[71]+z[31]*tmp1[72]+z[22]*tmp1[73]+z[13]*tmp1[74]+z[4]*tmp1[75]+z[92]*tmp1[76]+z[83]*tmp1[77]+z[74]*tmp1[78]+z[65]*tmp1[79]+z[56]*tmp1[80]
            +z[47]*tmp1[81]+z[38]*tmp1[82]+z[29]*tmp1[83]+z[20]*tmp1[84]+z[11]*tmp1[85]+z[2]*tmp1[86]+z[90]*tmp1[87]+z[81]*tmp1[88]+z[72]*tmp1[89]+z[63]*tmp1[90]
            +z[54]*tmp1[91]+z[45]*tmp1[92]+z[36]*tmp1[93]+z[27]*tmp1[94]+z[18]*tmp1[95]+z[9]*tmp1[96];
            tab[i+89*stg_first]=tab[i+89*stg_first]+z[0]*tmp1[0]
            +z[89]*tmp1[1]+z[81]*tmp1[2]+z[73]*tmp1[3]+z[65]*tmp1[4]+z[57]*tmp1[5]+z[49]*tmp1[6]+z[41]*tmp1[7]+z[33]*tmp1[8]+z[25]*tmp1[9]+z[17]*tmp1[10]
            +z[9]*tmp1[11]+z[1]*tmp1[12]+z[90]*tmp1[13]+z[82]*tmp1[14]+z[74]*tmp1[15]+z[66]*tmp1[16]+z[58]*tmp1[17]+z[50]*tmp1[18]+z[42]*tmp1[19]+z[34]*tmp1[20]
            +z[26]*tmp1[21]+z[18]*tmp1[22]+z[10]*tmp1[23]+z[2]*tmp1[24]+z[91]*tmp1[25]+z[83]*tmp1[26]+z[75]*tmp1[27]+z[67]*tmp1[28]+z[59]*tmp1[29]+z[51]*tmp1[30]
            +z[43]*tmp1[31]+z[35]*tmp1[32]+z[27]*tmp1[33]+z[19]*tmp1[34]+z[11]*tmp1[35]+z[3]*tmp1[36]+z[92]*tmp1[37]+z[84]*tmp1[38]+z[76]*tmp1[39]+z[68]*tmp1[40]
            +z[60]*tmp1[41]+z[52]*tmp1[42]+z[44]*tmp1[43]+z[36]*tmp1[44]+z[28]*tmp1[45]+z[20]*tmp1[46]+z[12]*tmp1[47]+z[4]*tmp1[48]+z[93]*tmp1[49]+z[85]*tmp1[50]
            +z[77]*tmp1[51]+z[69]*tmp1[52]+z[61]*tmp1[53]+z[53]*tmp1[54]+z[45]*tmp1[55]+z[37]*tmp1[56]+z[29]*tmp1[57]+z[21]*tmp1[58]+z[13]*tmp1[59]+z[5]*tmp1[60]
            +z[94]*tmp1[61]+z[86]*tmp1[62]+z[78]*tmp1[63]+z[70]*tmp1[64]+z[62]*tmp1[65]+z[54]*tmp1[66]+z[46]*tmp1[67]+z[38]*tmp1[68]+z[30]*tmp1[69]+z[22]*tmp1[70]
            +z[14]*tmp1[71]+z[6]*tmp1[72]+z[95]*tmp1[73]+z[87]*tmp1[74]+z[79]*tmp1[75]+z[71]*tmp1[76]+z[63]*tmp1[77]+z[55]*tmp1[78]+z[47]*tmp1[79]+z[39]*tmp1[80]
            +z[31]*tmp1[81]+z[23]*tmp1[82]+z[15]*tmp1[83]+z[7]*tmp1[84]+z[96]*tmp1[85]+z[88]*tmp1[86]+z[80]*tmp1[87]+z[72]*tmp1[88]+z[64]*tmp1[89]+z[56]*tmp1[90]
            +z[48]*tmp1[91]+z[40]*tmp1[92]+z[32]*tmp1[93]+z[24]*tmp1[94]+z[16]*tmp1[95]+z[8]*tmp1[96];
            tab[i+90*stg_first]=tab[i+90*stg_first]+z[0]*tmp1[0]
            +z[90]*tmp1[1]+z[83]*tmp1[2]+z[76]*tmp1[3]+z[69]*tmp1[4]+z[62]*tmp1[5]+z[55]*tmp1[6]+z[48]*tmp1[7]+z[41]*tmp1[8]+z[34]*tmp1[9]+z[27]*tmp1[10]
            +z[20]*tmp1[11]+z[13]*tmp1[12]+z[6]*tmp1[13]+z[96]*tmp1[14]+z[89]*tmp1[15]+z[82]*tmp1[16]+z[75]*tmp1[17]+z[68]*tmp1[18]+z[61]*tmp1[19]+z[54]*tmp1[20]
            +z[47]*tmp1[21]+z[40]*tmp1[22]+z[33]*tmp1[23]+z[26]*tmp1[24]+z[19]*tmp1[25]+z[12]*tmp1[26]+z[5]*tmp1[27]+z[95]*tmp1[28]+z[88]*tmp1[29]+z[81]*tmp1[30]
            +z[74]*tmp1[31]+z[67]*tmp1[32]+z[60]*tmp1[33]+z[53]*tmp1[34]+z[46]*tmp1[35]+z[39]*tmp1[36]+z[32]*tmp1[37]+z[25]*tmp1[38]+z[18]*tmp1[39]+z[11]*tmp1[40]
            +z[4]*tmp1[41]+z[94]*tmp1[42]+z[87]*tmp1[43]+z[80]*tmp1[44]+z[73]*tmp1[45]+z[66]*tmp1[46]+z[59]*tmp1[47]+z[52]*tmp1[48]+z[45]*tmp1[49]+z[38]*tmp1[50]
            +z[31]*tmp1[51]+z[24]*tmp1[52]+z[17]*tmp1[53]+z[10]*tmp1[54]+z[3]*tmp1[55]+z[93]*tmp1[56]+z[86]*tmp1[57]+z[79]*tmp1[58]+z[72]*tmp1[59]+z[65]*tmp1[60]
            +z[58]*tmp1[61]+z[51]*tmp1[62]+z[44]*tmp1[63]+z[37]*tmp1[64]+z[30]*tmp1[65]+z[23]*tmp1[66]+z[16]*tmp1[67]+z[9]*tmp1[68]+z[2]*tmp1[69]+z[92]*tmp1[70]
            +z[85]*tmp1[71]+z[78]*tmp1[72]+z[71]*tmp1[73]+z[64]*tmp1[74]+z[57]*tmp1[75]+z[50]*tmp1[76]+z[43]*tmp1[77]+z[36]*tmp1[78]+z[29]*tmp1[79]+z[22]*tmp1[80]
            +z[15]*tmp1[81]+z[8]*tmp1[82]+z[1]*tmp1[83]+z[91]*tmp1[84]+z[84]*tmp1[85]+z[77]*tmp1[86]+z[70]*tmp1[87]+z[63]*tmp1[88]+z[56]*tmp1[89]+z[49]*tmp1[90]
            +z[42]*tmp1[91]+z[35]*tmp1[92]+z[28]*tmp1[93]+z[21]*tmp1[94]+z[14]*tmp1[95]+z[7]*tmp1[96];
            tab[i+91*stg_first]=tab[i+91*stg_first]+z[0]*tmp1[0]
            +z[91]*tmp1[1]+z[85]*tmp1[2]+z[79]*tmp1[3]+z[73]*tmp1[4]+z[67]*tmp1[5]+z[61]*tmp1[6]+z[55]*tmp1[7]+z[49]*tmp1[8]+z[43]*tmp1[9]+z[37]*tmp1[10]
            +z[31]*tmp1[11]+z[25]*tmp1[12]+z[19]*tmp1[13]+z[13]*tmp1[14]+z[7]*tmp1[15]+z[1]*tmp1[16]+z[92]*tmp1[17]+z[86]*tmp1[18]+z[80]*tmp1[19]+z[74]*tmp1[20]
            +z[68]*tmp1[21]+z[62]*tmp1[22]+z[56]*tmp1[23]+z[50]*tmp1[24]+z[44]*tmp1[25]+z[38]*tmp1[26]+z[32]*tmp1[27]+z[26]*tmp1[28]+z[20]*tmp1[29]+z[14]*tmp1[30]
            +z[8]*tmp1[31]+z[2]*tmp1[32]+z[93]*tmp1[33]+z[87]*tmp1[34]+z[81]*tmp1[35]+z[75]*tmp1[36]+z[69]*tmp1[37]+z[63]*tmp1[38]+z[57]*tmp1[39]+z[51]*tmp1[40]
            +z[45]*tmp1[41]+z[39]*tmp1[42]+z[33]*tmp1[43]+z[27]*tmp1[44]+z[21]*tmp1[45]+z[15]*tmp1[46]+z[9]*tmp1[47]+z[3]*tmp1[48]+z[94]*tmp1[49]+z[88]*tmp1[50]
            +z[82]*tmp1[51]+z[76]*tmp1[52]+z[70]*tmp1[53]+z[64]*tmp1[54]+z[58]*tmp1[55]+z[52]*tmp1[56]+z[46]*tmp1[57]+z[40]*tmp1[58]+z[34]*tmp1[59]+z[28]*tmp1[60]
            +z[22]*tmp1[61]+z[16]*tmp1[62]+z[10]*tmp1[63]+z[4]*tmp1[64]+z[95]*tmp1[65]+z[89]*tmp1[66]+z[83]*tmp1[67]+z[77]*tmp1[68]+z[71]*tmp1[69]+z[65]*tmp1[70]
            +z[59]*tmp1[71]+z[53]*tmp1[72]+z[47]*tmp1[73]+z[41]*tmp1[74]+z[35]*tmp1[75]+z[29]*tmp1[76]+z[23]*tmp1[77]+z[17]*tmp1[78]+z[11]*tmp1[79]+z[5]*tmp1[80]
            +z[96]*tmp1[81]+z[90]*tmp1[82]+z[84]*tmp1[83]+z[78]*tmp1[84]+z[72]*tmp1[85]+z[66]*tmp1[86]+z[60]*tmp1[87]+z[54]*tmp1[88]+z[48]*tmp1[89]+z[42]*tmp1[90]
            +z[36]*tmp1[91]+z[30]*tmp1[92]+z[24]*tmp1[93]+z[18]*tmp1[94]+z[12]*tmp1[95]+z[6]*tmp1[96];
            tab[i+92*stg_first]=tab[i+92*stg_first]+z[0]*tmp1[0]
            +z[92]*tmp1[1]+z[87]*tmp1[2]+z[82]*tmp1[3]+z[77]*tmp1[4]+z[72]*tmp1[5]+z[67]*tmp1[6]+z[62]*tmp1[7]+z[57]*tmp1[8]+z[52]*tmp1[9]+z[47]*tmp1[10]
            +z[42]*tmp1[11]+z[37]*tmp1[12]+z[32]*tmp1[13]+z[27]*tmp1[14]+z[22]*tmp1[15]+z[17]*tmp1[16]+z[12]*tmp1[17]+z[7]*tmp1[18]+z[2]*tmp1[19]+z[94]*tmp1[20]
            +z[89]*tmp1[21]+z[84]*tmp1[22]+z[79]*tmp1[23]+z[74]*tmp1[24]+z[69]*tmp1[25]+z[64]*tmp1[26]+z[59]*tmp1[27]+z[54]*tmp1[28]+z[49]*tmp1[29]+z[44]*tmp1[30]
            +z[39]*tmp1[31]+z[34]*tmp1[32]+z[29]*tmp1[33]+z[24]*tmp1[34]+z[19]*tmp1[35]+z[14]*tmp1[36]+z[9]*tmp1[37]+z[4]*tmp1[38]+z[96]*tmp1[39]+z[91]*tmp1[40]
            +z[86]*tmp1[41]+z[81]*tmp1[42]+z[76]*tmp1[43]+z[71]*tmp1[44]+z[66]*tmp1[45]+z[61]*tmp1[46]+z[56]*tmp1[47]+z[51]*tmp1[48]+z[46]*tmp1[49]+z[41]*tmp1[50]
            +z[36]*tmp1[51]+z[31]*tmp1[52]+z[26]*tmp1[53]+z[21]*tmp1[54]+z[16]*tmp1[55]+z[11]*tmp1[56]+z[6]*tmp1[57]+z[1]*tmp1[58]+z[93]*tmp1[59]+z[88]*tmp1[60]
            +z[83]*tmp1[61]+z[78]*tmp1[62]+z[73]*tmp1[63]+z[68]*tmp1[64]+z[63]*tmp1[65]+z[58]*tmp1[66]+z[53]*tmp1[67]+z[48]*tmp1[68]+z[43]*tmp1[69]+z[38]*tmp1[70]
            +z[33]*tmp1[71]+z[28]*tmp1[72]+z[23]*tmp1[73]+z[18]*tmp1[74]+z[13]*tmp1[75]+z[8]*tmp1[76]+z[3]*tmp1[77]+z[95]*tmp1[78]+z[90]*tmp1[79]+z[85]*tmp1[80]
            +z[80]*tmp1[81]+z[75]*tmp1[82]+z[70]*tmp1[83]+z[65]*tmp1[84]+z[60]*tmp1[85]+z[55]*tmp1[86]+z[50]*tmp1[87]+z[45]*tmp1[88]+z[40]*tmp1[89]+z[35]*tmp1[90]
            +z[30]*tmp1[91]+z[25]*tmp1[92]+z[20]*tmp1[93]+z[15]*tmp1[94]+z[10]*tmp1[95]+z[5]*tmp1[96];
            tab[i+93*stg_first]=tab[i+93*stg_first]+z[0]*tmp1[0]
            +z[93]*tmp1[1]+z[89]*tmp1[2]+z[85]*tmp1[3]+z[81]*tmp1[4]+z[77]*tmp1[5]+z[73]*tmp1[6]+z[69]*tmp1[7]+z[65]*tmp1[8]+z[61]*tmp1[9]+z[57]*tmp1[10]
            +z[53]*tmp1[11]+z[49]*tmp1[12]+z[45]*tmp1[13]+z[41]*tmp1[14]+z[37]*tmp1[15]+z[33]*tmp1[16]+z[29]*tmp1[17]+z[25]*tmp1[18]+z[21]*tmp1[19]+z[17]*tmp1[20]
            +z[13]*tmp1[21]+z[9]*tmp1[22]+z[5]*tmp1[23]+z[1]*tmp1[24]+z[94]*tmp1[25]+z[90]*tmp1[26]+z[86]*tmp1[27]+z[82]*tmp1[28]+z[78]*tmp1[29]+z[74]*tmp1[30]
            +z[70]*tmp1[31]+z[66]*tmp1[32]+z[62]*tmp1[33]+z[58]*tmp1[34]+z[54]*tmp1[35]+z[50]*tmp1[36]+z[46]*tmp1[37]+z[42]*tmp1[38]+z[38]*tmp1[39]+z[34]*tmp1[40]
            +z[30]*tmp1[41]+z[26]*tmp1[42]+z[22]*tmp1[43]+z[18]*tmp1[44]+z[14]*tmp1[45]+z[10]*tmp1[46]+z[6]*tmp1[47]+z[2]*tmp1[48]+z[95]*tmp1[49]+z[91]*tmp1[50]
            +z[87]*tmp1[51]+z[83]*tmp1[52]+z[79]*tmp1[53]+z[75]*tmp1[54]+z[71]*tmp1[55]+z[67]*tmp1[56]+z[63]*tmp1[57]+z[59]*tmp1[58]+z[55]*tmp1[59]+z[51]*tmp1[60]
            +z[47]*tmp1[61]+z[43]*tmp1[62]+z[39]*tmp1[63]+z[35]*tmp1[64]+z[31]*tmp1[65]+z[27]*tmp1[66]+z[23]*tmp1[67]+z[19]*tmp1[68]+z[15]*tmp1[69]+z[11]*tmp1[70]
            +z[7]*tmp1[71]+z[3]*tmp1[72]+z[96]*tmp1[73]+z[92]*tmp1[74]+z[88]*tmp1[75]+z[84]*tmp1[76]+z[80]*tmp1[77]+z[76]*tmp1[78]+z[72]*tmp1[79]+z[68]*tmp1[80]
            +z[64]*tmp1[81]+z[60]*tmp1[82]+z[56]*tmp1[83]+z[52]*tmp1[84]+z[48]*tmp1[85]+z[44]*tmp1[86]+z[40]*tmp1[87]+z[36]*tmp1[88]+z[32]*tmp1[89]+z[28]*tmp1[90]
            +z[24]*tmp1[91]+z[20]*tmp1[92]+z[16]*tmp1[93]+z[12]*tmp1[94]+z[8]*tmp1[95]+z[4]*tmp1[96];
            tab[i+94*stg_first]=tab[i+94*stg_first]+z[0]*tmp1[0]
            +z[94]*tmp1[1]+z[91]*tmp1[2]+z[88]*tmp1[3]+z[85]*tmp1[4]+z[82]*tmp1[5]+z[79]*tmp1[6]+z[76]*tmp1[7]+z[73]*tmp1[8]+z[70]*tmp1[9]+z[67]*tmp1[10]
            +z[64]*tmp1[11]+z[61]*tmp1[12]+z[58]*tmp1[13]+z[55]*tmp1[14]+z[52]*tmp1[15]+z[49]*tmp1[16]+z[46]*tmp1[17]+z[43]*tmp1[18]+z[40]*tmp1[19]+z[37]*tmp1[20]
            +z[34]*tmp1[21]+z[31]*tmp1[22]+z[28]*tmp1[23]+z[25]*tmp1[24]+z[22]*tmp1[25]+z[19]*tmp1[26]+z[16]*tmp1[27]+z[13]*tmp1[28]+z[10]*tmp1[29]+z[7]*tmp1[30]
            +z[4]*tmp1[31]+z[1]*tmp1[32]+z[95]*tmp1[33]+z[92]*tmp1[34]+z[89]*tmp1[35]+z[86]*tmp1[36]+z[83]*tmp1[37]+z[80]*tmp1[38]+z[77]*tmp1[39]+z[74]*tmp1[40]
            +z[71]*tmp1[41]+z[68]*tmp1[42]+z[65]*tmp1[43]+z[62]*tmp1[44]+z[59]*tmp1[45]+z[56]*tmp1[46]+z[53]*tmp1[47]+z[50]*tmp1[48]+z[47]*tmp1[49]+z[44]*tmp1[50]
            +z[41]*tmp1[51]+z[38]*tmp1[52]+z[35]*tmp1[53]+z[32]*tmp1[54]+z[29]*tmp1[55]+z[26]*tmp1[56]+z[23]*tmp1[57]+z[20]*tmp1[58]+z[17]*tmp1[59]+z[14]*tmp1[60]
            +z[11]*tmp1[61]+z[8]*tmp1[62]+z[5]*tmp1[63]+z[2]*tmp1[64]+z[96]*tmp1[65]+z[93]*tmp1[66]+z[90]*tmp1[67]+z[87]*tmp1[68]+z[84]*tmp1[69]+z[81]*tmp1[70]
            +z[78]*tmp1[71]+z[75]*tmp1[72]+z[72]*tmp1[73]+z[69]*tmp1[74]+z[66]*tmp1[75]+z[63]*tmp1[76]+z[60]*tmp1[77]+z[57]*tmp1[78]+z[54]*tmp1[79]+z[51]*tmp1[80]
            +z[48]*tmp1[81]+z[45]*tmp1[82]+z[42]*tmp1[83]+z[39]*tmp1[84]+z[36]*tmp1[85]+z[33]*tmp1[86]+z[30]*tmp1[87]+z[27]*tmp1[88]+z[24]*tmp1[89]+z[21]*tmp1[90]
            +z[18]*tmp1[91]+z[15]*tmp1[92]+z[12]*tmp1[93]+z[9]*tmp1[94]+z[6]*tmp1[95]+z[3]*tmp1[96];
            tab[i+95*stg_first]=tab[i+95*stg_first]+z[0]*tmp1[0]
            +z[95]*tmp1[1]+z[93]*tmp1[2]+z[91]*tmp1[3]+z[89]*tmp1[4]+z[87]*tmp1[5]+z[85]*tmp1[6]+z[83]*tmp1[7]+z[81]*tmp1[8]+z[79]*tmp1[9]+z[77]*tmp1[10]
            +z[75]*tmp1[11]+z[73]*tmp1[12]+z[71]*tmp1[13]+z[69]*tmp1[14]+z[67]*tmp1[15]+z[65]*tmp1[16]+z[63]*tmp1[17]+z[61]*tmp1[18]+z[59]*tmp1[19]+z[57]*tmp1[20]
            +z[55]*tmp1[21]+z[53]*tmp1[22]+z[51]*tmp1[23]+z[49]*tmp1[24]+z[47]*tmp1[25]+z[45]*tmp1[26]+z[43]*tmp1[27]+z[41]*tmp1[28]+z[39]*tmp1[29]+z[37]*tmp1[30]
            +z[35]*tmp1[31]+z[33]*tmp1[32]+z[31]*tmp1[33]+z[29]*tmp1[34]+z[27]*tmp1[35]+z[25]*tmp1[36]+z[23]*tmp1[37]+z[21]*tmp1[38]+z[19]*tmp1[39]+z[17]*tmp1[40]
            +z[15]*tmp1[41]+z[13]*tmp1[42]+z[11]*tmp1[43]+z[9]*tmp1[44]+z[7]*tmp1[45]+z[5]*tmp1[46]+z[3]*tmp1[47]+z[1]*tmp1[48]+z[96]*tmp1[49]+z[94]*tmp1[50]
            +z[92]*tmp1[51]+z[90]*tmp1[52]+z[88]*tmp1[53]+z[86]*tmp1[54]+z[84]*tmp1[55]+z[82]*tmp1[56]+z[80]*tmp1[57]+z[78]*tmp1[58]+z[76]*tmp1[59]+z[74]*tmp1[60]
            +z[72]*tmp1[61]+z[70]*tmp1[62]+z[68]*tmp1[63]+z[66]*tmp1[64]+z[64]*tmp1[65]+z[62]*tmp1[66]+z[60]*tmp1[67]+z[58]*tmp1[68]+z[56]*tmp1[69]+z[54]*tmp1[70]
            +z[52]*tmp1[71]+z[50]*tmp1[72]+z[48]*tmp1[73]+z[46]*tmp1[74]+z[44]*tmp1[75]+z[42]*tmp1[76]+z[40]*tmp1[77]+z[38]*tmp1[78]+z[36]*tmp1[79]+z[34]*tmp1[80]
            +z[32]*tmp1[81]+z[30]*tmp1[82]+z[28]*tmp1[83]+z[26]*tmp1[84]+z[24]*tmp1[85]+z[22]*tmp1[86]+z[20]*tmp1[87]+z[18]*tmp1[88]+z[16]*tmp1[89]+z[14]*tmp1[90]
            +z[12]*tmp1[91]+z[10]*tmp1[92]+z[8]*tmp1[93]+z[6]*tmp1[94]+z[4]*tmp1[95]+z[2]*tmp1[96];
            tab[i+96*stg_first]=tab[i+96*stg_first]+z[0]*tmp1[0]
            +z[96]*tmp1[1]+z[95]*tmp1[2]+z[94]*tmp1[3]+z[93]*tmp1[4]+z[92]*tmp1[5]+z[91]*tmp1[6]+z[90]*tmp1[7]+z[89]*tmp1[8]+z[88]*tmp1[9]+z[87]*tmp1[10]
            +z[86]*tmp1[11]+z[85]*tmp1[12]+z[84]*tmp1[13]+z[83]*tmp1[14]+z[82]*tmp1[15]+z[81]*tmp1[16]+z[80]*tmp1[17]+z[79]*tmp1[18]+z[78]*tmp1[19]+z[77]*tmp1[20]
            +z[76]*tmp1[21]+z[75]*tmp1[22]+z[74]*tmp1[23]+z[73]*tmp1[24]+z[72]*tmp1[25]+z[71]*tmp1[26]+z[70]*tmp1[27]+z[69]*tmp1[28]+z[68]*tmp1[29]+z[67]*tmp1[30]
            +z[66]*tmp1[31]+z[65]*tmp1[32]+z[64]*tmp1[33]+z[63]*tmp1[34]+z[62]*tmp1[35]+z[61]*tmp1[36]+z[60]*tmp1[37]+z[59]*tmp1[38]+z[58]*tmp1[39]+z[57]*tmp1[40]
            +z[56]*tmp1[41]+z[55]*tmp1[42]+z[54]*tmp1[43]+z[53]*tmp1[44]+z[52]*tmp1[45]+z[51]*tmp1[46]+z[50]*tmp1[47]+z[49]*tmp1[48]+z[48]*tmp1[49]+z[47]*tmp1[50]
            +z[46]*tmp1[51]+z[45]*tmp1[52]+z[44]*tmp1[53]+z[43]*tmp1[54]+z[42]*tmp1[55]+z[41]*tmp1[56]+z[40]*tmp1[57]+z[39]*tmp1[58]+z[38]*tmp1[59]+z[37]*tmp1[60]
            +z[36]*tmp1[61]+z[35]*tmp1[62]+z[34]*tmp1[63]+z[33]*tmp1[64]+z[32]*tmp1[65]+z[31]*tmp1[66]+z[30]*tmp1[67]+z[29]*tmp1[68]+z[28]*tmp1[69]+z[27]*tmp1[70]
            +z[26]*tmp1[71]+z[25]*tmp1[72]+z[24]*tmp1[73]+z[23]*tmp1[74]+z[22]*tmp1[75]+z[21]*tmp1[76]+z[20]*tmp1[77]+z[19]*tmp1[78]+z[18]*tmp1[79]+z[17]*tmp1[80]
            +z[16]*tmp1[81]+z[15]*tmp1[82]+z[14]*tmp1[83]+z[13]*tmp1[84]+z[12]*tmp1[85]+z[11]*tmp1[86]+z[10]*tmp1[87]+z[9]*tmp1[88]+z[8]*tmp1[89]+z[7]*tmp1[90]
            +z[6]*tmp1[91]+z[5]*tmp1[92]+z[4]*tmp1[93]+z[3]*tmp1[94]+z[2]*tmp1[95]+z[1]*tmp1[96];

        }
    }
///////////////////////////////////////////////////

void fun_inverse_fourier_transform_FFT_mixed_radix(int N,std::complex<double> tab[])
{

//author marcin matysek (r)ewertyn.PL 2marcin56@gmail.com
//open-source
    const double pi=3.141592653589793238462;
    std::complex<double> *tab2 = new std::complex<double>[N];    // tab2[]==N

    double tmp2;
    double tmp3;
    double tmp6;
    double tmp7;
    double tmp8;
    double tmp9;
    double tmp10;
    double tmp11;
    double tmp97;

    double fi2=fi; //shift only for stage nr 1

    int i=0;
    int rx5=5,rx4=4,rx3=3,rx2=2,rx7=7,rx11=11,rx97=97;
    int stg[100]={};
    int nb_stages=0;
    int nb1,nb2,nb3,nb4,nb5_stg_previous,stg_first;

    tmp6=2*pi/(N/1);

    tmp2=2*pi/(2/1);
    tmp3=2*pi/(3/1);
    tmp10=2*pi/(4/1);
    tmp8=2*pi/(5/1);
    tmp7=2*pi/(7/1);
    tmp11=2*pi/(11/1);
    tmp97=2*pi/(97/1);



    std::complex<double>  z_rx2[2]={{1,0}};
    std::complex<double>  z_rx3[3]={{1,0}};
    std::complex<double>  z_rx4[2]={{1,0}};
    std::complex<double>  z_rx5[5]={{1,0}};
    std::complex<double>  z_rx7[7]={{1,0}};
    std::complex<double>  z_rx11[11]={{1,0}};
    std::complex<double>  z_rx97[10000]={{1,0}};

//radix 2 fundament
          z_rx2[0].real(cos(0*tmp2));
          z_rx2[0].imag(sin(0*tmp2));
          z_rx2[1].real(cos(1*tmp2));
          z_rx2[1].imag(sin(1*tmp2));
//radix 3 fundament
          z_rx3[0].real(cos(0*tmp3));
          z_rx3[0].imag(sin(0*tmp3));
          z_rx3[1].real(cos(1*tmp3));
          z_rx3[1].imag(sin(1*tmp3));
          z_rx3[2].real(cos(2*tmp3));
          z_rx3[2].imag(sin(2*tmp3));
//radix 4 fundament
          z_rx4[0].real(cos(0*tmp10));
          z_rx4[0].imag(sin(0*tmp10));
          z_rx4[1].real(cos(1*tmp10));
          z_rx4[1].imag(sin(1*tmp10));
//radix 5 fundament
          z_rx5[0].real(cos(0*tmp8));
          z_rx5[0].imag(sin(0*tmp8));
          z_rx5[1].real(cos(1*tmp8));
          z_rx5[1].imag(sin(1*tmp8));
          z_rx5[2].real(cos(2*tmp8));
          z_rx5[2].imag(sin(2*tmp8));
          z_rx5[3].real(cos(3*tmp8));
          z_rx5[3].imag(sin(3*tmp8));
          z_rx5[4].real(cos(4*tmp8));
          z_rx5[4].imag(sin(4*tmp8));
//radix 7 fundament
          z_rx7[0].real(cos(0*tmp7));
          z_rx7[0].imag(sin(0*tmp7));
          z_rx7[1].real(cos(1*tmp7));
          z_rx7[1].imag(sin(1*tmp7));
          z_rx7[2].real(cos(2*tmp7));
          z_rx7[2].imag(sin(2*tmp7));
          z_rx7[3].real(cos(3*tmp7));
          z_rx7[3].imag(sin(3*tmp7));
          z_rx7[4].real(cos(4*tmp7));
          z_rx7[4].imag(sin(4*tmp7));
          z_rx7[5].real(cos(5*tmp7));
          z_rx7[5].imag(sin(5*tmp7));
          z_rx7[6].real(cos(6*tmp7));
          z_rx7[6].imag(sin(6*tmp7));
//radix 11 fundament
          z_rx11[0].real(cos(0*tmp11));
          z_rx11[0].imag(sin(0*tmp11));
          z_rx11[1].real(cos(1*tmp11));
          z_rx11[1].imag(sin(1*tmp11));
          z_rx11[2].real(cos(2*tmp11));
          z_rx11[2].imag(sin(2*tmp11));
          z_rx11[3].real(cos(3*tmp11));
          z_rx11[3].imag(sin(3*tmp11));
          z_rx11[4].real(cos(4*tmp11));
          z_rx11[4].imag(sin(4*tmp11));
          z_rx11[5].real(cos(5*tmp11));
          z_rx11[5].imag(sin(5*tmp11));
          z_rx11[6].real(cos(6*tmp11));
          z_rx11[6].imag(sin(6*tmp11));
          z_rx11[7].real(cos(7*tmp11));
          z_rx11[7].imag(sin(7*tmp11));
          z_rx11[8].real(cos(8*tmp11));
          z_rx11[8].imag(sin(8*tmp11));
          z_rx11[9].real(cos(9*tmp11));
          z_rx11[9].imag(sin(9*tmp11));
          z_rx11[10].real(cos(10*tmp11));
          z_rx11[10].imag(sin(10*tmp11));


        for (int i2=0;i2<97;i2++)
        {
            for (int j2=0;j2<97;j2++)
            {
                z_rx97[i2*j2].real(cos(i2*j2*tmp97));
                z_rx97[i2*j2].imag(sin(i2*j2*tmp97));
            }
        }


/*
for(int j=0;j<102;j++)
{
    for(int i=0;i<102;i++)
    {
      if(((fabs(round(z_rx11[j].imag()*1000)/1000-round(z_rx11[i].imag()*1000)/1000)<0.001)
         &&(fabs(round(z_rx11[j].real()*1000)/1000-round(z_rx11[i].real()*1000)/1000)<0.001))
         ||((fabs(round(z_rx11[j].imag()*1000)/1000+round(z_rx11[i].imag()*1000)/1000)<0.001)
         &&(fabs(round(z_rx11[j].real()*1000)/1000+round(z_rx11[i].real()*1000)/1000)<0.001)))
         {
             cout<<j<<" "<<round(z_rx11[j].real()*1000)/1000<<" "<<round(z_rx11[j].imag()*1000)/1000<<"   ";
             cout<<i<<" "<<round(z_rx11[i].real()*1000)/1000<<" "<<round(z_rx11[i].imag()*1000)/1000<<endl;
         }
             //cout<<endl;
    }
    //system("pause");
}
*/

	nb_stages=radix_base(N,stg);

        if(nb_stages>=1){cout<<"N= "<<N<<endl;}
        for(int m=1;m<=nb_stages;m++)
        {
        cout <<"stage:"<<m<<" = radix-"<<stg[m] << endl;
        }
        cout << endl;



    if(nb_stages>=1)
    {
        stg_first=N/stg[1];
        if(stg[1]==2)
        {
            fun_inverse_fourier_transform_FFT_radix_2_stg_first(tab,stg_first,fi2,z_rx2);
        }
        else if(stg[1]==3)
        {
            fun_inverse_fourier_transform_FFT_radix_3_stg_first(tab,stg_first,fi2,z_rx3);
        }
        else if(stg[1]==4)
        {
            fun_inverse_fourier_transform_FFT_radix_4_stg_first(tab,stg_first,fi2,z_rx4);
        }
        else if(stg[1]==5)
        {
            fun_inverse_fourier_transform_FFT_radix_5_stg_first(tab,stg_first,fi2,z_rx5);
        }
        else if(stg[1]==7)
        {
            fun_inverse_fourier_transform_FFT_radix_7_stg_first(tab,stg_first,fi2,z_rx7);
        }
        else if(stg[1]==11)
        {
            fun_inverse_fourier_transform_FFT_radix_11_stg_first(tab,stg_first,fi2,z_rx11);
        }
        else if(stg[1]==97)
        {
            fun_inverse_fourier_transform_FFT_radix_97_stg_first(tab,stg_first,fi2,z_rx97);
        }
        else{}
        nb1=N;
        nb4=1;
        for(int i=0;i<nb_stages-1;i++)
        {
            nb1=nb1/stg[0+i];
            nb2=nb1/stg[1+i];
            nb3=nb2/stg[2+i];
            nb4=nb4*stg[0+i];
            nb5_stg_previous=stg[1+i];

            if(stg[i+2]==2)
            {
                fun_inverse_fourier_transform_FFT_radix_2_stg_rest(tab,nb1,nb2,nb3,nb4,nb5_stg_previous,tmp6,z_rx2);
            }
            else if(stg[i+2]==3)
            {
                fun_inverse_fourier_transform_FFT_radix_3_stg_rest(tab,nb1,nb2,nb3,nb4,nb5_stg_previous,tmp6,z_rx3);
            }
            else if(stg[i+2]==4)
            {
                fun_inverse_fourier_transform_FFT_radix_4_stg_rest(tab,nb1,nb2,nb3,nb4,nb5_stg_previous,tmp6,z_rx4);
            }
            else if(stg[i+2]==5)
            {
                fun_inverse_fourier_transform_FFT_radix_5_stg_rest(tab,nb1,nb2,nb3,nb4,nb5_stg_previous,tmp6,z_rx5);
            }
            else if(stg[i+2]==7)
            {
                fun_inverse_fourier_transform_FFT_radix_7_stg_rest(tab,nb1,nb2,nb3,nb4,nb5_stg_previous,tmp6,z_rx7);
            }
            else if(stg[i+2]==11)
            {
                fun_inverse_fourier_transform_FFT_radix_11_stg_rest(tab,nb1,nb2,nb3,nb4,nb5_stg_previous,tmp6,z_rx11);
            }
            else if(stg[i+2]==97)
            {
                fun_inverse_fourier_transform_FFT_radix_97_stg_rest(tab,nb1,nb2,nb3,nb4,nb5_stg_previous,tmp6,z_rx97);
            }
            else{}

        }
    }

//new:
    for(int j=0;j<N;j++)
    {
     tab[j].real(tab[j].real()*0.5);
     tab[j].imag(tab[j].imag()*0.5);
    }
    delete [] tab2;
}
///////////////////////////////////////////////////


  void fun_inverse_fourier_transform_FFT_radix_2_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp6,std::complex<double> z[])
    {
        std::complex<double> tmp1,tmp2,tmp40,tmp50;
        std::complex<double>  w[2]={{1,0}};
        double tmp100=0.0;
        double tmp200=0.0;
        double tmp300=0.0;
        int nb_tmp1=0.0;
        int nb_tmp2=0.0;
        int nb_tmp3=0.0;

        for(int b=0;b<nb5_stg_previous;b=b+1)
        {
            tmp300=nb4*b*tmp6;
            tmp100=nb3*tmp300;
            for(int j=0;j<nb3;j=j+1)
            {
                tmp200=j*tmp300;
                w[0].real( cos(0*tmp100+tmp200));
                w[0].imag(sin(0*tmp100+tmp200));
                w[1].real( cos(1*tmp100+tmp200));
                //w[1].imag(sin(nb4*b*(1*nb3+j)*tmp5));
                w[1].imag(sin(1*tmp100+tmp200));

                nb_tmp1=b*nb2+j;
                for(int i=0;i<nb4;i=i+1)
                {
                    nb_tmp2=i*nb1;
                    nb_tmp3=nb_tmp1+nb_tmp2;
                    tmp1=w[0]*tab[nb_tmp3+0*nb3];
                    tmp2=w[1]*tab[nb_tmp3+1*nb3];

                    tmp40=z[0]*(tmp1+tmp2);

                    tab[nb_tmp3+0*nb3]=tmp40;
                    tab[nb_tmp3+1*nb3]=z[0]*tmp1+z[1]*tmp2;
                }
            }
        }
    }
///////////////////////////////////////////////////////////////


    void fun_inverse_fourier_transform_FFT_radix_3_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp6,std::complex<double> z[])
    {
        std::complex<double> tmp1,tmp2,tmp3,tmp4,tmp5;
        std::complex<double>  w[3]={{1,0}};
        double tmp100=0.0;
        double tmp200=0.0;
        double tmp300=0.0;
        int nb_tmp1=0.0;
        int nb_tmp2=0.0;
        int nb_tmp3=0.0;

        for(int b=0;b<nb5_stg_previous;b=b+1)
        {
            tmp300=nb4*b*tmp6;
            tmp100=nb3*tmp300;
            for(int j=0;j<nb3;j=j+1)
            {
                tmp200=j*tmp300;
                w[0].real(cos(0*tmp100+tmp200));
                w[0].imag(sin(0*tmp100+tmp200));
                w[1].real(cos(1*tmp100+tmp200));
                w[1].imag(sin(1*tmp100+tmp200));
                //w[2].real(cos(nb4*b*(2*nb3+j)*tmp6));
                w[2].real(cos(2*tmp100+tmp200));
                w[2].imag(sin(2*tmp100+tmp200));

                nb_tmp1=b*nb2+j;
                for(int i=0;i<nb4;i=i+1)
                {
                    nb_tmp2=i*nb1;
                    nb_tmp3=nb_tmp1+nb_tmp2;
                  tmp1=w[0]*tab[nb_tmp3+0*nb3];
                  tmp2=w[1]*tab[nb_tmp3+1*nb3];
                  tmp3=w[2]*tab[nb_tmp3+2*nb3];

                  tmp4=z[0]*tmp1;
                  tmp5=z[0]*(tmp2+tmp3);

                 //radix-3
                  tab[nb_tmp3+0*nb3]  = tmp4+tmp5;
                  tab[nb_tmp3+1*nb3]  = tmp4+z[1]*tmp2+z[2]*tmp3;
                  tab[nb_tmp3+2*nb3]  = tmp4+z[2]*tmp2+z[1]*tmp3;
                }
            }
        }
    }
///////////////////////////////////////////////////////////////


    void fun_inverse_fourier_transform_FFT_radix_4_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp6,std::complex<double> z[])
    {
        std::complex<double> tmp1,tmp2,tmp3,tmp4;
        std::complex<double> tmp101,tmp102,tmp103,tmp104;
        std::complex<double> tmp111,tmp112,tmp113,tmp114;
        std::complex<double>  w[4]={{1,0}};
        double tmp100=0.0;
        double tmp200=0.0;
        double tmp300=0.0;
        int nb_tmp1=0.0;
        int nb_tmp2=0.0;
        int nb_tmp3=0.0;

        for(int b=0;b<nb5_stg_previous;b=b+1)
        {
            tmp300=nb4*b*tmp6;
            tmp100=nb3*tmp300;
            for(int j=0;j<nb3;j=j+1)
            {
              tmp200=j*tmp300;
              w[0].real(cos(0*tmp100+tmp200));
              w[0].imag(sin(0*tmp100+tmp200));
              w[1].real(cos(1*tmp100+tmp200));
              w[1].imag(sin(1*tmp100+tmp200));
              w[2].real(cos(2*tmp100+tmp200));
              w[2].imag(sin(2*tmp100+tmp200));
              //w[3].real(cos(nb4*b*(3*nb3+j)*tmp5));
              w[3].real(cos(3*tmp100+tmp200));
              w[3].imag(sin(3*tmp100+tmp200));

                nb_tmp1=b*nb2+j;
                for(int i=0;i<nb4;i=i+1)
                {
                    nb_tmp2=i*nb1;
                    nb_tmp3=nb_tmp1+nb_tmp2;
                    tmp1=w[0]*tab[nb_tmp3+0*nb3];
                    tmp2=w[1]*tab[nb_tmp3+1*nb3];
                    tmp3=w[2]*tab[nb_tmp3+2*nb3];
                    tmp4=w[3]*tab[nb_tmp3+3*nb3];


                    tmp101=tmp2-tmp4;
                    tmp102=tmp2+tmp4;
                    tmp103=tmp1-tmp3;
                    tmp104=tmp1+tmp3;

                    tmp111=z[0]*(tmp104+tmp102);
                    tmp112=z[0]*tmp103;
                    tmp113=z[0]*(tmp104-tmp102);
                    tmp114=z[0]*tmp103;

                    //radix-4
                    tab[nb_tmp3+0*nb3]   =tmp111;
                    tab[nb_tmp3+1*nb3]   =tmp112+z[1]*tmp101;
                    tab[nb_tmp3+2*nb3]   =tmp113;
                    tab[nb_tmp3+3*nb3]   =tmp114-z[1]*tmp101;
                }
            }
        }
    }
/////////////////////////////////////////////////////////////////

   void fun_inverse_fourier_transform_FFT_radix_5_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp6,std::complex<double> z[])
    {
        std::complex<double> tmp1,tmp2,tmp3,tmp4,tmp5,tmp10,tmp20;
        std::complex<double>  w[5]={{1,0}};
        double tmp100=0.0;
        double tmp200=0.0;
        double tmp300=0.0;
        int nb_tmp1=0.0;
        int nb_tmp2=0.0;
        int nb_tmp3=0.0;

        for(int b=0;b<nb5_stg_previous;b=b+1)
        {
            tmp300=nb4*b*tmp6;
            tmp100=nb3*tmp300;
            for(int j=0;j<nb3;j=j+1)
            {
                tmp200=j*tmp300;
                w[0].real(cos(0*tmp100+tmp200));
                w[0].imag(sin(0*tmp100+tmp200));
                w[1].real(cos(1*tmp100+tmp200));
                w[1].imag(sin(1*tmp100+tmp200));
                //w[2].real(cos(nb4*b*(2*nb3+j)*tmp6));
                w[2].real(cos(2*tmp100+tmp200));
                w[2].imag(sin(2*tmp100+tmp200));
                w[3].real(cos(3*tmp100+tmp200));
                w[3].imag(sin(3*tmp100+tmp200));
                w[4].real(cos(4*tmp100+tmp200));
                w[4].imag(sin(4*tmp100+tmp200));

                nb_tmp1=b*nb2+j;
                for(int i=0;i<nb4;i=i+1)
                {
                    nb_tmp2=i*nb1;
                    nb_tmp3=nb_tmp1+nb_tmp2;
                  tmp1=w[0]*tab[nb_tmp3+0*nb3];
                  tmp2=w[1]*tab[nb_tmp3+1*nb3];
                  tmp3=w[2]*tab[nb_tmp3+2*nb3];
                  tmp4=w[3]*tab[nb_tmp3+3*nb3];
                  tmp5=w[4]*tab[nb_tmp3+4*nb3];

                  tmp10=z[0]*tmp1;
                  tmp20=z[0]*(tmp1+tmp2+tmp3+tmp4+tmp5);

                 //radix-5
                  tab[nb_tmp3+0*nb3] =tmp20;
                  tab[nb_tmp3+1*nb3] =tmp10+z[1]*tmp2+z[2]*tmp3+z[3]*tmp4+z[4]*tmp5;
                  tab[nb_tmp3+2*nb3] =tmp10+z[2]*tmp2+z[4]*tmp3+z[1]*tmp4+z[3]*tmp5;
                  tab[nb_tmp3+3*nb3] =tmp10+z[3]*tmp2+z[1]*tmp3+z[4]*tmp4+z[2]*tmp5;
                  tab[nb_tmp3+4*nb3] =tmp10+z[4]*tmp2+z[3]*tmp3+z[2]*tmp4+z[1]*tmp5;
                }
            }
        }
    }
///////////////////////////////////////////////////////////////////////

   void fun_inverse_fourier_transform_FFT_radix_7_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp99,std::complex<double> z[])
    {
        std::complex<double> tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp10,tmp20;
        std::complex<double>  w[7]={{1,0}};
        double tmp100=0.0;
        double tmp200=0.0;
        double tmp300=0.0;
        int nb_tmp1=0.0;
        int nb_tmp2=0.0;
        int nb_tmp3=0.0;

        for(int b=0;b<nb5_stg_previous;b=b+1)
        {
            tmp300=nb4*b*tmp99;
            tmp100=nb3*tmp300;
            for(int j=0;j<nb3;j=j+1)
            {
                tmp200=j*tmp300;
                w[0].real(cos(0*tmp100+tmp200));
                w[0].imag(sin(0*tmp100+tmp200));
                w[1].real(cos(1*tmp100+tmp200));
                w[1].imag(sin(1*tmp100+tmp200));
                //w[2].real(cos(nb4*b*(2*nb3+j)*tmp99));
                w[2].real(cos(2*tmp100+tmp200));
                w[2].imag(sin(2*tmp100+tmp200));
                w[3].real(cos(3*tmp100+tmp200));
                w[3].imag(sin(3*tmp100+tmp200));
                w[4].real(cos(4*tmp100+tmp200));
                w[4].imag(sin(4*tmp100+tmp200));
                w[5].real(cos(5*tmp100+tmp200));
                w[5].imag(sin(5*tmp100+tmp200));
                w[6].real(cos(6*tmp100+tmp200));
                w[6].imag(sin(6*tmp100+tmp200));

                nb_tmp1=b*nb2+j;
                for(int i=0;i<nb4;i=i+1)
                {
                    nb_tmp2=i*nb1;
                    nb_tmp3=nb_tmp1+nb_tmp2;
                  tmp1=w[0]*tab[nb_tmp3+0*nb3];
                  tmp2=w[1]*tab[nb_tmp3+1*nb3];
                  tmp3=w[2]*tab[nb_tmp3+2*nb3];
                  tmp4=w[3]*tab[nb_tmp3+3*nb3];
                  tmp5=w[4]*tab[nb_tmp3+4*nb3];
                  tmp6=w[5]*tab[nb_tmp3+5*nb3];
                  tmp7=w[6]*tab[nb_tmp3+6*nb3];

                  tmp10=z[0]*tmp1;
                  tmp20=z[0]*(tmp1+tmp2+tmp3+tmp4+tmp5+tmp6+tmp7);

                 //radix-7
                  tab[nb_tmp3+0*nb3] =tmp20;
                  tab[nb_tmp3+1*nb3] =tmp10+z[1]*tmp2+z[2]*tmp3+z[3]*tmp4+z[4]*tmp5+z[5]*tmp6+z[6]*tmp7;
                  tab[nb_tmp3+2*nb3] =tmp10+z[2]*tmp2+z[4]*tmp3+z[6]*tmp4+z[1]*tmp5+z[3]*tmp6+z[5]*tmp7;
                  tab[nb_tmp3+3*nb3] =tmp10+z[3]*tmp2+z[6]*tmp3+z[2]*tmp4+z[5]*tmp5+z[1]*tmp6+z[4]*tmp7;
                  tab[nb_tmp3+4*nb3] =tmp10+z[4]*tmp2+z[1]*tmp3+z[5]*tmp4+z[2]*tmp5+z[6]*tmp6+z[3]*tmp7;
                  tab[nb_tmp3+5*nb3] =tmp10+z[5]*tmp2+z[3]*tmp3+z[1]*tmp4+z[6]*tmp5+z[4]*tmp6+z[2]*tmp7;
                  tab[nb_tmp3+6*nb3] =tmp10+z[6]*tmp2+z[5]*tmp3+z[4]*tmp4+z[3]*tmp5+z[2]*tmp6+z[1]*tmp7;
                }
            }
        }
    }
///////////////////////////////////////////////////////////////////////

   void fun_inverse_fourier_transform_FFT_radix_11_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp99,std::complex<double> z[])
    {
        std::complex<double> tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10,tmp11,tmp30,tmp20;
        std::complex<double>  w[11]={{1,0}};
        double tmp100=0.0;
        double tmp200=0.0;
        double tmp300=0.0;
        int nb_tmp1=0.0;
        int nb_tmp2=0.0;
        int nb_tmp3=0.0;

        for(int b=0;b<nb5_stg_previous;b=b+1)
        {
            tmp300=nb4*b*tmp99;
            tmp100=nb3*tmp300;
            for(int j=0;j<nb3;j=j+1)
            {
                tmp200=j*tmp300;
                w[0].real(cos(0*tmp100+tmp200));
                w[0].imag(sin(0*tmp100+tmp200));
                w[1].real(cos(1*tmp100+tmp200));
                w[1].imag(sin(1*tmp100+tmp200));
                //w[2].real(cos(nb4*b*(2*nb3+j)*tmp99));
                w[2].real(cos(2*tmp100+tmp200));
                w[2].imag(sin(2*tmp100+tmp200));
                w[3].real(cos(3*tmp100+tmp200));
                w[3].imag(sin(3*tmp100+tmp200));
                w[4].real(cos(4*tmp100+tmp200));
                w[4].imag(sin(4*tmp100+tmp200));
                w[5].real(cos(5*tmp100+tmp200));
                w[5].imag(sin(5*tmp100+tmp200));
                w[6].real(cos(6*tmp100+tmp200));
                w[6].imag(sin(6*tmp100+tmp200));
                w[7].real(cos(7*tmp100+tmp200));
                w[7].imag(sin(7*tmp100+tmp200));
                w[8].real(cos(8*tmp100+tmp200));
                w[8].imag(sin(8*tmp100+tmp200));
                w[9].real(cos(9*tmp100+tmp200));
                w[9].imag(sin(9*tmp100+tmp200));
                w[10].real(cos(10*tmp100+tmp200));
                w[10].imag(sin(10*tmp100+tmp200));

                nb_tmp1=b*nb2+j;
                for(int i=0;i<nb4;i=i+1)
                {
                    nb_tmp2=i*nb1;
                    nb_tmp3=nb_tmp1+nb_tmp2;
                  tmp1=w[0]*tab[nb_tmp3+0*nb3];
                  tmp2=w[1]*tab[nb_tmp3+1*nb3];
                  tmp3=w[2]*tab[nb_tmp3+2*nb3];
                  tmp4=w[3]*tab[nb_tmp3+3*nb3];
                  tmp5=w[4]*tab[nb_tmp3+4*nb3];
                  tmp6=w[5]*tab[nb_tmp3+5*nb3];
                  tmp7=w[6]*tab[nb_tmp3+6*nb3];
                  tmp8=w[7]*tab[nb_tmp3+7*nb3];
                  tmp9=w[8]*tab[nb_tmp3+8*nb3];
                  tmp10=w[9]*tab[nb_tmp3+9*nb3];
                  tmp11=w[10]*tab[nb_tmp3+10*nb3];

                  tmp30=z[0]*tmp1;
                  tmp20=z[0]*(tmp1+tmp2+tmp3+tmp4+tmp5+tmp6+tmp7+tmp8+tmp9+tmp10+tmp11);

                 //radix-11
                  tab[nb_tmp3+0*nb3] =tmp20;
                  tab[nb_tmp3+1*nb3] =tmp30+z[1]*tmp2+z[2]*tmp3+z[3]*tmp4+z[4]*tmp5+z[5]*tmp6+z[6]*tmp7+z[7]*tmp8+z[8]*tmp9+z[9]*tmp10+z[10]*tmp11;
                  tab[nb_tmp3+2*nb3] =tmp30+z[2]*tmp2+z[4]*tmp3+z[6]*tmp4+z[8]*tmp5+z[10]*tmp6+z[1]*tmp7+z[3]*tmp8+z[5]*tmp9+z[7]*tmp10+z[9]*tmp11;
                  tab[nb_tmp3+3*nb3] =tmp30+z[3]*tmp2+z[6]*tmp3+z[9]*tmp4+z[1]*tmp5+z[4]*tmp6+z[7]*tmp7+z[10]*tmp8+z[2]*tmp9+z[5]*tmp10+z[8]*tmp11;
                  tab[nb_tmp3+4*nb3] =tmp30+z[4]*tmp2+z[8]*tmp3+z[1]*tmp4+z[5]*tmp5+z[9]*tmp6+z[2]*tmp7+z[6]*tmp8+z[10]*tmp9+z[3]*tmp10+z[7]*tmp11;
                  tab[nb_tmp3+5*nb3] =tmp30+z[5]*tmp2+z[10]*tmp3+z[4]*tmp4+z[9]*tmp5+z[3]*tmp6+z[8]*tmp7+z[2]*tmp8+z[7]*tmp9+z[1]*tmp10+z[6]*tmp11;
                  tab[nb_tmp3+6*nb3] =tmp30+z[6]*tmp2+z[1]*tmp3+z[7]*tmp4+z[2]*tmp5+z[8]*tmp6+z[3]*tmp7+z[9]*tmp8+z[4]*tmp9+z[10]*tmp10+z[5]*tmp11;
                  tab[nb_tmp3+7*nb3] =tmp30+z[7]*tmp2+z[3]*tmp3+z[10]*tmp4+z[6]*tmp5+z[2]*tmp6+z[9]*tmp7+z[5]*tmp8+z[1]*tmp9+z[8]*tmp10+z[4]*tmp11;
                  tab[nb_tmp3+8*nb3] =tmp30+z[8]*tmp2+z[5]*tmp3+z[2]*tmp4+z[10]*tmp5+z[7]*tmp6+z[4]*tmp7+z[1]*tmp8+z[9]*tmp9+z[6]*tmp10+z[3]*tmp11;
                  tab[nb_tmp3+9*nb3] =tmp30+z[9]*tmp2+z[7]*tmp3+z[5]*tmp4+z[3]*tmp5+z[1]*tmp6+z[10]*tmp7+z[8]*tmp8+z[6]*tmp9+z[4]*tmp10+z[2]*tmp11;
                  tab[nb_tmp3+10*nb3] =tmp30+z[10]*tmp2+z[9]*tmp3+z[8]*tmp4+z[7]*tmp5+z[6]*tmp6+z[5]*tmp7+z[4]*tmp8+z[3]*tmp9+z[2]*tmp10+z[1]*tmp11;
                }

            }
        }
    }
///////////////////////////////////////////////////////////////////////

void fun_inverse_fourier_transform_FFT_radix_97_stg_rest(std::complex<double> tab[],int nb1,int nb2,int nb3,int nb4,int nb5_stg_previous,double tmp99,std::complex<double> z[])
    {
        std::complex<double> tmp1[97],tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10,tmp11,tmp30,tmp20;
        std::complex<double>  w[97]={{1,0}};
        double tmp100=0.0;
        double tmp200=0.0;
        double tmp300=0.0;
        int nb_tmp1=0.0;
        int nb_tmp2=0.0;
        int nb_tmp3=0.0;

        for(int b=0;b<nb5_stg_previous;b=b+1)
        {
            tmp300=nb4*b*tmp99;
            tmp100=nb3*tmp300;
            for(int j=0;j<nb3;j=j+1)
            {
                tmp200=j*tmp300;
                for(int j9=0;j9<97;j9++)
                {
                     w[j9].real(cos(j9*tmp100+tmp200));
                     w[j9].imag(sin(j9*tmp100+tmp200));
                }

                nb_tmp1=b*nb2+j;
                for(int i=0;i<nb4;i=i+1)
                {
                        nb_tmp2=i*nb1;
                        nb_tmp3=nb_tmp1+nb_tmp2;
                for(int j9=0;j9<97;j9++)
                {
                       tmp1[j9]=w[j9]*tab[nb_tmp3+j9*nb3];
                }


                         //radix-97
                   for(int i9=0;i9<97;i9++)
                  {
                    for(int j9=0;j9<97;j9++)
                    {
                        tab[nb_tmp3+i9*nb3].real(0);
                        tab[nb_tmp3+i9*nb3].imag(0);
                    }
                  }

                    tab[nb_tmp3+0*nb3]=tab[nb_tmp3+0*nb3]+z[0]*tmp1[0]
                    +z[0]*tmp1[1]+z[0]*tmp1[2]+z[0]*tmp1[3]+z[0]*tmp1[4]+z[0]*tmp1[5]+z[0]*tmp1[6]+z[0]*tmp1[7]+z[0]*tmp1[8]+z[0]*tmp1[9]+z[0]*tmp1[10]
                    +z[0]*tmp1[11]+z[0]*tmp1[12]+z[0]*tmp1[13]+z[0]*tmp1[14]+z[0]*tmp1[15]+z[0]*tmp1[16]+z[0]*tmp1[17]+z[0]*tmp1[18]+z[0]*tmp1[19]+z[0]*tmp1[20]
                    +z[0]*tmp1[21]+z[0]*tmp1[22]+z[0]*tmp1[23]+z[0]*tmp1[24]+z[0]*tmp1[25]+z[0]*tmp1[26]+z[0]*tmp1[27]+z[0]*tmp1[28]+z[0]*tmp1[29]+z[0]*tmp1[30]
                    +z[0]*tmp1[31]+z[0]*tmp1[32]+z[0]*tmp1[33]+z[0]*tmp1[34]+z[0]*tmp1[35]+z[0]*tmp1[36]+z[0]*tmp1[37]+z[0]*tmp1[38]+z[0]*tmp1[39]+z[0]*tmp1[40]
                    +z[0]*tmp1[41]+z[0]*tmp1[42]+z[0]*tmp1[43]+z[0]*tmp1[44]+z[0]*tmp1[45]+z[0]*tmp1[46]+z[0]*tmp1[47]+z[0]*tmp1[48]+z[0]*tmp1[49]+z[0]*tmp1[50]
                    +z[0]*tmp1[51]+z[0]*tmp1[52]+z[0]*tmp1[53]+z[0]*tmp1[54]+z[0]*tmp1[55]+z[0]*tmp1[56]+z[0]*tmp1[57]+z[0]*tmp1[58]+z[0]*tmp1[59]+z[0]*tmp1[60]
                    +z[0]*tmp1[61]+z[0]*tmp1[62]+z[0]*tmp1[63]+z[0]*tmp1[64]+z[0]*tmp1[65]+z[0]*tmp1[66]+z[0]*tmp1[67]+z[0]*tmp1[68]+z[0]*tmp1[69]+z[0]*tmp1[70]
                    +z[0]*tmp1[71]+z[0]*tmp1[72]+z[0]*tmp1[73]+z[0]*tmp1[74]+z[0]*tmp1[75]+z[0]*tmp1[76]+z[0]*tmp1[77]+z[0]*tmp1[78]+z[0]*tmp1[79]+z[0]*tmp1[80]
                    +z[0]*tmp1[81]+z[0]*tmp1[82]+z[0]*tmp1[83]+z[0]*tmp1[84]+z[0]*tmp1[85]+z[0]*tmp1[86]+z[0]*tmp1[87]+z[0]*tmp1[88]+z[0]*tmp1[89]+z[0]*tmp1[90]
                    +z[0]*tmp1[91]+z[0]*tmp1[92]+z[0]*tmp1[93]+z[0]*tmp1[94]+z[0]*tmp1[95]+z[0]*tmp1[96];
                    tab[nb_tmp3+1*nb3]=tab[nb_tmp3+1*nb3]+z[0]*tmp1[0]
                    +z[1]*tmp1[1]+z[2]*tmp1[2]+z[3]*tmp1[3]+z[4]*tmp1[4]+z[5]*tmp1[5]+z[6]*tmp1[6]+z[7]*tmp1[7]+z[8]*tmp1[8]+z[9]*tmp1[9]+z[10]*tmp1[10]
                    +z[11]*tmp1[11]+z[12]*tmp1[12]+z[13]*tmp1[13]+z[14]*tmp1[14]+z[15]*tmp1[15]+z[16]*tmp1[16]+z[17]*tmp1[17]+z[18]*tmp1[18]+z[19]*tmp1[19]+z[20]*tmp1[20]
                    +z[21]*tmp1[21]+z[22]*tmp1[22]+z[23]*tmp1[23]+z[24]*tmp1[24]+z[25]*tmp1[25]+z[26]*tmp1[26]+z[27]*tmp1[27]+z[28]*tmp1[28]+z[29]*tmp1[29]+z[30]*tmp1[30]
                    +z[31]*tmp1[31]+z[32]*tmp1[32]+z[33]*tmp1[33]+z[34]*tmp1[34]+z[35]*tmp1[35]+z[36]*tmp1[36]+z[37]*tmp1[37]+z[38]*tmp1[38]+z[39]*tmp1[39]+z[40]*tmp1[40]
                    +z[41]*tmp1[41]+z[42]*tmp1[42]+z[43]*tmp1[43]+z[44]*tmp1[44]+z[45]*tmp1[45]+z[46]*tmp1[46]+z[47]*tmp1[47]+z[48]*tmp1[48]+z[49]*tmp1[49]+z[50]*tmp1[50]
                    +z[51]*tmp1[51]+z[52]*tmp1[52]+z[53]*tmp1[53]+z[54]*tmp1[54]+z[55]*tmp1[55]+z[56]*tmp1[56]+z[57]*tmp1[57]+z[58]*tmp1[58]+z[59]*tmp1[59]+z[60]*tmp1[60]
                    +z[61]*tmp1[61]+z[62]*tmp1[62]+z[63]*tmp1[63]+z[64]*tmp1[64]+z[65]*tmp1[65]+z[66]*tmp1[66]+z[67]*tmp1[67]+z[68]*tmp1[68]+z[69]*tmp1[69]+z[70]*tmp1[70]
                    +z[71]*tmp1[71]+z[72]*tmp1[72]+z[73]*tmp1[73]+z[74]*tmp1[74]+z[75]*tmp1[75]+z[76]*tmp1[76]+z[77]*tmp1[77]+z[78]*tmp1[78]+z[79]*tmp1[79]+z[80]*tmp1[80]
                    +z[81]*tmp1[81]+z[82]*tmp1[82]+z[83]*tmp1[83]+z[84]*tmp1[84]+z[85]*tmp1[85]+z[86]*tmp1[86]+z[87]*tmp1[87]+z[88]*tmp1[88]+z[89]*tmp1[89]+z[90]*tmp1[90]
                    +z[91]*tmp1[91]+z[92]*tmp1[92]+z[93]*tmp1[93]+z[94]*tmp1[94]+z[95]*tmp1[95]+z[96]*tmp1[96];
                    tab[nb_tmp3+2*nb3]=tab[nb_tmp3+2*nb3]+z[0]*tmp1[0]
                    +z[2]*tmp1[1]+z[4]*tmp1[2]+z[6]*tmp1[3]+z[8]*tmp1[4]+z[10]*tmp1[5]+z[12]*tmp1[6]+z[14]*tmp1[7]+z[16]*tmp1[8]+z[18]*tmp1[9]+z[20]*tmp1[10]
                    +z[22]*tmp1[11]+z[24]*tmp1[12]+z[26]*tmp1[13]+z[28]*tmp1[14]+z[30]*tmp1[15]+z[32]*tmp1[16]+z[34]*tmp1[17]+z[36]*tmp1[18]+z[38]*tmp1[19]+z[40]*tmp1[20]
                    +z[42]*tmp1[21]+z[44]*tmp1[22]+z[46]*tmp1[23]+z[48]*tmp1[24]+z[50]*tmp1[25]+z[52]*tmp1[26]+z[54]*tmp1[27]+z[56]*tmp1[28]+z[58]*tmp1[29]+z[60]*tmp1[30]
                    +z[62]*tmp1[31]+z[64]*tmp1[32]+z[66]*tmp1[33]+z[68]*tmp1[34]+z[70]*tmp1[35]+z[72]*tmp1[36]+z[74]*tmp1[37]+z[76]*tmp1[38]+z[78]*tmp1[39]+z[80]*tmp1[40]
                    +z[82]*tmp1[41]+z[84]*tmp1[42]+z[86]*tmp1[43]+z[88]*tmp1[44]+z[90]*tmp1[45]+z[92]*tmp1[46]+z[94]*tmp1[47]+z[96]*tmp1[48]+z[1]*tmp1[49]+z[3]*tmp1[50]
                    +z[5]*tmp1[51]+z[7]*tmp1[52]+z[9]*tmp1[53]+z[11]*tmp1[54]+z[13]*tmp1[55]+z[15]*tmp1[56]+z[17]*tmp1[57]+z[19]*tmp1[58]+z[21]*tmp1[59]+z[23]*tmp1[60]
                    +z[25]*tmp1[61]+z[27]*tmp1[62]+z[29]*tmp1[63]+z[31]*tmp1[64]+z[33]*tmp1[65]+z[35]*tmp1[66]+z[37]*tmp1[67]+z[39]*tmp1[68]+z[41]*tmp1[69]+z[43]*tmp1[70]
                    +z[45]*tmp1[71]+z[47]*tmp1[72]+z[49]*tmp1[73]+z[51]*tmp1[74]+z[53]*tmp1[75]+z[55]*tmp1[76]+z[57]*tmp1[77]+z[59]*tmp1[78]+z[61]*tmp1[79]+z[63]*tmp1[80]
                    +z[65]*tmp1[81]+z[67]*tmp1[82]+z[69]*tmp1[83]+z[71]*tmp1[84]+z[73]*tmp1[85]+z[75]*tmp1[86]+z[77]*tmp1[87]+z[79]*tmp1[88]+z[81]*tmp1[89]+z[83]*tmp1[90]
                    +z[85]*tmp1[91]+z[87]*tmp1[92]+z[89]*tmp1[93]+z[91]*tmp1[94]+z[93]*tmp1[95]+z[95]*tmp1[96];
                    tab[nb_tmp3+3*nb3]=tab[nb_tmp3+3*nb3]+z[0]*tmp1[0]
                    +z[3]*tmp1[1]+z[6]*tmp1[2]+z[9]*tmp1[3]+z[12]*tmp1[4]+z[15]*tmp1[5]+z[18]*tmp1[6]+z[21]*tmp1[7]+z[24]*tmp1[8]+z[27]*tmp1[9]+z[30]*tmp1[10]
                    +z[33]*tmp1[11]+z[36]*tmp1[12]+z[39]*tmp1[13]+z[42]*tmp1[14]+z[45]*tmp1[15]+z[48]*tmp1[16]+z[51]*tmp1[17]+z[54]*tmp1[18]+z[57]*tmp1[19]+z[60]*tmp1[20]
                    +z[63]*tmp1[21]+z[66]*tmp1[22]+z[69]*tmp1[23]+z[72]*tmp1[24]+z[75]*tmp1[25]+z[78]*tmp1[26]+z[81]*tmp1[27]+z[84]*tmp1[28]+z[87]*tmp1[29]+z[90]*tmp1[30]
                    +z[93]*tmp1[31]+z[96]*tmp1[32]+z[2]*tmp1[33]+z[5]*tmp1[34]+z[8]*tmp1[35]+z[11]*tmp1[36]+z[14]*tmp1[37]+z[17]*tmp1[38]+z[20]*tmp1[39]+z[23]*tmp1[40]
                    +z[26]*tmp1[41]+z[29]*tmp1[42]+z[32]*tmp1[43]+z[35]*tmp1[44]+z[38]*tmp1[45]+z[41]*tmp1[46]+z[44]*tmp1[47]+z[47]*tmp1[48]+z[50]*tmp1[49]+z[53]*tmp1[50]
                    +z[56]*tmp1[51]+z[59]*tmp1[52]+z[62]*tmp1[53]+z[65]*tmp1[54]+z[68]*tmp1[55]+z[71]*tmp1[56]+z[74]*tmp1[57]+z[77]*tmp1[58]+z[80]*tmp1[59]+z[83]*tmp1[60]
                    +z[86]*tmp1[61]+z[89]*tmp1[62]+z[92]*tmp1[63]+z[95]*tmp1[64]+z[1]*tmp1[65]+z[4]*tmp1[66]+z[7]*tmp1[67]+z[10]*tmp1[68]+z[13]*tmp1[69]+z[16]*tmp1[70]
                    +z[19]*tmp1[71]+z[22]*tmp1[72]+z[25]*tmp1[73]+z[28]*tmp1[74]+z[31]*tmp1[75]+z[34]*tmp1[76]+z[37]*tmp1[77]+z[40]*tmp1[78]+z[43]*tmp1[79]+z[46]*tmp1[80]
                    +z[49]*tmp1[81]+z[52]*tmp1[82]+z[55]*tmp1[83]+z[58]*tmp1[84]+z[61]*tmp1[85]+z[64]*tmp1[86]+z[67]*tmp1[87]+z[70]*tmp1[88]+z[73]*tmp1[89]+z[76]*tmp1[90]
                    +z[79]*tmp1[91]+z[82]*tmp1[92]+z[85]*tmp1[93]+z[88]*tmp1[94]+z[91]*tmp1[95]+z[94]*tmp1[96];
                    tab[nb_tmp3+4*nb3]=tab[nb_tmp3+4*nb3]+z[0]*tmp1[0]
                    +z[4]*tmp1[1]+z[8]*tmp1[2]+z[12]*tmp1[3]+z[16]*tmp1[4]+z[20]*tmp1[5]+z[24]*tmp1[6]+z[28]*tmp1[7]+z[32]*tmp1[8]+z[36]*tmp1[9]+z[40]*tmp1[10]
                    +z[44]*tmp1[11]+z[48]*tmp1[12]+z[52]*tmp1[13]+z[56]*tmp1[14]+z[60]*tmp1[15]+z[64]*tmp1[16]+z[68]*tmp1[17]+z[72]*tmp1[18]+z[76]*tmp1[19]+z[80]*tmp1[20]
                    +z[84]*tmp1[21]+z[88]*tmp1[22]+z[92]*tmp1[23]+z[96]*tmp1[24]+z[3]*tmp1[25]+z[7]*tmp1[26]+z[11]*tmp1[27]+z[15]*tmp1[28]+z[19]*tmp1[29]+z[23]*tmp1[30]
                    +z[27]*tmp1[31]+z[31]*tmp1[32]+z[35]*tmp1[33]+z[39]*tmp1[34]+z[43]*tmp1[35]+z[47]*tmp1[36]+z[51]*tmp1[37]+z[55]*tmp1[38]+z[59]*tmp1[39]+z[63]*tmp1[40]
                    +z[67]*tmp1[41]+z[71]*tmp1[42]+z[75]*tmp1[43]+z[79]*tmp1[44]+z[83]*tmp1[45]+z[87]*tmp1[46]+z[91]*tmp1[47]+z[95]*tmp1[48]+z[2]*tmp1[49]+z[6]*tmp1[50]
                    +z[10]*tmp1[51]+z[14]*tmp1[52]+z[18]*tmp1[53]+z[22]*tmp1[54]+z[26]*tmp1[55]+z[30]*tmp1[56]+z[34]*tmp1[57]+z[38]*tmp1[58]+z[42]*tmp1[59]+z[46]*tmp1[60]
                    +z[50]*tmp1[61]+z[54]*tmp1[62]+z[58]*tmp1[63]+z[62]*tmp1[64]+z[66]*tmp1[65]+z[70]*tmp1[66]+z[74]*tmp1[67]+z[78]*tmp1[68]+z[82]*tmp1[69]+z[86]*tmp1[70]
                    +z[90]*tmp1[71]+z[94]*tmp1[72]+z[1]*tmp1[73]+z[5]*tmp1[74]+z[9]*tmp1[75]+z[13]*tmp1[76]+z[17]*tmp1[77]+z[21]*tmp1[78]+z[25]*tmp1[79]+z[29]*tmp1[80]
                    +z[33]*tmp1[81]+z[37]*tmp1[82]+z[41]*tmp1[83]+z[45]*tmp1[84]+z[49]*tmp1[85]+z[53]*tmp1[86]+z[57]*tmp1[87]+z[61]*tmp1[88]+z[65]*tmp1[89]+z[69]*tmp1[90]
                    +z[73]*tmp1[91]+z[77]*tmp1[92]+z[81]*tmp1[93]+z[85]*tmp1[94]+z[89]*tmp1[95]+z[93]*tmp1[96];
                    tab[nb_tmp3+5*nb3]=tab[nb_tmp3+5*nb3]+z[0]*tmp1[0]
                    +z[5]*tmp1[1]+z[10]*tmp1[2]+z[15]*tmp1[3]+z[20]*tmp1[4]+z[25]*tmp1[5]+z[30]*tmp1[6]+z[35]*tmp1[7]+z[40]*tmp1[8]+z[45]*tmp1[9]+z[50]*tmp1[10]
                    +z[55]*tmp1[11]+z[60]*tmp1[12]+z[65]*tmp1[13]+z[70]*tmp1[14]+z[75]*tmp1[15]+z[80]*tmp1[16]+z[85]*tmp1[17]+z[90]*tmp1[18]+z[95]*tmp1[19]+z[3]*tmp1[20]
                    +z[8]*tmp1[21]+z[13]*tmp1[22]+z[18]*tmp1[23]+z[23]*tmp1[24]+z[28]*tmp1[25]+z[33]*tmp1[26]+z[38]*tmp1[27]+z[43]*tmp1[28]+z[48]*tmp1[29]+z[53]*tmp1[30]
                    +z[58]*tmp1[31]+z[63]*tmp1[32]+z[68]*tmp1[33]+z[73]*tmp1[34]+z[78]*tmp1[35]+z[83]*tmp1[36]+z[88]*tmp1[37]+z[93]*tmp1[38]+z[1]*tmp1[39]+z[6]*tmp1[40]
                    +z[11]*tmp1[41]+z[16]*tmp1[42]+z[21]*tmp1[43]+z[26]*tmp1[44]+z[31]*tmp1[45]+z[36]*tmp1[46]+z[41]*tmp1[47]+z[46]*tmp1[48]+z[51]*tmp1[49]+z[56]*tmp1[50]
                    +z[61]*tmp1[51]+z[66]*tmp1[52]+z[71]*tmp1[53]+z[76]*tmp1[54]+z[81]*tmp1[55]+z[86]*tmp1[56]+z[91]*tmp1[57]+z[96]*tmp1[58]+z[4]*tmp1[59]+z[9]*tmp1[60]
                    +z[14]*tmp1[61]+z[19]*tmp1[62]+z[24]*tmp1[63]+z[29]*tmp1[64]+z[34]*tmp1[65]+z[39]*tmp1[66]+z[44]*tmp1[67]+z[49]*tmp1[68]+z[54]*tmp1[69]+z[59]*tmp1[70]
                    +z[64]*tmp1[71]+z[69]*tmp1[72]+z[74]*tmp1[73]+z[79]*tmp1[74]+z[84]*tmp1[75]+z[89]*tmp1[76]+z[94]*tmp1[77]+z[2]*tmp1[78]+z[7]*tmp1[79]+z[12]*tmp1[80]
                    +z[17]*tmp1[81]+z[22]*tmp1[82]+z[27]*tmp1[83]+z[32]*tmp1[84]+z[37]*tmp1[85]+z[42]*tmp1[86]+z[47]*tmp1[87]+z[52]*tmp1[88]+z[57]*tmp1[89]+z[62]*tmp1[90]
                    +z[67]*tmp1[91]+z[72]*tmp1[92]+z[77]*tmp1[93]+z[82]*tmp1[94]+z[87]*tmp1[95]+z[92]*tmp1[96];
                    tab[nb_tmp3+6*nb3]=tab[nb_tmp3+6*nb3]+z[0]*tmp1[0]
                    +z[6]*tmp1[1]+z[12]*tmp1[2]+z[18]*tmp1[3]+z[24]*tmp1[4]+z[30]*tmp1[5]+z[36]*tmp1[6]+z[42]*tmp1[7]+z[48]*tmp1[8]+z[54]*tmp1[9]+z[60]*tmp1[10]
                    +z[66]*tmp1[11]+z[72]*tmp1[12]+z[78]*tmp1[13]+z[84]*tmp1[14]+z[90]*tmp1[15]+z[96]*tmp1[16]+z[5]*tmp1[17]+z[11]*tmp1[18]+z[17]*tmp1[19]+z[23]*tmp1[20]
                    +z[29]*tmp1[21]+z[35]*tmp1[22]+z[41]*tmp1[23]+z[47]*tmp1[24]+z[53]*tmp1[25]+z[59]*tmp1[26]+z[65]*tmp1[27]+z[71]*tmp1[28]+z[77]*tmp1[29]+z[83]*tmp1[30]
                    +z[89]*tmp1[31]+z[95]*tmp1[32]+z[4]*tmp1[33]+z[10]*tmp1[34]+z[16]*tmp1[35]+z[22]*tmp1[36]+z[28]*tmp1[37]+z[34]*tmp1[38]+z[40]*tmp1[39]+z[46]*tmp1[40]
                    +z[52]*tmp1[41]+z[58]*tmp1[42]+z[64]*tmp1[43]+z[70]*tmp1[44]+z[76]*tmp1[45]+z[82]*tmp1[46]+z[88]*tmp1[47]+z[94]*tmp1[48]+z[3]*tmp1[49]+z[9]*tmp1[50]
                    +z[15]*tmp1[51]+z[21]*tmp1[52]+z[27]*tmp1[53]+z[33]*tmp1[54]+z[39]*tmp1[55]+z[45]*tmp1[56]+z[51]*tmp1[57]+z[57]*tmp1[58]+z[63]*tmp1[59]+z[69]*tmp1[60]
                    +z[75]*tmp1[61]+z[81]*tmp1[62]+z[87]*tmp1[63]+z[93]*tmp1[64]+z[2]*tmp1[65]+z[8]*tmp1[66]+z[14]*tmp1[67]+z[20]*tmp1[68]+z[26]*tmp1[69]+z[32]*tmp1[70]
                    +z[38]*tmp1[71]+z[44]*tmp1[72]+z[50]*tmp1[73]+z[56]*tmp1[74]+z[62]*tmp1[75]+z[68]*tmp1[76]+z[74]*tmp1[77]+z[80]*tmp1[78]+z[86]*tmp1[79]+z[92]*tmp1[80]
                    +z[1]*tmp1[81]+z[7]*tmp1[82]+z[13]*tmp1[83]+z[19]*tmp1[84]+z[25]*tmp1[85]+z[31]*tmp1[86]+z[37]*tmp1[87]+z[43]*tmp1[88]+z[49]*tmp1[89]+z[55]*tmp1[90]
                    +z[61]*tmp1[91]+z[67]*tmp1[92]+z[73]*tmp1[93]+z[79]*tmp1[94]+z[85]*tmp1[95]+z[91]*tmp1[96];
                    tab[nb_tmp3+7*nb3]=tab[nb_tmp3+7*nb3]+z[0]*tmp1[0]
                    +z[7]*tmp1[1]+z[14]*tmp1[2]+z[21]*tmp1[3]+z[28]*tmp1[4]+z[35]*tmp1[5]+z[42]*tmp1[6]+z[49]*tmp1[7]+z[56]*tmp1[8]+z[63]*tmp1[9]+z[70]*tmp1[10]
                    +z[77]*tmp1[11]+z[84]*tmp1[12]+z[91]*tmp1[13]+z[1]*tmp1[14]+z[8]*tmp1[15]+z[15]*tmp1[16]+z[22]*tmp1[17]+z[29]*tmp1[18]+z[36]*tmp1[19]+z[43]*tmp1[20]
                    +z[50]*tmp1[21]+z[57]*tmp1[22]+z[64]*tmp1[23]+z[71]*tmp1[24]+z[78]*tmp1[25]+z[85]*tmp1[26]+z[92]*tmp1[27]+z[2]*tmp1[28]+z[9]*tmp1[29]+z[16]*tmp1[30]
                    +z[23]*tmp1[31]+z[30]*tmp1[32]+z[37]*tmp1[33]+z[44]*tmp1[34]+z[51]*tmp1[35]+z[58]*tmp1[36]+z[65]*tmp1[37]+z[72]*tmp1[38]+z[79]*tmp1[39]+z[86]*tmp1[40]
                    +z[93]*tmp1[41]+z[3]*tmp1[42]+z[10]*tmp1[43]+z[17]*tmp1[44]+z[24]*tmp1[45]+z[31]*tmp1[46]+z[38]*tmp1[47]+z[45]*tmp1[48]+z[52]*tmp1[49]+z[59]*tmp1[50]
                    +z[66]*tmp1[51]+z[73]*tmp1[52]+z[80]*tmp1[53]+z[87]*tmp1[54]+z[94]*tmp1[55]+z[4]*tmp1[56]+z[11]*tmp1[57]+z[18]*tmp1[58]+z[25]*tmp1[59]+z[32]*tmp1[60]
                    +z[39]*tmp1[61]+z[46]*tmp1[62]+z[53]*tmp1[63]+z[60]*tmp1[64]+z[67]*tmp1[65]+z[74]*tmp1[66]+z[81]*tmp1[67]+z[88]*tmp1[68]+z[95]*tmp1[69]+z[5]*tmp1[70]
                    +z[12]*tmp1[71]+z[19]*tmp1[72]+z[26]*tmp1[73]+z[33]*tmp1[74]+z[40]*tmp1[75]+z[47]*tmp1[76]+z[54]*tmp1[77]+z[61]*tmp1[78]+z[68]*tmp1[79]+z[75]*tmp1[80]
                    +z[82]*tmp1[81]+z[89]*tmp1[82]+z[96]*tmp1[83]+z[6]*tmp1[84]+z[13]*tmp1[85]+z[20]*tmp1[86]+z[27]*tmp1[87]+z[34]*tmp1[88]+z[41]*tmp1[89]+z[48]*tmp1[90]
                    +z[55]*tmp1[91]+z[62]*tmp1[92]+z[69]*tmp1[93]+z[76]*tmp1[94]+z[83]*tmp1[95]+z[90]*tmp1[96];
                    tab[nb_tmp3+8*nb3]=tab[nb_tmp3+8*nb3]+z[0]*tmp1[0]
                    +z[8]*tmp1[1]+z[16]*tmp1[2]+z[24]*tmp1[3]+z[32]*tmp1[4]+z[40]*tmp1[5]+z[48]*tmp1[6]+z[56]*tmp1[7]+z[64]*tmp1[8]+z[72]*tmp1[9]+z[80]*tmp1[10]
                    +z[88]*tmp1[11]+z[96]*tmp1[12]+z[7]*tmp1[13]+z[15]*tmp1[14]+z[23]*tmp1[15]+z[31]*tmp1[16]+z[39]*tmp1[17]+z[47]*tmp1[18]+z[55]*tmp1[19]+z[63]*tmp1[20]
                    +z[71]*tmp1[21]+z[79]*tmp1[22]+z[87]*tmp1[23]+z[95]*tmp1[24]+z[6]*tmp1[25]+z[14]*tmp1[26]+z[22]*tmp1[27]+z[30]*tmp1[28]+z[38]*tmp1[29]+z[46]*tmp1[30]
                    +z[54]*tmp1[31]+z[62]*tmp1[32]+z[70]*tmp1[33]+z[78]*tmp1[34]+z[86]*tmp1[35]+z[94]*tmp1[36]+z[5]*tmp1[37]+z[13]*tmp1[38]+z[21]*tmp1[39]+z[29]*tmp1[40]
                    +z[37]*tmp1[41]+z[45]*tmp1[42]+z[53]*tmp1[43]+z[61]*tmp1[44]+z[69]*tmp1[45]+z[77]*tmp1[46]+z[85]*tmp1[47]+z[93]*tmp1[48]+z[4]*tmp1[49]+z[12]*tmp1[50]
                    +z[20]*tmp1[51]+z[28]*tmp1[52]+z[36]*tmp1[53]+z[44]*tmp1[54]+z[52]*tmp1[55]+z[60]*tmp1[56]+z[68]*tmp1[57]+z[76]*tmp1[58]+z[84]*tmp1[59]+z[92]*tmp1[60]
                    +z[3]*tmp1[61]+z[11]*tmp1[62]+z[19]*tmp1[63]+z[27]*tmp1[64]+z[35]*tmp1[65]+z[43]*tmp1[66]+z[51]*tmp1[67]+z[59]*tmp1[68]+z[67]*tmp1[69]+z[75]*tmp1[70]
                    +z[83]*tmp1[71]+z[91]*tmp1[72]+z[2]*tmp1[73]+z[10]*tmp1[74]+z[18]*tmp1[75]+z[26]*tmp1[76]+z[34]*tmp1[77]+z[42]*tmp1[78]+z[50]*tmp1[79]+z[58]*tmp1[80]
                    +z[66]*tmp1[81]+z[74]*tmp1[82]+z[82]*tmp1[83]+z[90]*tmp1[84]+z[1]*tmp1[85]+z[9]*tmp1[86]+z[17]*tmp1[87]+z[25]*tmp1[88]+z[33]*tmp1[89]+z[41]*tmp1[90]
                    +z[49]*tmp1[91]+z[57]*tmp1[92]+z[65]*tmp1[93]+z[73]*tmp1[94]+z[81]*tmp1[95]+z[89]*tmp1[96];
                    tab[nb_tmp3+9*nb3]=tab[nb_tmp3+9*nb3]+z[0]*tmp1[0]
                    +z[9]*tmp1[1]+z[18]*tmp1[2]+z[27]*tmp1[3]+z[36]*tmp1[4]+z[45]*tmp1[5]+z[54]*tmp1[6]+z[63]*tmp1[7]+z[72]*tmp1[8]+z[81]*tmp1[9]+z[90]*tmp1[10]
                    +z[2]*tmp1[11]+z[11]*tmp1[12]+z[20]*tmp1[13]+z[29]*tmp1[14]+z[38]*tmp1[15]+z[47]*tmp1[16]+z[56]*tmp1[17]+z[65]*tmp1[18]+z[74]*tmp1[19]+z[83]*tmp1[20]
                    +z[92]*tmp1[21]+z[4]*tmp1[22]+z[13]*tmp1[23]+z[22]*tmp1[24]+z[31]*tmp1[25]+z[40]*tmp1[26]+z[49]*tmp1[27]+z[58]*tmp1[28]+z[67]*tmp1[29]+z[76]*tmp1[30]
                    +z[85]*tmp1[31]+z[94]*tmp1[32]+z[6]*tmp1[33]+z[15]*tmp1[34]+z[24]*tmp1[35]+z[33]*tmp1[36]+z[42]*tmp1[37]+z[51]*tmp1[38]+z[60]*tmp1[39]+z[69]*tmp1[40]
                    +z[78]*tmp1[41]+z[87]*tmp1[42]+z[96]*tmp1[43]+z[8]*tmp1[44]+z[17]*tmp1[45]+z[26]*tmp1[46]+z[35]*tmp1[47]+z[44]*tmp1[48]+z[53]*tmp1[49]+z[62]*tmp1[50]
                    +z[71]*tmp1[51]+z[80]*tmp1[52]+z[89]*tmp1[53]+z[1]*tmp1[54]+z[10]*tmp1[55]+z[19]*tmp1[56]+z[28]*tmp1[57]+z[37]*tmp1[58]+z[46]*tmp1[59]+z[55]*tmp1[60]
                    +z[64]*tmp1[61]+z[73]*tmp1[62]+z[82]*tmp1[63]+z[91]*tmp1[64]+z[3]*tmp1[65]+z[12]*tmp1[66]+z[21]*tmp1[67]+z[30]*tmp1[68]+z[39]*tmp1[69]+z[48]*tmp1[70]
                    +z[57]*tmp1[71]+z[66]*tmp1[72]+z[75]*tmp1[73]+z[84]*tmp1[74]+z[93]*tmp1[75]+z[5]*tmp1[76]+z[14]*tmp1[77]+z[23]*tmp1[78]+z[32]*tmp1[79]+z[41]*tmp1[80]
                    +z[50]*tmp1[81]+z[59]*tmp1[82]+z[68]*tmp1[83]+z[77]*tmp1[84]+z[86]*tmp1[85]+z[95]*tmp1[86]+z[7]*tmp1[87]+z[16]*tmp1[88]+z[25]*tmp1[89]+z[34]*tmp1[90]
                    +z[43]*tmp1[91]+z[52]*tmp1[92]+z[61]*tmp1[93]+z[70]*tmp1[94]+z[79]*tmp1[95]+z[88]*tmp1[96];
                    tab[nb_tmp3+10*nb3]=tab[nb_tmp3+10*nb3]+z[0]*tmp1[0]
                    +z[10]*tmp1[1]+z[20]*tmp1[2]+z[30]*tmp1[3]+z[40]*tmp1[4]+z[50]*tmp1[5]+z[60]*tmp1[6]+z[70]*tmp1[7]+z[80]*tmp1[8]+z[90]*tmp1[9]+z[3]*tmp1[10]
                    +z[13]*tmp1[11]+z[23]*tmp1[12]+z[33]*tmp1[13]+z[43]*tmp1[14]+z[53]*tmp1[15]+z[63]*tmp1[16]+z[73]*tmp1[17]+z[83]*tmp1[18]+z[93]*tmp1[19]+z[6]*tmp1[20]
                    +z[16]*tmp1[21]+z[26]*tmp1[22]+z[36]*tmp1[23]+z[46]*tmp1[24]+z[56]*tmp1[25]+z[66]*tmp1[26]+z[76]*tmp1[27]+z[86]*tmp1[28]+z[96]*tmp1[29]+z[9]*tmp1[30]
                    +z[19]*tmp1[31]+z[29]*tmp1[32]+z[39]*tmp1[33]+z[49]*tmp1[34]+z[59]*tmp1[35]+z[69]*tmp1[36]+z[79]*tmp1[37]+z[89]*tmp1[38]+z[2]*tmp1[39]+z[12]*tmp1[40]
                    +z[22]*tmp1[41]+z[32]*tmp1[42]+z[42]*tmp1[43]+z[52]*tmp1[44]+z[62]*tmp1[45]+z[72]*tmp1[46]+z[82]*tmp1[47]+z[92]*tmp1[48]+z[5]*tmp1[49]+z[15]*tmp1[50]
                    +z[25]*tmp1[51]+z[35]*tmp1[52]+z[45]*tmp1[53]+z[55]*tmp1[54]+z[65]*tmp1[55]+z[75]*tmp1[56]+z[85]*tmp1[57]+z[95]*tmp1[58]+z[8]*tmp1[59]+z[18]*tmp1[60]
                    +z[28]*tmp1[61]+z[38]*tmp1[62]+z[48]*tmp1[63]+z[58]*tmp1[64]+z[68]*tmp1[65]+z[78]*tmp1[66]+z[88]*tmp1[67]+z[1]*tmp1[68]+z[11]*tmp1[69]+z[21]*tmp1[70]
                    +z[31]*tmp1[71]+z[41]*tmp1[72]+z[51]*tmp1[73]+z[61]*tmp1[74]+z[71]*tmp1[75]+z[81]*tmp1[76]+z[91]*tmp1[77]+z[4]*tmp1[78]+z[14]*tmp1[79]+z[24]*tmp1[80]
                    +z[34]*tmp1[81]+z[44]*tmp1[82]+z[54]*tmp1[83]+z[64]*tmp1[84]+z[74]*tmp1[85]+z[84]*tmp1[86]+z[94]*tmp1[87]+z[7]*tmp1[88]+z[17]*tmp1[89]+z[27]*tmp1[90]
                    +z[37]*tmp1[91]+z[47]*tmp1[92]+z[57]*tmp1[93]+z[67]*tmp1[94]+z[77]*tmp1[95]+z[87]*tmp1[96];
                    tab[nb_tmp3+11*nb3]=tab[nb_tmp3+11*nb3]+z[0]*tmp1[0]
                    +z[11]*tmp1[1]+z[22]*tmp1[2]+z[33]*tmp1[3]+z[44]*tmp1[4]+z[55]*tmp1[5]+z[66]*tmp1[6]+z[77]*tmp1[7]+z[88]*tmp1[8]+z[2]*tmp1[9]+z[13]*tmp1[10]
                    +z[24]*tmp1[11]+z[35]*tmp1[12]+z[46]*tmp1[13]+z[57]*tmp1[14]+z[68]*tmp1[15]+z[79]*tmp1[16]+z[90]*tmp1[17]+z[4]*tmp1[18]+z[15]*tmp1[19]+z[26]*tmp1[20]
                    +z[37]*tmp1[21]+z[48]*tmp1[22]+z[59]*tmp1[23]+z[70]*tmp1[24]+z[81]*tmp1[25]+z[92]*tmp1[26]+z[6]*tmp1[27]+z[17]*tmp1[28]+z[28]*tmp1[29]+z[39]*tmp1[30]
                    +z[50]*tmp1[31]+z[61]*tmp1[32]+z[72]*tmp1[33]+z[83]*tmp1[34]+z[94]*tmp1[35]+z[8]*tmp1[36]+z[19]*tmp1[37]+z[30]*tmp1[38]+z[41]*tmp1[39]+z[52]*tmp1[40]
                    +z[63]*tmp1[41]+z[74]*tmp1[42]+z[85]*tmp1[43]+z[96]*tmp1[44]+z[10]*tmp1[45]+z[21]*tmp1[46]+z[32]*tmp1[47]+z[43]*tmp1[48]+z[54]*tmp1[49]+z[65]*tmp1[50]
                    +z[76]*tmp1[51]+z[87]*tmp1[52]+z[1]*tmp1[53]+z[12]*tmp1[54]+z[23]*tmp1[55]+z[34]*tmp1[56]+z[45]*tmp1[57]+z[56]*tmp1[58]+z[67]*tmp1[59]+z[78]*tmp1[60]
                    +z[89]*tmp1[61]+z[3]*tmp1[62]+z[14]*tmp1[63]+z[25]*tmp1[64]+z[36]*tmp1[65]+z[47]*tmp1[66]+z[58]*tmp1[67]+z[69]*tmp1[68]+z[80]*tmp1[69]+z[91]*tmp1[70]
                    +z[5]*tmp1[71]+z[16]*tmp1[72]+z[27]*tmp1[73]+z[38]*tmp1[74]+z[49]*tmp1[75]+z[60]*tmp1[76]+z[71]*tmp1[77]+z[82]*tmp1[78]+z[93]*tmp1[79]+z[7]*tmp1[80]
                    +z[18]*tmp1[81]+z[29]*tmp1[82]+z[40]*tmp1[83]+z[51]*tmp1[84]+z[62]*tmp1[85]+z[73]*tmp1[86]+z[84]*tmp1[87]+z[95]*tmp1[88]+z[9]*tmp1[89]+z[20]*tmp1[90]
                    +z[31]*tmp1[91]+z[42]*tmp1[92]+z[53]*tmp1[93]+z[64]*tmp1[94]+z[75]*tmp1[95]+z[86]*tmp1[96];
                    tab[nb_tmp3+12*nb3]=tab[nb_tmp3+12*nb3]+z[0]*tmp1[0]
                    +z[12]*tmp1[1]+z[24]*tmp1[2]+z[36]*tmp1[3]+z[48]*tmp1[4]+z[60]*tmp1[5]+z[72]*tmp1[6]+z[84]*tmp1[7]+z[96]*tmp1[8]+z[11]*tmp1[9]+z[23]*tmp1[10]
                    +z[35]*tmp1[11]+z[47]*tmp1[12]+z[59]*tmp1[13]+z[71]*tmp1[14]+z[83]*tmp1[15]+z[95]*tmp1[16]+z[10]*tmp1[17]+z[22]*tmp1[18]+z[34]*tmp1[19]+z[46]*tmp1[20]
                    +z[58]*tmp1[21]+z[70]*tmp1[22]+z[82]*tmp1[23]+z[94]*tmp1[24]+z[9]*tmp1[25]+z[21]*tmp1[26]+z[33]*tmp1[27]+z[45]*tmp1[28]+z[57]*tmp1[29]+z[69]*tmp1[30]
                    +z[81]*tmp1[31]+z[93]*tmp1[32]+z[8]*tmp1[33]+z[20]*tmp1[34]+z[32]*tmp1[35]+z[44]*tmp1[36]+z[56]*tmp1[37]+z[68]*tmp1[38]+z[80]*tmp1[39]+z[92]*tmp1[40]
                    +z[7]*tmp1[41]+z[19]*tmp1[42]+z[31]*tmp1[43]+z[43]*tmp1[44]+z[55]*tmp1[45]+z[67]*tmp1[46]+z[79]*tmp1[47]+z[91]*tmp1[48]+z[6]*tmp1[49]+z[18]*tmp1[50]
                    +z[30]*tmp1[51]+z[42]*tmp1[52]+z[54]*tmp1[53]+z[66]*tmp1[54]+z[78]*tmp1[55]+z[90]*tmp1[56]+z[5]*tmp1[57]+z[17]*tmp1[58]+z[29]*tmp1[59]+z[41]*tmp1[60]
                    +z[53]*tmp1[61]+z[65]*tmp1[62]+z[77]*tmp1[63]+z[89]*tmp1[64]+z[4]*tmp1[65]+z[16]*tmp1[66]+z[28]*tmp1[67]+z[40]*tmp1[68]+z[52]*tmp1[69]+z[64]*tmp1[70]
                    +z[76]*tmp1[71]+z[88]*tmp1[72]+z[3]*tmp1[73]+z[15]*tmp1[74]+z[27]*tmp1[75]+z[39]*tmp1[76]+z[51]*tmp1[77]+z[63]*tmp1[78]+z[75]*tmp1[79]+z[87]*tmp1[80]
                    +z[2]*tmp1[81]+z[14]*tmp1[82]+z[26]*tmp1[83]+z[38]*tmp1[84]+z[50]*tmp1[85]+z[62]*tmp1[86]+z[74]*tmp1[87]+z[86]*tmp1[88]+z[1]*tmp1[89]+z[13]*tmp1[90]
                    +z[25]*tmp1[91]+z[37]*tmp1[92]+z[49]*tmp1[93]+z[61]*tmp1[94]+z[73]*tmp1[95]+z[85]*tmp1[96];
                    tab[nb_tmp3+13*nb3]=tab[nb_tmp3+13*nb3]+z[0]*tmp1[0]
                    +z[13]*tmp1[1]+z[26]*tmp1[2]+z[39]*tmp1[3]+z[52]*tmp1[4]+z[65]*tmp1[5]+z[78]*tmp1[6]+z[91]*tmp1[7]+z[7]*tmp1[8]+z[20]*tmp1[9]+z[33]*tmp1[10]
                    +z[46]*tmp1[11]+z[59]*tmp1[12]+z[72]*tmp1[13]+z[85]*tmp1[14]+z[1]*tmp1[15]+z[14]*tmp1[16]+z[27]*tmp1[17]+z[40]*tmp1[18]+z[53]*tmp1[19]+z[66]*tmp1[20]
                    +z[79]*tmp1[21]+z[92]*tmp1[22]+z[8]*tmp1[23]+z[21]*tmp1[24]+z[34]*tmp1[25]+z[47]*tmp1[26]+z[60]*tmp1[27]+z[73]*tmp1[28]+z[86]*tmp1[29]+z[2]*tmp1[30]
                    +z[15]*tmp1[31]+z[28]*tmp1[32]+z[41]*tmp1[33]+z[54]*tmp1[34]+z[67]*tmp1[35]+z[80]*tmp1[36]+z[93]*tmp1[37]+z[9]*tmp1[38]+z[22]*tmp1[39]+z[35]*tmp1[40]
                    +z[48]*tmp1[41]+z[61]*tmp1[42]+z[74]*tmp1[43]+z[87]*tmp1[44]+z[3]*tmp1[45]+z[16]*tmp1[46]+z[29]*tmp1[47]+z[42]*tmp1[48]+z[55]*tmp1[49]+z[68]*tmp1[50]
                    +z[81]*tmp1[51]+z[94]*tmp1[52]+z[10]*tmp1[53]+z[23]*tmp1[54]+z[36]*tmp1[55]+z[49]*tmp1[56]+z[62]*tmp1[57]+z[75]*tmp1[58]+z[88]*tmp1[59]+z[4]*tmp1[60]
                    +z[17]*tmp1[61]+z[30]*tmp1[62]+z[43]*tmp1[63]+z[56]*tmp1[64]+z[69]*tmp1[65]+z[82]*tmp1[66]+z[95]*tmp1[67]+z[11]*tmp1[68]+z[24]*tmp1[69]+z[37]*tmp1[70]
                    +z[50]*tmp1[71]+z[63]*tmp1[72]+z[76]*tmp1[73]+z[89]*tmp1[74]+z[5]*tmp1[75]+z[18]*tmp1[76]+z[31]*tmp1[77]+z[44]*tmp1[78]+z[57]*tmp1[79]+z[70]*tmp1[80]
                    +z[83]*tmp1[81]+z[96]*tmp1[82]+z[12]*tmp1[83]+z[25]*tmp1[84]+z[38]*tmp1[85]+z[51]*tmp1[86]+z[64]*tmp1[87]+z[77]*tmp1[88]+z[90]*tmp1[89]+z[6]*tmp1[90]
                    +z[19]*tmp1[91]+z[32]*tmp1[92]+z[45]*tmp1[93]+z[58]*tmp1[94]+z[71]*tmp1[95]+z[84]*tmp1[96];
                    tab[nb_tmp3+14*nb3]=tab[nb_tmp3+14*nb3]+z[0]*tmp1[0]
                    +z[14]*tmp1[1]+z[28]*tmp1[2]+z[42]*tmp1[3]+z[56]*tmp1[4]+z[70]*tmp1[5]+z[84]*tmp1[6]+z[1]*tmp1[7]+z[15]*tmp1[8]+z[29]*tmp1[9]+z[43]*tmp1[10]
                    +z[57]*tmp1[11]+z[71]*tmp1[12]+z[85]*tmp1[13]+z[2]*tmp1[14]+z[16]*tmp1[15]+z[30]*tmp1[16]+z[44]*tmp1[17]+z[58]*tmp1[18]+z[72]*tmp1[19]+z[86]*tmp1[20]
                    +z[3]*tmp1[21]+z[17]*tmp1[22]+z[31]*tmp1[23]+z[45]*tmp1[24]+z[59]*tmp1[25]+z[73]*tmp1[26]+z[87]*tmp1[27]+z[4]*tmp1[28]+z[18]*tmp1[29]+z[32]*tmp1[30]
                    +z[46]*tmp1[31]+z[60]*tmp1[32]+z[74]*tmp1[33]+z[88]*tmp1[34]+z[5]*tmp1[35]+z[19]*tmp1[36]+z[33]*tmp1[37]+z[47]*tmp1[38]+z[61]*tmp1[39]+z[75]*tmp1[40]
                    +z[89]*tmp1[41]+z[6]*tmp1[42]+z[20]*tmp1[43]+z[34]*tmp1[44]+z[48]*tmp1[45]+z[62]*tmp1[46]+z[76]*tmp1[47]+z[90]*tmp1[48]+z[7]*tmp1[49]+z[21]*tmp1[50]
                    +z[35]*tmp1[51]+z[49]*tmp1[52]+z[63]*tmp1[53]+z[77]*tmp1[54]+z[91]*tmp1[55]+z[8]*tmp1[56]+z[22]*tmp1[57]+z[36]*tmp1[58]+z[50]*tmp1[59]+z[64]*tmp1[60]
                    +z[78]*tmp1[61]+z[92]*tmp1[62]+z[9]*tmp1[63]+z[23]*tmp1[64]+z[37]*tmp1[65]+z[51]*tmp1[66]+z[65]*tmp1[67]+z[79]*tmp1[68]+z[93]*tmp1[69]+z[10]*tmp1[70]
                    +z[24]*tmp1[71]+z[38]*tmp1[72]+z[52]*tmp1[73]+z[66]*tmp1[74]+z[80]*tmp1[75]+z[94]*tmp1[76]+z[11]*tmp1[77]+z[25]*tmp1[78]+z[39]*tmp1[79]+z[53]*tmp1[80]
                    +z[67]*tmp1[81]+z[81]*tmp1[82]+z[95]*tmp1[83]+z[12]*tmp1[84]+z[26]*tmp1[85]+z[40]*tmp1[86]+z[54]*tmp1[87]+z[68]*tmp1[88]+z[82]*tmp1[89]+z[96]*tmp1[90]
                    +z[13]*tmp1[91]+z[27]*tmp1[92]+z[41]*tmp1[93]+z[55]*tmp1[94]+z[69]*tmp1[95]+z[83]*tmp1[96];
                    tab[nb_tmp3+15*nb3]=tab[nb_tmp3+15*nb3]+z[0]*tmp1[0]
                    +z[15]*tmp1[1]+z[30]*tmp1[2]+z[45]*tmp1[3]+z[60]*tmp1[4]+z[75]*tmp1[5]+z[90]*tmp1[6]+z[8]*tmp1[7]+z[23]*tmp1[8]+z[38]*tmp1[9]+z[53]*tmp1[10]
                    +z[68]*tmp1[11]+z[83]*tmp1[12]+z[1]*tmp1[13]+z[16]*tmp1[14]+z[31]*tmp1[15]+z[46]*tmp1[16]+z[61]*tmp1[17]+z[76]*tmp1[18]+z[91]*tmp1[19]+z[9]*tmp1[20]
                    +z[24]*tmp1[21]+z[39]*tmp1[22]+z[54]*tmp1[23]+z[69]*tmp1[24]+z[84]*tmp1[25]+z[2]*tmp1[26]+z[17]*tmp1[27]+z[32]*tmp1[28]+z[47]*tmp1[29]+z[62]*tmp1[30]
                    +z[77]*tmp1[31]+z[92]*tmp1[32]+z[10]*tmp1[33]+z[25]*tmp1[34]+z[40]*tmp1[35]+z[55]*tmp1[36]+z[70]*tmp1[37]+z[85]*tmp1[38]+z[3]*tmp1[39]+z[18]*tmp1[40]
                    +z[33]*tmp1[41]+z[48]*tmp1[42]+z[63]*tmp1[43]+z[78]*tmp1[44]+z[93]*tmp1[45]+z[11]*tmp1[46]+z[26]*tmp1[47]+z[41]*tmp1[48]+z[56]*tmp1[49]+z[71]*tmp1[50]
                    +z[86]*tmp1[51]+z[4]*tmp1[52]+z[19]*tmp1[53]+z[34]*tmp1[54]+z[49]*tmp1[55]+z[64]*tmp1[56]+z[79]*tmp1[57]+z[94]*tmp1[58]+z[12]*tmp1[59]+z[27]*tmp1[60]
                    +z[42]*tmp1[61]+z[57]*tmp1[62]+z[72]*tmp1[63]+z[87]*tmp1[64]+z[5]*tmp1[65]+z[20]*tmp1[66]+z[35]*tmp1[67]+z[50]*tmp1[68]+z[65]*tmp1[69]+z[80]*tmp1[70]
                    +z[95]*tmp1[71]+z[13]*tmp1[72]+z[28]*tmp1[73]+z[43]*tmp1[74]+z[58]*tmp1[75]+z[73]*tmp1[76]+z[88]*tmp1[77]+z[6]*tmp1[78]+z[21]*tmp1[79]+z[36]*tmp1[80]
                    +z[51]*tmp1[81]+z[66]*tmp1[82]+z[81]*tmp1[83]+z[96]*tmp1[84]+z[14]*tmp1[85]+z[29]*tmp1[86]+z[44]*tmp1[87]+z[59]*tmp1[88]+z[74]*tmp1[89]+z[89]*tmp1[90]
                    +z[7]*tmp1[91]+z[22]*tmp1[92]+z[37]*tmp1[93]+z[52]*tmp1[94]+z[67]*tmp1[95]+z[82]*tmp1[96];
                    tab[nb_tmp3+16*nb3]=tab[nb_tmp3+16*nb3]+z[0]*tmp1[0]
                    +z[16]*tmp1[1]+z[32]*tmp1[2]+z[48]*tmp1[3]+z[64]*tmp1[4]+z[80]*tmp1[5]+z[96]*tmp1[6]+z[15]*tmp1[7]+z[31]*tmp1[8]+z[47]*tmp1[9]+z[63]*tmp1[10]
                    +z[79]*tmp1[11]+z[95]*tmp1[12]+z[14]*tmp1[13]+z[30]*tmp1[14]+z[46]*tmp1[15]+z[62]*tmp1[16]+z[78]*tmp1[17]+z[94]*tmp1[18]+z[13]*tmp1[19]+z[29]*tmp1[20]
                    +z[45]*tmp1[21]+z[61]*tmp1[22]+z[77]*tmp1[23]+z[93]*tmp1[24]+z[12]*tmp1[25]+z[28]*tmp1[26]+z[44]*tmp1[27]+z[60]*tmp1[28]+z[76]*tmp1[29]+z[92]*tmp1[30]
                    +z[11]*tmp1[31]+z[27]*tmp1[32]+z[43]*tmp1[33]+z[59]*tmp1[34]+z[75]*tmp1[35]+z[91]*tmp1[36]+z[10]*tmp1[37]+z[26]*tmp1[38]+z[42]*tmp1[39]+z[58]*tmp1[40]
                    +z[74]*tmp1[41]+z[90]*tmp1[42]+z[9]*tmp1[43]+z[25]*tmp1[44]+z[41]*tmp1[45]+z[57]*tmp1[46]+z[73]*tmp1[47]+z[89]*tmp1[48]+z[8]*tmp1[49]+z[24]*tmp1[50]
                    +z[40]*tmp1[51]+z[56]*tmp1[52]+z[72]*tmp1[53]+z[88]*tmp1[54]+z[7]*tmp1[55]+z[23]*tmp1[56]+z[39]*tmp1[57]+z[55]*tmp1[58]+z[71]*tmp1[59]+z[87]*tmp1[60]
                    +z[6]*tmp1[61]+z[22]*tmp1[62]+z[38]*tmp1[63]+z[54]*tmp1[64]+z[70]*tmp1[65]+z[86]*tmp1[66]+z[5]*tmp1[67]+z[21]*tmp1[68]+z[37]*tmp1[69]+z[53]*tmp1[70]
                    +z[69]*tmp1[71]+z[85]*tmp1[72]+z[4]*tmp1[73]+z[20]*tmp1[74]+z[36]*tmp1[75]+z[52]*tmp1[76]+z[68]*tmp1[77]+z[84]*tmp1[78]+z[3]*tmp1[79]+z[19]*tmp1[80]
                    +z[35]*tmp1[81]+z[51]*tmp1[82]+z[67]*tmp1[83]+z[83]*tmp1[84]+z[2]*tmp1[85]+z[18]*tmp1[86]+z[34]*tmp1[87]+z[50]*tmp1[88]+z[66]*tmp1[89]+z[82]*tmp1[90]
                    +z[1]*tmp1[91]+z[17]*tmp1[92]+z[33]*tmp1[93]+z[49]*tmp1[94]+z[65]*tmp1[95]+z[81]*tmp1[96];
                    tab[nb_tmp3+17*nb3]=tab[nb_tmp3+17*nb3]+z[0]*tmp1[0]
                    +z[17]*tmp1[1]+z[34]*tmp1[2]+z[51]*tmp1[3]+z[68]*tmp1[4]+z[85]*tmp1[5]+z[5]*tmp1[6]+z[22]*tmp1[7]+z[39]*tmp1[8]+z[56]*tmp1[9]+z[73]*tmp1[10]
                    +z[90]*tmp1[11]+z[10]*tmp1[12]+z[27]*tmp1[13]+z[44]*tmp1[14]+z[61]*tmp1[15]+z[78]*tmp1[16]+z[95]*tmp1[17]+z[15]*tmp1[18]+z[32]*tmp1[19]+z[49]*tmp1[20]
                    +z[66]*tmp1[21]+z[83]*tmp1[22]+z[3]*tmp1[23]+z[20]*tmp1[24]+z[37]*tmp1[25]+z[54]*tmp1[26]+z[71]*tmp1[27]+z[88]*tmp1[28]+z[8]*tmp1[29]+z[25]*tmp1[30]
                    +z[42]*tmp1[31]+z[59]*tmp1[32]+z[76]*tmp1[33]+z[93]*tmp1[34]+z[13]*tmp1[35]+z[30]*tmp1[36]+z[47]*tmp1[37]+z[64]*tmp1[38]+z[81]*tmp1[39]+z[1]*tmp1[40]
                    +z[18]*tmp1[41]+z[35]*tmp1[42]+z[52]*tmp1[43]+z[69]*tmp1[44]+z[86]*tmp1[45]+z[6]*tmp1[46]+z[23]*tmp1[47]+z[40]*tmp1[48]+z[57]*tmp1[49]+z[74]*tmp1[50]
                    +z[91]*tmp1[51]+z[11]*tmp1[52]+z[28]*tmp1[53]+z[45]*tmp1[54]+z[62]*tmp1[55]+z[79]*tmp1[56]+z[96]*tmp1[57]+z[16]*tmp1[58]+z[33]*tmp1[59]+z[50]*tmp1[60]
                    +z[67]*tmp1[61]+z[84]*tmp1[62]+z[4]*tmp1[63]+z[21]*tmp1[64]+z[38]*tmp1[65]+z[55]*tmp1[66]+z[72]*tmp1[67]+z[89]*tmp1[68]+z[9]*tmp1[69]+z[26]*tmp1[70]
                    +z[43]*tmp1[71]+z[60]*tmp1[72]+z[77]*tmp1[73]+z[94]*tmp1[74]+z[14]*tmp1[75]+z[31]*tmp1[76]+z[48]*tmp1[77]+z[65]*tmp1[78]+z[82]*tmp1[79]+z[2]*tmp1[80]
                    +z[19]*tmp1[81]+z[36]*tmp1[82]+z[53]*tmp1[83]+z[70]*tmp1[84]+z[87]*tmp1[85]+z[7]*tmp1[86]+z[24]*tmp1[87]+z[41]*tmp1[88]+z[58]*tmp1[89]+z[75]*tmp1[90]
                    +z[92]*tmp1[91]+z[12]*tmp1[92]+z[29]*tmp1[93]+z[46]*tmp1[94]+z[63]*tmp1[95]+z[80]*tmp1[96];
                    tab[nb_tmp3+18*nb3]=tab[nb_tmp3+18*nb3]+z[0]*tmp1[0]
                    +z[18]*tmp1[1]+z[36]*tmp1[2]+z[54]*tmp1[3]+z[72]*tmp1[4]+z[90]*tmp1[5]+z[11]*tmp1[6]+z[29]*tmp1[7]+z[47]*tmp1[8]+z[65]*tmp1[9]+z[83]*tmp1[10]
                    +z[4]*tmp1[11]+z[22]*tmp1[12]+z[40]*tmp1[13]+z[58]*tmp1[14]+z[76]*tmp1[15]+z[94]*tmp1[16]+z[15]*tmp1[17]+z[33]*tmp1[18]+z[51]*tmp1[19]+z[69]*tmp1[20]
                    +z[87]*tmp1[21]+z[8]*tmp1[22]+z[26]*tmp1[23]+z[44]*tmp1[24]+z[62]*tmp1[25]+z[80]*tmp1[26]+z[1]*tmp1[27]+z[19]*tmp1[28]+z[37]*tmp1[29]+z[55]*tmp1[30]
                    +z[73]*tmp1[31]+z[91]*tmp1[32]+z[12]*tmp1[33]+z[30]*tmp1[34]+z[48]*tmp1[35]+z[66]*tmp1[36]+z[84]*tmp1[37]+z[5]*tmp1[38]+z[23]*tmp1[39]+z[41]*tmp1[40]
                    +z[59]*tmp1[41]+z[77]*tmp1[42]+z[95]*tmp1[43]+z[16]*tmp1[44]+z[34]*tmp1[45]+z[52]*tmp1[46]+z[70]*tmp1[47]+z[88]*tmp1[48]+z[9]*tmp1[49]+z[27]*tmp1[50]
                    +z[45]*tmp1[51]+z[63]*tmp1[52]+z[81]*tmp1[53]+z[2]*tmp1[54]+z[20]*tmp1[55]+z[38]*tmp1[56]+z[56]*tmp1[57]+z[74]*tmp1[58]+z[92]*tmp1[59]+z[13]*tmp1[60]
                    +z[31]*tmp1[61]+z[49]*tmp1[62]+z[67]*tmp1[63]+z[85]*tmp1[64]+z[6]*tmp1[65]+z[24]*tmp1[66]+z[42]*tmp1[67]+z[60]*tmp1[68]+z[78]*tmp1[69]+z[96]*tmp1[70]
                    +z[17]*tmp1[71]+z[35]*tmp1[72]+z[53]*tmp1[73]+z[71]*tmp1[74]+z[89]*tmp1[75]+z[10]*tmp1[76]+z[28]*tmp1[77]+z[46]*tmp1[78]+z[64]*tmp1[79]+z[82]*tmp1[80]
                    +z[3]*tmp1[81]+z[21]*tmp1[82]+z[39]*tmp1[83]+z[57]*tmp1[84]+z[75]*tmp1[85]+z[93]*tmp1[86]+z[14]*tmp1[87]+z[32]*tmp1[88]+z[50]*tmp1[89]+z[68]*tmp1[90]
                    +z[86]*tmp1[91]+z[7]*tmp1[92]+z[25]*tmp1[93]+z[43]*tmp1[94]+z[61]*tmp1[95]+z[79]*tmp1[96];
                    tab[nb_tmp3+19*nb3]=tab[nb_tmp3+19*nb3]+z[0]*tmp1[0]
                    +z[19]*tmp1[1]+z[38]*tmp1[2]+z[57]*tmp1[3]+z[76]*tmp1[4]+z[95]*tmp1[5]+z[17]*tmp1[6]+z[36]*tmp1[7]+z[55]*tmp1[8]+z[74]*tmp1[9]+z[93]*tmp1[10]
                    +z[15]*tmp1[11]+z[34]*tmp1[12]+z[53]*tmp1[13]+z[72]*tmp1[14]+z[91]*tmp1[15]+z[13]*tmp1[16]+z[32]*tmp1[17]+z[51]*tmp1[18]+z[70]*tmp1[19]+z[89]*tmp1[20]
                    +z[11]*tmp1[21]+z[30]*tmp1[22]+z[49]*tmp1[23]+z[68]*tmp1[24]+z[87]*tmp1[25]+z[9]*tmp1[26]+z[28]*tmp1[27]+z[47]*tmp1[28]+z[66]*tmp1[29]+z[85]*tmp1[30]
                    +z[7]*tmp1[31]+z[26]*tmp1[32]+z[45]*tmp1[33]+z[64]*tmp1[34]+z[83]*tmp1[35]+z[5]*tmp1[36]+z[24]*tmp1[37]+z[43]*tmp1[38]+z[62]*tmp1[39]+z[81]*tmp1[40]
                    +z[3]*tmp1[41]+z[22]*tmp1[42]+z[41]*tmp1[43]+z[60]*tmp1[44]+z[79]*tmp1[45]+z[1]*tmp1[46]+z[20]*tmp1[47]+z[39]*tmp1[48]+z[58]*tmp1[49]+z[77]*tmp1[50]
                    +z[96]*tmp1[51]+z[18]*tmp1[52]+z[37]*tmp1[53]+z[56]*tmp1[54]+z[75]*tmp1[55]+z[94]*tmp1[56]+z[16]*tmp1[57]+z[35]*tmp1[58]+z[54]*tmp1[59]+z[73]*tmp1[60]
                    +z[92]*tmp1[61]+z[14]*tmp1[62]+z[33]*tmp1[63]+z[52]*tmp1[64]+z[71]*tmp1[65]+z[90]*tmp1[66]+z[12]*tmp1[67]+z[31]*tmp1[68]+z[50]*tmp1[69]+z[69]*tmp1[70]
                    +z[88]*tmp1[71]+z[10]*tmp1[72]+z[29]*tmp1[73]+z[48]*tmp1[74]+z[67]*tmp1[75]+z[86]*tmp1[76]+z[8]*tmp1[77]+z[27]*tmp1[78]+z[46]*tmp1[79]+z[65]*tmp1[80]
                    +z[84]*tmp1[81]+z[6]*tmp1[82]+z[25]*tmp1[83]+z[44]*tmp1[84]+z[63]*tmp1[85]+z[82]*tmp1[86]+z[4]*tmp1[87]+z[23]*tmp1[88]+z[42]*tmp1[89]+z[61]*tmp1[90]
                    +z[80]*tmp1[91]+z[2]*tmp1[92]+z[21]*tmp1[93]+z[40]*tmp1[94]+z[59]*tmp1[95]+z[78]*tmp1[96];
                    tab[nb_tmp3+20*nb3]=tab[nb_tmp3+20*nb3]+z[0]*tmp1[0]
                    +z[20]*tmp1[1]+z[40]*tmp1[2]+z[60]*tmp1[3]+z[80]*tmp1[4]+z[3]*tmp1[5]+z[23]*tmp1[6]+z[43]*tmp1[7]+z[63]*tmp1[8]+z[83]*tmp1[9]+z[6]*tmp1[10]
                    +z[26]*tmp1[11]+z[46]*tmp1[12]+z[66]*tmp1[13]+z[86]*tmp1[14]+z[9]*tmp1[15]+z[29]*tmp1[16]+z[49]*tmp1[17]+z[69]*tmp1[18]+z[89]*tmp1[19]+z[12]*tmp1[20]
                    +z[32]*tmp1[21]+z[52]*tmp1[22]+z[72]*tmp1[23]+z[92]*tmp1[24]+z[15]*tmp1[25]+z[35]*tmp1[26]+z[55]*tmp1[27]+z[75]*tmp1[28]+z[95]*tmp1[29]+z[18]*tmp1[30]
                    +z[38]*tmp1[31]+z[58]*tmp1[32]+z[78]*tmp1[33]+z[1]*tmp1[34]+z[21]*tmp1[35]+z[41]*tmp1[36]+z[61]*tmp1[37]+z[81]*tmp1[38]+z[4]*tmp1[39]+z[24]*tmp1[40]
                    +z[44]*tmp1[41]+z[64]*tmp1[42]+z[84]*tmp1[43]+z[7]*tmp1[44]+z[27]*tmp1[45]+z[47]*tmp1[46]+z[67]*tmp1[47]+z[87]*tmp1[48]+z[10]*tmp1[49]+z[30]*tmp1[50]
                    +z[50]*tmp1[51]+z[70]*tmp1[52]+z[90]*tmp1[53]+z[13]*tmp1[54]+z[33]*tmp1[55]+z[53]*tmp1[56]+z[73]*tmp1[57]+z[93]*tmp1[58]+z[16]*tmp1[59]+z[36]*tmp1[60]
                    +z[56]*tmp1[61]+z[76]*tmp1[62]+z[96]*tmp1[63]+z[19]*tmp1[64]+z[39]*tmp1[65]+z[59]*tmp1[66]+z[79]*tmp1[67]+z[2]*tmp1[68]+z[22]*tmp1[69]+z[42]*tmp1[70]
                    +z[62]*tmp1[71]+z[82]*tmp1[72]+z[5]*tmp1[73]+z[25]*tmp1[74]+z[45]*tmp1[75]+z[65]*tmp1[76]+z[85]*tmp1[77]+z[8]*tmp1[78]+z[28]*tmp1[79]+z[48]*tmp1[80]
                    +z[68]*tmp1[81]+z[88]*tmp1[82]+z[11]*tmp1[83]+z[31]*tmp1[84]+z[51]*tmp1[85]+z[71]*tmp1[86]+z[91]*tmp1[87]+z[14]*tmp1[88]+z[34]*tmp1[89]+z[54]*tmp1[90]
                    +z[74]*tmp1[91]+z[94]*tmp1[92]+z[17]*tmp1[93]+z[37]*tmp1[94]+z[57]*tmp1[95]+z[77]*tmp1[96];
                    tab[nb_tmp3+21*nb3]=tab[nb_tmp3+21*nb3]+z[0]*tmp1[0]
                    +z[21]*tmp1[1]+z[42]*tmp1[2]+z[63]*tmp1[3]+z[84]*tmp1[4]+z[8]*tmp1[5]+z[29]*tmp1[6]+z[50]*tmp1[7]+z[71]*tmp1[8]+z[92]*tmp1[9]+z[16]*tmp1[10]
                    +z[37]*tmp1[11]+z[58]*tmp1[12]+z[79]*tmp1[13]+z[3]*tmp1[14]+z[24]*tmp1[15]+z[45]*tmp1[16]+z[66]*tmp1[17]+z[87]*tmp1[18]+z[11]*tmp1[19]+z[32]*tmp1[20]
                    +z[53]*tmp1[21]+z[74]*tmp1[22]+z[95]*tmp1[23]+z[19]*tmp1[24]+z[40]*tmp1[25]+z[61]*tmp1[26]+z[82]*tmp1[27]+z[6]*tmp1[28]+z[27]*tmp1[29]+z[48]*tmp1[30]
                    +z[69]*tmp1[31]+z[90]*tmp1[32]+z[14]*tmp1[33]+z[35]*tmp1[34]+z[56]*tmp1[35]+z[77]*tmp1[36]+z[1]*tmp1[37]+z[22]*tmp1[38]+z[43]*tmp1[39]+z[64]*tmp1[40]
                    +z[85]*tmp1[41]+z[9]*tmp1[42]+z[30]*tmp1[43]+z[51]*tmp1[44]+z[72]*tmp1[45]+z[93]*tmp1[46]+z[17]*tmp1[47]+z[38]*tmp1[48]+z[59]*tmp1[49]+z[80]*tmp1[50]
                    +z[4]*tmp1[51]+z[25]*tmp1[52]+z[46]*tmp1[53]+z[67]*tmp1[54]+z[88]*tmp1[55]+z[12]*tmp1[56]+z[33]*tmp1[57]+z[54]*tmp1[58]+z[75]*tmp1[59]+z[96]*tmp1[60]
                    +z[20]*tmp1[61]+z[41]*tmp1[62]+z[62]*tmp1[63]+z[83]*tmp1[64]+z[7]*tmp1[65]+z[28]*tmp1[66]+z[49]*tmp1[67]+z[70]*tmp1[68]+z[91]*tmp1[69]+z[15]*tmp1[70]
                    +z[36]*tmp1[71]+z[57]*tmp1[72]+z[78]*tmp1[73]+z[2]*tmp1[74]+z[23]*tmp1[75]+z[44]*tmp1[76]+z[65]*tmp1[77]+z[86]*tmp1[78]+z[10]*tmp1[79]+z[31]*tmp1[80]
                    +z[52]*tmp1[81]+z[73]*tmp1[82]+z[94]*tmp1[83]+z[18]*tmp1[84]+z[39]*tmp1[85]+z[60]*tmp1[86]+z[81]*tmp1[87]+z[5]*tmp1[88]+z[26]*tmp1[89]+z[47]*tmp1[90]
                    +z[68]*tmp1[91]+z[89]*tmp1[92]+z[13]*tmp1[93]+z[34]*tmp1[94]+z[55]*tmp1[95]+z[76]*tmp1[96];
                    tab[nb_tmp3+22*nb3]=tab[nb_tmp3+22*nb3]+z[0]*tmp1[0]
                    +z[22]*tmp1[1]+z[44]*tmp1[2]+z[66]*tmp1[3]+z[88]*tmp1[4]+z[13]*tmp1[5]+z[35]*tmp1[6]+z[57]*tmp1[7]+z[79]*tmp1[8]+z[4]*tmp1[9]+z[26]*tmp1[10]
                    +z[48]*tmp1[11]+z[70]*tmp1[12]+z[92]*tmp1[13]+z[17]*tmp1[14]+z[39]*tmp1[15]+z[61]*tmp1[16]+z[83]*tmp1[17]+z[8]*tmp1[18]+z[30]*tmp1[19]+z[52]*tmp1[20]
                    +z[74]*tmp1[21]+z[96]*tmp1[22]+z[21]*tmp1[23]+z[43]*tmp1[24]+z[65]*tmp1[25]+z[87]*tmp1[26]+z[12]*tmp1[27]+z[34]*tmp1[28]+z[56]*tmp1[29]+z[78]*tmp1[30]
                    +z[3]*tmp1[31]+z[25]*tmp1[32]+z[47]*tmp1[33]+z[69]*tmp1[34]+z[91]*tmp1[35]+z[16]*tmp1[36]+z[38]*tmp1[37]+z[60]*tmp1[38]+z[82]*tmp1[39]+z[7]*tmp1[40]
                    +z[29]*tmp1[41]+z[51]*tmp1[42]+z[73]*tmp1[43]+z[95]*tmp1[44]+z[20]*tmp1[45]+z[42]*tmp1[46]+z[64]*tmp1[47]+z[86]*tmp1[48]+z[11]*tmp1[49]+z[33]*tmp1[50]
                    +z[55]*tmp1[51]+z[77]*tmp1[52]+z[2]*tmp1[53]+z[24]*tmp1[54]+z[46]*tmp1[55]+z[68]*tmp1[56]+z[90]*tmp1[57]+z[15]*tmp1[58]+z[37]*tmp1[59]+z[59]*tmp1[60]
                    +z[81]*tmp1[61]+z[6]*tmp1[62]+z[28]*tmp1[63]+z[50]*tmp1[64]+z[72]*tmp1[65]+z[94]*tmp1[66]+z[19]*tmp1[67]+z[41]*tmp1[68]+z[63]*tmp1[69]+z[85]*tmp1[70]
                    +z[10]*tmp1[71]+z[32]*tmp1[72]+z[54]*tmp1[73]+z[76]*tmp1[74]+z[1]*tmp1[75]+z[23]*tmp1[76]+z[45]*tmp1[77]+z[67]*tmp1[78]+z[89]*tmp1[79]+z[14]*tmp1[80]
                    +z[36]*tmp1[81]+z[58]*tmp1[82]+z[80]*tmp1[83]+z[5]*tmp1[84]+z[27]*tmp1[85]+z[49]*tmp1[86]+z[71]*tmp1[87]+z[93]*tmp1[88]+z[18]*tmp1[89]+z[40]*tmp1[90]
                    +z[62]*tmp1[91]+z[84]*tmp1[92]+z[9]*tmp1[93]+z[31]*tmp1[94]+z[53]*tmp1[95]+z[75]*tmp1[96];
                    tab[nb_tmp3+23*nb3]=tab[nb_tmp3+23*nb3]+z[0]*tmp1[0]
                    +z[23]*tmp1[1]+z[46]*tmp1[2]+z[69]*tmp1[3]+z[92]*tmp1[4]+z[18]*tmp1[5]+z[41]*tmp1[6]+z[64]*tmp1[7]+z[87]*tmp1[8]+z[13]*tmp1[9]+z[36]*tmp1[10]
                    +z[59]*tmp1[11]+z[82]*tmp1[12]+z[8]*tmp1[13]+z[31]*tmp1[14]+z[54]*tmp1[15]+z[77]*tmp1[16]+z[3]*tmp1[17]+z[26]*tmp1[18]+z[49]*tmp1[19]+z[72]*tmp1[20]
                    +z[95]*tmp1[21]+z[21]*tmp1[22]+z[44]*tmp1[23]+z[67]*tmp1[24]+z[90]*tmp1[25]+z[16]*tmp1[26]+z[39]*tmp1[27]+z[62]*tmp1[28]+z[85]*tmp1[29]+z[11]*tmp1[30]
                    +z[34]*tmp1[31]+z[57]*tmp1[32]+z[80]*tmp1[33]+z[6]*tmp1[34]+z[29]*tmp1[35]+z[52]*tmp1[36]+z[75]*tmp1[37]+z[1]*tmp1[38]+z[24]*tmp1[39]+z[47]*tmp1[40]
                    +z[70]*tmp1[41]+z[93]*tmp1[42]+z[19]*tmp1[43]+z[42]*tmp1[44]+z[65]*tmp1[45]+z[88]*tmp1[46]+z[14]*tmp1[47]+z[37]*tmp1[48]+z[60]*tmp1[49]+z[83]*tmp1[50]
                    +z[9]*tmp1[51]+z[32]*tmp1[52]+z[55]*tmp1[53]+z[78]*tmp1[54]+z[4]*tmp1[55]+z[27]*tmp1[56]+z[50]*tmp1[57]+z[73]*tmp1[58]+z[96]*tmp1[59]+z[22]*tmp1[60]
                    +z[45]*tmp1[61]+z[68]*tmp1[62]+z[91]*tmp1[63]+z[17]*tmp1[64]+z[40]*tmp1[65]+z[63]*tmp1[66]+z[86]*tmp1[67]+z[12]*tmp1[68]+z[35]*tmp1[69]+z[58]*tmp1[70]
                    +z[81]*tmp1[71]+z[7]*tmp1[72]+z[30]*tmp1[73]+z[53]*tmp1[74]+z[76]*tmp1[75]+z[2]*tmp1[76]+z[25]*tmp1[77]+z[48]*tmp1[78]+z[71]*tmp1[79]+z[94]*tmp1[80]
                    +z[20]*tmp1[81]+z[43]*tmp1[82]+z[66]*tmp1[83]+z[89]*tmp1[84]+z[15]*tmp1[85]+z[38]*tmp1[86]+z[61]*tmp1[87]+z[84]*tmp1[88]+z[10]*tmp1[89]+z[33]*tmp1[90]
                    +z[56]*tmp1[91]+z[79]*tmp1[92]+z[5]*tmp1[93]+z[28]*tmp1[94]+z[51]*tmp1[95]+z[74]*tmp1[96];
                    tab[nb_tmp3+24*nb3]=tab[nb_tmp3+24*nb3]+z[0]*tmp1[0]
                    +z[24]*tmp1[1]+z[48]*tmp1[2]+z[72]*tmp1[3]+z[96]*tmp1[4]+z[23]*tmp1[5]+z[47]*tmp1[6]+z[71]*tmp1[7]+z[95]*tmp1[8]+z[22]*tmp1[9]+z[46]*tmp1[10]
                    +z[70]*tmp1[11]+z[94]*tmp1[12]+z[21]*tmp1[13]+z[45]*tmp1[14]+z[69]*tmp1[15]+z[93]*tmp1[16]+z[20]*tmp1[17]+z[44]*tmp1[18]+z[68]*tmp1[19]+z[92]*tmp1[20]
                    +z[19]*tmp1[21]+z[43]*tmp1[22]+z[67]*tmp1[23]+z[91]*tmp1[24]+z[18]*tmp1[25]+z[42]*tmp1[26]+z[66]*tmp1[27]+z[90]*tmp1[28]+z[17]*tmp1[29]+z[41]*tmp1[30]
                    +z[65]*tmp1[31]+z[89]*tmp1[32]+z[16]*tmp1[33]+z[40]*tmp1[34]+z[64]*tmp1[35]+z[88]*tmp1[36]+z[15]*tmp1[37]+z[39]*tmp1[38]+z[63]*tmp1[39]+z[87]*tmp1[40]
                    +z[14]*tmp1[41]+z[38]*tmp1[42]+z[62]*tmp1[43]+z[86]*tmp1[44]+z[13]*tmp1[45]+z[37]*tmp1[46]+z[61]*tmp1[47]+z[85]*tmp1[48]+z[12]*tmp1[49]+z[36]*tmp1[50]
                    +z[60]*tmp1[51]+z[84]*tmp1[52]+z[11]*tmp1[53]+z[35]*tmp1[54]+z[59]*tmp1[55]+z[83]*tmp1[56]+z[10]*tmp1[57]+z[34]*tmp1[58]+z[58]*tmp1[59]+z[82]*tmp1[60]
                    +z[9]*tmp1[61]+z[33]*tmp1[62]+z[57]*tmp1[63]+z[81]*tmp1[64]+z[8]*tmp1[65]+z[32]*tmp1[66]+z[56]*tmp1[67]+z[80]*tmp1[68]+z[7]*tmp1[69]+z[31]*tmp1[70]
                    +z[55]*tmp1[71]+z[79]*tmp1[72]+z[6]*tmp1[73]+z[30]*tmp1[74]+z[54]*tmp1[75]+z[78]*tmp1[76]+z[5]*tmp1[77]+z[29]*tmp1[78]+z[53]*tmp1[79]+z[77]*tmp1[80]
                    +z[4]*tmp1[81]+z[28]*tmp1[82]+z[52]*tmp1[83]+z[76]*tmp1[84]+z[3]*tmp1[85]+z[27]*tmp1[86]+z[51]*tmp1[87]+z[75]*tmp1[88]+z[2]*tmp1[89]+z[26]*tmp1[90]
                    +z[50]*tmp1[91]+z[74]*tmp1[92]+z[1]*tmp1[93]+z[25]*tmp1[94]+z[49]*tmp1[95]+z[73]*tmp1[96];
                    tab[nb_tmp3+25*nb3]=tab[nb_tmp3+25*nb3]+z[0]*tmp1[0]
                    +z[25]*tmp1[1]+z[50]*tmp1[2]+z[75]*tmp1[3]+z[3]*tmp1[4]+z[28]*tmp1[5]+z[53]*tmp1[6]+z[78]*tmp1[7]+z[6]*tmp1[8]+z[31]*tmp1[9]+z[56]*tmp1[10]
                    +z[81]*tmp1[11]+z[9]*tmp1[12]+z[34]*tmp1[13]+z[59]*tmp1[14]+z[84]*tmp1[15]+z[12]*tmp1[16]+z[37]*tmp1[17]+z[62]*tmp1[18]+z[87]*tmp1[19]+z[15]*tmp1[20]
                    +z[40]*tmp1[21]+z[65]*tmp1[22]+z[90]*tmp1[23]+z[18]*tmp1[24]+z[43]*tmp1[25]+z[68]*tmp1[26]+z[93]*tmp1[27]+z[21]*tmp1[28]+z[46]*tmp1[29]+z[71]*tmp1[30]
                    +z[96]*tmp1[31]+z[24]*tmp1[32]+z[49]*tmp1[33]+z[74]*tmp1[34]+z[2]*tmp1[35]+z[27]*tmp1[36]+z[52]*tmp1[37]+z[77]*tmp1[38]+z[5]*tmp1[39]+z[30]*tmp1[40]
                    +z[55]*tmp1[41]+z[80]*tmp1[42]+z[8]*tmp1[43]+z[33]*tmp1[44]+z[58]*tmp1[45]+z[83]*tmp1[46]+z[11]*tmp1[47]+z[36]*tmp1[48]+z[61]*tmp1[49]+z[86]*tmp1[50]
                    +z[14]*tmp1[51]+z[39]*tmp1[52]+z[64]*tmp1[53]+z[89]*tmp1[54]+z[17]*tmp1[55]+z[42]*tmp1[56]+z[67]*tmp1[57]+z[92]*tmp1[58]+z[20]*tmp1[59]+z[45]*tmp1[60]
                    +z[70]*tmp1[61]+z[95]*tmp1[62]+z[23]*tmp1[63]+z[48]*tmp1[64]+z[73]*tmp1[65]+z[1]*tmp1[66]+z[26]*tmp1[67]+z[51]*tmp1[68]+z[76]*tmp1[69]+z[4]*tmp1[70]
                    +z[29]*tmp1[71]+z[54]*tmp1[72]+z[79]*tmp1[73]+z[7]*tmp1[74]+z[32]*tmp1[75]+z[57]*tmp1[76]+z[82]*tmp1[77]+z[10]*tmp1[78]+z[35]*tmp1[79]+z[60]*tmp1[80]
                    +z[85]*tmp1[81]+z[13]*tmp1[82]+z[38]*tmp1[83]+z[63]*tmp1[84]+z[88]*tmp1[85]+z[16]*tmp1[86]+z[41]*tmp1[87]+z[66]*tmp1[88]+z[91]*tmp1[89]+z[19]*tmp1[90]
                    +z[44]*tmp1[91]+z[69]*tmp1[92]+z[94]*tmp1[93]+z[22]*tmp1[94]+z[47]*tmp1[95]+z[72]*tmp1[96];
                    tab[nb_tmp3+26*nb3]=tab[nb_tmp3+26*nb3]+z[0]*tmp1[0]
                    +z[26]*tmp1[1]+z[52]*tmp1[2]+z[78]*tmp1[3]+z[7]*tmp1[4]+z[33]*tmp1[5]+z[59]*tmp1[6]+z[85]*tmp1[7]+z[14]*tmp1[8]+z[40]*tmp1[9]+z[66]*tmp1[10]
                    +z[92]*tmp1[11]+z[21]*tmp1[12]+z[47]*tmp1[13]+z[73]*tmp1[14]+z[2]*tmp1[15]+z[28]*tmp1[16]+z[54]*tmp1[17]+z[80]*tmp1[18]+z[9]*tmp1[19]+z[35]*tmp1[20]
                    +z[61]*tmp1[21]+z[87]*tmp1[22]+z[16]*tmp1[23]+z[42]*tmp1[24]+z[68]*tmp1[25]+z[94]*tmp1[26]+z[23]*tmp1[27]+z[49]*tmp1[28]+z[75]*tmp1[29]+z[4]*tmp1[30]
                    +z[30]*tmp1[31]+z[56]*tmp1[32]+z[82]*tmp1[33]+z[11]*tmp1[34]+z[37]*tmp1[35]+z[63]*tmp1[36]+z[89]*tmp1[37]+z[18]*tmp1[38]+z[44]*tmp1[39]+z[70]*tmp1[40]
                    +z[96]*tmp1[41]+z[25]*tmp1[42]+z[51]*tmp1[43]+z[77]*tmp1[44]+z[6]*tmp1[45]+z[32]*tmp1[46]+z[58]*tmp1[47]+z[84]*tmp1[48]+z[13]*tmp1[49]+z[39]*tmp1[50]
                    +z[65]*tmp1[51]+z[91]*tmp1[52]+z[20]*tmp1[53]+z[46]*tmp1[54]+z[72]*tmp1[55]+z[1]*tmp1[56]+z[27]*tmp1[57]+z[53]*tmp1[58]+z[79]*tmp1[59]+z[8]*tmp1[60]
                    +z[34]*tmp1[61]+z[60]*tmp1[62]+z[86]*tmp1[63]+z[15]*tmp1[64]+z[41]*tmp1[65]+z[67]*tmp1[66]+z[93]*tmp1[67]+z[22]*tmp1[68]+z[48]*tmp1[69]+z[74]*tmp1[70]
                    +z[3]*tmp1[71]+z[29]*tmp1[72]+z[55]*tmp1[73]+z[81]*tmp1[74]+z[10]*tmp1[75]+z[36]*tmp1[76]+z[62]*tmp1[77]+z[88]*tmp1[78]+z[17]*tmp1[79]+z[43]*tmp1[80]
                    +z[69]*tmp1[81]+z[95]*tmp1[82]+z[24]*tmp1[83]+z[50]*tmp1[84]+z[76]*tmp1[85]+z[5]*tmp1[86]+z[31]*tmp1[87]+z[57]*tmp1[88]+z[83]*tmp1[89]+z[12]*tmp1[90]
                    +z[38]*tmp1[91]+z[64]*tmp1[92]+z[90]*tmp1[93]+z[19]*tmp1[94]+z[45]*tmp1[95]+z[71]*tmp1[96];
                    tab[nb_tmp3+27*nb3]=tab[nb_tmp3+27*nb3]+z[0]*tmp1[0]
                    +z[27]*tmp1[1]+z[54]*tmp1[2]+z[81]*tmp1[3]+z[11]*tmp1[4]+z[38]*tmp1[5]+z[65]*tmp1[6]+z[92]*tmp1[7]+z[22]*tmp1[8]+z[49]*tmp1[9]+z[76]*tmp1[10]
                    +z[6]*tmp1[11]+z[33]*tmp1[12]+z[60]*tmp1[13]+z[87]*tmp1[14]+z[17]*tmp1[15]+z[44]*tmp1[16]+z[71]*tmp1[17]+z[1]*tmp1[18]+z[28]*tmp1[19]+z[55]*tmp1[20]
                    +z[82]*tmp1[21]+z[12]*tmp1[22]+z[39]*tmp1[23]+z[66]*tmp1[24]+z[93]*tmp1[25]+z[23]*tmp1[26]+z[50]*tmp1[27]+z[77]*tmp1[28]+z[7]*tmp1[29]+z[34]*tmp1[30]
                    +z[61]*tmp1[31]+z[88]*tmp1[32]+z[18]*tmp1[33]+z[45]*tmp1[34]+z[72]*tmp1[35]+z[2]*tmp1[36]+z[29]*tmp1[37]+z[56]*tmp1[38]+z[83]*tmp1[39]+z[13]*tmp1[40]
                    +z[40]*tmp1[41]+z[67]*tmp1[42]+z[94]*tmp1[43]+z[24]*tmp1[44]+z[51]*tmp1[45]+z[78]*tmp1[46]+z[8]*tmp1[47]+z[35]*tmp1[48]+z[62]*tmp1[49]+z[89]*tmp1[50]
                    +z[19]*tmp1[51]+z[46]*tmp1[52]+z[73]*tmp1[53]+z[3]*tmp1[54]+z[30]*tmp1[55]+z[57]*tmp1[56]+z[84]*tmp1[57]+z[14]*tmp1[58]+z[41]*tmp1[59]+z[68]*tmp1[60]
                    +z[95]*tmp1[61]+z[25]*tmp1[62]+z[52]*tmp1[63]+z[79]*tmp1[64]+z[9]*tmp1[65]+z[36]*tmp1[66]+z[63]*tmp1[67]+z[90]*tmp1[68]+z[20]*tmp1[69]+z[47]*tmp1[70]
                    +z[74]*tmp1[71]+z[4]*tmp1[72]+z[31]*tmp1[73]+z[58]*tmp1[74]+z[85]*tmp1[75]+z[15]*tmp1[76]+z[42]*tmp1[77]+z[69]*tmp1[78]+z[96]*tmp1[79]+z[26]*tmp1[80]
                    +z[53]*tmp1[81]+z[80]*tmp1[82]+z[10]*tmp1[83]+z[37]*tmp1[84]+z[64]*tmp1[85]+z[91]*tmp1[86]+z[21]*tmp1[87]+z[48]*tmp1[88]+z[75]*tmp1[89]+z[5]*tmp1[90]
                    +z[32]*tmp1[91]+z[59]*tmp1[92]+z[86]*tmp1[93]+z[16]*tmp1[94]+z[43]*tmp1[95]+z[70]*tmp1[96];
                    tab[nb_tmp3+28*nb3]=tab[nb_tmp3+28*nb3]+z[0]*tmp1[0]
                    +z[28]*tmp1[1]+z[56]*tmp1[2]+z[84]*tmp1[3]+z[15]*tmp1[4]+z[43]*tmp1[5]+z[71]*tmp1[6]+z[2]*tmp1[7]+z[30]*tmp1[8]+z[58]*tmp1[9]+z[86]*tmp1[10]
                    +z[17]*tmp1[11]+z[45]*tmp1[12]+z[73]*tmp1[13]+z[4]*tmp1[14]+z[32]*tmp1[15]+z[60]*tmp1[16]+z[88]*tmp1[17]+z[19]*tmp1[18]+z[47]*tmp1[19]+z[75]*tmp1[20]
                    +z[6]*tmp1[21]+z[34]*tmp1[22]+z[62]*tmp1[23]+z[90]*tmp1[24]+z[21]*tmp1[25]+z[49]*tmp1[26]+z[77]*tmp1[27]+z[8]*tmp1[28]+z[36]*tmp1[29]+z[64]*tmp1[30]
                    +z[92]*tmp1[31]+z[23]*tmp1[32]+z[51]*tmp1[33]+z[79]*tmp1[34]+z[10]*tmp1[35]+z[38]*tmp1[36]+z[66]*tmp1[37]+z[94]*tmp1[38]+z[25]*tmp1[39]+z[53]*tmp1[40]
                    +z[81]*tmp1[41]+z[12]*tmp1[42]+z[40]*tmp1[43]+z[68]*tmp1[44]+z[96]*tmp1[45]+z[27]*tmp1[46]+z[55]*tmp1[47]+z[83]*tmp1[48]+z[14]*tmp1[49]+z[42]*tmp1[50]
                    +z[70]*tmp1[51]+z[1]*tmp1[52]+z[29]*tmp1[53]+z[57]*tmp1[54]+z[85]*tmp1[55]+z[16]*tmp1[56]+z[44]*tmp1[57]+z[72]*tmp1[58]+z[3]*tmp1[59]+z[31]*tmp1[60]
                    +z[59]*tmp1[61]+z[87]*tmp1[62]+z[18]*tmp1[63]+z[46]*tmp1[64]+z[74]*tmp1[65]+z[5]*tmp1[66]+z[33]*tmp1[67]+z[61]*tmp1[68]+z[89]*tmp1[69]+z[20]*tmp1[70]
                    +z[48]*tmp1[71]+z[76]*tmp1[72]+z[7]*tmp1[73]+z[35]*tmp1[74]+z[63]*tmp1[75]+z[91]*tmp1[76]+z[22]*tmp1[77]+z[50]*tmp1[78]+z[78]*tmp1[79]+z[9]*tmp1[80]
                    +z[37]*tmp1[81]+z[65]*tmp1[82]+z[93]*tmp1[83]+z[24]*tmp1[84]+z[52]*tmp1[85]+z[80]*tmp1[86]+z[11]*tmp1[87]+z[39]*tmp1[88]+z[67]*tmp1[89]+z[95]*tmp1[90]
                    +z[26]*tmp1[91]+z[54]*tmp1[92]+z[82]*tmp1[93]+z[13]*tmp1[94]+z[41]*tmp1[95]+z[69]*tmp1[96];
                    tab[nb_tmp3+29*nb3]=tab[nb_tmp3+29*nb3]+z[0]*tmp1[0]
                    +z[29]*tmp1[1]+z[58]*tmp1[2]+z[87]*tmp1[3]+z[19]*tmp1[4]+z[48]*tmp1[5]+z[77]*tmp1[6]+z[9]*tmp1[7]+z[38]*tmp1[8]+z[67]*tmp1[9]+z[96]*tmp1[10]
                    +z[28]*tmp1[11]+z[57]*tmp1[12]+z[86]*tmp1[13]+z[18]*tmp1[14]+z[47]*tmp1[15]+z[76]*tmp1[16]+z[8]*tmp1[17]+z[37]*tmp1[18]+z[66]*tmp1[19]+z[95]*tmp1[20]
                    +z[27]*tmp1[21]+z[56]*tmp1[22]+z[85]*tmp1[23]+z[17]*tmp1[24]+z[46]*tmp1[25]+z[75]*tmp1[26]+z[7]*tmp1[27]+z[36]*tmp1[28]+z[65]*tmp1[29]+z[94]*tmp1[30]
                    +z[26]*tmp1[31]+z[55]*tmp1[32]+z[84]*tmp1[33]+z[16]*tmp1[34]+z[45]*tmp1[35]+z[74]*tmp1[36]+z[6]*tmp1[37]+z[35]*tmp1[38]+z[64]*tmp1[39]+z[93]*tmp1[40]
                    +z[25]*tmp1[41]+z[54]*tmp1[42]+z[83]*tmp1[43]+z[15]*tmp1[44]+z[44]*tmp1[45]+z[73]*tmp1[46]+z[5]*tmp1[47]+z[34]*tmp1[48]+z[63]*tmp1[49]+z[92]*tmp1[50]
                    +z[24]*tmp1[51]+z[53]*tmp1[52]+z[82]*tmp1[53]+z[14]*tmp1[54]+z[43]*tmp1[55]+z[72]*tmp1[56]+z[4]*tmp1[57]+z[33]*tmp1[58]+z[62]*tmp1[59]+z[91]*tmp1[60]
                    +z[23]*tmp1[61]+z[52]*tmp1[62]+z[81]*tmp1[63]+z[13]*tmp1[64]+z[42]*tmp1[65]+z[71]*tmp1[66]+z[3]*tmp1[67]+z[32]*tmp1[68]+z[61]*tmp1[69]+z[90]*tmp1[70]
                    +z[22]*tmp1[71]+z[51]*tmp1[72]+z[80]*tmp1[73]+z[12]*tmp1[74]+z[41]*tmp1[75]+z[70]*tmp1[76]+z[2]*tmp1[77]+z[31]*tmp1[78]+z[60]*tmp1[79]+z[89]*tmp1[80]
                    +z[21]*tmp1[81]+z[50]*tmp1[82]+z[79]*tmp1[83]+z[11]*tmp1[84]+z[40]*tmp1[85]+z[69]*tmp1[86]+z[1]*tmp1[87]+z[30]*tmp1[88]+z[59]*tmp1[89]+z[88]*tmp1[90]
                    +z[20]*tmp1[91]+z[49]*tmp1[92]+z[78]*tmp1[93]+z[10]*tmp1[94]+z[39]*tmp1[95]+z[68]*tmp1[96];
                    tab[nb_tmp3+30*nb3]=tab[nb_tmp3+30*nb3]+z[0]*tmp1[0]
                    +z[30]*tmp1[1]+z[60]*tmp1[2]+z[90]*tmp1[3]+z[23]*tmp1[4]+z[53]*tmp1[5]+z[83]*tmp1[6]+z[16]*tmp1[7]+z[46]*tmp1[8]+z[76]*tmp1[9]+z[9]*tmp1[10]
                    +z[39]*tmp1[11]+z[69]*tmp1[12]+z[2]*tmp1[13]+z[32]*tmp1[14]+z[62]*tmp1[15]+z[92]*tmp1[16]+z[25]*tmp1[17]+z[55]*tmp1[18]+z[85]*tmp1[19]+z[18]*tmp1[20]
                    +z[48]*tmp1[21]+z[78]*tmp1[22]+z[11]*tmp1[23]+z[41]*tmp1[24]+z[71]*tmp1[25]+z[4]*tmp1[26]+z[34]*tmp1[27]+z[64]*tmp1[28]+z[94]*tmp1[29]+z[27]*tmp1[30]
                    +z[57]*tmp1[31]+z[87]*tmp1[32]+z[20]*tmp1[33]+z[50]*tmp1[34]+z[80]*tmp1[35]+z[13]*tmp1[36]+z[43]*tmp1[37]+z[73]*tmp1[38]+z[6]*tmp1[39]+z[36]*tmp1[40]
                    +z[66]*tmp1[41]+z[96]*tmp1[42]+z[29]*tmp1[43]+z[59]*tmp1[44]+z[89]*tmp1[45]+z[22]*tmp1[46]+z[52]*tmp1[47]+z[82]*tmp1[48]+z[15]*tmp1[49]+z[45]*tmp1[50]
                    +z[75]*tmp1[51]+z[8]*tmp1[52]+z[38]*tmp1[53]+z[68]*tmp1[54]+z[1]*tmp1[55]+z[31]*tmp1[56]+z[61]*tmp1[57]+z[91]*tmp1[58]+z[24]*tmp1[59]+z[54]*tmp1[60]
                    +z[84]*tmp1[61]+z[17]*tmp1[62]+z[47]*tmp1[63]+z[77]*tmp1[64]+z[10]*tmp1[65]+z[40]*tmp1[66]+z[70]*tmp1[67]+z[3]*tmp1[68]+z[33]*tmp1[69]+z[63]*tmp1[70]
                    +z[93]*tmp1[71]+z[26]*tmp1[72]+z[56]*tmp1[73]+z[86]*tmp1[74]+z[19]*tmp1[75]+z[49]*tmp1[76]+z[79]*tmp1[77]+z[12]*tmp1[78]+z[42]*tmp1[79]+z[72]*tmp1[80]
                    +z[5]*tmp1[81]+z[35]*tmp1[82]+z[65]*tmp1[83]+z[95]*tmp1[84]+z[28]*tmp1[85]+z[58]*tmp1[86]+z[88]*tmp1[87]+z[21]*tmp1[88]+z[51]*tmp1[89]+z[81]*tmp1[90]
                    +z[14]*tmp1[91]+z[44]*tmp1[92]+z[74]*tmp1[93]+z[7]*tmp1[94]+z[37]*tmp1[95]+z[67]*tmp1[96];
                    tab[nb_tmp3+31*nb3]=tab[nb_tmp3+31*nb3]+z[0]*tmp1[0]
                    +z[31]*tmp1[1]+z[62]*tmp1[2]+z[93]*tmp1[3]+z[27]*tmp1[4]+z[58]*tmp1[5]+z[89]*tmp1[6]+z[23]*tmp1[7]+z[54]*tmp1[8]+z[85]*tmp1[9]+z[19]*tmp1[10]
                    +z[50]*tmp1[11]+z[81]*tmp1[12]+z[15]*tmp1[13]+z[46]*tmp1[14]+z[77]*tmp1[15]+z[11]*tmp1[16]+z[42]*tmp1[17]+z[73]*tmp1[18]+z[7]*tmp1[19]+z[38]*tmp1[20]
                    +z[69]*tmp1[21]+z[3]*tmp1[22]+z[34]*tmp1[23]+z[65]*tmp1[24]+z[96]*tmp1[25]+z[30]*tmp1[26]+z[61]*tmp1[27]+z[92]*tmp1[28]+z[26]*tmp1[29]+z[57]*tmp1[30]
                    +z[88]*tmp1[31]+z[22]*tmp1[32]+z[53]*tmp1[33]+z[84]*tmp1[34]+z[18]*tmp1[35]+z[49]*tmp1[36]+z[80]*tmp1[37]+z[14]*tmp1[38]+z[45]*tmp1[39]+z[76]*tmp1[40]
                    +z[10]*tmp1[41]+z[41]*tmp1[42]+z[72]*tmp1[43]+z[6]*tmp1[44]+z[37]*tmp1[45]+z[68]*tmp1[46]+z[2]*tmp1[47]+z[33]*tmp1[48]+z[64]*tmp1[49]+z[95]*tmp1[50]
                    +z[29]*tmp1[51]+z[60]*tmp1[52]+z[91]*tmp1[53]+z[25]*tmp1[54]+z[56]*tmp1[55]+z[87]*tmp1[56]+z[21]*tmp1[57]+z[52]*tmp1[58]+z[83]*tmp1[59]+z[17]*tmp1[60]
                    +z[48]*tmp1[61]+z[79]*tmp1[62]+z[13]*tmp1[63]+z[44]*tmp1[64]+z[75]*tmp1[65]+z[9]*tmp1[66]+z[40]*tmp1[67]+z[71]*tmp1[68]+z[5]*tmp1[69]+z[36]*tmp1[70]
                    +z[67]*tmp1[71]+z[1]*tmp1[72]+z[32]*tmp1[73]+z[63]*tmp1[74]+z[94]*tmp1[75]+z[28]*tmp1[76]+z[59]*tmp1[77]+z[90]*tmp1[78]+z[24]*tmp1[79]+z[55]*tmp1[80]
                    +z[86]*tmp1[81]+z[20]*tmp1[82]+z[51]*tmp1[83]+z[82]*tmp1[84]+z[16]*tmp1[85]+z[47]*tmp1[86]+z[78]*tmp1[87]+z[12]*tmp1[88]+z[43]*tmp1[89]+z[74]*tmp1[90]
                    +z[8]*tmp1[91]+z[39]*tmp1[92]+z[70]*tmp1[93]+z[4]*tmp1[94]+z[35]*tmp1[95]+z[66]*tmp1[96];
                    tab[nb_tmp3+32*nb3]=tab[nb_tmp3+32*nb3]+z[0]*tmp1[0]
                    +z[32]*tmp1[1]+z[64]*tmp1[2]+z[96]*tmp1[3]+z[31]*tmp1[4]+z[63]*tmp1[5]+z[95]*tmp1[6]+z[30]*tmp1[7]+z[62]*tmp1[8]+z[94]*tmp1[9]+z[29]*tmp1[10]
                    +z[61]*tmp1[11]+z[93]*tmp1[12]+z[28]*tmp1[13]+z[60]*tmp1[14]+z[92]*tmp1[15]+z[27]*tmp1[16]+z[59]*tmp1[17]+z[91]*tmp1[18]+z[26]*tmp1[19]+z[58]*tmp1[20]
                    +z[90]*tmp1[21]+z[25]*tmp1[22]+z[57]*tmp1[23]+z[89]*tmp1[24]+z[24]*tmp1[25]+z[56]*tmp1[26]+z[88]*tmp1[27]+z[23]*tmp1[28]+z[55]*tmp1[29]+z[87]*tmp1[30]
                    +z[22]*tmp1[31]+z[54]*tmp1[32]+z[86]*tmp1[33]+z[21]*tmp1[34]+z[53]*tmp1[35]+z[85]*tmp1[36]+z[20]*tmp1[37]+z[52]*tmp1[38]+z[84]*tmp1[39]+z[19]*tmp1[40]
                    +z[51]*tmp1[41]+z[83]*tmp1[42]+z[18]*tmp1[43]+z[50]*tmp1[44]+z[82]*tmp1[45]+z[17]*tmp1[46]+z[49]*tmp1[47]+z[81]*tmp1[48]+z[16]*tmp1[49]+z[48]*tmp1[50]
                    +z[80]*tmp1[51]+z[15]*tmp1[52]+z[47]*tmp1[53]+z[79]*tmp1[54]+z[14]*tmp1[55]+z[46]*tmp1[56]+z[78]*tmp1[57]+z[13]*tmp1[58]+z[45]*tmp1[59]+z[77]*tmp1[60]
                    +z[12]*tmp1[61]+z[44]*tmp1[62]+z[76]*tmp1[63]+z[11]*tmp1[64]+z[43]*tmp1[65]+z[75]*tmp1[66]+z[10]*tmp1[67]+z[42]*tmp1[68]+z[74]*tmp1[69]+z[9]*tmp1[70]
                    +z[41]*tmp1[71]+z[73]*tmp1[72]+z[8]*tmp1[73]+z[40]*tmp1[74]+z[72]*tmp1[75]+z[7]*tmp1[76]+z[39]*tmp1[77]+z[71]*tmp1[78]+z[6]*tmp1[79]+z[38]*tmp1[80]
                    +z[70]*tmp1[81]+z[5]*tmp1[82]+z[37]*tmp1[83]+z[69]*tmp1[84]+z[4]*tmp1[85]+z[36]*tmp1[86]+z[68]*tmp1[87]+z[3]*tmp1[88]+z[35]*tmp1[89]+z[67]*tmp1[90]
                    +z[2]*tmp1[91]+z[34]*tmp1[92]+z[66]*tmp1[93]+z[1]*tmp1[94]+z[33]*tmp1[95]+z[65]*tmp1[96];
                    tab[nb_tmp3+33*nb3]=tab[nb_tmp3+33*nb3]+z[0]*tmp1[0]
                    +z[33]*tmp1[1]+z[66]*tmp1[2]+z[2]*tmp1[3]+z[35]*tmp1[4]+z[68]*tmp1[5]+z[4]*tmp1[6]+z[37]*tmp1[7]+z[70]*tmp1[8]+z[6]*tmp1[9]+z[39]*tmp1[10]
                    +z[72]*tmp1[11]+z[8]*tmp1[12]+z[41]*tmp1[13]+z[74]*tmp1[14]+z[10]*tmp1[15]+z[43]*tmp1[16]+z[76]*tmp1[17]+z[12]*tmp1[18]+z[45]*tmp1[19]+z[78]*tmp1[20]
                    +z[14]*tmp1[21]+z[47]*tmp1[22]+z[80]*tmp1[23]+z[16]*tmp1[24]+z[49]*tmp1[25]+z[82]*tmp1[26]+z[18]*tmp1[27]+z[51]*tmp1[28]+z[84]*tmp1[29]+z[20]*tmp1[30]
                    +z[53]*tmp1[31]+z[86]*tmp1[32]+z[22]*tmp1[33]+z[55]*tmp1[34]+z[88]*tmp1[35]+z[24]*tmp1[36]+z[57]*tmp1[37]+z[90]*tmp1[38]+z[26]*tmp1[39]+z[59]*tmp1[40]
                    +z[92]*tmp1[41]+z[28]*tmp1[42]+z[61]*tmp1[43]+z[94]*tmp1[44]+z[30]*tmp1[45]+z[63]*tmp1[46]+z[96]*tmp1[47]+z[32]*tmp1[48]+z[65]*tmp1[49]+z[1]*tmp1[50]
                    +z[34]*tmp1[51]+z[67]*tmp1[52]+z[3]*tmp1[53]+z[36]*tmp1[54]+z[69]*tmp1[55]+z[5]*tmp1[56]+z[38]*tmp1[57]+z[71]*tmp1[58]+z[7]*tmp1[59]+z[40]*tmp1[60]
                    +z[73]*tmp1[61]+z[9]*tmp1[62]+z[42]*tmp1[63]+z[75]*tmp1[64]+z[11]*tmp1[65]+z[44]*tmp1[66]+z[77]*tmp1[67]+z[13]*tmp1[68]+z[46]*tmp1[69]+z[79]*tmp1[70]
                    +z[15]*tmp1[71]+z[48]*tmp1[72]+z[81]*tmp1[73]+z[17]*tmp1[74]+z[50]*tmp1[75]+z[83]*tmp1[76]+z[19]*tmp1[77]+z[52]*tmp1[78]+z[85]*tmp1[79]+z[21]*tmp1[80]
                    +z[54]*tmp1[81]+z[87]*tmp1[82]+z[23]*tmp1[83]+z[56]*tmp1[84]+z[89]*tmp1[85]+z[25]*tmp1[86]+z[58]*tmp1[87]+z[91]*tmp1[88]+z[27]*tmp1[89]+z[60]*tmp1[90]
                    +z[93]*tmp1[91]+z[29]*tmp1[92]+z[62]*tmp1[93]+z[95]*tmp1[94]+z[31]*tmp1[95]+z[64]*tmp1[96];
                    tab[nb_tmp3+34*nb3]=tab[nb_tmp3+34*nb3]+z[0]*tmp1[0]
                    +z[34]*tmp1[1]+z[68]*tmp1[2]+z[5]*tmp1[3]+z[39]*tmp1[4]+z[73]*tmp1[5]+z[10]*tmp1[6]+z[44]*tmp1[7]+z[78]*tmp1[8]+z[15]*tmp1[9]+z[49]*tmp1[10]
                    +z[83]*tmp1[11]+z[20]*tmp1[12]+z[54]*tmp1[13]+z[88]*tmp1[14]+z[25]*tmp1[15]+z[59]*tmp1[16]+z[93]*tmp1[17]+z[30]*tmp1[18]+z[64]*tmp1[19]+z[1]*tmp1[20]
                    +z[35]*tmp1[21]+z[69]*tmp1[22]+z[6]*tmp1[23]+z[40]*tmp1[24]+z[74]*tmp1[25]+z[11]*tmp1[26]+z[45]*tmp1[27]+z[79]*tmp1[28]+z[16]*tmp1[29]+z[50]*tmp1[30]
                    +z[84]*tmp1[31]+z[21]*tmp1[32]+z[55]*tmp1[33]+z[89]*tmp1[34]+z[26]*tmp1[35]+z[60]*tmp1[36]+z[94]*tmp1[37]+z[31]*tmp1[38]+z[65]*tmp1[39]+z[2]*tmp1[40]
                    +z[36]*tmp1[41]+z[70]*tmp1[42]+z[7]*tmp1[43]+z[41]*tmp1[44]+z[75]*tmp1[45]+z[12]*tmp1[46]+z[46]*tmp1[47]+z[80]*tmp1[48]+z[17]*tmp1[49]+z[51]*tmp1[50]
                    +z[85]*tmp1[51]+z[22]*tmp1[52]+z[56]*tmp1[53]+z[90]*tmp1[54]+z[27]*tmp1[55]+z[61]*tmp1[56]+z[95]*tmp1[57]+z[32]*tmp1[58]+z[66]*tmp1[59]+z[3]*tmp1[60]
                    +z[37]*tmp1[61]+z[71]*tmp1[62]+z[8]*tmp1[63]+z[42]*tmp1[64]+z[76]*tmp1[65]+z[13]*tmp1[66]+z[47]*tmp1[67]+z[81]*tmp1[68]+z[18]*tmp1[69]+z[52]*tmp1[70]
                    +z[86]*tmp1[71]+z[23]*tmp1[72]+z[57]*tmp1[73]+z[91]*tmp1[74]+z[28]*tmp1[75]+z[62]*tmp1[76]+z[96]*tmp1[77]+z[33]*tmp1[78]+z[67]*tmp1[79]+z[4]*tmp1[80]
                    +z[38]*tmp1[81]+z[72]*tmp1[82]+z[9]*tmp1[83]+z[43]*tmp1[84]+z[77]*tmp1[85]+z[14]*tmp1[86]+z[48]*tmp1[87]+z[82]*tmp1[88]+z[19]*tmp1[89]+z[53]*tmp1[90]
                    +z[87]*tmp1[91]+z[24]*tmp1[92]+z[58]*tmp1[93]+z[92]*tmp1[94]+z[29]*tmp1[95]+z[63]*tmp1[96];
                    tab[nb_tmp3+35*nb3]=tab[nb_tmp3+35*nb3]+z[0]*tmp1[0]
                    +z[35]*tmp1[1]+z[70]*tmp1[2]+z[8]*tmp1[3]+z[43]*tmp1[4]+z[78]*tmp1[5]+z[16]*tmp1[6]+z[51]*tmp1[7]+z[86]*tmp1[8]+z[24]*tmp1[9]+z[59]*tmp1[10]
                    +z[94]*tmp1[11]+z[32]*tmp1[12]+z[67]*tmp1[13]+z[5]*tmp1[14]+z[40]*tmp1[15]+z[75]*tmp1[16]+z[13]*tmp1[17]+z[48]*tmp1[18]+z[83]*tmp1[19]+z[21]*tmp1[20]
                    +z[56]*tmp1[21]+z[91]*tmp1[22]+z[29]*tmp1[23]+z[64]*tmp1[24]+z[2]*tmp1[25]+z[37]*tmp1[26]+z[72]*tmp1[27]+z[10]*tmp1[28]+z[45]*tmp1[29]+z[80]*tmp1[30]
                    +z[18]*tmp1[31]+z[53]*tmp1[32]+z[88]*tmp1[33]+z[26]*tmp1[34]+z[61]*tmp1[35]+z[96]*tmp1[36]+z[34]*tmp1[37]+z[69]*tmp1[38]+z[7]*tmp1[39]+z[42]*tmp1[40]
                    +z[77]*tmp1[41]+z[15]*tmp1[42]+z[50]*tmp1[43]+z[85]*tmp1[44]+z[23]*tmp1[45]+z[58]*tmp1[46]+z[93]*tmp1[47]+z[31]*tmp1[48]+z[66]*tmp1[49]+z[4]*tmp1[50]
                    +z[39]*tmp1[51]+z[74]*tmp1[52]+z[12]*tmp1[53]+z[47]*tmp1[54]+z[82]*tmp1[55]+z[20]*tmp1[56]+z[55]*tmp1[57]+z[90]*tmp1[58]+z[28]*tmp1[59]+z[63]*tmp1[60]
                    +z[1]*tmp1[61]+z[36]*tmp1[62]+z[71]*tmp1[63]+z[9]*tmp1[64]+z[44]*tmp1[65]+z[79]*tmp1[66]+z[17]*tmp1[67]+z[52]*tmp1[68]+z[87]*tmp1[69]+z[25]*tmp1[70]
                    +z[60]*tmp1[71]+z[95]*tmp1[72]+z[33]*tmp1[73]+z[68]*tmp1[74]+z[6]*tmp1[75]+z[41]*tmp1[76]+z[76]*tmp1[77]+z[14]*tmp1[78]+z[49]*tmp1[79]+z[84]*tmp1[80]
                    +z[22]*tmp1[81]+z[57]*tmp1[82]+z[92]*tmp1[83]+z[30]*tmp1[84]+z[65]*tmp1[85]+z[3]*tmp1[86]+z[38]*tmp1[87]+z[73]*tmp1[88]+z[11]*tmp1[89]+z[46]*tmp1[90]
                    +z[81]*tmp1[91]+z[19]*tmp1[92]+z[54]*tmp1[93]+z[89]*tmp1[94]+z[27]*tmp1[95]+z[62]*tmp1[96];
                    tab[nb_tmp3+36*nb3]=tab[nb_tmp3+36*nb3]+z[0]*tmp1[0]
                    +z[36]*tmp1[1]+z[72]*tmp1[2]+z[11]*tmp1[3]+z[47]*tmp1[4]+z[83]*tmp1[5]+z[22]*tmp1[6]+z[58]*tmp1[7]+z[94]*tmp1[8]+z[33]*tmp1[9]+z[69]*tmp1[10]
                    +z[8]*tmp1[11]+z[44]*tmp1[12]+z[80]*tmp1[13]+z[19]*tmp1[14]+z[55]*tmp1[15]+z[91]*tmp1[16]+z[30]*tmp1[17]+z[66]*tmp1[18]+z[5]*tmp1[19]+z[41]*tmp1[20]
                    +z[77]*tmp1[21]+z[16]*tmp1[22]+z[52]*tmp1[23]+z[88]*tmp1[24]+z[27]*tmp1[25]+z[63]*tmp1[26]+z[2]*tmp1[27]+z[38]*tmp1[28]+z[74]*tmp1[29]+z[13]*tmp1[30]
                    +z[49]*tmp1[31]+z[85]*tmp1[32]+z[24]*tmp1[33]+z[60]*tmp1[34]+z[96]*tmp1[35]+z[35]*tmp1[36]+z[71]*tmp1[37]+z[10]*tmp1[38]+z[46]*tmp1[39]+z[82]*tmp1[40]
                    +z[21]*tmp1[41]+z[57]*tmp1[42]+z[93]*tmp1[43]+z[32]*tmp1[44]+z[68]*tmp1[45]+z[7]*tmp1[46]+z[43]*tmp1[47]+z[79]*tmp1[48]+z[18]*tmp1[49]+z[54]*tmp1[50]
                    +z[90]*tmp1[51]+z[29]*tmp1[52]+z[65]*tmp1[53]+z[4]*tmp1[54]+z[40]*tmp1[55]+z[76]*tmp1[56]+z[15]*tmp1[57]+z[51]*tmp1[58]+z[87]*tmp1[59]+z[26]*tmp1[60]
                    +z[62]*tmp1[61]+z[1]*tmp1[62]+z[37]*tmp1[63]+z[73]*tmp1[64]+z[12]*tmp1[65]+z[48]*tmp1[66]+z[84]*tmp1[67]+z[23]*tmp1[68]+z[59]*tmp1[69]+z[95]*tmp1[70]
                    +z[34]*tmp1[71]+z[70]*tmp1[72]+z[9]*tmp1[73]+z[45]*tmp1[74]+z[81]*tmp1[75]+z[20]*tmp1[76]+z[56]*tmp1[77]+z[92]*tmp1[78]+z[31]*tmp1[79]+z[67]*tmp1[80]
                    +z[6]*tmp1[81]+z[42]*tmp1[82]+z[78]*tmp1[83]+z[17]*tmp1[84]+z[53]*tmp1[85]+z[89]*tmp1[86]+z[28]*tmp1[87]+z[64]*tmp1[88]+z[3]*tmp1[89]+z[39]*tmp1[90]
                    +z[75]*tmp1[91]+z[14]*tmp1[92]+z[50]*tmp1[93]+z[86]*tmp1[94]+z[25]*tmp1[95]+z[61]*tmp1[96];
                    tab[nb_tmp3+37*nb3]=tab[nb_tmp3+37*nb3]+z[0]*tmp1[0]
                    +z[37]*tmp1[1]+z[74]*tmp1[2]+z[14]*tmp1[3]+z[51]*tmp1[4]+z[88]*tmp1[5]+z[28]*tmp1[6]+z[65]*tmp1[7]+z[5]*tmp1[8]+z[42]*tmp1[9]+z[79]*tmp1[10]
                    +z[19]*tmp1[11]+z[56]*tmp1[12]+z[93]*tmp1[13]+z[33]*tmp1[14]+z[70]*tmp1[15]+z[10]*tmp1[16]+z[47]*tmp1[17]+z[84]*tmp1[18]+z[24]*tmp1[19]+z[61]*tmp1[20]
                    +z[1]*tmp1[21]+z[38]*tmp1[22]+z[75]*tmp1[23]+z[15]*tmp1[24]+z[52]*tmp1[25]+z[89]*tmp1[26]+z[29]*tmp1[27]+z[66]*tmp1[28]+z[6]*tmp1[29]+z[43]*tmp1[30]
                    +z[80]*tmp1[31]+z[20]*tmp1[32]+z[57]*tmp1[33]+z[94]*tmp1[34]+z[34]*tmp1[35]+z[71]*tmp1[36]+z[11]*tmp1[37]+z[48]*tmp1[38]+z[85]*tmp1[39]+z[25]*tmp1[40]
                    +z[62]*tmp1[41]+z[2]*tmp1[42]+z[39]*tmp1[43]+z[76]*tmp1[44]+z[16]*tmp1[45]+z[53]*tmp1[46]+z[90]*tmp1[47]+z[30]*tmp1[48]+z[67]*tmp1[49]+z[7]*tmp1[50]
                    +z[44]*tmp1[51]+z[81]*tmp1[52]+z[21]*tmp1[53]+z[58]*tmp1[54]+z[95]*tmp1[55]+z[35]*tmp1[56]+z[72]*tmp1[57]+z[12]*tmp1[58]+z[49]*tmp1[59]+z[86]*tmp1[60]
                    +z[26]*tmp1[61]+z[63]*tmp1[62]+z[3]*tmp1[63]+z[40]*tmp1[64]+z[77]*tmp1[65]+z[17]*tmp1[66]+z[54]*tmp1[67]+z[91]*tmp1[68]+z[31]*tmp1[69]+z[68]*tmp1[70]
                    +z[8]*tmp1[71]+z[45]*tmp1[72]+z[82]*tmp1[73]+z[22]*tmp1[74]+z[59]*tmp1[75]+z[96]*tmp1[76]+z[36]*tmp1[77]+z[73]*tmp1[78]+z[13]*tmp1[79]+z[50]*tmp1[80]
                    +z[87]*tmp1[81]+z[27]*tmp1[82]+z[64]*tmp1[83]+z[4]*tmp1[84]+z[41]*tmp1[85]+z[78]*tmp1[86]+z[18]*tmp1[87]+z[55]*tmp1[88]+z[92]*tmp1[89]+z[32]*tmp1[90]
                    +z[69]*tmp1[91]+z[9]*tmp1[92]+z[46]*tmp1[93]+z[83]*tmp1[94]+z[23]*tmp1[95]+z[60]*tmp1[96];
                    tab[nb_tmp3+38*nb3]=tab[nb_tmp3+38*nb3]+z[0]*tmp1[0]
                    +z[38]*tmp1[1]+z[76]*tmp1[2]+z[17]*tmp1[3]+z[55]*tmp1[4]+z[93]*tmp1[5]+z[34]*tmp1[6]+z[72]*tmp1[7]+z[13]*tmp1[8]+z[51]*tmp1[9]+z[89]*tmp1[10]
                    +z[30]*tmp1[11]+z[68]*tmp1[12]+z[9]*tmp1[13]+z[47]*tmp1[14]+z[85]*tmp1[15]+z[26]*tmp1[16]+z[64]*tmp1[17]+z[5]*tmp1[18]+z[43]*tmp1[19]+z[81]*tmp1[20]
                    +z[22]*tmp1[21]+z[60]*tmp1[22]+z[1]*tmp1[23]+z[39]*tmp1[24]+z[77]*tmp1[25]+z[18]*tmp1[26]+z[56]*tmp1[27]+z[94]*tmp1[28]+z[35]*tmp1[29]+z[73]*tmp1[30]
                    +z[14]*tmp1[31]+z[52]*tmp1[32]+z[90]*tmp1[33]+z[31]*tmp1[34]+z[69]*tmp1[35]+z[10]*tmp1[36]+z[48]*tmp1[37]+z[86]*tmp1[38]+z[27]*tmp1[39]+z[65]*tmp1[40]
                    +z[6]*tmp1[41]+z[44]*tmp1[42]+z[82]*tmp1[43]+z[23]*tmp1[44]+z[61]*tmp1[45]+z[2]*tmp1[46]+z[40]*tmp1[47]+z[78]*tmp1[48]+z[19]*tmp1[49]+z[57]*tmp1[50]
                    +z[95]*tmp1[51]+z[36]*tmp1[52]+z[74]*tmp1[53]+z[15]*tmp1[54]+z[53]*tmp1[55]+z[91]*tmp1[56]+z[32]*tmp1[57]+z[70]*tmp1[58]+z[11]*tmp1[59]+z[49]*tmp1[60]
                    +z[87]*tmp1[61]+z[28]*tmp1[62]+z[66]*tmp1[63]+z[7]*tmp1[64]+z[45]*tmp1[65]+z[83]*tmp1[66]+z[24]*tmp1[67]+z[62]*tmp1[68]+z[3]*tmp1[69]+z[41]*tmp1[70]
                    +z[79]*tmp1[71]+z[20]*tmp1[72]+z[58]*tmp1[73]+z[96]*tmp1[74]+z[37]*tmp1[75]+z[75]*tmp1[76]+z[16]*tmp1[77]+z[54]*tmp1[78]+z[92]*tmp1[79]+z[33]*tmp1[80]
                    +z[71]*tmp1[81]+z[12]*tmp1[82]+z[50]*tmp1[83]+z[88]*tmp1[84]+z[29]*tmp1[85]+z[67]*tmp1[86]+z[8]*tmp1[87]+z[46]*tmp1[88]+z[84]*tmp1[89]+z[25]*tmp1[90]
                    +z[63]*tmp1[91]+z[4]*tmp1[92]+z[42]*tmp1[93]+z[80]*tmp1[94]+z[21]*tmp1[95]+z[59]*tmp1[96];
                    tab[nb_tmp3+39*nb3]=tab[nb_tmp3+39*nb3]+z[0]*tmp1[0]
                    +z[39]*tmp1[1]+z[78]*tmp1[2]+z[20]*tmp1[3]+z[59]*tmp1[4]+z[1]*tmp1[5]+z[40]*tmp1[6]+z[79]*tmp1[7]+z[21]*tmp1[8]+z[60]*tmp1[9]+z[2]*tmp1[10]
                    +z[41]*tmp1[11]+z[80]*tmp1[12]+z[22]*tmp1[13]+z[61]*tmp1[14]+z[3]*tmp1[15]+z[42]*tmp1[16]+z[81]*tmp1[17]+z[23]*tmp1[18]+z[62]*tmp1[19]+z[4]*tmp1[20]
                    +z[43]*tmp1[21]+z[82]*tmp1[22]+z[24]*tmp1[23]+z[63]*tmp1[24]+z[5]*tmp1[25]+z[44]*tmp1[26]+z[83]*tmp1[27]+z[25]*tmp1[28]+z[64]*tmp1[29]+z[6]*tmp1[30]
                    +z[45]*tmp1[31]+z[84]*tmp1[32]+z[26]*tmp1[33]+z[65]*tmp1[34]+z[7]*tmp1[35]+z[46]*tmp1[36]+z[85]*tmp1[37]+z[27]*tmp1[38]+z[66]*tmp1[39]+z[8]*tmp1[40]
                    +z[47]*tmp1[41]+z[86]*tmp1[42]+z[28]*tmp1[43]+z[67]*tmp1[44]+z[9]*tmp1[45]+z[48]*tmp1[46]+z[87]*tmp1[47]+z[29]*tmp1[48]+z[68]*tmp1[49]+z[10]*tmp1[50]
                    +z[49]*tmp1[51]+z[88]*tmp1[52]+z[30]*tmp1[53]+z[69]*tmp1[54]+z[11]*tmp1[55]+z[50]*tmp1[56]+z[89]*tmp1[57]+z[31]*tmp1[58]+z[70]*tmp1[59]+z[12]*tmp1[60]
                    +z[51]*tmp1[61]+z[90]*tmp1[62]+z[32]*tmp1[63]+z[71]*tmp1[64]+z[13]*tmp1[65]+z[52]*tmp1[66]+z[91]*tmp1[67]+z[33]*tmp1[68]+z[72]*tmp1[69]+z[14]*tmp1[70]
                    +z[53]*tmp1[71]+z[92]*tmp1[72]+z[34]*tmp1[73]+z[73]*tmp1[74]+z[15]*tmp1[75]+z[54]*tmp1[76]+z[93]*tmp1[77]+z[35]*tmp1[78]+z[74]*tmp1[79]+z[16]*tmp1[80]
                    +z[55]*tmp1[81]+z[94]*tmp1[82]+z[36]*tmp1[83]+z[75]*tmp1[84]+z[17]*tmp1[85]+z[56]*tmp1[86]+z[95]*tmp1[87]+z[37]*tmp1[88]+z[76]*tmp1[89]+z[18]*tmp1[90]
                    +z[57]*tmp1[91]+z[96]*tmp1[92]+z[38]*tmp1[93]+z[77]*tmp1[94]+z[19]*tmp1[95]+z[58]*tmp1[96];
                    tab[nb_tmp3+40*nb3]=tab[nb_tmp3+40*nb3]+z[0]*tmp1[0]
                    +z[40]*tmp1[1]+z[80]*tmp1[2]+z[23]*tmp1[3]+z[63]*tmp1[4]+z[6]*tmp1[5]+z[46]*tmp1[6]+z[86]*tmp1[7]+z[29]*tmp1[8]+z[69]*tmp1[9]+z[12]*tmp1[10]
                    +z[52]*tmp1[11]+z[92]*tmp1[12]+z[35]*tmp1[13]+z[75]*tmp1[14]+z[18]*tmp1[15]+z[58]*tmp1[16]+z[1]*tmp1[17]+z[41]*tmp1[18]+z[81]*tmp1[19]+z[24]*tmp1[20]
                    +z[64]*tmp1[21]+z[7]*tmp1[22]+z[47]*tmp1[23]+z[87]*tmp1[24]+z[30]*tmp1[25]+z[70]*tmp1[26]+z[13]*tmp1[27]+z[53]*tmp1[28]+z[93]*tmp1[29]+z[36]*tmp1[30]
                    +z[76]*tmp1[31]+z[19]*tmp1[32]+z[59]*tmp1[33]+z[2]*tmp1[34]+z[42]*tmp1[35]+z[82]*tmp1[36]+z[25]*tmp1[37]+z[65]*tmp1[38]+z[8]*tmp1[39]+z[48]*tmp1[40]
                    +z[88]*tmp1[41]+z[31]*tmp1[42]+z[71]*tmp1[43]+z[14]*tmp1[44]+z[54]*tmp1[45]+z[94]*tmp1[46]+z[37]*tmp1[47]+z[77]*tmp1[48]+z[20]*tmp1[49]+z[60]*tmp1[50]
                    +z[3]*tmp1[51]+z[43]*tmp1[52]+z[83]*tmp1[53]+z[26]*tmp1[54]+z[66]*tmp1[55]+z[9]*tmp1[56]+z[49]*tmp1[57]+z[89]*tmp1[58]+z[32]*tmp1[59]+z[72]*tmp1[60]
                    +z[15]*tmp1[61]+z[55]*tmp1[62]+z[95]*tmp1[63]+z[38]*tmp1[64]+z[78]*tmp1[65]+z[21]*tmp1[66]+z[61]*tmp1[67]+z[4]*tmp1[68]+z[44]*tmp1[69]+z[84]*tmp1[70]
                    +z[27]*tmp1[71]+z[67]*tmp1[72]+z[10]*tmp1[73]+z[50]*tmp1[74]+z[90]*tmp1[75]+z[33]*tmp1[76]+z[73]*tmp1[77]+z[16]*tmp1[78]+z[56]*tmp1[79]+z[96]*tmp1[80]
                    +z[39]*tmp1[81]+z[79]*tmp1[82]+z[22]*tmp1[83]+z[62]*tmp1[84]+z[5]*tmp1[85]+z[45]*tmp1[86]+z[85]*tmp1[87]+z[28]*tmp1[88]+z[68]*tmp1[89]+z[11]*tmp1[90]
                    +z[51]*tmp1[91]+z[91]*tmp1[92]+z[34]*tmp1[93]+z[74]*tmp1[94]+z[17]*tmp1[95]+z[57]*tmp1[96];
                    tab[nb_tmp3+41*nb3]=tab[nb_tmp3+41*nb3]+z[0]*tmp1[0]
                    +z[41]*tmp1[1]+z[82]*tmp1[2]+z[26]*tmp1[3]+z[67]*tmp1[4]+z[11]*tmp1[5]+z[52]*tmp1[6]+z[93]*tmp1[7]+z[37]*tmp1[8]+z[78]*tmp1[9]+z[22]*tmp1[10]
                    +z[63]*tmp1[11]+z[7]*tmp1[12]+z[48]*tmp1[13]+z[89]*tmp1[14]+z[33]*tmp1[15]+z[74]*tmp1[16]+z[18]*tmp1[17]+z[59]*tmp1[18]+z[3]*tmp1[19]+z[44]*tmp1[20]
                    +z[85]*tmp1[21]+z[29]*tmp1[22]+z[70]*tmp1[23]+z[14]*tmp1[24]+z[55]*tmp1[25]+z[96]*tmp1[26]+z[40]*tmp1[27]+z[81]*tmp1[28]+z[25]*tmp1[29]+z[66]*tmp1[30]
                    +z[10]*tmp1[31]+z[51]*tmp1[32]+z[92]*tmp1[33]+z[36]*tmp1[34]+z[77]*tmp1[35]+z[21]*tmp1[36]+z[62]*tmp1[37]+z[6]*tmp1[38]+z[47]*tmp1[39]+z[88]*tmp1[40]
                    +z[32]*tmp1[41]+z[73]*tmp1[42]+z[17]*tmp1[43]+z[58]*tmp1[44]+z[2]*tmp1[45]+z[43]*tmp1[46]+z[84]*tmp1[47]+z[28]*tmp1[48]+z[69]*tmp1[49]+z[13]*tmp1[50]
                    +z[54]*tmp1[51]+z[95]*tmp1[52]+z[39]*tmp1[53]+z[80]*tmp1[54]+z[24]*tmp1[55]+z[65]*tmp1[56]+z[9]*tmp1[57]+z[50]*tmp1[58]+z[91]*tmp1[59]+z[35]*tmp1[60]
                    +z[76]*tmp1[61]+z[20]*tmp1[62]+z[61]*tmp1[63]+z[5]*tmp1[64]+z[46]*tmp1[65]+z[87]*tmp1[66]+z[31]*tmp1[67]+z[72]*tmp1[68]+z[16]*tmp1[69]+z[57]*tmp1[70]
                    +z[1]*tmp1[71]+z[42]*tmp1[72]+z[83]*tmp1[73]+z[27]*tmp1[74]+z[68]*tmp1[75]+z[12]*tmp1[76]+z[53]*tmp1[77]+z[94]*tmp1[78]+z[38]*tmp1[79]+z[79]*tmp1[80]
                    +z[23]*tmp1[81]+z[64]*tmp1[82]+z[8]*tmp1[83]+z[49]*tmp1[84]+z[90]*tmp1[85]+z[34]*tmp1[86]+z[75]*tmp1[87]+z[19]*tmp1[88]+z[60]*tmp1[89]+z[4]*tmp1[90]
                    +z[45]*tmp1[91]+z[86]*tmp1[92]+z[30]*tmp1[93]+z[71]*tmp1[94]+z[15]*tmp1[95]+z[56]*tmp1[96];
                    tab[nb_tmp3+42*nb3]=tab[nb_tmp3+42*nb3]+z[0]*tmp1[0]
                    +z[42]*tmp1[1]+z[84]*tmp1[2]+z[29]*tmp1[3]+z[71]*tmp1[4]+z[16]*tmp1[5]+z[58]*tmp1[6]+z[3]*tmp1[7]+z[45]*tmp1[8]+z[87]*tmp1[9]+z[32]*tmp1[10]
                    +z[74]*tmp1[11]+z[19]*tmp1[12]+z[61]*tmp1[13]+z[6]*tmp1[14]+z[48]*tmp1[15]+z[90]*tmp1[16]+z[35]*tmp1[17]+z[77]*tmp1[18]+z[22]*tmp1[19]+z[64]*tmp1[20]
                    +z[9]*tmp1[21]+z[51]*tmp1[22]+z[93]*tmp1[23]+z[38]*tmp1[24]+z[80]*tmp1[25]+z[25]*tmp1[26]+z[67]*tmp1[27]+z[12]*tmp1[28]+z[54]*tmp1[29]+z[96]*tmp1[30]
                    +z[41]*tmp1[31]+z[83]*tmp1[32]+z[28]*tmp1[33]+z[70]*tmp1[34]+z[15]*tmp1[35]+z[57]*tmp1[36]+z[2]*tmp1[37]+z[44]*tmp1[38]+z[86]*tmp1[39]+z[31]*tmp1[40]
                    +z[73]*tmp1[41]+z[18]*tmp1[42]+z[60]*tmp1[43]+z[5]*tmp1[44]+z[47]*tmp1[45]+z[89]*tmp1[46]+z[34]*tmp1[47]+z[76]*tmp1[48]+z[21]*tmp1[49]+z[63]*tmp1[50]
                    +z[8]*tmp1[51]+z[50]*tmp1[52]+z[92]*tmp1[53]+z[37]*tmp1[54]+z[79]*tmp1[55]+z[24]*tmp1[56]+z[66]*tmp1[57]+z[11]*tmp1[58]+z[53]*tmp1[59]+z[95]*tmp1[60]
                    +z[40]*tmp1[61]+z[82]*tmp1[62]+z[27]*tmp1[63]+z[69]*tmp1[64]+z[14]*tmp1[65]+z[56]*tmp1[66]+z[1]*tmp1[67]+z[43]*tmp1[68]+z[85]*tmp1[69]+z[30]*tmp1[70]
                    +z[72]*tmp1[71]+z[17]*tmp1[72]+z[59]*tmp1[73]+z[4]*tmp1[74]+z[46]*tmp1[75]+z[88]*tmp1[76]+z[33]*tmp1[77]+z[75]*tmp1[78]+z[20]*tmp1[79]+z[62]*tmp1[80]
                    +z[7]*tmp1[81]+z[49]*tmp1[82]+z[91]*tmp1[83]+z[36]*tmp1[84]+z[78]*tmp1[85]+z[23]*tmp1[86]+z[65]*tmp1[87]+z[10]*tmp1[88]+z[52]*tmp1[89]+z[94]*tmp1[90]
                    +z[39]*tmp1[91]+z[81]*tmp1[92]+z[26]*tmp1[93]+z[68]*tmp1[94]+z[13]*tmp1[95]+z[55]*tmp1[96];
                    tab[nb_tmp3+43*nb3]=tab[nb_tmp3+43*nb3]+z[0]*tmp1[0]
                    +z[43]*tmp1[1]+z[86]*tmp1[2]+z[32]*tmp1[3]+z[75]*tmp1[4]+z[21]*tmp1[5]+z[64]*tmp1[6]+z[10]*tmp1[7]+z[53]*tmp1[8]+z[96]*tmp1[9]+z[42]*tmp1[10]
                    +z[85]*tmp1[11]+z[31]*tmp1[12]+z[74]*tmp1[13]+z[20]*tmp1[14]+z[63]*tmp1[15]+z[9]*tmp1[16]+z[52]*tmp1[17]+z[95]*tmp1[18]+z[41]*tmp1[19]+z[84]*tmp1[20]
                    +z[30]*tmp1[21]+z[73]*tmp1[22]+z[19]*tmp1[23]+z[62]*tmp1[24]+z[8]*tmp1[25]+z[51]*tmp1[26]+z[94]*tmp1[27]+z[40]*tmp1[28]+z[83]*tmp1[29]+z[29]*tmp1[30]
                    +z[72]*tmp1[31]+z[18]*tmp1[32]+z[61]*tmp1[33]+z[7]*tmp1[34]+z[50]*tmp1[35]+z[93]*tmp1[36]+z[39]*tmp1[37]+z[82]*tmp1[38]+z[28]*tmp1[39]+z[71]*tmp1[40]
                    +z[17]*tmp1[41]+z[60]*tmp1[42]+z[6]*tmp1[43]+z[49]*tmp1[44]+z[92]*tmp1[45]+z[38]*tmp1[46]+z[81]*tmp1[47]+z[27]*tmp1[48]+z[70]*tmp1[49]+z[16]*tmp1[50]
                    +z[59]*tmp1[51]+z[5]*tmp1[52]+z[48]*tmp1[53]+z[91]*tmp1[54]+z[37]*tmp1[55]+z[80]*tmp1[56]+z[26]*tmp1[57]+z[69]*tmp1[58]+z[15]*tmp1[59]+z[58]*tmp1[60]
                    +z[4]*tmp1[61]+z[47]*tmp1[62]+z[90]*tmp1[63]+z[36]*tmp1[64]+z[79]*tmp1[65]+z[25]*tmp1[66]+z[68]*tmp1[67]+z[14]*tmp1[68]+z[57]*tmp1[69]+z[3]*tmp1[70]
                    +z[46]*tmp1[71]+z[89]*tmp1[72]+z[35]*tmp1[73]+z[78]*tmp1[74]+z[24]*tmp1[75]+z[67]*tmp1[76]+z[13]*tmp1[77]+z[56]*tmp1[78]+z[2]*tmp1[79]+z[45]*tmp1[80]
                    +z[88]*tmp1[81]+z[34]*tmp1[82]+z[77]*tmp1[83]+z[23]*tmp1[84]+z[66]*tmp1[85]+z[12]*tmp1[86]+z[55]*tmp1[87]+z[1]*tmp1[88]+z[44]*tmp1[89]+z[87]*tmp1[90]
                    +z[33]*tmp1[91]+z[76]*tmp1[92]+z[22]*tmp1[93]+z[65]*tmp1[94]+z[11]*tmp1[95]+z[54]*tmp1[96];
                    tab[nb_tmp3+44*nb3]=tab[nb_tmp3+44*nb3]+z[0]*tmp1[0]
                    +z[44]*tmp1[1]+z[88]*tmp1[2]+z[35]*tmp1[3]+z[79]*tmp1[4]+z[26]*tmp1[5]+z[70]*tmp1[6]+z[17]*tmp1[7]+z[61]*tmp1[8]+z[8]*tmp1[9]+z[52]*tmp1[10]
                    +z[96]*tmp1[11]+z[43]*tmp1[12]+z[87]*tmp1[13]+z[34]*tmp1[14]+z[78]*tmp1[15]+z[25]*tmp1[16]+z[69]*tmp1[17]+z[16]*tmp1[18]+z[60]*tmp1[19]+z[7]*tmp1[20]
                    +z[51]*tmp1[21]+z[95]*tmp1[22]+z[42]*tmp1[23]+z[86]*tmp1[24]+z[33]*tmp1[25]+z[77]*tmp1[26]+z[24]*tmp1[27]+z[68]*tmp1[28]+z[15]*tmp1[29]+z[59]*tmp1[30]
                    +z[6]*tmp1[31]+z[50]*tmp1[32]+z[94]*tmp1[33]+z[41]*tmp1[34]+z[85]*tmp1[35]+z[32]*tmp1[36]+z[76]*tmp1[37]+z[23]*tmp1[38]+z[67]*tmp1[39]+z[14]*tmp1[40]
                    +z[58]*tmp1[41]+z[5]*tmp1[42]+z[49]*tmp1[43]+z[93]*tmp1[44]+z[40]*tmp1[45]+z[84]*tmp1[46]+z[31]*tmp1[47]+z[75]*tmp1[48]+z[22]*tmp1[49]+z[66]*tmp1[50]
                    +z[13]*tmp1[51]+z[57]*tmp1[52]+z[4]*tmp1[53]+z[48]*tmp1[54]+z[92]*tmp1[55]+z[39]*tmp1[56]+z[83]*tmp1[57]+z[30]*tmp1[58]+z[74]*tmp1[59]+z[21]*tmp1[60]
                    +z[65]*tmp1[61]+z[12]*tmp1[62]+z[56]*tmp1[63]+z[3]*tmp1[64]+z[47]*tmp1[65]+z[91]*tmp1[66]+z[38]*tmp1[67]+z[82]*tmp1[68]+z[29]*tmp1[69]+z[73]*tmp1[70]
                    +z[20]*tmp1[71]+z[64]*tmp1[72]+z[11]*tmp1[73]+z[55]*tmp1[74]+z[2]*tmp1[75]+z[46]*tmp1[76]+z[90]*tmp1[77]+z[37]*tmp1[78]+z[81]*tmp1[79]+z[28]*tmp1[80]
                    +z[72]*tmp1[81]+z[19]*tmp1[82]+z[63]*tmp1[83]+z[10]*tmp1[84]+z[54]*tmp1[85]+z[1]*tmp1[86]+z[45]*tmp1[87]+z[89]*tmp1[88]+z[36]*tmp1[89]+z[80]*tmp1[90]
                    +z[27]*tmp1[91]+z[71]*tmp1[92]+z[18]*tmp1[93]+z[62]*tmp1[94]+z[9]*tmp1[95]+z[53]*tmp1[96];
                    tab[nb_tmp3+45*nb3]=tab[nb_tmp3+45*nb3]+z[0]*tmp1[0]
                    +z[45]*tmp1[1]+z[90]*tmp1[2]+z[38]*tmp1[3]+z[83]*tmp1[4]+z[31]*tmp1[5]+z[76]*tmp1[6]+z[24]*tmp1[7]+z[69]*tmp1[8]+z[17]*tmp1[9]+z[62]*tmp1[10]
                    +z[10]*tmp1[11]+z[55]*tmp1[12]+z[3]*tmp1[13]+z[48]*tmp1[14]+z[93]*tmp1[15]+z[41]*tmp1[16]+z[86]*tmp1[17]+z[34]*tmp1[18]+z[79]*tmp1[19]+z[27]*tmp1[20]
                    +z[72]*tmp1[21]+z[20]*tmp1[22]+z[65]*tmp1[23]+z[13]*tmp1[24]+z[58]*tmp1[25]+z[6]*tmp1[26]+z[51]*tmp1[27]+z[96]*tmp1[28]+z[44]*tmp1[29]+z[89]*tmp1[30]
                    +z[37]*tmp1[31]+z[82]*tmp1[32]+z[30]*tmp1[33]+z[75]*tmp1[34]+z[23]*tmp1[35]+z[68]*tmp1[36]+z[16]*tmp1[37]+z[61]*tmp1[38]+z[9]*tmp1[39]+z[54]*tmp1[40]
                    +z[2]*tmp1[41]+z[47]*tmp1[42]+z[92]*tmp1[43]+z[40]*tmp1[44]+z[85]*tmp1[45]+z[33]*tmp1[46]+z[78]*tmp1[47]+z[26]*tmp1[48]+z[71]*tmp1[49]+z[19]*tmp1[50]
                    +z[64]*tmp1[51]+z[12]*tmp1[52]+z[57]*tmp1[53]+z[5]*tmp1[54]+z[50]*tmp1[55]+z[95]*tmp1[56]+z[43]*tmp1[57]+z[88]*tmp1[58]+z[36]*tmp1[59]+z[81]*tmp1[60]
                    +z[29]*tmp1[61]+z[74]*tmp1[62]+z[22]*tmp1[63]+z[67]*tmp1[64]+z[15]*tmp1[65]+z[60]*tmp1[66]+z[8]*tmp1[67]+z[53]*tmp1[68]+z[1]*tmp1[69]+z[46]*tmp1[70]
                    +z[91]*tmp1[71]+z[39]*tmp1[72]+z[84]*tmp1[73]+z[32]*tmp1[74]+z[77]*tmp1[75]+z[25]*tmp1[76]+z[70]*tmp1[77]+z[18]*tmp1[78]+z[63]*tmp1[79]+z[11]*tmp1[80]
                    +z[56]*tmp1[81]+z[4]*tmp1[82]+z[49]*tmp1[83]+z[94]*tmp1[84]+z[42]*tmp1[85]+z[87]*tmp1[86]+z[35]*tmp1[87]+z[80]*tmp1[88]+z[28]*tmp1[89]+z[73]*tmp1[90]
                    +z[21]*tmp1[91]+z[66]*tmp1[92]+z[14]*tmp1[93]+z[59]*tmp1[94]+z[7]*tmp1[95]+z[52]*tmp1[96];
                    tab[nb_tmp3+46*nb3]=tab[nb_tmp3+46*nb3]+z[0]*tmp1[0]
                    +z[46]*tmp1[1]+z[92]*tmp1[2]+z[41]*tmp1[3]+z[87]*tmp1[4]+z[36]*tmp1[5]+z[82]*tmp1[6]+z[31]*tmp1[7]+z[77]*tmp1[8]+z[26]*tmp1[9]+z[72]*tmp1[10]
                    +z[21]*tmp1[11]+z[67]*tmp1[12]+z[16]*tmp1[13]+z[62]*tmp1[14]+z[11]*tmp1[15]+z[57]*tmp1[16]+z[6]*tmp1[17]+z[52]*tmp1[18]+z[1]*tmp1[19]+z[47]*tmp1[20]
                    +z[93]*tmp1[21]+z[42]*tmp1[22]+z[88]*tmp1[23]+z[37]*tmp1[24]+z[83]*tmp1[25]+z[32]*tmp1[26]+z[78]*tmp1[27]+z[27]*tmp1[28]+z[73]*tmp1[29]+z[22]*tmp1[30]
                    +z[68]*tmp1[31]+z[17]*tmp1[32]+z[63]*tmp1[33]+z[12]*tmp1[34]+z[58]*tmp1[35]+z[7]*tmp1[36]+z[53]*tmp1[37]+z[2]*tmp1[38]+z[48]*tmp1[39]+z[94]*tmp1[40]
                    +z[43]*tmp1[41]+z[89]*tmp1[42]+z[38]*tmp1[43]+z[84]*tmp1[44]+z[33]*tmp1[45]+z[79]*tmp1[46]+z[28]*tmp1[47]+z[74]*tmp1[48]+z[23]*tmp1[49]+z[69]*tmp1[50]
                    +z[18]*tmp1[51]+z[64]*tmp1[52]+z[13]*tmp1[53]+z[59]*tmp1[54]+z[8]*tmp1[55]+z[54]*tmp1[56]+z[3]*tmp1[57]+z[49]*tmp1[58]+z[95]*tmp1[59]+z[44]*tmp1[60]
                    +z[90]*tmp1[61]+z[39]*tmp1[62]+z[85]*tmp1[63]+z[34]*tmp1[64]+z[80]*tmp1[65]+z[29]*tmp1[66]+z[75]*tmp1[67]+z[24]*tmp1[68]+z[70]*tmp1[69]+z[19]*tmp1[70]
                    +z[65]*tmp1[71]+z[14]*tmp1[72]+z[60]*tmp1[73]+z[9]*tmp1[74]+z[55]*tmp1[75]+z[4]*tmp1[76]+z[50]*tmp1[77]+z[96]*tmp1[78]+z[45]*tmp1[79]+z[91]*tmp1[80]
                    +z[40]*tmp1[81]+z[86]*tmp1[82]+z[35]*tmp1[83]+z[81]*tmp1[84]+z[30]*tmp1[85]+z[76]*tmp1[86]+z[25]*tmp1[87]+z[71]*tmp1[88]+z[20]*tmp1[89]+z[66]*tmp1[90]
                    +z[15]*tmp1[91]+z[61]*tmp1[92]+z[10]*tmp1[93]+z[56]*tmp1[94]+z[5]*tmp1[95]+z[51]*tmp1[96];
                    tab[nb_tmp3+47*nb3]=tab[nb_tmp3+47*nb3]+z[0]*tmp1[0]
                    +z[47]*tmp1[1]+z[94]*tmp1[2]+z[44]*tmp1[3]+z[91]*tmp1[4]+z[41]*tmp1[5]+z[88]*tmp1[6]+z[38]*tmp1[7]+z[85]*tmp1[8]+z[35]*tmp1[9]+z[82]*tmp1[10]
                    +z[32]*tmp1[11]+z[79]*tmp1[12]+z[29]*tmp1[13]+z[76]*tmp1[14]+z[26]*tmp1[15]+z[73]*tmp1[16]+z[23]*tmp1[17]+z[70]*tmp1[18]+z[20]*tmp1[19]+z[67]*tmp1[20]
                    +z[17]*tmp1[21]+z[64]*tmp1[22]+z[14]*tmp1[23]+z[61]*tmp1[24]+z[11]*tmp1[25]+z[58]*tmp1[26]+z[8]*tmp1[27]+z[55]*tmp1[28]+z[5]*tmp1[29]+z[52]*tmp1[30]
                    +z[2]*tmp1[31]+z[49]*tmp1[32]+z[96]*tmp1[33]+z[46]*tmp1[34]+z[93]*tmp1[35]+z[43]*tmp1[36]+z[90]*tmp1[37]+z[40]*tmp1[38]+z[87]*tmp1[39]+z[37]*tmp1[40]
                    +z[84]*tmp1[41]+z[34]*tmp1[42]+z[81]*tmp1[43]+z[31]*tmp1[44]+z[78]*tmp1[45]+z[28]*tmp1[46]+z[75]*tmp1[47]+z[25]*tmp1[48]+z[72]*tmp1[49]+z[22]*tmp1[50]
                    +z[69]*tmp1[51]+z[19]*tmp1[52]+z[66]*tmp1[53]+z[16]*tmp1[54]+z[63]*tmp1[55]+z[13]*tmp1[56]+z[60]*tmp1[57]+z[10]*tmp1[58]+z[57]*tmp1[59]+z[7]*tmp1[60]
                    +z[54]*tmp1[61]+z[4]*tmp1[62]+z[51]*tmp1[63]+z[1]*tmp1[64]+z[48]*tmp1[65]+z[95]*tmp1[66]+z[45]*tmp1[67]+z[92]*tmp1[68]+z[42]*tmp1[69]+z[89]*tmp1[70]
                    +z[39]*tmp1[71]+z[86]*tmp1[72]+z[36]*tmp1[73]+z[83]*tmp1[74]+z[33]*tmp1[75]+z[80]*tmp1[76]+z[30]*tmp1[77]+z[77]*tmp1[78]+z[27]*tmp1[79]+z[74]*tmp1[80]
                    +z[24]*tmp1[81]+z[71]*tmp1[82]+z[21]*tmp1[83]+z[68]*tmp1[84]+z[18]*tmp1[85]+z[65]*tmp1[86]+z[15]*tmp1[87]+z[62]*tmp1[88]+z[12]*tmp1[89]+z[59]*tmp1[90]
                    +z[9]*tmp1[91]+z[56]*tmp1[92]+z[6]*tmp1[93]+z[53]*tmp1[94]+z[3]*tmp1[95]+z[50]*tmp1[96];
                    tab[nb_tmp3+48*nb3]=tab[nb_tmp3+48*nb3]+z[0]*tmp1[0]
                    +z[48]*tmp1[1]+z[96]*tmp1[2]+z[47]*tmp1[3]+z[95]*tmp1[4]+z[46]*tmp1[5]+z[94]*tmp1[6]+z[45]*tmp1[7]+z[93]*tmp1[8]+z[44]*tmp1[9]+z[92]*tmp1[10]
                    +z[43]*tmp1[11]+z[91]*tmp1[12]+z[42]*tmp1[13]+z[90]*tmp1[14]+z[41]*tmp1[15]+z[89]*tmp1[16]+z[40]*tmp1[17]+z[88]*tmp1[18]+z[39]*tmp1[19]+z[87]*tmp1[20]
                    +z[38]*tmp1[21]+z[86]*tmp1[22]+z[37]*tmp1[23]+z[85]*tmp1[24]+z[36]*tmp1[25]+z[84]*tmp1[26]+z[35]*tmp1[27]+z[83]*tmp1[28]+z[34]*tmp1[29]+z[82]*tmp1[30]
                    +z[33]*tmp1[31]+z[81]*tmp1[32]+z[32]*tmp1[33]+z[80]*tmp1[34]+z[31]*tmp1[35]+z[79]*tmp1[36]+z[30]*tmp1[37]+z[78]*tmp1[38]+z[29]*tmp1[39]+z[77]*tmp1[40]
                    +z[28]*tmp1[41]+z[76]*tmp1[42]+z[27]*tmp1[43]+z[75]*tmp1[44]+z[26]*tmp1[45]+z[74]*tmp1[46]+z[25]*tmp1[47]+z[73]*tmp1[48]+z[24]*tmp1[49]+z[72]*tmp1[50]
                    +z[23]*tmp1[51]+z[71]*tmp1[52]+z[22]*tmp1[53]+z[70]*tmp1[54]+z[21]*tmp1[55]+z[69]*tmp1[56]+z[20]*tmp1[57]+z[68]*tmp1[58]+z[19]*tmp1[59]+z[67]*tmp1[60]
                    +z[18]*tmp1[61]+z[66]*tmp1[62]+z[17]*tmp1[63]+z[65]*tmp1[64]+z[16]*tmp1[65]+z[64]*tmp1[66]+z[15]*tmp1[67]+z[63]*tmp1[68]+z[14]*tmp1[69]+z[62]*tmp1[70]
                    +z[13]*tmp1[71]+z[61]*tmp1[72]+z[12]*tmp1[73]+z[60]*tmp1[74]+z[11]*tmp1[75]+z[59]*tmp1[76]+z[10]*tmp1[77]+z[58]*tmp1[78]+z[9]*tmp1[79]+z[57]*tmp1[80]
                    +z[8]*tmp1[81]+z[56]*tmp1[82]+z[7]*tmp1[83]+z[55]*tmp1[84]+z[6]*tmp1[85]+z[54]*tmp1[86]+z[5]*tmp1[87]+z[53]*tmp1[88]+z[4]*tmp1[89]+z[52]*tmp1[90]
                    +z[3]*tmp1[91]+z[51]*tmp1[92]+z[2]*tmp1[93]+z[50]*tmp1[94]+z[1]*tmp1[95]+z[49]*tmp1[96];
                    tab[nb_tmp3+49*nb3]=tab[nb_tmp3+49*nb3]+z[0]*tmp1[0]
                    +z[49]*tmp1[1]+z[1]*tmp1[2]+z[50]*tmp1[3]+z[2]*tmp1[4]+z[51]*tmp1[5]+z[3]*tmp1[6]+z[52]*tmp1[7]+z[4]*tmp1[8]+z[53]*tmp1[9]+z[5]*tmp1[10]
                    +z[54]*tmp1[11]+z[6]*tmp1[12]+z[55]*tmp1[13]+z[7]*tmp1[14]+z[56]*tmp1[15]+z[8]*tmp1[16]+z[57]*tmp1[17]+z[9]*tmp1[18]+z[58]*tmp1[19]+z[10]*tmp1[20]
                    +z[59]*tmp1[21]+z[11]*tmp1[22]+z[60]*tmp1[23]+z[12]*tmp1[24]+z[61]*tmp1[25]+z[13]*tmp1[26]+z[62]*tmp1[27]+z[14]*tmp1[28]+z[63]*tmp1[29]+z[15]*tmp1[30]
                    +z[64]*tmp1[31]+z[16]*tmp1[32]+z[65]*tmp1[33]+z[17]*tmp1[34]+z[66]*tmp1[35]+z[18]*tmp1[36]+z[67]*tmp1[37]+z[19]*tmp1[38]+z[68]*tmp1[39]+z[20]*tmp1[40]
                    +z[69]*tmp1[41]+z[21]*tmp1[42]+z[70]*tmp1[43]+z[22]*tmp1[44]+z[71]*tmp1[45]+z[23]*tmp1[46]+z[72]*tmp1[47]+z[24]*tmp1[48]+z[73]*tmp1[49]+z[25]*tmp1[50]
                    +z[74]*tmp1[51]+z[26]*tmp1[52]+z[75]*tmp1[53]+z[27]*tmp1[54]+z[76]*tmp1[55]+z[28]*tmp1[56]+z[77]*tmp1[57]+z[29]*tmp1[58]+z[78]*tmp1[59]+z[30]*tmp1[60]
                    +z[79]*tmp1[61]+z[31]*tmp1[62]+z[80]*tmp1[63]+z[32]*tmp1[64]+z[81]*tmp1[65]+z[33]*tmp1[66]+z[82]*tmp1[67]+z[34]*tmp1[68]+z[83]*tmp1[69]+z[35]*tmp1[70]
                    +z[84]*tmp1[71]+z[36]*tmp1[72]+z[85]*tmp1[73]+z[37]*tmp1[74]+z[86]*tmp1[75]+z[38]*tmp1[76]+z[87]*tmp1[77]+z[39]*tmp1[78]+z[88]*tmp1[79]+z[40]*tmp1[80]
                    +z[89]*tmp1[81]+z[41]*tmp1[82]+z[90]*tmp1[83]+z[42]*tmp1[84]+z[91]*tmp1[85]+z[43]*tmp1[86]+z[92]*tmp1[87]+z[44]*tmp1[88]+z[93]*tmp1[89]+z[45]*tmp1[90]
                    +z[94]*tmp1[91]+z[46]*tmp1[92]+z[95]*tmp1[93]+z[47]*tmp1[94]+z[96]*tmp1[95]+z[48]*tmp1[96];
                    tab[nb_tmp3+50*nb3]=tab[nb_tmp3+50*nb3]+z[0]*tmp1[0]
                    +z[50]*tmp1[1]+z[3]*tmp1[2]+z[53]*tmp1[3]+z[6]*tmp1[4]+z[56]*tmp1[5]+z[9]*tmp1[6]+z[59]*tmp1[7]+z[12]*tmp1[8]+z[62]*tmp1[9]+z[15]*tmp1[10]
                    +z[65]*tmp1[11]+z[18]*tmp1[12]+z[68]*tmp1[13]+z[21]*tmp1[14]+z[71]*tmp1[15]+z[24]*tmp1[16]+z[74]*tmp1[17]+z[27]*tmp1[18]+z[77]*tmp1[19]+z[30]*tmp1[20]
                    +z[80]*tmp1[21]+z[33]*tmp1[22]+z[83]*tmp1[23]+z[36]*tmp1[24]+z[86]*tmp1[25]+z[39]*tmp1[26]+z[89]*tmp1[27]+z[42]*tmp1[28]+z[92]*tmp1[29]+z[45]*tmp1[30]
                    +z[95]*tmp1[31]+z[48]*tmp1[32]+z[1]*tmp1[33]+z[51]*tmp1[34]+z[4]*tmp1[35]+z[54]*tmp1[36]+z[7]*tmp1[37]+z[57]*tmp1[38]+z[10]*tmp1[39]+z[60]*tmp1[40]
                    +z[13]*tmp1[41]+z[63]*tmp1[42]+z[16]*tmp1[43]+z[66]*tmp1[44]+z[19]*tmp1[45]+z[69]*tmp1[46]+z[22]*tmp1[47]+z[72]*tmp1[48]+z[25]*tmp1[49]+z[75]*tmp1[50]
                    +z[28]*tmp1[51]+z[78]*tmp1[52]+z[31]*tmp1[53]+z[81]*tmp1[54]+z[34]*tmp1[55]+z[84]*tmp1[56]+z[37]*tmp1[57]+z[87]*tmp1[58]+z[40]*tmp1[59]+z[90]*tmp1[60]
                    +z[43]*tmp1[61]+z[93]*tmp1[62]+z[46]*tmp1[63]+z[96]*tmp1[64]+z[49]*tmp1[65]+z[2]*tmp1[66]+z[52]*tmp1[67]+z[5]*tmp1[68]+z[55]*tmp1[69]+z[8]*tmp1[70]
                    +z[58]*tmp1[71]+z[11]*tmp1[72]+z[61]*tmp1[73]+z[14]*tmp1[74]+z[64]*tmp1[75]+z[17]*tmp1[76]+z[67]*tmp1[77]+z[20]*tmp1[78]+z[70]*tmp1[79]+z[23]*tmp1[80]
                    +z[73]*tmp1[81]+z[26]*tmp1[82]+z[76]*tmp1[83]+z[29]*tmp1[84]+z[79]*tmp1[85]+z[32]*tmp1[86]+z[82]*tmp1[87]+z[35]*tmp1[88]+z[85]*tmp1[89]+z[38]*tmp1[90]
                    +z[88]*tmp1[91]+z[41]*tmp1[92]+z[91]*tmp1[93]+z[44]*tmp1[94]+z[94]*tmp1[95]+z[47]*tmp1[96];
                    tab[nb_tmp3+51*nb3]=tab[nb_tmp3+51*nb3]+z[0]*tmp1[0]
                    +z[51]*tmp1[1]+z[5]*tmp1[2]+z[56]*tmp1[3]+z[10]*tmp1[4]+z[61]*tmp1[5]+z[15]*tmp1[6]+z[66]*tmp1[7]+z[20]*tmp1[8]+z[71]*tmp1[9]+z[25]*tmp1[10]
                    +z[76]*tmp1[11]+z[30]*tmp1[12]+z[81]*tmp1[13]+z[35]*tmp1[14]+z[86]*tmp1[15]+z[40]*tmp1[16]+z[91]*tmp1[17]+z[45]*tmp1[18]+z[96]*tmp1[19]+z[50]*tmp1[20]
                    +z[4]*tmp1[21]+z[55]*tmp1[22]+z[9]*tmp1[23]+z[60]*tmp1[24]+z[14]*tmp1[25]+z[65]*tmp1[26]+z[19]*tmp1[27]+z[70]*tmp1[28]+z[24]*tmp1[29]+z[75]*tmp1[30]
                    +z[29]*tmp1[31]+z[80]*tmp1[32]+z[34]*tmp1[33]+z[85]*tmp1[34]+z[39]*tmp1[35]+z[90]*tmp1[36]+z[44]*tmp1[37]+z[95]*tmp1[38]+z[49]*tmp1[39]+z[3]*tmp1[40]
                    +z[54]*tmp1[41]+z[8]*tmp1[42]+z[59]*tmp1[43]+z[13]*tmp1[44]+z[64]*tmp1[45]+z[18]*tmp1[46]+z[69]*tmp1[47]+z[23]*tmp1[48]+z[74]*tmp1[49]+z[28]*tmp1[50]
                    +z[79]*tmp1[51]+z[33]*tmp1[52]+z[84]*tmp1[53]+z[38]*tmp1[54]+z[89]*tmp1[55]+z[43]*tmp1[56]+z[94]*tmp1[57]+z[48]*tmp1[58]+z[2]*tmp1[59]+z[53]*tmp1[60]
                    +z[7]*tmp1[61]+z[58]*tmp1[62]+z[12]*tmp1[63]+z[63]*tmp1[64]+z[17]*tmp1[65]+z[68]*tmp1[66]+z[22]*tmp1[67]+z[73]*tmp1[68]+z[27]*tmp1[69]+z[78]*tmp1[70]
                    +z[32]*tmp1[71]+z[83]*tmp1[72]+z[37]*tmp1[73]+z[88]*tmp1[74]+z[42]*tmp1[75]+z[93]*tmp1[76]+z[47]*tmp1[77]+z[1]*tmp1[78]+z[52]*tmp1[79]+z[6]*tmp1[80]
                    +z[57]*tmp1[81]+z[11]*tmp1[82]+z[62]*tmp1[83]+z[16]*tmp1[84]+z[67]*tmp1[85]+z[21]*tmp1[86]+z[72]*tmp1[87]+z[26]*tmp1[88]+z[77]*tmp1[89]+z[31]*tmp1[90]
                    +z[82]*tmp1[91]+z[36]*tmp1[92]+z[87]*tmp1[93]+z[41]*tmp1[94]+z[92]*tmp1[95]+z[46]*tmp1[96];
                    tab[nb_tmp3+52*nb3]=tab[nb_tmp3+52*nb3]+z[0]*tmp1[0]
                    +z[52]*tmp1[1]+z[7]*tmp1[2]+z[59]*tmp1[3]+z[14]*tmp1[4]+z[66]*tmp1[5]+z[21]*tmp1[6]+z[73]*tmp1[7]+z[28]*tmp1[8]+z[80]*tmp1[9]+z[35]*tmp1[10]
                    +z[87]*tmp1[11]+z[42]*tmp1[12]+z[94]*tmp1[13]+z[49]*tmp1[14]+z[4]*tmp1[15]+z[56]*tmp1[16]+z[11]*tmp1[17]+z[63]*tmp1[18]+z[18]*tmp1[19]+z[70]*tmp1[20]
                    +z[25]*tmp1[21]+z[77]*tmp1[22]+z[32]*tmp1[23]+z[84]*tmp1[24]+z[39]*tmp1[25]+z[91]*tmp1[26]+z[46]*tmp1[27]+z[1]*tmp1[28]+z[53]*tmp1[29]+z[8]*tmp1[30]
                    +z[60]*tmp1[31]+z[15]*tmp1[32]+z[67]*tmp1[33]+z[22]*tmp1[34]+z[74]*tmp1[35]+z[29]*tmp1[36]+z[81]*tmp1[37]+z[36]*tmp1[38]+z[88]*tmp1[39]+z[43]*tmp1[40]
                    +z[95]*tmp1[41]+z[50]*tmp1[42]+z[5]*tmp1[43]+z[57]*tmp1[44]+z[12]*tmp1[45]+z[64]*tmp1[46]+z[19]*tmp1[47]+z[71]*tmp1[48]+z[26]*tmp1[49]+z[78]*tmp1[50]
                    +z[33]*tmp1[51]+z[85]*tmp1[52]+z[40]*tmp1[53]+z[92]*tmp1[54]+z[47]*tmp1[55]+z[2]*tmp1[56]+z[54]*tmp1[57]+z[9]*tmp1[58]+z[61]*tmp1[59]+z[16]*tmp1[60]
                    +z[68]*tmp1[61]+z[23]*tmp1[62]+z[75]*tmp1[63]+z[30]*tmp1[64]+z[82]*tmp1[65]+z[37]*tmp1[66]+z[89]*tmp1[67]+z[44]*tmp1[68]+z[96]*tmp1[69]+z[51]*tmp1[70]
                    +z[6]*tmp1[71]+z[58]*tmp1[72]+z[13]*tmp1[73]+z[65]*tmp1[74]+z[20]*tmp1[75]+z[72]*tmp1[76]+z[27]*tmp1[77]+z[79]*tmp1[78]+z[34]*tmp1[79]+z[86]*tmp1[80]
                    +z[41]*tmp1[81]+z[93]*tmp1[82]+z[48]*tmp1[83]+z[3]*tmp1[84]+z[55]*tmp1[85]+z[10]*tmp1[86]+z[62]*tmp1[87]+z[17]*tmp1[88]+z[69]*tmp1[89]+z[24]*tmp1[90]
                    +z[76]*tmp1[91]+z[31]*tmp1[92]+z[83]*tmp1[93]+z[38]*tmp1[94]+z[90]*tmp1[95]+z[45]*tmp1[96];
                    tab[nb_tmp3+53*nb3]=tab[nb_tmp3+53*nb3]+z[0]*tmp1[0]
                    +z[53]*tmp1[1]+z[9]*tmp1[2]+z[62]*tmp1[3]+z[18]*tmp1[4]+z[71]*tmp1[5]+z[27]*tmp1[6]+z[80]*tmp1[7]+z[36]*tmp1[8]+z[89]*tmp1[9]+z[45]*tmp1[10]
                    +z[1]*tmp1[11]+z[54]*tmp1[12]+z[10]*tmp1[13]+z[63]*tmp1[14]+z[19]*tmp1[15]+z[72]*tmp1[16]+z[28]*tmp1[17]+z[81]*tmp1[18]+z[37]*tmp1[19]+z[90]*tmp1[20]
                    +z[46]*tmp1[21]+z[2]*tmp1[22]+z[55]*tmp1[23]+z[11]*tmp1[24]+z[64]*tmp1[25]+z[20]*tmp1[26]+z[73]*tmp1[27]+z[29]*tmp1[28]+z[82]*tmp1[29]+z[38]*tmp1[30]
                    +z[91]*tmp1[31]+z[47]*tmp1[32]+z[3]*tmp1[33]+z[56]*tmp1[34]+z[12]*tmp1[35]+z[65]*tmp1[36]+z[21]*tmp1[37]+z[74]*tmp1[38]+z[30]*tmp1[39]+z[83]*tmp1[40]
                    +z[39]*tmp1[41]+z[92]*tmp1[42]+z[48]*tmp1[43]+z[4]*tmp1[44]+z[57]*tmp1[45]+z[13]*tmp1[46]+z[66]*tmp1[47]+z[22]*tmp1[48]+z[75]*tmp1[49]+z[31]*tmp1[50]
                    +z[84]*tmp1[51]+z[40]*tmp1[52]+z[93]*tmp1[53]+z[49]*tmp1[54]+z[5]*tmp1[55]+z[58]*tmp1[56]+z[14]*tmp1[57]+z[67]*tmp1[58]+z[23]*tmp1[59]+z[76]*tmp1[60]
                    +z[32]*tmp1[61]+z[85]*tmp1[62]+z[41]*tmp1[63]+z[94]*tmp1[64]+z[50]*tmp1[65]+z[6]*tmp1[66]+z[59]*tmp1[67]+z[15]*tmp1[68]+z[68]*tmp1[69]+z[24]*tmp1[70]
                    +z[77]*tmp1[71]+z[33]*tmp1[72]+z[86]*tmp1[73]+z[42]*tmp1[74]+z[95]*tmp1[75]+z[51]*tmp1[76]+z[7]*tmp1[77]+z[60]*tmp1[78]+z[16]*tmp1[79]+z[69]*tmp1[80]
                    +z[25]*tmp1[81]+z[78]*tmp1[82]+z[34]*tmp1[83]+z[87]*tmp1[84]+z[43]*tmp1[85]+z[96]*tmp1[86]+z[52]*tmp1[87]+z[8]*tmp1[88]+z[61]*tmp1[89]+z[17]*tmp1[90]
                    +z[70]*tmp1[91]+z[26]*tmp1[92]+z[79]*tmp1[93]+z[35]*tmp1[94]+z[88]*tmp1[95]+z[44]*tmp1[96];
                    tab[nb_tmp3+54*nb3]=tab[nb_tmp3+54*nb3]+z[0]*tmp1[0]
                    +z[54]*tmp1[1]+z[11]*tmp1[2]+z[65]*tmp1[3]+z[22]*tmp1[4]+z[76]*tmp1[5]+z[33]*tmp1[6]+z[87]*tmp1[7]+z[44]*tmp1[8]+z[1]*tmp1[9]+z[55]*tmp1[10]
                    +z[12]*tmp1[11]+z[66]*tmp1[12]+z[23]*tmp1[13]+z[77]*tmp1[14]+z[34]*tmp1[15]+z[88]*tmp1[16]+z[45]*tmp1[17]+z[2]*tmp1[18]+z[56]*tmp1[19]+z[13]*tmp1[20]
                    +z[67]*tmp1[21]+z[24]*tmp1[22]+z[78]*tmp1[23]+z[35]*tmp1[24]+z[89]*tmp1[25]+z[46]*tmp1[26]+z[3]*tmp1[27]+z[57]*tmp1[28]+z[14]*tmp1[29]+z[68]*tmp1[30]
                    +z[25]*tmp1[31]+z[79]*tmp1[32]+z[36]*tmp1[33]+z[90]*tmp1[34]+z[47]*tmp1[35]+z[4]*tmp1[36]+z[58]*tmp1[37]+z[15]*tmp1[38]+z[69]*tmp1[39]+z[26]*tmp1[40]
                    +z[80]*tmp1[41]+z[37]*tmp1[42]+z[91]*tmp1[43]+z[48]*tmp1[44]+z[5]*tmp1[45]+z[59]*tmp1[46]+z[16]*tmp1[47]+z[70]*tmp1[48]+z[27]*tmp1[49]+z[81]*tmp1[50]
                    +z[38]*tmp1[51]+z[92]*tmp1[52]+z[49]*tmp1[53]+z[6]*tmp1[54]+z[60]*tmp1[55]+z[17]*tmp1[56]+z[71]*tmp1[57]+z[28]*tmp1[58]+z[82]*tmp1[59]+z[39]*tmp1[60]
                    +z[93]*tmp1[61]+z[50]*tmp1[62]+z[7]*tmp1[63]+z[61]*tmp1[64]+z[18]*tmp1[65]+z[72]*tmp1[66]+z[29]*tmp1[67]+z[83]*tmp1[68]+z[40]*tmp1[69]+z[94]*tmp1[70]
                    +z[51]*tmp1[71]+z[8]*tmp1[72]+z[62]*tmp1[73]+z[19]*tmp1[74]+z[73]*tmp1[75]+z[30]*tmp1[76]+z[84]*tmp1[77]+z[41]*tmp1[78]+z[95]*tmp1[79]+z[52]*tmp1[80]
                    +z[9]*tmp1[81]+z[63]*tmp1[82]+z[20]*tmp1[83]+z[74]*tmp1[84]+z[31]*tmp1[85]+z[85]*tmp1[86]+z[42]*tmp1[87]+z[96]*tmp1[88]+z[53]*tmp1[89]+z[10]*tmp1[90]
                    +z[64]*tmp1[91]+z[21]*tmp1[92]+z[75]*tmp1[93]+z[32]*tmp1[94]+z[86]*tmp1[95]+z[43]*tmp1[96];
                    tab[nb_tmp3+55*nb3]=tab[nb_tmp3+55*nb3]+z[0]*tmp1[0]
                    +z[55]*tmp1[1]+z[13]*tmp1[2]+z[68]*tmp1[3]+z[26]*tmp1[4]+z[81]*tmp1[5]+z[39]*tmp1[6]+z[94]*tmp1[7]+z[52]*tmp1[8]+z[10]*tmp1[9]+z[65]*tmp1[10]
                    +z[23]*tmp1[11]+z[78]*tmp1[12]+z[36]*tmp1[13]+z[91]*tmp1[14]+z[49]*tmp1[15]+z[7]*tmp1[16]+z[62]*tmp1[17]+z[20]*tmp1[18]+z[75]*tmp1[19]+z[33]*tmp1[20]
                    +z[88]*tmp1[21]+z[46]*tmp1[22]+z[4]*tmp1[23]+z[59]*tmp1[24]+z[17]*tmp1[25]+z[72]*tmp1[26]+z[30]*tmp1[27]+z[85]*tmp1[28]+z[43]*tmp1[29]+z[1]*tmp1[30]
                    +z[56]*tmp1[31]+z[14]*tmp1[32]+z[69]*tmp1[33]+z[27]*tmp1[34]+z[82]*tmp1[35]+z[40]*tmp1[36]+z[95]*tmp1[37]+z[53]*tmp1[38]+z[11]*tmp1[39]+z[66]*tmp1[40]
                    +z[24]*tmp1[41]+z[79]*tmp1[42]+z[37]*tmp1[43]+z[92]*tmp1[44]+z[50]*tmp1[45]+z[8]*tmp1[46]+z[63]*tmp1[47]+z[21]*tmp1[48]+z[76]*tmp1[49]+z[34]*tmp1[50]
                    +z[89]*tmp1[51]+z[47]*tmp1[52]+z[5]*tmp1[53]+z[60]*tmp1[54]+z[18]*tmp1[55]+z[73]*tmp1[56]+z[31]*tmp1[57]+z[86]*tmp1[58]+z[44]*tmp1[59]+z[2]*tmp1[60]
                    +z[57]*tmp1[61]+z[15]*tmp1[62]+z[70]*tmp1[63]+z[28]*tmp1[64]+z[83]*tmp1[65]+z[41]*tmp1[66]+z[96]*tmp1[67]+z[54]*tmp1[68]+z[12]*tmp1[69]+z[67]*tmp1[70]
                    +z[25]*tmp1[71]+z[80]*tmp1[72]+z[38]*tmp1[73]+z[93]*tmp1[74]+z[51]*tmp1[75]+z[9]*tmp1[76]+z[64]*tmp1[77]+z[22]*tmp1[78]+z[77]*tmp1[79]+z[35]*tmp1[80]
                    +z[90]*tmp1[81]+z[48]*tmp1[82]+z[6]*tmp1[83]+z[61]*tmp1[84]+z[19]*tmp1[85]+z[74]*tmp1[86]+z[32]*tmp1[87]+z[87]*tmp1[88]+z[45]*tmp1[89]+z[3]*tmp1[90]
                    +z[58]*tmp1[91]+z[16]*tmp1[92]+z[71]*tmp1[93]+z[29]*tmp1[94]+z[84]*tmp1[95]+z[42]*tmp1[96];
                    tab[nb_tmp3+56*nb3]=tab[nb_tmp3+56*nb3]+z[0]*tmp1[0]
                    +z[56]*tmp1[1]+z[15]*tmp1[2]+z[71]*tmp1[3]+z[30]*tmp1[4]+z[86]*tmp1[5]+z[45]*tmp1[6]+z[4]*tmp1[7]+z[60]*tmp1[8]+z[19]*tmp1[9]+z[75]*tmp1[10]
                    +z[34]*tmp1[11]+z[90]*tmp1[12]+z[49]*tmp1[13]+z[8]*tmp1[14]+z[64]*tmp1[15]+z[23]*tmp1[16]+z[79]*tmp1[17]+z[38]*tmp1[18]+z[94]*tmp1[19]+z[53]*tmp1[20]
                    +z[12]*tmp1[21]+z[68]*tmp1[22]+z[27]*tmp1[23]+z[83]*tmp1[24]+z[42]*tmp1[25]+z[1]*tmp1[26]+z[57]*tmp1[27]+z[16]*tmp1[28]+z[72]*tmp1[29]+z[31]*tmp1[30]
                    +z[87]*tmp1[31]+z[46]*tmp1[32]+z[5]*tmp1[33]+z[61]*tmp1[34]+z[20]*tmp1[35]+z[76]*tmp1[36]+z[35]*tmp1[37]+z[91]*tmp1[38]+z[50]*tmp1[39]+z[9]*tmp1[40]
                    +z[65]*tmp1[41]+z[24]*tmp1[42]+z[80]*tmp1[43]+z[39]*tmp1[44]+z[95]*tmp1[45]+z[54]*tmp1[46]+z[13]*tmp1[47]+z[69]*tmp1[48]+z[28]*tmp1[49]+z[84]*tmp1[50]
                    +z[43]*tmp1[51]+z[2]*tmp1[52]+z[58]*tmp1[53]+z[17]*tmp1[54]+z[73]*tmp1[55]+z[32]*tmp1[56]+z[88]*tmp1[57]+z[47]*tmp1[58]+z[6]*tmp1[59]+z[62]*tmp1[60]
                    +z[21]*tmp1[61]+z[77]*tmp1[62]+z[36]*tmp1[63]+z[92]*tmp1[64]+z[51]*tmp1[65]+z[10]*tmp1[66]+z[66]*tmp1[67]+z[25]*tmp1[68]+z[81]*tmp1[69]+z[40]*tmp1[70]
                    +z[96]*tmp1[71]+z[55]*tmp1[72]+z[14]*tmp1[73]+z[70]*tmp1[74]+z[29]*tmp1[75]+z[85]*tmp1[76]+z[44]*tmp1[77]+z[3]*tmp1[78]+z[59]*tmp1[79]+z[18]*tmp1[80]
                    +z[74]*tmp1[81]+z[33]*tmp1[82]+z[89]*tmp1[83]+z[48]*tmp1[84]+z[7]*tmp1[85]+z[63]*tmp1[86]+z[22]*tmp1[87]+z[78]*tmp1[88]+z[37]*tmp1[89]+z[93]*tmp1[90]
                    +z[52]*tmp1[91]+z[11]*tmp1[92]+z[67]*tmp1[93]+z[26]*tmp1[94]+z[82]*tmp1[95]+z[41]*tmp1[96];
                    tab[nb_tmp3+57*nb3]=tab[nb_tmp3+57*nb3]+z[0]*tmp1[0]
                    +z[57]*tmp1[1]+z[17]*tmp1[2]+z[74]*tmp1[3]+z[34]*tmp1[4]+z[91]*tmp1[5]+z[51]*tmp1[6]+z[11]*tmp1[7]+z[68]*tmp1[8]+z[28]*tmp1[9]+z[85]*tmp1[10]
                    +z[45]*tmp1[11]+z[5]*tmp1[12]+z[62]*tmp1[13]+z[22]*tmp1[14]+z[79]*tmp1[15]+z[39]*tmp1[16]+z[96]*tmp1[17]+z[56]*tmp1[18]+z[16]*tmp1[19]+z[73]*tmp1[20]
                    +z[33]*tmp1[21]+z[90]*tmp1[22]+z[50]*tmp1[23]+z[10]*tmp1[24]+z[67]*tmp1[25]+z[27]*tmp1[26]+z[84]*tmp1[27]+z[44]*tmp1[28]+z[4]*tmp1[29]+z[61]*tmp1[30]
                    +z[21]*tmp1[31]+z[78]*tmp1[32]+z[38]*tmp1[33]+z[95]*tmp1[34]+z[55]*tmp1[35]+z[15]*tmp1[36]+z[72]*tmp1[37]+z[32]*tmp1[38]+z[89]*tmp1[39]+z[49]*tmp1[40]
                    +z[9]*tmp1[41]+z[66]*tmp1[42]+z[26]*tmp1[43]+z[83]*tmp1[44]+z[43]*tmp1[45]+z[3]*tmp1[46]+z[60]*tmp1[47]+z[20]*tmp1[48]+z[77]*tmp1[49]+z[37]*tmp1[50]
                    +z[94]*tmp1[51]+z[54]*tmp1[52]+z[14]*tmp1[53]+z[71]*tmp1[54]+z[31]*tmp1[55]+z[88]*tmp1[56]+z[48]*tmp1[57]+z[8]*tmp1[58]+z[65]*tmp1[59]+z[25]*tmp1[60]
                    +z[82]*tmp1[61]+z[42]*tmp1[62]+z[2]*tmp1[63]+z[59]*tmp1[64]+z[19]*tmp1[65]+z[76]*tmp1[66]+z[36]*tmp1[67]+z[93]*tmp1[68]+z[53]*tmp1[69]+z[13]*tmp1[70]
                    +z[70]*tmp1[71]+z[30]*tmp1[72]+z[87]*tmp1[73]+z[47]*tmp1[74]+z[7]*tmp1[75]+z[64]*tmp1[76]+z[24]*tmp1[77]+z[81]*tmp1[78]+z[41]*tmp1[79]+z[1]*tmp1[80]
                    +z[58]*tmp1[81]+z[18]*tmp1[82]+z[75]*tmp1[83]+z[35]*tmp1[84]+z[92]*tmp1[85]+z[52]*tmp1[86]+z[12]*tmp1[87]+z[69]*tmp1[88]+z[29]*tmp1[89]+z[86]*tmp1[90]
                    +z[46]*tmp1[91]+z[6]*tmp1[92]+z[63]*tmp1[93]+z[23]*tmp1[94]+z[80]*tmp1[95]+z[40]*tmp1[96];
                    tab[nb_tmp3+58*nb3]=tab[nb_tmp3+58*nb3]+z[0]*tmp1[0]
                    +z[58]*tmp1[1]+z[19]*tmp1[2]+z[77]*tmp1[3]+z[38]*tmp1[4]+z[96]*tmp1[5]+z[57]*tmp1[6]+z[18]*tmp1[7]+z[76]*tmp1[8]+z[37]*tmp1[9]+z[95]*tmp1[10]
                    +z[56]*tmp1[11]+z[17]*tmp1[12]+z[75]*tmp1[13]+z[36]*tmp1[14]+z[94]*tmp1[15]+z[55]*tmp1[16]+z[16]*tmp1[17]+z[74]*tmp1[18]+z[35]*tmp1[19]+z[93]*tmp1[20]
                    +z[54]*tmp1[21]+z[15]*tmp1[22]+z[73]*tmp1[23]+z[34]*tmp1[24]+z[92]*tmp1[25]+z[53]*tmp1[26]+z[14]*tmp1[27]+z[72]*tmp1[28]+z[33]*tmp1[29]+z[91]*tmp1[30]
                    +z[52]*tmp1[31]+z[13]*tmp1[32]+z[71]*tmp1[33]+z[32]*tmp1[34]+z[90]*tmp1[35]+z[51]*tmp1[36]+z[12]*tmp1[37]+z[70]*tmp1[38]+z[31]*tmp1[39]+z[89]*tmp1[40]
                    +z[50]*tmp1[41]+z[11]*tmp1[42]+z[69]*tmp1[43]+z[30]*tmp1[44]+z[88]*tmp1[45]+z[49]*tmp1[46]+z[10]*tmp1[47]+z[68]*tmp1[48]+z[29]*tmp1[49]+z[87]*tmp1[50]
                    +z[48]*tmp1[51]+z[9]*tmp1[52]+z[67]*tmp1[53]+z[28]*tmp1[54]+z[86]*tmp1[55]+z[47]*tmp1[56]+z[8]*tmp1[57]+z[66]*tmp1[58]+z[27]*tmp1[59]+z[85]*tmp1[60]
                    +z[46]*tmp1[61]+z[7]*tmp1[62]+z[65]*tmp1[63]+z[26]*tmp1[64]+z[84]*tmp1[65]+z[45]*tmp1[66]+z[6]*tmp1[67]+z[64]*tmp1[68]+z[25]*tmp1[69]+z[83]*tmp1[70]
                    +z[44]*tmp1[71]+z[5]*tmp1[72]+z[63]*tmp1[73]+z[24]*tmp1[74]+z[82]*tmp1[75]+z[43]*tmp1[76]+z[4]*tmp1[77]+z[62]*tmp1[78]+z[23]*tmp1[79]+z[81]*tmp1[80]
                    +z[42]*tmp1[81]+z[3]*tmp1[82]+z[61]*tmp1[83]+z[22]*tmp1[84]+z[80]*tmp1[85]+z[41]*tmp1[86]+z[2]*tmp1[87]+z[60]*tmp1[88]+z[21]*tmp1[89]+z[79]*tmp1[90]
                    +z[40]*tmp1[91]+z[1]*tmp1[92]+z[59]*tmp1[93]+z[20]*tmp1[94]+z[78]*tmp1[95]+z[39]*tmp1[96];
                    tab[nb_tmp3+59*nb3]=tab[nb_tmp3+59*nb3]+z[0]*tmp1[0]
                    +z[59]*tmp1[1]+z[21]*tmp1[2]+z[80]*tmp1[3]+z[42]*tmp1[4]+z[4]*tmp1[5]+z[63]*tmp1[6]+z[25]*tmp1[7]+z[84]*tmp1[8]+z[46]*tmp1[9]+z[8]*tmp1[10]
                    +z[67]*tmp1[11]+z[29]*tmp1[12]+z[88]*tmp1[13]+z[50]*tmp1[14]+z[12]*tmp1[15]+z[71]*tmp1[16]+z[33]*tmp1[17]+z[92]*tmp1[18]+z[54]*tmp1[19]+z[16]*tmp1[20]
                    +z[75]*tmp1[21]+z[37]*tmp1[22]+z[96]*tmp1[23]+z[58]*tmp1[24]+z[20]*tmp1[25]+z[79]*tmp1[26]+z[41]*tmp1[27]+z[3]*tmp1[28]+z[62]*tmp1[29]+z[24]*tmp1[30]
                    +z[83]*tmp1[31]+z[45]*tmp1[32]+z[7]*tmp1[33]+z[66]*tmp1[34]+z[28]*tmp1[35]+z[87]*tmp1[36]+z[49]*tmp1[37]+z[11]*tmp1[38]+z[70]*tmp1[39]+z[32]*tmp1[40]
                    +z[91]*tmp1[41]+z[53]*tmp1[42]+z[15]*tmp1[43]+z[74]*tmp1[44]+z[36]*tmp1[45]+z[95]*tmp1[46]+z[57]*tmp1[47]+z[19]*tmp1[48]+z[78]*tmp1[49]+z[40]*tmp1[50]
                    +z[2]*tmp1[51]+z[61]*tmp1[52]+z[23]*tmp1[53]+z[82]*tmp1[54]+z[44]*tmp1[55]+z[6]*tmp1[56]+z[65]*tmp1[57]+z[27]*tmp1[58]+z[86]*tmp1[59]+z[48]*tmp1[60]
                    +z[10]*tmp1[61]+z[69]*tmp1[62]+z[31]*tmp1[63]+z[90]*tmp1[64]+z[52]*tmp1[65]+z[14]*tmp1[66]+z[73]*tmp1[67]+z[35]*tmp1[68]+z[94]*tmp1[69]+z[56]*tmp1[70]
                    +z[18]*tmp1[71]+z[77]*tmp1[72]+z[39]*tmp1[73]+z[1]*tmp1[74]+z[60]*tmp1[75]+z[22]*tmp1[76]+z[81]*tmp1[77]+z[43]*tmp1[78]+z[5]*tmp1[79]+z[64]*tmp1[80]
                    +z[26]*tmp1[81]+z[85]*tmp1[82]+z[47]*tmp1[83]+z[9]*tmp1[84]+z[68]*tmp1[85]+z[30]*tmp1[86]+z[89]*tmp1[87]+z[51]*tmp1[88]+z[13]*tmp1[89]+z[72]*tmp1[90]
                    +z[34]*tmp1[91]+z[93]*tmp1[92]+z[55]*tmp1[93]+z[17]*tmp1[94]+z[76]*tmp1[95]+z[38]*tmp1[96];
                    tab[nb_tmp3+60*nb3]=tab[nb_tmp3+60*nb3]+z[0]*tmp1[0]
                    +z[60]*tmp1[1]+z[23]*tmp1[2]+z[83]*tmp1[3]+z[46]*tmp1[4]+z[9]*tmp1[5]+z[69]*tmp1[6]+z[32]*tmp1[7]+z[92]*tmp1[8]+z[55]*tmp1[9]+z[18]*tmp1[10]
                    +z[78]*tmp1[11]+z[41]*tmp1[12]+z[4]*tmp1[13]+z[64]*tmp1[14]+z[27]*tmp1[15]+z[87]*tmp1[16]+z[50]*tmp1[17]+z[13]*tmp1[18]+z[73]*tmp1[19]+z[36]*tmp1[20]
                    +z[96]*tmp1[21]+z[59]*tmp1[22]+z[22]*tmp1[23]+z[82]*tmp1[24]+z[45]*tmp1[25]+z[8]*tmp1[26]+z[68]*tmp1[27]+z[31]*tmp1[28]+z[91]*tmp1[29]+z[54]*tmp1[30]
                    +z[17]*tmp1[31]+z[77]*tmp1[32]+z[40]*tmp1[33]+z[3]*tmp1[34]+z[63]*tmp1[35]+z[26]*tmp1[36]+z[86]*tmp1[37]+z[49]*tmp1[38]+z[12]*tmp1[39]+z[72]*tmp1[40]
                    +z[35]*tmp1[41]+z[95]*tmp1[42]+z[58]*tmp1[43]+z[21]*tmp1[44]+z[81]*tmp1[45]+z[44]*tmp1[46]+z[7]*tmp1[47]+z[67]*tmp1[48]+z[30]*tmp1[49]+z[90]*tmp1[50]
                    +z[53]*tmp1[51]+z[16]*tmp1[52]+z[76]*tmp1[53]+z[39]*tmp1[54]+z[2]*tmp1[55]+z[62]*tmp1[56]+z[25]*tmp1[57]+z[85]*tmp1[58]+z[48]*tmp1[59]+z[11]*tmp1[60]
                    +z[71]*tmp1[61]+z[34]*tmp1[62]+z[94]*tmp1[63]+z[57]*tmp1[64]+z[20]*tmp1[65]+z[80]*tmp1[66]+z[43]*tmp1[67]+z[6]*tmp1[68]+z[66]*tmp1[69]+z[29]*tmp1[70]
                    +z[89]*tmp1[71]+z[52]*tmp1[72]+z[15]*tmp1[73]+z[75]*tmp1[74]+z[38]*tmp1[75]+z[1]*tmp1[76]+z[61]*tmp1[77]+z[24]*tmp1[78]+z[84]*tmp1[79]+z[47]*tmp1[80]
                    +z[10]*tmp1[81]+z[70]*tmp1[82]+z[33]*tmp1[83]+z[93]*tmp1[84]+z[56]*tmp1[85]+z[19]*tmp1[86]+z[79]*tmp1[87]+z[42]*tmp1[88]+z[5]*tmp1[89]+z[65]*tmp1[90]
                    +z[28]*tmp1[91]+z[88]*tmp1[92]+z[51]*tmp1[93]+z[14]*tmp1[94]+z[74]*tmp1[95]+z[37]*tmp1[96];
                    tab[nb_tmp3+61*nb3]=tab[nb_tmp3+61*nb3]+z[0]*tmp1[0]
                    +z[61]*tmp1[1]+z[25]*tmp1[2]+z[86]*tmp1[3]+z[50]*tmp1[4]+z[14]*tmp1[5]+z[75]*tmp1[6]+z[39]*tmp1[7]+z[3]*tmp1[8]+z[64]*tmp1[9]+z[28]*tmp1[10]
                    +z[89]*tmp1[11]+z[53]*tmp1[12]+z[17]*tmp1[13]+z[78]*tmp1[14]+z[42]*tmp1[15]+z[6]*tmp1[16]+z[67]*tmp1[17]+z[31]*tmp1[18]+z[92]*tmp1[19]+z[56]*tmp1[20]
                    +z[20]*tmp1[21]+z[81]*tmp1[22]+z[45]*tmp1[23]+z[9]*tmp1[24]+z[70]*tmp1[25]+z[34]*tmp1[26]+z[95]*tmp1[27]+z[59]*tmp1[28]+z[23]*tmp1[29]+z[84]*tmp1[30]
                    +z[48]*tmp1[31]+z[12]*tmp1[32]+z[73]*tmp1[33]+z[37]*tmp1[34]+z[1]*tmp1[35]+z[62]*tmp1[36]+z[26]*tmp1[37]+z[87]*tmp1[38]+z[51]*tmp1[39]+z[15]*tmp1[40]
                    +z[76]*tmp1[41]+z[40]*tmp1[42]+z[4]*tmp1[43]+z[65]*tmp1[44]+z[29]*tmp1[45]+z[90]*tmp1[46]+z[54]*tmp1[47]+z[18]*tmp1[48]+z[79]*tmp1[49]+z[43]*tmp1[50]
                    +z[7]*tmp1[51]+z[68]*tmp1[52]+z[32]*tmp1[53]+z[93]*tmp1[54]+z[57]*tmp1[55]+z[21]*tmp1[56]+z[82]*tmp1[57]+z[46]*tmp1[58]+z[10]*tmp1[59]+z[71]*tmp1[60]
                    +z[35]*tmp1[61]+z[96]*tmp1[62]+z[60]*tmp1[63]+z[24]*tmp1[64]+z[85]*tmp1[65]+z[49]*tmp1[66]+z[13]*tmp1[67]+z[74]*tmp1[68]+z[38]*tmp1[69]+z[2]*tmp1[70]
                    +z[63]*tmp1[71]+z[27]*tmp1[72]+z[88]*tmp1[73]+z[52]*tmp1[74]+z[16]*tmp1[75]+z[77]*tmp1[76]+z[41]*tmp1[77]+z[5]*tmp1[78]+z[66]*tmp1[79]+z[30]*tmp1[80]
                    +z[91]*tmp1[81]+z[55]*tmp1[82]+z[19]*tmp1[83]+z[80]*tmp1[84]+z[44]*tmp1[85]+z[8]*tmp1[86]+z[69]*tmp1[87]+z[33]*tmp1[88]+z[94]*tmp1[89]+z[58]*tmp1[90]
                    +z[22]*tmp1[91]+z[83]*tmp1[92]+z[47]*tmp1[93]+z[11]*tmp1[94]+z[72]*tmp1[95]+z[36]*tmp1[96];
                    tab[nb_tmp3+62*nb3]=tab[nb_tmp3+62*nb3]+z[0]*tmp1[0]
                    +z[62]*tmp1[1]+z[27]*tmp1[2]+z[89]*tmp1[3]+z[54]*tmp1[4]+z[19]*tmp1[5]+z[81]*tmp1[6]+z[46]*tmp1[7]+z[11]*tmp1[8]+z[73]*tmp1[9]+z[38]*tmp1[10]
                    +z[3]*tmp1[11]+z[65]*tmp1[12]+z[30]*tmp1[13]+z[92]*tmp1[14]+z[57]*tmp1[15]+z[22]*tmp1[16]+z[84]*tmp1[17]+z[49]*tmp1[18]+z[14]*tmp1[19]+z[76]*tmp1[20]
                    +z[41]*tmp1[21]+z[6]*tmp1[22]+z[68]*tmp1[23]+z[33]*tmp1[24]+z[95]*tmp1[25]+z[60]*tmp1[26]+z[25]*tmp1[27]+z[87]*tmp1[28]+z[52]*tmp1[29]+z[17]*tmp1[30]
                    +z[79]*tmp1[31]+z[44]*tmp1[32]+z[9]*tmp1[33]+z[71]*tmp1[34]+z[36]*tmp1[35]+z[1]*tmp1[36]+z[63]*tmp1[37]+z[28]*tmp1[38]+z[90]*tmp1[39]+z[55]*tmp1[40]
                    +z[20]*tmp1[41]+z[82]*tmp1[42]+z[47]*tmp1[43]+z[12]*tmp1[44]+z[74]*tmp1[45]+z[39]*tmp1[46]+z[4]*tmp1[47]+z[66]*tmp1[48]+z[31]*tmp1[49]+z[93]*tmp1[50]
                    +z[58]*tmp1[51]+z[23]*tmp1[52]+z[85]*tmp1[53]+z[50]*tmp1[54]+z[15]*tmp1[55]+z[77]*tmp1[56]+z[42]*tmp1[57]+z[7]*tmp1[58]+z[69]*tmp1[59]+z[34]*tmp1[60]
                    +z[96]*tmp1[61]+z[61]*tmp1[62]+z[26]*tmp1[63]+z[88]*tmp1[64]+z[53]*tmp1[65]+z[18]*tmp1[66]+z[80]*tmp1[67]+z[45]*tmp1[68]+z[10]*tmp1[69]+z[72]*tmp1[70]
                    +z[37]*tmp1[71]+z[2]*tmp1[72]+z[64]*tmp1[73]+z[29]*tmp1[74]+z[91]*tmp1[75]+z[56]*tmp1[76]+z[21]*tmp1[77]+z[83]*tmp1[78]+z[48]*tmp1[79]+z[13]*tmp1[80]
                    +z[75]*tmp1[81]+z[40]*tmp1[82]+z[5]*tmp1[83]+z[67]*tmp1[84]+z[32]*tmp1[85]+z[94]*tmp1[86]+z[59]*tmp1[87]+z[24]*tmp1[88]+z[86]*tmp1[89]+z[51]*tmp1[90]
                    +z[16]*tmp1[91]+z[78]*tmp1[92]+z[43]*tmp1[93]+z[8]*tmp1[94]+z[70]*tmp1[95]+z[35]*tmp1[96];
                    tab[nb_tmp3+63*nb3]=tab[nb_tmp3+63*nb3]+z[0]*tmp1[0]
                    +z[63]*tmp1[1]+z[29]*tmp1[2]+z[92]*tmp1[3]+z[58]*tmp1[4]+z[24]*tmp1[5]+z[87]*tmp1[6]+z[53]*tmp1[7]+z[19]*tmp1[8]+z[82]*tmp1[9]+z[48]*tmp1[10]
                    +z[14]*tmp1[11]+z[77]*tmp1[12]+z[43]*tmp1[13]+z[9]*tmp1[14]+z[72]*tmp1[15]+z[38]*tmp1[16]+z[4]*tmp1[17]+z[67]*tmp1[18]+z[33]*tmp1[19]+z[96]*tmp1[20]
                    +z[62]*tmp1[21]+z[28]*tmp1[22]+z[91]*tmp1[23]+z[57]*tmp1[24]+z[23]*tmp1[25]+z[86]*tmp1[26]+z[52]*tmp1[27]+z[18]*tmp1[28]+z[81]*tmp1[29]+z[47]*tmp1[30]
                    +z[13]*tmp1[31]+z[76]*tmp1[32]+z[42]*tmp1[33]+z[8]*tmp1[34]+z[71]*tmp1[35]+z[37]*tmp1[36]+z[3]*tmp1[37]+z[66]*tmp1[38]+z[32]*tmp1[39]+z[95]*tmp1[40]
                    +z[61]*tmp1[41]+z[27]*tmp1[42]+z[90]*tmp1[43]+z[56]*tmp1[44]+z[22]*tmp1[45]+z[85]*tmp1[46]+z[51]*tmp1[47]+z[17]*tmp1[48]+z[80]*tmp1[49]+z[46]*tmp1[50]
                    +z[12]*tmp1[51]+z[75]*tmp1[52]+z[41]*tmp1[53]+z[7]*tmp1[54]+z[70]*tmp1[55]+z[36]*tmp1[56]+z[2]*tmp1[57]+z[65]*tmp1[58]+z[31]*tmp1[59]+z[94]*tmp1[60]
                    +z[60]*tmp1[61]+z[26]*tmp1[62]+z[89]*tmp1[63]+z[55]*tmp1[64]+z[21]*tmp1[65]+z[84]*tmp1[66]+z[50]*tmp1[67]+z[16]*tmp1[68]+z[79]*tmp1[69]+z[45]*tmp1[70]
                    +z[11]*tmp1[71]+z[74]*tmp1[72]+z[40]*tmp1[73]+z[6]*tmp1[74]+z[69]*tmp1[75]+z[35]*tmp1[76]+z[1]*tmp1[77]+z[64]*tmp1[78]+z[30]*tmp1[79]+z[93]*tmp1[80]
                    +z[59]*tmp1[81]+z[25]*tmp1[82]+z[88]*tmp1[83]+z[54]*tmp1[84]+z[20]*tmp1[85]+z[83]*tmp1[86]+z[49]*tmp1[87]+z[15]*tmp1[88]+z[78]*tmp1[89]+z[44]*tmp1[90]
                    +z[10]*tmp1[91]+z[73]*tmp1[92]+z[39]*tmp1[93]+z[5]*tmp1[94]+z[68]*tmp1[95]+z[34]*tmp1[96];
                    tab[nb_tmp3+64*nb3]=tab[nb_tmp3+64*nb3]+z[0]*tmp1[0]
                    +z[64]*tmp1[1]+z[31]*tmp1[2]+z[95]*tmp1[3]+z[62]*tmp1[4]+z[29]*tmp1[5]+z[93]*tmp1[6]+z[60]*tmp1[7]+z[27]*tmp1[8]+z[91]*tmp1[9]+z[58]*tmp1[10]
                    +z[25]*tmp1[11]+z[89]*tmp1[12]+z[56]*tmp1[13]+z[23]*tmp1[14]+z[87]*tmp1[15]+z[54]*tmp1[16]+z[21]*tmp1[17]+z[85]*tmp1[18]+z[52]*tmp1[19]+z[19]*tmp1[20]
                    +z[83]*tmp1[21]+z[50]*tmp1[22]+z[17]*tmp1[23]+z[81]*tmp1[24]+z[48]*tmp1[25]+z[15]*tmp1[26]+z[79]*tmp1[27]+z[46]*tmp1[28]+z[13]*tmp1[29]+z[77]*tmp1[30]
                    +z[44]*tmp1[31]+z[11]*tmp1[32]+z[75]*tmp1[33]+z[42]*tmp1[34]+z[9]*tmp1[35]+z[73]*tmp1[36]+z[40]*tmp1[37]+z[7]*tmp1[38]+z[71]*tmp1[39]+z[38]*tmp1[40]
                    +z[5]*tmp1[41]+z[69]*tmp1[42]+z[36]*tmp1[43]+z[3]*tmp1[44]+z[67]*tmp1[45]+z[34]*tmp1[46]+z[1]*tmp1[47]+z[65]*tmp1[48]+z[32]*tmp1[49]+z[96]*tmp1[50]
                    +z[63]*tmp1[51]+z[30]*tmp1[52]+z[94]*tmp1[53]+z[61]*tmp1[54]+z[28]*tmp1[55]+z[92]*tmp1[56]+z[59]*tmp1[57]+z[26]*tmp1[58]+z[90]*tmp1[59]+z[57]*tmp1[60]
                    +z[24]*tmp1[61]+z[88]*tmp1[62]+z[55]*tmp1[63]+z[22]*tmp1[64]+z[86]*tmp1[65]+z[53]*tmp1[66]+z[20]*tmp1[67]+z[84]*tmp1[68]+z[51]*tmp1[69]+z[18]*tmp1[70]
                    +z[82]*tmp1[71]+z[49]*tmp1[72]+z[16]*tmp1[73]+z[80]*tmp1[74]+z[47]*tmp1[75]+z[14]*tmp1[76]+z[78]*tmp1[77]+z[45]*tmp1[78]+z[12]*tmp1[79]+z[76]*tmp1[80]
                    +z[43]*tmp1[81]+z[10]*tmp1[82]+z[74]*tmp1[83]+z[41]*tmp1[84]+z[8]*tmp1[85]+z[72]*tmp1[86]+z[39]*tmp1[87]+z[6]*tmp1[88]+z[70]*tmp1[89]+z[37]*tmp1[90]
                    +z[4]*tmp1[91]+z[68]*tmp1[92]+z[35]*tmp1[93]+z[2]*tmp1[94]+z[66]*tmp1[95]+z[33]*tmp1[96];
                    tab[nb_tmp3+65*nb3]=tab[nb_tmp3+65*nb3]+z[0]*tmp1[0]
                    +z[65]*tmp1[1]+z[33]*tmp1[2]+z[1]*tmp1[3]+z[66]*tmp1[4]+z[34]*tmp1[5]+z[2]*tmp1[6]+z[67]*tmp1[7]+z[35]*tmp1[8]+z[3]*tmp1[9]+z[68]*tmp1[10]
                    +z[36]*tmp1[11]+z[4]*tmp1[12]+z[69]*tmp1[13]+z[37]*tmp1[14]+z[5]*tmp1[15]+z[70]*tmp1[16]+z[38]*tmp1[17]+z[6]*tmp1[18]+z[71]*tmp1[19]+z[39]*tmp1[20]
                    +z[7]*tmp1[21]+z[72]*tmp1[22]+z[40]*tmp1[23]+z[8]*tmp1[24]+z[73]*tmp1[25]+z[41]*tmp1[26]+z[9]*tmp1[27]+z[74]*tmp1[28]+z[42]*tmp1[29]+z[10]*tmp1[30]
                    +z[75]*tmp1[31]+z[43]*tmp1[32]+z[11]*tmp1[33]+z[76]*tmp1[34]+z[44]*tmp1[35]+z[12]*tmp1[36]+z[77]*tmp1[37]+z[45]*tmp1[38]+z[13]*tmp1[39]+z[78]*tmp1[40]
                    +z[46]*tmp1[41]+z[14]*tmp1[42]+z[79]*tmp1[43]+z[47]*tmp1[44]+z[15]*tmp1[45]+z[80]*tmp1[46]+z[48]*tmp1[47]+z[16]*tmp1[48]+z[81]*tmp1[49]+z[49]*tmp1[50]
                    +z[17]*tmp1[51]+z[82]*tmp1[52]+z[50]*tmp1[53]+z[18]*tmp1[54]+z[83]*tmp1[55]+z[51]*tmp1[56]+z[19]*tmp1[57]+z[84]*tmp1[58]+z[52]*tmp1[59]+z[20]*tmp1[60]
                    +z[85]*tmp1[61]+z[53]*tmp1[62]+z[21]*tmp1[63]+z[86]*tmp1[64]+z[54]*tmp1[65]+z[22]*tmp1[66]+z[87]*tmp1[67]+z[55]*tmp1[68]+z[23]*tmp1[69]+z[88]*tmp1[70]
                    +z[56]*tmp1[71]+z[24]*tmp1[72]+z[89]*tmp1[73]+z[57]*tmp1[74]+z[25]*tmp1[75]+z[90]*tmp1[76]+z[58]*tmp1[77]+z[26]*tmp1[78]+z[91]*tmp1[79]+z[59]*tmp1[80]
                    +z[27]*tmp1[81]+z[92]*tmp1[82]+z[60]*tmp1[83]+z[28]*tmp1[84]+z[93]*tmp1[85]+z[61]*tmp1[86]+z[29]*tmp1[87]+z[94]*tmp1[88]+z[62]*tmp1[89]+z[30]*tmp1[90]
                    +z[95]*tmp1[91]+z[63]*tmp1[92]+z[31]*tmp1[93]+z[96]*tmp1[94]+z[64]*tmp1[95]+z[32]*tmp1[96];
                    tab[nb_tmp3+66*nb3]=tab[nb_tmp3+66*nb3]+z[0]*tmp1[0]
                    +z[66]*tmp1[1]+z[35]*tmp1[2]+z[4]*tmp1[3]+z[70]*tmp1[4]+z[39]*tmp1[5]+z[8]*tmp1[6]+z[74]*tmp1[7]+z[43]*tmp1[8]+z[12]*tmp1[9]+z[78]*tmp1[10]
                    +z[47]*tmp1[11]+z[16]*tmp1[12]+z[82]*tmp1[13]+z[51]*tmp1[14]+z[20]*tmp1[15]+z[86]*tmp1[16]+z[55]*tmp1[17]+z[24]*tmp1[18]+z[90]*tmp1[19]+z[59]*tmp1[20]
                    +z[28]*tmp1[21]+z[94]*tmp1[22]+z[63]*tmp1[23]+z[32]*tmp1[24]+z[1]*tmp1[25]+z[67]*tmp1[26]+z[36]*tmp1[27]+z[5]*tmp1[28]+z[71]*tmp1[29]+z[40]*tmp1[30]
                    +z[9]*tmp1[31]+z[75]*tmp1[32]+z[44]*tmp1[33]+z[13]*tmp1[34]+z[79]*tmp1[35]+z[48]*tmp1[36]+z[17]*tmp1[37]+z[83]*tmp1[38]+z[52]*tmp1[39]+z[21]*tmp1[40]
                    +z[87]*tmp1[41]+z[56]*tmp1[42]+z[25]*tmp1[43]+z[91]*tmp1[44]+z[60]*tmp1[45]+z[29]*tmp1[46]+z[95]*tmp1[47]+z[64]*tmp1[48]+z[33]*tmp1[49]+z[2]*tmp1[50]
                    +z[68]*tmp1[51]+z[37]*tmp1[52]+z[6]*tmp1[53]+z[72]*tmp1[54]+z[41]*tmp1[55]+z[10]*tmp1[56]+z[76]*tmp1[57]+z[45]*tmp1[58]+z[14]*tmp1[59]+z[80]*tmp1[60]
                    +z[49]*tmp1[61]+z[18]*tmp1[62]+z[84]*tmp1[63]+z[53]*tmp1[64]+z[22]*tmp1[65]+z[88]*tmp1[66]+z[57]*tmp1[67]+z[26]*tmp1[68]+z[92]*tmp1[69]+z[61]*tmp1[70]
                    +z[30]*tmp1[71]+z[96]*tmp1[72]+z[65]*tmp1[73]+z[34]*tmp1[74]+z[3]*tmp1[75]+z[69]*tmp1[76]+z[38]*tmp1[77]+z[7]*tmp1[78]+z[73]*tmp1[79]+z[42]*tmp1[80]
                    +z[11]*tmp1[81]+z[77]*tmp1[82]+z[46]*tmp1[83]+z[15]*tmp1[84]+z[81]*tmp1[85]+z[50]*tmp1[86]+z[19]*tmp1[87]+z[85]*tmp1[88]+z[54]*tmp1[89]+z[23]*tmp1[90]
                    +z[89]*tmp1[91]+z[58]*tmp1[92]+z[27]*tmp1[93]+z[93]*tmp1[94]+z[62]*tmp1[95]+z[31]*tmp1[96];
                    tab[nb_tmp3+67*nb3]=tab[nb_tmp3+67*nb3]+z[0]*tmp1[0]
                    +z[67]*tmp1[1]+z[37]*tmp1[2]+z[7]*tmp1[3]+z[74]*tmp1[4]+z[44]*tmp1[5]+z[14]*tmp1[6]+z[81]*tmp1[7]+z[51]*tmp1[8]+z[21]*tmp1[9]+z[88]*tmp1[10]
                    +z[58]*tmp1[11]+z[28]*tmp1[12]+z[95]*tmp1[13]+z[65]*tmp1[14]+z[35]*tmp1[15]+z[5]*tmp1[16]+z[72]*tmp1[17]+z[42]*tmp1[18]+z[12]*tmp1[19]+z[79]*tmp1[20]
                    +z[49]*tmp1[21]+z[19]*tmp1[22]+z[86]*tmp1[23]+z[56]*tmp1[24]+z[26]*tmp1[25]+z[93]*tmp1[26]+z[63]*tmp1[27]+z[33]*tmp1[28]+z[3]*tmp1[29]+z[70]*tmp1[30]
                    +z[40]*tmp1[31]+z[10]*tmp1[32]+z[77]*tmp1[33]+z[47]*tmp1[34]+z[17]*tmp1[35]+z[84]*tmp1[36]+z[54]*tmp1[37]+z[24]*tmp1[38]+z[91]*tmp1[39]+z[61]*tmp1[40]
                    +z[31]*tmp1[41]+z[1]*tmp1[42]+z[68]*tmp1[43]+z[38]*tmp1[44]+z[8]*tmp1[45]+z[75]*tmp1[46]+z[45]*tmp1[47]+z[15]*tmp1[48]+z[82]*tmp1[49]+z[52]*tmp1[50]
                    +z[22]*tmp1[51]+z[89]*tmp1[52]+z[59]*tmp1[53]+z[29]*tmp1[54]+z[96]*tmp1[55]+z[66]*tmp1[56]+z[36]*tmp1[57]+z[6]*tmp1[58]+z[73]*tmp1[59]+z[43]*tmp1[60]
                    +z[13]*tmp1[61]+z[80]*tmp1[62]+z[50]*tmp1[63]+z[20]*tmp1[64]+z[87]*tmp1[65]+z[57]*tmp1[66]+z[27]*tmp1[67]+z[94]*tmp1[68]+z[64]*tmp1[69]+z[34]*tmp1[70]
                    +z[4]*tmp1[71]+z[71]*tmp1[72]+z[41]*tmp1[73]+z[11]*tmp1[74]+z[78]*tmp1[75]+z[48]*tmp1[76]+z[18]*tmp1[77]+z[85]*tmp1[78]+z[55]*tmp1[79]+z[25]*tmp1[80]
                    +z[92]*tmp1[81]+z[62]*tmp1[82]+z[32]*tmp1[83]+z[2]*tmp1[84]+z[69]*tmp1[85]+z[39]*tmp1[86]+z[9]*tmp1[87]+z[76]*tmp1[88]+z[46]*tmp1[89]+z[16]*tmp1[90]
                    +z[83]*tmp1[91]+z[53]*tmp1[92]+z[23]*tmp1[93]+z[90]*tmp1[94]+z[60]*tmp1[95]+z[30]*tmp1[96];
                    tab[nb_tmp3+68*nb3]=tab[nb_tmp3+68*nb3]+z[0]*tmp1[0]
                    +z[68]*tmp1[1]+z[39]*tmp1[2]+z[10]*tmp1[3]+z[78]*tmp1[4]+z[49]*tmp1[5]+z[20]*tmp1[6]+z[88]*tmp1[7]+z[59]*tmp1[8]+z[30]*tmp1[9]+z[1]*tmp1[10]
                    +z[69]*tmp1[11]+z[40]*tmp1[12]+z[11]*tmp1[13]+z[79]*tmp1[14]+z[50]*tmp1[15]+z[21]*tmp1[16]+z[89]*tmp1[17]+z[60]*tmp1[18]+z[31]*tmp1[19]+z[2]*tmp1[20]
                    +z[70]*tmp1[21]+z[41]*tmp1[22]+z[12]*tmp1[23]+z[80]*tmp1[24]+z[51]*tmp1[25]+z[22]*tmp1[26]+z[90]*tmp1[27]+z[61]*tmp1[28]+z[32]*tmp1[29]+z[3]*tmp1[30]
                    +z[71]*tmp1[31]+z[42]*tmp1[32]+z[13]*tmp1[33]+z[81]*tmp1[34]+z[52]*tmp1[35]+z[23]*tmp1[36]+z[91]*tmp1[37]+z[62]*tmp1[38]+z[33]*tmp1[39]+z[4]*tmp1[40]
                    +z[72]*tmp1[41]+z[43]*tmp1[42]+z[14]*tmp1[43]+z[82]*tmp1[44]+z[53]*tmp1[45]+z[24]*tmp1[46]+z[92]*tmp1[47]+z[63]*tmp1[48]+z[34]*tmp1[49]+z[5]*tmp1[50]
                    +z[73]*tmp1[51]+z[44]*tmp1[52]+z[15]*tmp1[53]+z[83]*tmp1[54]+z[54]*tmp1[55]+z[25]*tmp1[56]+z[93]*tmp1[57]+z[64]*tmp1[58]+z[35]*tmp1[59]+z[6]*tmp1[60]
                    +z[74]*tmp1[61]+z[45]*tmp1[62]+z[16]*tmp1[63]+z[84]*tmp1[64]+z[55]*tmp1[65]+z[26]*tmp1[66]+z[94]*tmp1[67]+z[65]*tmp1[68]+z[36]*tmp1[69]+z[7]*tmp1[70]
                    +z[75]*tmp1[71]+z[46]*tmp1[72]+z[17]*tmp1[73]+z[85]*tmp1[74]+z[56]*tmp1[75]+z[27]*tmp1[76]+z[95]*tmp1[77]+z[66]*tmp1[78]+z[37]*tmp1[79]+z[8]*tmp1[80]
                    +z[76]*tmp1[81]+z[47]*tmp1[82]+z[18]*tmp1[83]+z[86]*tmp1[84]+z[57]*tmp1[85]+z[28]*tmp1[86]+z[96]*tmp1[87]+z[67]*tmp1[88]+z[38]*tmp1[89]+z[9]*tmp1[90]
                    +z[77]*tmp1[91]+z[48]*tmp1[92]+z[19]*tmp1[93]+z[87]*tmp1[94]+z[58]*tmp1[95]+z[29]*tmp1[96];
                    tab[nb_tmp3+69*nb3]=tab[nb_tmp3+69*nb3]+z[0]*tmp1[0]
                    +z[69]*tmp1[1]+z[41]*tmp1[2]+z[13]*tmp1[3]+z[82]*tmp1[4]+z[54]*tmp1[5]+z[26]*tmp1[6]+z[95]*tmp1[7]+z[67]*tmp1[8]+z[39]*tmp1[9]+z[11]*tmp1[10]
                    +z[80]*tmp1[11]+z[52]*tmp1[12]+z[24]*tmp1[13]+z[93]*tmp1[14]+z[65]*tmp1[15]+z[37]*tmp1[16]+z[9]*tmp1[17]+z[78]*tmp1[18]+z[50]*tmp1[19]+z[22]*tmp1[20]
                    +z[91]*tmp1[21]+z[63]*tmp1[22]+z[35]*tmp1[23]+z[7]*tmp1[24]+z[76]*tmp1[25]+z[48]*tmp1[26]+z[20]*tmp1[27]+z[89]*tmp1[28]+z[61]*tmp1[29]+z[33]*tmp1[30]
                    +z[5]*tmp1[31]+z[74]*tmp1[32]+z[46]*tmp1[33]+z[18]*tmp1[34]+z[87]*tmp1[35]+z[59]*tmp1[36]+z[31]*tmp1[37]+z[3]*tmp1[38]+z[72]*tmp1[39]+z[44]*tmp1[40]
                    +z[16]*tmp1[41]+z[85]*tmp1[42]+z[57]*tmp1[43]+z[29]*tmp1[44]+z[1]*tmp1[45]+z[70]*tmp1[46]+z[42]*tmp1[47]+z[14]*tmp1[48]+z[83]*tmp1[49]+z[55]*tmp1[50]
                    +z[27]*tmp1[51]+z[96]*tmp1[52]+z[68]*tmp1[53]+z[40]*tmp1[54]+z[12]*tmp1[55]+z[81]*tmp1[56]+z[53]*tmp1[57]+z[25]*tmp1[58]+z[94]*tmp1[59]+z[66]*tmp1[60]
                    +z[38]*tmp1[61]+z[10]*tmp1[62]+z[79]*tmp1[63]+z[51]*tmp1[64]+z[23]*tmp1[65]+z[92]*tmp1[66]+z[64]*tmp1[67]+z[36]*tmp1[68]+z[8]*tmp1[69]+z[77]*tmp1[70]
                    +z[49]*tmp1[71]+z[21]*tmp1[72]+z[90]*tmp1[73]+z[62]*tmp1[74]+z[34]*tmp1[75]+z[6]*tmp1[76]+z[75]*tmp1[77]+z[47]*tmp1[78]+z[19]*tmp1[79]+z[88]*tmp1[80]
                    +z[60]*tmp1[81]+z[32]*tmp1[82]+z[4]*tmp1[83]+z[73]*tmp1[84]+z[45]*tmp1[85]+z[17]*tmp1[86]+z[86]*tmp1[87]+z[58]*tmp1[88]+z[30]*tmp1[89]+z[2]*tmp1[90]
                    +z[71]*tmp1[91]+z[43]*tmp1[92]+z[15]*tmp1[93]+z[84]*tmp1[94]+z[56]*tmp1[95]+z[28]*tmp1[96];
                    tab[nb_tmp3+70*nb3]=tab[nb_tmp3+70*nb3]+z[0]*tmp1[0]
                    +z[70]*tmp1[1]+z[43]*tmp1[2]+z[16]*tmp1[3]+z[86]*tmp1[4]+z[59]*tmp1[5]+z[32]*tmp1[6]+z[5]*tmp1[7]+z[75]*tmp1[8]+z[48]*tmp1[9]+z[21]*tmp1[10]
                    +z[91]*tmp1[11]+z[64]*tmp1[12]+z[37]*tmp1[13]+z[10]*tmp1[14]+z[80]*tmp1[15]+z[53]*tmp1[16]+z[26]*tmp1[17]+z[96]*tmp1[18]+z[69]*tmp1[19]+z[42]*tmp1[20]
                    +z[15]*tmp1[21]+z[85]*tmp1[22]+z[58]*tmp1[23]+z[31]*tmp1[24]+z[4]*tmp1[25]+z[74]*tmp1[26]+z[47]*tmp1[27]+z[20]*tmp1[28]+z[90]*tmp1[29]+z[63]*tmp1[30]
                    +z[36]*tmp1[31]+z[9]*tmp1[32]+z[79]*tmp1[33]+z[52]*tmp1[34]+z[25]*tmp1[35]+z[95]*tmp1[36]+z[68]*tmp1[37]+z[41]*tmp1[38]+z[14]*tmp1[39]+z[84]*tmp1[40]
                    +z[57]*tmp1[41]+z[30]*tmp1[42]+z[3]*tmp1[43]+z[73]*tmp1[44]+z[46]*tmp1[45]+z[19]*tmp1[46]+z[89]*tmp1[47]+z[62]*tmp1[48]+z[35]*tmp1[49]+z[8]*tmp1[50]
                    +z[78]*tmp1[51]+z[51]*tmp1[52]+z[24]*tmp1[53]+z[94]*tmp1[54]+z[67]*tmp1[55]+z[40]*tmp1[56]+z[13]*tmp1[57]+z[83]*tmp1[58]+z[56]*tmp1[59]+z[29]*tmp1[60]
                    +z[2]*tmp1[61]+z[72]*tmp1[62]+z[45]*tmp1[63]+z[18]*tmp1[64]+z[88]*tmp1[65]+z[61]*tmp1[66]+z[34]*tmp1[67]+z[7]*tmp1[68]+z[77]*tmp1[69]+z[50]*tmp1[70]
                    +z[23]*tmp1[71]+z[93]*tmp1[72]+z[66]*tmp1[73]+z[39]*tmp1[74]+z[12]*tmp1[75]+z[82]*tmp1[76]+z[55]*tmp1[77]+z[28]*tmp1[78]+z[1]*tmp1[79]+z[71]*tmp1[80]
                    +z[44]*tmp1[81]+z[17]*tmp1[82]+z[87]*tmp1[83]+z[60]*tmp1[84]+z[33]*tmp1[85]+z[6]*tmp1[86]+z[76]*tmp1[87]+z[49]*tmp1[88]+z[22]*tmp1[89]+z[92]*tmp1[90]
                    +z[65]*tmp1[91]+z[38]*tmp1[92]+z[11]*tmp1[93]+z[81]*tmp1[94]+z[54]*tmp1[95]+z[27]*tmp1[96];
                    tab[nb_tmp3+71*nb3]=tab[nb_tmp3+71*nb3]+z[0]*tmp1[0]
                    +z[71]*tmp1[1]+z[45]*tmp1[2]+z[19]*tmp1[3]+z[90]*tmp1[4]+z[64]*tmp1[5]+z[38]*tmp1[6]+z[12]*tmp1[7]+z[83]*tmp1[8]+z[57]*tmp1[9]+z[31]*tmp1[10]
                    +z[5]*tmp1[11]+z[76]*tmp1[12]+z[50]*tmp1[13]+z[24]*tmp1[14]+z[95]*tmp1[15]+z[69]*tmp1[16]+z[43]*tmp1[17]+z[17]*tmp1[18]+z[88]*tmp1[19]+z[62]*tmp1[20]
                    +z[36]*tmp1[21]+z[10]*tmp1[22]+z[81]*tmp1[23]+z[55]*tmp1[24]+z[29]*tmp1[25]+z[3]*tmp1[26]+z[74]*tmp1[27]+z[48]*tmp1[28]+z[22]*tmp1[29]+z[93]*tmp1[30]
                    +z[67]*tmp1[31]+z[41]*tmp1[32]+z[15]*tmp1[33]+z[86]*tmp1[34]+z[60]*tmp1[35]+z[34]*tmp1[36]+z[8]*tmp1[37]+z[79]*tmp1[38]+z[53]*tmp1[39]+z[27]*tmp1[40]
                    +z[1]*tmp1[41]+z[72]*tmp1[42]+z[46]*tmp1[43]+z[20]*tmp1[44]+z[91]*tmp1[45]+z[65]*tmp1[46]+z[39]*tmp1[47]+z[13]*tmp1[48]+z[84]*tmp1[49]+z[58]*tmp1[50]
                    +z[32]*tmp1[51]+z[6]*tmp1[52]+z[77]*tmp1[53]+z[51]*tmp1[54]+z[25]*tmp1[55]+z[96]*tmp1[56]+z[70]*tmp1[57]+z[44]*tmp1[58]+z[18]*tmp1[59]+z[89]*tmp1[60]
                    +z[63]*tmp1[61]+z[37]*tmp1[62]+z[11]*tmp1[63]+z[82]*tmp1[64]+z[56]*tmp1[65]+z[30]*tmp1[66]+z[4]*tmp1[67]+z[75]*tmp1[68]+z[49]*tmp1[69]+z[23]*tmp1[70]
                    +z[94]*tmp1[71]+z[68]*tmp1[72]+z[42]*tmp1[73]+z[16]*tmp1[74]+z[87]*tmp1[75]+z[61]*tmp1[76]+z[35]*tmp1[77]+z[9]*tmp1[78]+z[80]*tmp1[79]+z[54]*tmp1[80]
                    +z[28]*tmp1[81]+z[2]*tmp1[82]+z[73]*tmp1[83]+z[47]*tmp1[84]+z[21]*tmp1[85]+z[92]*tmp1[86]+z[66]*tmp1[87]+z[40]*tmp1[88]+z[14]*tmp1[89]+z[85]*tmp1[90]
                    +z[59]*tmp1[91]+z[33]*tmp1[92]+z[7]*tmp1[93]+z[78]*tmp1[94]+z[52]*tmp1[95]+z[26]*tmp1[96];
                    tab[nb_tmp3+72*nb3]=tab[nb_tmp3+72*nb3]+z[0]*tmp1[0]
                    +z[72]*tmp1[1]+z[47]*tmp1[2]+z[22]*tmp1[3]+z[94]*tmp1[4]+z[69]*tmp1[5]+z[44]*tmp1[6]+z[19]*tmp1[7]+z[91]*tmp1[8]+z[66]*tmp1[9]+z[41]*tmp1[10]
                    +z[16]*tmp1[11]+z[88]*tmp1[12]+z[63]*tmp1[13]+z[38]*tmp1[14]+z[13]*tmp1[15]+z[85]*tmp1[16]+z[60]*tmp1[17]+z[35]*tmp1[18]+z[10]*tmp1[19]+z[82]*tmp1[20]
                    +z[57]*tmp1[21]+z[32]*tmp1[22]+z[7]*tmp1[23]+z[79]*tmp1[24]+z[54]*tmp1[25]+z[29]*tmp1[26]+z[4]*tmp1[27]+z[76]*tmp1[28]+z[51]*tmp1[29]+z[26]*tmp1[30]
                    +z[1]*tmp1[31]+z[73]*tmp1[32]+z[48]*tmp1[33]+z[23]*tmp1[34]+z[95]*tmp1[35]+z[70]*tmp1[36]+z[45]*tmp1[37]+z[20]*tmp1[38]+z[92]*tmp1[39]+z[67]*tmp1[40]
                    +z[42]*tmp1[41]+z[17]*tmp1[42]+z[89]*tmp1[43]+z[64]*tmp1[44]+z[39]*tmp1[45]+z[14]*tmp1[46]+z[86]*tmp1[47]+z[61]*tmp1[48]+z[36]*tmp1[49]+z[11]*tmp1[50]
                    +z[83]*tmp1[51]+z[58]*tmp1[52]+z[33]*tmp1[53]+z[8]*tmp1[54]+z[80]*tmp1[55]+z[55]*tmp1[56]+z[30]*tmp1[57]+z[5]*tmp1[58]+z[77]*tmp1[59]+z[52]*tmp1[60]
                    +z[27]*tmp1[61]+z[2]*tmp1[62]+z[74]*tmp1[63]+z[49]*tmp1[64]+z[24]*tmp1[65]+z[96]*tmp1[66]+z[71]*tmp1[67]+z[46]*tmp1[68]+z[21]*tmp1[69]+z[93]*tmp1[70]
                    +z[68]*tmp1[71]+z[43]*tmp1[72]+z[18]*tmp1[73]+z[90]*tmp1[74]+z[65]*tmp1[75]+z[40]*tmp1[76]+z[15]*tmp1[77]+z[87]*tmp1[78]+z[62]*tmp1[79]+z[37]*tmp1[80]
                    +z[12]*tmp1[81]+z[84]*tmp1[82]+z[59]*tmp1[83]+z[34]*tmp1[84]+z[9]*tmp1[85]+z[81]*tmp1[86]+z[56]*tmp1[87]+z[31]*tmp1[88]+z[6]*tmp1[89]+z[78]*tmp1[90]
                    +z[53]*tmp1[91]+z[28]*tmp1[92]+z[3]*tmp1[93]+z[75]*tmp1[94]+z[50]*tmp1[95]+z[25]*tmp1[96];
                    tab[nb_tmp3+73*nb3]=tab[nb_tmp3+73*nb3]+z[0]*tmp1[0]
                    +z[73]*tmp1[1]+z[49]*tmp1[2]+z[25]*tmp1[3]+z[1]*tmp1[4]+z[74]*tmp1[5]+z[50]*tmp1[6]+z[26]*tmp1[7]+z[2]*tmp1[8]+z[75]*tmp1[9]+z[51]*tmp1[10]
                    +z[27]*tmp1[11]+z[3]*tmp1[12]+z[76]*tmp1[13]+z[52]*tmp1[14]+z[28]*tmp1[15]+z[4]*tmp1[16]+z[77]*tmp1[17]+z[53]*tmp1[18]+z[29]*tmp1[19]+z[5]*tmp1[20]
                    +z[78]*tmp1[21]+z[54]*tmp1[22]+z[30]*tmp1[23]+z[6]*tmp1[24]+z[79]*tmp1[25]+z[55]*tmp1[26]+z[31]*tmp1[27]+z[7]*tmp1[28]+z[80]*tmp1[29]+z[56]*tmp1[30]
                    +z[32]*tmp1[31]+z[8]*tmp1[32]+z[81]*tmp1[33]+z[57]*tmp1[34]+z[33]*tmp1[35]+z[9]*tmp1[36]+z[82]*tmp1[37]+z[58]*tmp1[38]+z[34]*tmp1[39]+z[10]*tmp1[40]
                    +z[83]*tmp1[41]+z[59]*tmp1[42]+z[35]*tmp1[43]+z[11]*tmp1[44]+z[84]*tmp1[45]+z[60]*tmp1[46]+z[36]*tmp1[47]+z[12]*tmp1[48]+z[85]*tmp1[49]+z[61]*tmp1[50]
                    +z[37]*tmp1[51]+z[13]*tmp1[52]+z[86]*tmp1[53]+z[62]*tmp1[54]+z[38]*tmp1[55]+z[14]*tmp1[56]+z[87]*tmp1[57]+z[63]*tmp1[58]+z[39]*tmp1[59]+z[15]*tmp1[60]
                    +z[88]*tmp1[61]+z[64]*tmp1[62]+z[40]*tmp1[63]+z[16]*tmp1[64]+z[89]*tmp1[65]+z[65]*tmp1[66]+z[41]*tmp1[67]+z[17]*tmp1[68]+z[90]*tmp1[69]+z[66]*tmp1[70]
                    +z[42]*tmp1[71]+z[18]*tmp1[72]+z[91]*tmp1[73]+z[67]*tmp1[74]+z[43]*tmp1[75]+z[19]*tmp1[76]+z[92]*tmp1[77]+z[68]*tmp1[78]+z[44]*tmp1[79]+z[20]*tmp1[80]
                    +z[93]*tmp1[81]+z[69]*tmp1[82]+z[45]*tmp1[83]+z[21]*tmp1[84]+z[94]*tmp1[85]+z[70]*tmp1[86]+z[46]*tmp1[87]+z[22]*tmp1[88]+z[95]*tmp1[89]+z[71]*tmp1[90]
                    +z[47]*tmp1[91]+z[23]*tmp1[92]+z[96]*tmp1[93]+z[72]*tmp1[94]+z[48]*tmp1[95]+z[24]*tmp1[96];
                    tab[nb_tmp3+74*nb3]=tab[nb_tmp3+74*nb3]+z[0]*tmp1[0]
                    +z[74]*tmp1[1]+z[51]*tmp1[2]+z[28]*tmp1[3]+z[5]*tmp1[4]+z[79]*tmp1[5]+z[56]*tmp1[6]+z[33]*tmp1[7]+z[10]*tmp1[8]+z[84]*tmp1[9]+z[61]*tmp1[10]
                    +z[38]*tmp1[11]+z[15]*tmp1[12]+z[89]*tmp1[13]+z[66]*tmp1[14]+z[43]*tmp1[15]+z[20]*tmp1[16]+z[94]*tmp1[17]+z[71]*tmp1[18]+z[48]*tmp1[19]+z[25]*tmp1[20]
                    +z[2]*tmp1[21]+z[76]*tmp1[22]+z[53]*tmp1[23]+z[30]*tmp1[24]+z[7]*tmp1[25]+z[81]*tmp1[26]+z[58]*tmp1[27]+z[35]*tmp1[28]+z[12]*tmp1[29]+z[86]*tmp1[30]
                    +z[63]*tmp1[31]+z[40]*tmp1[32]+z[17]*tmp1[33]+z[91]*tmp1[34]+z[68]*tmp1[35]+z[45]*tmp1[36]+z[22]*tmp1[37]+z[96]*tmp1[38]+z[73]*tmp1[39]+z[50]*tmp1[40]
                    +z[27]*tmp1[41]+z[4]*tmp1[42]+z[78]*tmp1[43]+z[55]*tmp1[44]+z[32]*tmp1[45]+z[9]*tmp1[46]+z[83]*tmp1[47]+z[60]*tmp1[48]+z[37]*tmp1[49]+z[14]*tmp1[50]
                    +z[88]*tmp1[51]+z[65]*tmp1[52]+z[42]*tmp1[53]+z[19]*tmp1[54]+z[93]*tmp1[55]+z[70]*tmp1[56]+z[47]*tmp1[57]+z[24]*tmp1[58]+z[1]*tmp1[59]+z[75]*tmp1[60]
                    +z[52]*tmp1[61]+z[29]*tmp1[62]+z[6]*tmp1[63]+z[80]*tmp1[64]+z[57]*tmp1[65]+z[34]*tmp1[66]+z[11]*tmp1[67]+z[85]*tmp1[68]+z[62]*tmp1[69]+z[39]*tmp1[70]
                    +z[16]*tmp1[71]+z[90]*tmp1[72]+z[67]*tmp1[73]+z[44]*tmp1[74]+z[21]*tmp1[75]+z[95]*tmp1[76]+z[72]*tmp1[77]+z[49]*tmp1[78]+z[26]*tmp1[79]+z[3]*tmp1[80]
                    +z[77]*tmp1[81]+z[54]*tmp1[82]+z[31]*tmp1[83]+z[8]*tmp1[84]+z[82]*tmp1[85]+z[59]*tmp1[86]+z[36]*tmp1[87]+z[13]*tmp1[88]+z[87]*tmp1[89]+z[64]*tmp1[90]
                    +z[41]*tmp1[91]+z[18]*tmp1[92]+z[92]*tmp1[93]+z[69]*tmp1[94]+z[46]*tmp1[95]+z[23]*tmp1[96];
                    tab[nb_tmp3+75*nb3]=tab[nb_tmp3+75*nb3]+z[0]*tmp1[0]
                    +z[75]*tmp1[1]+z[53]*tmp1[2]+z[31]*tmp1[3]+z[9]*tmp1[4]+z[84]*tmp1[5]+z[62]*tmp1[6]+z[40]*tmp1[7]+z[18]*tmp1[8]+z[93]*tmp1[9]+z[71]*tmp1[10]
                    +z[49]*tmp1[11]+z[27]*tmp1[12]+z[5]*tmp1[13]+z[80]*tmp1[14]+z[58]*tmp1[15]+z[36]*tmp1[16]+z[14]*tmp1[17]+z[89]*tmp1[18]+z[67]*tmp1[19]+z[45]*tmp1[20]
                    +z[23]*tmp1[21]+z[1]*tmp1[22]+z[76]*tmp1[23]+z[54]*tmp1[24]+z[32]*tmp1[25]+z[10]*tmp1[26]+z[85]*tmp1[27]+z[63]*tmp1[28]+z[41]*tmp1[29]+z[19]*tmp1[30]
                    +z[94]*tmp1[31]+z[72]*tmp1[32]+z[50]*tmp1[33]+z[28]*tmp1[34]+z[6]*tmp1[35]+z[81]*tmp1[36]+z[59]*tmp1[37]+z[37]*tmp1[38]+z[15]*tmp1[39]+z[90]*tmp1[40]
                    +z[68]*tmp1[41]+z[46]*tmp1[42]+z[24]*tmp1[43]+z[2]*tmp1[44]+z[77]*tmp1[45]+z[55]*tmp1[46]+z[33]*tmp1[47]+z[11]*tmp1[48]+z[86]*tmp1[49]+z[64]*tmp1[50]
                    +z[42]*tmp1[51]+z[20]*tmp1[52]+z[95]*tmp1[53]+z[73]*tmp1[54]+z[51]*tmp1[55]+z[29]*tmp1[56]+z[7]*tmp1[57]+z[82]*tmp1[58]+z[60]*tmp1[59]+z[38]*tmp1[60]
                    +z[16]*tmp1[61]+z[91]*tmp1[62]+z[69]*tmp1[63]+z[47]*tmp1[64]+z[25]*tmp1[65]+z[3]*tmp1[66]+z[78]*tmp1[67]+z[56]*tmp1[68]+z[34]*tmp1[69]+z[12]*tmp1[70]
                    +z[87]*tmp1[71]+z[65]*tmp1[72]+z[43]*tmp1[73]+z[21]*tmp1[74]+z[96]*tmp1[75]+z[74]*tmp1[76]+z[52]*tmp1[77]+z[30]*tmp1[78]+z[8]*tmp1[79]+z[83]*tmp1[80]
                    +z[61]*tmp1[81]+z[39]*tmp1[82]+z[17]*tmp1[83]+z[92]*tmp1[84]+z[70]*tmp1[85]+z[48]*tmp1[86]+z[26]*tmp1[87]+z[4]*tmp1[88]+z[79]*tmp1[89]+z[57]*tmp1[90]
                    +z[35]*tmp1[91]+z[13]*tmp1[92]+z[88]*tmp1[93]+z[66]*tmp1[94]+z[44]*tmp1[95]+z[22]*tmp1[96];
                    tab[nb_tmp3+76*nb3]=tab[nb_tmp3+76*nb3]+z[0]*tmp1[0]
                    +z[76]*tmp1[1]+z[55]*tmp1[2]+z[34]*tmp1[3]+z[13]*tmp1[4]+z[89]*tmp1[5]+z[68]*tmp1[6]+z[47]*tmp1[7]+z[26]*tmp1[8]+z[5]*tmp1[9]+z[81]*tmp1[10]
                    +z[60]*tmp1[11]+z[39]*tmp1[12]+z[18]*tmp1[13]+z[94]*tmp1[14]+z[73]*tmp1[15]+z[52]*tmp1[16]+z[31]*tmp1[17]+z[10]*tmp1[18]+z[86]*tmp1[19]+z[65]*tmp1[20]
                    +z[44]*tmp1[21]+z[23]*tmp1[22]+z[2]*tmp1[23]+z[78]*tmp1[24]+z[57]*tmp1[25]+z[36]*tmp1[26]+z[15]*tmp1[27]+z[91]*tmp1[28]+z[70]*tmp1[29]+z[49]*tmp1[30]
                    +z[28]*tmp1[31]+z[7]*tmp1[32]+z[83]*tmp1[33]+z[62]*tmp1[34]+z[41]*tmp1[35]+z[20]*tmp1[36]+z[96]*tmp1[37]+z[75]*tmp1[38]+z[54]*tmp1[39]+z[33]*tmp1[40]
                    +z[12]*tmp1[41]+z[88]*tmp1[42]+z[67]*tmp1[43]+z[46]*tmp1[44]+z[25]*tmp1[45]+z[4]*tmp1[46]+z[80]*tmp1[47]+z[59]*tmp1[48]+z[38]*tmp1[49]+z[17]*tmp1[50]
                    +z[93]*tmp1[51]+z[72]*tmp1[52]+z[51]*tmp1[53]+z[30]*tmp1[54]+z[9]*tmp1[55]+z[85]*tmp1[56]+z[64]*tmp1[57]+z[43]*tmp1[58]+z[22]*tmp1[59]+z[1]*tmp1[60]
                    +z[77]*tmp1[61]+z[56]*tmp1[62]+z[35]*tmp1[63]+z[14]*tmp1[64]+z[90]*tmp1[65]+z[69]*tmp1[66]+z[48]*tmp1[67]+z[27]*tmp1[68]+z[6]*tmp1[69]+z[82]*tmp1[70]
                    +z[61]*tmp1[71]+z[40]*tmp1[72]+z[19]*tmp1[73]+z[95]*tmp1[74]+z[74]*tmp1[75]+z[53]*tmp1[76]+z[32]*tmp1[77]+z[11]*tmp1[78]+z[87]*tmp1[79]+z[66]*tmp1[80]
                    +z[45]*tmp1[81]+z[24]*tmp1[82]+z[3]*tmp1[83]+z[79]*tmp1[84]+z[58]*tmp1[85]+z[37]*tmp1[86]+z[16]*tmp1[87]+z[92]*tmp1[88]+z[71]*tmp1[89]+z[50]*tmp1[90]
                    +z[29]*tmp1[91]+z[8]*tmp1[92]+z[84]*tmp1[93]+z[63]*tmp1[94]+z[42]*tmp1[95]+z[21]*tmp1[96];
                    tab[nb_tmp3+77*nb3]=tab[nb_tmp3+77*nb3]+z[0]*tmp1[0]
                    +z[77]*tmp1[1]+z[57]*tmp1[2]+z[37]*tmp1[3]+z[17]*tmp1[4]+z[94]*tmp1[5]+z[74]*tmp1[6]+z[54]*tmp1[7]+z[34]*tmp1[8]+z[14]*tmp1[9]+z[91]*tmp1[10]
                    +z[71]*tmp1[11]+z[51]*tmp1[12]+z[31]*tmp1[13]+z[11]*tmp1[14]+z[88]*tmp1[15]+z[68]*tmp1[16]+z[48]*tmp1[17]+z[28]*tmp1[18]+z[8]*tmp1[19]+z[85]*tmp1[20]
                    +z[65]*tmp1[21]+z[45]*tmp1[22]+z[25]*tmp1[23]+z[5]*tmp1[24]+z[82]*tmp1[25]+z[62]*tmp1[26]+z[42]*tmp1[27]+z[22]*tmp1[28]+z[2]*tmp1[29]+z[79]*tmp1[30]
                    +z[59]*tmp1[31]+z[39]*tmp1[32]+z[19]*tmp1[33]+z[96]*tmp1[34]+z[76]*tmp1[35]+z[56]*tmp1[36]+z[36]*tmp1[37]+z[16]*tmp1[38]+z[93]*tmp1[39]+z[73]*tmp1[40]
                    +z[53]*tmp1[41]+z[33]*tmp1[42]+z[13]*tmp1[43]+z[90]*tmp1[44]+z[70]*tmp1[45]+z[50]*tmp1[46]+z[30]*tmp1[47]+z[10]*tmp1[48]+z[87]*tmp1[49]+z[67]*tmp1[50]
                    +z[47]*tmp1[51]+z[27]*tmp1[52]+z[7]*tmp1[53]+z[84]*tmp1[54]+z[64]*tmp1[55]+z[44]*tmp1[56]+z[24]*tmp1[57]+z[4]*tmp1[58]+z[81]*tmp1[59]+z[61]*tmp1[60]
                    +z[41]*tmp1[61]+z[21]*tmp1[62]+z[1]*tmp1[63]+z[78]*tmp1[64]+z[58]*tmp1[65]+z[38]*tmp1[66]+z[18]*tmp1[67]+z[95]*tmp1[68]+z[75]*tmp1[69]+z[55]*tmp1[70]
                    +z[35]*tmp1[71]+z[15]*tmp1[72]+z[92]*tmp1[73]+z[72]*tmp1[74]+z[52]*tmp1[75]+z[32]*tmp1[76]+z[12]*tmp1[77]+z[89]*tmp1[78]+z[69]*tmp1[79]+z[49]*tmp1[80]
                    +z[29]*tmp1[81]+z[9]*tmp1[82]+z[86]*tmp1[83]+z[66]*tmp1[84]+z[46]*tmp1[85]+z[26]*tmp1[86]+z[6]*tmp1[87]+z[83]*tmp1[88]+z[63]*tmp1[89]+z[43]*tmp1[90]
                    +z[23]*tmp1[91]+z[3]*tmp1[92]+z[80]*tmp1[93]+z[60]*tmp1[94]+z[40]*tmp1[95]+z[20]*tmp1[96];
                    tab[nb_tmp3+78*nb3]=tab[nb_tmp3+78*nb3]+z[0]*tmp1[0]
                    +z[78]*tmp1[1]+z[59]*tmp1[2]+z[40]*tmp1[3]+z[21]*tmp1[4]+z[2]*tmp1[5]+z[80]*tmp1[6]+z[61]*tmp1[7]+z[42]*tmp1[8]+z[23]*tmp1[9]+z[4]*tmp1[10]
                    +z[82]*tmp1[11]+z[63]*tmp1[12]+z[44]*tmp1[13]+z[25]*tmp1[14]+z[6]*tmp1[15]+z[84]*tmp1[16]+z[65]*tmp1[17]+z[46]*tmp1[18]+z[27]*tmp1[19]+z[8]*tmp1[20]
                    +z[86]*tmp1[21]+z[67]*tmp1[22]+z[48]*tmp1[23]+z[29]*tmp1[24]+z[10]*tmp1[25]+z[88]*tmp1[26]+z[69]*tmp1[27]+z[50]*tmp1[28]+z[31]*tmp1[29]+z[12]*tmp1[30]
                    +z[90]*tmp1[31]+z[71]*tmp1[32]+z[52]*tmp1[33]+z[33]*tmp1[34]+z[14]*tmp1[35]+z[92]*tmp1[36]+z[73]*tmp1[37]+z[54]*tmp1[38]+z[35]*tmp1[39]+z[16]*tmp1[40]
                    +z[94]*tmp1[41]+z[75]*tmp1[42]+z[56]*tmp1[43]+z[37]*tmp1[44]+z[18]*tmp1[45]+z[96]*tmp1[46]+z[77]*tmp1[47]+z[58]*tmp1[48]+z[39]*tmp1[49]+z[20]*tmp1[50]
                    +z[1]*tmp1[51]+z[79]*tmp1[52]+z[60]*tmp1[53]+z[41]*tmp1[54]+z[22]*tmp1[55]+z[3]*tmp1[56]+z[81]*tmp1[57]+z[62]*tmp1[58]+z[43]*tmp1[59]+z[24]*tmp1[60]
                    +z[5]*tmp1[61]+z[83]*tmp1[62]+z[64]*tmp1[63]+z[45]*tmp1[64]+z[26]*tmp1[65]+z[7]*tmp1[66]+z[85]*tmp1[67]+z[66]*tmp1[68]+z[47]*tmp1[69]+z[28]*tmp1[70]
                    +z[9]*tmp1[71]+z[87]*tmp1[72]+z[68]*tmp1[73]+z[49]*tmp1[74]+z[30]*tmp1[75]+z[11]*tmp1[76]+z[89]*tmp1[77]+z[70]*tmp1[78]+z[51]*tmp1[79]+z[32]*tmp1[80]
                    +z[13]*tmp1[81]+z[91]*tmp1[82]+z[72]*tmp1[83]+z[53]*tmp1[84]+z[34]*tmp1[85]+z[15]*tmp1[86]+z[93]*tmp1[87]+z[74]*tmp1[88]+z[55]*tmp1[89]+z[36]*tmp1[90]
                    +z[17]*tmp1[91]+z[95]*tmp1[92]+z[76]*tmp1[93]+z[57]*tmp1[94]+z[38]*tmp1[95]+z[19]*tmp1[96];
                    tab[nb_tmp3+79*nb3]=tab[nb_tmp3+79*nb3]+z[0]*tmp1[0]
                    +z[79]*tmp1[1]+z[61]*tmp1[2]+z[43]*tmp1[3]+z[25]*tmp1[4]+z[7]*tmp1[5]+z[86]*tmp1[6]+z[68]*tmp1[7]+z[50]*tmp1[8]+z[32]*tmp1[9]+z[14]*tmp1[10]
                    +z[93]*tmp1[11]+z[75]*tmp1[12]+z[57]*tmp1[13]+z[39]*tmp1[14]+z[21]*tmp1[15]+z[3]*tmp1[16]+z[82]*tmp1[17]+z[64]*tmp1[18]+z[46]*tmp1[19]+z[28]*tmp1[20]
                    +z[10]*tmp1[21]+z[89]*tmp1[22]+z[71]*tmp1[23]+z[53]*tmp1[24]+z[35]*tmp1[25]+z[17]*tmp1[26]+z[96]*tmp1[27]+z[78]*tmp1[28]+z[60]*tmp1[29]+z[42]*tmp1[30]
                    +z[24]*tmp1[31]+z[6]*tmp1[32]+z[85]*tmp1[33]+z[67]*tmp1[34]+z[49]*tmp1[35]+z[31]*tmp1[36]+z[13]*tmp1[37]+z[92]*tmp1[38]+z[74]*tmp1[39]+z[56]*tmp1[40]
                    +z[38]*tmp1[41]+z[20]*tmp1[42]+z[2]*tmp1[43]+z[81]*tmp1[44]+z[63]*tmp1[45]+z[45]*tmp1[46]+z[27]*tmp1[47]+z[9]*tmp1[48]+z[88]*tmp1[49]+z[70]*tmp1[50]
                    +z[52]*tmp1[51]+z[34]*tmp1[52]+z[16]*tmp1[53]+z[95]*tmp1[54]+z[77]*tmp1[55]+z[59]*tmp1[56]+z[41]*tmp1[57]+z[23]*tmp1[58]+z[5]*tmp1[59]+z[84]*tmp1[60]
                    +z[66]*tmp1[61]+z[48]*tmp1[62]+z[30]*tmp1[63]+z[12]*tmp1[64]+z[91]*tmp1[65]+z[73]*tmp1[66]+z[55]*tmp1[67]+z[37]*tmp1[68]+z[19]*tmp1[69]+z[1]*tmp1[70]
                    +z[80]*tmp1[71]+z[62]*tmp1[72]+z[44]*tmp1[73]+z[26]*tmp1[74]+z[8]*tmp1[75]+z[87]*tmp1[76]+z[69]*tmp1[77]+z[51]*tmp1[78]+z[33]*tmp1[79]+z[15]*tmp1[80]
                    +z[94]*tmp1[81]+z[76]*tmp1[82]+z[58]*tmp1[83]+z[40]*tmp1[84]+z[22]*tmp1[85]+z[4]*tmp1[86]+z[83]*tmp1[87]+z[65]*tmp1[88]+z[47]*tmp1[89]+z[29]*tmp1[90]
                    +z[11]*tmp1[91]+z[90]*tmp1[92]+z[72]*tmp1[93]+z[54]*tmp1[94]+z[36]*tmp1[95]+z[18]*tmp1[96];
                    tab[nb_tmp3+80*nb3]=tab[nb_tmp3+80*nb3]+z[0]*tmp1[0]
                    +z[80]*tmp1[1]+z[63]*tmp1[2]+z[46]*tmp1[3]+z[29]*tmp1[4]+z[12]*tmp1[5]+z[92]*tmp1[6]+z[75]*tmp1[7]+z[58]*tmp1[8]+z[41]*tmp1[9]+z[24]*tmp1[10]
                    +z[7]*tmp1[11]+z[87]*tmp1[12]+z[70]*tmp1[13]+z[53]*tmp1[14]+z[36]*tmp1[15]+z[19]*tmp1[16]+z[2]*tmp1[17]+z[82]*tmp1[18]+z[65]*tmp1[19]+z[48]*tmp1[20]
                    +z[31]*tmp1[21]+z[14]*tmp1[22]+z[94]*tmp1[23]+z[77]*tmp1[24]+z[60]*tmp1[25]+z[43]*tmp1[26]+z[26]*tmp1[27]+z[9]*tmp1[28]+z[89]*tmp1[29]+z[72]*tmp1[30]
                    +z[55]*tmp1[31]+z[38]*tmp1[32]+z[21]*tmp1[33]+z[4]*tmp1[34]+z[84]*tmp1[35]+z[67]*tmp1[36]+z[50]*tmp1[37]+z[33]*tmp1[38]+z[16]*tmp1[39]+z[96]*tmp1[40]
                    +z[79]*tmp1[41]+z[62]*tmp1[42]+z[45]*tmp1[43]+z[28]*tmp1[44]+z[11]*tmp1[45]+z[91]*tmp1[46]+z[74]*tmp1[47]+z[57]*tmp1[48]+z[40]*tmp1[49]+z[23]*tmp1[50]
                    +z[6]*tmp1[51]+z[86]*tmp1[52]+z[69]*tmp1[53]+z[52]*tmp1[54]+z[35]*tmp1[55]+z[18]*tmp1[56]+z[1]*tmp1[57]+z[81]*tmp1[58]+z[64]*tmp1[59]+z[47]*tmp1[60]
                    +z[30]*tmp1[61]+z[13]*tmp1[62]+z[93]*tmp1[63]+z[76]*tmp1[64]+z[59]*tmp1[65]+z[42]*tmp1[66]+z[25]*tmp1[67]+z[8]*tmp1[68]+z[88]*tmp1[69]+z[71]*tmp1[70]
                    +z[54]*tmp1[71]+z[37]*tmp1[72]+z[20]*tmp1[73]+z[3]*tmp1[74]+z[83]*tmp1[75]+z[66]*tmp1[76]+z[49]*tmp1[77]+z[32]*tmp1[78]+z[15]*tmp1[79]+z[95]*tmp1[80]
                    +z[78]*tmp1[81]+z[61]*tmp1[82]+z[44]*tmp1[83]+z[27]*tmp1[84]+z[10]*tmp1[85]+z[90]*tmp1[86]+z[73]*tmp1[87]+z[56]*tmp1[88]+z[39]*tmp1[89]+z[22]*tmp1[90]
                    +z[5]*tmp1[91]+z[85]*tmp1[92]+z[68]*tmp1[93]+z[51]*tmp1[94]+z[34]*tmp1[95]+z[17]*tmp1[96];
                    tab[nb_tmp3+81*nb3]=tab[nb_tmp3+81*nb3]+z[0]*tmp1[0]
                    +z[81]*tmp1[1]+z[65]*tmp1[2]+z[49]*tmp1[3]+z[33]*tmp1[4]+z[17]*tmp1[5]+z[1]*tmp1[6]+z[82]*tmp1[7]+z[66]*tmp1[8]+z[50]*tmp1[9]+z[34]*tmp1[10]
                    +z[18]*tmp1[11]+z[2]*tmp1[12]+z[83]*tmp1[13]+z[67]*tmp1[14]+z[51]*tmp1[15]+z[35]*tmp1[16]+z[19]*tmp1[17]+z[3]*tmp1[18]+z[84]*tmp1[19]+z[68]*tmp1[20]
                    +z[52]*tmp1[21]+z[36]*tmp1[22]+z[20]*tmp1[23]+z[4]*tmp1[24]+z[85]*tmp1[25]+z[69]*tmp1[26]+z[53]*tmp1[27]+z[37]*tmp1[28]+z[21]*tmp1[29]+z[5]*tmp1[30]
                    +z[86]*tmp1[31]+z[70]*tmp1[32]+z[54]*tmp1[33]+z[38]*tmp1[34]+z[22]*tmp1[35]+z[6]*tmp1[36]+z[87]*tmp1[37]+z[71]*tmp1[38]+z[55]*tmp1[39]+z[39]*tmp1[40]
                    +z[23]*tmp1[41]+z[7]*tmp1[42]+z[88]*tmp1[43]+z[72]*tmp1[44]+z[56]*tmp1[45]+z[40]*tmp1[46]+z[24]*tmp1[47]+z[8]*tmp1[48]+z[89]*tmp1[49]+z[73]*tmp1[50]
                    +z[57]*tmp1[51]+z[41]*tmp1[52]+z[25]*tmp1[53]+z[9]*tmp1[54]+z[90]*tmp1[55]+z[74]*tmp1[56]+z[58]*tmp1[57]+z[42]*tmp1[58]+z[26]*tmp1[59]+z[10]*tmp1[60]
                    +z[91]*tmp1[61]+z[75]*tmp1[62]+z[59]*tmp1[63]+z[43]*tmp1[64]+z[27]*tmp1[65]+z[11]*tmp1[66]+z[92]*tmp1[67]+z[76]*tmp1[68]+z[60]*tmp1[69]+z[44]*tmp1[70]
                    +z[28]*tmp1[71]+z[12]*tmp1[72]+z[93]*tmp1[73]+z[77]*tmp1[74]+z[61]*tmp1[75]+z[45]*tmp1[76]+z[29]*tmp1[77]+z[13]*tmp1[78]+z[94]*tmp1[79]+z[78]*tmp1[80]
                    +z[62]*tmp1[81]+z[46]*tmp1[82]+z[30]*tmp1[83]+z[14]*tmp1[84]+z[95]*tmp1[85]+z[79]*tmp1[86]+z[63]*tmp1[87]+z[47]*tmp1[88]+z[31]*tmp1[89]+z[15]*tmp1[90]
                    +z[96]*tmp1[91]+z[80]*tmp1[92]+z[64]*tmp1[93]+z[48]*tmp1[94]+z[32]*tmp1[95]+z[16]*tmp1[96];
                    tab[nb_tmp3+82*nb3]=tab[nb_tmp3+82*nb3]+z[0]*tmp1[0]
                    +z[82]*tmp1[1]+z[67]*tmp1[2]+z[52]*tmp1[3]+z[37]*tmp1[4]+z[22]*tmp1[5]+z[7]*tmp1[6]+z[89]*tmp1[7]+z[74]*tmp1[8]+z[59]*tmp1[9]+z[44]*tmp1[10]
                    +z[29]*tmp1[11]+z[14]*tmp1[12]+z[96]*tmp1[13]+z[81]*tmp1[14]+z[66]*tmp1[15]+z[51]*tmp1[16]+z[36]*tmp1[17]+z[21]*tmp1[18]+z[6]*tmp1[19]+z[88]*tmp1[20]
                    +z[73]*tmp1[21]+z[58]*tmp1[22]+z[43]*tmp1[23]+z[28]*tmp1[24]+z[13]*tmp1[25]+z[95]*tmp1[26]+z[80]*tmp1[27]+z[65]*tmp1[28]+z[50]*tmp1[29]+z[35]*tmp1[30]
                    +z[20]*tmp1[31]+z[5]*tmp1[32]+z[87]*tmp1[33]+z[72]*tmp1[34]+z[57]*tmp1[35]+z[42]*tmp1[36]+z[27]*tmp1[37]+z[12]*tmp1[38]+z[94]*tmp1[39]+z[79]*tmp1[40]
                    +z[64]*tmp1[41]+z[49]*tmp1[42]+z[34]*tmp1[43]+z[19]*tmp1[44]+z[4]*tmp1[45]+z[86]*tmp1[46]+z[71]*tmp1[47]+z[56]*tmp1[48]+z[41]*tmp1[49]+z[26]*tmp1[50]
                    +z[11]*tmp1[51]+z[93]*tmp1[52]+z[78]*tmp1[53]+z[63]*tmp1[54]+z[48]*tmp1[55]+z[33]*tmp1[56]+z[18]*tmp1[57]+z[3]*tmp1[58]+z[85]*tmp1[59]+z[70]*tmp1[60]
                    +z[55]*tmp1[61]+z[40]*tmp1[62]+z[25]*tmp1[63]+z[10]*tmp1[64]+z[92]*tmp1[65]+z[77]*tmp1[66]+z[62]*tmp1[67]+z[47]*tmp1[68]+z[32]*tmp1[69]+z[17]*tmp1[70]
                    +z[2]*tmp1[71]+z[84]*tmp1[72]+z[69]*tmp1[73]+z[54]*tmp1[74]+z[39]*tmp1[75]+z[24]*tmp1[76]+z[9]*tmp1[77]+z[91]*tmp1[78]+z[76]*tmp1[79]+z[61]*tmp1[80]
                    +z[46]*tmp1[81]+z[31]*tmp1[82]+z[16]*tmp1[83]+z[1]*tmp1[84]+z[83]*tmp1[85]+z[68]*tmp1[86]+z[53]*tmp1[87]+z[38]*tmp1[88]+z[23]*tmp1[89]+z[8]*tmp1[90]
                    +z[90]*tmp1[91]+z[75]*tmp1[92]+z[60]*tmp1[93]+z[45]*tmp1[94]+z[30]*tmp1[95]+z[15]*tmp1[96];
                    tab[nb_tmp3+83*nb3]=tab[nb_tmp3+83*nb3]+z[0]*tmp1[0]
                    +z[83]*tmp1[1]+z[69]*tmp1[2]+z[55]*tmp1[3]+z[41]*tmp1[4]+z[27]*tmp1[5]+z[13]*tmp1[6]+z[96]*tmp1[7]+z[82]*tmp1[8]+z[68]*tmp1[9]+z[54]*tmp1[10]
                    +z[40]*tmp1[11]+z[26]*tmp1[12]+z[12]*tmp1[13]+z[95]*tmp1[14]+z[81]*tmp1[15]+z[67]*tmp1[16]+z[53]*tmp1[17]+z[39]*tmp1[18]+z[25]*tmp1[19]+z[11]*tmp1[20]
                    +z[94]*tmp1[21]+z[80]*tmp1[22]+z[66]*tmp1[23]+z[52]*tmp1[24]+z[38]*tmp1[25]+z[24]*tmp1[26]+z[10]*tmp1[27]+z[93]*tmp1[28]+z[79]*tmp1[29]+z[65]*tmp1[30]
                    +z[51]*tmp1[31]+z[37]*tmp1[32]+z[23]*tmp1[33]+z[9]*tmp1[34]+z[92]*tmp1[35]+z[78]*tmp1[36]+z[64]*tmp1[37]+z[50]*tmp1[38]+z[36]*tmp1[39]+z[22]*tmp1[40]
                    +z[8]*tmp1[41]+z[91]*tmp1[42]+z[77]*tmp1[43]+z[63]*tmp1[44]+z[49]*tmp1[45]+z[35]*tmp1[46]+z[21]*tmp1[47]+z[7]*tmp1[48]+z[90]*tmp1[49]+z[76]*tmp1[50]
                    +z[62]*tmp1[51]+z[48]*tmp1[52]+z[34]*tmp1[53]+z[20]*tmp1[54]+z[6]*tmp1[55]+z[89]*tmp1[56]+z[75]*tmp1[57]+z[61]*tmp1[58]+z[47]*tmp1[59]+z[33]*tmp1[60]
                    +z[19]*tmp1[61]+z[5]*tmp1[62]+z[88]*tmp1[63]+z[74]*tmp1[64]+z[60]*tmp1[65]+z[46]*tmp1[66]+z[32]*tmp1[67]+z[18]*tmp1[68]+z[4]*tmp1[69]+z[87]*tmp1[70]
                    +z[73]*tmp1[71]+z[59]*tmp1[72]+z[45]*tmp1[73]+z[31]*tmp1[74]+z[17]*tmp1[75]+z[3]*tmp1[76]+z[86]*tmp1[77]+z[72]*tmp1[78]+z[58]*tmp1[79]+z[44]*tmp1[80]
                    +z[30]*tmp1[81]+z[16]*tmp1[82]+z[2]*tmp1[83]+z[85]*tmp1[84]+z[71]*tmp1[85]+z[57]*tmp1[86]+z[43]*tmp1[87]+z[29]*tmp1[88]+z[15]*tmp1[89]+z[1]*tmp1[90]
                    +z[84]*tmp1[91]+z[70]*tmp1[92]+z[56]*tmp1[93]+z[42]*tmp1[94]+z[28]*tmp1[95]+z[14]*tmp1[96];
                    tab[nb_tmp3+84*nb3]=tab[nb_tmp3+84*nb3]+z[0]*tmp1[0]
                    +z[84]*tmp1[1]+z[71]*tmp1[2]+z[58]*tmp1[3]+z[45]*tmp1[4]+z[32]*tmp1[5]+z[19]*tmp1[6]+z[6]*tmp1[7]+z[90]*tmp1[8]+z[77]*tmp1[9]+z[64]*tmp1[10]
                    +z[51]*tmp1[11]+z[38]*tmp1[12]+z[25]*tmp1[13]+z[12]*tmp1[14]+z[96]*tmp1[15]+z[83]*tmp1[16]+z[70]*tmp1[17]+z[57]*tmp1[18]+z[44]*tmp1[19]+z[31]*tmp1[20]
                    +z[18]*tmp1[21]+z[5]*tmp1[22]+z[89]*tmp1[23]+z[76]*tmp1[24]+z[63]*tmp1[25]+z[50]*tmp1[26]+z[37]*tmp1[27]+z[24]*tmp1[28]+z[11]*tmp1[29]+z[95]*tmp1[30]
                    +z[82]*tmp1[31]+z[69]*tmp1[32]+z[56]*tmp1[33]+z[43]*tmp1[34]+z[30]*tmp1[35]+z[17]*tmp1[36]+z[4]*tmp1[37]+z[88]*tmp1[38]+z[75]*tmp1[39]+z[62]*tmp1[40]
                    +z[49]*tmp1[41]+z[36]*tmp1[42]+z[23]*tmp1[43]+z[10]*tmp1[44]+z[94]*tmp1[45]+z[81]*tmp1[46]+z[68]*tmp1[47]+z[55]*tmp1[48]+z[42]*tmp1[49]+z[29]*tmp1[50]
                    +z[16]*tmp1[51]+z[3]*tmp1[52]+z[87]*tmp1[53]+z[74]*tmp1[54]+z[61]*tmp1[55]+z[48]*tmp1[56]+z[35]*tmp1[57]+z[22]*tmp1[58]+z[9]*tmp1[59]+z[93]*tmp1[60]
                    +z[80]*tmp1[61]+z[67]*tmp1[62]+z[54]*tmp1[63]+z[41]*tmp1[64]+z[28]*tmp1[65]+z[15]*tmp1[66]+z[2]*tmp1[67]+z[86]*tmp1[68]+z[73]*tmp1[69]+z[60]*tmp1[70]
                    +z[47]*tmp1[71]+z[34]*tmp1[72]+z[21]*tmp1[73]+z[8]*tmp1[74]+z[92]*tmp1[75]+z[79]*tmp1[76]+z[66]*tmp1[77]+z[53]*tmp1[78]+z[40]*tmp1[79]+z[27]*tmp1[80]
                    +z[14]*tmp1[81]+z[1]*tmp1[82]+z[85]*tmp1[83]+z[72]*tmp1[84]+z[59]*tmp1[85]+z[46]*tmp1[86]+z[33]*tmp1[87]+z[20]*tmp1[88]+z[7]*tmp1[89]+z[91]*tmp1[90]
                    +z[78]*tmp1[91]+z[65]*tmp1[92]+z[52]*tmp1[93]+z[39]*tmp1[94]+z[26]*tmp1[95]+z[13]*tmp1[96];
                    tab[nb_tmp3+85*nb3]=tab[nb_tmp3+85*nb3]+z[0]*tmp1[0]
                    +z[85]*tmp1[1]+z[73]*tmp1[2]+z[61]*tmp1[3]+z[49]*tmp1[4]+z[37]*tmp1[5]+z[25]*tmp1[6]+z[13]*tmp1[7]+z[1]*tmp1[8]+z[86]*tmp1[9]+z[74]*tmp1[10]
                    +z[62]*tmp1[11]+z[50]*tmp1[12]+z[38]*tmp1[13]+z[26]*tmp1[14]+z[14]*tmp1[15]+z[2]*tmp1[16]+z[87]*tmp1[17]+z[75]*tmp1[18]+z[63]*tmp1[19]+z[51]*tmp1[20]
                    +z[39]*tmp1[21]+z[27]*tmp1[22]+z[15]*tmp1[23]+z[3]*tmp1[24]+z[88]*tmp1[25]+z[76]*tmp1[26]+z[64]*tmp1[27]+z[52]*tmp1[28]+z[40]*tmp1[29]+z[28]*tmp1[30]
                    +z[16]*tmp1[31]+z[4]*tmp1[32]+z[89]*tmp1[33]+z[77]*tmp1[34]+z[65]*tmp1[35]+z[53]*tmp1[36]+z[41]*tmp1[37]+z[29]*tmp1[38]+z[17]*tmp1[39]+z[5]*tmp1[40]
                    +z[90]*tmp1[41]+z[78]*tmp1[42]+z[66]*tmp1[43]+z[54]*tmp1[44]+z[42]*tmp1[45]+z[30]*tmp1[46]+z[18]*tmp1[47]+z[6]*tmp1[48]+z[91]*tmp1[49]+z[79]*tmp1[50]
                    +z[67]*tmp1[51]+z[55]*tmp1[52]+z[43]*tmp1[53]+z[31]*tmp1[54]+z[19]*tmp1[55]+z[7]*tmp1[56]+z[92]*tmp1[57]+z[80]*tmp1[58]+z[68]*tmp1[59]+z[56]*tmp1[60]
                    +z[44]*tmp1[61]+z[32]*tmp1[62]+z[20]*tmp1[63]+z[8]*tmp1[64]+z[93]*tmp1[65]+z[81]*tmp1[66]+z[69]*tmp1[67]+z[57]*tmp1[68]+z[45]*tmp1[69]+z[33]*tmp1[70]
                    +z[21]*tmp1[71]+z[9]*tmp1[72]+z[94]*tmp1[73]+z[82]*tmp1[74]+z[70]*tmp1[75]+z[58]*tmp1[76]+z[46]*tmp1[77]+z[34]*tmp1[78]+z[22]*tmp1[79]+z[10]*tmp1[80]
                    +z[95]*tmp1[81]+z[83]*tmp1[82]+z[71]*tmp1[83]+z[59]*tmp1[84]+z[47]*tmp1[85]+z[35]*tmp1[86]+z[23]*tmp1[87]+z[11]*tmp1[88]+z[96]*tmp1[89]+z[84]*tmp1[90]
                    +z[72]*tmp1[91]+z[60]*tmp1[92]+z[48]*tmp1[93]+z[36]*tmp1[94]+z[24]*tmp1[95]+z[12]*tmp1[96];
                    tab[nb_tmp3+86*nb3]=tab[nb_tmp3+86*nb3]+z[0]*tmp1[0]
                    +z[86]*tmp1[1]+z[75]*tmp1[2]+z[64]*tmp1[3]+z[53]*tmp1[4]+z[42]*tmp1[5]+z[31]*tmp1[6]+z[20]*tmp1[7]+z[9]*tmp1[8]+z[95]*tmp1[9]+z[84]*tmp1[10]
                    +z[73]*tmp1[11]+z[62]*tmp1[12]+z[51]*tmp1[13]+z[40]*tmp1[14]+z[29]*tmp1[15]+z[18]*tmp1[16]+z[7]*tmp1[17]+z[93]*tmp1[18]+z[82]*tmp1[19]+z[71]*tmp1[20]
                    +z[60]*tmp1[21]+z[49]*tmp1[22]+z[38]*tmp1[23]+z[27]*tmp1[24]+z[16]*tmp1[25]+z[5]*tmp1[26]+z[91]*tmp1[27]+z[80]*tmp1[28]+z[69]*tmp1[29]+z[58]*tmp1[30]
                    +z[47]*tmp1[31]+z[36]*tmp1[32]+z[25]*tmp1[33]+z[14]*tmp1[34]+z[3]*tmp1[35]+z[89]*tmp1[36]+z[78]*tmp1[37]+z[67]*tmp1[38]+z[56]*tmp1[39]+z[45]*tmp1[40]
                    +z[34]*tmp1[41]+z[23]*tmp1[42]+z[12]*tmp1[43]+z[1]*tmp1[44]+z[87]*tmp1[45]+z[76]*tmp1[46]+z[65]*tmp1[47]+z[54]*tmp1[48]+z[43]*tmp1[49]+z[32]*tmp1[50]
                    +z[21]*tmp1[51]+z[10]*tmp1[52]+z[96]*tmp1[53]+z[85]*tmp1[54]+z[74]*tmp1[55]+z[63]*tmp1[56]+z[52]*tmp1[57]+z[41]*tmp1[58]+z[30]*tmp1[59]+z[19]*tmp1[60]
                    +z[8]*tmp1[61]+z[94]*tmp1[62]+z[83]*tmp1[63]+z[72]*tmp1[64]+z[61]*tmp1[65]+z[50]*tmp1[66]+z[39]*tmp1[67]+z[28]*tmp1[68]+z[17]*tmp1[69]+z[6]*tmp1[70]
                    +z[92]*tmp1[71]+z[81]*tmp1[72]+z[70]*tmp1[73]+z[59]*tmp1[74]+z[48]*tmp1[75]+z[37]*tmp1[76]+z[26]*tmp1[77]+z[15]*tmp1[78]+z[4]*tmp1[79]+z[90]*tmp1[80]
                    +z[79]*tmp1[81]+z[68]*tmp1[82]+z[57]*tmp1[83]+z[46]*tmp1[84]+z[35]*tmp1[85]+z[24]*tmp1[86]+z[13]*tmp1[87]+z[2]*tmp1[88]+z[88]*tmp1[89]+z[77]*tmp1[90]
                    +z[66]*tmp1[91]+z[55]*tmp1[92]+z[44]*tmp1[93]+z[33]*tmp1[94]+z[22]*tmp1[95]+z[11]*tmp1[96];
                    tab[nb_tmp3+87*nb3]=tab[nb_tmp3+87*nb3]+z[0]*tmp1[0]
                    +z[87]*tmp1[1]+z[77]*tmp1[2]+z[67]*tmp1[3]+z[57]*tmp1[4]+z[47]*tmp1[5]+z[37]*tmp1[6]+z[27]*tmp1[7]+z[17]*tmp1[8]+z[7]*tmp1[9]+z[94]*tmp1[10]
                    +z[84]*tmp1[11]+z[74]*tmp1[12]+z[64]*tmp1[13]+z[54]*tmp1[14]+z[44]*tmp1[15]+z[34]*tmp1[16]+z[24]*tmp1[17]+z[14]*tmp1[18]+z[4]*tmp1[19]+z[91]*tmp1[20]
                    +z[81]*tmp1[21]+z[71]*tmp1[22]+z[61]*tmp1[23]+z[51]*tmp1[24]+z[41]*tmp1[25]+z[31]*tmp1[26]+z[21]*tmp1[27]+z[11]*tmp1[28]+z[1]*tmp1[29]+z[88]*tmp1[30]
                    +z[78]*tmp1[31]+z[68]*tmp1[32]+z[58]*tmp1[33]+z[48]*tmp1[34]+z[38]*tmp1[35]+z[28]*tmp1[36]+z[18]*tmp1[37]+z[8]*tmp1[38]+z[95]*tmp1[39]+z[85]*tmp1[40]
                    +z[75]*tmp1[41]+z[65]*tmp1[42]+z[55]*tmp1[43]+z[45]*tmp1[44]+z[35]*tmp1[45]+z[25]*tmp1[46]+z[15]*tmp1[47]+z[5]*tmp1[48]+z[92]*tmp1[49]+z[82]*tmp1[50]
                    +z[72]*tmp1[51]+z[62]*tmp1[52]+z[52]*tmp1[53]+z[42]*tmp1[54]+z[32]*tmp1[55]+z[22]*tmp1[56]+z[12]*tmp1[57]+z[2]*tmp1[58]+z[89]*tmp1[59]+z[79]*tmp1[60]
                    +z[69]*tmp1[61]+z[59]*tmp1[62]+z[49]*tmp1[63]+z[39]*tmp1[64]+z[29]*tmp1[65]+z[19]*tmp1[66]+z[9]*tmp1[67]+z[96]*tmp1[68]+z[86]*tmp1[69]+z[76]*tmp1[70]
                    +z[66]*tmp1[71]+z[56]*tmp1[72]+z[46]*tmp1[73]+z[36]*tmp1[74]+z[26]*tmp1[75]+z[16]*tmp1[76]+z[6]*tmp1[77]+z[93]*tmp1[78]+z[83]*tmp1[79]+z[73]*tmp1[80]
                    +z[63]*tmp1[81]+z[53]*tmp1[82]+z[43]*tmp1[83]+z[33]*tmp1[84]+z[23]*tmp1[85]+z[13]*tmp1[86]+z[3]*tmp1[87]+z[90]*tmp1[88]+z[80]*tmp1[89]+z[70]*tmp1[90]
                    +z[60]*tmp1[91]+z[50]*tmp1[92]+z[40]*tmp1[93]+z[30]*tmp1[94]+z[20]*tmp1[95]+z[10]*tmp1[96];
                    tab[nb_tmp3+88*nb3]=tab[nb_tmp3+88*nb3]+z[0]*tmp1[0]
                    +z[88]*tmp1[1]+z[79]*tmp1[2]+z[70]*tmp1[3]+z[61]*tmp1[4]+z[52]*tmp1[5]+z[43]*tmp1[6]+z[34]*tmp1[7]+z[25]*tmp1[8]+z[16]*tmp1[9]+z[7]*tmp1[10]
                    +z[95]*tmp1[11]+z[86]*tmp1[12]+z[77]*tmp1[13]+z[68]*tmp1[14]+z[59]*tmp1[15]+z[50]*tmp1[16]+z[41]*tmp1[17]+z[32]*tmp1[18]+z[23]*tmp1[19]+z[14]*tmp1[20]
                    +z[5]*tmp1[21]+z[93]*tmp1[22]+z[84]*tmp1[23]+z[75]*tmp1[24]+z[66]*tmp1[25]+z[57]*tmp1[26]+z[48]*tmp1[27]+z[39]*tmp1[28]+z[30]*tmp1[29]+z[21]*tmp1[30]
                    +z[12]*tmp1[31]+z[3]*tmp1[32]+z[91]*tmp1[33]+z[82]*tmp1[34]+z[73]*tmp1[35]+z[64]*tmp1[36]+z[55]*tmp1[37]+z[46]*tmp1[38]+z[37]*tmp1[39]+z[28]*tmp1[40]
                    +z[19]*tmp1[41]+z[10]*tmp1[42]+z[1]*tmp1[43]+z[89]*tmp1[44]+z[80]*tmp1[45]+z[71]*tmp1[46]+z[62]*tmp1[47]+z[53]*tmp1[48]+z[44]*tmp1[49]+z[35]*tmp1[50]
                    +z[26]*tmp1[51]+z[17]*tmp1[52]+z[8]*tmp1[53]+z[96]*tmp1[54]+z[87]*tmp1[55]+z[78]*tmp1[56]+z[69]*tmp1[57]+z[60]*tmp1[58]+z[51]*tmp1[59]+z[42]*tmp1[60]
                    +z[33]*tmp1[61]+z[24]*tmp1[62]+z[15]*tmp1[63]+z[6]*tmp1[64]+z[94]*tmp1[65]+z[85]*tmp1[66]+z[76]*tmp1[67]+z[67]*tmp1[68]+z[58]*tmp1[69]+z[49]*tmp1[70]
                    +z[40]*tmp1[71]+z[31]*tmp1[72]+z[22]*tmp1[73]+z[13]*tmp1[74]+z[4]*tmp1[75]+z[92]*tmp1[76]+z[83]*tmp1[77]+z[74]*tmp1[78]+z[65]*tmp1[79]+z[56]*tmp1[80]
                    +z[47]*tmp1[81]+z[38]*tmp1[82]+z[29]*tmp1[83]+z[20]*tmp1[84]+z[11]*tmp1[85]+z[2]*tmp1[86]+z[90]*tmp1[87]+z[81]*tmp1[88]+z[72]*tmp1[89]+z[63]*tmp1[90]
                    +z[54]*tmp1[91]+z[45]*tmp1[92]+z[36]*tmp1[93]+z[27]*tmp1[94]+z[18]*tmp1[95]+z[9]*tmp1[96];
                    tab[nb_tmp3+89*nb3]=tab[nb_tmp3+89*nb3]+z[0]*tmp1[0]
                    +z[89]*tmp1[1]+z[81]*tmp1[2]+z[73]*tmp1[3]+z[65]*tmp1[4]+z[57]*tmp1[5]+z[49]*tmp1[6]+z[41]*tmp1[7]+z[33]*tmp1[8]+z[25]*tmp1[9]+z[17]*tmp1[10]
                    +z[9]*tmp1[11]+z[1]*tmp1[12]+z[90]*tmp1[13]+z[82]*tmp1[14]+z[74]*tmp1[15]+z[66]*tmp1[16]+z[58]*tmp1[17]+z[50]*tmp1[18]+z[42]*tmp1[19]+z[34]*tmp1[20]
                    +z[26]*tmp1[21]+z[18]*tmp1[22]+z[10]*tmp1[23]+z[2]*tmp1[24]+z[91]*tmp1[25]+z[83]*tmp1[26]+z[75]*tmp1[27]+z[67]*tmp1[28]+z[59]*tmp1[29]+z[51]*tmp1[30]
                    +z[43]*tmp1[31]+z[35]*tmp1[32]+z[27]*tmp1[33]+z[19]*tmp1[34]+z[11]*tmp1[35]+z[3]*tmp1[36]+z[92]*tmp1[37]+z[84]*tmp1[38]+z[76]*tmp1[39]+z[68]*tmp1[40]
                    +z[60]*tmp1[41]+z[52]*tmp1[42]+z[44]*tmp1[43]+z[36]*tmp1[44]+z[28]*tmp1[45]+z[20]*tmp1[46]+z[12]*tmp1[47]+z[4]*tmp1[48]+z[93]*tmp1[49]+z[85]*tmp1[50]
                    +z[77]*tmp1[51]+z[69]*tmp1[52]+z[61]*tmp1[53]+z[53]*tmp1[54]+z[45]*tmp1[55]+z[37]*tmp1[56]+z[29]*tmp1[57]+z[21]*tmp1[58]+z[13]*tmp1[59]+z[5]*tmp1[60]
                    +z[94]*tmp1[61]+z[86]*tmp1[62]+z[78]*tmp1[63]+z[70]*tmp1[64]+z[62]*tmp1[65]+z[54]*tmp1[66]+z[46]*tmp1[67]+z[38]*tmp1[68]+z[30]*tmp1[69]+z[22]*tmp1[70]
                    +z[14]*tmp1[71]+z[6]*tmp1[72]+z[95]*tmp1[73]+z[87]*tmp1[74]+z[79]*tmp1[75]+z[71]*tmp1[76]+z[63]*tmp1[77]+z[55]*tmp1[78]+z[47]*tmp1[79]+z[39]*tmp1[80]
                    +z[31]*tmp1[81]+z[23]*tmp1[82]+z[15]*tmp1[83]+z[7]*tmp1[84]+z[96]*tmp1[85]+z[88]*tmp1[86]+z[80]*tmp1[87]+z[72]*tmp1[88]+z[64]*tmp1[89]+z[56]*tmp1[90]
                    +z[48]*tmp1[91]+z[40]*tmp1[92]+z[32]*tmp1[93]+z[24]*tmp1[94]+z[16]*tmp1[95]+z[8]*tmp1[96];
                    tab[nb_tmp3+90*nb3]=tab[nb_tmp3+90*nb3]+z[0]*tmp1[0]
                    +z[90]*tmp1[1]+z[83]*tmp1[2]+z[76]*tmp1[3]+z[69]*tmp1[4]+z[62]*tmp1[5]+z[55]*tmp1[6]+z[48]*tmp1[7]+z[41]*tmp1[8]+z[34]*tmp1[9]+z[27]*tmp1[10]
                    +z[20]*tmp1[11]+z[13]*tmp1[12]+z[6]*tmp1[13]+z[96]*tmp1[14]+z[89]*tmp1[15]+z[82]*tmp1[16]+z[75]*tmp1[17]+z[68]*tmp1[18]+z[61]*tmp1[19]+z[54]*tmp1[20]
                    +z[47]*tmp1[21]+z[40]*tmp1[22]+z[33]*tmp1[23]+z[26]*tmp1[24]+z[19]*tmp1[25]+z[12]*tmp1[26]+z[5]*tmp1[27]+z[95]*tmp1[28]+z[88]*tmp1[29]+z[81]*tmp1[30]
                    +z[74]*tmp1[31]+z[67]*tmp1[32]+z[60]*tmp1[33]+z[53]*tmp1[34]+z[46]*tmp1[35]+z[39]*tmp1[36]+z[32]*tmp1[37]+z[25]*tmp1[38]+z[18]*tmp1[39]+z[11]*tmp1[40]
                    +z[4]*tmp1[41]+z[94]*tmp1[42]+z[87]*tmp1[43]+z[80]*tmp1[44]+z[73]*tmp1[45]+z[66]*tmp1[46]+z[59]*tmp1[47]+z[52]*tmp1[48]+z[45]*tmp1[49]+z[38]*tmp1[50]
                    +z[31]*tmp1[51]+z[24]*tmp1[52]+z[17]*tmp1[53]+z[10]*tmp1[54]+z[3]*tmp1[55]+z[93]*tmp1[56]+z[86]*tmp1[57]+z[79]*tmp1[58]+z[72]*tmp1[59]+z[65]*tmp1[60]
                    +z[58]*tmp1[61]+z[51]*tmp1[62]+z[44]*tmp1[63]+z[37]*tmp1[64]+z[30]*tmp1[65]+z[23]*tmp1[66]+z[16]*tmp1[67]+z[9]*tmp1[68]+z[2]*tmp1[69]+z[92]*tmp1[70]
                    +z[85]*tmp1[71]+z[78]*tmp1[72]+z[71]*tmp1[73]+z[64]*tmp1[74]+z[57]*tmp1[75]+z[50]*tmp1[76]+z[43]*tmp1[77]+z[36]*tmp1[78]+z[29]*tmp1[79]+z[22]*tmp1[80]
                    +z[15]*tmp1[81]+z[8]*tmp1[82]+z[1]*tmp1[83]+z[91]*tmp1[84]+z[84]*tmp1[85]+z[77]*tmp1[86]+z[70]*tmp1[87]+z[63]*tmp1[88]+z[56]*tmp1[89]+z[49]*tmp1[90]
                    +z[42]*tmp1[91]+z[35]*tmp1[92]+z[28]*tmp1[93]+z[21]*tmp1[94]+z[14]*tmp1[95]+z[7]*tmp1[96];
                    tab[nb_tmp3+91*nb3]=tab[nb_tmp3+91*nb3]+z[0]*tmp1[0]
                    +z[91]*tmp1[1]+z[85]*tmp1[2]+z[79]*tmp1[3]+z[73]*tmp1[4]+z[67]*tmp1[5]+z[61]*tmp1[6]+z[55]*tmp1[7]+z[49]*tmp1[8]+z[43]*tmp1[9]+z[37]*tmp1[10]
                    +z[31]*tmp1[11]+z[25]*tmp1[12]+z[19]*tmp1[13]+z[13]*tmp1[14]+z[7]*tmp1[15]+z[1]*tmp1[16]+z[92]*tmp1[17]+z[86]*tmp1[18]+z[80]*tmp1[19]+z[74]*tmp1[20]
                    +z[68]*tmp1[21]+z[62]*tmp1[22]+z[56]*tmp1[23]+z[50]*tmp1[24]+z[44]*tmp1[25]+z[38]*tmp1[26]+z[32]*tmp1[27]+z[26]*tmp1[28]+z[20]*tmp1[29]+z[14]*tmp1[30]
                    +z[8]*tmp1[31]+z[2]*tmp1[32]+z[93]*tmp1[33]+z[87]*tmp1[34]+z[81]*tmp1[35]+z[75]*tmp1[36]+z[69]*tmp1[37]+z[63]*tmp1[38]+z[57]*tmp1[39]+z[51]*tmp1[40]
                    +z[45]*tmp1[41]+z[39]*tmp1[42]+z[33]*tmp1[43]+z[27]*tmp1[44]+z[21]*tmp1[45]+z[15]*tmp1[46]+z[9]*tmp1[47]+z[3]*tmp1[48]+z[94]*tmp1[49]+z[88]*tmp1[50]
                    +z[82]*tmp1[51]+z[76]*tmp1[52]+z[70]*tmp1[53]+z[64]*tmp1[54]+z[58]*tmp1[55]+z[52]*tmp1[56]+z[46]*tmp1[57]+z[40]*tmp1[58]+z[34]*tmp1[59]+z[28]*tmp1[60]
                    +z[22]*tmp1[61]+z[16]*tmp1[62]+z[10]*tmp1[63]+z[4]*tmp1[64]+z[95]*tmp1[65]+z[89]*tmp1[66]+z[83]*tmp1[67]+z[77]*tmp1[68]+z[71]*tmp1[69]+z[65]*tmp1[70]
                    +z[59]*tmp1[71]+z[53]*tmp1[72]+z[47]*tmp1[73]+z[41]*tmp1[74]+z[35]*tmp1[75]+z[29]*tmp1[76]+z[23]*tmp1[77]+z[17]*tmp1[78]+z[11]*tmp1[79]+z[5]*tmp1[80]
                    +z[96]*tmp1[81]+z[90]*tmp1[82]+z[84]*tmp1[83]+z[78]*tmp1[84]+z[72]*tmp1[85]+z[66]*tmp1[86]+z[60]*tmp1[87]+z[54]*tmp1[88]+z[48]*tmp1[89]+z[42]*tmp1[90]
                    +z[36]*tmp1[91]+z[30]*tmp1[92]+z[24]*tmp1[93]+z[18]*tmp1[94]+z[12]*tmp1[95]+z[6]*tmp1[96];
                    tab[nb_tmp3+92*nb3]=tab[nb_tmp3+92*nb3]+z[0]*tmp1[0]
                    +z[92]*tmp1[1]+z[87]*tmp1[2]+z[82]*tmp1[3]+z[77]*tmp1[4]+z[72]*tmp1[5]+z[67]*tmp1[6]+z[62]*tmp1[7]+z[57]*tmp1[8]+z[52]*tmp1[9]+z[47]*tmp1[10]
                    +z[42]*tmp1[11]+z[37]*tmp1[12]+z[32]*tmp1[13]+z[27]*tmp1[14]+z[22]*tmp1[15]+z[17]*tmp1[16]+z[12]*tmp1[17]+z[7]*tmp1[18]+z[2]*tmp1[19]+z[94]*tmp1[20]
                    +z[89]*tmp1[21]+z[84]*tmp1[22]+z[79]*tmp1[23]+z[74]*tmp1[24]+z[69]*tmp1[25]+z[64]*tmp1[26]+z[59]*tmp1[27]+z[54]*tmp1[28]+z[49]*tmp1[29]+z[44]*tmp1[30]
                    +z[39]*tmp1[31]+z[34]*tmp1[32]+z[29]*tmp1[33]+z[24]*tmp1[34]+z[19]*tmp1[35]+z[14]*tmp1[36]+z[9]*tmp1[37]+z[4]*tmp1[38]+z[96]*tmp1[39]+z[91]*tmp1[40]
                    +z[86]*tmp1[41]+z[81]*tmp1[42]+z[76]*tmp1[43]+z[71]*tmp1[44]+z[66]*tmp1[45]+z[61]*tmp1[46]+z[56]*tmp1[47]+z[51]*tmp1[48]+z[46]*tmp1[49]+z[41]*tmp1[50]
                    +z[36]*tmp1[51]+z[31]*tmp1[52]+z[26]*tmp1[53]+z[21]*tmp1[54]+z[16]*tmp1[55]+z[11]*tmp1[56]+z[6]*tmp1[57]+z[1]*tmp1[58]+z[93]*tmp1[59]+z[88]*tmp1[60]
                    +z[83]*tmp1[61]+z[78]*tmp1[62]+z[73]*tmp1[63]+z[68]*tmp1[64]+z[63]*tmp1[65]+z[58]*tmp1[66]+z[53]*tmp1[67]+z[48]*tmp1[68]+z[43]*tmp1[69]+z[38]*tmp1[70]
                    +z[33]*tmp1[71]+z[28]*tmp1[72]+z[23]*tmp1[73]+z[18]*tmp1[74]+z[13]*tmp1[75]+z[8]*tmp1[76]+z[3]*tmp1[77]+z[95]*tmp1[78]+z[90]*tmp1[79]+z[85]*tmp1[80]
                    +z[80]*tmp1[81]+z[75]*tmp1[82]+z[70]*tmp1[83]+z[65]*tmp1[84]+z[60]*tmp1[85]+z[55]*tmp1[86]+z[50]*tmp1[87]+z[45]*tmp1[88]+z[40]*tmp1[89]+z[35]*tmp1[90]
                    +z[30]*tmp1[91]+z[25]*tmp1[92]+z[20]*tmp1[93]+z[15]*tmp1[94]+z[10]*tmp1[95]+z[5]*tmp1[96];
                    tab[nb_tmp3+93*nb3]=tab[nb_tmp3+93*nb3]+z[0]*tmp1[0]
                    +z[93]*tmp1[1]+z[89]*tmp1[2]+z[85]*tmp1[3]+z[81]*tmp1[4]+z[77]*tmp1[5]+z[73]*tmp1[6]+z[69]*tmp1[7]+z[65]*tmp1[8]+z[61]*tmp1[9]+z[57]*tmp1[10]
                    +z[53]*tmp1[11]+z[49]*tmp1[12]+z[45]*tmp1[13]+z[41]*tmp1[14]+z[37]*tmp1[15]+z[33]*tmp1[16]+z[29]*tmp1[17]+z[25]*tmp1[18]+z[21]*tmp1[19]+z[17]*tmp1[20]
                    +z[13]*tmp1[21]+z[9]*tmp1[22]+z[5]*tmp1[23]+z[1]*tmp1[24]+z[94]*tmp1[25]+z[90]*tmp1[26]+z[86]*tmp1[27]+z[82]*tmp1[28]+z[78]*tmp1[29]+z[74]*tmp1[30]
                    +z[70]*tmp1[31]+z[66]*tmp1[32]+z[62]*tmp1[33]+z[58]*tmp1[34]+z[54]*tmp1[35]+z[50]*tmp1[36]+z[46]*tmp1[37]+z[42]*tmp1[38]+z[38]*tmp1[39]+z[34]*tmp1[40]
                    +z[30]*tmp1[41]+z[26]*tmp1[42]+z[22]*tmp1[43]+z[18]*tmp1[44]+z[14]*tmp1[45]+z[10]*tmp1[46]+z[6]*tmp1[47]+z[2]*tmp1[48]+z[95]*tmp1[49]+z[91]*tmp1[50]
                    +z[87]*tmp1[51]+z[83]*tmp1[52]+z[79]*tmp1[53]+z[75]*tmp1[54]+z[71]*tmp1[55]+z[67]*tmp1[56]+z[63]*tmp1[57]+z[59]*tmp1[58]+z[55]*tmp1[59]+z[51]*tmp1[60]
                    +z[47]*tmp1[61]+z[43]*tmp1[62]+z[39]*tmp1[63]+z[35]*tmp1[64]+z[31]*tmp1[65]+z[27]*tmp1[66]+z[23]*tmp1[67]+z[19]*tmp1[68]+z[15]*tmp1[69]+z[11]*tmp1[70]
                    +z[7]*tmp1[71]+z[3]*tmp1[72]+z[96]*tmp1[73]+z[92]*tmp1[74]+z[88]*tmp1[75]+z[84]*tmp1[76]+z[80]*tmp1[77]+z[76]*tmp1[78]+z[72]*tmp1[79]+z[68]*tmp1[80]
                    +z[64]*tmp1[81]+z[60]*tmp1[82]+z[56]*tmp1[83]+z[52]*tmp1[84]+z[48]*tmp1[85]+z[44]*tmp1[86]+z[40]*tmp1[87]+z[36]*tmp1[88]+z[32]*tmp1[89]+z[28]*tmp1[90]
                    +z[24]*tmp1[91]+z[20]*tmp1[92]+z[16]*tmp1[93]+z[12]*tmp1[94]+z[8]*tmp1[95]+z[4]*tmp1[96];
                    tab[nb_tmp3+94*nb3]=tab[nb_tmp3+94*nb3]+z[0]*tmp1[0]
                    +z[94]*tmp1[1]+z[91]*tmp1[2]+z[88]*tmp1[3]+z[85]*tmp1[4]+z[82]*tmp1[5]+z[79]*tmp1[6]+z[76]*tmp1[7]+z[73]*tmp1[8]+z[70]*tmp1[9]+z[67]*tmp1[10]
                    +z[64]*tmp1[11]+z[61]*tmp1[12]+z[58]*tmp1[13]+z[55]*tmp1[14]+z[52]*tmp1[15]+z[49]*tmp1[16]+z[46]*tmp1[17]+z[43]*tmp1[18]+z[40]*tmp1[19]+z[37]*tmp1[20]
                    +z[34]*tmp1[21]+z[31]*tmp1[22]+z[28]*tmp1[23]+z[25]*tmp1[24]+z[22]*tmp1[25]+z[19]*tmp1[26]+z[16]*tmp1[27]+z[13]*tmp1[28]+z[10]*tmp1[29]+z[7]*tmp1[30]
                    +z[4]*tmp1[31]+z[1]*tmp1[32]+z[95]*tmp1[33]+z[92]*tmp1[34]+z[89]*tmp1[35]+z[86]*tmp1[36]+z[83]*tmp1[37]+z[80]*tmp1[38]+z[77]*tmp1[39]+z[74]*tmp1[40]
                    +z[71]*tmp1[41]+z[68]*tmp1[42]+z[65]*tmp1[43]+z[62]*tmp1[44]+z[59]*tmp1[45]+z[56]*tmp1[46]+z[53]*tmp1[47]+z[50]*tmp1[48]+z[47]*tmp1[49]+z[44]*tmp1[50]
                    +z[41]*tmp1[51]+z[38]*tmp1[52]+z[35]*tmp1[53]+z[32]*tmp1[54]+z[29]*tmp1[55]+z[26]*tmp1[56]+z[23]*tmp1[57]+z[20]*tmp1[58]+z[17]*tmp1[59]+z[14]*tmp1[60]
                    +z[11]*tmp1[61]+z[8]*tmp1[62]+z[5]*tmp1[63]+z[2]*tmp1[64]+z[96]*tmp1[65]+z[93]*tmp1[66]+z[90]*tmp1[67]+z[87]*tmp1[68]+z[84]*tmp1[69]+z[81]*tmp1[70]
                    +z[78]*tmp1[71]+z[75]*tmp1[72]+z[72]*tmp1[73]+z[69]*tmp1[74]+z[66]*tmp1[75]+z[63]*tmp1[76]+z[60]*tmp1[77]+z[57]*tmp1[78]+z[54]*tmp1[79]+z[51]*tmp1[80]
                    +z[48]*tmp1[81]+z[45]*tmp1[82]+z[42]*tmp1[83]+z[39]*tmp1[84]+z[36]*tmp1[85]+z[33]*tmp1[86]+z[30]*tmp1[87]+z[27]*tmp1[88]+z[24]*tmp1[89]+z[21]*tmp1[90]
                    +z[18]*tmp1[91]+z[15]*tmp1[92]+z[12]*tmp1[93]+z[9]*tmp1[94]+z[6]*tmp1[95]+z[3]*tmp1[96];
                    tab[nb_tmp3+95*nb3]=tab[nb_tmp3+95*nb3]+z[0]*tmp1[0]
                    +z[95]*tmp1[1]+z[93]*tmp1[2]+z[91]*tmp1[3]+z[89]*tmp1[4]+z[87]*tmp1[5]+z[85]*tmp1[6]+z[83]*tmp1[7]+z[81]*tmp1[8]+z[79]*tmp1[9]+z[77]*tmp1[10]
                    +z[75]*tmp1[11]+z[73]*tmp1[12]+z[71]*tmp1[13]+z[69]*tmp1[14]+z[67]*tmp1[15]+z[65]*tmp1[16]+z[63]*tmp1[17]+z[61]*tmp1[18]+z[59]*tmp1[19]+z[57]*tmp1[20]
                    +z[55]*tmp1[21]+z[53]*tmp1[22]+z[51]*tmp1[23]+z[49]*tmp1[24]+z[47]*tmp1[25]+z[45]*tmp1[26]+z[43]*tmp1[27]+z[41]*tmp1[28]+z[39]*tmp1[29]+z[37]*tmp1[30]
                    +z[35]*tmp1[31]+z[33]*tmp1[32]+z[31]*tmp1[33]+z[29]*tmp1[34]+z[27]*tmp1[35]+z[25]*tmp1[36]+z[23]*tmp1[37]+z[21]*tmp1[38]+z[19]*tmp1[39]+z[17]*tmp1[40]
                    +z[15]*tmp1[41]+z[13]*tmp1[42]+z[11]*tmp1[43]+z[9]*tmp1[44]+z[7]*tmp1[45]+z[5]*tmp1[46]+z[3]*tmp1[47]+z[1]*tmp1[48]+z[96]*tmp1[49]+z[94]*tmp1[50]
                    +z[92]*tmp1[51]+z[90]*tmp1[52]+z[88]*tmp1[53]+z[86]*tmp1[54]+z[84]*tmp1[55]+z[82]*tmp1[56]+z[80]*tmp1[57]+z[78]*tmp1[58]+z[76]*tmp1[59]+z[74]*tmp1[60]
                    +z[72]*tmp1[61]+z[70]*tmp1[62]+z[68]*tmp1[63]+z[66]*tmp1[64]+z[64]*tmp1[65]+z[62]*tmp1[66]+z[60]*tmp1[67]+z[58]*tmp1[68]+z[56]*tmp1[69]+z[54]*tmp1[70]
                    +z[52]*tmp1[71]+z[50]*tmp1[72]+z[48]*tmp1[73]+z[46]*tmp1[74]+z[44]*tmp1[75]+z[42]*tmp1[76]+z[40]*tmp1[77]+z[38]*tmp1[78]+z[36]*tmp1[79]+z[34]*tmp1[80]
                    +z[32]*tmp1[81]+z[30]*tmp1[82]+z[28]*tmp1[83]+z[26]*tmp1[84]+z[24]*tmp1[85]+z[22]*tmp1[86]+z[20]*tmp1[87]+z[18]*tmp1[88]+z[16]*tmp1[89]+z[14]*tmp1[90]
                    +z[12]*tmp1[91]+z[10]*tmp1[92]+z[8]*tmp1[93]+z[6]*tmp1[94]+z[4]*tmp1[95]+z[2]*tmp1[96];
                    tab[nb_tmp3+96*nb3]=tab[nb_tmp3+96*nb3]+z[0]*tmp1[0]
                    +z[96]*tmp1[1]+z[95]*tmp1[2]+z[94]*tmp1[3]+z[93]*tmp1[4]+z[92]*tmp1[5]+z[91]*tmp1[6]+z[90]*tmp1[7]+z[89]*tmp1[8]+z[88]*tmp1[9]+z[87]*tmp1[10]
                    +z[86]*tmp1[11]+z[85]*tmp1[12]+z[84]*tmp1[13]+z[83]*tmp1[14]+z[82]*tmp1[15]+z[81]*tmp1[16]+z[80]*tmp1[17]+z[79]*tmp1[18]+z[78]*tmp1[19]+z[77]*tmp1[20]
                    +z[76]*tmp1[21]+z[75]*tmp1[22]+z[74]*tmp1[23]+z[73]*tmp1[24]+z[72]*tmp1[25]+z[71]*tmp1[26]+z[70]*tmp1[27]+z[69]*tmp1[28]+z[68]*tmp1[29]+z[67]*tmp1[30]
                    +z[66]*tmp1[31]+z[65]*tmp1[32]+z[64]*tmp1[33]+z[63]*tmp1[34]+z[62]*tmp1[35]+z[61]*tmp1[36]+z[60]*tmp1[37]+z[59]*tmp1[38]+z[58]*tmp1[39]+z[57]*tmp1[40]
                    +z[56]*tmp1[41]+z[55]*tmp1[42]+z[54]*tmp1[43]+z[53]*tmp1[44]+z[52]*tmp1[45]+z[51]*tmp1[46]+z[50]*tmp1[47]+z[49]*tmp1[48]+z[48]*tmp1[49]+z[47]*tmp1[50]
                    +z[46]*tmp1[51]+z[45]*tmp1[52]+z[44]*tmp1[53]+z[43]*tmp1[54]+z[42]*tmp1[55]+z[41]*tmp1[56]+z[40]*tmp1[57]+z[39]*tmp1[58]+z[38]*tmp1[59]+z[37]*tmp1[60]
                    +z[36]*tmp1[61]+z[35]*tmp1[62]+z[34]*tmp1[63]+z[33]*tmp1[64]+z[32]*tmp1[65]+z[31]*tmp1[66]+z[30]*tmp1[67]+z[29]*tmp1[68]+z[28]*tmp1[69]+z[27]*tmp1[70]
                    +z[26]*tmp1[71]+z[25]*tmp1[72]+z[24]*tmp1[73]+z[23]*tmp1[74]+z[22]*tmp1[75]+z[21]*tmp1[76]+z[20]*tmp1[77]+z[19]*tmp1[78]+z[18]*tmp1[79]+z[17]*tmp1[80]
                    +z[16]*tmp1[81]+z[15]*tmp1[82]+z[14]*tmp1[83]+z[13]*tmp1[84]+z[12]*tmp1[85]+z[11]*tmp1[86]+z[10]*tmp1[87]+z[9]*tmp1[88]+z[8]*tmp1[89]+z[7]*tmp1[90]
                    +z[6]*tmp1[91]+z[5]*tmp1[92]+z[4]*tmp1[93]+z[3]*tmp1[94]+z[2]*tmp1[95]+z[1]*tmp1[96];



                }

            }
        }
    }
///////////////////////////////////////////////////////////////////////


	void fun_inverse_fourier_transform_FFT_radix_2_stg_first(std::complex<double> tab[],int stg_first,double fi2,std::complex<double> z[])
    {
        std::complex<double> tmp1,tmp2,tmp40,tmp50;
        std::complex<double>  w[2]={{1,0}};

        w[0].real(cos(0+fi2));
        w[0].imag(sin(0+fi2));
        w[1].real(cos(0+fi2));
        w[1].imag(sin(0+fi2));

        for(int i=0;i<stg_first;i=i+1)
        {
            tmp1=w[0]*tab[i+0*stg_first];
            tmp2=w[1]*tab[i+1*stg_first];

            tmp40=z[0]*(tmp1+tmp2);

            tab[i+0*stg_first]=tmp40;
            tab[i+1*stg_first]=z[0]*tmp1+z[1]*tmp2;
        }
    }
///////////////////////////////////////////////////////////////


   void fun_inverse_fourier_transform_FFT_radix_3_stg_first(std::complex<double> tab[],int stg_first,double fi2,std::complex<double> z[])
    {
        std::complex<double> tmp1,tmp2,tmp3,tmp4,tmp5;
        std::complex<double>  w[3]={{1,0}};

        w[0].real(cos(0+fi2));
        w[0].imag(sin(0+fi2));
        w[1].real(cos(0+fi2));
        w[1].imag(sin(0+fi2));
        w[2].real(cos(0+fi2));
        w[2].imag(sin(0+fi2));

        for(int i=0;i<stg_first;i=i+1)
        {
          tmp1=w[0]*tab[i+0*stg_first];
          tmp2=w[1]*tab[i+1*stg_first];
          tmp3=w[2]*tab[i+2*stg_first];

          tmp4=z[0]*tmp1;
          tmp5=z[0]*(tmp2+tmp3);

         //radix-3
          tab[i+0*stg_first]  = tmp4+tmp5;
          tab[i+1*stg_first]  = tmp4+z[1]*tmp2+z[2]*tmp3;
          tab[i+2*stg_first]  = tmp4+z[2]*tmp2+z[1]*tmp3;
        }
    }
///////////////////////////////////////////////////////////////


   void fun_inverse_fourier_transform_FFT_radix_4_stg_first(std::complex<double> tab[],int stg_first,double fi2,std::complex<double> z[])
    {
        std::complex<double> tmp1,tmp2,tmp3,tmp4;
        std::complex<double> tmp101,tmp102,tmp103,tmp104;
        std::complex<double> tmp111,tmp112,tmp113,tmp114;
        std::complex<double>  w[4]={{1,0}};

        w[0].real(cos(0+fi2));
        w[0].imag(sin(0+fi2));
        w[1].real(cos(0+fi2));
        w[1].imag(sin(0+fi2));
        w[2].real(cos(0+fi2));
        w[2].imag(sin(0+fi2));
        w[3].real(cos(0+fi2));
        w[3].imag(sin(0+fi2));
        for(int i=0;i<stg_first;i=i+1)
        {
            tmp1=w[0]*tab[i+0*stg_first];
            tmp2=w[1]*tab[i+1*stg_first];
            tmp3=w[2]*tab[i+2*stg_first];
            tmp4=w[3]*tab[i+3*stg_first];


            tmp101=tmp2-tmp4;
            tmp102=tmp2+tmp4;
            tmp103=tmp1-tmp3;
            tmp104=tmp1+tmp3;

            tmp111=z[0]*(tmp104+tmp102);
            tmp112=z[0]*tmp103;
            tmp113=z[0]*(tmp104-tmp102);
            tmp114=z[0]*tmp103;

            //radix-4
            tab[i+0*stg_first]   =tmp111;
            tab[i+1*stg_first]   =tmp112+z[1]*tmp101;
            tab[i+2*stg_first]   =tmp113;
            tab[i+3*stg_first]   =tmp114-z[1]*tmp101;
        }
    }
/////////////////////////////////////////////////////////////////


void fun_inverse_fourier_transform_FFT_radix_5_stg_first(std::complex<double> tab[],int stg_first,double fi2,std::complex<double> z_rx5[])
  {
        std::complex<double> tmp1,tmp2,tmp3,tmp4,tmp5,tmp33,tmp34;
        std::complex<double>  w[5]={{1,0}};

        w[0].real(cos(0+fi2));
        w[0].imag(sin(0+fi2));
        w[1].real(cos(0+fi2));
        w[1].imag(sin(0+fi2));
        w[2].real(cos(0+fi2));
        w[2].imag(sin(0+fi2));
        w[3].real(cos(0+fi2));
        w[3].imag(sin(0+fi2));
        w[4].real(cos(0+fi2));
        w[4].imag(sin(0+fi2));

        for(int i=0;i<stg_first;i++)
        {
          tmp1=w[0]*tab[i+0*stg_first];
          tmp2=w[1]*tab[i+1*stg_first];
          tmp3=w[2]*tab[i+2*stg_first];
          tmp4=w[4]*tab[i+3*stg_first];
          tmp5=w[4]*tab[i+4*stg_first];

          tmp33=z_rx5[0]*tmp1;
          tmp34=z_rx5[0]*(tmp1+tmp2+tmp3+tmp4+tmp5);
         //radix-5
          tab[i+0*stg_first] =tmp34;
          tab[i+1*stg_first] =tmp33+z_rx5[1]*tmp2+z_rx5[2]*tmp3+z_rx5[3]*tmp4+z_rx5[4]*tmp5;
          tab[i+2*stg_first] =tmp33+z_rx5[2]*tmp2+z_rx5[4]*tmp3+z_rx5[1]*tmp4+z_rx5[3]*tmp5;
          tab[i+3*stg_first] =tmp33+z_rx5[3]*tmp2+z_rx5[1]*tmp3+z_rx5[4]*tmp4+z_rx5[2]*tmp5;
          tab[i+4*stg_first] =tmp33+z_rx5[4]*tmp2+z_rx5[3]*tmp3+z_rx5[2]*tmp4+z_rx5[1]*tmp5;
        }
  }
/////////////////////////////////////////


   void fun_inverse_fourier_transform_FFT_radix_7_stg_first(std::complex<double> tab[],int stg_first,double fi2,std::complex<double> z[])
  {
        std::complex<double> tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp10,tmp20;
        std::complex<double>  w[7]={{1,0}};


        w[0].real(cos(0+fi2));
        w[0].imag(sin(0+fi2));
        w[1].real(cos(0+fi2));
        w[1].imag(sin(0+fi2));
        w[2].real(cos(0+fi2));
        w[2].imag(sin(0+fi2));
        w[3].real(cos(0+fi2));
        w[3].imag(sin(0+fi2));
        w[4].real(cos(0+fi2));
        w[4].imag(sin(0+fi2));
        w[5].real(cos(0+fi2));
        w[5].imag(sin(0+fi2));
        w[6].real(cos(0+fi2));
        w[6].imag(sin(0+fi2));

        for(int i=0;i<stg_first;i=i+1)
        {
          tmp1=w[0]*tab[i+0*stg_first];
          tmp2=w[1]*tab[i+1*stg_first];
          tmp3=w[2]*tab[i+2*stg_first];
          tmp4=w[3]*tab[i+3*stg_first];
          tmp5=w[4]*tab[i+4*stg_first];
          tmp6=w[5]*tab[i+5*stg_first];
          tmp7=w[6]*tab[i+6*stg_first];

          tmp10=z[0]*tmp1;
          tmp20=z[0]*(tmp1+tmp2+tmp3+tmp4+tmp5+tmp6+tmp7);

         //radix-7
          tab[i+0*stg_first] =tmp20;
          tab[i+1*stg_first] =tmp10+z[1]*tmp2+z[2]*tmp3+z[3]*tmp4+z[4]*tmp5+z[5]*tmp6+z[6]*tmp7;
          tab[i+2*stg_first] =tmp10+z[2]*tmp2+z[4]*tmp3+z[6]*tmp4+z[1]*tmp5+z[3]*tmp6+z[5]*tmp7;
          tab[i+3*stg_first] =tmp10+z[3]*tmp2+z[6]*tmp3+z[2]*tmp4+z[5]*tmp5+z[1]*tmp6+z[4]*tmp7;
          tab[i+4*stg_first] =tmp10+z[4]*tmp2+z[1]*tmp3+z[5]*tmp4+z[2]*tmp5+z[6]*tmp6+z[3]*tmp7;
          tab[i+5*stg_first] =tmp10+z[5]*tmp2+z[3]*tmp3+z[1]*tmp4+z[6]*tmp5+z[4]*tmp6+z[2]*tmp7;
          tab[i+6*stg_first] =tmp10+z[6]*tmp2+z[5]*tmp3+z[4]*tmp4+z[3]*tmp5+z[2]*tmp6+z[1]*tmp7;
        }

    }
///////////////////////////////////////////////////////////////////////


 void fun_inverse_fourier_transform_FFT_radix_11_stg_first(std::complex<double> tab[],int stg_first,double fi2,std::complex<double> z[])
  {
        std::complex<double> tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,tmp7,tmp8,tmp9,tmp10,tmp11,tmp33,tmp34;
        std::complex<double>  w[11]={{1,0}};

        w[0].real(cos(0+fi2));
        w[0].imag(sin(0+fi2));
        w[1].real(cos(0+fi2));
        w[1].imag(sin(0+fi2));
        w[2].real(cos(0+fi2));
        w[2].imag(sin(0+fi2));
        w[3].real(cos(0+fi2));
        w[3].imag(sin(0+fi2));
        w[4].real(cos(0+fi2));
        w[4].imag(sin(0+fi2));
        w[5].real(cos(0+fi2));
        w[5].imag(sin(0+fi2));
        w[6].real(cos(0+fi2));
        w[6].imag(sin(0+fi2));
        w[7].real(cos(0+fi2));
        w[7].imag(sin(0+fi2));
        w[8].real(cos(0+fi2));
        w[8].imag(sin(0+fi2));
        w[9].real(cos(0+fi2));
        w[9].imag(sin(0+fi2));
        w[10].real(cos(0+fi2));
        w[10].imag(sin(0+fi2));


        for(int i=0;i<stg_first;i=i+1)
        {
          tmp1=w[0]*tab[i+0*stg_first];
          tmp2=w[1]*tab[i+1*stg_first];
          tmp3=w[2]*tab[i+2*stg_first];
          tmp4=w[3]*tab[i+3*stg_first];
          tmp5=w[4]*tab[i+4*stg_first];
          tmp6=w[5]*tab[i+5*stg_first];
          tmp7=w[6]*tab[i+6*stg_first];
          tmp8=w[7]*tab[i+7*stg_first];
          tmp9=w[8]*tab[i+8*stg_first];
          tmp10=w[9]*tab[i+9*stg_first];
          tmp11=w[10]*tab[i+10*stg_first];

          tmp33=z[0]*tmp1;
          tmp34=z[0]*(tmp1+tmp2+tmp3+tmp4+tmp5+tmp6+tmp7+tmp8+tmp9+tmp10+tmp11);

         //radix-11
          tab[i+0*stg_first] =tmp34;
          tab[i+1*stg_first] =tmp33+z[1]*tmp2+z[2]*tmp3+z[3]*tmp4+z[4]*tmp5+z[5]*tmp6+z[6]*tmp7+z[7]*tmp8+z[8]*tmp9+z[9]*tmp10+z[10]*tmp11;
          tab[i+2*stg_first] =tmp33+z[2]*tmp2+z[4]*tmp3+z[6]*tmp4+z[8]*tmp5+z[10]*tmp6+z[1]*tmp7+z[3]*tmp8+z[5]*tmp9+z[7]*tmp10+z[9]*tmp11;
          tab[i+3*stg_first] =tmp33+z[3]*tmp2+z[6]*tmp3+z[9]*tmp4+z[1]*tmp5+z[4]*tmp6+z[7]*tmp7+z[10]*tmp8+z[2]*tmp9+z[5]*tmp10+z[8]*tmp11;
          tab[i+4*stg_first] =tmp33+z[4]*tmp2+z[8]*tmp3+z[1]*tmp4+z[5]*tmp5+z[9]*tmp6+z[2]*tmp7+z[6]*tmp8+z[10]*tmp9+z[3]*tmp10+z[7]*tmp11;
          tab[i+5*stg_first] =tmp33+z[5]*tmp2+z[10]*tmp3+z[4]*tmp4+z[9]*tmp5+z[3]*tmp6+z[8]*tmp7+z[2]*tmp8+z[7]*tmp9+z[1]*tmp10+z[6]*tmp11;
          tab[i+6*stg_first] =tmp33+z[6]*tmp2+z[1]*tmp3+z[7]*tmp4+z[2]*tmp5+z[8]*tmp6+z[3]*tmp7+z[9]*tmp8+z[4]*tmp9+z[10]*tmp10+z[5]*tmp11;
          tab[i+7*stg_first] =tmp33+z[7]*tmp2+z[3]*tmp3+z[10]*tmp4+z[6]*tmp5+z[2]*tmp6+z[9]*tmp7+z[5]*tmp8+z[1]*tmp9+z[8]*tmp10+z[4]*tmp11;
          tab[i+8*stg_first] =tmp33+z[8]*tmp2+z[5]*tmp3+z[2]*tmp4+z[10]*tmp5+z[7]*tmp6+z[4]*tmp7+z[1]*tmp8+z[9]*tmp9+z[6]*tmp10+z[3]*tmp11;
          tab[i+9*stg_first] =tmp33+z[9]*tmp2+z[7]*tmp3+z[5]*tmp4+z[3]*tmp5+z[1]*tmp6+z[10]*tmp7+z[8]*tmp8+z[6]*tmp9+z[4]*tmp10+z[2]*tmp11;
          tab[i+10*stg_first] =tmp33+z[10]*tmp2+z[9]*tmp3+z[8]*tmp4+z[7]*tmp5+z[6]*tmp6+z[5]*tmp7+z[4]*tmp8+z[3]*tmp9+z[2]*tmp10+z[1]*tmp11;
        }
    }
///////////////////////////////////////////////////


 void fun_inverse_fourier_transform_FFT_radix_97_stg_first(std::complex<double> tab[],int stg_first,double fi2,std::complex<double> z[])
  {
        std::complex<double> tmp1[97];
        std::complex<double>  w[97]={{1,0}};


	for(int j9=0;j9<97;j9++)
	{
        w[j9].real(cos(0+fi2));
        w[j9].imag(sin(0+fi2));
	}

        for(int i=0;i<stg_first;i=i+1)
        {
          for(int j9=0;j9<97;j9++)
          {
                tmp1[j9]=w[j9]*tab[i+j9*stg_first];
          }

          //radix-97
          for(int i9=0;i9<97;i9++)
          {
            for(int j9=0;j9<97;j9++)
            {
                tab[i+i9*stg_first].real(0);
                tab[i+i9*stg_first].imag(0);
            }
          }

            tab[i+0*stg_first]=tab[i+0*stg_first]+z[0]*tmp1[0]
            +z[0]*tmp1[1]+z[0]*tmp1[2]+z[0]*tmp1[3]+z[0]*tmp1[4]+z[0]*tmp1[5]+z[0]*tmp1[6]+z[0]*tmp1[7]+z[0]*tmp1[8]+z[0]*tmp1[9]+z[0]*tmp1[10]
            +z[0]*tmp1[11]+z[0]*tmp1[12]+z[0]*tmp1[13]+z[0]*tmp1[14]+z[0]*tmp1[15]+z[0]*tmp1[16]+z[0]*tmp1[17]+z[0]*tmp1[18]+z[0]*tmp1[19]+z[0]*tmp1[20]
            +z[0]*tmp1[21]+z[0]*tmp1[22]+z[0]*tmp1[23]+z[0]*tmp1[24]+z[0]*tmp1[25]+z[0]*tmp1[26]+z[0]*tmp1[27]+z[0]*tmp1[28]+z[0]*tmp1[29]+z[0]*tmp1[30]
            +z[0]*tmp1[31]+z[0]*tmp1[32]+z[0]*tmp1[33]+z[0]*tmp1[34]+z[0]*tmp1[35]+z[0]*tmp1[36]+z[0]*tmp1[37]+z[0]*tmp1[38]+z[0]*tmp1[39]+z[0]*tmp1[40]
            +z[0]*tmp1[41]+z[0]*tmp1[42]+z[0]*tmp1[43]+z[0]*tmp1[44]+z[0]*tmp1[45]+z[0]*tmp1[46]+z[0]*tmp1[47]+z[0]*tmp1[48]+z[0]*tmp1[49]+z[0]*tmp1[50]
            +z[0]*tmp1[51]+z[0]*tmp1[52]+z[0]*tmp1[53]+z[0]*tmp1[54]+z[0]*tmp1[55]+z[0]*tmp1[56]+z[0]*tmp1[57]+z[0]*tmp1[58]+z[0]*tmp1[59]+z[0]*tmp1[60]
            +z[0]*tmp1[61]+z[0]*tmp1[62]+z[0]*tmp1[63]+z[0]*tmp1[64]+z[0]*tmp1[65]+z[0]*tmp1[66]+z[0]*tmp1[67]+z[0]*tmp1[68]+z[0]*tmp1[69]+z[0]*tmp1[70]
            +z[0]*tmp1[71]+z[0]*tmp1[72]+z[0]*tmp1[73]+z[0]*tmp1[74]+z[0]*tmp1[75]+z[0]*tmp1[76]+z[0]*tmp1[77]+z[0]*tmp1[78]+z[0]*tmp1[79]+z[0]*tmp1[80]
            +z[0]*tmp1[81]+z[0]*tmp1[82]+z[0]*tmp1[83]+z[0]*tmp1[84]+z[0]*tmp1[85]+z[0]*tmp1[86]+z[0]*tmp1[87]+z[0]*tmp1[88]+z[0]*tmp1[89]+z[0]*tmp1[90]
            +z[0]*tmp1[91]+z[0]*tmp1[92]+z[0]*tmp1[93]+z[0]*tmp1[94]+z[0]*tmp1[95]+z[0]*tmp1[96];
            tab[i+1*stg_first]=tab[i+1*stg_first]+z[0]*tmp1[0]
            +z[1]*tmp1[1]+z[2]*tmp1[2]+z[3]*tmp1[3]+z[4]*tmp1[4]+z[5]*tmp1[5]+z[6]*tmp1[6]+z[7]*tmp1[7]+z[8]*tmp1[8]+z[9]*tmp1[9]+z[10]*tmp1[10]
            +z[11]*tmp1[11]+z[12]*tmp1[12]+z[13]*tmp1[13]+z[14]*tmp1[14]+z[15]*tmp1[15]+z[16]*tmp1[16]+z[17]*tmp1[17]+z[18]*tmp1[18]+z[19]*tmp1[19]+z[20]*tmp1[20]
            +z[21]*tmp1[21]+z[22]*tmp1[22]+z[23]*tmp1[23]+z[24]*tmp1[24]+z[25]*tmp1[25]+z[26]*tmp1[26]+z[27]*tmp1[27]+z[28]*tmp1[28]+z[29]*tmp1[29]+z[30]*tmp1[30]
            +z[31]*tmp1[31]+z[32]*tmp1[32]+z[33]*tmp1[33]+z[34]*tmp1[34]+z[35]*tmp1[35]+z[36]*tmp1[36]+z[37]*tmp1[37]+z[38]*tmp1[38]+z[39]*tmp1[39]+z[40]*tmp1[40]
            +z[41]*tmp1[41]+z[42]*tmp1[42]+z[43]*tmp1[43]+z[44]*tmp1[44]+z[45]*tmp1[45]+z[46]*tmp1[46]+z[47]*tmp1[47]+z[48]*tmp1[48]+z[49]*tmp1[49]+z[50]*tmp1[50]
            +z[51]*tmp1[51]+z[52]*tmp1[52]+z[53]*tmp1[53]+z[54]*tmp1[54]+z[55]*tmp1[55]+z[56]*tmp1[56]+z[57]*tmp1[57]+z[58]*tmp1[58]+z[59]*tmp1[59]+z[60]*tmp1[60]
            +z[61]*tmp1[61]+z[62]*tmp1[62]+z[63]*tmp1[63]+z[64]*tmp1[64]+z[65]*tmp1[65]+z[66]*tmp1[66]+z[67]*tmp1[67]+z[68]*tmp1[68]+z[69]*tmp1[69]+z[70]*tmp1[70]
            +z[71]*tmp1[71]+z[72]*tmp1[72]+z[73]*tmp1[73]+z[74]*tmp1[74]+z[75]*tmp1[75]+z[76]*tmp1[76]+z[77]*tmp1[77]+z[78]*tmp1[78]+z[79]*tmp1[79]+z[80]*tmp1[80]
            +z[81]*tmp1[81]+z[82]*tmp1[82]+z[83]*tmp1[83]+z[84]*tmp1[84]+z[85]*tmp1[85]+z[86]*tmp1[86]+z[87]*tmp1[87]+z[88]*tmp1[88]+z[89]*tmp1[89]+z[90]*tmp1[90]
            +z[91]*tmp1[91]+z[92]*tmp1[92]+z[93]*tmp1[93]+z[94]*tmp1[94]+z[95]*tmp1[95]+z[96]*tmp1[96];
            tab[i+2*stg_first]=tab[i+2*stg_first]+z[0]*tmp1[0]
            +z[2]*tmp1[1]+z[4]*tmp1[2]+z[6]*tmp1[3]+z[8]*tmp1[4]+z[10]*tmp1[5]+z[12]*tmp1[6]+z[14]*tmp1[7]+z[16]*tmp1[8]+z[18]*tmp1[9]+z[20]*tmp1[10]
            +z[22]*tmp1[11]+z[24]*tmp1[12]+z[26]*tmp1[13]+z[28]*tmp1[14]+z[30]*tmp1[15]+z[32]*tmp1[16]+z[34]*tmp1[17]+z[36]*tmp1[18]+z[38]*tmp1[19]+z[40]*tmp1[20]
            +z[42]*tmp1[21]+z[44]*tmp1[22]+z[46]*tmp1[23]+z[48]*tmp1[24]+z[50]*tmp1[25]+z[52]*tmp1[26]+z[54]*tmp1[27]+z[56]*tmp1[28]+z[58]*tmp1[29]+z[60]*tmp1[30]
            +z[62]*tmp1[31]+z[64]*tmp1[32]+z[66]*tmp1[33]+z[68]*tmp1[34]+z[70]*tmp1[35]+z[72]*tmp1[36]+z[74]*tmp1[37]+z[76]*tmp1[38]+z[78]*tmp1[39]+z[80]*tmp1[40]
            +z[82]*tmp1[41]+z[84]*tmp1[42]+z[86]*tmp1[43]+z[88]*tmp1[44]+z[90]*tmp1[45]+z[92]*tmp1[46]+z[94]*tmp1[47]+z[96]*tmp1[48]+z[1]*tmp1[49]+z[3]*tmp1[50]
            +z[5]*tmp1[51]+z[7]*tmp1[52]+z[9]*tmp1[53]+z[11]*tmp1[54]+z[13]*tmp1[55]+z[15]*tmp1[56]+z[17]*tmp1[57]+z[19]*tmp1[58]+z[21]*tmp1[59]+z[23]*tmp1[60]
            +z[25]*tmp1[61]+z[27]*tmp1[62]+z[29]*tmp1[63]+z[31]*tmp1[64]+z[33]*tmp1[65]+z[35]*tmp1[66]+z[37]*tmp1[67]+z[39]*tmp1[68]+z[41]*tmp1[69]+z[43]*tmp1[70]
            +z[45]*tmp1[71]+z[47]*tmp1[72]+z[49]*tmp1[73]+z[51]*tmp1[74]+z[53]*tmp1[75]+z[55]*tmp1[76]+z[57]*tmp1[77]+z[59]*tmp1[78]+z[61]*tmp1[79]+z[63]*tmp1[80]
            +z[65]*tmp1[81]+z[67]*tmp1[82]+z[69]*tmp1[83]+z[71]*tmp1[84]+z[73]*tmp1[85]+z[75]*tmp1[86]+z[77]*tmp1[87]+z[79]*tmp1[88]+z[81]*tmp1[89]+z[83]*tmp1[90]
            +z[85]*tmp1[91]+z[87]*tmp1[92]+z[89]*tmp1[93]+z[91]*tmp1[94]+z[93]*tmp1[95]+z[95]*tmp1[96];
            tab[i+3*stg_first]=tab[i+3*stg_first]+z[0]*tmp1[0]
            +z[3]*tmp1[1]+z[6]*tmp1[2]+z[9]*tmp1[3]+z[12]*tmp1[4]+z[15]*tmp1[5]+z[18]*tmp1[6]+z[21]*tmp1[7]+z[24]*tmp1[8]+z[27]*tmp1[9]+z[30]*tmp1[10]
            +z[33]*tmp1[11]+z[36]*tmp1[12]+z[39]*tmp1[13]+z[42]*tmp1[14]+z[45]*tmp1[15]+z[48]*tmp1[16]+z[51]*tmp1[17]+z[54]*tmp1[18]+z[57]*tmp1[19]+z[60]*tmp1[20]
            +z[63]*tmp1[21]+z[66]*tmp1[22]+z[69]*tmp1[23]+z[72]*tmp1[24]+z[75]*tmp1[25]+z[78]*tmp1[26]+z[81]*tmp1[27]+z[84]*tmp1[28]+z[87]*tmp1[29]+z[90]*tmp1[30]
            +z[93]*tmp1[31]+z[96]*tmp1[32]+z[2]*tmp1[33]+z[5]*tmp1[34]+z[8]*tmp1[35]+z[11]*tmp1[36]+z[14]*tmp1[37]+z[17]*tmp1[38]+z[20]*tmp1[39]+z[23]*tmp1[40]
            +z[26]*tmp1[41]+z[29]*tmp1[42]+z[32]*tmp1[43]+z[35]*tmp1[44]+z[38]*tmp1[45]+z[41]*tmp1[46]+z[44]*tmp1[47]+z[47]*tmp1[48]+z[50]*tmp1[49]+z[53]*tmp1[50]
            +z[56]*tmp1[51]+z[59]*tmp1[52]+z[62]*tmp1[53]+z[65]*tmp1[54]+z[68]*tmp1[55]+z[71]*tmp1[56]+z[74]*tmp1[57]+z[77]*tmp1[58]+z[80]*tmp1[59]+z[83]*tmp1[60]
            +z[86]*tmp1[61]+z[89]*tmp1[62]+z[92]*tmp1[63]+z[95]*tmp1[64]+z[1]*tmp1[65]+z[4]*tmp1[66]+z[7]*tmp1[67]+z[10]*tmp1[68]+z[13]*tmp1[69]+z[16]*tmp1[70]
            +z[19]*tmp1[71]+z[22]*tmp1[72]+z[25]*tmp1[73]+z[28]*tmp1[74]+z[31]*tmp1[75]+z[34]*tmp1[76]+z[37]*tmp1[77]+z[40]*tmp1[78]+z[43]*tmp1[79]+z[46]*tmp1[80]
            +z[49]*tmp1[81]+z[52]*tmp1[82]+z[55]*tmp1[83]+z[58]*tmp1[84]+z[61]*tmp1[85]+z[64]*tmp1[86]+z[67]*tmp1[87]+z[70]*tmp1[88]+z[73]*tmp1[89]+z[76]*tmp1[90]
            +z[79]*tmp1[91]+z[82]*tmp1[92]+z[85]*tmp1[93]+z[88]*tmp1[94]+z[91]*tmp1[95]+z[94]*tmp1[96];
            tab[i+4*stg_first]=tab[i+4*stg_first]+z[0]*tmp1[0]
            +z[4]*tmp1[1]+z[8]*tmp1[2]+z[12]*tmp1[3]+z[16]*tmp1[4]+z[20]*tmp1[5]+z[24]*tmp1[6]+z[28]*tmp1[7]+z[32]*tmp1[8]+z[36]*tmp1[9]+z[40]*tmp1[10]
            +z[44]*tmp1[11]+z[48]*tmp1[12]+z[52]*tmp1[13]+z[56]*tmp1[14]+z[60]*tmp1[15]+z[64]*tmp1[16]+z[68]*tmp1[17]+z[72]*tmp1[18]+z[76]*tmp1[19]+z[80]*tmp1[20]
            +z[84]*tmp1[21]+z[88]*tmp1[22]+z[92]*tmp1[23]+z[96]*tmp1[24]+z[3]*tmp1[25]+z[7]*tmp1[26]+z[11]*tmp1[27]+z[15]*tmp1[28]+z[19]*tmp1[29]+z[23]*tmp1[30]
            +z[27]*tmp1[31]+z[31]*tmp1[32]+z[35]*tmp1[33]+z[39]*tmp1[34]+z[43]*tmp1[35]+z[47]*tmp1[36]+z[51]*tmp1[37]+z[55]*tmp1[38]+z[59]*tmp1[39]+z[63]*tmp1[40]
            +z[67]*tmp1[41]+z[71]*tmp1[42]+z[75]*tmp1[43]+z[79]*tmp1[44]+z[83]*tmp1[45]+z[87]*tmp1[46]+z[91]*tmp1[47]+z[95]*tmp1[48]+z[2]*tmp1[49]+z[6]*tmp1[50]
            +z[10]*tmp1[51]+z[14]*tmp1[52]+z[18]*tmp1[53]+z[22]*tmp1[54]+z[26]*tmp1[55]+z[30]*tmp1[56]+z[34]*tmp1[57]+z[38]*tmp1[58]+z[42]*tmp1[59]+z[46]*tmp1[60]
            +z[50]*tmp1[61]+z[54]*tmp1[62]+z[58]*tmp1[63]+z[62]*tmp1[64]+z[66]*tmp1[65]+z[70]*tmp1[66]+z[74]*tmp1[67]+z[78]*tmp1[68]+z[82]*tmp1[69]+z[86]*tmp1[70]
            +z[90]*tmp1[71]+z[94]*tmp1[72]+z[1]*tmp1[73]+z[5]*tmp1[74]+z[9]*tmp1[75]+z[13]*tmp1[76]+z[17]*tmp1[77]+z[21]*tmp1[78]+z[25]*tmp1[79]+z[29]*tmp1[80]
            +z[33]*tmp1[81]+z[37]*tmp1[82]+z[41]*tmp1[83]+z[45]*tmp1[84]+z[49]*tmp1[85]+z[53]*tmp1[86]+z[57]*tmp1[87]+z[61]*tmp1[88]+z[65]*tmp1[89]+z[69]*tmp1[90]
            +z[73]*tmp1[91]+z[77]*tmp1[92]+z[81]*tmp1[93]+z[85]*tmp1[94]+z[89]*tmp1[95]+z[93]*tmp1[96];
            tab[i+5*stg_first]=tab[i+5*stg_first]+z[0]*tmp1[0]
            +z[5]*tmp1[1]+z[10]*tmp1[2]+z[15]*tmp1[3]+z[20]*tmp1[4]+z[25]*tmp1[5]+z[30]*tmp1[6]+z[35]*tmp1[7]+z[40]*tmp1[8]+z[45]*tmp1[9]+z[50]*tmp1[10]
            +z[55]*tmp1[11]+z[60]*tmp1[12]+z[65]*tmp1[13]+z[70]*tmp1[14]+z[75]*tmp1[15]+z[80]*tmp1[16]+z[85]*tmp1[17]+z[90]*tmp1[18]+z[95]*tmp1[19]+z[3]*tmp1[20]
            +z[8]*tmp1[21]+z[13]*tmp1[22]+z[18]*tmp1[23]+z[23]*tmp1[24]+z[28]*tmp1[25]+z[33]*tmp1[26]+z[38]*tmp1[27]+z[43]*tmp1[28]+z[48]*tmp1[29]+z[53]*tmp1[30]
            +z[58]*tmp1[31]+z[63]*tmp1[32]+z[68]*tmp1[33]+z[73]*tmp1[34]+z[78]*tmp1[35]+z[83]*tmp1[36]+z[88]*tmp1[37]+z[93]*tmp1[38]+z[1]*tmp1[39]+z[6]*tmp1[40]
            +z[11]*tmp1[41]+z[16]*tmp1[42]+z[21]*tmp1[43]+z[26]*tmp1[44]+z[31]*tmp1[45]+z[36]*tmp1[46]+z[41]*tmp1[47]+z[46]*tmp1[48]+z[51]*tmp1[49]+z[56]*tmp1[50]
            +z[61]*tmp1[51]+z[66]*tmp1[52]+z[71]*tmp1[53]+z[76]*tmp1[54]+z[81]*tmp1[55]+z[86]*tmp1[56]+z[91]*tmp1[57]+z[96]*tmp1[58]+z[4]*tmp1[59]+z[9]*tmp1[60]
            +z[14]*tmp1[61]+z[19]*tmp1[62]+z[24]*tmp1[63]+z[29]*tmp1[64]+z[34]*tmp1[65]+z[39]*tmp1[66]+z[44]*tmp1[67]+z[49]*tmp1[68]+z[54]*tmp1[69]+z[59]*tmp1[70]
            +z[64]*tmp1[71]+z[69]*tmp1[72]+z[74]*tmp1[73]+z[79]*tmp1[74]+z[84]*tmp1[75]+z[89]*tmp1[76]+z[94]*tmp1[77]+z[2]*tmp1[78]+z[7]*tmp1[79]+z[12]*tmp1[80]
            +z[17]*tmp1[81]+z[22]*tmp1[82]+z[27]*tmp1[83]+z[32]*tmp1[84]+z[37]*tmp1[85]+z[42]*tmp1[86]+z[47]*tmp1[87]+z[52]*tmp1[88]+z[57]*tmp1[89]+z[62]*tmp1[90]
            +z[67]*tmp1[91]+z[72]*tmp1[92]+z[77]*tmp1[93]+z[82]*tmp1[94]+z[87]*tmp1[95]+z[92]*tmp1[96];
            tab[i+6*stg_first]=tab[i+6*stg_first]+z[0]*tmp1[0]
            +z[6]*tmp1[1]+z[12]*tmp1[2]+z[18]*tmp1[3]+z[24]*tmp1[4]+z[30]*tmp1[5]+z[36]*tmp1[6]+z[42]*tmp1[7]+z[48]*tmp1[8]+z[54]*tmp1[9]+z[60]*tmp1[10]
            +z[66]*tmp1[11]+z[72]*tmp1[12]+z[78]*tmp1[13]+z[84]*tmp1[14]+z[90]*tmp1[15]+z[96]*tmp1[16]+z[5]*tmp1[17]+z[11]*tmp1[18]+z[17]*tmp1[19]+z[23]*tmp1[20]
            +z[29]*tmp1[21]+z[35]*tmp1[22]+z[41]*tmp1[23]+z[47]*tmp1[24]+z[53]*tmp1[25]+z[59]*tmp1[26]+z[65]*tmp1[27]+z[71]*tmp1[28]+z[77]*tmp1[29]+z[83]*tmp1[30]
            +z[89]*tmp1[31]+z[95]*tmp1[32]+z[4]*tmp1[33]+z[10]*tmp1[34]+z[16]*tmp1[35]+z[22]*tmp1[36]+z[28]*tmp1[37]+z[34]*tmp1[38]+z[40]*tmp1[39]+z[46]*tmp1[40]
            +z[52]*tmp1[41]+z[58]*tmp1[42]+z[64]*tmp1[43]+z[70]*tmp1[44]+z[76]*tmp1[45]+z[82]*tmp1[46]+z[88]*tmp1[47]+z[94]*tmp1[48]+z[3]*tmp1[49]+z[9]*tmp1[50]
            +z[15]*tmp1[51]+z[21]*tmp1[52]+z[27]*tmp1[53]+z[33]*tmp1[54]+z[39]*tmp1[55]+z[45]*tmp1[56]+z[51]*tmp1[57]+z[57]*tmp1[58]+z[63]*tmp1[59]+z[69]*tmp1[60]
            +z[75]*tmp1[61]+z[81]*tmp1[62]+z[87]*tmp1[63]+z[93]*tmp1[64]+z[2]*tmp1[65]+z[8]*tmp1[66]+z[14]*tmp1[67]+z[20]*tmp1[68]+z[26]*tmp1[69]+z[32]*tmp1[70]
            +z[38]*tmp1[71]+z[44]*tmp1[72]+z[50]*tmp1[73]+z[56]*tmp1[74]+z[62]*tmp1[75]+z[68]*tmp1[76]+z[74]*tmp1[77]+z[80]*tmp1[78]+z[86]*tmp1[79]+z[92]*tmp1[80]
            +z[1]*tmp1[81]+z[7]*tmp1[82]+z[13]*tmp1[83]+z[19]*tmp1[84]+z[25]*tmp1[85]+z[31]*tmp1[86]+z[37]*tmp1[87]+z[43]*tmp1[88]+z[49]*tmp1[89]+z[55]*tmp1[90]
            +z[61]*tmp1[91]+z[67]*tmp1[92]+z[73]*tmp1[93]+z[79]*tmp1[94]+z[85]*tmp1[95]+z[91]*tmp1[96];
            tab[i+7*stg_first]=tab[i+7*stg_first]+z[0]*tmp1[0]
            +z[7]*tmp1[1]+z[14]*tmp1[2]+z[21]*tmp1[3]+z[28]*tmp1[4]+z[35]*tmp1[5]+z[42]*tmp1[6]+z[49]*tmp1[7]+z[56]*tmp1[8]+z[63]*tmp1[9]+z[70]*tmp1[10]
            +z[77]*tmp1[11]+z[84]*tmp1[12]+z[91]*tmp1[13]+z[1]*tmp1[14]+z[8]*tmp1[15]+z[15]*tmp1[16]+z[22]*tmp1[17]+z[29]*tmp1[18]+z[36]*tmp1[19]+z[43]*tmp1[20]
            +z[50]*tmp1[21]+z[57]*tmp1[22]+z[64]*tmp1[23]+z[71]*tmp1[24]+z[78]*tmp1[25]+z[85]*tmp1[26]+z[92]*tmp1[27]+z[2]*tmp1[28]+z[9]*tmp1[29]+z[16]*tmp1[30]
            +z[23]*tmp1[31]+z[30]*tmp1[32]+z[37]*tmp1[33]+z[44]*tmp1[34]+z[51]*tmp1[35]+z[58]*tmp1[36]+z[65]*tmp1[37]+z[72]*tmp1[38]+z[79]*tmp1[39]+z[86]*tmp1[40]
            +z[93]*tmp1[41]+z[3]*tmp1[42]+z[10]*tmp1[43]+z[17]*tmp1[44]+z[24]*tmp1[45]+z[31]*tmp1[46]+z[38]*tmp1[47]+z[45]*tmp1[48]+z[52]*tmp1[49]+z[59]*tmp1[50]
            +z[66]*tmp1[51]+z[73]*tmp1[52]+z[80]*tmp1[53]+z[87]*tmp1[54]+z[94]*tmp1[55]+z[4]*tmp1[56]+z[11]*tmp1[57]+z[18]*tmp1[58]+z[25]*tmp1[59]+z[32]*tmp1[60]
            +z[39]*tmp1[61]+z[46]*tmp1[62]+z[53]*tmp1[63]+z[60]*tmp1[64]+z[67]*tmp1[65]+z[74]*tmp1[66]+z[81]*tmp1[67]+z[88]*tmp1[68]+z[95]*tmp1[69]+z[5]*tmp1[70]
            +z[12]*tmp1[71]+z[19]*tmp1[72]+z[26]*tmp1[73]+z[33]*tmp1[74]+z[40]*tmp1[75]+z[47]*tmp1[76]+z[54]*tmp1[77]+z[61]*tmp1[78]+z[68]*tmp1[79]+z[75]*tmp1[80]
            +z[82]*tmp1[81]+z[89]*tmp1[82]+z[96]*tmp1[83]+z[6]*tmp1[84]+z[13]*tmp1[85]+z[20]*tmp1[86]+z[27]*tmp1[87]+z[34]*tmp1[88]+z[41]*tmp1[89]+z[48]*tmp1[90]
            +z[55]*tmp1[91]+z[62]*tmp1[92]+z[69]*tmp1[93]+z[76]*tmp1[94]+z[83]*tmp1[95]+z[90]*tmp1[96];
            tab[i+8*stg_first]=tab[i+8*stg_first]+z[0]*tmp1[0]
            +z[8]*tmp1[1]+z[16]*tmp1[2]+z[24]*tmp1[3]+z[32]*tmp1[4]+z[40]*tmp1[5]+z[48]*tmp1[6]+z[56]*tmp1[7]+z[64]*tmp1[8]+z[72]*tmp1[9]+z[80]*tmp1[10]
            +z[88]*tmp1[11]+z[96]*tmp1[12]+z[7]*tmp1[13]+z[15]*tmp1[14]+z[23]*tmp1[15]+z[31]*tmp1[16]+z[39]*tmp1[17]+z[47]*tmp1[18]+z[55]*tmp1[19]+z[63]*tmp1[20]
            +z[71]*tmp1[21]+z[79]*tmp1[22]+z[87]*tmp1[23]+z[95]*tmp1[24]+z[6]*tmp1[25]+z[14]*tmp1[26]+z[22]*tmp1[27]+z[30]*tmp1[28]+z[38]*tmp1[29]+z[46]*tmp1[30]
            +z[54]*tmp1[31]+z[62]*tmp1[32]+z[70]*tmp1[33]+z[78]*tmp1[34]+z[86]*tmp1[35]+z[94]*tmp1[36]+z[5]*tmp1[37]+z[13]*tmp1[38]+z[21]*tmp1[39]+z[29]*tmp1[40]
            +z[37]*tmp1[41]+z[45]*tmp1[42]+z[53]*tmp1[43]+z[61]*tmp1[44]+z[69]*tmp1[45]+z[77]*tmp1[46]+z[85]*tmp1[47]+z[93]*tmp1[48]+z[4]*tmp1[49]+z[12]*tmp1[50]
            +z[20]*tmp1[51]+z[28]*tmp1[52]+z[36]*tmp1[53]+z[44]*tmp1[54]+z[52]*tmp1[55]+z[60]*tmp1[56]+z[68]*tmp1[57]+z[76]*tmp1[58]+z[84]*tmp1[59]+z[92]*tmp1[60]
            +z[3]*tmp1[61]+z[11]*tmp1[62]+z[19]*tmp1[63]+z[27]*tmp1[64]+z[35]*tmp1[65]+z[43]*tmp1[66]+z[51]*tmp1[67]+z[59]*tmp1[68]+z[67]*tmp1[69]+z[75]*tmp1[70]
            +z[83]*tmp1[71]+z[91]*tmp1[72]+z[2]*tmp1[73]+z[10]*tmp1[74]+z[18]*tmp1[75]+z[26]*tmp1[76]+z[34]*tmp1[77]+z[42]*tmp1[78]+z[50]*tmp1[79]+z[58]*tmp1[80]
            +z[66]*tmp1[81]+z[74]*tmp1[82]+z[82]*tmp1[83]+z[90]*tmp1[84]+z[1]*tmp1[85]+z[9]*tmp1[86]+z[17]*tmp1[87]+z[25]*tmp1[88]+z[33]*tmp1[89]+z[41]*tmp1[90]
            +z[49]*tmp1[91]+z[57]*tmp1[92]+z[65]*tmp1[93]+z[73]*tmp1[94]+z[81]*tmp1[95]+z[89]*tmp1[96];
            tab[i+9*stg_first]=tab[i+9*stg_first]+z[0]*tmp1[0]
            +z[9]*tmp1[1]+z[18]*tmp1[2]+z[27]*tmp1[3]+z[36]*tmp1[4]+z[45]*tmp1[5]+z[54]*tmp1[6]+z[63]*tmp1[7]+z[72]*tmp1[8]+z[81]*tmp1[9]+z[90]*tmp1[10]
            +z[2]*tmp1[11]+z[11]*tmp1[12]+z[20]*tmp1[13]+z[29]*tmp1[14]+z[38]*tmp1[15]+z[47]*tmp1[16]+z[56]*tmp1[17]+z[65]*tmp1[18]+z[74]*tmp1[19]+z[83]*tmp1[20]
            +z[92]*tmp1[21]+z[4]*tmp1[22]+z[13]*tmp1[23]+z[22]*tmp1[24]+z[31]*tmp1[25]+z[40]*tmp1[26]+z[49]*tmp1[27]+z[58]*tmp1[28]+z[67]*tmp1[29]+z[76]*tmp1[30]
            +z[85]*tmp1[31]+z[94]*tmp1[32]+z[6]*tmp1[33]+z[15]*tmp1[34]+z[24]*tmp1[35]+z[33]*tmp1[36]+z[42]*tmp1[37]+z[51]*tmp1[38]+z[60]*tmp1[39]+z[69]*tmp1[40]
            +z[78]*tmp1[41]+z[87]*tmp1[42]+z[96]*tmp1[43]+z[8]*tmp1[44]+z[17]*tmp1[45]+z[26]*tmp1[46]+z[35]*tmp1[47]+z[44]*tmp1[48]+z[53]*tmp1[49]+z[62]*tmp1[50]
            +z[71]*tmp1[51]+z[80]*tmp1[52]+z[89]*tmp1[53]+z[1]*tmp1[54]+z[10]*tmp1[55]+z[19]*tmp1[56]+z[28]*tmp1[57]+z[37]*tmp1[58]+z[46]*tmp1[59]+z[55]*tmp1[60]
            +z[64]*tmp1[61]+z[73]*tmp1[62]+z[82]*tmp1[63]+z[91]*tmp1[64]+z[3]*tmp1[65]+z[12]*tmp1[66]+z[21]*tmp1[67]+z[30]*tmp1[68]+z[39]*tmp1[69]+z[48]*tmp1[70]
            +z[57]*tmp1[71]+z[66]*tmp1[72]+z[75]*tmp1[73]+z[84]*tmp1[74]+z[93]*tmp1[75]+z[5]*tmp1[76]+z[14]*tmp1[77]+z[23]*tmp1[78]+z[32]*tmp1[79]+z[41]*tmp1[80]
            +z[50]*tmp1[81]+z[59]*tmp1[82]+z[68]*tmp1[83]+z[77]*tmp1[84]+z[86]*tmp1[85]+z[95]*tmp1[86]+z[7]*tmp1[87]+z[16]*tmp1[88]+z[25]*tmp1[89]+z[34]*tmp1[90]
            +z[43]*tmp1[91]+z[52]*tmp1[92]+z[61]*tmp1[93]+z[70]*tmp1[94]+z[79]*tmp1[95]+z[88]*tmp1[96];
            tab[i+10*stg_first]=tab[i+10*stg_first]+z[0]*tmp1[0]
            +z[10]*tmp1[1]+z[20]*tmp1[2]+z[30]*tmp1[3]+z[40]*tmp1[4]+z[50]*tmp1[5]+z[60]*tmp1[6]+z[70]*tmp1[7]+z[80]*tmp1[8]+z[90]*tmp1[9]+z[3]*tmp1[10]
            +z[13]*tmp1[11]+z[23]*tmp1[12]+z[33]*tmp1[13]+z[43]*tmp1[14]+z[53]*tmp1[15]+z[63]*tmp1[16]+z[73]*tmp1[17]+z[83]*tmp1[18]+z[93]*tmp1[19]+z[6]*tmp1[20]
            +z[16]*tmp1[21]+z[26]*tmp1[22]+z[36]*tmp1[23]+z[46]*tmp1[24]+z[56]*tmp1[25]+z[66]*tmp1[26]+z[76]*tmp1[27]+z[86]*tmp1[28]+z[96]*tmp1[29]+z[9]*tmp1[30]
            +z[19]*tmp1[31]+z[29]*tmp1[32]+z[39]*tmp1[33]+z[49]*tmp1[34]+z[59]*tmp1[35]+z[69]*tmp1[36]+z[79]*tmp1[37]+z[89]*tmp1[38]+z[2]*tmp1[39]+z[12]*tmp1[40]
            +z[22]*tmp1[41]+z[32]*tmp1[42]+z[42]*tmp1[43]+z[52]*tmp1[44]+z[62]*tmp1[45]+z[72]*tmp1[46]+z[82]*tmp1[47]+z[92]*tmp1[48]+z[5]*tmp1[49]+z[15]*tmp1[50]
            +z[25]*tmp1[51]+z[35]*tmp1[52]+z[45]*tmp1[53]+z[55]*tmp1[54]+z[65]*tmp1[55]+z[75]*tmp1[56]+z[85]*tmp1[57]+z[95]*tmp1[58]+z[8]*tmp1[59]+z[18]*tmp1[60]
            +z[28]*tmp1[61]+z[38]*tmp1[62]+z[48]*tmp1[63]+z[58]*tmp1[64]+z[68]*tmp1[65]+z[78]*tmp1[66]+z[88]*tmp1[67]+z[1]*tmp1[68]+z[11]*tmp1[69]+z[21]*tmp1[70]
            +z[31]*tmp1[71]+z[41]*tmp1[72]+z[51]*tmp1[73]+z[61]*tmp1[74]+z[71]*tmp1[75]+z[81]*tmp1[76]+z[91]*tmp1[77]+z[4]*tmp1[78]+z[14]*tmp1[79]+z[24]*tmp1[80]
            +z[34]*tmp1[81]+z[44]*tmp1[82]+z[54]*tmp1[83]+z[64]*tmp1[84]+z[74]*tmp1[85]+z[84]*tmp1[86]+z[94]*tmp1[87]+z[7]*tmp1[88]+z[17]*tmp1[89]+z[27]*tmp1[90]
            +z[37]*tmp1[91]+z[47]*tmp1[92]+z[57]*tmp1[93]+z[67]*tmp1[94]+z[77]*tmp1[95]+z[87]*tmp1[96];
            tab[i+11*stg_first]=tab[i+11*stg_first]+z[0]*tmp1[0]
            +z[11]*tmp1[1]+z[22]*tmp1[2]+z[33]*tmp1[3]+z[44]*tmp1[4]+z[55]*tmp1[5]+z[66]*tmp1[6]+z[77]*tmp1[7]+z[88]*tmp1[8]+z[2]*tmp1[9]+z[13]*tmp1[10]
            +z[24]*tmp1[11]+z[35]*tmp1[12]+z[46]*tmp1[13]+z[57]*tmp1[14]+z[68]*tmp1[15]+z[79]*tmp1[16]+z[90]*tmp1[17]+z[4]*tmp1[18]+z[15]*tmp1[19]+z[26]*tmp1[20]
            +z[37]*tmp1[21]+z[48]*tmp1[22]+z[59]*tmp1[23]+z[70]*tmp1[24]+z[81]*tmp1[25]+z[92]*tmp1[26]+z[6]*tmp1[27]+z[17]*tmp1[28]+z[28]*tmp1[29]+z[39]*tmp1[30]
            +z[50]*tmp1[31]+z[61]*tmp1[32]+z[72]*tmp1[33]+z[83]*tmp1[34]+z[94]*tmp1[35]+z[8]*tmp1[36]+z[19]*tmp1[37]+z[30]*tmp1[38]+z[41]*tmp1[39]+z[52]*tmp1[40]
            +z[63]*tmp1[41]+z[74]*tmp1[42]+z[85]*tmp1[43]+z[96]*tmp1[44]+z[10]*tmp1[45]+z[21]*tmp1[46]+z[32]*tmp1[47]+z[43]*tmp1[48]+z[54]*tmp1[49]+z[65]*tmp1[50]
            +z[76]*tmp1[51]+z[87]*tmp1[52]+z[1]*tmp1[53]+z[12]*tmp1[54]+z[23]*tmp1[55]+z[34]*tmp1[56]+z[45]*tmp1[57]+z[56]*tmp1[58]+z[67]*tmp1[59]+z[78]*tmp1[60]
            +z[89]*tmp1[61]+z[3]*tmp1[62]+z[14]*tmp1[63]+z[25]*tmp1[64]+z[36]*tmp1[65]+z[47]*tmp1[66]+z[58]*tmp1[67]+z[69]*tmp1[68]+z[80]*tmp1[69]+z[91]*tmp1[70]
            +z[5]*tmp1[71]+z[16]*tmp1[72]+z[27]*tmp1[73]+z[38]*tmp1[74]+z[49]*tmp1[75]+z[60]*tmp1[76]+z[71]*tmp1[77]+z[82]*tmp1[78]+z[93]*tmp1[79]+z[7]*tmp1[80]
            +z[18]*tmp1[81]+z[29]*tmp1[82]+z[40]*tmp1[83]+z[51]*tmp1[84]+z[62]*tmp1[85]+z[73]*tmp1[86]+z[84]*tmp1[87]+z[95]*tmp1[88]+z[9]*tmp1[89]+z[20]*tmp1[90]
            +z[31]*tmp1[91]+z[42]*tmp1[92]+z[53]*tmp1[93]+z[64]*tmp1[94]+z[75]*tmp1[95]+z[86]*tmp1[96];
            tab[i+12*stg_first]=tab[i+12*stg_first]+z[0]*tmp1[0]
            +z[12]*tmp1[1]+z[24]*tmp1[2]+z[36]*tmp1[3]+z[48]*tmp1[4]+z[60]*tmp1[5]+z[72]*tmp1[6]+z[84]*tmp1[7]+z[96]*tmp1[8]+z[11]*tmp1[9]+z[23]*tmp1[10]
            +z[35]*tmp1[11]+z[47]*tmp1[12]+z[59]*tmp1[13]+z[71]*tmp1[14]+z[83]*tmp1[15]+z[95]*tmp1[16]+z[10]*tmp1[17]+z[22]*tmp1[18]+z[34]*tmp1[19]+z[46]*tmp1[20]
            +z[58]*tmp1[21]+z[70]*tmp1[22]+z[82]*tmp1[23]+z[94]*tmp1[24]+z[9]*tmp1[25]+z[21]*tmp1[26]+z[33]*tmp1[27]+z[45]*tmp1[28]+z[57]*tmp1[29]+z[69]*tmp1[30]
            +z[81]*tmp1[31]+z[93]*tmp1[32]+z[8]*tmp1[33]+z[20]*tmp1[34]+z[32]*tmp1[35]+z[44]*tmp1[36]+z[56]*tmp1[37]+z[68]*tmp1[38]+z[80]*tmp1[39]+z[92]*tmp1[40]
            +z[7]*tmp1[41]+z[19]*tmp1[42]+z[31]*tmp1[43]+z[43]*tmp1[44]+z[55]*tmp1[45]+z[67]*tmp1[46]+z[79]*tmp1[47]+z[91]*tmp1[48]+z[6]*tmp1[49]+z[18]*tmp1[50]
            +z[30]*tmp1[51]+z[42]*tmp1[52]+z[54]*tmp1[53]+z[66]*tmp1[54]+z[78]*tmp1[55]+z[90]*tmp1[56]+z[5]*tmp1[57]+z[17]*tmp1[58]+z[29]*tmp1[59]+z[41]*tmp1[60]
            +z[53]*tmp1[61]+z[65]*tmp1[62]+z[77]*tmp1[63]+z[89]*tmp1[64]+z[4]*tmp1[65]+z[16]*tmp1[66]+z[28]*tmp1[67]+z[40]*tmp1[68]+z[52]*tmp1[69]+z[64]*tmp1[70]
            +z[76]*tmp1[71]+z[88]*tmp1[72]+z[3]*tmp1[73]+z[15]*tmp1[74]+z[27]*tmp1[75]+z[39]*tmp1[76]+z[51]*tmp1[77]+z[63]*tmp1[78]+z[75]*tmp1[79]+z[87]*tmp1[80]
            +z[2]*tmp1[81]+z[14]*tmp1[82]+z[26]*tmp1[83]+z[38]*tmp1[84]+z[50]*tmp1[85]+z[62]*tmp1[86]+z[74]*tmp1[87]+z[86]*tmp1[88]+z[1]*tmp1[89]+z[13]*tmp1[90]
            +z[25]*tmp1[91]+z[37]*tmp1[92]+z[49]*tmp1[93]+z[61]*tmp1[94]+z[73]*tmp1[95]+z[85]*tmp1[96];
            tab[i+13*stg_first]=tab[i+13*stg_first]+z[0]*tmp1[0]
            +z[13]*tmp1[1]+z[26]*tmp1[2]+z[39]*tmp1[3]+z[52]*tmp1[4]+z[65]*tmp1[5]+z[78]*tmp1[6]+z[91]*tmp1[7]+z[7]*tmp1[8]+z[20]*tmp1[9]+z[33]*tmp1[10]
            +z[46]*tmp1[11]+z[59]*tmp1[12]+z[72]*tmp1[13]+z[85]*tmp1[14]+z[1]*tmp1[15]+z[14]*tmp1[16]+z[27]*tmp1[17]+z[40]*tmp1[18]+z[53]*tmp1[19]+z[66]*tmp1[20]
            +z[79]*tmp1[21]+z[92]*tmp1[22]+z[8]*tmp1[23]+z[21]*tmp1[24]+z[34]*tmp1[25]+z[47]*tmp1[26]+z[60]*tmp1[27]+z[73]*tmp1[28]+z[86]*tmp1[29]+z[2]*tmp1[30]
            +z[15]*tmp1[31]+z[28]*tmp1[32]+z[41]*tmp1[33]+z[54]*tmp1[34]+z[67]*tmp1[35]+z[80]*tmp1[36]+z[93]*tmp1[37]+z[9]*tmp1[38]+z[22]*tmp1[39]+z[35]*tmp1[40]
            +z[48]*tmp1[41]+z[61]*tmp1[42]+z[74]*tmp1[43]+z[87]*tmp1[44]+z[3]*tmp1[45]+z[16]*tmp1[46]+z[29]*tmp1[47]+z[42]*tmp1[48]+z[55]*tmp1[49]+z[68]*tmp1[50]
            +z[81]*tmp1[51]+z[94]*tmp1[52]+z[10]*tmp1[53]+z[23]*tmp1[54]+z[36]*tmp1[55]+z[49]*tmp1[56]+z[62]*tmp1[57]+z[75]*tmp1[58]+z[88]*tmp1[59]+z[4]*tmp1[60]
            +z[17]*tmp1[61]+z[30]*tmp1[62]+z[43]*tmp1[63]+z[56]*tmp1[64]+z[69]*tmp1[65]+z[82]*tmp1[66]+z[95]*tmp1[67]+z[11]*tmp1[68]+z[24]*tmp1[69]+z[37]*tmp1[70]
            +z[50]*tmp1[71]+z[63]*tmp1[72]+z[76]*tmp1[73]+z[89]*tmp1[74]+z[5]*tmp1[75]+z[18]*tmp1[76]+z[31]*tmp1[77]+z[44]*tmp1[78]+z[57]*tmp1[79]+z[70]*tmp1[80]
            +z[83]*tmp1[81]+z[96]*tmp1[82]+z[12]*tmp1[83]+z[25]*tmp1[84]+z[38]*tmp1[85]+z[51]*tmp1[86]+z[64]*tmp1[87]+z[77]*tmp1[88]+z[90]*tmp1[89]+z[6]*tmp1[90]
            +z[19]*tmp1[91]+z[32]*tmp1[92]+z[45]*tmp1[93]+z[58]*tmp1[94]+z[71]*tmp1[95]+z[84]*tmp1[96];
            tab[i+14*stg_first]=tab[i+14*stg_first]+z[0]*tmp1[0]
            +z[14]*tmp1[1]+z[28]*tmp1[2]+z[42]*tmp1[3]+z[56]*tmp1[4]+z[70]*tmp1[5]+z[84]*tmp1[6]+z[1]*tmp1[7]+z[15]*tmp1[8]+z[29]*tmp1[9]+z[43]*tmp1[10]
            +z[57]*tmp1[11]+z[71]*tmp1[12]+z[85]*tmp1[13]+z[2]*tmp1[14]+z[16]*tmp1[15]+z[30]*tmp1[16]+z[44]*tmp1[17]+z[58]*tmp1[18]+z[72]*tmp1[19]+z[86]*tmp1[20]
            +z[3]*tmp1[21]+z[17]*tmp1[22]+z[31]*tmp1[23]+z[45]*tmp1[24]+z[59]*tmp1[25]+z[73]*tmp1[26]+z[87]*tmp1[27]+z[4]*tmp1[28]+z[18]*tmp1[29]+z[32]*tmp1[30]
            +z[46]*tmp1[31]+z[60]*tmp1[32]+z[74]*tmp1[33]+z[88]*tmp1[34]+z[5]*tmp1[35]+z[19]*tmp1[36]+z[33]*tmp1[37]+z[47]*tmp1[38]+z[61]*tmp1[39]+z[75]*tmp1[40]
            +z[89]*tmp1[41]+z[6]*tmp1[42]+z[20]*tmp1[43]+z[34]*tmp1[44]+z[48]*tmp1[45]+z[62]*tmp1[46]+z[76]*tmp1[47]+z[90]*tmp1[48]+z[7]*tmp1[49]+z[21]*tmp1[50]
            +z[35]*tmp1[51]+z[49]*tmp1[52]+z[63]*tmp1[53]+z[77]*tmp1[54]+z[91]*tmp1[55]+z[8]*tmp1[56]+z[22]*tmp1[57]+z[36]*tmp1[58]+z[50]*tmp1[59]+z[64]*tmp1[60]
            +z[78]*tmp1[61]+z[92]*tmp1[62]+z[9]*tmp1[63]+z[23]*tmp1[64]+z[37]*tmp1[65]+z[51]*tmp1[66]+z[65]*tmp1[67]+z[79]*tmp1[68]+z[93]*tmp1[69]+z[10]*tmp1[70]
            +z[24]*tmp1[71]+z[38]*tmp1[72]+z[52]*tmp1[73]+z[66]*tmp1[74]+z[80]*tmp1[75]+z[94]*tmp1[76]+z[11]*tmp1[77]+z[25]*tmp1[78]+z[39]*tmp1[79]+z[53]*tmp1[80]
            +z[67]*tmp1[81]+z[81]*tmp1[82]+z[95]*tmp1[83]+z[12]*tmp1[84]+z[26]*tmp1[85]+z[40]*tmp1[86]+z[54]*tmp1[87]+z[68]*tmp1[88]+z[82]*tmp1[89]+z[96]*tmp1[90]
            +z[13]*tmp1[91]+z[27]*tmp1[92]+z[41]*tmp1[93]+z[55]*tmp1[94]+z[69]*tmp1[95]+z[83]*tmp1[96];
            tab[i+15*stg_first]=tab[i+15*stg_first]+z[0]*tmp1[0]
            +z[15]*tmp1[1]+z[30]*tmp1[2]+z[45]*tmp1[3]+z[60]*tmp1[4]+z[75]*tmp1[5]+z[90]*tmp1[6]+z[8]*tmp1[7]+z[23]*tmp1[8]+z[38]*tmp1[9]+z[53]*tmp1[10]
            +z[68]*tmp1[11]+z[83]*tmp1[12]+z[1]*tmp1[13]+z[16]*tmp1[14]+z[31]*tmp1[15]+z[46]*tmp1[16]+z[61]*tmp1[17]+z[76]*tmp1[18]+z[91]*tmp1[19]+z[9]*tmp1[20]
            +z[24]*tmp1[21]+z[39]*tmp1[22]+z[54]*tmp1[23]+z[69]*tmp1[24]+z[84]*tmp1[25]+z[2]*tmp1[26]+z[17]*tmp1[27]+z[32]*tmp1[28]+z[47]*tmp1[29]+z[62]*tmp1[30]
            +z[77]*tmp1[31]+z[92]*tmp1[32]+z[10]*tmp1[33]+z[25]*tmp1[34]+z[40]*tmp1[35]+z[55]*tmp1[36]+z[70]*tmp1[37]+z[85]*tmp1[38]+z[3]*tmp1[39]+z[18]*tmp1[40]
            +z[33]*tmp1[41]+z[48]*tmp1[42]+z[63]*tmp1[43]+z[78]*tmp1[44]+z[93]*tmp1[45]+z[11]*tmp1[46]+z[26]*tmp1[47]+z[41]*tmp1[48]+z[56]*tmp1[49]+z[71]*tmp1[50]
            +z[86]*tmp1[51]+z[4]*tmp1[52]+z[19]*tmp1[53]+z[34]*tmp1[54]+z[49]*tmp1[55]+z[64]*tmp1[56]+z[79]*tmp1[57]+z[94]*tmp1[58]+z[12]*tmp1[59]+z[27]*tmp1[60]
            +z[42]*tmp1[61]+z[57]*tmp1[62]+z[72]*tmp1[63]+z[87]*tmp1[64]+z[5]*tmp1[65]+z[20]*tmp1[66]+z[35]*tmp1[67]+z[50]*tmp1[68]+z[65]*tmp1[69]+z[80]*tmp1[70]
            +z[95]*tmp1[71]+z[13]*tmp1[72]+z[28]*tmp1[73]+z[43]*tmp1[74]+z[58]*tmp1[75]+z[73]*tmp1[76]+z[88]*tmp1[77]+z[6]*tmp1[78]+z[21]*tmp1[79]+z[36]*tmp1[80]
            +z[51]*tmp1[81]+z[66]*tmp1[82]+z[81]*tmp1[83]+z[96]*tmp1[84]+z[14]*tmp1[85]+z[29]*tmp1[86]+z[44]*tmp1[87]+z[59]*tmp1[88]+z[74]*tmp1[89]+z[89]*tmp1[90]
            +z[7]*tmp1[91]+z[22]*tmp1[92]+z[37]*tmp1[93]+z[52]*tmp1[94]+z[67]*tmp1[95]+z[82]*tmp1[96];
            tab[i+16*stg_first]=tab[i+16*stg_first]+z[0]*tmp1[0]
            +z[16]*tmp1[1]+z[32]*tmp1[2]+z[48]*tmp1[3]+z[64]*tmp1[4]+z[80]*tmp1[5]+z[96]*tmp1[6]+z[15]*tmp1[7]+z[31]*tmp1[8]+z[47]*tmp1[9]+z[63]*tmp1[10]
            +z[79]*tmp1[11]+z[95]*tmp1[12]+z[14]*tmp1[13]+z[30]*tmp1[14]+z[46]*tmp1[15]+z[62]*tmp1[16]+z[78]*tmp1[17]+z[94]*tmp1[18]+z[13]*tmp1[19]+z[29]*tmp1[20]
            +z[45]*tmp1[21]+z[61]*tmp1[22]+z[77]*tmp1[23]+z[93]*tmp1[24]+z[12]*tmp1[25]+z[28]*tmp1[26]+z[44]*tmp1[27]+z[60]*tmp1[28]+z[76]*tmp1[29]+z[92]*tmp1[30]
            +z[11]*tmp1[31]+z[27]*tmp1[32]+z[43]*tmp1[33]+z[59]*tmp1[34]+z[75]*tmp1[35]+z[91]*tmp1[36]+z[10]*tmp1[37]+z[26]*tmp1[38]+z[42]*tmp1[39]+z[58]*tmp1[40]
            +z[74]*tmp1[41]+z[90]*tmp1[42]+z[9]*tmp1[43]+z[25]*tmp1[44]+z[41]*tmp1[45]+z[57]*tmp1[46]+z[73]*tmp1[47]+z[89]*tmp1[48]+z[8]*tmp1[49]+z[24]*tmp1[50]
            +z[40]*tmp1[51]+z[56]*tmp1[52]+z[72]*tmp1[53]+z[88]*tmp1[54]+z[7]*tmp1[55]+z[23]*tmp1[56]+z[39]*tmp1[57]+z[55]*tmp1[58]+z[71]*tmp1[59]+z[87]*tmp1[60]
            +z[6]*tmp1[61]+z[22]*tmp1[62]+z[38]*tmp1[63]+z[54]*tmp1[64]+z[70]*tmp1[65]+z[86]*tmp1[66]+z[5]*tmp1[67]+z[21]*tmp1[68]+z[37]*tmp1[69]+z[53]*tmp1[70]
            +z[69]*tmp1[71]+z[85]*tmp1[72]+z[4]*tmp1[73]+z[20]*tmp1[74]+z[36]*tmp1[75]+z[52]*tmp1[76]+z[68]*tmp1[77]+z[84]*tmp1[78]+z[3]*tmp1[79]+z[19]*tmp1[80]
            +z[35]*tmp1[81]+z[51]*tmp1[82]+z[67]*tmp1[83]+z[83]*tmp1[84]+z[2]*tmp1[85]+z[18]*tmp1[86]+z[34]*tmp1[87]+z[50]*tmp1[88]+z[66]*tmp1[89]+z[82]*tmp1[90]
            +z[1]*tmp1[91]+z[17]*tmp1[92]+z[33]*tmp1[93]+z[49]*tmp1[94]+z[65]*tmp1[95]+z[81]*tmp1[96];
            tab[i+17*stg_first]=tab[i+17*stg_first]+z[0]*tmp1[0]
            +z[17]*tmp1[1]+z[34]*tmp1[2]+z[51]*tmp1[3]+z[68]*tmp1[4]+z[85]*tmp1[5]+z[5]*tmp1[6]+z[22]*tmp1[7]+z[39]*tmp1[8]+z[56]*tmp1[9]+z[73]*tmp1[10]
            +z[90]*tmp1[11]+z[10]*tmp1[12]+z[27]*tmp1[13]+z[44]*tmp1[14]+z[61]*tmp1[15]+z[78]*tmp1[16]+z[95]*tmp1[17]+z[15]*tmp1[18]+z[32]*tmp1[19]+z[49]*tmp1[20]
            +z[66]*tmp1[21]+z[83]*tmp1[22]+z[3]*tmp1[23]+z[20]*tmp1[24]+z[37]*tmp1[25]+z[54]*tmp1[26]+z[71]*tmp1[27]+z[88]*tmp1[28]+z[8]*tmp1[29]+z[25]*tmp1[30]
            +z[42]*tmp1[31]+z[59]*tmp1[32]+z[76]*tmp1[33]+z[93]*tmp1[34]+z[13]*tmp1[35]+z[30]*tmp1[36]+z[47]*tmp1[37]+z[64]*tmp1[38]+z[81]*tmp1[39]+z[1]*tmp1[40]
            +z[18]*tmp1[41]+z[35]*tmp1[42]+z[52]*tmp1[43]+z[69]*tmp1[44]+z[86]*tmp1[45]+z[6]*tmp1[46]+z[23]*tmp1[47]+z[40]*tmp1[48]+z[57]*tmp1[49]+z[74]*tmp1[50]
            +z[91]*tmp1[51]+z[11]*tmp1[52]+z[28]*tmp1[53]+z[45]*tmp1[54]+z[62]*tmp1[55]+z[79]*tmp1[56]+z[96]*tmp1[57]+z[16]*tmp1[58]+z[33]*tmp1[59]+z[50]*tmp1[60]
            +z[67]*tmp1[61]+z[84]*tmp1[62]+z[4]*tmp1[63]+z[21]*tmp1[64]+z[38]*tmp1[65]+z[55]*tmp1[66]+z[72]*tmp1[67]+z[89]*tmp1[68]+z[9]*tmp1[69]+z[26]*tmp1[70]
            +z[43]*tmp1[71]+z[60]*tmp1[72]+z[77]*tmp1[73]+z[94]*tmp1[74]+z[14]*tmp1[75]+z[31]*tmp1[76]+z[48]*tmp1[77]+z[65]*tmp1[78]+z[82]*tmp1[79]+z[2]*tmp1[80]
            +z[19]*tmp1[81]+z[36]*tmp1[82]+z[53]*tmp1[83]+z[70]*tmp1[84]+z[87]*tmp1[85]+z[7]*tmp1[86]+z[24]*tmp1[87]+z[41]*tmp1[88]+z[58]*tmp1[89]+z[75]*tmp1[90]
            +z[92]*tmp1[91]+z[12]*tmp1[92]+z[29]*tmp1[93]+z[46]*tmp1[94]+z[63]*tmp1[95]+z[80]*tmp1[96];
            tab[i+18*stg_first]=tab[i+18*stg_first]+z[0]*tmp1[0]
            +z[18]*tmp1[1]+z[36]*tmp1[2]+z[54]*tmp1[3]+z[72]*tmp1[4]+z[90]*tmp1[5]+z[11]*tmp1[6]+z[29]*tmp1[7]+z[47]*tmp1[8]+z[65]*tmp1[9]+z[83]*tmp1[10]
            +z[4]*tmp1[11]+z[22]*tmp1[12]+z[40]*tmp1[13]+z[58]*tmp1[14]+z[76]*tmp1[15]+z[94]*tmp1[16]+z[15]*tmp1[17]+z[33]*tmp1[18]+z[51]*tmp1[19]+z[69]*tmp1[20]
            +z[87]*tmp1[21]+z[8]*tmp1[22]+z[26]*tmp1[23]+z[44]*tmp1[24]+z[62]*tmp1[25]+z[80]*tmp1[26]+z[1]*tmp1[27]+z[19]*tmp1[28]+z[37]*tmp1[29]+z[55]*tmp1[30]
            +z[73]*tmp1[31]+z[91]*tmp1[32]+z[12]*tmp1[33]+z[30]*tmp1[34]+z[48]*tmp1[35]+z[66]*tmp1[36]+z[84]*tmp1[37]+z[5]*tmp1[38]+z[23]*tmp1[39]+z[41]*tmp1[40]
            +z[59]*tmp1[41]+z[77]*tmp1[42]+z[95]*tmp1[43]+z[16]*tmp1[44]+z[34]*tmp1[45]+z[52]*tmp1[46]+z[70]*tmp1[47]+z[88]*tmp1[48]+z[9]*tmp1[49]+z[27]*tmp1[50]
            +z[45]*tmp1[51]+z[63]*tmp1[52]+z[81]*tmp1[53]+z[2]*tmp1[54]+z[20]*tmp1[55]+z[38]*tmp1[56]+z[56]*tmp1[57]+z[74]*tmp1[58]+z[92]*tmp1[59]+z[13]*tmp1[60]
            +z[31]*tmp1[61]+z[49]*tmp1[62]+z[67]*tmp1[63]+z[85]*tmp1[64]+z[6]*tmp1[65]+z[24]*tmp1[66]+z[42]*tmp1[67]+z[60]*tmp1[68]+z[78]*tmp1[69]+z[96]*tmp1[70]
            +z[17]*tmp1[71]+z[35]*tmp1[72]+z[53]*tmp1[73]+z[71]*tmp1[74]+z[89]*tmp1[75]+z[10]*tmp1[76]+z[28]*tmp1[77]+z[46]*tmp1[78]+z[64]*tmp1[79]+z[82]*tmp1[80]
            +z[3]*tmp1[81]+z[21]*tmp1[82]+z[39]*tmp1[83]+z[57]*tmp1[84]+z[75]*tmp1[85]+z[93]*tmp1[86]+z[14]*tmp1[87]+z[32]*tmp1[88]+z[50]*tmp1[89]+z[68]*tmp1[90]
            +z[86]*tmp1[91]+z[7]*tmp1[92]+z[25]*tmp1[93]+z[43]*tmp1[94]+z[61]*tmp1[95]+z[79]*tmp1[96];
            tab[i+19*stg_first]=tab[i+19*stg_first]+z[0]*tmp1[0]
            +z[19]*tmp1[1]+z[38]*tmp1[2]+z[57]*tmp1[3]+z[76]*tmp1[4]+z[95]*tmp1[5]+z[17]*tmp1[6]+z[36]*tmp1[7]+z[55]*tmp1[8]+z[74]*tmp1[9]+z[93]*tmp1[10]
            +z[15]*tmp1[11]+z[34]*tmp1[12]+z[53]*tmp1[13]+z[72]*tmp1[14]+z[91]*tmp1[15]+z[13]*tmp1[16]+z[32]*tmp1[17]+z[51]*tmp1[18]+z[70]*tmp1[19]+z[89]*tmp1[20]
            +z[11]*tmp1[21]+z[30]*tmp1[22]+z[49]*tmp1[23]+z[68]*tmp1[24]+z[87]*tmp1[25]+z[9]*tmp1[26]+z[28]*tmp1[27]+z[47]*tmp1[28]+z[66]*tmp1[29]+z[85]*tmp1[30]
            +z[7]*tmp1[31]+z[26]*tmp1[32]+z[45]*tmp1[33]+z[64]*tmp1[34]+z[83]*tmp1[35]+z[5]*tmp1[36]+z[24]*tmp1[37]+z[43]*tmp1[38]+z[62]*tmp1[39]+z[81]*tmp1[40]
            +z[3]*tmp1[41]+z[22]*tmp1[42]+z[41]*tmp1[43]+z[60]*tmp1[44]+z[79]*tmp1[45]+z[1]*tmp1[46]+z[20]*tmp1[47]+z[39]*tmp1[48]+z[58]*tmp1[49]+z[77]*tmp1[50]
            +z[96]*tmp1[51]+z[18]*tmp1[52]+z[37]*tmp1[53]+z[56]*tmp1[54]+z[75]*tmp1[55]+z[94]*tmp1[56]+z[16]*tmp1[57]+z[35]*tmp1[58]+z[54]*tmp1[59]+z[73]*tmp1[60]
            +z[92]*tmp1[61]+z[14]*tmp1[62]+z[33]*tmp1[63]+z[52]*tmp1[64]+z[71]*tmp1[65]+z[90]*tmp1[66]+z[12]*tmp1[67]+z[31]*tmp1[68]+z[50]*tmp1[69]+z[69]*tmp1[70]
            +z[88]*tmp1[71]+z[10]*tmp1[72]+z[29]*tmp1[73]+z[48]*tmp1[74]+z[67]*tmp1[75]+z[86]*tmp1[76]+z[8]*tmp1[77]+z[27]*tmp1[78]+z[46]*tmp1[79]+z[65]*tmp1[80]
            +z[84]*tmp1[81]+z[6]*tmp1[82]+z[25]*tmp1[83]+z[44]*tmp1[84]+z[63]*tmp1[85]+z[82]*tmp1[86]+z[4]*tmp1[87]+z[23]*tmp1[88]+z[42]*tmp1[89]+z[61]*tmp1[90]
            +z[80]*tmp1[91]+z[2]*tmp1[92]+z[21]*tmp1[93]+z[40]*tmp1[94]+z[59]*tmp1[95]+z[78]*tmp1[96];
            tab[i+20*stg_first]=tab[i+20*stg_first]+z[0]*tmp1[0]
            +z[20]*tmp1[1]+z[40]*tmp1[2]+z[60]*tmp1[3]+z[80]*tmp1[4]+z[3]*tmp1[5]+z[23]*tmp1[6]+z[43]*tmp1[7]+z[63]*tmp1[8]+z[83]*tmp1[9]+z[6]*tmp1[10]
            +z[26]*tmp1[11]+z[46]*tmp1[12]+z[66]*tmp1[13]+z[86]*tmp1[14]+z[9]*tmp1[15]+z[29]*tmp1[16]+z[49]*tmp1[17]+z[69]*tmp1[18]+z[89]*tmp1[19]+z[12]*tmp1[20]
            +z[32]*tmp1[21]+z[52]*tmp1[22]+z[72]*tmp1[23]+z[92]*tmp1[24]+z[15]*tmp1[25]+z[35]*tmp1[26]+z[55]*tmp1[27]+z[75]*tmp1[28]+z[95]*tmp1[29]+z[18]*tmp1[30]
            +z[38]*tmp1[31]+z[58]*tmp1[32]+z[78]*tmp1[33]+z[1]*tmp1[34]+z[21]*tmp1[35]+z[41]*tmp1[36]+z[61]*tmp1[37]+z[81]*tmp1[38]+z[4]*tmp1[39]+z[24]*tmp1[40]
            +z[44]*tmp1[41]+z[64]*tmp1[42]+z[84]*tmp1[43]+z[7]*tmp1[44]+z[27]*tmp1[45]+z[47]*tmp1[46]+z[67]*tmp1[47]+z[87]*tmp1[48]+z[10]*tmp1[49]+z[30]*tmp1[50]
            +z[50]*tmp1[51]+z[70]*tmp1[52]+z[90]*tmp1[53]+z[13]*tmp1[54]+z[33]*tmp1[55]+z[53]*tmp1[56]+z[73]*tmp1[57]+z[93]*tmp1[58]+z[16]*tmp1[59]+z[36]*tmp1[60]
            +z[56]*tmp1[61]+z[76]*tmp1[62]+z[96]*tmp1[63]+z[19]*tmp1[64]+z[39]*tmp1[65]+z[59]*tmp1[66]+z[79]*tmp1[67]+z[2]*tmp1[68]+z[22]*tmp1[69]+z[42]*tmp1[70]
            +z[62]*tmp1[71]+z[82]*tmp1[72]+z[5]*tmp1[73]+z[25]*tmp1[74]+z[45]*tmp1[75]+z[65]*tmp1[76]+z[85]*tmp1[77]+z[8]*tmp1[78]+z[28]*tmp1[79]+z[48]*tmp1[80]
            +z[68]*tmp1[81]+z[88]*tmp1[82]+z[11]*tmp1[83]+z[31]*tmp1[84]+z[51]*tmp1[85]+z[71]*tmp1[86]+z[91]*tmp1[87]+z[14]*tmp1[88]+z[34]*tmp1[89]+z[54]*tmp1[90]
            +z[74]*tmp1[91]+z[94]*tmp1[92]+z[17]*tmp1[93]+z[37]*tmp1[94]+z[57]*tmp1[95]+z[77]*tmp1[96];
            tab[i+21*stg_first]=tab[i+21*stg_first]+z[0]*tmp1[0]
            +z[21]*tmp1[1]+z[42]*tmp1[2]+z[63]*tmp1[3]+z[84]*tmp1[4]+z[8]*tmp1[5]+z[29]*tmp1[6]+z[50]*tmp1[7]+z[71]*tmp1[8]+z[92]*tmp1[9]+z[16]*tmp1[10]
            +z[37]*tmp1[11]+z[58]*tmp1[12]+z[79]*tmp1[13]+z[3]*tmp1[14]+z[24]*tmp1[15]+z[45]*tmp1[16]+z[66]*tmp1[17]+z[87]*tmp1[18]+z[11]*tmp1[19]+z[32]*tmp1[20]
            +z[53]*tmp1[21]+z[74]*tmp1[22]+z[95]*tmp1[23]+z[19]*tmp1[24]+z[40]*tmp1[25]+z[61]*tmp1[26]+z[82]*tmp1[27]+z[6]*tmp1[28]+z[27]*tmp1[29]+z[48]*tmp1[30]
            +z[69]*tmp1[31]+z[90]*tmp1[32]+z[14]*tmp1[33]+z[35]*tmp1[34]+z[56]*tmp1[35]+z[77]*tmp1[36]+z[1]*tmp1[37]+z[22]*tmp1[38]+z[43]*tmp1[39]+z[64]*tmp1[40]
            +z[85]*tmp1[41]+z[9]*tmp1[42]+z[30]*tmp1[43]+z[51]*tmp1[44]+z[72]*tmp1[45]+z[93]*tmp1[46]+z[17]*tmp1[47]+z[38]*tmp1[48]+z[59]*tmp1[49]+z[80]*tmp1[50]
            +z[4]*tmp1[51]+z[25]*tmp1[52]+z[46]*tmp1[53]+z[67]*tmp1[54]+z[88]*tmp1[55]+z[12]*tmp1[56]+z[33]*tmp1[57]+z[54]*tmp1[58]+z[75]*tmp1[59]+z[96]*tmp1[60]
            +z[20]*tmp1[61]+z[41]*tmp1[62]+z[62]*tmp1[63]+z[83]*tmp1[64]+z[7]*tmp1[65]+z[28]*tmp1[66]+z[49]*tmp1[67]+z[70]*tmp1[68]+z[91]*tmp1[69]+z[15]*tmp1[70]
            +z[36]*tmp1[71]+z[57]*tmp1[72]+z[78]*tmp1[73]+z[2]*tmp1[74]+z[23]*tmp1[75]+z[44]*tmp1[76]+z[65]*tmp1[77]+z[86]*tmp1[78]+z[10]*tmp1[79]+z[31]*tmp1[80]
            +z[52]*tmp1[81]+z[73]*tmp1[82]+z[94]*tmp1[83]+z[18]*tmp1[84]+z[39]*tmp1[85]+z[60]*tmp1[86]+z[81]*tmp1[87]+z[5]*tmp1[88]+z[26]*tmp1[89]+z[47]*tmp1[90]
            +z[68]*tmp1[91]+z[89]*tmp1[92]+z[13]*tmp1[93]+z[34]*tmp1[94]+z[55]*tmp1[95]+z[76]*tmp1[96];
            tab[i+22*stg_first]=tab[i+22*stg_first]+z[0]*tmp1[0]
            +z[22]*tmp1[1]+z[44]*tmp1[2]+z[66]*tmp1[3]+z[88]*tmp1[4]+z[13]*tmp1[5]+z[35]*tmp1[6]+z[57]*tmp1[7]+z[79]*tmp1[8]+z[4]*tmp1[9]+z[26]*tmp1[10]
            +z[48]*tmp1[11]+z[70]*tmp1[12]+z[92]*tmp1[13]+z[17]*tmp1[14]+z[39]*tmp1[15]+z[61]*tmp1[16]+z[83]*tmp1[17]+z[8]*tmp1[18]+z[30]*tmp1[19]+z[52]*tmp1[20]
            +z[74]*tmp1[21]+z[96]*tmp1[22]+z[21]*tmp1[23]+z[43]*tmp1[24]+z[65]*tmp1[25]+z[87]*tmp1[26]+z[12]*tmp1[27]+z[34]*tmp1[28]+z[56]*tmp1[29]+z[78]*tmp1[30]
            +z[3]*tmp1[31]+z[25]*tmp1[32]+z[47]*tmp1[33]+z[69]*tmp1[34]+z[91]*tmp1[35]+z[16]*tmp1[36]+z[38]*tmp1[37]+z[60]*tmp1[38]+z[82]*tmp1[39]+z[7]*tmp1[40]
            +z[29]*tmp1[41]+z[51]*tmp1[42]+z[73]*tmp1[43]+z[95]*tmp1[44]+z[20]*tmp1[45]+z[42]*tmp1[46]+z[64]*tmp1[47]+z[86]*tmp1[48]+z[11]*tmp1[49]+z[33]*tmp1[50]
            +z[55]*tmp1[51]+z[77]*tmp1[52]+z[2]*tmp1[53]+z[24]*tmp1[54]+z[46]*tmp1[55]+z[68]*tmp1[56]+z[90]*tmp1[57]+z[15]*tmp1[58]+z[37]*tmp1[59]+z[59]*tmp1[60]
            +z[81]*tmp1[61]+z[6]*tmp1[62]+z[28]*tmp1[63]+z[50]*tmp1[64]+z[72]*tmp1[65]+z[94]*tmp1[66]+z[19]*tmp1[67]+z[41]*tmp1[68]+z[63]*tmp1[69]+z[85]*tmp1[70]
            +z[10]*tmp1[71]+z[32]*tmp1[72]+z[54]*tmp1[73]+z[76]*tmp1[74]+z[1]*tmp1[75]+z[23]*tmp1[76]+z[45]*tmp1[77]+z[67]*tmp1[78]+z[89]*tmp1[79]+z[14]*tmp1[80]
            +z[36]*tmp1[81]+z[58]*tmp1[82]+z[80]*tmp1[83]+z[5]*tmp1[84]+z[27]*tmp1[85]+z[49]*tmp1[86]+z[71]*tmp1[87]+z[93]*tmp1[88]+z[18]*tmp1[89]+z[40]*tmp1[90]
            +z[62]*tmp1[91]+z[84]*tmp1[92]+z[9]*tmp1[93]+z[31]*tmp1[94]+z[53]*tmp1[95]+z[75]*tmp1[96];
            tab[i+23*stg_first]=tab[i+23*stg_first]+z[0]*tmp1[0]
            +z[23]*tmp1[1]+z[46]*tmp1[2]+z[69]*tmp1[3]+z[92]*tmp1[4]+z[18]*tmp1[5]+z[41]*tmp1[6]+z[64]*tmp1[7]+z[87]*tmp1[8]+z[13]*tmp1[9]+z[36]*tmp1[10]
            +z[59]*tmp1[11]+z[82]*tmp1[12]+z[8]*tmp1[13]+z[31]*tmp1[14]+z[54]*tmp1[15]+z[77]*tmp1[16]+z[3]*tmp1[17]+z[26]*tmp1[18]+z[49]*tmp1[19]+z[72]*tmp1[20]
            +z[95]*tmp1[21]+z[21]*tmp1[22]+z[44]*tmp1[23]+z[67]*tmp1[24]+z[90]*tmp1[25]+z[16]*tmp1[26]+z[39]*tmp1[27]+z[62]*tmp1[28]+z[85]*tmp1[29]+z[11]*tmp1[30]
            +z[34]*tmp1[31]+z[57]*tmp1[32]+z[80]*tmp1[33]+z[6]*tmp1[34]+z[29]*tmp1[35]+z[52]*tmp1[36]+z[75]*tmp1[37]+z[1]*tmp1[38]+z[24]*tmp1[39]+z[47]*tmp1[40]
            +z[70]*tmp1[41]+z[93]*tmp1[42]+z[19]*tmp1[43]+z[42]*tmp1[44]+z[65]*tmp1[45]+z[88]*tmp1[46]+z[14]*tmp1[47]+z[37]*tmp1[48]+z[60]*tmp1[49]+z[83]*tmp1[50]
            +z[9]*tmp1[51]+z[32]*tmp1[52]+z[55]*tmp1[53]+z[78]*tmp1[54]+z[4]*tmp1[55]+z[27]*tmp1[56]+z[50]*tmp1[57]+z[73]*tmp1[58]+z[96]*tmp1[59]+z[22]*tmp1[60]
            +z[45]*tmp1[61]+z[68]*tmp1[62]+z[91]*tmp1[63]+z[17]*tmp1[64]+z[40]*tmp1[65]+z[63]*tmp1[66]+z[86]*tmp1[67]+z[12]*tmp1[68]+z[35]*tmp1[69]+z[58]*tmp1[70]
            +z[81]*tmp1[71]+z[7]*tmp1[72]+z[30]*tmp1[73]+z[53]*tmp1[74]+z[76]*tmp1[75]+z[2]*tmp1[76]+z[25]*tmp1[77]+z[48]*tmp1[78]+z[71]*tmp1[79]+z[94]*tmp1[80]
            +z[20]*tmp1[81]+z[43]*tmp1[82]+z[66]*tmp1[83]+z[89]*tmp1[84]+z[15]*tmp1[85]+z[38]*tmp1[86]+z[61]*tmp1[87]+z[84]*tmp1[88]+z[10]*tmp1[89]+z[33]*tmp1[90]
            +z[56]*tmp1[91]+z[79]*tmp1[92]+z[5]*tmp1[93]+z[28]*tmp1[94]+z[51]*tmp1[95]+z[74]*tmp1[96];
            tab[i+24*stg_first]=tab[i+24*stg_first]+z[0]*tmp1[0]
            +z[24]*tmp1[1]+z[48]*tmp1[2]+z[72]*tmp1[3]+z[96]*tmp1[4]+z[23]*tmp1[5]+z[47]*tmp1[6]+z[71]*tmp1[7]+z[95]*tmp1[8]+z[22]*tmp1[9]+z[46]*tmp1[10]
            +z[70]*tmp1[11]+z[94]*tmp1[12]+z[21]*tmp1[13]+z[45]*tmp1[14]+z[69]*tmp1[15]+z[93]*tmp1[16]+z[20]*tmp1[17]+z[44]*tmp1[18]+z[68]*tmp1[19]+z[92]*tmp1[20]
            +z[19]*tmp1[21]+z[43]*tmp1[22]+z[67]*tmp1[23]+z[91]*tmp1[24]+z[18]*tmp1[25]+z[42]*tmp1[26]+z[66]*tmp1[27]+z[90]*tmp1[28]+z[17]*tmp1[29]+z[41]*tmp1[30]
            +z[65]*tmp1[31]+z[89]*tmp1[32]+z[16]*tmp1[33]+z[40]*tmp1[34]+z[64]*tmp1[35]+z[88]*tmp1[36]+z[15]*tmp1[37]+z[39]*tmp1[38]+z[63]*tmp1[39]+z[87]*tmp1[40]
            +z[14]*tmp1[41]+z[38]*tmp1[42]+z[62]*tmp1[43]+z[86]*tmp1[44]+z[13]*tmp1[45]+z[37]*tmp1[46]+z[61]*tmp1[47]+z[85]*tmp1[48]+z[12]*tmp1[49]+z[36]*tmp1[50]
            +z[60]*tmp1[51]+z[84]*tmp1[52]+z[11]*tmp1[53]+z[35]*tmp1[54]+z[59]*tmp1[55]+z[83]*tmp1[56]+z[10]*tmp1[57]+z[34]*tmp1[58]+z[58]*tmp1[59]+z[82]*tmp1[60]
            +z[9]*tmp1[61]+z[33]*tmp1[62]+z[57]*tmp1[63]+z[81]*tmp1[64]+z[8]*tmp1[65]+z[32]*tmp1[66]+z[56]*tmp1[67]+z[80]*tmp1[68]+z[7]*tmp1[69]+z[31]*tmp1[70]
            +z[55]*tmp1[71]+z[79]*tmp1[72]+z[6]*tmp1[73]+z[30]*tmp1[74]+z[54]*tmp1[75]+z[78]*tmp1[76]+z[5]*tmp1[77]+z[29]*tmp1[78]+z[53]*tmp1[79]+z[77]*tmp1[80]
            +z[4]*tmp1[81]+z[28]*tmp1[82]+z[52]*tmp1[83]+z[76]*tmp1[84]+z[3]*tmp1[85]+z[27]*tmp1[86]+z[51]*tmp1[87]+z[75]*tmp1[88]+z[2]*tmp1[89]+z[26]*tmp1[90]
            +z[50]*tmp1[91]+z[74]*tmp1[92]+z[1]*tmp1[93]+z[25]*tmp1[94]+z[49]*tmp1[95]+z[73]*tmp1[96];
            tab[i+25*stg_first]=tab[i+25*stg_first]+z[0]*tmp1[0]
            +z[25]*tmp1[1]+z[50]*tmp1[2]+z[75]*tmp1[3]+z[3]*tmp1[4]+z[28]*tmp1[5]+z[53]*tmp1[6]+z[78]*tmp1[7]+z[6]*tmp1[8]+z[31]*tmp1[9]+z[56]*tmp1[10]
            +z[81]*tmp1[11]+z[9]*tmp1[12]+z[34]*tmp1[13]+z[59]*tmp1[14]+z[84]*tmp1[15]+z[12]*tmp1[16]+z[37]*tmp1[17]+z[62]*tmp1[18]+z[87]*tmp1[19]+z[15]*tmp1[20]
            +z[40]*tmp1[21]+z[65]*tmp1[22]+z[90]*tmp1[23]+z[18]*tmp1[24]+z[43]*tmp1[25]+z[68]*tmp1[26]+z[93]*tmp1[27]+z[21]*tmp1[28]+z[46]*tmp1[29]+z[71]*tmp1[30]
            +z[96]*tmp1[31]+z[24]*tmp1[32]+z[49]*tmp1[33]+z[74]*tmp1[34]+z[2]*tmp1[35]+z[27]*tmp1[36]+z[52]*tmp1[37]+z[77]*tmp1[38]+z[5]*tmp1[39]+z[30]*tmp1[40]
            +z[55]*tmp1[41]+z[80]*tmp1[42]+z[8]*tmp1[43]+z[33]*tmp1[44]+z[58]*tmp1[45]+z[83]*tmp1[46]+z[11]*tmp1[47]+z[36]*tmp1[48]+z[61]*tmp1[49]+z[86]*tmp1[50]
            +z[14]*tmp1[51]+z[39]*tmp1[52]+z[64]*tmp1[53]+z[89]*tmp1[54]+z[17]*tmp1[55]+z[42]*tmp1[56]+z[67]*tmp1[57]+z[92]*tmp1[58]+z[20]*tmp1[59]+z[45]*tmp1[60]
            +z[70]*tmp1[61]+z[95]*tmp1[62]+z[23]*tmp1[63]+z[48]*tmp1[64]+z[73]*tmp1[65]+z[1]*tmp1[66]+z[26]*tmp1[67]+z[51]*tmp1[68]+z[76]*tmp1[69]+z[4]*tmp1[70]
            +z[29]*tmp1[71]+z[54]*tmp1[72]+z[79]*tmp1[73]+z[7]*tmp1[74]+z[32]*tmp1[75]+z[57]*tmp1[76]+z[82]*tmp1[77]+z[10]*tmp1[78]+z[35]*tmp1[79]+z[60]*tmp1[80]
            +z[85]*tmp1[81]+z[13]*tmp1[82]+z[38]*tmp1[83]+z[63]*tmp1[84]+z[88]*tmp1[85]+z[16]*tmp1[86]+z[41]*tmp1[87]+z[66]*tmp1[88]+z[91]*tmp1[89]+z[19]*tmp1[90]
            +z[44]*tmp1[91]+z[69]*tmp1[92]+z[94]*tmp1[93]+z[22]*tmp1[94]+z[47]*tmp1[95]+z[72]*tmp1[96];
            tab[i+26*stg_first]=tab[i+26*stg_first]+z[0]*tmp1[0]
            +z[26]*tmp1[1]+z[52]*tmp1[2]+z[78]*tmp1[3]+z[7]*tmp1[4]+z[33]*tmp1[5]+z[59]*tmp1[6]+z[85]*tmp1[7]+z[14]*tmp1[8]+z[40]*tmp1[9]+z[66]*tmp1[10]
            +z[92]*tmp1[11]+z[21]*tmp1[12]+z[47]*tmp1[13]+z[73]*tmp1[14]+z[2]*tmp1[15]+z[28]*tmp1[16]+z[54]*tmp1[17]+z[80]*tmp1[18]+z[9]*tmp1[19]+z[35]*tmp1[20]
            +z[61]*tmp1[21]+z[87]*tmp1[22]+z[16]*tmp1[23]+z[42]*tmp1[24]+z[68]*tmp1[25]+z[94]*tmp1[26]+z[23]*tmp1[27]+z[49]*tmp1[28]+z[75]*tmp1[29]+z[4]*tmp1[30]
            +z[30]*tmp1[31]+z[56]*tmp1[32]+z[82]*tmp1[33]+z[11]*tmp1[34]+z[37]*tmp1[35]+z[63]*tmp1[36]+z[89]*tmp1[37]+z[18]*tmp1[38]+z[44]*tmp1[39]+z[70]*tmp1[40]
            +z[96]*tmp1[41]+z[25]*tmp1[42]+z[51]*tmp1[43]+z[77]*tmp1[44]+z[6]*tmp1[45]+z[32]*tmp1[46]+z[58]*tmp1[47]+z[84]*tmp1[48]+z[13]*tmp1[49]+z[39]*tmp1[50]
            +z[65]*tmp1[51]+z[91]*tmp1[52]+z[20]*tmp1[53]+z[46]*tmp1[54]+z[72]*tmp1[55]+z[1]*tmp1[56]+z[27]*tmp1[57]+z[53]*tmp1[58]+z[79]*tmp1[59]+z[8]*tmp1[60]
            +z[34]*tmp1[61]+z[60]*tmp1[62]+z[86]*tmp1[63]+z[15]*tmp1[64]+z[41]*tmp1[65]+z[67]*tmp1[66]+z[93]*tmp1[67]+z[22]*tmp1[68]+z[48]*tmp1[69]+z[74]*tmp1[70]
            +z[3]*tmp1[71]+z[29]*tmp1[72]+z[55]*tmp1[73]+z[81]*tmp1[74]+z[10]*tmp1[75]+z[36]*tmp1[76]+z[62]*tmp1[77]+z[88]*tmp1[78]+z[17]*tmp1[79]+z[43]*tmp1[80]
            +z[69]*tmp1[81]+z[95]*tmp1[82]+z[24]*tmp1[83]+z[50]*tmp1[84]+z[76]*tmp1[85]+z[5]*tmp1[86]+z[31]*tmp1[87]+z[57]*tmp1[88]+z[83]*tmp1[89]+z[12]*tmp1[90]
            +z[38]*tmp1[91]+z[64]*tmp1[92]+z[90]*tmp1[93]+z[19]*tmp1[94]+z[45]*tmp1[95]+z[71]*tmp1[96];
            tab[i+27*stg_first]=tab[i+27*stg_first]+z[0]*tmp1[0]
            +z[27]*tmp1[1]+z[54]*tmp1[2]+z[81]*tmp1[3]+z[11]*tmp1[4]+z[38]*tmp1[5]+z[65]*tmp1[6]+z[92]*tmp1[7]+z[22]*tmp1[8]+z[49]*tmp1[9]+z[76]*tmp1[10]
            +z[6]*tmp1[11]+z[33]*tmp1[12]+z[60]*tmp1[13]+z[87]*tmp1[14]+z[17]*tmp1[15]+z[44]*tmp1[16]+z[71]*tmp1[17]+z[1]*tmp1[18]+z[28]*tmp1[19]+z[55]*tmp1[20]
            +z[82]*tmp1[21]+z[12]*tmp1[22]+z[39]*tmp1[23]+z[66]*tmp1[24]+z[93]*tmp1[25]+z[23]*tmp1[26]+z[50]*tmp1[27]+z[77]*tmp1[28]+z[7]*tmp1[29]+z[34]*tmp1[30]
            +z[61]*tmp1[31]+z[88]*tmp1[32]+z[18]*tmp1[33]+z[45]*tmp1[34]+z[72]*tmp1[35]+z[2]*tmp1[36]+z[29]*tmp1[37]+z[56]*tmp1[38]+z[83]*tmp1[39]+z[13]*tmp1[40]
            +z[40]*tmp1[41]+z[67]*tmp1[42]+z[94]*tmp1[43]+z[24]*tmp1[44]+z[51]*tmp1[45]+z[78]*tmp1[46]+z[8]*tmp1[47]+z[35]*tmp1[48]+z[62]*tmp1[49]+z[89]*tmp1[50]
            +z[19]*tmp1[51]+z[46]*tmp1[52]+z[73]*tmp1[53]+z[3]*tmp1[54]+z[30]*tmp1[55]+z[57]*tmp1[56]+z[84]*tmp1[57]+z[14]*tmp1[58]+z[41]*tmp1[59]+z[68]*tmp1[60]
            +z[95]*tmp1[61]+z[25]*tmp1[62]+z[52]*tmp1[63]+z[79]*tmp1[64]+z[9]*tmp1[65]+z[36]*tmp1[66]+z[63]*tmp1[67]+z[90]*tmp1[68]+z[20]*tmp1[69]+z[47]*tmp1[70]
            +z[74]*tmp1[71]+z[4]*tmp1[72]+z[31]*tmp1[73]+z[58]*tmp1[74]+z[85]*tmp1[75]+z[15]*tmp1[76]+z[42]*tmp1[77]+z[69]*tmp1[78]+z[96]*tmp1[79]+z[26]*tmp1[80]
            +z[53]*tmp1[81]+z[80]*tmp1[82]+z[10]*tmp1[83]+z[37]*tmp1[84]+z[64]*tmp1[85]+z[91]*tmp1[86]+z[21]*tmp1[87]+z[48]*tmp1[88]+z[75]*tmp1[89]+z[5]*tmp1[90]
            +z[32]*tmp1[91]+z[59]*tmp1[92]+z[86]*tmp1[93]+z[16]*tmp1[94]+z[43]*tmp1[95]+z[70]*tmp1[96];
            tab[i+28*stg_first]=tab[i+28*stg_first]+z[0]*tmp1[0]
            +z[28]*tmp1[1]+z[56]*tmp1[2]+z[84]*tmp1[3]+z[15]*tmp1[4]+z[43]*tmp1[5]+z[71]*tmp1[6]+z[2]*tmp1[7]+z[30]*tmp1[8]+z[58]*tmp1[9]+z[86]*tmp1[10]
            +z[17]*tmp1[11]+z[45]*tmp1[12]+z[73]*tmp1[13]+z[4]*tmp1[14]+z[32]*tmp1[15]+z[60]*tmp1[16]+z[88]*tmp1[17]+z[19]*tmp1[18]+z[47]*tmp1[19]+z[75]*tmp1[20]
            +z[6]*tmp1[21]+z[34]*tmp1[22]+z[62]*tmp1[23]+z[90]*tmp1[24]+z[21]*tmp1[25]+z[49]*tmp1[26]+z[77]*tmp1[27]+z[8]*tmp1[28]+z[36]*tmp1[29]+z[64]*tmp1[30]
            +z[92]*tmp1[31]+z[23]*tmp1[32]+z[51]*tmp1[33]+z[79]*tmp1[34]+z[10]*tmp1[35]+z[38]*tmp1[36]+z[66]*tmp1[37]+z[94]*tmp1[38]+z[25]*tmp1[39]+z[53]*tmp1[40]
            +z[81]*tmp1[41]+z[12]*tmp1[42]+z[40]*tmp1[43]+z[68]*tmp1[44]+z[96]*tmp1[45]+z[27]*tmp1[46]+z[55]*tmp1[47]+z[83]*tmp1[48]+z[14]*tmp1[49]+z[42]*tmp1[50]
            +z[70]*tmp1[51]+z[1]*tmp1[52]+z[29]*tmp1[53]+z[57]*tmp1[54]+z[85]*tmp1[55]+z[16]*tmp1[56]+z[44]*tmp1[57]+z[72]*tmp1[58]+z[3]*tmp1[59]+z[31]*tmp1[60]
            +z[59]*tmp1[61]+z[87]*tmp1[62]+z[18]*tmp1[63]+z[46]*tmp1[64]+z[74]*tmp1[65]+z[5]*tmp1[66]+z[33]*tmp1[67]+z[61]*tmp1[68]+z[89]*tmp1[69]+z[20]*tmp1[70]
            +z[48]*tmp1[71]+z[76]*tmp1[72]+z[7]*tmp1[73]+z[35]*tmp1[74]+z[63]*tmp1[75]+z[91]*tmp1[76]+z[22]*tmp1[77]+z[50]*tmp1[78]+z[78]*tmp1[79]+z[9]*tmp1[80]
            +z[37]*tmp1[81]+z[65]*tmp1[82]+z[93]*tmp1[83]+z[24]*tmp1[84]+z[52]*tmp1[85]+z[80]*tmp1[86]+z[11]*tmp1[87]+z[39]*tmp1[88]+z[67]*tmp1[89]+z[95]*tmp1[90]
            +z[26]*tmp1[91]+z[54]*tmp1[92]+z[82]*tmp1[93]+z[13]*tmp1[94]+z[41]*tmp1[95]+z[69]*tmp1[96];
            tab[i+29*stg_first]=tab[i+29*stg_first]+z[0]*tmp1[0]
            +z[29]*tmp1[1]+z[58]*tmp1[2]+z[87]*tmp1[3]+z[19]*tmp1[4]+z[48]*tmp1[5]+z[77]*tmp1[6]+z[9]*tmp1[7]+z[38]*tmp1[8]+z[67]*tmp1[9]+z[96]*tmp1[10]
            +z[28]*tmp1[11]+z[57]*tmp1[12]+z[86]*tmp1[13]+z[18]*tmp1[14]+z[47]*tmp1[15]+z[76]*tmp1[16]+z[8]*tmp1[17]+z[37]*tmp1[18]+z[66]*tmp1[19]+z[95]*tmp1[20]
            +z[27]*tmp1[21]+z[56]*tmp1[22]+z[85]*tmp1[23]+z[17]*tmp1[24]+z[46]*tmp1[25]+z[75]*tmp1[26]+z[7]*tmp1[27]+z[36]*tmp1[28]+z[65]*tmp1[29]+z[94]*tmp1[30]
            +z[26]*tmp1[31]+z[55]*tmp1[32]+z[84]*tmp1[33]+z[16]*tmp1[34]+z[45]*tmp1[35]+z[74]*tmp1[36]+z[6]*tmp1[37]+z[35]*tmp1[38]+z[64]*tmp1[39]+z[93]*tmp1[40]
            +z[25]*tmp1[41]+z[54]*tmp1[42]+z[83]*tmp1[43]+z[15]*tmp1[44]+z[44]*tmp1[45]+z[73]*tmp1[46]+z[5]*tmp1[47]+z[34]*tmp1[48]+z[63]*tmp1[49]+z[92]*tmp1[50]
            +z[24]*tmp1[51]+z[53]*tmp1[52]+z[82]*tmp1[53]+z[14]*tmp1[54]+z[43]*tmp1[55]+z[72]*tmp1[56]+z[4]*tmp1[57]+z[33]*tmp1[58]+z[62]*tmp1[59]+z[91]*tmp1[60]
            +z[23]*tmp1[61]+z[52]*tmp1[62]+z[81]*tmp1[63]+z[13]*tmp1[64]+z[42]*tmp1[65]+z[71]*tmp1[66]+z[3]*tmp1[67]+z[32]*tmp1[68]+z[61]*tmp1[69]+z[90]*tmp1[70]
            +z[22]*tmp1[71]+z[51]*tmp1[72]+z[80]*tmp1[73]+z[12]*tmp1[74]+z[41]*tmp1[75]+z[70]*tmp1[76]+z[2]*tmp1[77]+z[31]*tmp1[78]+z[60]*tmp1[79]+z[89]*tmp1[80]
            +z[21]*tmp1[81]+z[50]*tmp1[82]+z[79]*tmp1[83]+z[11]*tmp1[84]+z[40]*tmp1[85]+z[69]*tmp1[86]+z[1]*tmp1[87]+z[30]*tmp1[88]+z[59]*tmp1[89]+z[88]*tmp1[90]
            +z[20]*tmp1[91]+z[49]*tmp1[92]+z[78]*tmp1[93]+z[10]*tmp1[94]+z[39]*tmp1[95]+z[68]*tmp1[96];
            tab[i+30*stg_first]=tab[i+30*stg_first]+z[0]*tmp1[0]
            +z[30]*tmp1[1]+z[60]*tmp1[2]+z[90]*tmp1[3]+z[23]*tmp1[4]+z[53]*tmp1[5]+z[83]*tmp1[6]+z[16]*tmp1[7]+z[46]*tmp1[8]+z[76]*tmp1[9]+z[9]*tmp1[10]
            +z[39]*tmp1[11]+z[69]*tmp1[12]+z[2]*tmp1[13]+z[32]*tmp1[14]+z[62]*tmp1[15]+z[92]*tmp1[16]+z[25]*tmp1[17]+z[55]*tmp1[18]+z[85]*tmp1[19]+z[18]*tmp1[20]
            +z[48]*tmp1[21]+z[78]*tmp1[22]+z[11]*tmp1[23]+z[41]*tmp1[24]+z[71]*tmp1[25]+z[4]*tmp1[26]+z[34]*tmp1[27]+z[64]*tmp1[28]+z[94]*tmp1[29]+z[27]*tmp1[30]
            +z[57]*tmp1[31]+z[87]*tmp1[32]+z[20]*tmp1[33]+z[50]*tmp1[34]+z[80]*tmp1[35]+z[13]*tmp1[36]+z[43]*tmp1[37]+z[73]*tmp1[38]+z[6]*tmp1[39]+z[36]*tmp1[40]
            +z[66]*tmp1[41]+z[96]*tmp1[42]+z[29]*tmp1[43]+z[59]*tmp1[44]+z[89]*tmp1[45]+z[22]*tmp1[46]+z[52]*tmp1[47]+z[82]*tmp1[48]+z[15]*tmp1[49]+z[45]*tmp1[50]
            +z[75]*tmp1[51]+z[8]*tmp1[52]+z[38]*tmp1[53]+z[68]*tmp1[54]+z[1]*tmp1[55]+z[31]*tmp1[56]+z[61]*tmp1[57]+z[91]*tmp1[58]+z[24]*tmp1[59]+z[54]*tmp1[60]
            +z[84]*tmp1[61]+z[17]*tmp1[62]+z[47]*tmp1[63]+z[77]*tmp1[64]+z[10]*tmp1[65]+z[40]*tmp1[66]+z[70]*tmp1[67]+z[3]*tmp1[68]+z[33]*tmp1[69]+z[63]*tmp1[70]
            +z[93]*tmp1[71]+z[26]*tmp1[72]+z[56]*tmp1[73]+z[86]*tmp1[74]+z[19]*tmp1[75]+z[49]*tmp1[76]+z[79]*tmp1[77]+z[12]*tmp1[78]+z[42]*tmp1[79]+z[72]*tmp1[80]
            +z[5]*tmp1[81]+z[35]*tmp1[82]+z[65]*tmp1[83]+z[95]*tmp1[84]+z[28]*tmp1[85]+z[58]*tmp1[86]+z[88]*tmp1[87]+z[21]*tmp1[88]+z[51]*tmp1[89]+z[81]*tmp1[90]
            +z[14]*tmp1[91]+z[44]*tmp1[92]+z[74]*tmp1[93]+z[7]*tmp1[94]+z[37]*tmp1[95]+z[67]*tmp1[96];
            tab[i+31*stg_first]=tab[i+31*stg_first]+z[0]*tmp1[0]
            +z[31]*tmp1[1]+z[62]*tmp1[2]+z[93]*tmp1[3]+z[27]*tmp1[4]+z[58]*tmp1[5]+z[89]*tmp1[6]+z[23]*tmp1[7]+z[54]*tmp1[8]+z[85]*tmp1[9]+z[19]*tmp1[10]
            +z[50]*tmp1[11]+z[81]*tmp1[12]+z[15]*tmp1[13]+z[46]*tmp1[14]+z[77]*tmp1[15]+z[11]*tmp1[16]+z[42]*tmp1[17]+z[73]*tmp1[18]+z[7]*tmp1[19]+z[38]*tmp1[20]
            +z[69]*tmp1[21]+z[3]*tmp1[22]+z[34]*tmp1[23]+z[65]*tmp1[24]+z[96]*tmp1[25]+z[30]*tmp1[26]+z[61]*tmp1[27]+z[92]*tmp1[28]+z[26]*tmp1[29]+z[57]*tmp1[30]
            +z[88]*tmp1[31]+z[22]*tmp1[32]+z[53]*tmp1[33]+z[84]*tmp1[34]+z[18]*tmp1[35]+z[49]*tmp1[36]+z[80]*tmp1[37]+z[14]*tmp1[38]+z[45]*tmp1[39]+z[76]*tmp1[40]
            +z[10]*tmp1[41]+z[41]*tmp1[42]+z[72]*tmp1[43]+z[6]*tmp1[44]+z[37]*tmp1[45]+z[68]*tmp1[46]+z[2]*tmp1[47]+z[33]*tmp1[48]+z[64]*tmp1[49]+z[95]*tmp1[50]
            +z[29]*tmp1[51]+z[60]*tmp1[52]+z[91]*tmp1[53]+z[25]*tmp1[54]+z[56]*tmp1[55]+z[87]*tmp1[56]+z[21]*tmp1[57]+z[52]*tmp1[58]+z[83]*tmp1[59]+z[17]*tmp1[60]
            +z[48]*tmp1[61]+z[79]*tmp1[62]+z[13]*tmp1[63]+z[44]*tmp1[64]+z[75]*tmp1[65]+z[9]*tmp1[66]+z[40]*tmp1[67]+z[71]*tmp1[68]+z[5]*tmp1[69]+z[36]*tmp1[70]
            +z[67]*tmp1[71]+z[1]*tmp1[72]+z[32]*tmp1[73]+z[63]*tmp1[74]+z[94]*tmp1[75]+z[28]*tmp1[76]+z[59]*tmp1[77]+z[90]*tmp1[78]+z[24]*tmp1[79]+z[55]*tmp1[80]
            +z[86]*tmp1[81]+z[20]*tmp1[82]+z[51]*tmp1[83]+z[82]*tmp1[84]+z[16]*tmp1[85]+z[47]*tmp1[86]+z[78]*tmp1[87]+z[12]*tmp1[88]+z[43]*tmp1[89]+z[74]*tmp1[90]
            +z[8]*tmp1[91]+z[39]*tmp1[92]+z[70]*tmp1[93]+z[4]*tmp1[94]+z[35]*tmp1[95]+z[66]*tmp1[96];
            tab[i+32*stg_first]=tab[i+32*stg_first]+z[0]*tmp1[0]
            +z[32]*tmp1[1]+z[64]*tmp1[2]+z[96]*tmp1[3]+z[31]*tmp1[4]+z[63]*tmp1[5]+z[95]*tmp1[6]+z[30]*tmp1[7]+z[62]*tmp1[8]+z[94]*tmp1[9]+z[29]*tmp1[10]
            +z[61]*tmp1[11]+z[93]*tmp1[12]+z[28]*tmp1[13]+z[60]*tmp1[14]+z[92]*tmp1[15]+z[27]*tmp1[16]+z[59]*tmp1[17]+z[91]*tmp1[18]+z[26]*tmp1[19]+z[58]*tmp1[20]
            +z[90]*tmp1[21]+z[25]*tmp1[22]+z[57]*tmp1[23]+z[89]*tmp1[24]+z[24]*tmp1[25]+z[56]*tmp1[26]+z[88]*tmp1[27]+z[23]*tmp1[28]+z[55]*tmp1[29]+z[87]*tmp1[30]
            +z[22]*tmp1[31]+z[54]*tmp1[32]+z[86]*tmp1[33]+z[21]*tmp1[34]+z[53]*tmp1[35]+z[85]*tmp1[36]+z[20]*tmp1[37]+z[52]*tmp1[38]+z[84]*tmp1[39]+z[19]*tmp1[40]
            +z[51]*tmp1[41]+z[83]*tmp1[42]+z[18]*tmp1[43]+z[50]*tmp1[44]+z[82]*tmp1[45]+z[17]*tmp1[46]+z[49]*tmp1[47]+z[81]*tmp1[48]+z[16]*tmp1[49]+z[48]*tmp1[50]
            +z[80]*tmp1[51]+z[15]*tmp1[52]+z[47]*tmp1[53]+z[79]*tmp1[54]+z[14]*tmp1[55]+z[46]*tmp1[56]+z[78]*tmp1[57]+z[13]*tmp1[58]+z[45]*tmp1[59]+z[77]*tmp1[60]
            +z[12]*tmp1[61]+z[44]*tmp1[62]+z[76]*tmp1[63]+z[11]*tmp1[64]+z[43]*tmp1[65]+z[75]*tmp1[66]+z[10]*tmp1[67]+z[42]*tmp1[68]+z[74]*tmp1[69]+z[9]*tmp1[70]
            +z[41]*tmp1[71]+z[73]*tmp1[72]+z[8]*tmp1[73]+z[40]*tmp1[74]+z[72]*tmp1[75]+z[7]*tmp1[76]+z[39]*tmp1[77]+z[71]*tmp1[78]+z[6]*tmp1[79]+z[38]*tmp1[80]
            +z[70]*tmp1[81]+z[5]*tmp1[82]+z[37]*tmp1[83]+z[69]*tmp1[84]+z[4]*tmp1[85]+z[36]*tmp1[86]+z[68]*tmp1[87]+z[3]*tmp1[88]+z[35]*tmp1[89]+z[67]*tmp1[90]
            +z[2]*tmp1[91]+z[34]*tmp1[92]+z[66]*tmp1[93]+z[1]*tmp1[94]+z[33]*tmp1[95]+z[65]*tmp1[96];
            tab[i+33*stg_first]=tab[i+33*stg_first]+z[0]*tmp1[0]
            +z[33]*tmp1[1]+z[66]*tmp1[2]+z[2]*tmp1[3]+z[35]*tmp1[4]+z[68]*tmp1[5]+z[4]*tmp1[6]+z[37]*tmp1[7]+z[70]*tmp1[8]+z[6]*tmp1[9]+z[39]*tmp1[10]
            +z[72]*tmp1[11]+z[8]*tmp1[12]+z[41]*tmp1[13]+z[74]*tmp1[14]+z[10]*tmp1[15]+z[43]*tmp1[16]+z[76]*tmp1[17]+z[12]*tmp1[18]+z[45]*tmp1[19]+z[78]*tmp1[20]
            +z[14]*tmp1[21]+z[47]*tmp1[22]+z[80]*tmp1[23]+z[16]*tmp1[24]+z[49]*tmp1[25]+z[82]*tmp1[26]+z[18]*tmp1[27]+z[51]*tmp1[28]+z[84]*tmp1[29]+z[20]*tmp1[30]
            +z[53]*tmp1[31]+z[86]*tmp1[32]+z[22]*tmp1[33]+z[55]*tmp1[34]+z[88]*tmp1[35]+z[24]*tmp1[36]+z[57]*tmp1[37]+z[90]*tmp1[38]+z[26]*tmp1[39]+z[59]*tmp1[40]
            +z[92]*tmp1[41]+z[28]*tmp1[42]+z[61]*tmp1[43]+z[94]*tmp1[44]+z[30]*tmp1[45]+z[63]*tmp1[46]+z[96]*tmp1[47]+z[32]*tmp1[48]+z[65]*tmp1[49]+z[1]*tmp1[50]
            +z[34]*tmp1[51]+z[67]*tmp1[52]+z[3]*tmp1[53]+z[36]*tmp1[54]+z[69]*tmp1[55]+z[5]*tmp1[56]+z[38]*tmp1[57]+z[71]*tmp1[58]+z[7]*tmp1[59]+z[40]*tmp1[60]
            +z[73]*tmp1[61]+z[9]*tmp1[62]+z[42]*tmp1[63]+z[75]*tmp1[64]+z[11]*tmp1[65]+z[44]*tmp1[66]+z[77]*tmp1[67]+z[13]*tmp1[68]+z[46]*tmp1[69]+z[79]*tmp1[70]
            +z[15]*tmp1[71]+z[48]*tmp1[72]+z[81]*tmp1[73]+z[17]*tmp1[74]+z[50]*tmp1[75]+z[83]*tmp1[76]+z[19]*tmp1[77]+z[52]*tmp1[78]+z[85]*tmp1[79]+z[21]*tmp1[80]
            +z[54]*tmp1[81]+z[87]*tmp1[82]+z[23]*tmp1[83]+z[56]*tmp1[84]+z[89]*tmp1[85]+z[25]*tmp1[86]+z[58]*tmp1[87]+z[91]*tmp1[88]+z[27]*tmp1[89]+z[60]*tmp1[90]
            +z[93]*tmp1[91]+z[29]*tmp1[92]+z[62]*tmp1[93]+z[95]*tmp1[94]+z[31]*tmp1[95]+z[64]*tmp1[96];
            tab[i+34*stg_first]=tab[i+34*stg_first]+z[0]*tmp1[0]
            +z[34]*tmp1[1]+z[68]*tmp1[2]+z[5]*tmp1[3]+z[39]*tmp1[4]+z[73]*tmp1[5]+z[10]*tmp1[6]+z[44]*tmp1[7]+z[78]*tmp1[8]+z[15]*tmp1[9]+z[49]*tmp1[10]
            +z[83]*tmp1[11]+z[20]*tmp1[12]+z[54]*tmp1[13]+z[88]*tmp1[14]+z[25]*tmp1[15]+z[59]*tmp1[16]+z[93]*tmp1[17]+z[30]*tmp1[18]+z[64]*tmp1[19]+z[1]*tmp1[20]
            +z[35]*tmp1[21]+z[69]*tmp1[22]+z[6]*tmp1[23]+z[40]*tmp1[24]+z[74]*tmp1[25]+z[11]*tmp1[26]+z[45]*tmp1[27]+z[79]*tmp1[28]+z[16]*tmp1[29]+z[50]*tmp1[30]
            +z[84]*tmp1[31]+z[21]*tmp1[32]+z[55]*tmp1[33]+z[89]*tmp1[34]+z[26]*tmp1[35]+z[60]*tmp1[36]+z[94]*tmp1[37]+z[31]*tmp1[38]+z[65]*tmp1[39]+z[2]*tmp1[40]
            +z[36]*tmp1[41]+z[70]*tmp1[42]+z[7]*tmp1[43]+z[41]*tmp1[44]+z[75]*tmp1[45]+z[12]*tmp1[46]+z[46]*tmp1[47]+z[80]*tmp1[48]+z[17]*tmp1[49]+z[51]*tmp1[50]
            +z[85]*tmp1[51]+z[22]*tmp1[52]+z[56]*tmp1[53]+z[90]*tmp1[54]+z[27]*tmp1[55]+z[61]*tmp1[56]+z[95]*tmp1[57]+z[32]*tmp1[58]+z[66]*tmp1[59]+z[3]*tmp1[60]
            +z[37]*tmp1[61]+z[71]*tmp1[62]+z[8]*tmp1[63]+z[42]*tmp1[64]+z[76]*tmp1[65]+z[13]*tmp1[66]+z[47]*tmp1[67]+z[81]*tmp1[68]+z[18]*tmp1[69]+z[52]*tmp1[70]
            +z[86]*tmp1[71]+z[23]*tmp1[72]+z[57]*tmp1[73]+z[91]*tmp1[74]+z[28]*tmp1[75]+z[62]*tmp1[76]+z[96]*tmp1[77]+z[33]*tmp1[78]+z[67]*tmp1[79]+z[4]*tmp1[80]
            +z[38]*tmp1[81]+z[72]*tmp1[82]+z[9]*tmp1[83]+z[43]*tmp1[84]+z[77]*tmp1[85]+z[14]*tmp1[86]+z[48]*tmp1[87]+z[82]*tmp1[88]+z[19]*tmp1[89]+z[53]*tmp1[90]
            +z[87]*tmp1[91]+z[24]*tmp1[92]+z[58]*tmp1[93]+z[92]*tmp1[94]+z[29]*tmp1[95]+z[63]*tmp1[96];
            tab[i+35*stg_first]=tab[i+35*stg_first]+z[0]*tmp1[0]
            +z[35]*tmp1[1]+z[70]*tmp1[2]+z[8]*tmp1[3]+z[43]*tmp1[4]+z[78]*tmp1[5]+z[16]*tmp1[6]+z[51]*tmp1[7]+z[86]*tmp1[8]+z[24]*tmp1[9]+z[59]*tmp1[10]
            +z[94]*tmp1[11]+z[32]*tmp1[12]+z[67]*tmp1[13]+z[5]*tmp1[14]+z[40]*tmp1[15]+z[75]*tmp1[16]+z[13]*tmp1[17]+z[48]*tmp1[18]+z[83]*tmp1[19]+z[21]*tmp1[20]
            +z[56]*tmp1[21]+z[91]*tmp1[22]+z[29]*tmp1[23]+z[64]*tmp1[24]+z[2]*tmp1[25]+z[37]*tmp1[26]+z[72]*tmp1[27]+z[10]*tmp1[28]+z[45]*tmp1[29]+z[80]*tmp1[30]
            +z[18]*tmp1[31]+z[53]*tmp1[32]+z[88]*tmp1[33]+z[26]*tmp1[34]+z[61]*tmp1[35]+z[96]*tmp1[36]+z[34]*tmp1[37]+z[69]*tmp1[38]+z[7]*tmp1[39]+z[42]*tmp1[40]
            +z[77]*tmp1[41]+z[15]*tmp1[42]+z[50]*tmp1[43]+z[85]*tmp1[44]+z[23]*tmp1[45]+z[58]*tmp1[46]+z[93]*tmp1[47]+z[31]*tmp1[48]+z[66]*tmp1[49]+z[4]*tmp1[50]
            +z[39]*tmp1[51]+z[74]*tmp1[52]+z[12]*tmp1[53]+z[47]*tmp1[54]+z[82]*tmp1[55]+z[20]*tmp1[56]+z[55]*tmp1[57]+z[90]*tmp1[58]+z[28]*tmp1[59]+z[63]*tmp1[60]
            +z[1]*tmp1[61]+z[36]*tmp1[62]+z[71]*tmp1[63]+z[9]*tmp1[64]+z[44]*tmp1[65]+z[79]*tmp1[66]+z[17]*tmp1[67]+z[52]*tmp1[68]+z[87]*tmp1[69]+z[25]*tmp1[70]
            +z[60]*tmp1[71]+z[95]*tmp1[72]+z[33]*tmp1[73]+z[68]*tmp1[74]+z[6]*tmp1[75]+z[41]*tmp1[76]+z[76]*tmp1[77]+z[14]*tmp1[78]+z[49]*tmp1[79]+z[84]*tmp1[80]
            +z[22]*tmp1[81]+z[57]*tmp1[82]+z[92]*tmp1[83]+z[30]*tmp1[84]+z[65]*tmp1[85]+z[3]*tmp1[86]+z[38]*tmp1[87]+z[73]*tmp1[88]+z[11]*tmp1[89]+z[46]*tmp1[90]
            +z[81]*tmp1[91]+z[19]*tmp1[92]+z[54]*tmp1[93]+z[89]*tmp1[94]+z[27]*tmp1[95]+z[62]*tmp1[96];
            tab[i+36*stg_first]=tab[i+36*stg_first]+z[0]*tmp1[0]
            +z[36]*tmp1[1]+z[72]*tmp1[2]+z[11]*tmp1[3]+z[47]*tmp1[4]+z[83]*tmp1[5]+z[22]*tmp1[6]+z[58]*tmp1[7]+z[94]*tmp1[8]+z[33]*tmp1[9]+z[69]*tmp1[10]
            +z[8]*tmp1[11]+z[44]*tmp1[12]+z[80]*tmp1[13]+z[19]*tmp1[14]+z[55]*tmp1[15]+z[91]*tmp1[16]+z[30]*tmp1[17]+z[66]*tmp1[18]+z[5]*tmp1[19]+z[41]*tmp1[20]
            +z[77]*tmp1[21]+z[16]*tmp1[22]+z[52]*tmp1[23]+z[88]*tmp1[24]+z[27]*tmp1[25]+z[63]*tmp1[26]+z[2]*tmp1[27]+z[38]*tmp1[28]+z[74]*tmp1[29]+z[13]*tmp1[30]
            +z[49]*tmp1[31]+z[85]*tmp1[32]+z[24]*tmp1[33]+z[60]*tmp1[34]+z[96]*tmp1[35]+z[35]*tmp1[36]+z[71]*tmp1[37]+z[10]*tmp1[38]+z[46]*tmp1[39]+z[82]*tmp1[40]
            +z[21]*tmp1[41]+z[57]*tmp1[42]+z[93]*tmp1[43]+z[32]*tmp1[44]+z[68]*tmp1[45]+z[7]*tmp1[46]+z[43]*tmp1[47]+z[79]*tmp1[48]+z[18]*tmp1[49]+z[54]*tmp1[50]
            +z[90]*tmp1[51]+z[29]*tmp1[52]+z[65]*tmp1[53]+z[4]*tmp1[54]+z[40]*tmp1[55]+z[76]*tmp1[56]+z[15]*tmp1[57]+z[51]*tmp1[58]+z[87]*tmp1[59]+z[26]*tmp1[60]
            +z[62]*tmp1[61]+z[1]*tmp1[62]+z[37]*tmp1[63]+z[73]*tmp1[64]+z[12]*tmp1[65]+z[48]*tmp1[66]+z[84]*tmp1[67]+z[23]*tmp1[68]+z[59]*tmp1[69]+z[95]*tmp1[70]
            +z[34]*tmp1[71]+z[70]*tmp1[72]+z[9]*tmp1[73]+z[45]*tmp1[74]+z[81]*tmp1[75]+z[20]*tmp1[76]+z[56]*tmp1[77]+z[92]*tmp1[78]+z[31]*tmp1[79]+z[67]*tmp1[80]
            +z[6]*tmp1[81]+z[42]*tmp1[82]+z[78]*tmp1[83]+z[17]*tmp1[84]+z[53]*tmp1[85]+z[89]*tmp1[86]+z[28]*tmp1[87]+z[64]*tmp1[88]+z[3]*tmp1[89]+z[39]*tmp1[90]
            +z[75]*tmp1[91]+z[14]*tmp1[92]+z[50]*tmp1[93]+z[86]*tmp1[94]+z[25]*tmp1[95]+z[61]*tmp1[96];
            tab[i+37*stg_first]=tab[i+37*stg_first]+z[0]*tmp1[0]
            +z[37]*tmp1[1]+z[74]*tmp1[2]+z[14]*tmp1[3]+z[51]*tmp1[4]+z[88]*tmp1[5]+z[28]*tmp1[6]+z[65]*tmp1[7]+z[5]*tmp1[8]+z[42]*tmp1[9]+z[79]*tmp1[10]
            +z[19]*tmp1[11]+z[56]*tmp1[12]+z[93]*tmp1[13]+z[33]*tmp1[14]+z[70]*tmp1[15]+z[10]*tmp1[16]+z[47]*tmp1[17]+z[84]*tmp1[18]+z[24]*tmp1[19]+z[61]*tmp1[20]
            +z[1]*tmp1[21]+z[38]*tmp1[22]+z[75]*tmp1[23]+z[15]*tmp1[24]+z[52]*tmp1[25]+z[89]*tmp1[26]+z[29]*tmp1[27]+z[66]*tmp1[28]+z[6]*tmp1[29]+z[43]*tmp1[30]
            +z[80]*tmp1[31]+z[20]*tmp1[32]+z[57]*tmp1[33]+z[94]*tmp1[34]+z[34]*tmp1[35]+z[71]*tmp1[36]+z[11]*tmp1[37]+z[48]*tmp1[38]+z[85]*tmp1[39]+z[25]*tmp1[40]
            +z[62]*tmp1[41]+z[2]*tmp1[42]+z[39]*tmp1[43]+z[76]*tmp1[44]+z[16]*tmp1[45]+z[53]*tmp1[46]+z[90]*tmp1[47]+z[30]*tmp1[48]+z[67]*tmp1[49]+z[7]*tmp1[50]
            +z[44]*tmp1[51]+z[81]*tmp1[52]+z[21]*tmp1[53]+z[58]*tmp1[54]+z[95]*tmp1[55]+z[35]*tmp1[56]+z[72]*tmp1[57]+z[12]*tmp1[58]+z[49]*tmp1[59]+z[86]*tmp1[60]
            +z[26]*tmp1[61]+z[63]*tmp1[62]+z[3]*tmp1[63]+z[40]*tmp1[64]+z[77]*tmp1[65]+z[17]*tmp1[66]+z[54]*tmp1[67]+z[91]*tmp1[68]+z[31]*tmp1[69]+z[68]*tmp1[70]
            +z[8]*tmp1[71]+z[45]*tmp1[72]+z[82]*tmp1[73]+z[22]*tmp1[74]+z[59]*tmp1[75]+z[96]*tmp1[76]+z[36]*tmp1[77]+z[73]*tmp1[78]+z[13]*tmp1[79]+z[50]*tmp1[80]
            +z[87]*tmp1[81]+z[27]*tmp1[82]+z[64]*tmp1[83]+z[4]*tmp1[84]+z[41]*tmp1[85]+z[78]*tmp1[86]+z[18]*tmp1[87]+z[55]*tmp1[88]+z[92]*tmp1[89]+z[32]*tmp1[90]
            +z[69]*tmp1[91]+z[9]*tmp1[92]+z[46]*tmp1[93]+z[83]*tmp1[94]+z[23]*tmp1[95]+z[60]*tmp1[96];
            tab[i+38*stg_first]=tab[i+38*stg_first]+z[0]*tmp1[0]
            +z[38]*tmp1[1]+z[76]*tmp1[2]+z[17]*tmp1[3]+z[55]*tmp1[4]+z[93]*tmp1[5]+z[34]*tmp1[6]+z[72]*tmp1[7]+z[13]*tmp1[8]+z[51]*tmp1[9]+z[89]*tmp1[10]
            +z[30]*tmp1[11]+z[68]*tmp1[12]+z[9]*tmp1[13]+z[47]*tmp1[14]+z[85]*tmp1[15]+z[26]*tmp1[16]+z[64]*tmp1[17]+z[5]*tmp1[18]+z[43]*tmp1[19]+z[81]*tmp1[20]
            +z[22]*tmp1[21]+z[60]*tmp1[22]+z[1]*tmp1[23]+z[39]*tmp1[24]+z[77]*tmp1[25]+z[18]*tmp1[26]+z[56]*tmp1[27]+z[94]*tmp1[28]+z[35]*tmp1[29]+z[73]*tmp1[30]
            +z[14]*tmp1[31]+z[52]*tmp1[32]+z[90]*tmp1[33]+z[31]*tmp1[34]+z[69]*tmp1[35]+z[10]*tmp1[36]+z[48]*tmp1[37]+z[86]*tmp1[38]+z[27]*tmp1[39]+z[65]*tmp1[40]
            +z[6]*tmp1[41]+z[44]*tmp1[42]+z[82]*tmp1[43]+z[23]*tmp1[44]+z[61]*tmp1[45]+z[2]*tmp1[46]+z[40]*tmp1[47]+z[78]*tmp1[48]+z[19]*tmp1[49]+z[57]*tmp1[50]
            +z[95]*tmp1[51]+z[36]*tmp1[52]+z[74]*tmp1[53]+z[15]*tmp1[54]+z[53]*tmp1[55]+z[91]*tmp1[56]+z[32]*tmp1[57]+z[70]*tmp1[58]+z[11]*tmp1[59]+z[49]*tmp1[60]
            +z[87]*tmp1[61]+z[28]*tmp1[62]+z[66]*tmp1[63]+z[7]*tmp1[64]+z[45]*tmp1[65]+z[83]*tmp1[66]+z[24]*tmp1[67]+z[62]*tmp1[68]+z[3]*tmp1[69]+z[41]*tmp1[70]
            +z[79]*tmp1[71]+z[20]*tmp1[72]+z[58]*tmp1[73]+z[96]*tmp1[74]+z[37]*tmp1[75]+z[75]*tmp1[76]+z[16]*tmp1[77]+z[54]*tmp1[78]+z[92]*tmp1[79]+z[33]*tmp1[80]
            +z[71]*tmp1[81]+z[12]*tmp1[82]+z[50]*tmp1[83]+z[88]*tmp1[84]+z[29]*tmp1[85]+z[67]*tmp1[86]+z[8]*tmp1[87]+z[46]*tmp1[88]+z[84]*tmp1[89]+z[25]*tmp1[90]
            +z[63]*tmp1[91]+z[4]*tmp1[92]+z[42]*tmp1[93]+z[80]*tmp1[94]+z[21]*tmp1[95]+z[59]*tmp1[96];
            tab[i+39*stg_first]=tab[i+39*stg_first]+z[0]*tmp1[0]
            +z[39]*tmp1[1]+z[78]*tmp1[2]+z[20]*tmp1[3]+z[59]*tmp1[4]+z[1]*tmp1[5]+z[40]*tmp1[6]+z[79]*tmp1[7]+z[21]*tmp1[8]+z[60]*tmp1[9]+z[2]*tmp1[10]
            +z[41]*tmp1[11]+z[80]*tmp1[12]+z[22]*tmp1[13]+z[61]*tmp1[14]+z[3]*tmp1[15]+z[42]*tmp1[16]+z[81]*tmp1[17]+z[23]*tmp1[18]+z[62]*tmp1[19]+z[4]*tmp1[20]
            +z[43]*tmp1[21]+z[82]*tmp1[22]+z[24]*tmp1[23]+z[63]*tmp1[24]+z[5]*tmp1[25]+z[44]*tmp1[26]+z[83]*tmp1[27]+z[25]*tmp1[28]+z[64]*tmp1[29]+z[6]*tmp1[30]
            +z[45]*tmp1[31]+z[84]*tmp1[32]+z[26]*tmp1[33]+z[65]*tmp1[34]+z[7]*tmp1[35]+z[46]*tmp1[36]+z[85]*tmp1[37]+z[27]*tmp1[38]+z[66]*tmp1[39]+z[8]*tmp1[40]
            +z[47]*tmp1[41]+z[86]*tmp1[42]+z[28]*tmp1[43]+z[67]*tmp1[44]+z[9]*tmp1[45]+z[48]*tmp1[46]+z[87]*tmp1[47]+z[29]*tmp1[48]+z[68]*tmp1[49]+z[10]*tmp1[50]
            +z[49]*tmp1[51]+z[88]*tmp1[52]+z[30]*tmp1[53]+z[69]*tmp1[54]+z[11]*tmp1[55]+z[50]*tmp1[56]+z[89]*tmp1[57]+z[31]*tmp1[58]+z[70]*tmp1[59]+z[12]*tmp1[60]
            +z[51]*tmp1[61]+z[90]*tmp1[62]+z[32]*tmp1[63]+z[71]*tmp1[64]+z[13]*tmp1[65]+z[52]*tmp1[66]+z[91]*tmp1[67]+z[33]*tmp1[68]+z[72]*tmp1[69]+z[14]*tmp1[70]
            +z[53]*tmp1[71]+z[92]*tmp1[72]+z[34]*tmp1[73]+z[73]*tmp1[74]+z[15]*tmp1[75]+z[54]*tmp1[76]+z[93]*tmp1[77]+z[35]*tmp1[78]+z[74]*tmp1[79]+z[16]*tmp1[80]
            +z[55]*tmp1[81]+z[94]*tmp1[82]+z[36]*tmp1[83]+z[75]*tmp1[84]+z[17]*tmp1[85]+z[56]*tmp1[86]+z[95]*tmp1[87]+z[37]*tmp1[88]+z[76]*tmp1[89]+z[18]*tmp1[90]
            +z[57]*tmp1[91]+z[96]*tmp1[92]+z[38]*tmp1[93]+z[77]*tmp1[94]+z[19]*tmp1[95]+z[58]*tmp1[96];
            tab[i+40*stg_first]=tab[i+40*stg_first]+z[0]*tmp1[0]
            +z[40]*tmp1[1]+z[80]*tmp1[2]+z[23]*tmp1[3]+z[63]*tmp1[4]+z[6]*tmp1[5]+z[46]*tmp1[6]+z[86]*tmp1[7]+z[29]*tmp1[8]+z[69]*tmp1[9]+z[12]*tmp1[10]
            +z[52]*tmp1[11]+z[92]*tmp1[12]+z[35]*tmp1[13]+z[75]*tmp1[14]+z[18]*tmp1[15]+z[58]*tmp1[16]+z[1]*tmp1[17]+z[41]*tmp1[18]+z[81]*tmp1[19]+z[24]*tmp1[20]
            +z[64]*tmp1[21]+z[7]*tmp1[22]+z[47]*tmp1[23]+z[87]*tmp1[24]+z[30]*tmp1[25]+z[70]*tmp1[26]+z[13]*tmp1[27]+z[53]*tmp1[28]+z[93]*tmp1[29]+z[36]*tmp1[30]
            +z[76]*tmp1[31]+z[19]*tmp1[32]+z[59]*tmp1[33]+z[2]*tmp1[34]+z[42]*tmp1[35]+z[82]*tmp1[36]+z[25]*tmp1[37]+z[65]*tmp1[38]+z[8]*tmp1[39]+z[48]*tmp1[40]
            +z[88]*tmp1[41]+z[31]*tmp1[42]+z[71]*tmp1[43]+z[14]*tmp1[44]+z[54]*tmp1[45]+z[94]*tmp1[46]+z[37]*tmp1[47]+z[77]*tmp1[48]+z[20]*tmp1[49]+z[60]*tmp1[50]
            +z[3]*tmp1[51]+z[43]*tmp1[52]+z[83]*tmp1[53]+z[26]*tmp1[54]+z[66]*tmp1[55]+z[9]*tmp1[56]+z[49]*tmp1[57]+z[89]*tmp1[58]+z[32]*tmp1[59]+z[72]*tmp1[60]
            +z[15]*tmp1[61]+z[55]*tmp1[62]+z[95]*tmp1[63]+z[38]*tmp1[64]+z[78]*tmp1[65]+z[21]*tmp1[66]+z[61]*tmp1[67]+z[4]*tmp1[68]+z[44]*tmp1[69]+z[84]*tmp1[70]
            +z[27]*tmp1[71]+z[67]*tmp1[72]+z[10]*tmp1[73]+z[50]*tmp1[74]+z[90]*tmp1[75]+z[33]*tmp1[76]+z[73]*tmp1[77]+z[16]*tmp1[78]+z[56]*tmp1[79]+z[96]*tmp1[80]
            +z[39]*tmp1[81]+z[79]*tmp1[82]+z[22]*tmp1[83]+z[62]*tmp1[84]+z[5]*tmp1[85]+z[45]*tmp1[86]+z[85]*tmp1[87]+z[28]*tmp1[88]+z[68]*tmp1[89]+z[11]*tmp1[90]
            +z[51]*tmp1[91]+z[91]*tmp1[92]+z[34]*tmp1[93]+z[74]*tmp1[94]+z[17]*tmp1[95]+z[57]*tmp1[96];
            tab[i+41*stg_first]=tab[i+41*stg_first]+z[0]*tmp1[0]
            +z[41]*tmp1[1]+z[82]*tmp1[2]+z[26]*tmp1[3]+z[67]*tmp1[4]+z[11]*tmp1[5]+z[52]*tmp1[6]+z[93]*tmp1[7]+z[37]*tmp1[8]+z[78]*tmp1[9]+z[22]*tmp1[10]
            +z[63]*tmp1[11]+z[7]*tmp1[12]+z[48]*tmp1[13]+z[89]*tmp1[14]+z[33]*tmp1[15]+z[74]*tmp1[16]+z[18]*tmp1[17]+z[59]*tmp1[18]+z[3]*tmp1[19]+z[44]*tmp1[20]
            +z[85]*tmp1[21]+z[29]*tmp1[22]+z[70]*tmp1[23]+z[14]*tmp1[24]+z[55]*tmp1[25]+z[96]*tmp1[26]+z[40]*tmp1[27]+z[81]*tmp1[28]+z[25]*tmp1[29]+z[66]*tmp1[30]
            +z[10]*tmp1[31]+z[51]*tmp1[32]+z[92]*tmp1[33]+z[36]*tmp1[34]+z[77]*tmp1[35]+z[21]*tmp1[36]+z[62]*tmp1[37]+z[6]*tmp1[38]+z[47]*tmp1[39]+z[88]*tmp1[40]
            +z[32]*tmp1[41]+z[73]*tmp1[42]+z[17]*tmp1[43]+z[58]*tmp1[44]+z[2]*tmp1[45]+z[43]*tmp1[46]+z[84]*tmp1[47]+z[28]*tmp1[48]+z[69]*tmp1[49]+z[13]*tmp1[50]
            +z[54]*tmp1[51]+z[95]*tmp1[52]+z[39]*tmp1[53]+z[80]*tmp1[54]+z[24]*tmp1[55]+z[65]*tmp1[56]+z[9]*tmp1[57]+z[50]*tmp1[58]+z[91]*tmp1[59]+z[35]*tmp1[60]
            +z[76]*tmp1[61]+z[20]*tmp1[62]+z[61]*tmp1[63]+z[5]*tmp1[64]+z[46]*tmp1[65]+z[87]*tmp1[66]+z[31]*tmp1[67]+z[72]*tmp1[68]+z[16]*tmp1[69]+z[57]*tmp1[70]
            +z[1]*tmp1[71]+z[42]*tmp1[72]+z[83]*tmp1[73]+z[27]*tmp1[74]+z[68]*tmp1[75]+z[12]*tmp1[76]+z[53]*tmp1[77]+z[94]*tmp1[78]+z[38]*tmp1[79]+z[79]*tmp1[80]
            +z[23]*tmp1[81]+z[64]*tmp1[82]+z[8]*tmp1[83]+z[49]*tmp1[84]+z[90]*tmp1[85]+z[34]*tmp1[86]+z[75]*tmp1[87]+z[19]*tmp1[88]+z[60]*tmp1[89]+z[4]*tmp1[90]
            +z[45]*tmp1[91]+z[86]*tmp1[92]+z[30]*tmp1[93]+z[71]*tmp1[94]+z[15]*tmp1[95]+z[56]*tmp1[96];
            tab[i+42*stg_first]=tab[i+42*stg_first]+z[0]*tmp1[0]
            +z[42]*tmp1[1]+z[84]*tmp1[2]+z[29]*tmp1[3]+z[71]*tmp1[4]+z[16]*tmp1[5]+z[58]*tmp1[6]+z[3]*tmp1[7]+z[45]*tmp1[8]+z[87]*tmp1[9]+z[32]*tmp1[10]
            +z[74]*tmp1[11]+z[19]*tmp1[12]+z[61]*tmp1[13]+z[6]*tmp1[14]+z[48]*tmp1[15]+z[90]*tmp1[16]+z[35]*tmp1[17]+z[77]*tmp1[18]+z[22]*tmp1[19]+z[64]*tmp1[20]
            +z[9]*tmp1[21]+z[51]*tmp1[22]+z[93]*tmp1[23]+z[38]*tmp1[24]+z[80]*tmp1[25]+z[25]*tmp1[26]+z[67]*tmp1[27]+z[12]*tmp1[28]+z[54]*tmp1[29]+z[96]*tmp1[30]
            +z[41]*tmp1[31]+z[83]*tmp1[32]+z[28]*tmp1[33]+z[70]*tmp1[34]+z[15]*tmp1[35]+z[57]*tmp1[36]+z[2]*tmp1[37]+z[44]*tmp1[38]+z[86]*tmp1[39]+z[31]*tmp1[40]
            +z[73]*tmp1[41]+z[18]*tmp1[42]+z[60]*tmp1[43]+z[5]*tmp1[44]+z[47]*tmp1[45]+z[89]*tmp1[46]+z[34]*tmp1[47]+z[76]*tmp1[48]+z[21]*tmp1[49]+z[63]*tmp1[50]
            +z[8]*tmp1[51]+z[50]*tmp1[52]+z[92]*tmp1[53]+z[37]*tmp1[54]+z[79]*tmp1[55]+z[24]*tmp1[56]+z[66]*tmp1[57]+z[11]*tmp1[58]+z[53]*tmp1[59]+z[95]*tmp1[60]
            +z[40]*tmp1[61]+z[82]*tmp1[62]+z[27]*tmp1[63]+z[69]*tmp1[64]+z[14]*tmp1[65]+z[56]*tmp1[66]+z[1]*tmp1[67]+z[43]*tmp1[68]+z[85]*tmp1[69]+z[30]*tmp1[70]
            +z[72]*tmp1[71]+z[17]*tmp1[72]+z[59]*tmp1[73]+z[4]*tmp1[74]+z[46]*tmp1[75]+z[88]*tmp1[76]+z[33]*tmp1[77]+z[75]*tmp1[78]+z[20]*tmp1[79]+z[62]*tmp1[80]
            +z[7]*tmp1[81]+z[49]*tmp1[82]+z[91]*tmp1[83]+z[36]*tmp1[84]+z[78]*tmp1[85]+z[23]*tmp1[86]+z[65]*tmp1[87]+z[10]*tmp1[88]+z[52]*tmp1[89]+z[94]*tmp1[90]
            +z[39]*tmp1[91]+z[81]*tmp1[92]+z[26]*tmp1[93]+z[68]*tmp1[94]+z[13]*tmp1[95]+z[55]*tmp1[96];
            tab[i+43*stg_first]=tab[i+43*stg_first]+z[0]*tmp1[0]
            +z[43]*tmp1[1]+z[86]*tmp1[2]+z[32]*tmp1[3]+z[75]*tmp1[4]+z[21]*tmp1[5]+z[64]*tmp1[6]+z[10]*tmp1[7]+z[53]*tmp1[8]+z[96]*tmp1[9]+z[42]*tmp1[10]
            +z[85]*tmp1[11]+z[31]*tmp1[12]+z[74]*tmp1[13]+z[20]*tmp1[14]+z[63]*tmp1[15]+z[9]*tmp1[16]+z[52]*tmp1[17]+z[95]*tmp1[18]+z[41]*tmp1[19]+z[84]*tmp1[20]
            +z[30]*tmp1[21]+z[73]*tmp1[22]+z[19]*tmp1[23]+z[62]*tmp1[24]+z[8]*tmp1[25]+z[51]*tmp1[26]+z[94]*tmp1[27]+z[40]*tmp1[28]+z[83]*tmp1[29]+z[29]*tmp1[30]
            +z[72]*tmp1[31]+z[18]*tmp1[32]+z[61]*tmp1[33]+z[7]*tmp1[34]+z[50]*tmp1[35]+z[93]*tmp1[36]+z[39]*tmp1[37]+z[82]*tmp1[38]+z[28]*tmp1[39]+z[71]*tmp1[40]
            +z[17]*tmp1[41]+z[60]*tmp1[42]+z[6]*tmp1[43]+z[49]*tmp1[44]+z[92]*tmp1[45]+z[38]*tmp1[46]+z[81]*tmp1[47]+z[27]*tmp1[48]+z[70]*tmp1[49]+z[16]*tmp1[50]
            +z[59]*tmp1[51]+z[5]*tmp1[52]+z[48]*tmp1[53]+z[91]*tmp1[54]+z[37]*tmp1[55]+z[80]*tmp1[56]+z[26]*tmp1[57]+z[69]*tmp1[58]+z[15]*tmp1[59]+z[58]*tmp1[60]
            +z[4]*tmp1[61]+z[47]*tmp1[62]+z[90]*tmp1[63]+z[36]*tmp1[64]+z[79]*tmp1[65]+z[25]*tmp1[66]+z[68]*tmp1[67]+z[14]*tmp1[68]+z[57]*tmp1[69]+z[3]*tmp1[70]
            +z[46]*tmp1[71]+z[89]*tmp1[72]+z[35]*tmp1[73]+z[78]*tmp1[74]+z[24]*tmp1[75]+z[67]*tmp1[76]+z[13]*tmp1[77]+z[56]*tmp1[78]+z[2]*tmp1[79]+z[45]*tmp1[80]
            +z[88]*tmp1[81]+z[34]*tmp1[82]+z[77]*tmp1[83]+z[23]*tmp1[84]+z[66]*tmp1[85]+z[12]*tmp1[86]+z[55]*tmp1[87]+z[1]*tmp1[88]+z[44]*tmp1[89]+z[87]*tmp1[90]
            +z[33]*tmp1[91]+z[76]*tmp1[92]+z[22]*tmp1[93]+z[65]*tmp1[94]+z[11]*tmp1[95]+z[54]*tmp1[96];
            tab[i+44*stg_first]=tab[i+44*stg_first]+z[0]*tmp1[0]
            +z[44]*tmp1[1]+z[88]*tmp1[2]+z[35]*tmp1[3]+z[79]*tmp1[4]+z[26]*tmp1[5]+z[70]*tmp1[6]+z[17]*tmp1[7]+z[61]*tmp1[8]+z[8]*tmp1[9]+z[52]*tmp1[10]
            +z[96]*tmp1[11]+z[43]*tmp1[12]+z[87]*tmp1[13]+z[34]*tmp1[14]+z[78]*tmp1[15]+z[25]*tmp1[16]+z[69]*tmp1[17]+z[16]*tmp1[18]+z[60]*tmp1[19]+z[7]*tmp1[20]
            +z[51]*tmp1[21]+z[95]*tmp1[22]+z[42]*tmp1[23]+z[86]*tmp1[24]+z[33]*tmp1[25]+z[77]*tmp1[26]+z[24]*tmp1[27]+z[68]*tmp1[28]+z[15]*tmp1[29]+z[59]*tmp1[30]
            +z[6]*tmp1[31]+z[50]*tmp1[32]+z[94]*tmp1[33]+z[41]*tmp1[34]+z[85]*tmp1[35]+z[32]*tmp1[36]+z[76]*tmp1[37]+z[23]*tmp1[38]+z[67]*tmp1[39]+z[14]*tmp1[40]
            +z[58]*tmp1[41]+z[5]*tmp1[42]+z[49]*tmp1[43]+z[93]*tmp1[44]+z[40]*tmp1[45]+z[84]*tmp1[46]+z[31]*tmp1[47]+z[75]*tmp1[48]+z[22]*tmp1[49]+z[66]*tmp1[50]
            +z[13]*tmp1[51]+z[57]*tmp1[52]+z[4]*tmp1[53]+z[48]*tmp1[54]+z[92]*tmp1[55]+z[39]*tmp1[56]+z[83]*tmp1[57]+z[30]*tmp1[58]+z[74]*tmp1[59]+z[21]*tmp1[60]
            +z[65]*tmp1[61]+z[12]*tmp1[62]+z[56]*tmp1[63]+z[3]*tmp1[64]+z[47]*tmp1[65]+z[91]*tmp1[66]+z[38]*tmp1[67]+z[82]*tmp1[68]+z[29]*tmp1[69]+z[73]*tmp1[70]
            +z[20]*tmp1[71]+z[64]*tmp1[72]+z[11]*tmp1[73]+z[55]*tmp1[74]+z[2]*tmp1[75]+z[46]*tmp1[76]+z[90]*tmp1[77]+z[37]*tmp1[78]+z[81]*tmp1[79]+z[28]*tmp1[80]
            +z[72]*tmp1[81]+z[19]*tmp1[82]+z[63]*tmp1[83]+z[10]*tmp1[84]+z[54]*tmp1[85]+z[1]*tmp1[86]+z[45]*tmp1[87]+z[89]*tmp1[88]+z[36]*tmp1[89]+z[80]*tmp1[90]
            +z[27]*tmp1[91]+z[71]*tmp1[92]+z[18]*tmp1[93]+z[62]*tmp1[94]+z[9]*tmp1[95]+z[53]*tmp1[96];
            tab[i+45*stg_first]=tab[i+45*stg_first]+z[0]*tmp1[0]
            +z[45]*tmp1[1]+z[90]*tmp1[2]+z[38]*tmp1[3]+z[83]*tmp1[4]+z[31]*tmp1[5]+z[76]*tmp1[6]+z[24]*tmp1[7]+z[69]*tmp1[8]+z[17]*tmp1[9]+z[62]*tmp1[10]
            +z[10]*tmp1[11]+z[55]*tmp1[12]+z[3]*tmp1[13]+z[48]*tmp1[14]+z[93]*tmp1[15]+z[41]*tmp1[16]+z[86]*tmp1[17]+z[34]*tmp1[18]+z[79]*tmp1[19]+z[27]*tmp1[20]
            +z[72]*tmp1[21]+z[20]*tmp1[22]+z[65]*tmp1[23]+z[13]*tmp1[24]+z[58]*tmp1[25]+z[6]*tmp1[26]+z[51]*tmp1[27]+z[96]*tmp1[28]+z[44]*tmp1[29]+z[89]*tmp1[30]
            +z[37]*tmp1[31]+z[82]*tmp1[32]+z[30]*tmp1[33]+z[75]*tmp1[34]+z[23]*tmp1[35]+z[68]*tmp1[36]+z[16]*tmp1[37]+z[61]*tmp1[38]+z[9]*tmp1[39]+z[54]*tmp1[40]
            +z[2]*tmp1[41]+z[47]*tmp1[42]+z[92]*tmp1[43]+z[40]*tmp1[44]+z[85]*tmp1[45]+z[33]*tmp1[46]+z[78]*tmp1[47]+z[26]*tmp1[48]+z[71]*tmp1[49]+z[19]*tmp1[50]
            +z[64]*tmp1[51]+z[12]*tmp1[52]+z[57]*tmp1[53]+z[5]*tmp1[54]+z[50]*tmp1[55]+z[95]*tmp1[56]+z[43]*tmp1[57]+z[88]*tmp1[58]+z[36]*tmp1[59]+z[81]*tmp1[60]
            +z[29]*tmp1[61]+z[74]*tmp1[62]+z[22]*tmp1[63]+z[67]*tmp1[64]+z[15]*tmp1[65]+z[60]*tmp1[66]+z[8]*tmp1[67]+z[53]*tmp1[68]+z[1]*tmp1[69]+z[46]*tmp1[70]
            +z[91]*tmp1[71]+z[39]*tmp1[72]+z[84]*tmp1[73]+z[32]*tmp1[74]+z[77]*tmp1[75]+z[25]*tmp1[76]+z[70]*tmp1[77]+z[18]*tmp1[78]+z[63]*tmp1[79]+z[11]*tmp1[80]
            +z[56]*tmp1[81]+z[4]*tmp1[82]+z[49]*tmp1[83]+z[94]*tmp1[84]+z[42]*tmp1[85]+z[87]*tmp1[86]+z[35]*tmp1[87]+z[80]*tmp1[88]+z[28]*tmp1[89]+z[73]*tmp1[90]
            +z[21]*tmp1[91]+z[66]*tmp1[92]+z[14]*tmp1[93]+z[59]*tmp1[94]+z[7]*tmp1[95]+z[52]*tmp1[96];
            tab[i+46*stg_first]=tab[i+46*stg_first]+z[0]*tmp1[0]
            +z[46]*tmp1[1]+z[92]*tmp1[2]+z[41]*tmp1[3]+z[87]*tmp1[4]+z[36]*tmp1[5]+z[82]*tmp1[6]+z[31]*tmp1[7]+z[77]*tmp1[8]+z[26]*tmp1[9]+z[72]*tmp1[10]
            +z[21]*tmp1[11]+z[67]*tmp1[12]+z[16]*tmp1[13]+z[62]*tmp1[14]+z[11]*tmp1[15]+z[57]*tmp1[16]+z[6]*tmp1[17]+z[52]*tmp1[18]+z[1]*tmp1[19]+z[47]*tmp1[20]
            +z[93]*tmp1[21]+z[42]*tmp1[22]+z[88]*tmp1[23]+z[37]*tmp1[24]+z[83]*tmp1[25]+z[32]*tmp1[26]+z[78]*tmp1[27]+z[27]*tmp1[28]+z[73]*tmp1[29]+z[22]*tmp1[30]
            +z[68]*tmp1[31]+z[17]*tmp1[32]+z[63]*tmp1[33]+z[12]*tmp1[34]+z[58]*tmp1[35]+z[7]*tmp1[36]+z[53]*tmp1[37]+z[2]*tmp1[38]+z[48]*tmp1[39]+z[94]*tmp1[40]
            +z[43]*tmp1[41]+z[89]*tmp1[42]+z[38]*tmp1[43]+z[84]*tmp1[44]+z[33]*tmp1[45]+z[79]*tmp1[46]+z[28]*tmp1[47]+z[74]*tmp1[48]+z[23]*tmp1[49]+z[69]*tmp1[50]
            +z[18]*tmp1[51]+z[64]*tmp1[52]+z[13]*tmp1[53]+z[59]*tmp1[54]+z[8]*tmp1[55]+z[54]*tmp1[56]+z[3]*tmp1[57]+z[49]*tmp1[58]+z[95]*tmp1[59]+z[44]*tmp1[60]
            +z[90]*tmp1[61]+z[39]*tmp1[62]+z[85]*tmp1[63]+z[34]*tmp1[64]+z[80]*tmp1[65]+z[29]*tmp1[66]+z[75]*tmp1[67]+z[24]*tmp1[68]+z[70]*tmp1[69]+z[19]*tmp1[70]
            +z[65]*tmp1[71]+z[14]*tmp1[72]+z[60]*tmp1[73]+z[9]*tmp1[74]+z[55]*tmp1[75]+z[4]*tmp1[76]+z[50]*tmp1[77]+z[96]*tmp1[78]+z[45]*tmp1[79]+z[91]*tmp1[80]
            +z[40]*tmp1[81]+z[86]*tmp1[82]+z[35]*tmp1[83]+z[81]*tmp1[84]+z[30]*tmp1[85]+z[76]*tmp1[86]+z[25]*tmp1[87]+z[71]*tmp1[88]+z[20]*tmp1[89]+z[66]*tmp1[90]
            +z[15]*tmp1[91]+z[61]*tmp1[92]+z[10]*tmp1[93]+z[56]*tmp1[94]+z[5]*tmp1[95]+z[51]*tmp1[96];
            tab[i+47*stg_first]=tab[i+47*stg_first]+z[0]*tmp1[0]
            +z[47]*tmp1[1]+z[94]*tmp1[2]+z[44]*tmp1[3]+z[91]*tmp1[4]+z[41]*tmp1[5]+z[88]*tmp1[6]+z[38]*tmp1[7]+z[85]*tmp1[8]+z[35]*tmp1[9]+z[82]*tmp1[10]
            +z[32]*tmp1[11]+z[79]*tmp1[12]+z[29]*tmp1[13]+z[76]*tmp1[14]+z[26]*tmp1[15]+z[73]*tmp1[16]+z[23]*tmp1[17]+z[70]*tmp1[18]+z[20]*tmp1[19]+z[67]*tmp1[20]
            +z[17]*tmp1[21]+z[64]*tmp1[22]+z[14]*tmp1[23]+z[61]*tmp1[24]+z[11]*tmp1[25]+z[58]*tmp1[26]+z[8]*tmp1[27]+z[55]*tmp1[28]+z[5]*tmp1[29]+z[52]*tmp1[30]
            +z[2]*tmp1[31]+z[49]*tmp1[32]+z[96]*tmp1[33]+z[46]*tmp1[34]+z[93]*tmp1[35]+z[43]*tmp1[36]+z[90]*tmp1[37]+z[40]*tmp1[38]+z[87]*tmp1[39]+z[37]*tmp1[40]
            +z[84]*tmp1[41]+z[34]*tmp1[42]+z[81]*tmp1[43]+z[31]*tmp1[44]+z[78]*tmp1[45]+z[28]*tmp1[46]+z[75]*tmp1[47]+z[25]*tmp1[48]+z[72]*tmp1[49]+z[22]*tmp1[50]
            +z[69]*tmp1[51]+z[19]*tmp1[52]+z[66]*tmp1[53]+z[16]*tmp1[54]+z[63]*tmp1[55]+z[13]*tmp1[56]+z[60]*tmp1[57]+z[10]*tmp1[58]+z[57]*tmp1[59]+z[7]*tmp1[60]
            +z[54]*tmp1[61]+z[4]*tmp1[62]+z[51]*tmp1[63]+z[1]*tmp1[64]+z[48]*tmp1[65]+z[95]*tmp1[66]+z[45]*tmp1[67]+z[92]*tmp1[68]+z[42]*tmp1[69]+z[89]*tmp1[70]
            +z[39]*tmp1[71]+z[86]*tmp1[72]+z[36]*tmp1[73]+z[83]*tmp1[74]+z[33]*tmp1[75]+z[80]*tmp1[76]+z[30]*tmp1[77]+z[77]*tmp1[78]+z[27]*tmp1[79]+z[74]*tmp1[80]
            +z[24]*tmp1[81]+z[71]*tmp1[82]+z[21]*tmp1[83]+z[68]*tmp1[84]+z[18]*tmp1[85]+z[65]*tmp1[86]+z[15]*tmp1[87]+z[62]*tmp1[88]+z[12]*tmp1[89]+z[59]*tmp1[90]
            +z[9]*tmp1[91]+z[56]*tmp1[92]+z[6]*tmp1[93]+z[53]*tmp1[94]+z[3]*tmp1[95]+z[50]*tmp1[96];
            tab[i+48*stg_first]=tab[i+48*stg_first]+z[0]*tmp1[0]
            +z[48]*tmp1[1]+z[96]*tmp1[2]+z[47]*tmp1[3]+z[95]*tmp1[4]+z[46]*tmp1[5]+z[94]*tmp1[6]+z[45]*tmp1[7]+z[93]*tmp1[8]+z[44]*tmp1[9]+z[92]*tmp1[10]
            +z[43]*tmp1[11]+z[91]*tmp1[12]+z[42]*tmp1[13]+z[90]*tmp1[14]+z[41]*tmp1[15]+z[89]*tmp1[16]+z[40]*tmp1[17]+z[88]*tmp1[18]+z[39]*tmp1[19]+z[87]*tmp1[20]
            +z[38]*tmp1[21]+z[86]*tmp1[22]+z[37]*tmp1[23]+z[85]*tmp1[24]+z[36]*tmp1[25]+z[84]*tmp1[26]+z[35]*tmp1[27]+z[83]*tmp1[28]+z[34]*tmp1[29]+z[82]*tmp1[30]
            +z[33]*tmp1[31]+z[81]*tmp1[32]+z[32]*tmp1[33]+z[80]*tmp1[34]+z[31]*tmp1[35]+z[79]*tmp1[36]+z[30]*tmp1[37]+z[78]*tmp1[38]+z[29]*tmp1[39]+z[77]*tmp1[40]
            +z[28]*tmp1[41]+z[76]*tmp1[42]+z[27]*tmp1[43]+z[75]*tmp1[44]+z[26]*tmp1[45]+z[74]*tmp1[46]+z[25]*tmp1[47]+z[73]*tmp1[48]+z[24]*tmp1[49]+z[72]*tmp1[50]
            +z[23]*tmp1[51]+z[71]*tmp1[52]+z[22]*tmp1[53]+z[70]*tmp1[54]+z[21]*tmp1[55]+z[69]*tmp1[56]+z[20]*tmp1[57]+z[68]*tmp1[58]+z[19]*tmp1[59]+z[67]*tmp1[60]
            +z[18]*tmp1[61]+z[66]*tmp1[62]+z[17]*tmp1[63]+z[65]*tmp1[64]+z[16]*tmp1[65]+z[64]*tmp1[66]+z[15]*tmp1[67]+z[63]*tmp1[68]+z[14]*tmp1[69]+z[62]*tmp1[70]
            +z[13]*tmp1[71]+z[61]*tmp1[72]+z[12]*tmp1[73]+z[60]*tmp1[74]+z[11]*tmp1[75]+z[59]*tmp1[76]+z[10]*tmp1[77]+z[58]*tmp1[78]+z[9]*tmp1[79]+z[57]*tmp1[80]
            +z[8]*tmp1[81]+z[56]*tmp1[82]+z[7]*tmp1[83]+z[55]*tmp1[84]+z[6]*tmp1[85]+z[54]*tmp1[86]+z[5]*tmp1[87]+z[53]*tmp1[88]+z[4]*tmp1[89]+z[52]*tmp1[90]
            +z[3]*tmp1[91]+z[51]*tmp1[92]+z[2]*tmp1[93]+z[50]*tmp1[94]+z[1]*tmp1[95]+z[49]*tmp1[96];
            tab[i+49*stg_first]=tab[i+49*stg_first]+z[0]*tmp1[0]
            +z[49]*tmp1[1]+z[1]*tmp1[2]+z[50]*tmp1[3]+z[2]*tmp1[4]+z[51]*tmp1[5]+z[3]*tmp1[6]+z[52]*tmp1[7]+z[4]*tmp1[8]+z[53]*tmp1[9]+z[5]*tmp1[10]
            +z[54]*tmp1[11]+z[6]*tmp1[12]+z[55]*tmp1[13]+z[7]*tmp1[14]+z[56]*tmp1[15]+z[8]*tmp1[16]+z[57]*tmp1[17]+z[9]*tmp1[18]+z[58]*tmp1[19]+z[10]*tmp1[20]
            +z[59]*tmp1[21]+z[11]*tmp1[22]+z[60]*tmp1[23]+z[12]*tmp1[24]+z[61]*tmp1[25]+z[13]*tmp1[26]+z[62]*tmp1[27]+z[14]*tmp1[28]+z[63]*tmp1[29]+z[15]*tmp1[30]
            +z[64]*tmp1[31]+z[16]*tmp1[32]+z[65]*tmp1[33]+z[17]*tmp1[34]+z[66]*tmp1[35]+z[18]*tmp1[36]+z[67]*tmp1[37]+z[19]*tmp1[38]+z[68]*tmp1[39]+z[20]*tmp1[40]
            +z[69]*tmp1[41]+z[21]*tmp1[42]+z[70]*tmp1[43]+z[22]*tmp1[44]+z[71]*tmp1[45]+z[23]*tmp1[46]+z[72]*tmp1[47]+z[24]*tmp1[48]+z[73]*tmp1[49]+z[25]*tmp1[50]
            +z[74]*tmp1[51]+z[26]*tmp1[52]+z[75]*tmp1[53]+z[27]*tmp1[54]+z[76]*tmp1[55]+z[28]*tmp1[56]+z[77]*tmp1[57]+z[29]*tmp1[58]+z[78]*tmp1[59]+z[30]*tmp1[60]
            +z[79]*tmp1[61]+z[31]*tmp1[62]+z[80]*tmp1[63]+z[32]*tmp1[64]+z[81]*tmp1[65]+z[33]*tmp1[66]+z[82]*tmp1[67]+z[34]*tmp1[68]+z[83]*tmp1[69]+z[35]*tmp1[70]
            +z[84]*tmp1[71]+z[36]*tmp1[72]+z[85]*tmp1[73]+z[37]*tmp1[74]+z[86]*tmp1[75]+z[38]*tmp1[76]+z[87]*tmp1[77]+z[39]*tmp1[78]+z[88]*tmp1[79]+z[40]*tmp1[80]
            +z[89]*tmp1[81]+z[41]*tmp1[82]+z[90]*tmp1[83]+z[42]*tmp1[84]+z[91]*tmp1[85]+z[43]*tmp1[86]+z[92]*tmp1[87]+z[44]*tmp1[88]+z[93]*tmp1[89]+z[45]*tmp1[90]
            +z[94]*tmp1[91]+z[46]*tmp1[92]+z[95]*tmp1[93]+z[47]*tmp1[94]+z[96]*tmp1[95]+z[48]*tmp1[96];
            tab[i+50*stg_first]=tab[i+50*stg_first]+z[0]*tmp1[0]
            +z[50]*tmp1[1]+z[3]*tmp1[2]+z[53]*tmp1[3]+z[6]*tmp1[4]+z[56]*tmp1[5]+z[9]*tmp1[6]+z[59]*tmp1[7]+z[12]*tmp1[8]+z[62]*tmp1[9]+z[15]*tmp1[10]
            +z[65]*tmp1[11]+z[18]*tmp1[12]+z[68]*tmp1[13]+z[21]*tmp1[14]+z[71]*tmp1[15]+z[24]*tmp1[16]+z[74]*tmp1[17]+z[27]*tmp1[18]+z[77]*tmp1[19]+z[30]*tmp1[20]
            +z[80]*tmp1[21]+z[33]*tmp1[22]+z[83]*tmp1[23]+z[36]*tmp1[24]+z[86]*tmp1[25]+z[39]*tmp1[26]+z[89]*tmp1[27]+z[42]*tmp1[28]+z[92]*tmp1[29]+z[45]*tmp1[30]
            +z[95]*tmp1[31]+z[48]*tmp1[32]+z[1]*tmp1[33]+z[51]*tmp1[34]+z[4]*tmp1[35]+z[54]*tmp1[36]+z[7]*tmp1[37]+z[57]*tmp1[38]+z[10]*tmp1[39]+z[60]*tmp1[40]
            +z[13]*tmp1[41]+z[63]*tmp1[42]+z[16]*tmp1[43]+z[66]*tmp1[44]+z[19]*tmp1[45]+z[69]*tmp1[46]+z[22]*tmp1[47]+z[72]*tmp1[48]+z[25]*tmp1[49]+z[75]*tmp1[50]
            +z[28]*tmp1[51]+z[78]*tmp1[52]+z[31]*tmp1[53]+z[81]*tmp1[54]+z[34]*tmp1[55]+z[84]*tmp1[56]+z[37]*tmp1[57]+z[87]*tmp1[58]+z[40]*tmp1[59]+z[90]*tmp1[60]
            +z[43]*tmp1[61]+z[93]*tmp1[62]+z[46]*tmp1[63]+z[96]*tmp1[64]+z[49]*tmp1[65]+z[2]*tmp1[66]+z[52]*tmp1[67]+z[5]*tmp1[68]+z[55]*tmp1[69]+z[8]*tmp1[70]
            +z[58]*tmp1[71]+z[11]*tmp1[72]+z[61]*tmp1[73]+z[14]*tmp1[74]+z[64]*tmp1[75]+z[17]*tmp1[76]+z[67]*tmp1[77]+z[20]*tmp1[78]+z[70]*tmp1[79]+z[23]*tmp1[80]
            +z[73]*tmp1[81]+z[26]*tmp1[82]+z[76]*tmp1[83]+z[29]*tmp1[84]+z[79]*tmp1[85]+z[32]*tmp1[86]+z[82]*tmp1[87]+z[35]*tmp1[88]+z[85]*tmp1[89]+z[38]*tmp1[90]
            +z[88]*tmp1[91]+z[41]*tmp1[92]+z[91]*tmp1[93]+z[44]*tmp1[94]+z[94]*tmp1[95]+z[47]*tmp1[96];
            tab[i+51*stg_first]=tab[i+51*stg_first]+z[0]*tmp1[0]
            +z[51]*tmp1[1]+z[5]*tmp1[2]+z[56]*tmp1[3]+z[10]*tmp1[4]+z[61]*tmp1[5]+z[15]*tmp1[6]+z[66]*tmp1[7]+z[20]*tmp1[8]+z[71]*tmp1[9]+z[25]*tmp1[10]
            +z[76]*tmp1[11]+z[30]*tmp1[12]+z[81]*tmp1[13]+z[35]*tmp1[14]+z[86]*tmp1[15]+z[40]*tmp1[16]+z[91]*tmp1[17]+z[45]*tmp1[18]+z[96]*tmp1[19]+z[50]*tmp1[20]
            +z[4]*tmp1[21]+z[55]*tmp1[22]+z[9]*tmp1[23]+z[60]*tmp1[24]+z[14]*tmp1[25]+z[65]*tmp1[26]+z[19]*tmp1[27]+z[70]*tmp1[28]+z[24]*tmp1[29]+z[75]*tmp1[30]
            +z[29]*tmp1[31]+z[80]*tmp1[32]+z[34]*tmp1[33]+z[85]*tmp1[34]+z[39]*tmp1[35]+z[90]*tmp1[36]+z[44]*tmp1[37]+z[95]*tmp1[38]+z[49]*tmp1[39]+z[3]*tmp1[40]
            +z[54]*tmp1[41]+z[8]*tmp1[42]+z[59]*tmp1[43]+z[13]*tmp1[44]+z[64]*tmp1[45]+z[18]*tmp1[46]+z[69]*tmp1[47]+z[23]*tmp1[48]+z[74]*tmp1[49]+z[28]*tmp1[50]
            +z[79]*tmp1[51]+z[33]*tmp1[52]+z[84]*tmp1[53]+z[38]*tmp1[54]+z[89]*tmp1[55]+z[43]*tmp1[56]+z[94]*tmp1[57]+z[48]*tmp1[58]+z[2]*tmp1[59]+z[53]*tmp1[60]
            +z[7]*tmp1[61]+z[58]*tmp1[62]+z[12]*tmp1[63]+z[63]*tmp1[64]+z[17]*tmp1[65]+z[68]*tmp1[66]+z[22]*tmp1[67]+z[73]*tmp1[68]+z[27]*tmp1[69]+z[78]*tmp1[70]
            +z[32]*tmp1[71]+z[83]*tmp1[72]+z[37]*tmp1[73]+z[88]*tmp1[74]+z[42]*tmp1[75]+z[93]*tmp1[76]+z[47]*tmp1[77]+z[1]*tmp1[78]+z[52]*tmp1[79]+z[6]*tmp1[80]
            +z[57]*tmp1[81]+z[11]*tmp1[82]+z[62]*tmp1[83]+z[16]*tmp1[84]+z[67]*tmp1[85]+z[21]*tmp1[86]+z[72]*tmp1[87]+z[26]*tmp1[88]+z[77]*tmp1[89]+z[31]*tmp1[90]
            +z[82]*tmp1[91]+z[36]*tmp1[92]+z[87]*tmp1[93]+z[41]*tmp1[94]+z[92]*tmp1[95]+z[46]*tmp1[96];
            tab[i+52*stg_first]=tab[i+52*stg_first]+z[0]*tmp1[0]
            +z[52]*tmp1[1]+z[7]*tmp1[2]+z[59]*tmp1[3]+z[14]*tmp1[4]+z[66]*tmp1[5]+z[21]*tmp1[6]+z[73]*tmp1[7]+z[28]*tmp1[8]+z[80]*tmp1[9]+z[35]*tmp1[10]
            +z[87]*tmp1[11]+z[42]*tmp1[12]+z[94]*tmp1[13]+z[49]*tmp1[14]+z[4]*tmp1[15]+z[56]*tmp1[16]+z[11]*tmp1[17]+z[63]*tmp1[18]+z[18]*tmp1[19]+z[70]*tmp1[20]
            +z[25]*tmp1[21]+z[77]*tmp1[22]+z[32]*tmp1[23]+z[84]*tmp1[24]+z[39]*tmp1[25]+z[91]*tmp1[26]+z[46]*tmp1[27]+z[1]*tmp1[28]+z[53]*tmp1[29]+z[8]*tmp1[30]
            +z[60]*tmp1[31]+z[15]*tmp1[32]+z[67]*tmp1[33]+z[22]*tmp1[34]+z[74]*tmp1[35]+z[29]*tmp1[36]+z[81]*tmp1[37]+z[36]*tmp1[38]+z[88]*tmp1[39]+z[43]*tmp1[40]
            +z[95]*tmp1[41]+z[50]*tmp1[42]+z[5]*tmp1[43]+z[57]*tmp1[44]+z[12]*tmp1[45]+z[64]*tmp1[46]+z[19]*tmp1[47]+z[71]*tmp1[48]+z[26]*tmp1[49]+z[78]*tmp1[50]
            +z[33]*tmp1[51]+z[85]*tmp1[52]+z[40]*tmp1[53]+z[92]*tmp1[54]+z[47]*tmp1[55]+z[2]*tmp1[56]+z[54]*tmp1[57]+z[9]*tmp1[58]+z[61]*tmp1[59]+z[16]*tmp1[60]
            +z[68]*tmp1[61]+z[23]*tmp1[62]+z[75]*tmp1[63]+z[30]*tmp1[64]+z[82]*tmp1[65]+z[37]*tmp1[66]+z[89]*tmp1[67]+z[44]*tmp1[68]+z[96]*tmp1[69]+z[51]*tmp1[70]
            +z[6]*tmp1[71]+z[58]*tmp1[72]+z[13]*tmp1[73]+z[65]*tmp1[74]+z[20]*tmp1[75]+z[72]*tmp1[76]+z[27]*tmp1[77]+z[79]*tmp1[78]+z[34]*tmp1[79]+z[86]*tmp1[80]
            +z[41]*tmp1[81]+z[93]*tmp1[82]+z[48]*tmp1[83]+z[3]*tmp1[84]+z[55]*tmp1[85]+z[10]*tmp1[86]+z[62]*tmp1[87]+z[17]*tmp1[88]+z[69]*tmp1[89]+z[24]*tmp1[90]
            +z[76]*tmp1[91]+z[31]*tmp1[92]+z[83]*tmp1[93]+z[38]*tmp1[94]+z[90]*tmp1[95]+z[45]*tmp1[96];
            tab[i+53*stg_first]=tab[i+53*stg_first]+z[0]*tmp1[0]
            +z[53]*tmp1[1]+z[9]*tmp1[2]+z[62]*tmp1[3]+z[18]*tmp1[4]+z[71]*tmp1[5]+z[27]*tmp1[6]+z[80]*tmp1[7]+z[36]*tmp1[8]+z[89]*tmp1[9]+z[45]*tmp1[10]
            +z[1]*tmp1[11]+z[54]*tmp1[12]+z[10]*tmp1[13]+z[63]*tmp1[14]+z[19]*tmp1[15]+z[72]*tmp1[16]+z[28]*tmp1[17]+z[81]*tmp1[18]+z[37]*tmp1[19]+z[90]*tmp1[20]
            +z[46]*tmp1[21]+z[2]*tmp1[22]+z[55]*tmp1[23]+z[11]*tmp1[24]+z[64]*tmp1[25]+z[20]*tmp1[26]+z[73]*tmp1[27]+z[29]*tmp1[28]+z[82]*tmp1[29]+z[38]*tmp1[30]
            +z[91]*tmp1[31]+z[47]*tmp1[32]+z[3]*tmp1[33]+z[56]*tmp1[34]+z[12]*tmp1[35]+z[65]*tmp1[36]+z[21]*tmp1[37]+z[74]*tmp1[38]+z[30]*tmp1[39]+z[83]*tmp1[40]
            +z[39]*tmp1[41]+z[92]*tmp1[42]+z[48]*tmp1[43]+z[4]*tmp1[44]+z[57]*tmp1[45]+z[13]*tmp1[46]+z[66]*tmp1[47]+z[22]*tmp1[48]+z[75]*tmp1[49]+z[31]*tmp1[50]
            +z[84]*tmp1[51]+z[40]*tmp1[52]+z[93]*tmp1[53]+z[49]*tmp1[54]+z[5]*tmp1[55]+z[58]*tmp1[56]+z[14]*tmp1[57]+z[67]*tmp1[58]+z[23]*tmp1[59]+z[76]*tmp1[60]
            +z[32]*tmp1[61]+z[85]*tmp1[62]+z[41]*tmp1[63]+z[94]*tmp1[64]+z[50]*tmp1[65]+z[6]*tmp1[66]+z[59]*tmp1[67]+z[15]*tmp1[68]+z[68]*tmp1[69]+z[24]*tmp1[70]
            +z[77]*tmp1[71]+z[33]*tmp1[72]+z[86]*tmp1[73]+z[42]*tmp1[74]+z[95]*tmp1[75]+z[51]*tmp1[76]+z[7]*tmp1[77]+z[60]*tmp1[78]+z[16]*tmp1[79]+z[69]*tmp1[80]
            +z[25]*tmp1[81]+z[78]*tmp1[82]+z[34]*tmp1[83]+z[87]*tmp1[84]+z[43]*tmp1[85]+z[96]*tmp1[86]+z[52]*tmp1[87]+z[8]*tmp1[88]+z[61]*tmp1[89]+z[17]*tmp1[90]
            +z[70]*tmp1[91]+z[26]*tmp1[92]+z[79]*tmp1[93]+z[35]*tmp1[94]+z[88]*tmp1[95]+z[44]*tmp1[96];
            tab[i+54*stg_first]=tab[i+54*stg_first]+z[0]*tmp1[0]
            +z[54]*tmp1[1]+z[11]*tmp1[2]+z[65]*tmp1[3]+z[22]*tmp1[4]+z[76]*tmp1[5]+z[33]*tmp1[6]+z[87]*tmp1[7]+z[44]*tmp1[8]+z[1]*tmp1[9]+z[55]*tmp1[10]
            +z[12]*tmp1[11]+z[66]*tmp1[12]+z[23]*tmp1[13]+z[77]*tmp1[14]+z[34]*tmp1[15]+z[88]*tmp1[16]+z[45]*tmp1[17]+z[2]*tmp1[18]+z[56]*tmp1[19]+z[13]*tmp1[20]
            +z[67]*tmp1[21]+z[24]*tmp1[22]+z[78]*tmp1[23]+z[35]*tmp1[24]+z[89]*tmp1[25]+z[46]*tmp1[26]+z[3]*tmp1[27]+z[57]*tmp1[28]+z[14]*tmp1[29]+z[68]*tmp1[30]
            +z[25]*tmp1[31]+z[79]*tmp1[32]+z[36]*tmp1[33]+z[90]*tmp1[34]+z[47]*tmp1[35]+z[4]*tmp1[36]+z[58]*tmp1[37]+z[15]*tmp1[38]+z[69]*tmp1[39]+z[26]*tmp1[40]
            +z[80]*tmp1[41]+z[37]*tmp1[42]+z[91]*tmp1[43]+z[48]*tmp1[44]+z[5]*tmp1[45]+z[59]*tmp1[46]+z[16]*tmp1[47]+z[70]*tmp1[48]+z[27]*tmp1[49]+z[81]*tmp1[50]
            +z[38]*tmp1[51]+z[92]*tmp1[52]+z[49]*tmp1[53]+z[6]*tmp1[54]+z[60]*tmp1[55]+z[17]*tmp1[56]+z[71]*tmp1[57]+z[28]*tmp1[58]+z[82]*tmp1[59]+z[39]*tmp1[60]
            +z[93]*tmp1[61]+z[50]*tmp1[62]+z[7]*tmp1[63]+z[61]*tmp1[64]+z[18]*tmp1[65]+z[72]*tmp1[66]+z[29]*tmp1[67]+z[83]*tmp1[68]+z[40]*tmp1[69]+z[94]*tmp1[70]
            +z[51]*tmp1[71]+z[8]*tmp1[72]+z[62]*tmp1[73]+z[19]*tmp1[74]+z[73]*tmp1[75]+z[30]*tmp1[76]+z[84]*tmp1[77]+z[41]*tmp1[78]+z[95]*tmp1[79]+z[52]*tmp1[80]
            +z[9]*tmp1[81]+z[63]*tmp1[82]+z[20]*tmp1[83]+z[74]*tmp1[84]+z[31]*tmp1[85]+z[85]*tmp1[86]+z[42]*tmp1[87]+z[96]*tmp1[88]+z[53]*tmp1[89]+z[10]*tmp1[90]
            +z[64]*tmp1[91]+z[21]*tmp1[92]+z[75]*tmp1[93]+z[32]*tmp1[94]+z[86]*tmp1[95]+z[43]*tmp1[96];
            tab[i+55*stg_first]=tab[i+55*stg_first]+z[0]*tmp1[0]
            +z[55]*tmp1[1]+z[13]*tmp1[2]+z[68]*tmp1[3]+z[26]*tmp1[4]+z[81]*tmp1[5]+z[39]*tmp1[6]+z[94]*tmp1[7]+z[52]*tmp1[8]+z[10]*tmp1[9]+z[65]*tmp1[10]
            +z[23]*tmp1[11]+z[78]*tmp1[12]+z[36]*tmp1[13]+z[91]*tmp1[14]+z[49]*tmp1[15]+z[7]*tmp1[16]+z[62]*tmp1[17]+z[20]*tmp1[18]+z[75]*tmp1[19]+z[33]*tmp1[20]
            +z[88]*tmp1[21]+z[46]*tmp1[22]+z[4]*tmp1[23]+z[59]*tmp1[24]+z[17]*tmp1[25]+z[72]*tmp1[26]+z[30]*tmp1[27]+z[85]*tmp1[28]+z[43]*tmp1[29]+z[1]*tmp1[30]
            +z[56]*tmp1[31]+z[14]*tmp1[32]+z[69]*tmp1[33]+z[27]*tmp1[34]+z[82]*tmp1[35]+z[40]*tmp1[36]+z[95]*tmp1[37]+z[53]*tmp1[38]+z[11]*tmp1[39]+z[66]*tmp1[40]
            +z[24]*tmp1[41]+z[79]*tmp1[42]+z[37]*tmp1[43]+z[92]*tmp1[44]+z[50]*tmp1[45]+z[8]*tmp1[46]+z[63]*tmp1[47]+z[21]*tmp1[48]+z[76]*tmp1[49]+z[34]*tmp1[50]
            +z[89]*tmp1[51]+z[47]*tmp1[52]+z[5]*tmp1[53]+z[60]*tmp1[54]+z[18]*tmp1[55]+z[73]*tmp1[56]+z[31]*tmp1[57]+z[86]*tmp1[58]+z[44]*tmp1[59]+z[2]*tmp1[60]
            +z[57]*tmp1[61]+z[15]*tmp1[62]+z[70]*tmp1[63]+z[28]*tmp1[64]+z[83]*tmp1[65]+z[41]*tmp1[66]+z[96]*tmp1[67]+z[54]*tmp1[68]+z[12]*tmp1[69]+z[67]*tmp1[70]
            +z[25]*tmp1[71]+z[80]*tmp1[72]+z[38]*tmp1[73]+z[93]*tmp1[74]+z[51]*tmp1[75]+z[9]*tmp1[76]+z[64]*tmp1[77]+z[22]*tmp1[78]+z[77]*tmp1[79]+z[35]*tmp1[80]
            +z[90]*tmp1[81]+z[48]*tmp1[82]+z[6]*tmp1[83]+z[61]*tmp1[84]+z[19]*tmp1[85]+z[74]*tmp1[86]+z[32]*tmp1[87]+z[87]*tmp1[88]+z[45]*tmp1[89]+z[3]*tmp1[90]
            +z[58]*tmp1[91]+z[16]*tmp1[92]+z[71]*tmp1[93]+z[29]*tmp1[94]+z[84]*tmp1[95]+z[42]*tmp1[96];
            tab[i+56*stg_first]=tab[i+56*stg_first]+z[0]*tmp1[0]
            +z[56]*tmp1[1]+z[15]*tmp1[2]+z[71]*tmp1[3]+z[30]*tmp1[4]+z[86]*tmp1[5]+z[45]*tmp1[6]+z[4]*tmp1[7]+z[60]*tmp1[8]+z[19]*tmp1[9]+z[75]*tmp1[10]
            +z[34]*tmp1[11]+z[90]*tmp1[12]+z[49]*tmp1[13]+z[8]*tmp1[14]+z[64]*tmp1[15]+z[23]*tmp1[16]+z[79]*tmp1[17]+z[38]*tmp1[18]+z[94]*tmp1[19]+z[53]*tmp1[20]
            +z[12]*tmp1[21]+z[68]*tmp1[22]+z[27]*tmp1[23]+z[83]*tmp1[24]+z[42]*tmp1[25]+z[1]*tmp1[26]+z[57]*tmp1[27]+z[16]*tmp1[28]+z[72]*tmp1[29]+z[31]*tmp1[30]
            +z[87]*tmp1[31]+z[46]*tmp1[32]+z[5]*tmp1[33]+z[61]*tmp1[34]+z[20]*tmp1[35]+z[76]*tmp1[36]+z[35]*tmp1[37]+z[91]*tmp1[38]+z[50]*tmp1[39]+z[9]*tmp1[40]
            +z[65]*tmp1[41]+z[24]*tmp1[42]+z[80]*tmp1[43]+z[39]*tmp1[44]+z[95]*tmp1[45]+z[54]*tmp1[46]+z[13]*tmp1[47]+z[69]*tmp1[48]+z[28]*tmp1[49]+z[84]*tmp1[50]
            +z[43]*tmp1[51]+z[2]*tmp1[52]+z[58]*tmp1[53]+z[17]*tmp1[54]+z[73]*tmp1[55]+z[32]*tmp1[56]+z[88]*tmp1[57]+z[47]*tmp1[58]+z[6]*tmp1[59]+z[62]*tmp1[60]
            +z[21]*tmp1[61]+z[77]*tmp1[62]+z[36]*tmp1[63]+z[92]*tmp1[64]+z[51]*tmp1[65]+z[10]*tmp1[66]+z[66]*tmp1[67]+z[25]*tmp1[68]+z[81]*tmp1[69]+z[40]*tmp1[70]
            +z[96]*tmp1[71]+z[55]*tmp1[72]+z[14]*tmp1[73]+z[70]*tmp1[74]+z[29]*tmp1[75]+z[85]*tmp1[76]+z[44]*tmp1[77]+z[3]*tmp1[78]+z[59]*tmp1[79]+z[18]*tmp1[80]
            +z[74]*tmp1[81]+z[33]*tmp1[82]+z[89]*tmp1[83]+z[48]*tmp1[84]+z[7]*tmp1[85]+z[63]*tmp1[86]+z[22]*tmp1[87]+z[78]*tmp1[88]+z[37]*tmp1[89]+z[93]*tmp1[90]
            +z[52]*tmp1[91]+z[11]*tmp1[92]+z[67]*tmp1[93]+z[26]*tmp1[94]+z[82]*tmp1[95]+z[41]*tmp1[96];
            tab[i+57*stg_first]=tab[i+57*stg_first]+z[0]*tmp1[0]
            +z[57]*tmp1[1]+z[17]*tmp1[2]+z[74]*tmp1[3]+z[34]*tmp1[4]+z[91]*tmp1[5]+z[51]*tmp1[6]+z[11]*tmp1[7]+z[68]*tmp1[8]+z[28]*tmp1[9]+z[85]*tmp1[10]
            +z[45]*tmp1[11]+z[5]*tmp1[12]+z[62]*tmp1[13]+z[22]*tmp1[14]+z[79]*tmp1[15]+z[39]*tmp1[16]+z[96]*tmp1[17]+z[56]*tmp1[18]+z[16]*tmp1[19]+z[73]*tmp1[20]
            +z[33]*tmp1[21]+z[90]*tmp1[22]+z[50]*tmp1[23]+z[10]*tmp1[24]+z[67]*tmp1[25]+z[27]*tmp1[26]+z[84]*tmp1[27]+z[44]*tmp1[28]+z[4]*tmp1[29]+z[61]*tmp1[30]
            +z[21]*tmp1[31]+z[78]*tmp1[32]+z[38]*tmp1[33]+z[95]*tmp1[34]+z[55]*tmp1[35]+z[15]*tmp1[36]+z[72]*tmp1[37]+z[32]*tmp1[38]+z[89]*tmp1[39]+z[49]*tmp1[40]
            +z[9]*tmp1[41]+z[66]*tmp1[42]+z[26]*tmp1[43]+z[83]*tmp1[44]+z[43]*tmp1[45]+z[3]*tmp1[46]+z[60]*tmp1[47]+z[20]*tmp1[48]+z[77]*tmp1[49]+z[37]*tmp1[50]
            +z[94]*tmp1[51]+z[54]*tmp1[52]+z[14]*tmp1[53]+z[71]*tmp1[54]+z[31]*tmp1[55]+z[88]*tmp1[56]+z[48]*tmp1[57]+z[8]*tmp1[58]+z[65]*tmp1[59]+z[25]*tmp1[60]
            +z[82]*tmp1[61]+z[42]*tmp1[62]+z[2]*tmp1[63]+z[59]*tmp1[64]+z[19]*tmp1[65]+z[76]*tmp1[66]+z[36]*tmp1[67]+z[93]*tmp1[68]+z[53]*tmp1[69]+z[13]*tmp1[70]
            +z[70]*tmp1[71]+z[30]*tmp1[72]+z[87]*tmp1[73]+z[47]*tmp1[74]+z[7]*tmp1[75]+z[64]*tmp1[76]+z[24]*tmp1[77]+z[81]*tmp1[78]+z[41]*tmp1[79]+z[1]*tmp1[80]
            +z[58]*tmp1[81]+z[18]*tmp1[82]+z[75]*tmp1[83]+z[35]*tmp1[84]+z[92]*tmp1[85]+z[52]*tmp1[86]+z[12]*tmp1[87]+z[69]*tmp1[88]+z[29]*tmp1[89]+z[86]*tmp1[90]
            +z[46]*tmp1[91]+z[6]*tmp1[92]+z[63]*tmp1[93]+z[23]*tmp1[94]+z[80]*tmp1[95]+z[40]*tmp1[96];
            tab[i+58*stg_first]=tab[i+58*stg_first]+z[0]*tmp1[0]
            +z[58]*tmp1[1]+z[19]*tmp1[2]+z[77]*tmp1[3]+z[38]*tmp1[4]+z[96]*tmp1[5]+z[57]*tmp1[6]+z[18]*tmp1[7]+z[76]*tmp1[8]+z[37]*tmp1[9]+z[95]*tmp1[10]
            +z[56]*tmp1[11]+z[17]*tmp1[12]+z[75]*tmp1[13]+z[36]*tmp1[14]+z[94]*tmp1[15]+z[55]*tmp1[16]+z[16]*tmp1[17]+z[74]*tmp1[18]+z[35]*tmp1[19]+z[93]*tmp1[20]
            +z[54]*tmp1[21]+z[15]*tmp1[22]+z[73]*tmp1[23]+z[34]*tmp1[24]+z[92]*tmp1[25]+z[53]*tmp1[26]+z[14]*tmp1[27]+z[72]*tmp1[28]+z[33]*tmp1[29]+z[91]*tmp1[30]
            +z[52]*tmp1[31]+z[13]*tmp1[32]+z[71]*tmp1[33]+z[32]*tmp1[34]+z[90]*tmp1[35]+z[51]*tmp1[36]+z[12]*tmp1[37]+z[70]*tmp1[38]+z[31]*tmp1[39]+z[89]*tmp1[40]
            +z[50]*tmp1[41]+z[11]*tmp1[42]+z[69]*tmp1[43]+z[30]*tmp1[44]+z[88]*tmp1[45]+z[49]*tmp1[46]+z[10]*tmp1[47]+z[68]*tmp1[48]+z[29]*tmp1[49]+z[87]*tmp1[50]
            +z[48]*tmp1[51]+z[9]*tmp1[52]+z[67]*tmp1[53]+z[28]*tmp1[54]+z[86]*tmp1[55]+z[47]*tmp1[56]+z[8]*tmp1[57]+z[66]*tmp1[58]+z[27]*tmp1[59]+z[85]*tmp1[60]
            +z[46]*tmp1[61]+z[7]*tmp1[62]+z[65]*tmp1[63]+z[26]*tmp1[64]+z[84]*tmp1[65]+z[45]*tmp1[66]+z[6]*tmp1[67]+z[64]*tmp1[68]+z[25]*tmp1[69]+z[83]*tmp1[70]
            +z[44]*tmp1[71]+z[5]*tmp1[72]+z[63]*tmp1[73]+z[24]*tmp1[74]+z[82]*tmp1[75]+z[43]*tmp1[76]+z[4]*tmp1[77]+z[62]*tmp1[78]+z[23]*tmp1[79]+z[81]*tmp1[80]
            +z[42]*tmp1[81]+z[3]*tmp1[82]+z[61]*tmp1[83]+z[22]*tmp1[84]+z[80]*tmp1[85]+z[41]*tmp1[86]+z[2]*tmp1[87]+z[60]*tmp1[88]+z[21]*tmp1[89]+z[79]*tmp1[90]
            +z[40]*tmp1[91]+z[1]*tmp1[92]+z[59]*tmp1[93]+z[20]*tmp1[94]+z[78]*tmp1[95]+z[39]*tmp1[96];
            tab[i+59*stg_first]=tab[i+59*stg_first]+z[0]*tmp1[0]
            +z[59]*tmp1[1]+z[21]*tmp1[2]+z[80]*tmp1[3]+z[42]*tmp1[4]+z[4]*tmp1[5]+z[63]*tmp1[6]+z[25]*tmp1[7]+z[84]*tmp1[8]+z[46]*tmp1[9]+z[8]*tmp1[10]
            +z[67]*tmp1[11]+z[29]*tmp1[12]+z[88]*tmp1[13]+z[50]*tmp1[14]+z[12]*tmp1[15]+z[71]*tmp1[16]+z[33]*tmp1[17]+z[92]*tmp1[18]+z[54]*tmp1[19]+z[16]*tmp1[20]
            +z[75]*tmp1[21]+z[37]*tmp1[22]+z[96]*tmp1[23]+z[58]*tmp1[24]+z[20]*tmp1[25]+z[79]*tmp1[26]+z[41]*tmp1[27]+z[3]*tmp1[28]+z[62]*tmp1[29]+z[24]*tmp1[30]
            +z[83]*tmp1[31]+z[45]*tmp1[32]+z[7]*tmp1[33]+z[66]*tmp1[34]+z[28]*tmp1[35]+z[87]*tmp1[36]+z[49]*tmp1[37]+z[11]*tmp1[38]+z[70]*tmp1[39]+z[32]*tmp1[40]
            +z[91]*tmp1[41]+z[53]*tmp1[42]+z[15]*tmp1[43]+z[74]*tmp1[44]+z[36]*tmp1[45]+z[95]*tmp1[46]+z[57]*tmp1[47]+z[19]*tmp1[48]+z[78]*tmp1[49]+z[40]*tmp1[50]
            +z[2]*tmp1[51]+z[61]*tmp1[52]+z[23]*tmp1[53]+z[82]*tmp1[54]+z[44]*tmp1[55]+z[6]*tmp1[56]+z[65]*tmp1[57]+z[27]*tmp1[58]+z[86]*tmp1[59]+z[48]*tmp1[60]
            +z[10]*tmp1[61]+z[69]*tmp1[62]+z[31]*tmp1[63]+z[90]*tmp1[64]+z[52]*tmp1[65]+z[14]*tmp1[66]+z[73]*tmp1[67]+z[35]*tmp1[68]+z[94]*tmp1[69]+z[56]*tmp1[70]
            +z[18]*tmp1[71]+z[77]*tmp1[72]+z[39]*tmp1[73]+z[1]*tmp1[74]+z[60]*tmp1[75]+z[22]*tmp1[76]+z[81]*tmp1[77]+z[43]*tmp1[78]+z[5]*tmp1[79]+z[64]*tmp1[80]
            +z[26]*tmp1[81]+z[85]*tmp1[82]+z[47]*tmp1[83]+z[9]*tmp1[84]+z[68]*tmp1[85]+z[30]*tmp1[86]+z[89]*tmp1[87]+z[51]*tmp1[88]+z[13]*tmp1[89]+z[72]*tmp1[90]
            +z[34]*tmp1[91]+z[93]*tmp1[92]+z[55]*tmp1[93]+z[17]*tmp1[94]+z[76]*tmp1[95]+z[38]*tmp1[96];
            tab[i+60*stg_first]=tab[i+60*stg_first]+z[0]*tmp1[0]
            +z[60]*tmp1[1]+z[23]*tmp1[2]+z[83]*tmp1[3]+z[46]*tmp1[4]+z[9]*tmp1[5]+z[69]*tmp1[6]+z[32]*tmp1[7]+z[92]*tmp1[8]+z[55]*tmp1[9]+z[18]*tmp1[10]
            +z[78]*tmp1[11]+z[41]*tmp1[12]+z[4]*tmp1[13]+z[64]*tmp1[14]+z[27]*tmp1[15]+z[87]*tmp1[16]+z[50]*tmp1[17]+z[13]*tmp1[18]+z[73]*tmp1[19]+z[36]*tmp1[20]
            +z[96]*tmp1[21]+z[59]*tmp1[22]+z[22]*tmp1[23]+z[82]*tmp1[24]+z[45]*tmp1[25]+z[8]*tmp1[26]+z[68]*tmp1[27]+z[31]*tmp1[28]+z[91]*tmp1[29]+z[54]*tmp1[30]
            +z[17]*tmp1[31]+z[77]*tmp1[32]+z[40]*tmp1[33]+z[3]*tmp1[34]+z[63]*tmp1[35]+z[26]*tmp1[36]+z[86]*tmp1[37]+z[49]*tmp1[38]+z[12]*tmp1[39]+z[72]*tmp1[40]
            +z[35]*tmp1[41]+z[95]*tmp1[42]+z[58]*tmp1[43]+z[21]*tmp1[44]+z[81]*tmp1[45]+z[44]*tmp1[46]+z[7]*tmp1[47]+z[67]*tmp1[48]+z[30]*tmp1[49]+z[90]*tmp1[50]
            +z[53]*tmp1[51]+z[16]*tmp1[52]+z[76]*tmp1[53]+z[39]*tmp1[54]+z[2]*tmp1[55]+z[62]*tmp1[56]+z[25]*tmp1[57]+z[85]*tmp1[58]+z[48]*tmp1[59]+z[11]*tmp1[60]
            +z[71]*tmp1[61]+z[34]*tmp1[62]+z[94]*tmp1[63]+z[57]*tmp1[64]+z[20]*tmp1[65]+z[80]*tmp1[66]+z[43]*tmp1[67]+z[6]*tmp1[68]+z[66]*tmp1[69]+z[29]*tmp1[70]
            +z[89]*tmp1[71]+z[52]*tmp1[72]+z[15]*tmp1[73]+z[75]*tmp1[74]+z[38]*tmp1[75]+z[1]*tmp1[76]+z[61]*tmp1[77]+z[24]*tmp1[78]+z[84]*tmp1[79]+z[47]*tmp1[80]
            +z[10]*tmp1[81]+z[70]*tmp1[82]+z[33]*tmp1[83]+z[93]*tmp1[84]+z[56]*tmp1[85]+z[19]*tmp1[86]+z[79]*tmp1[87]+z[42]*tmp1[88]+z[5]*tmp1[89]+z[65]*tmp1[90]
            +z[28]*tmp1[91]+z[88]*tmp1[92]+z[51]*tmp1[93]+z[14]*tmp1[94]+z[74]*tmp1[95]+z[37]*tmp1[96];
            tab[i+61*stg_first]=tab[i+61*stg_first]+z[0]*tmp1[0]
            +z[61]*tmp1[1]+z[25]*tmp1[2]+z[86]*tmp1[3]+z[50]*tmp1[4]+z[14]*tmp1[5]+z[75]*tmp1[6]+z[39]*tmp1[7]+z[3]*tmp1[8]+z[64]*tmp1[9]+z[28]*tmp1[10]
            +z[89]*tmp1[11]+z[53]*tmp1[12]+z[17]*tmp1[13]+z[78]*tmp1[14]+z[42]*tmp1[15]+z[6]*tmp1[16]+z[67]*tmp1[17]+z[31]*tmp1[18]+z[92]*tmp1[19]+z[56]*tmp1[20]
            +z[20]*tmp1[21]+z[81]*tmp1[22]+z[45]*tmp1[23]+z[9]*tmp1[24]+z[70]*tmp1[25]+z[34]*tmp1[26]+z[95]*tmp1[27]+z[59]*tmp1[28]+z[23]*tmp1[29]+z[84]*tmp1[30]
            +z[48]*tmp1[31]+z[12]*tmp1[32]+z[73]*tmp1[33]+z[37]*tmp1[34]+z[1]*tmp1[35]+z[62]*tmp1[36]+z[26]*tmp1[37]+z[87]*tmp1[38]+z[51]*tmp1[39]+z[15]*tmp1[40]
            +z[76]*tmp1[41]+z[40]*tmp1[42]+z[4]*tmp1[43]+z[65]*tmp1[44]+z[29]*tmp1[45]+z[90]*tmp1[46]+z[54]*tmp1[47]+z[18]*tmp1[48]+z[79]*tmp1[49]+z[43]*tmp1[50]
            +z[7]*tmp1[51]+z[68]*tmp1[52]+z[32]*tmp1[53]+z[93]*tmp1[54]+z[57]*tmp1[55]+z[21]*tmp1[56]+z[82]*tmp1[57]+z[46]*tmp1[58]+z[10]*tmp1[59]+z[71]*tmp1[60]
            +z[35]*tmp1[61]+z[96]*tmp1[62]+z[60]*tmp1[63]+z[24]*tmp1[64]+z[85]*tmp1[65]+z[49]*tmp1[66]+z[13]*tmp1[67]+z[74]*tmp1[68]+z[38]*tmp1[69]+z[2]*tmp1[70]
            +z[63]*tmp1[71]+z[27]*tmp1[72]+z[88]*tmp1[73]+z[52]*tmp1[74]+z[16]*tmp1[75]+z[77]*tmp1[76]+z[41]*tmp1[77]+z[5]*tmp1[78]+z[66]*tmp1[79]+z[30]*tmp1[80]
            +z[91]*tmp1[81]+z[55]*tmp1[82]+z[19]*tmp1[83]+z[80]*tmp1[84]+z[44]*tmp1[85]+z[8]*tmp1[86]+z[69]*tmp1[87]+z[33]*tmp1[88]+z[94]*tmp1[89]+z[58]*tmp1[90]
            +z[22]*tmp1[91]+z[83]*tmp1[92]+z[47]*tmp1[93]+z[11]*tmp1[94]+z[72]*tmp1[95]+z[36]*tmp1[96];
            tab[i+62*stg_first]=tab[i+62*stg_first]+z[0]*tmp1[0]
            +z[62]*tmp1[1]+z[27]*tmp1[2]+z[89]*tmp1[3]+z[54]*tmp1[4]+z[19]*tmp1[5]+z[81]*tmp1[6]+z[46]*tmp1[7]+z[11]*tmp1[8]+z[73]*tmp1[9]+z[38]*tmp1[10]
            +z[3]*tmp1[11]+z[65]*tmp1[12]+z[30]*tmp1[13]+z[92]*tmp1[14]+z[57]*tmp1[15]+z[22]*tmp1[16]+z[84]*tmp1[17]+z[49]*tmp1[18]+z[14]*tmp1[19]+z[76]*tmp1[20]
            +z[41]*tmp1[21]+z[6]*tmp1[22]+z[68]*tmp1[23]+z[33]*tmp1[24]+z[95]*tmp1[25]+z[60]*tmp1[26]+z[25]*tmp1[27]+z[87]*tmp1[28]+z[52]*tmp1[29]+z[17]*tmp1[30]
            +z[79]*tmp1[31]+z[44]*tmp1[32]+z[9]*tmp1[33]+z[71]*tmp1[34]+z[36]*tmp1[35]+z[1]*tmp1[36]+z[63]*tmp1[37]+z[28]*tmp1[38]+z[90]*tmp1[39]+z[55]*tmp1[40]
            +z[20]*tmp1[41]+z[82]*tmp1[42]+z[47]*tmp1[43]+z[12]*tmp1[44]+z[74]*tmp1[45]+z[39]*tmp1[46]+z[4]*tmp1[47]+z[66]*tmp1[48]+z[31]*tmp1[49]+z[93]*tmp1[50]
            +z[58]*tmp1[51]+z[23]*tmp1[52]+z[85]*tmp1[53]+z[50]*tmp1[54]+z[15]*tmp1[55]+z[77]*tmp1[56]+z[42]*tmp1[57]+z[7]*tmp1[58]+z[69]*tmp1[59]+z[34]*tmp1[60]
            +z[96]*tmp1[61]+z[61]*tmp1[62]+z[26]*tmp1[63]+z[88]*tmp1[64]+z[53]*tmp1[65]+z[18]*tmp1[66]+z[80]*tmp1[67]+z[45]*tmp1[68]+z[10]*tmp1[69]+z[72]*tmp1[70]
            +z[37]*tmp1[71]+z[2]*tmp1[72]+z[64]*tmp1[73]+z[29]*tmp1[74]+z[91]*tmp1[75]+z[56]*tmp1[76]+z[21]*tmp1[77]+z[83]*tmp1[78]+z[48]*tmp1[79]+z[13]*tmp1[80]
            +z[75]*tmp1[81]+z[40]*tmp1[82]+z[5]*tmp1[83]+z[67]*tmp1[84]+z[32]*tmp1[85]+z[94]*tmp1[86]+z[59]*tmp1[87]+z[24]*tmp1[88]+z[86]*tmp1[89]+z[51]*tmp1[90]
            +z[16]*tmp1[91]+z[78]*tmp1[92]+z[43]*tmp1[93]+z[8]*tmp1[94]+z[70]*tmp1[95]+z[35]*tmp1[96];
            tab[i+63*stg_first]=tab[i+63*stg_first]+z[0]*tmp1[0]
            +z[63]*tmp1[1]+z[29]*tmp1[2]+z[92]*tmp1[3]+z[58]*tmp1[4]+z[24]*tmp1[5]+z[87]*tmp1[6]+z[53]*tmp1[7]+z[19]*tmp1[8]+z[82]*tmp1[9]+z[48]*tmp1[10]
            +z[14]*tmp1[11]+z[77]*tmp1[12]+z[43]*tmp1[13]+z[9]*tmp1[14]+z[72]*tmp1[15]+z[38]*tmp1[16]+z[4]*tmp1[17]+z[67]*tmp1[18]+z[33]*tmp1[19]+z[96]*tmp1[20]
            +z[62]*tmp1[21]+z[28]*tmp1[22]+z[91]*tmp1[23]+z[57]*tmp1[24]+z[23]*tmp1[25]+z[86]*tmp1[26]+z[52]*tmp1[27]+z[18]*tmp1[28]+z[81]*tmp1[29]+z[47]*tmp1[30]
            +z[13]*tmp1[31]+z[76]*tmp1[32]+z[42]*tmp1[33]+z[8]*tmp1[34]+z[71]*tmp1[35]+z[37]*tmp1[36]+z[3]*tmp1[37]+z[66]*tmp1[38]+z[32]*tmp1[39]+z[95]*tmp1[40]
            +z[61]*tmp1[41]+z[27]*tmp1[42]+z[90]*tmp1[43]+z[56]*tmp1[44]+z[22]*tmp1[45]+z[85]*tmp1[46]+z[51]*tmp1[47]+z[17]*tmp1[48]+z[80]*tmp1[49]+z[46]*tmp1[50]
            +z[12]*tmp1[51]+z[75]*tmp1[52]+z[41]*tmp1[53]+z[7]*tmp1[54]+z[70]*tmp1[55]+z[36]*tmp1[56]+z[2]*tmp1[57]+z[65]*tmp1[58]+z[31]*tmp1[59]+z[94]*tmp1[60]
            +z[60]*tmp1[61]+z[26]*tmp1[62]+z[89]*tmp1[63]+z[55]*tmp1[64]+z[21]*tmp1[65]+z[84]*tmp1[66]+z[50]*tmp1[67]+z[16]*tmp1[68]+z[79]*tmp1[69]+z[45]*tmp1[70]
            +z[11]*tmp1[71]+z[74]*tmp1[72]+z[40]*tmp1[73]+z[6]*tmp1[74]+z[69]*tmp1[75]+z[35]*tmp1[76]+z[1]*tmp1[77]+z[64]*tmp1[78]+z[30]*tmp1[79]+z[93]*tmp1[80]
            +z[59]*tmp1[81]+z[25]*tmp1[82]+z[88]*tmp1[83]+z[54]*tmp1[84]+z[20]*tmp1[85]+z[83]*tmp1[86]+z[49]*tmp1[87]+z[15]*tmp1[88]+z[78]*tmp1[89]+z[44]*tmp1[90]
            +z[10]*tmp1[91]+z[73]*tmp1[92]+z[39]*tmp1[93]+z[5]*tmp1[94]+z[68]*tmp1[95]+z[34]*tmp1[96];
            tab[i+64*stg_first]=tab[i+64*stg_first]+z[0]*tmp1[0]
            +z[64]*tmp1[1]+z[31]*tmp1[2]+z[95]*tmp1[3]+z[62]*tmp1[4]+z[29]*tmp1[5]+z[93]*tmp1[6]+z[60]*tmp1[7]+z[27]*tmp1[8]+z[91]*tmp1[9]+z[58]*tmp1[10]
            +z[25]*tmp1[11]+z[89]*tmp1[12]+z[56]*tmp1[13]+z[23]*tmp1[14]+z[87]*tmp1[15]+z[54]*tmp1[16]+z[21]*tmp1[17]+z[85]*tmp1[18]+z[52]*tmp1[19]+z[19]*tmp1[20]
            +z[83]*tmp1[21]+z[50]*tmp1[22]+z[17]*tmp1[23]+z[81]*tmp1[24]+z[48]*tmp1[25]+z[15]*tmp1[26]+z[79]*tmp1[27]+z[46]*tmp1[28]+z[13]*tmp1[29]+z[77]*tmp1[30]
            +z[44]*tmp1[31]+z[11]*tmp1[32]+z[75]*tmp1[33]+z[42]*tmp1[34]+z[9]*tmp1[35]+z[73]*tmp1[36]+z[40]*tmp1[37]+z[7]*tmp1[38]+z[71]*tmp1[39]+z[38]*tmp1[40]
            +z[5]*tmp1[41]+z[69]*tmp1[42]+z[36]*tmp1[43]+z[3]*tmp1[44]+z[67]*tmp1[45]+z[34]*tmp1[46]+z[1]*tmp1[47]+z[65]*tmp1[48]+z[32]*tmp1[49]+z[96]*tmp1[50]
            +z[63]*tmp1[51]+z[30]*tmp1[52]+z[94]*tmp1[53]+z[61]*tmp1[54]+z[28]*tmp1[55]+z[92]*tmp1[56]+z[59]*tmp1[57]+z[26]*tmp1[58]+z[90]*tmp1[59]+z[57]*tmp1[60]
            +z[24]*tmp1[61]+z[88]*tmp1[62]+z[55]*tmp1[63]+z[22]*tmp1[64]+z[86]*tmp1[65]+z[53]*tmp1[66]+z[20]*tmp1[67]+z[84]*tmp1[68]+z[51]*tmp1[69]+z[18]*tmp1[70]
            +z[82]*tmp1[71]+z[49]*tmp1[72]+z[16]*tmp1[73]+z[80]*tmp1[74]+z[47]*tmp1[75]+z[14]*tmp1[76]+z[78]*tmp1[77]+z[45]*tmp1[78]+z[12]*tmp1[79]+z[76]*tmp1[80]
            +z[43]*tmp1[81]+z[10]*tmp1[82]+z[74]*tmp1[83]+z[41]*tmp1[84]+z[8]*tmp1[85]+z[72]*tmp1[86]+z[39]*tmp1[87]+z[6]*tmp1[88]+z[70]*tmp1[89]+z[37]*tmp1[90]
            +z[4]*tmp1[91]+z[68]*tmp1[92]+z[35]*tmp1[93]+z[2]*tmp1[94]+z[66]*tmp1[95]+z[33]*tmp1[96];
            tab[i+65*stg_first]=tab[i+65*stg_first]+z[0]*tmp1[0]
            +z[65]*tmp1[1]+z[33]*tmp1[2]+z[1]*tmp1[3]+z[66]*tmp1[4]+z[34]*tmp1[5]+z[2]*tmp1[6]+z[67]*tmp1[7]+z[35]*tmp1[8]+z[3]*tmp1[9]+z[68]*tmp1[10]
            +z[36]*tmp1[11]+z[4]*tmp1[12]+z[69]*tmp1[13]+z[37]*tmp1[14]+z[5]*tmp1[15]+z[70]*tmp1[16]+z[38]*tmp1[17]+z[6]*tmp1[18]+z[71]*tmp1[19]+z[39]*tmp1[20]
            +z[7]*tmp1[21]+z[72]*tmp1[22]+z[40]*tmp1[23]+z[8]*tmp1[24]+z[73]*tmp1[25]+z[41]*tmp1[26]+z[9]*tmp1[27]+z[74]*tmp1[28]+z[42]*tmp1[29]+z[10]*tmp1[30]
            +z[75]*tmp1[31]+z[43]*tmp1[32]+z[11]*tmp1[33]+z[76]*tmp1[34]+z[44]*tmp1[35]+z[12]*tmp1[36]+z[77]*tmp1[37]+z[45]*tmp1[38]+z[13]*tmp1[39]+z[78]*tmp1[40]
            +z[46]*tmp1[41]+z[14]*tmp1[42]+z[79]*tmp1[43]+z[47]*tmp1[44]+z[15]*tmp1[45]+z[80]*tmp1[46]+z[48]*tmp1[47]+z[16]*tmp1[48]+z[81]*tmp1[49]+z[49]*tmp1[50]
            +z[17]*tmp1[51]+z[82]*tmp1[52]+z[50]*tmp1[53]+z[18]*tmp1[54]+z[83]*tmp1[55]+z[51]*tmp1[56]+z[19]*tmp1[57]+z[84]*tmp1[58]+z[52]*tmp1[59]+z[20]*tmp1[60]
            +z[85]*tmp1[61]+z[53]*tmp1[62]+z[21]*tmp1[63]+z[86]*tmp1[64]+z[54]*tmp1[65]+z[22]*tmp1[66]+z[87]*tmp1[67]+z[55]*tmp1[68]+z[23]*tmp1[69]+z[88]*tmp1[70]
            +z[56]*tmp1[71]+z[24]*tmp1[72]+z[89]*tmp1[73]+z[57]*tmp1[74]+z[25]*tmp1[75]+z[90]*tmp1[76]+z[58]*tmp1[77]+z[26]*tmp1[78]+z[91]*tmp1[79]+z[59]*tmp1[80]
            +z[27]*tmp1[81]+z[92]*tmp1[82]+z[60]*tmp1[83]+z[28]*tmp1[84]+z[93]*tmp1[85]+z[61]*tmp1[86]+z[29]*tmp1[87]+z[94]*tmp1[88]+z[62]*tmp1[89]+z[30]*tmp1[90]
            +z[95]*tmp1[91]+z[63]*tmp1[92]+z[31]*tmp1[93]+z[96]*tmp1[94]+z[64]*tmp1[95]+z[32]*tmp1[96];
            tab[i+66*stg_first]=tab[i+66*stg_first]+z[0]*tmp1[0]
            +z[66]*tmp1[1]+z[35]*tmp1[2]+z[4]*tmp1[3]+z[70]*tmp1[4]+z[39]*tmp1[5]+z[8]*tmp1[6]+z[74]*tmp1[7]+z[43]*tmp1[8]+z[12]*tmp1[9]+z[78]*tmp1[10]
            +z[47]*tmp1[11]+z[16]*tmp1[12]+z[82]*tmp1[13]+z[51]*tmp1[14]+z[20]*tmp1[15]+z[86]*tmp1[16]+z[55]*tmp1[17]+z[24]*tmp1[18]+z[90]*tmp1[19]+z[59]*tmp1[20]
            +z[28]*tmp1[21]+z[94]*tmp1[22]+z[63]*tmp1[23]+z[32]*tmp1[24]+z[1]*tmp1[25]+z[67]*tmp1[26]+z[36]*tmp1[27]+z[5]*tmp1[28]+z[71]*tmp1[29]+z[40]*tmp1[30]
            +z[9]*tmp1[31]+z[75]*tmp1[32]+z[44]*tmp1[33]+z[13]*tmp1[34]+z[79]*tmp1[35]+z[48]*tmp1[36]+z[17]*tmp1[37]+z[83]*tmp1[38]+z[52]*tmp1[39]+z[21]*tmp1[40]
            +z[87]*tmp1[41]+z[56]*tmp1[42]+z[25]*tmp1[43]+z[91]*tmp1[44]+z[60]*tmp1[45]+z[29]*tmp1[46]+z[95]*tmp1[47]+z[64]*tmp1[48]+z[33]*tmp1[49]+z[2]*tmp1[50]
            +z[68]*tmp1[51]+z[37]*tmp1[52]+z[6]*tmp1[53]+z[72]*tmp1[54]+z[41]*tmp1[55]+z[10]*tmp1[56]+z[76]*tmp1[57]+z[45]*tmp1[58]+z[14]*tmp1[59]+z[80]*tmp1[60]
            +z[49]*tmp1[61]+z[18]*tmp1[62]+z[84]*tmp1[63]+z[53]*tmp1[64]+z[22]*tmp1[65]+z[88]*tmp1[66]+z[57]*tmp1[67]+z[26]*tmp1[68]+z[92]*tmp1[69]+z[61]*tmp1[70]
            +z[30]*tmp1[71]+z[96]*tmp1[72]+z[65]*tmp1[73]+z[34]*tmp1[74]+z[3]*tmp1[75]+z[69]*tmp1[76]+z[38]*tmp1[77]+z[7]*tmp1[78]+z[73]*tmp1[79]+z[42]*tmp1[80]
            +z[11]*tmp1[81]+z[77]*tmp1[82]+z[46]*tmp1[83]+z[15]*tmp1[84]+z[81]*tmp1[85]+z[50]*tmp1[86]+z[19]*tmp1[87]+z[85]*tmp1[88]+z[54]*tmp1[89]+z[23]*tmp1[90]
            +z[89]*tmp1[91]+z[58]*tmp1[92]+z[27]*tmp1[93]+z[93]*tmp1[94]+z[62]*tmp1[95]+z[31]*tmp1[96];
            tab[i+67*stg_first]=tab[i+67*stg_first]+z[0]*tmp1[0]
            +z[67]*tmp1[1]+z[37]*tmp1[2]+z[7]*tmp1[3]+z[74]*tmp1[4]+z[44]*tmp1[5]+z[14]*tmp1[6]+z[81]*tmp1[7]+z[51]*tmp1[8]+z[21]*tmp1[9]+z[88]*tmp1[10]
            +z[58]*tmp1[11]+z[28]*tmp1[12]+z[95]*tmp1[13]+z[65]*tmp1[14]+z[35]*tmp1[15]+z[5]*tmp1[16]+z[72]*tmp1[17]+z[42]*tmp1[18]+z[12]*tmp1[19]+z[79]*tmp1[20]
            +z[49]*tmp1[21]+z[19]*tmp1[22]+z[86]*tmp1[23]+z[56]*tmp1[24]+z[26]*tmp1[25]+z[93]*tmp1[26]+z[63]*tmp1[27]+z[33]*tmp1[28]+z[3]*tmp1[29]+z[70]*tmp1[30]
            +z[40]*tmp1[31]+z[10]*tmp1[32]+z[77]*tmp1[33]+z[47]*tmp1[34]+z[17]*tmp1[35]+z[84]*tmp1[36]+z[54]*tmp1[37]+z[24]*tmp1[38]+z[91]*tmp1[39]+z[61]*tmp1[40]
            +z[31]*tmp1[41]+z[1]*tmp1[42]+z[68]*tmp1[43]+z[38]*tmp1[44]+z[8]*tmp1[45]+z[75]*tmp1[46]+z[45]*tmp1[47]+z[15]*tmp1[48]+z[82]*tmp1[49]+z[52]*tmp1[50]
            +z[22]*tmp1[51]+z[89]*tmp1[52]+z[59]*tmp1[53]+z[29]*tmp1[54]+z[96]*tmp1[55]+z[66]*tmp1[56]+z[36]*tmp1[57]+z[6]*tmp1[58]+z[73]*tmp1[59]+z[43]*tmp1[60]
            +z[13]*tmp1[61]+z[80]*tmp1[62]+z[50]*tmp1[63]+z[20]*tmp1[64]+z[87]*tmp1[65]+z[57]*tmp1[66]+z[27]*tmp1[67]+z[94]*tmp1[68]+z[64]*tmp1[69]+z[34]*tmp1[70]
            +z[4]*tmp1[71]+z[71]*tmp1[72]+z[41]*tmp1[73]+z[11]*tmp1[74]+z[78]*tmp1[75]+z[48]*tmp1[76]+z[18]*tmp1[77]+z[85]*tmp1[78]+z[55]*tmp1[79]+z[25]*tmp1[80]
            +z[92]*tmp1[81]+z[62]*tmp1[82]+z[32]*tmp1[83]+z[2]*tmp1[84]+z[69]*tmp1[85]+z[39]*tmp1[86]+z[9]*tmp1[87]+z[76]*tmp1[88]+z[46]*tmp1[89]+z[16]*tmp1[90]
            +z[83]*tmp1[91]+z[53]*tmp1[92]+z[23]*tmp1[93]+z[90]*tmp1[94]+z[60]*tmp1[95]+z[30]*tmp1[96];
            tab[i+68*stg_first]=tab[i+68*stg_first]+z[0]*tmp1[0]
            +z[68]*tmp1[1]+z[39]*tmp1[2]+z[10]*tmp1[3]+z[78]*tmp1[4]+z[49]*tmp1[5]+z[20]*tmp1[6]+z[88]*tmp1[7]+z[59]*tmp1[8]+z[30]*tmp1[9]+z[1]*tmp1[10]
            +z[69]*tmp1[11]+z[40]*tmp1[12]+z[11]*tmp1[13]+z[79]*tmp1[14]+z[50]*tmp1[15]+z[21]*tmp1[16]+z[89]*tmp1[17]+z[60]*tmp1[18]+z[31]*tmp1[19]+z[2]*tmp1[20]
            +z[70]*tmp1[21]+z[41]*tmp1[22]+z[12]*tmp1[23]+z[80]*tmp1[24]+z[51]*tmp1[25]+z[22]*tmp1[26]+z[90]*tmp1[27]+z[61]*tmp1[28]+z[32]*tmp1[29]+z[3]*tmp1[30]
            +z[71]*tmp1[31]+z[42]*tmp1[32]+z[13]*tmp1[33]+z[81]*tmp1[34]+z[52]*tmp1[35]+z[23]*tmp1[36]+z[91]*tmp1[37]+z[62]*tmp1[38]+z[33]*tmp1[39]+z[4]*tmp1[40]
            +z[72]*tmp1[41]+z[43]*tmp1[42]+z[14]*tmp1[43]+z[82]*tmp1[44]+z[53]*tmp1[45]+z[24]*tmp1[46]+z[92]*tmp1[47]+z[63]*tmp1[48]+z[34]*tmp1[49]+z[5]*tmp1[50]
            +z[73]*tmp1[51]+z[44]*tmp1[52]+z[15]*tmp1[53]+z[83]*tmp1[54]+z[54]*tmp1[55]+z[25]*tmp1[56]+z[93]*tmp1[57]+z[64]*tmp1[58]+z[35]*tmp1[59]+z[6]*tmp1[60]
            +z[74]*tmp1[61]+z[45]*tmp1[62]+z[16]*tmp1[63]+z[84]*tmp1[64]+z[55]*tmp1[65]+z[26]*tmp1[66]+z[94]*tmp1[67]+z[65]*tmp1[68]+z[36]*tmp1[69]+z[7]*tmp1[70]
            +z[75]*tmp1[71]+z[46]*tmp1[72]+z[17]*tmp1[73]+z[85]*tmp1[74]+z[56]*tmp1[75]+z[27]*tmp1[76]+z[95]*tmp1[77]+z[66]*tmp1[78]+z[37]*tmp1[79]+z[8]*tmp1[80]
            +z[76]*tmp1[81]+z[47]*tmp1[82]+z[18]*tmp1[83]+z[86]*tmp1[84]+z[57]*tmp1[85]+z[28]*tmp1[86]+z[96]*tmp1[87]+z[67]*tmp1[88]+z[38]*tmp1[89]+z[9]*tmp1[90]
            +z[77]*tmp1[91]+z[48]*tmp1[92]+z[19]*tmp1[93]+z[87]*tmp1[94]+z[58]*tmp1[95]+z[29]*tmp1[96];
            tab[i+69*stg_first]=tab[i+69*stg_first]+z[0]*tmp1[0]
            +z[69]*tmp1[1]+z[41]*tmp1[2]+z[13]*tmp1[3]+z[82]*tmp1[4]+z[54]*tmp1[5]+z[26]*tmp1[6]+z[95]*tmp1[7]+z[67]*tmp1[8]+z[39]*tmp1[9]+z[11]*tmp1[10]
            +z[80]*tmp1[11]+z[52]*tmp1[12]+z[24]*tmp1[13]+z[93]*tmp1[14]+z[65]*tmp1[15]+z[37]*tmp1[16]+z[9]*tmp1[17]+z[78]*tmp1[18]+z[50]*tmp1[19]+z[22]*tmp1[20]
            +z[91]*tmp1[21]+z[63]*tmp1[22]+z[35]*tmp1[23]+z[7]*tmp1[24]+z[76]*tmp1[25]+z[48]*tmp1[26]+z[20]*tmp1[27]+z[89]*tmp1[28]+z[61]*tmp1[29]+z[33]*tmp1[30]
            +z[5]*tmp1[31]+z[74]*tmp1[32]+z[46]*tmp1[33]+z[18]*tmp1[34]+z[87]*tmp1[35]+z[59]*tmp1[36]+z[31]*tmp1[37]+z[3]*tmp1[38]+z[72]*tmp1[39]+z[44]*tmp1[40]
            +z[16]*tmp1[41]+z[85]*tmp1[42]+z[57]*tmp1[43]+z[29]*tmp1[44]+z[1]*tmp1[45]+z[70]*tmp1[46]+z[42]*tmp1[47]+z[14]*tmp1[48]+z[83]*tmp1[49]+z[55]*tmp1[50]
            +z[27]*tmp1[51]+z[96]*tmp1[52]+z[68]*tmp1[53]+z[40]*tmp1[54]+z[12]*tmp1[55]+z[81]*tmp1[56]+z[53]*tmp1[57]+z[25]*tmp1[58]+z[94]*tmp1[59]+z[66]*tmp1[60]
            +z[38]*tmp1[61]+z[10]*tmp1[62]+z[79]*tmp1[63]+z[51]*tmp1[64]+z[23]*tmp1[65]+z[92]*tmp1[66]+z[64]*tmp1[67]+z[36]*tmp1[68]+z[8]*tmp1[69]+z[77]*tmp1[70]
            +z[49]*tmp1[71]+z[21]*tmp1[72]+z[90]*tmp1[73]+z[62]*tmp1[74]+z[34]*tmp1[75]+z[6]*tmp1[76]+z[75]*tmp1[77]+z[47]*tmp1[78]+z[19]*tmp1[79]+z[88]*tmp1[80]
            +z[60]*tmp1[81]+z[32]*tmp1[82]+z[4]*tmp1[83]+z[73]*tmp1[84]+z[45]*tmp1[85]+z[17]*tmp1[86]+z[86]*tmp1[87]+z[58]*tmp1[88]+z[30]*tmp1[89]+z[2]*tmp1[90]
            +z[71]*tmp1[91]+z[43]*tmp1[92]+z[15]*tmp1[93]+z[84]*tmp1[94]+z[56]*tmp1[95]+z[28]*tmp1[96];
            tab[i+70*stg_first]=tab[i+70*stg_first]+z[0]*tmp1[0]
            +z[70]*tmp1[1]+z[43]*tmp1[2]+z[16]*tmp1[3]+z[86]*tmp1[4]+z[59]*tmp1[5]+z[32]*tmp1[6]+z[5]*tmp1[7]+z[75]*tmp1[8]+z[48]*tmp1[9]+z[21]*tmp1[10]
            +z[91]*tmp1[11]+z[64]*tmp1[12]+z[37]*tmp1[13]+z[10]*tmp1[14]+z[80]*tmp1[15]+z[53]*tmp1[16]+z[26]*tmp1[17]+z[96]*tmp1[18]+z[69]*tmp1[19]+z[42]*tmp1[20]
            +z[15]*tmp1[21]+z[85]*tmp1[22]+z[58]*tmp1[23]+z[31]*tmp1[24]+z[4]*tmp1[25]+z[74]*tmp1[26]+z[47]*tmp1[27]+z[20]*tmp1[28]+z[90]*tmp1[29]+z[63]*tmp1[30]
            +z[36]*tmp1[31]+z[9]*tmp1[32]+z[79]*tmp1[33]+z[52]*tmp1[34]+z[25]*tmp1[35]+z[95]*tmp1[36]+z[68]*tmp1[37]+z[41]*tmp1[38]+z[14]*tmp1[39]+z[84]*tmp1[40]
            +z[57]*tmp1[41]+z[30]*tmp1[42]+z[3]*tmp1[43]+z[73]*tmp1[44]+z[46]*tmp1[45]+z[19]*tmp1[46]+z[89]*tmp1[47]+z[62]*tmp1[48]+z[35]*tmp1[49]+z[8]*tmp1[50]
            +z[78]*tmp1[51]+z[51]*tmp1[52]+z[24]*tmp1[53]+z[94]*tmp1[54]+z[67]*tmp1[55]+z[40]*tmp1[56]+z[13]*tmp1[57]+z[83]*tmp1[58]+z[56]*tmp1[59]+z[29]*tmp1[60]
            +z[2]*tmp1[61]+z[72]*tmp1[62]+z[45]*tmp1[63]+z[18]*tmp1[64]+z[88]*tmp1[65]+z[61]*tmp1[66]+z[34]*tmp1[67]+z[7]*tmp1[68]+z[77]*tmp1[69]+z[50]*tmp1[70]
            +z[23]*tmp1[71]+z[93]*tmp1[72]+z[66]*tmp1[73]+z[39]*tmp1[74]+z[12]*tmp1[75]+z[82]*tmp1[76]+z[55]*tmp1[77]+z[28]*tmp1[78]+z[1]*tmp1[79]+z[71]*tmp1[80]
            +z[44]*tmp1[81]+z[17]*tmp1[82]+z[87]*tmp1[83]+z[60]*tmp1[84]+z[33]*tmp1[85]+z[6]*tmp1[86]+z[76]*tmp1[87]+z[49]*tmp1[88]+z[22]*tmp1[89]+z[92]*tmp1[90]
            +z[65]*tmp1[91]+z[38]*tmp1[92]+z[11]*tmp1[93]+z[81]*tmp1[94]+z[54]*tmp1[95]+z[27]*tmp1[96];
            tab[i+71*stg_first]=tab[i+71*stg_first]+z[0]*tmp1[0]
            +z[71]*tmp1[1]+z[45]*tmp1[2]+z[19]*tmp1[3]+z[90]*tmp1[4]+z[64]*tmp1[5]+z[38]*tmp1[6]+z[12]*tmp1[7]+z[83]*tmp1[8]+z[57]*tmp1[9]+z[31]*tmp1[10]
            +z[5]*tmp1[11]+z[76]*tmp1[12]+z[50]*tmp1[13]+z[24]*tmp1[14]+z[95]*tmp1[15]+z[69]*tmp1[16]+z[43]*tmp1[17]+z[17]*tmp1[18]+z[88]*tmp1[19]+z[62]*tmp1[20]
            +z[36]*tmp1[21]+z[10]*tmp1[22]+z[81]*tmp1[23]+z[55]*tmp1[24]+z[29]*tmp1[25]+z[3]*tmp1[26]+z[74]*tmp1[27]+z[48]*tmp1[28]+z[22]*tmp1[29]+z[93]*tmp1[30]
            +z[67]*tmp1[31]+z[41]*tmp1[32]+z[15]*tmp1[33]+z[86]*tmp1[34]+z[60]*tmp1[35]+z[34]*tmp1[36]+z[8]*tmp1[37]+z[79]*tmp1[38]+z[53]*tmp1[39]+z[27]*tmp1[40]
            +z[1]*tmp1[41]+z[72]*tmp1[42]+z[46]*tmp1[43]+z[20]*tmp1[44]+z[91]*tmp1[45]+z[65]*tmp1[46]+z[39]*tmp1[47]+z[13]*tmp1[48]+z[84]*tmp1[49]+z[58]*tmp1[50]
            +z[32]*tmp1[51]+z[6]*tmp1[52]+z[77]*tmp1[53]+z[51]*tmp1[54]+z[25]*tmp1[55]+z[96]*tmp1[56]+z[70]*tmp1[57]+z[44]*tmp1[58]+z[18]*tmp1[59]+z[89]*tmp1[60]
            +z[63]*tmp1[61]+z[37]*tmp1[62]+z[11]*tmp1[63]+z[82]*tmp1[64]+z[56]*tmp1[65]+z[30]*tmp1[66]+z[4]*tmp1[67]+z[75]*tmp1[68]+z[49]*tmp1[69]+z[23]*tmp1[70]
            +z[94]*tmp1[71]+z[68]*tmp1[72]+z[42]*tmp1[73]+z[16]*tmp1[74]+z[87]*tmp1[75]+z[61]*tmp1[76]+z[35]*tmp1[77]+z[9]*tmp1[78]+z[80]*tmp1[79]+z[54]*tmp1[80]
            +z[28]*tmp1[81]+z[2]*tmp1[82]+z[73]*tmp1[83]+z[47]*tmp1[84]+z[21]*tmp1[85]+z[92]*tmp1[86]+z[66]*tmp1[87]+z[40]*tmp1[88]+z[14]*tmp1[89]+z[85]*tmp1[90]
            +z[59]*tmp1[91]+z[33]*tmp1[92]+z[7]*tmp1[93]+z[78]*tmp1[94]+z[52]*tmp1[95]+z[26]*tmp1[96];
            tab[i+72*stg_first]=tab[i+72*stg_first]+z[0]*tmp1[0]
            +z[72]*tmp1[1]+z[47]*tmp1[2]+z[22]*tmp1[3]+z[94]*tmp1[4]+z[69]*tmp1[5]+z[44]*tmp1[6]+z[19]*tmp1[7]+z[91]*tmp1[8]+z[66]*tmp1[9]+z[41]*tmp1[10]
            +z[16]*tmp1[11]+z[88]*tmp1[12]+z[63]*tmp1[13]+z[38]*tmp1[14]+z[13]*tmp1[15]+z[85]*tmp1[16]+z[60]*tmp1[17]+z[35]*tmp1[18]+z[10]*tmp1[19]+z[82]*tmp1[20]
            +z[57]*tmp1[21]+z[32]*tmp1[22]+z[7]*tmp1[23]+z[79]*tmp1[24]+z[54]*tmp1[25]+z[29]*tmp1[26]+z[4]*tmp1[27]+z[76]*tmp1[28]+z[51]*tmp1[29]+z[26]*tmp1[30]
            +z[1]*tmp1[31]+z[73]*tmp1[32]+z[48]*tmp1[33]+z[23]*tmp1[34]+z[95]*tmp1[35]+z[70]*tmp1[36]+z[45]*tmp1[37]+z[20]*tmp1[38]+z[92]*tmp1[39]+z[67]*tmp1[40]
            +z[42]*tmp1[41]+z[17]*tmp1[42]+z[89]*tmp1[43]+z[64]*tmp1[44]+z[39]*tmp1[45]+z[14]*tmp1[46]+z[86]*tmp1[47]+z[61]*tmp1[48]+z[36]*tmp1[49]+z[11]*tmp1[50]
            +z[83]*tmp1[51]+z[58]*tmp1[52]+z[33]*tmp1[53]+z[8]*tmp1[54]+z[80]*tmp1[55]+z[55]*tmp1[56]+z[30]*tmp1[57]+z[5]*tmp1[58]+z[77]*tmp1[59]+z[52]*tmp1[60]
            +z[27]*tmp1[61]+z[2]*tmp1[62]+z[74]*tmp1[63]+z[49]*tmp1[64]+z[24]*tmp1[65]+z[96]*tmp1[66]+z[71]*tmp1[67]+z[46]*tmp1[68]+z[21]*tmp1[69]+z[93]*tmp1[70]
            +z[68]*tmp1[71]+z[43]*tmp1[72]+z[18]*tmp1[73]+z[90]*tmp1[74]+z[65]*tmp1[75]+z[40]*tmp1[76]+z[15]*tmp1[77]+z[87]*tmp1[78]+z[62]*tmp1[79]+z[37]*tmp1[80]
            +z[12]*tmp1[81]+z[84]*tmp1[82]+z[59]*tmp1[83]+z[34]*tmp1[84]+z[9]*tmp1[85]+z[81]*tmp1[86]+z[56]*tmp1[87]+z[31]*tmp1[88]+z[6]*tmp1[89]+z[78]*tmp1[90]
            +z[53]*tmp1[91]+z[28]*tmp1[92]+z[3]*tmp1[93]+z[75]*tmp1[94]+z[50]*tmp1[95]+z[25]*tmp1[96];
            tab[i+73*stg_first]=tab[i+73*stg_first]+z[0]*tmp1[0]
            +z[73]*tmp1[1]+z[49]*tmp1[2]+z[25]*tmp1[3]+z[1]*tmp1[4]+z[74]*tmp1[5]+z[50]*tmp1[6]+z[26]*tmp1[7]+z[2]*tmp1[8]+z[75]*tmp1[9]+z[51]*tmp1[10]
            +z[27]*tmp1[11]+z[3]*tmp1[12]+z[76]*tmp1[13]+z[52]*tmp1[14]+z[28]*tmp1[15]+z[4]*tmp1[16]+z[77]*tmp1[17]+z[53]*tmp1[18]+z[29]*tmp1[19]+z[5]*tmp1[20]
            +z[78]*tmp1[21]+z[54]*tmp1[22]+z[30]*tmp1[23]+z[6]*tmp1[24]+z[79]*tmp1[25]+z[55]*tmp1[26]+z[31]*tmp1[27]+z[7]*tmp1[28]+z[80]*tmp1[29]+z[56]*tmp1[30]
            +z[32]*tmp1[31]+z[8]*tmp1[32]+z[81]*tmp1[33]+z[57]*tmp1[34]+z[33]*tmp1[35]+z[9]*tmp1[36]+z[82]*tmp1[37]+z[58]*tmp1[38]+z[34]*tmp1[39]+z[10]*tmp1[40]
            +z[83]*tmp1[41]+z[59]*tmp1[42]+z[35]*tmp1[43]+z[11]*tmp1[44]+z[84]*tmp1[45]+z[60]*tmp1[46]+z[36]*tmp1[47]+z[12]*tmp1[48]+z[85]*tmp1[49]+z[61]*tmp1[50]
            +z[37]*tmp1[51]+z[13]*tmp1[52]+z[86]*tmp1[53]+z[62]*tmp1[54]+z[38]*tmp1[55]+z[14]*tmp1[56]+z[87]*tmp1[57]+z[63]*tmp1[58]+z[39]*tmp1[59]+z[15]*tmp1[60]
            +z[88]*tmp1[61]+z[64]*tmp1[62]+z[40]*tmp1[63]+z[16]*tmp1[64]+z[89]*tmp1[65]+z[65]*tmp1[66]+z[41]*tmp1[67]+z[17]*tmp1[68]+z[90]*tmp1[69]+z[66]*tmp1[70]
            +z[42]*tmp1[71]+z[18]*tmp1[72]+z[91]*tmp1[73]+z[67]*tmp1[74]+z[43]*tmp1[75]+z[19]*tmp1[76]+z[92]*tmp1[77]+z[68]*tmp1[78]+z[44]*tmp1[79]+z[20]*tmp1[80]
            +z[93]*tmp1[81]+z[69]*tmp1[82]+z[45]*tmp1[83]+z[21]*tmp1[84]+z[94]*tmp1[85]+z[70]*tmp1[86]+z[46]*tmp1[87]+z[22]*tmp1[88]+z[95]*tmp1[89]+z[71]*tmp1[90]
            +z[47]*tmp1[91]+z[23]*tmp1[92]+z[96]*tmp1[93]+z[72]*tmp1[94]+z[48]*tmp1[95]+z[24]*tmp1[96];
            tab[i+74*stg_first]=tab[i+74*stg_first]+z[0]*tmp1[0]
            +z[74]*tmp1[1]+z[51]*tmp1[2]+z[28]*tmp1[3]+z[5]*tmp1[4]+z[79]*tmp1[5]+z[56]*tmp1[6]+z[33]*tmp1[7]+z[10]*tmp1[8]+z[84]*tmp1[9]+z[61]*tmp1[10]
            +z[38]*tmp1[11]+z[15]*tmp1[12]+z[89]*tmp1[13]+z[66]*tmp1[14]+z[43]*tmp1[15]+z[20]*tmp1[16]+z[94]*tmp1[17]+z[71]*tmp1[18]+z[48]*tmp1[19]+z[25]*tmp1[20]
            +z[2]*tmp1[21]+z[76]*tmp1[22]+z[53]*tmp1[23]+z[30]*tmp1[24]+z[7]*tmp1[25]+z[81]*tmp1[26]+z[58]*tmp1[27]+z[35]*tmp1[28]+z[12]*tmp1[29]+z[86]*tmp1[30]
            +z[63]*tmp1[31]+z[40]*tmp1[32]+z[17]*tmp1[33]+z[91]*tmp1[34]+z[68]*tmp1[35]+z[45]*tmp1[36]+z[22]*tmp1[37]+z[96]*tmp1[38]+z[73]*tmp1[39]+z[50]*tmp1[40]
            +z[27]*tmp1[41]+z[4]*tmp1[42]+z[78]*tmp1[43]+z[55]*tmp1[44]+z[32]*tmp1[45]+z[9]*tmp1[46]+z[83]*tmp1[47]+z[60]*tmp1[48]+z[37]*tmp1[49]+z[14]*tmp1[50]
            +z[88]*tmp1[51]+z[65]*tmp1[52]+z[42]*tmp1[53]+z[19]*tmp1[54]+z[93]*tmp1[55]+z[70]*tmp1[56]+z[47]*tmp1[57]+z[24]*tmp1[58]+z[1]*tmp1[59]+z[75]*tmp1[60]
            +z[52]*tmp1[61]+z[29]*tmp1[62]+z[6]*tmp1[63]+z[80]*tmp1[64]+z[57]*tmp1[65]+z[34]*tmp1[66]+z[11]*tmp1[67]+z[85]*tmp1[68]+z[62]*tmp1[69]+z[39]*tmp1[70]
            +z[16]*tmp1[71]+z[90]*tmp1[72]+z[67]*tmp1[73]+z[44]*tmp1[74]+z[21]*tmp1[75]+z[95]*tmp1[76]+z[72]*tmp1[77]+z[49]*tmp1[78]+z[26]*tmp1[79]+z[3]*tmp1[80]
            +z[77]*tmp1[81]+z[54]*tmp1[82]+z[31]*tmp1[83]+z[8]*tmp1[84]+z[82]*tmp1[85]+z[59]*tmp1[86]+z[36]*tmp1[87]+z[13]*tmp1[88]+z[87]*tmp1[89]+z[64]*tmp1[90]
            +z[41]*tmp1[91]+z[18]*tmp1[92]+z[92]*tmp1[93]+z[69]*tmp1[94]+z[46]*tmp1[95]+z[23]*tmp1[96];
            tab[i+75*stg_first]=tab[i+75*stg_first]+z[0]*tmp1[0]
            +z[75]*tmp1[1]+z[53]*tmp1[2]+z[31]*tmp1[3]+z[9]*tmp1[4]+z[84]*tmp1[5]+z[62]*tmp1[6]+z[40]*tmp1[7]+z[18]*tmp1[8]+z[93]*tmp1[9]+z[71]*tmp1[10]
            +z[49]*tmp1[11]+z[27]*tmp1[12]+z[5]*tmp1[13]+z[80]*tmp1[14]+z[58]*tmp1[15]+z[36]*tmp1[16]+z[14]*tmp1[17]+z[89]*tmp1[18]+z[67]*tmp1[19]+z[45]*tmp1[20]
            +z[23]*tmp1[21]+z[1]*tmp1[22]+z[76]*tmp1[23]+z[54]*tmp1[24]+z[32]*tmp1[25]+z[10]*tmp1[26]+z[85]*tmp1[27]+z[63]*tmp1[28]+z[41]*tmp1[29]+z[19]*tmp1[30]
            +z[94]*tmp1[31]+z[72]*tmp1[32]+z[50]*tmp1[33]+z[28]*tmp1[34]+z[6]*tmp1[35]+z[81]*tmp1[36]+z[59]*tmp1[37]+z[37]*tmp1[38]+z[15]*tmp1[39]+z[90]*tmp1[40]
            +z[68]*tmp1[41]+z[46]*tmp1[42]+z[24]*tmp1[43]+z[2]*tmp1[44]+z[77]*tmp1[45]+z[55]*tmp1[46]+z[33]*tmp1[47]+z[11]*tmp1[48]+z[86]*tmp1[49]+z[64]*tmp1[50]
            +z[42]*tmp1[51]+z[20]*tmp1[52]+z[95]*tmp1[53]+z[73]*tmp1[54]+z[51]*tmp1[55]+z[29]*tmp1[56]+z[7]*tmp1[57]+z[82]*tmp1[58]+z[60]*tmp1[59]+z[38]*tmp1[60]
            +z[16]*tmp1[61]+z[91]*tmp1[62]+z[69]*tmp1[63]+z[47]*tmp1[64]+z[25]*tmp1[65]+z[3]*tmp1[66]+z[78]*tmp1[67]+z[56]*tmp1[68]+z[34]*tmp1[69]+z[12]*tmp1[70]
            +z[87]*tmp1[71]+z[65]*tmp1[72]+z[43]*tmp1[73]+z[21]*tmp1[74]+z[96]*tmp1[75]+z[74]*tmp1[76]+z[52]*tmp1[77]+z[30]*tmp1[78]+z[8]*tmp1[79]+z[83]*tmp1[80]
            +z[61]*tmp1[81]+z[39]*tmp1[82]+z[17]*tmp1[83]+z[92]*tmp1[84]+z[70]*tmp1[85]+z[48]*tmp1[86]+z[26]*tmp1[87]+z[4]*tmp1[88]+z[79]*tmp1[89]+z[57]*tmp1[90]
            +z[35]*tmp1[91]+z[13]*tmp1[92]+z[88]*tmp1[93]+z[66]*tmp1[94]+z[44]*tmp1[95]+z[22]*tmp1[96];
            tab[i+76*stg_first]=tab[i+76*stg_first]+z[0]*tmp1[0]
            +z[76]*tmp1[1]+z[55]*tmp1[2]+z[34]*tmp1[3]+z[13]*tmp1[4]+z[89]*tmp1[5]+z[68]*tmp1[6]+z[47]*tmp1[7]+z[26]*tmp1[8]+z[5]*tmp1[9]+z[81]*tmp1[10]
            +z[60]*tmp1[11]+z[39]*tmp1[12]+z[18]*tmp1[13]+z[94]*tmp1[14]+z[73]*tmp1[15]+z[52]*tmp1[16]+z[31]*tmp1[17]+z[10]*tmp1[18]+z[86]*tmp1[19]+z[65]*tmp1[20]
            +z[44]*tmp1[21]+z[23]*tmp1[22]+z[2]*tmp1[23]+z[78]*tmp1[24]+z[57]*tmp1[25]+z[36]*tmp1[26]+z[15]*tmp1[27]+z[91]*tmp1[28]+z[70]*tmp1[29]+z[49]*tmp1[30]
            +z[28]*tmp1[31]+z[7]*tmp1[32]+z[83]*tmp1[33]+z[62]*tmp1[34]+z[41]*tmp1[35]+z[20]*tmp1[36]+z[96]*tmp1[37]+z[75]*tmp1[38]+z[54]*tmp1[39]+z[33]*tmp1[40]
            +z[12]*tmp1[41]+z[88]*tmp1[42]+z[67]*tmp1[43]+z[46]*tmp1[44]+z[25]*tmp1[45]+z[4]*tmp1[46]+z[80]*tmp1[47]+z[59]*tmp1[48]+z[38]*tmp1[49]+z[17]*tmp1[50]
            +z[93]*tmp1[51]+z[72]*tmp1[52]+z[51]*tmp1[53]+z[30]*tmp1[54]+z[9]*tmp1[55]+z[85]*tmp1[56]+z[64]*tmp1[57]+z[43]*tmp1[58]+z[22]*tmp1[59]+z[1]*tmp1[60]
            +z[77]*tmp1[61]+z[56]*tmp1[62]+z[35]*tmp1[63]+z[14]*tmp1[64]+z[90]*tmp1[65]+z[69]*tmp1[66]+z[48]*tmp1[67]+z[27]*tmp1[68]+z[6]*tmp1[69]+z[82]*tmp1[70]
            +z[61]*tmp1[71]+z[40]*tmp1[72]+z[19]*tmp1[73]+z[95]*tmp1[74]+z[74]*tmp1[75]+z[53]*tmp1[76]+z[32]*tmp1[77]+z[11]*tmp1[78]+z[87]*tmp1[79]+z[66]*tmp1[80]
            +z[45]*tmp1[81]+z[24]*tmp1[82]+z[3]*tmp1[83]+z[79]*tmp1[84]+z[58]*tmp1[85]+z[37]*tmp1[86]+z[16]*tmp1[87]+z[92]*tmp1[88]+z[71]*tmp1[89]+z[50]*tmp1[90]
            +z[29]*tmp1[91]+z[8]*tmp1[92]+z[84]*tmp1[93]+z[63]*tmp1[94]+z[42]*tmp1[95]+z[21]*tmp1[96];
            tab[i+77*stg_first]=tab[i+77*stg_first]+z[0]*tmp1[0]
            +z[77]*tmp1[1]+z[57]*tmp1[2]+z[37]*tmp1[3]+z[17]*tmp1[4]+z[94]*tmp1[5]+z[74]*tmp1[6]+z[54]*tmp1[7]+z[34]*tmp1[8]+z[14]*tmp1[9]+z[91]*tmp1[10]
            +z[71]*tmp1[11]+z[51]*tmp1[12]+z[31]*tmp1[13]+z[11]*tmp1[14]+z[88]*tmp1[15]+z[68]*tmp1[16]+z[48]*tmp1[17]+z[28]*tmp1[18]+z[8]*tmp1[19]+z[85]*tmp1[20]
            +z[65]*tmp1[21]+z[45]*tmp1[22]+z[25]*tmp1[23]+z[5]*tmp1[24]+z[82]*tmp1[25]+z[62]*tmp1[26]+z[42]*tmp1[27]+z[22]*tmp1[28]+z[2]*tmp1[29]+z[79]*tmp1[30]
            +z[59]*tmp1[31]+z[39]*tmp1[32]+z[19]*tmp1[33]+z[96]*tmp1[34]+z[76]*tmp1[35]+z[56]*tmp1[36]+z[36]*tmp1[37]+z[16]*tmp1[38]+z[93]*tmp1[39]+z[73]*tmp1[40]
            +z[53]*tmp1[41]+z[33]*tmp1[42]+z[13]*tmp1[43]+z[90]*tmp1[44]+z[70]*tmp1[45]+z[50]*tmp1[46]+z[30]*tmp1[47]+z[10]*tmp1[48]+z[87]*tmp1[49]+z[67]*tmp1[50]
            +z[47]*tmp1[51]+z[27]*tmp1[52]+z[7]*tmp1[53]+z[84]*tmp1[54]+z[64]*tmp1[55]+z[44]*tmp1[56]+z[24]*tmp1[57]+z[4]*tmp1[58]+z[81]*tmp1[59]+z[61]*tmp1[60]
            +z[41]*tmp1[61]+z[21]*tmp1[62]+z[1]*tmp1[63]+z[78]*tmp1[64]+z[58]*tmp1[65]+z[38]*tmp1[66]+z[18]*tmp1[67]+z[95]*tmp1[68]+z[75]*tmp1[69]+z[55]*tmp1[70]
            +z[35]*tmp1[71]+z[15]*tmp1[72]+z[92]*tmp1[73]+z[72]*tmp1[74]+z[52]*tmp1[75]+z[32]*tmp1[76]+z[12]*tmp1[77]+z[89]*tmp1[78]+z[69]*tmp1[79]+z[49]*tmp1[80]
            +z[29]*tmp1[81]+z[9]*tmp1[82]+z[86]*tmp1[83]+z[66]*tmp1[84]+z[46]*tmp1[85]+z[26]*tmp1[86]+z[6]*tmp1[87]+z[83]*tmp1[88]+z[63]*tmp1[89]+z[43]*tmp1[90]
            +z[23]*tmp1[91]+z[3]*tmp1[92]+z[80]*tmp1[93]+z[60]*tmp1[94]+z[40]*tmp1[95]+z[20]*tmp1[96];
            tab[i+78*stg_first]=tab[i+78*stg_first]+z[0]*tmp1[0]
            +z[78]*tmp1[1]+z[59]*tmp1[2]+z[40]*tmp1[3]+z[21]*tmp1[4]+z[2]*tmp1[5]+z[80]*tmp1[6]+z[61]*tmp1[7]+z[42]*tmp1[8]+z[23]*tmp1[9]+z[4]*tmp1[10]
            +z[82]*tmp1[11]+z[63]*tmp1[12]+z[44]*tmp1[13]+z[25]*tmp1[14]+z[6]*tmp1[15]+z[84]*tmp1[16]+z[65]*tmp1[17]+z[46]*tmp1[18]+z[27]*tmp1[19]+z[8]*tmp1[20]
            +z[86]*tmp1[21]+z[67]*tmp1[22]+z[48]*tmp1[23]+z[29]*tmp1[24]+z[10]*tmp1[25]+z[88]*tmp1[26]+z[69]*tmp1[27]+z[50]*tmp1[28]+z[31]*tmp1[29]+z[12]*tmp1[30]
            +z[90]*tmp1[31]+z[71]*tmp1[32]+z[52]*tmp1[33]+z[33]*tmp1[34]+z[14]*tmp1[35]+z[92]*tmp1[36]+z[73]*tmp1[37]+z[54]*tmp1[38]+z[35]*tmp1[39]+z[16]*tmp1[40]
            +z[94]*tmp1[41]+z[75]*tmp1[42]+z[56]*tmp1[43]+z[37]*tmp1[44]+z[18]*tmp1[45]+z[96]*tmp1[46]+z[77]*tmp1[47]+z[58]*tmp1[48]+z[39]*tmp1[49]+z[20]*tmp1[50]
            +z[1]*tmp1[51]+z[79]*tmp1[52]+z[60]*tmp1[53]+z[41]*tmp1[54]+z[22]*tmp1[55]+z[3]*tmp1[56]+z[81]*tmp1[57]+z[62]*tmp1[58]+z[43]*tmp1[59]+z[24]*tmp1[60]
            +z[5]*tmp1[61]+z[83]*tmp1[62]+z[64]*tmp1[63]+z[45]*tmp1[64]+z[26]*tmp1[65]+z[7]*tmp1[66]+z[85]*tmp1[67]+z[66]*tmp1[68]+z[47]*tmp1[69]+z[28]*tmp1[70]
            +z[9]*tmp1[71]+z[87]*tmp1[72]+z[68]*tmp1[73]+z[49]*tmp1[74]+z[30]*tmp1[75]+z[11]*tmp1[76]+z[89]*tmp1[77]+z[70]*tmp1[78]+z[51]*tmp1[79]+z[32]*tmp1[80]
            +z[13]*tmp1[81]+z[91]*tmp1[82]+z[72]*tmp1[83]+z[53]*tmp1[84]+z[34]*tmp1[85]+z[15]*tmp1[86]+z[93]*tmp1[87]+z[74]*tmp1[88]+z[55]*tmp1[89]+z[36]*tmp1[90]
            +z[17]*tmp1[91]+z[95]*tmp1[92]+z[76]*tmp1[93]+z[57]*tmp1[94]+z[38]*tmp1[95]+z[19]*tmp1[96];
            tab[i+79*stg_first]=tab[i+79*stg_first]+z[0]*tmp1[0]
            +z[79]*tmp1[1]+z[61]*tmp1[2]+z[43]*tmp1[3]+z[25]*tmp1[4]+z[7]*tmp1[5]+z[86]*tmp1[6]+z[68]*tmp1[7]+z[50]*tmp1[8]+z[32]*tmp1[9]+z[14]*tmp1[10]
            +z[93]*tmp1[11]+z[75]*tmp1[12]+z[57]*tmp1[13]+z[39]*tmp1[14]+z[21]*tmp1[15]+z[3]*tmp1[16]+z[82]*tmp1[17]+z[64]*tmp1[18]+z[46]*tmp1[19]+z[28]*tmp1[20]
            +z[10]*tmp1[21]+z[89]*tmp1[22]+z[71]*tmp1[23]+z[53]*tmp1[24]+z[35]*tmp1[25]+z[17]*tmp1[26]+z[96]*tmp1[27]+z[78]*tmp1[28]+z[60]*tmp1[29]+z[42]*tmp1[30]
            +z[24]*tmp1[31]+z[6]*tmp1[32]+z[85]*tmp1[33]+z[67]*tmp1[34]+z[49]*tmp1[35]+z[31]*tmp1[36]+z[13]*tmp1[37]+z[92]*tmp1[38]+z[74]*tmp1[39]+z[56]*tmp1[40]
            +z[38]*tmp1[41]+z[20]*tmp1[42]+z[2]*tmp1[43]+z[81]*tmp1[44]+z[63]*tmp1[45]+z[45]*tmp1[46]+z[27]*tmp1[47]+z[9]*tmp1[48]+z[88]*tmp1[49]+z[70]*tmp1[50]
            +z[52]*tmp1[51]+z[34]*tmp1[52]+z[16]*tmp1[53]+z[95]*tmp1[54]+z[77]*tmp1[55]+z[59]*tmp1[56]+z[41]*tmp1[57]+z[23]*tmp1[58]+z[5]*tmp1[59]+z[84]*tmp1[60]
            +z[66]*tmp1[61]+z[48]*tmp1[62]+z[30]*tmp1[63]+z[12]*tmp1[64]+z[91]*tmp1[65]+z[73]*tmp1[66]+z[55]*tmp1[67]+z[37]*tmp1[68]+z[19]*tmp1[69]+z[1]*tmp1[70]
            +z[80]*tmp1[71]+z[62]*tmp1[72]+z[44]*tmp1[73]+z[26]*tmp1[74]+z[8]*tmp1[75]+z[87]*tmp1[76]+z[69]*tmp1[77]+z[51]*tmp1[78]+z[33]*tmp1[79]+z[15]*tmp1[80]
            +z[94]*tmp1[81]+z[76]*tmp1[82]+z[58]*tmp1[83]+z[40]*tmp1[84]+z[22]*tmp1[85]+z[4]*tmp1[86]+z[83]*tmp1[87]+z[65]*tmp1[88]+z[47]*tmp1[89]+z[29]*tmp1[90]
            +z[11]*tmp1[91]+z[90]*tmp1[92]+z[72]*tmp1[93]+z[54]*tmp1[94]+z[36]*tmp1[95]+z[18]*tmp1[96];
            tab[i+80*stg_first]=tab[i+80*stg_first]+z[0]*tmp1[0]
            +z[80]*tmp1[1]+z[63]*tmp1[2]+z[46]*tmp1[3]+z[29]*tmp1[4]+z[12]*tmp1[5]+z[92]*tmp1[6]+z[75]*tmp1[7]+z[58]*tmp1[8]+z[41]*tmp1[9]+z[24]*tmp1[10]
            +z[7]*tmp1[11]+z[87]*tmp1[12]+z[70]*tmp1[13]+z[53]*tmp1[14]+z[36]*tmp1[15]+z[19]*tmp1[16]+z[2]*tmp1[17]+z[82]*tmp1[18]+z[65]*tmp1[19]+z[48]*tmp1[20]
            +z[31]*tmp1[21]+z[14]*tmp1[22]+z[94]*tmp1[23]+z[77]*tmp1[24]+z[60]*tmp1[25]+z[43]*tmp1[26]+z[26]*tmp1[27]+z[9]*tmp1[28]+z[89]*tmp1[29]+z[72]*tmp1[30]
            +z[55]*tmp1[31]+z[38]*tmp1[32]+z[21]*tmp1[33]+z[4]*tmp1[34]+z[84]*tmp1[35]+z[67]*tmp1[36]+z[50]*tmp1[37]+z[33]*tmp1[38]+z[16]*tmp1[39]+z[96]*tmp1[40]
            +z[79]*tmp1[41]+z[62]*tmp1[42]+z[45]*tmp1[43]+z[28]*tmp1[44]+z[11]*tmp1[45]+z[91]*tmp1[46]+z[74]*tmp1[47]+z[57]*tmp1[48]+z[40]*tmp1[49]+z[23]*tmp1[50]
            +z[6]*tmp1[51]+z[86]*tmp1[52]+z[69]*tmp1[53]+z[52]*tmp1[54]+z[35]*tmp1[55]+z[18]*tmp1[56]+z[1]*tmp1[57]+z[81]*tmp1[58]+z[64]*tmp1[59]+z[47]*tmp1[60]
            +z[30]*tmp1[61]+z[13]*tmp1[62]+z[93]*tmp1[63]+z[76]*tmp1[64]+z[59]*tmp1[65]+z[42]*tmp1[66]+z[25]*tmp1[67]+z[8]*tmp1[68]+z[88]*tmp1[69]+z[71]*tmp1[70]
            +z[54]*tmp1[71]+z[37]*tmp1[72]+z[20]*tmp1[73]+z[3]*tmp1[74]+z[83]*tmp1[75]+z[66]*tmp1[76]+z[49]*tmp1[77]+z[32]*tmp1[78]+z[15]*tmp1[79]+z[95]*tmp1[80]
            +z[78]*tmp1[81]+z[61]*tmp1[82]+z[44]*tmp1[83]+z[27]*tmp1[84]+z[10]*tmp1[85]+z[90]*tmp1[86]+z[73]*tmp1[87]+z[56]*tmp1[88]+z[39]*tmp1[89]+z[22]*tmp1[90]
            +z[5]*tmp1[91]+z[85]*tmp1[92]+z[68]*tmp1[93]+z[51]*tmp1[94]+z[34]*tmp1[95]+z[17]*tmp1[96];
            tab[i+81*stg_first]=tab[i+81*stg_first]+z[0]*tmp1[0]
            +z[81]*tmp1[1]+z[65]*tmp1[2]+z[49]*tmp1[3]+z[33]*tmp1[4]+z[17]*tmp1[5]+z[1]*tmp1[6]+z[82]*tmp1[7]+z[66]*tmp1[8]+z[50]*tmp1[9]+z[34]*tmp1[10]
            +z[18]*tmp1[11]+z[2]*tmp1[12]+z[83]*tmp1[13]+z[67]*tmp1[14]+z[51]*tmp1[15]+z[35]*tmp1[16]+z[19]*tmp1[17]+z[3]*tmp1[18]+z[84]*tmp1[19]+z[68]*tmp1[20]
            +z[52]*tmp1[21]+z[36]*tmp1[22]+z[20]*tmp1[23]+z[4]*tmp1[24]+z[85]*tmp1[25]+z[69]*tmp1[26]+z[53]*tmp1[27]+z[37]*tmp1[28]+z[21]*tmp1[29]+z[5]*tmp1[30]
            +z[86]*tmp1[31]+z[70]*tmp1[32]+z[54]*tmp1[33]+z[38]*tmp1[34]+z[22]*tmp1[35]+z[6]*tmp1[36]+z[87]*tmp1[37]+z[71]*tmp1[38]+z[55]*tmp1[39]+z[39]*tmp1[40]
            +z[23]*tmp1[41]+z[7]*tmp1[42]+z[88]*tmp1[43]+z[72]*tmp1[44]+z[56]*tmp1[45]+z[40]*tmp1[46]+z[24]*tmp1[47]+z[8]*tmp1[48]+z[89]*tmp1[49]+z[73]*tmp1[50]
            +z[57]*tmp1[51]+z[41]*tmp1[52]+z[25]*tmp1[53]+z[9]*tmp1[54]+z[90]*tmp1[55]+z[74]*tmp1[56]+z[58]*tmp1[57]+z[42]*tmp1[58]+z[26]*tmp1[59]+z[10]*tmp1[60]
            +z[91]*tmp1[61]+z[75]*tmp1[62]+z[59]*tmp1[63]+z[43]*tmp1[64]+z[27]*tmp1[65]+z[11]*tmp1[66]+z[92]*tmp1[67]+z[76]*tmp1[68]+z[60]*tmp1[69]+z[44]*tmp1[70]
            +z[28]*tmp1[71]+z[12]*tmp1[72]+z[93]*tmp1[73]+z[77]*tmp1[74]+z[61]*tmp1[75]+z[45]*tmp1[76]+z[29]*tmp1[77]+z[13]*tmp1[78]+z[94]*tmp1[79]+z[78]*tmp1[80]
            +z[62]*tmp1[81]+z[46]*tmp1[82]+z[30]*tmp1[83]+z[14]*tmp1[84]+z[95]*tmp1[85]+z[79]*tmp1[86]+z[63]*tmp1[87]+z[47]*tmp1[88]+z[31]*tmp1[89]+z[15]*tmp1[90]
            +z[96]*tmp1[91]+z[80]*tmp1[92]+z[64]*tmp1[93]+z[48]*tmp1[94]+z[32]*tmp1[95]+z[16]*tmp1[96];
            tab[i+82*stg_first]=tab[i+82*stg_first]+z[0]*tmp1[0]
            +z[82]*tmp1[1]+z[67]*tmp1[2]+z[52]*tmp1[3]+z[37]*tmp1[4]+z[22]*tmp1[5]+z[7]*tmp1[6]+z[89]*tmp1[7]+z[74]*tmp1[8]+z[59]*tmp1[9]+z[44]*tmp1[10]
            +z[29]*tmp1[11]+z[14]*tmp1[12]+z[96]*tmp1[13]+z[81]*tmp1[14]+z[66]*tmp1[15]+z[51]*tmp1[16]+z[36]*tmp1[17]+z[21]*tmp1[18]+z[6]*tmp1[19]+z[88]*tmp1[20]
            +z[73]*tmp1[21]+z[58]*tmp1[22]+z[43]*tmp1[23]+z[28]*tmp1[24]+z[13]*tmp1[25]+z[95]*tmp1[26]+z[80]*tmp1[27]+z[65]*tmp1[28]+z[50]*tmp1[29]+z[35]*tmp1[30]
            +z[20]*tmp1[31]+z[5]*tmp1[32]+z[87]*tmp1[33]+z[72]*tmp1[34]+z[57]*tmp1[35]+z[42]*tmp1[36]+z[27]*tmp1[37]+z[12]*tmp1[38]+z[94]*tmp1[39]+z[79]*tmp1[40]
            +z[64]*tmp1[41]+z[49]*tmp1[42]+z[34]*tmp1[43]+z[19]*tmp1[44]+z[4]*tmp1[45]+z[86]*tmp1[46]+z[71]*tmp1[47]+z[56]*tmp1[48]+z[41]*tmp1[49]+z[26]*tmp1[50]
            +z[11]*tmp1[51]+z[93]*tmp1[52]+z[78]*tmp1[53]+z[63]*tmp1[54]+z[48]*tmp1[55]+z[33]*tmp1[56]+z[18]*tmp1[57]+z[3]*tmp1[58]+z[85]*tmp1[59]+z[70]*tmp1[60]
            +z[55]*tmp1[61]+z[40]*tmp1[62]+z[25]*tmp1[63]+z[10]*tmp1[64]+z[92]*tmp1[65]+z[77]*tmp1[66]+z[62]*tmp1[67]+z[47]*tmp1[68]+z[32]*tmp1[69]+z[17]*tmp1[70]
            +z[2]*tmp1[71]+z[84]*tmp1[72]+z[69]*tmp1[73]+z[54]*tmp1[74]+z[39]*tmp1[75]+z[24]*tmp1[76]+z[9]*tmp1[77]+z[91]*tmp1[78]+z[76]*tmp1[79]+z[61]*tmp1[80]
            +z[46]*tmp1[81]+z[31]*tmp1[82]+z[16]*tmp1[83]+z[1]*tmp1[84]+z[83]*tmp1[85]+z[68]*tmp1[86]+z[53]*tmp1[87]+z[38]*tmp1[88]+z[23]*tmp1[89]+z[8]*tmp1[90]
            +z[90]*tmp1[91]+z[75]*tmp1[92]+z[60]*tmp1[93]+z[45]*tmp1[94]+z[30]*tmp1[95]+z[15]*tmp1[96];
            tab[i+83*stg_first]=tab[i+83*stg_first]+z[0]*tmp1[0]
            +z[83]*tmp1[1]+z[69]*tmp1[2]+z[55]*tmp1[3]+z[41]*tmp1[4]+z[27]*tmp1[5]+z[13]*tmp1[6]+z[96]*tmp1[7]+z[82]*tmp1[8]+z[68]*tmp1[9]+z[54]*tmp1[10]
            +z[40]*tmp1[11]+z[26]*tmp1[12]+z[12]*tmp1[13]+z[95]*tmp1[14]+z[81]*tmp1[15]+z[67]*tmp1[16]+z[53]*tmp1[17]+z[39]*tmp1[18]+z[25]*tmp1[19]+z[11]*tmp1[20]
            +z[94]*tmp1[21]+z[80]*tmp1[22]+z[66]*tmp1[23]+z[52]*tmp1[24]+z[38]*tmp1[25]+z[24]*tmp1[26]+z[10]*tmp1[27]+z[93]*tmp1[28]+z[79]*tmp1[29]+z[65]*tmp1[30]
            +z[51]*tmp1[31]+z[37]*tmp1[32]+z[23]*tmp1[33]+z[9]*tmp1[34]+z[92]*tmp1[35]+z[78]*tmp1[36]+z[64]*tmp1[37]+z[50]*tmp1[38]+z[36]*tmp1[39]+z[22]*tmp1[40]
            +z[8]*tmp1[41]+z[91]*tmp1[42]+z[77]*tmp1[43]+z[63]*tmp1[44]+z[49]*tmp1[45]+z[35]*tmp1[46]+z[21]*tmp1[47]+z[7]*tmp1[48]+z[90]*tmp1[49]+z[76]*tmp1[50]
            +z[62]*tmp1[51]+z[48]*tmp1[52]+z[34]*tmp1[53]+z[20]*tmp1[54]+z[6]*tmp1[55]+z[89]*tmp1[56]+z[75]*tmp1[57]+z[61]*tmp1[58]+z[47]*tmp1[59]+z[33]*tmp1[60]
            +z[19]*tmp1[61]+z[5]*tmp1[62]+z[88]*tmp1[63]+z[74]*tmp1[64]+z[60]*tmp1[65]+z[46]*tmp1[66]+z[32]*tmp1[67]+z[18]*tmp1[68]+z[4]*tmp1[69]+z[87]*tmp1[70]
            +z[73]*tmp1[71]+z[59]*tmp1[72]+z[45]*tmp1[73]+z[31]*tmp1[74]+z[17]*tmp1[75]+z[3]*tmp1[76]+z[86]*tmp1[77]+z[72]*tmp1[78]+z[58]*tmp1[79]+z[44]*tmp1[80]
            +z[30]*tmp1[81]+z[16]*tmp1[82]+z[2]*tmp1[83]+z[85]*tmp1[84]+z[71]*tmp1[85]+z[57]*tmp1[86]+z[43]*tmp1[87]+z[29]*tmp1[88]+z[15]*tmp1[89]+z[1]*tmp1[90]
            +z[84]*tmp1[91]+z[70]*tmp1[92]+z[56]*tmp1[93]+z[42]*tmp1[94]+z[28]*tmp1[95]+z[14]*tmp1[96];
            tab[i+84*stg_first]=tab[i+84*stg_first]+z[0]*tmp1[0]
            +z[84]*tmp1[1]+z[71]*tmp1[2]+z[58]*tmp1[3]+z[45]*tmp1[4]+z[32]*tmp1[5]+z[19]*tmp1[6]+z[6]*tmp1[7]+z[90]*tmp1[8]+z[77]*tmp1[9]+z[64]*tmp1[10]
            +z[51]*tmp1[11]+z[38]*tmp1[12]+z[25]*tmp1[13]+z[12]*tmp1[14]+z[96]*tmp1[15]+z[83]*tmp1[16]+z[70]*tmp1[17]+z[57]*tmp1[18]+z[44]*tmp1[19]+z[31]*tmp1[20]
            +z[18]*tmp1[21]+z[5]*tmp1[22]+z[89]*tmp1[23]+z[76]*tmp1[24]+z[63]*tmp1[25]+z[50]*tmp1[26]+z[37]*tmp1[27]+z[24]*tmp1[28]+z[11]*tmp1[29]+z[95]*tmp1[30]
            +z[82]*tmp1[31]+z[69]*tmp1[32]+z[56]*tmp1[33]+z[43]*tmp1[34]+z[30]*tmp1[35]+z[17]*tmp1[36]+z[4]*tmp1[37]+z[88]*tmp1[38]+z[75]*tmp1[39]+z[62]*tmp1[40]
            +z[49]*tmp1[41]+z[36]*tmp1[42]+z[23]*tmp1[43]+z[10]*tmp1[44]+z[94]*tmp1[45]+z[81]*tmp1[46]+z[68]*tmp1[47]+z[55]*tmp1[48]+z[42]*tmp1[49]+z[29]*tmp1[50]
            +z[16]*tmp1[51]+z[3]*tmp1[52]+z[87]*tmp1[53]+z[74]*tmp1[54]+z[61]*tmp1[55]+z[48]*tmp1[56]+z[35]*tmp1[57]+z[22]*tmp1[58]+z[9]*tmp1[59]+z[93]*tmp1[60]
            +z[80]*tmp1[61]+z[67]*tmp1[62]+z[54]*tmp1[63]+z[41]*tmp1[64]+z[28]*tmp1[65]+z[15]*tmp1[66]+z[2]*tmp1[67]+z[86]*tmp1[68]+z[73]*tmp1[69]+z[60]*tmp1[70]
            +z[47]*tmp1[71]+z[34]*tmp1[72]+z[21]*tmp1[73]+z[8]*tmp1[74]+z[92]*tmp1[75]+z[79]*tmp1[76]+z[66]*tmp1[77]+z[53]*tmp1[78]+z[40]*tmp1[79]+z[27]*tmp1[80]
            +z[14]*tmp1[81]+z[1]*tmp1[82]+z[85]*tmp1[83]+z[72]*tmp1[84]+z[59]*tmp1[85]+z[46]*tmp1[86]+z[33]*tmp1[87]+z[20]*tmp1[88]+z[7]*tmp1[89]+z[91]*tmp1[90]
            +z[78]*tmp1[91]+z[65]*tmp1[92]+z[52]*tmp1[93]+z[39]*tmp1[94]+z[26]*tmp1[95]+z[13]*tmp1[96];
            tab[i+85*stg_first]=tab[i+85*stg_first]+z[0]*tmp1[0]
            +z[85]*tmp1[1]+z[73]*tmp1[2]+z[61]*tmp1[3]+z[49]*tmp1[4]+z[37]*tmp1[5]+z[25]*tmp1[6]+z[13]*tmp1[7]+z[1]*tmp1[8]+z[86]*tmp1[9]+z[74]*tmp1[10]
            +z[62]*tmp1[11]+z[50]*tmp1[12]+z[38]*tmp1[13]+z[26]*tmp1[14]+z[14]*tmp1[15]+z[2]*tmp1[16]+z[87]*tmp1[17]+z[75]*tmp1[18]+z[63]*tmp1[19]+z[51]*tmp1[20]
            +z[39]*tmp1[21]+z[27]*tmp1[22]+z[15]*tmp1[23]+z[3]*tmp1[24]+z[88]*tmp1[25]+z[76]*tmp1[26]+z[64]*tmp1[27]+z[52]*tmp1[28]+z[40]*tmp1[29]+z[28]*tmp1[30]
            +z[16]*tmp1[31]+z[4]*tmp1[32]+z[89]*tmp1[33]+z[77]*tmp1[34]+z[65]*tmp1[35]+z[53]*tmp1[36]+z[41]*tmp1[37]+z[29]*tmp1[38]+z[17]*tmp1[39]+z[5]*tmp1[40]
            +z[90]*tmp1[41]+z[78]*tmp1[42]+z[66]*tmp1[43]+z[54]*tmp1[44]+z[42]*tmp1[45]+z[30]*tmp1[46]+z[18]*tmp1[47]+z[6]*tmp1[48]+z[91]*tmp1[49]+z[79]*tmp1[50]
            +z[67]*tmp1[51]+z[55]*tmp1[52]+z[43]*tmp1[53]+z[31]*tmp1[54]+z[19]*tmp1[55]+z[7]*tmp1[56]+z[92]*tmp1[57]+z[80]*tmp1[58]+z[68]*tmp1[59]+z[56]*tmp1[60]
            +z[44]*tmp1[61]+z[32]*tmp1[62]+z[20]*tmp1[63]+z[8]*tmp1[64]+z[93]*tmp1[65]+z[81]*tmp1[66]+z[69]*tmp1[67]+z[57]*tmp1[68]+z[45]*tmp1[69]+z[33]*tmp1[70]
            +z[21]*tmp1[71]+z[9]*tmp1[72]+z[94]*tmp1[73]+z[82]*tmp1[74]+z[70]*tmp1[75]+z[58]*tmp1[76]+z[46]*tmp1[77]+z[34]*tmp1[78]+z[22]*tmp1[79]+z[10]*tmp1[80]
            +z[95]*tmp1[81]+z[83]*tmp1[82]+z[71]*tmp1[83]+z[59]*tmp1[84]+z[47]*tmp1[85]+z[35]*tmp1[86]+z[23]*tmp1[87]+z[11]*tmp1[88]+z[96]*tmp1[89]+z[84]*tmp1[90]
            +z[72]*tmp1[91]+z[60]*tmp1[92]+z[48]*tmp1[93]+z[36]*tmp1[94]+z[24]*tmp1[95]+z[12]*tmp1[96];
            tab[i+86*stg_first]=tab[i+86*stg_first]+z[0]*tmp1[0]
            +z[86]*tmp1[1]+z[75]*tmp1[2]+z[64]*tmp1[3]+z[53]*tmp1[4]+z[42]*tmp1[5]+z[31]*tmp1[6]+z[20]*tmp1[7]+z[9]*tmp1[8]+z[95]*tmp1[9]+z[84]*tmp1[10]
            +z[73]*tmp1[11]+z[62]*tmp1[12]+z[51]*tmp1[13]+z[40]*tmp1[14]+z[29]*tmp1[15]+z[18]*tmp1[16]+z[7]*tmp1[17]+z[93]*tmp1[18]+z[82]*tmp1[19]+z[71]*tmp1[20]
            +z[60]*tmp1[21]+z[49]*tmp1[22]+z[38]*tmp1[23]+z[27]*tmp1[24]+z[16]*tmp1[25]+z[5]*tmp1[26]+z[91]*tmp1[27]+z[80]*tmp1[28]+z[69]*tmp1[29]+z[58]*tmp1[30]
            +z[47]*tmp1[31]+z[36]*tmp1[32]+z[25]*tmp1[33]+z[14]*tmp1[34]+z[3]*tmp1[35]+z[89]*tmp1[36]+z[78]*tmp1[37]+z[67]*tmp1[38]+z[56]*tmp1[39]+z[45]*tmp1[40]
            +z[34]*tmp1[41]+z[23]*tmp1[42]+z[12]*tmp1[43]+z[1]*tmp1[44]+z[87]*tmp1[45]+z[76]*tmp1[46]+z[65]*tmp1[47]+z[54]*tmp1[48]+z[43]*tmp1[49]+z[32]*tmp1[50]
            +z[21]*tmp1[51]+z[10]*tmp1[52]+z[96]*tmp1[53]+z[85]*tmp1[54]+z[74]*tmp1[55]+z[63]*tmp1[56]+z[52]*tmp1[57]+z[41]*tmp1[58]+z[30]*tmp1[59]+z[19]*tmp1[60]
            +z[8]*tmp1[61]+z[94]*tmp1[62]+z[83]*tmp1[63]+z[72]*tmp1[64]+z[61]*tmp1[65]+z[50]*tmp1[66]+z[39]*tmp1[67]+z[28]*tmp1[68]+z[17]*tmp1[69]+z[6]*tmp1[70]
            +z[92]*tmp1[71]+z[81]*tmp1[72]+z[70]*tmp1[73]+z[59]*tmp1[74]+z[48]*tmp1[75]+z[37]*tmp1[76]+z[26]*tmp1[77]+z[15]*tmp1[78]+z[4]*tmp1[79]+z[90]*tmp1[80]
            +z[79]*tmp1[81]+z[68]*tmp1[82]+z[57]*tmp1[83]+z[46]*tmp1[84]+z[35]*tmp1[85]+z[24]*tmp1[86]+z[13]*tmp1[87]+z[2]*tmp1[88]+z[88]*tmp1[89]+z[77]*tmp1[90]
            +z[66]*tmp1[91]+z[55]*tmp1[92]+z[44]*tmp1[93]+z[33]*tmp1[94]+z[22]*tmp1[95]+z[11]*tmp1[96];
            tab[i+87*stg_first]=tab[i+87*stg_first]+z[0]*tmp1[0]
            +z[87]*tmp1[1]+z[77]*tmp1[2]+z[67]*tmp1[3]+z[57]*tmp1[4]+z[47]*tmp1[5]+z[37]*tmp1[6]+z[27]*tmp1[7]+z[17]*tmp1[8]+z[7]*tmp1[9]+z[94]*tmp1[10]
            +z[84]*tmp1[11]+z[74]*tmp1[12]+z[64]*tmp1[13]+z[54]*tmp1[14]+z[44]*tmp1[15]+z[34]*tmp1[16]+z[24]*tmp1[17]+z[14]*tmp1[18]+z[4]*tmp1[19]+z[91]*tmp1[20]
            +z[81]*tmp1[21]+z[71]*tmp1[22]+z[61]*tmp1[23]+z[51]*tmp1[24]+z[41]*tmp1[25]+z[31]*tmp1[26]+z[21]*tmp1[27]+z[11]*tmp1[28]+z[1]*tmp1[29]+z[88]*tmp1[30]
            +z[78]*tmp1[31]+z[68]*tmp1[32]+z[58]*tmp1[33]+z[48]*tmp1[34]+z[38]*tmp1[35]+z[28]*tmp1[36]+z[18]*tmp1[37]+z[8]*tmp1[38]+z[95]*tmp1[39]+z[85]*tmp1[40]
            +z[75]*tmp1[41]+z[65]*tmp1[42]+z[55]*tmp1[43]+z[45]*tmp1[44]+z[35]*tmp1[45]+z[25]*tmp1[46]+z[15]*tmp1[47]+z[5]*tmp1[48]+z[92]*tmp1[49]+z[82]*tmp1[50]
            +z[72]*tmp1[51]+z[62]*tmp1[52]+z[52]*tmp1[53]+z[42]*tmp1[54]+z[32]*tmp1[55]+z[22]*tmp1[56]+z[12]*tmp1[57]+z[2]*tmp1[58]+z[89]*tmp1[59]+z[79]*tmp1[60]
            +z[69]*tmp1[61]+z[59]*tmp1[62]+z[49]*tmp1[63]+z[39]*tmp1[64]+z[29]*tmp1[65]+z[19]*tmp1[66]+z[9]*tmp1[67]+z[96]*tmp1[68]+z[86]*tmp1[69]+z[76]*tmp1[70]
            +z[66]*tmp1[71]+z[56]*tmp1[72]+z[46]*tmp1[73]+z[36]*tmp1[74]+z[26]*tmp1[75]+z[16]*tmp1[76]+z[6]*tmp1[77]+z[93]*tmp1[78]+z[83]*tmp1[79]+z[73]*tmp1[80]
            +z[63]*tmp1[81]+z[53]*tmp1[82]+z[43]*tmp1[83]+z[33]*tmp1[84]+z[23]*tmp1[85]+z[13]*tmp1[86]+z[3]*tmp1[87]+z[90]*tmp1[88]+z[80]*tmp1[89]+z[70]*tmp1[90]
            +z[60]*tmp1[91]+z[50]*tmp1[92]+z[40]*tmp1[93]+z[30]*tmp1[94]+z[20]*tmp1[95]+z[10]*tmp1[96];
            tab[i+88*stg_first]=tab[i+88*stg_first]+z[0]*tmp1[0]
            +z[88]*tmp1[1]+z[79]*tmp1[2]+z[70]*tmp1[3]+z[61]*tmp1[4]+z[52]*tmp1[5]+z[43]*tmp1[6]+z[34]*tmp1[7]+z[25]*tmp1[8]+z[16]*tmp1[9]+z[7]*tmp1[10]
            +z[95]*tmp1[11]+z[86]*tmp1[12]+z[77]*tmp1[13]+z[68]*tmp1[14]+z[59]*tmp1[15]+z[50]*tmp1[16]+z[41]*tmp1[17]+z[32]*tmp1[18]+z[23]*tmp1[19]+z[14]*tmp1[20]
            +z[5]*tmp1[21]+z[93]*tmp1[22]+z[84]*tmp1[23]+z[75]*tmp1[24]+z[66]*tmp1[25]+z[57]*tmp1[26]+z[48]*tmp1[27]+z[39]*tmp1[28]+z[30]*tmp1[29]+z[21]*tmp1[30]
            +z[12]*tmp1[31]+z[3]*tmp1[32]+z[91]*tmp1[33]+z[82]*tmp1[34]+z[73]*tmp1[35]+z[64]*tmp1[36]+z[55]*tmp1[37]+z[46]*tmp1[38]+z[37]*tmp1[39]+z[28]*tmp1[40]
            +z[19]*tmp1[41]+z[10]*tmp1[42]+z[1]*tmp1[43]+z[89]*tmp1[44]+z[80]*tmp1[45]+z[71]*tmp1[46]+z[62]*tmp1[47]+z[53]*tmp1[48]+z[44]*tmp1[49]+z[35]*tmp1[50]
            +z[26]*tmp1[51]+z[17]*tmp1[52]+z[8]*tmp1[53]+z[96]*tmp1[54]+z[87]*tmp1[55]+z[78]*tmp1[56]+z[69]*tmp1[57]+z[60]*tmp1[58]+z[51]*tmp1[59]+z[42]*tmp1[60]
            +z[33]*tmp1[61]+z[24]*tmp1[62]+z[15]*tmp1[63]+z[6]*tmp1[64]+z[94]*tmp1[65]+z[85]*tmp1[66]+z[76]*tmp1[67]+z[67]*tmp1[68]+z[58]*tmp1[69]+z[49]*tmp1[70]
            +z[40]*tmp1[71]+z[31]*tmp1[72]+z[22]*tmp1[73]+z[13]*tmp1[74]+z[4]*tmp1[75]+z[92]*tmp1[76]+z[83]*tmp1[77]+z[74]*tmp1[78]+z[65]*tmp1[79]+z[56]*tmp1[80]
            +z[47]*tmp1[81]+z[38]*tmp1[82]+z[29]*tmp1[83]+z[20]*tmp1[84]+z[11]*tmp1[85]+z[2]*tmp1[86]+z[90]*tmp1[87]+z[81]*tmp1[88]+z[72]*tmp1[89]+z[63]*tmp1[90]
            +z[54]*tmp1[91]+z[45]*tmp1[92]+z[36]*tmp1[93]+z[27]*tmp1[94]+z[18]*tmp1[95]+z[9]*tmp1[96];
            tab[i+89*stg_first]=tab[i+89*stg_first]+z[0]*tmp1[0]
            +z[89]*tmp1[1]+z[81]*tmp1[2]+z[73]*tmp1[3]+z[65]*tmp1[4]+z[57]*tmp1[5]+z[49]*tmp1[6]+z[41]*tmp1[7]+z[33]*tmp1[8]+z[25]*tmp1[9]+z[17]*tmp1[10]
            +z[9]*tmp1[11]+z[1]*tmp1[12]+z[90]*tmp1[13]+z[82]*tmp1[14]+z[74]*tmp1[15]+z[66]*tmp1[16]+z[58]*tmp1[17]+z[50]*tmp1[18]+z[42]*tmp1[19]+z[34]*tmp1[20]
            +z[26]*tmp1[21]+z[18]*tmp1[22]+z[10]*tmp1[23]+z[2]*tmp1[24]+z[91]*tmp1[25]+z[83]*tmp1[26]+z[75]*tmp1[27]+z[67]*tmp1[28]+z[59]*tmp1[29]+z[51]*tmp1[30]
            +z[43]*tmp1[31]+z[35]*tmp1[32]+z[27]*tmp1[33]+z[19]*tmp1[34]+z[11]*tmp1[35]+z[3]*tmp1[36]+z[92]*tmp1[37]+z[84]*tmp1[38]+z[76]*tmp1[39]+z[68]*tmp1[40]
            +z[60]*tmp1[41]+z[52]*tmp1[42]+z[44]*tmp1[43]+z[36]*tmp1[44]+z[28]*tmp1[45]+z[20]*tmp1[46]+z[12]*tmp1[47]+z[4]*tmp1[48]+z[93]*tmp1[49]+z[85]*tmp1[50]
            +z[77]*tmp1[51]+z[69]*tmp1[52]+z[61]*tmp1[53]+z[53]*tmp1[54]+z[45]*tmp1[55]+z[37]*tmp1[56]+z[29]*tmp1[57]+z[21]*tmp1[58]+z[13]*tmp1[59]+z[5]*tmp1[60]
            +z[94]*tmp1[61]+z[86]*tmp1[62]+z[78]*tmp1[63]+z[70]*tmp1[64]+z[62]*tmp1[65]+z[54]*tmp1[66]+z[46]*tmp1[67]+z[38]*tmp1[68]+z[30]*tmp1[69]+z[22]*tmp1[70]
            +z[14]*tmp1[71]+z[6]*tmp1[72]+z[95]*tmp1[73]+z[87]*tmp1[74]+z[79]*tmp1[75]+z[71]*tmp1[76]+z[63]*tmp1[77]+z[55]*tmp1[78]+z[47]*tmp1[79]+z[39]*tmp1[80]
            +z[31]*tmp1[81]+z[23]*tmp1[82]+z[15]*tmp1[83]+z[7]*tmp1[84]+z[96]*tmp1[85]+z[88]*tmp1[86]+z[80]*tmp1[87]+z[72]*tmp1[88]+z[64]*tmp1[89]+z[56]*tmp1[90]
            +z[48]*tmp1[91]+z[40]*tmp1[92]+z[32]*tmp1[93]+z[24]*tmp1[94]+z[16]*tmp1[95]+z[8]*tmp1[96];
            tab[i+90*stg_first]=tab[i+90*stg_first]+z[0]*tmp1[0]
            +z[90]*tmp1[1]+z[83]*tmp1[2]+z[76]*tmp1[3]+z[69]*tmp1[4]+z[62]*tmp1[5]+z[55]*tmp1[6]+z[48]*tmp1[7]+z[41]*tmp1[8]+z[34]*tmp1[9]+z[27]*tmp1[10]
            +z[20]*tmp1[11]+z[13]*tmp1[12]+z[6]*tmp1[13]+z[96]*tmp1[14]+z[89]*tmp1[15]+z[82]*tmp1[16]+z[75]*tmp1[17]+z[68]*tmp1[18]+z[61]*tmp1[19]+z[54]*tmp1[20]
            +z[47]*tmp1[21]+z[40]*tmp1[22]+z[33]*tmp1[23]+z[26]*tmp1[24]+z[19]*tmp1[25]+z[12]*tmp1[26]+z[5]*tmp1[27]+z[95]*tmp1[28]+z[88]*tmp1[29]+z[81]*tmp1[30]
            +z[74]*tmp1[31]+z[67]*tmp1[32]+z[60]*tmp1[33]+z[53]*tmp1[34]+z[46]*tmp1[35]+z[39]*tmp1[36]+z[32]*tmp1[37]+z[25]*tmp1[38]+z[18]*tmp1[39]+z[11]*tmp1[40]
            +z[4]*tmp1[41]+z[94]*tmp1[42]+z[87]*tmp1[43]+z[80]*tmp1[44]+z[73]*tmp1[45]+z[66]*tmp1[46]+z[59]*tmp1[47]+z[52]*tmp1[48]+z[45]*tmp1[49]+z[38]*tmp1[50]
            +z[31]*tmp1[51]+z[24]*tmp1[52]+z[17]*tmp1[53]+z[10]*tmp1[54]+z[3]*tmp1[55]+z[93]*tmp1[56]+z[86]*tmp1[57]+z[79]*tmp1[58]+z[72]*tmp1[59]+z[65]*tmp1[60]
            +z[58]*tmp1[61]+z[51]*tmp1[62]+z[44]*tmp1[63]+z[37]*tmp1[64]+z[30]*tmp1[65]+z[23]*tmp1[66]+z[16]*tmp1[67]+z[9]*tmp1[68]+z[2]*tmp1[69]+z[92]*tmp1[70]
            +z[85]*tmp1[71]+z[78]*tmp1[72]+z[71]*tmp1[73]+z[64]*tmp1[74]+z[57]*tmp1[75]+z[50]*tmp1[76]+z[43]*tmp1[77]+z[36]*tmp1[78]+z[29]*tmp1[79]+z[22]*tmp1[80]
            +z[15]*tmp1[81]+z[8]*tmp1[82]+z[1]*tmp1[83]+z[91]*tmp1[84]+z[84]*tmp1[85]+z[77]*tmp1[86]+z[70]*tmp1[87]+z[63]*tmp1[88]+z[56]*tmp1[89]+z[49]*tmp1[90]
            +z[42]*tmp1[91]+z[35]*tmp1[92]+z[28]*tmp1[93]+z[21]*tmp1[94]+z[14]*tmp1[95]+z[7]*tmp1[96];
            tab[i+91*stg_first]=tab[i+91*stg_first]+z[0]*tmp1[0]
            +z[91]*tmp1[1]+z[85]*tmp1[2]+z[79]*tmp1[3]+z[73]*tmp1[4]+z[67]*tmp1[5]+z[61]*tmp1[6]+z[55]*tmp1[7]+z[49]*tmp1[8]+z[43]*tmp1[9]+z[37]*tmp1[10]
            +z[31]*tmp1[11]+z[25]*tmp1[12]+z[19]*tmp1[13]+z[13]*tmp1[14]+z[7]*tmp1[15]+z[1]*tmp1[16]+z[92]*tmp1[17]+z[86]*tmp1[18]+z[80]*tmp1[19]+z[74]*tmp1[20]
            +z[68]*tmp1[21]+z[62]*tmp1[22]+z[56]*tmp1[23]+z[50]*tmp1[24]+z[44]*tmp1[25]+z[38]*tmp1[26]+z[32]*tmp1[27]+z[26]*tmp1[28]+z[20]*tmp1[29]+z[14]*tmp1[30]
            +z[8]*tmp1[31]+z[2]*tmp1[32]+z[93]*tmp1[33]+z[87]*tmp1[34]+z[81]*tmp1[35]+z[75]*tmp1[36]+z[69]*tmp1[37]+z[63]*tmp1[38]+z[57]*tmp1[39]+z[51]*tmp1[40]
            +z[45]*tmp1[41]+z[39]*tmp1[42]+z[33]*tmp1[43]+z[27]*tmp1[44]+z[21]*tmp1[45]+z[15]*tmp1[46]+z[9]*tmp1[47]+z[3]*tmp1[48]+z[94]*tmp1[49]+z[88]*tmp1[50]
            +z[82]*tmp1[51]+z[76]*tmp1[52]+z[70]*tmp1[53]+z[64]*tmp1[54]+z[58]*tmp1[55]+z[52]*tmp1[56]+z[46]*tmp1[57]+z[40]*tmp1[58]+z[34]*tmp1[59]+z[28]*tmp1[60]
            +z[22]*tmp1[61]+z[16]*tmp1[62]+z[10]*tmp1[63]+z[4]*tmp1[64]+z[95]*tmp1[65]+z[89]*tmp1[66]+z[83]*tmp1[67]+z[77]*tmp1[68]+z[71]*tmp1[69]+z[65]*tmp1[70]
            +z[59]*tmp1[71]+z[53]*tmp1[72]+z[47]*tmp1[73]+z[41]*tmp1[74]+z[35]*tmp1[75]+z[29]*tmp1[76]+z[23]*tmp1[77]+z[17]*tmp1[78]+z[11]*tmp1[79]+z[5]*tmp1[80]
            +z[96]*tmp1[81]+z[90]*tmp1[82]+z[84]*tmp1[83]+z[78]*tmp1[84]+z[72]*tmp1[85]+z[66]*tmp1[86]+z[60]*tmp1[87]+z[54]*tmp1[88]+z[48]*tmp1[89]+z[42]*tmp1[90]
            +z[36]*tmp1[91]+z[30]*tmp1[92]+z[24]*tmp1[93]+z[18]*tmp1[94]+z[12]*tmp1[95]+z[6]*tmp1[96];
            tab[i+92*stg_first]=tab[i+92*stg_first]+z[0]*tmp1[0]
            +z[92]*tmp1[1]+z[87]*tmp1[2]+z[82]*tmp1[3]+z[77]*tmp1[4]+z[72]*tmp1[5]+z[67]*tmp1[6]+z[62]*tmp1[7]+z[57]*tmp1[8]+z[52]*tmp1[9]+z[47]*tmp1[10]
            +z[42]*tmp1[11]+z[37]*tmp1[12]+z[32]*tmp1[13]+z[27]*tmp1[14]+z[22]*tmp1[15]+z[17]*tmp1[16]+z[12]*tmp1[17]+z[7]*tmp1[18]+z[2]*tmp1[19]+z[94]*tmp1[20]
            +z[89]*tmp1[21]+z[84]*tmp1[22]+z[79]*tmp1[23]+z[74]*tmp1[24]+z[69]*tmp1[25]+z[64]*tmp1[26]+z[59]*tmp1[27]+z[54]*tmp1[28]+z[49]*tmp1[29]+z[44]*tmp1[30]
            +z[39]*tmp1[31]+z[34]*tmp1[32]+z[29]*tmp1[33]+z[24]*tmp1[34]+z[19]*tmp1[35]+z[14]*tmp1[36]+z[9]*tmp1[37]+z[4]*tmp1[38]+z[96]*tmp1[39]+z[91]*tmp1[40]
            +z[86]*tmp1[41]+z[81]*tmp1[42]+z[76]*tmp1[43]+z[71]*tmp1[44]+z[66]*tmp1[45]+z[61]*tmp1[46]+z[56]*tmp1[47]+z[51]*tmp1[48]+z[46]*tmp1[49]+z[41]*tmp1[50]
            +z[36]*tmp1[51]+z[31]*tmp1[52]+z[26]*tmp1[53]+z[21]*tmp1[54]+z[16]*tmp1[55]+z[11]*tmp1[56]+z[6]*tmp1[57]+z[1]*tmp1[58]+z[93]*tmp1[59]+z[88]*tmp1[60]
            +z[83]*tmp1[61]+z[78]*tmp1[62]+z[73]*tmp1[63]+z[68]*tmp1[64]+z[63]*tmp1[65]+z[58]*tmp1[66]+z[53]*tmp1[67]+z[48]*tmp1[68]+z[43]*tmp1[69]+z[38]*tmp1[70]
            +z[33]*tmp1[71]+z[28]*tmp1[72]+z[23]*tmp1[73]+z[18]*tmp1[74]+z[13]*tmp1[75]+z[8]*tmp1[76]+z[3]*tmp1[77]+z[95]*tmp1[78]+z[90]*tmp1[79]+z[85]*tmp1[80]
            +z[80]*tmp1[81]+z[75]*tmp1[82]+z[70]*tmp1[83]+z[65]*tmp1[84]+z[60]*tmp1[85]+z[55]*tmp1[86]+z[50]*tmp1[87]+z[45]*tmp1[88]+z[40]*tmp1[89]+z[35]*tmp1[90]
            +z[30]*tmp1[91]+z[25]*tmp1[92]+z[20]*tmp1[93]+z[15]*tmp1[94]+z[10]*tmp1[95]+z[5]*tmp1[96];
            tab[i+93*stg_first]=tab[i+93*stg_first]+z[0]*tmp1[0]
            +z[93]*tmp1[1]+z[89]*tmp1[2]+z[85]*tmp1[3]+z[81]*tmp1[4]+z[77]*tmp1[5]+z[73]*tmp1[6]+z[69]*tmp1[7]+z[65]*tmp1[8]+z[61]*tmp1[9]+z[57]*tmp1[10]
            +z[53]*tmp1[11]+z[49]*tmp1[12]+z[45]*tmp1[13]+z[41]*tmp1[14]+z[37]*tmp1[15]+z[33]*tmp1[16]+z[29]*tmp1[17]+z[25]*tmp1[18]+z[21]*tmp1[19]+z[17]*tmp1[20]
            +z[13]*tmp1[21]+z[9]*tmp1[22]+z[5]*tmp1[23]+z[1]*tmp1[24]+z[94]*tmp1[25]+z[90]*tmp1[26]+z[86]*tmp1[27]+z[82]*tmp1[28]+z[78]*tmp1[29]+z[74]*tmp1[30]
            +z[70]*tmp1[31]+z[66]*tmp1[32]+z[62]*tmp1[33]+z[58]*tmp1[34]+z[54]*tmp1[35]+z[50]*tmp1[36]+z[46]*tmp1[37]+z[42]*tmp1[38]+z[38]*tmp1[39]+z[34]*tmp1[40]
            +z[30]*tmp1[41]+z[26]*tmp1[42]+z[22]*tmp1[43]+z[18]*tmp1[44]+z[14]*tmp1[45]+z[10]*tmp1[46]+z[6]*tmp1[47]+z[2]*tmp1[48]+z[95]*tmp1[49]+z[91]*tmp1[50]
            +z[87]*tmp1[51]+z[83]*tmp1[52]+z[79]*tmp1[53]+z[75]*tmp1[54]+z[71]*tmp1[55]+z[67]*tmp1[56]+z[63]*tmp1[57]+z[59]*tmp1[58]+z[55]*tmp1[59]+z[51]*tmp1[60]
            +z[47]*tmp1[61]+z[43]*tmp1[62]+z[39]*tmp1[63]+z[35]*tmp1[64]+z[31]*tmp1[65]+z[27]*tmp1[66]+z[23]*tmp1[67]+z[19]*tmp1[68]+z[15]*tmp1[69]+z[11]*tmp1[70]
            +z[7]*tmp1[71]+z[3]*tmp1[72]+z[96]*tmp1[73]+z[92]*tmp1[74]+z[88]*tmp1[75]+z[84]*tmp1[76]+z[80]*tmp1[77]+z[76]*tmp1[78]+z[72]*tmp1[79]+z[68]*tmp1[80]
            +z[64]*tmp1[81]+z[60]*tmp1[82]+z[56]*tmp1[83]+z[52]*tmp1[84]+z[48]*tmp1[85]+z[44]*tmp1[86]+z[40]*tmp1[87]+z[36]*tmp1[88]+z[32]*tmp1[89]+z[28]*tmp1[90]
            +z[24]*tmp1[91]+z[20]*tmp1[92]+z[16]*tmp1[93]+z[12]*tmp1[94]+z[8]*tmp1[95]+z[4]*tmp1[96];
            tab[i+94*stg_first]=tab[i+94*stg_first]+z[0]*tmp1[0]
            +z[94]*tmp1[1]+z[91]*tmp1[2]+z[88]*tmp1[3]+z[85]*tmp1[4]+z[82]*tmp1[5]+z[79]*tmp1[6]+z[76]*tmp1[7]+z[73]*tmp1[8]+z[70]*tmp1[9]+z[67]*tmp1[10]
            +z[64]*tmp1[11]+z[61]*tmp1[12]+z[58]*tmp1[13]+z[55]*tmp1[14]+z[52]*tmp1[15]+z[49]*tmp1[16]+z[46]*tmp1[17]+z[43]*tmp1[18]+z[40]*tmp1[19]+z[37]*tmp1[20]
            +z[34]*tmp1[21]+z[31]*tmp1[22]+z[28]*tmp1[23]+z[25]*tmp1[24]+z[22]*tmp1[25]+z[19]*tmp1[26]+z[16]*tmp1[27]+z[13]*tmp1[28]+z[10]*tmp1[29]+z[7]*tmp1[30]
            +z[4]*tmp1[31]+z[1]*tmp1[32]+z[95]*tmp1[33]+z[92]*tmp1[34]+z[89]*tmp1[35]+z[86]*tmp1[36]+z[83]*tmp1[37]+z[80]*tmp1[38]+z[77]*tmp1[39]+z[74]*tmp1[40]
            +z[71]*tmp1[41]+z[68]*tmp1[42]+z[65]*tmp1[43]+z[62]*tmp1[44]+z[59]*tmp1[45]+z[56]*tmp1[46]+z[53]*tmp1[47]+z[50]*tmp1[48]+z[47]*tmp1[49]+z[44]*tmp1[50]
            +z[41]*tmp1[51]+z[38]*tmp1[52]+z[35]*tmp1[53]+z[32]*tmp1[54]+z[29]*tmp1[55]+z[26]*tmp1[56]+z[23]*tmp1[57]+z[20]*tmp1[58]+z[17]*tmp1[59]+z[14]*tmp1[60]
            +z[11]*tmp1[61]+z[8]*tmp1[62]+z[5]*tmp1[63]+z[2]*tmp1[64]+z[96]*tmp1[65]+z[93]*tmp1[66]+z[90]*tmp1[67]+z[87]*tmp1[68]+z[84]*tmp1[69]+z[81]*tmp1[70]
            +z[78]*tmp1[71]+z[75]*tmp1[72]+z[72]*tmp1[73]+z[69]*tmp1[74]+z[66]*tmp1[75]+z[63]*tmp1[76]+z[60]*tmp1[77]+z[57]*tmp1[78]+z[54]*tmp1[79]+z[51]*tmp1[80]
            +z[48]*tmp1[81]+z[45]*tmp1[82]+z[42]*tmp1[83]+z[39]*tmp1[84]+z[36]*tmp1[85]+z[33]*tmp1[86]+z[30]*tmp1[87]+z[27]*tmp1[88]+z[24]*tmp1[89]+z[21]*tmp1[90]
            +z[18]*tmp1[91]+z[15]*tmp1[92]+z[12]*tmp1[93]+z[9]*tmp1[94]+z[6]*tmp1[95]+z[3]*tmp1[96];
            tab[i+95*stg_first]=tab[i+95*stg_first]+z[0]*tmp1[0]
            +z[95]*tmp1[1]+z[93]*tmp1[2]+z[91]*tmp1[3]+z[89]*tmp1[4]+z[87]*tmp1[5]+z[85]*tmp1[6]+z[83]*tmp1[7]+z[81]*tmp1[8]+z[79]*tmp1[9]+z[77]*tmp1[10]
            +z[75]*tmp1[11]+z[73]*tmp1[12]+z[71]*tmp1[13]+z[69]*tmp1[14]+z[67]*tmp1[15]+z[65]*tmp1[16]+z[63]*tmp1[17]+z[61]*tmp1[18]+z[59]*tmp1[19]+z[57]*tmp1[20]
            +z[55]*tmp1[21]+z[53]*tmp1[22]+z[51]*tmp1[23]+z[49]*tmp1[24]+z[47]*tmp1[25]+z[45]*tmp1[26]+z[43]*tmp1[27]+z[41]*tmp1[28]+z[39]*tmp1[29]+z[37]*tmp1[30]
            +z[35]*tmp1[31]+z[33]*tmp1[32]+z[31]*tmp1[33]+z[29]*tmp1[34]+z[27]*tmp1[35]+z[25]*tmp1[36]+z[23]*tmp1[37]+z[21]*tmp1[38]+z[19]*tmp1[39]+z[17]*tmp1[40]
            +z[15]*tmp1[41]+z[13]*tmp1[42]+z[11]*tmp1[43]+z[9]*tmp1[44]+z[7]*tmp1[45]+z[5]*tmp1[46]+z[3]*tmp1[47]+z[1]*tmp1[48]+z[96]*tmp1[49]+z[94]*tmp1[50]
            +z[92]*tmp1[51]+z[90]*tmp1[52]+z[88]*tmp1[53]+z[86]*tmp1[54]+z[84]*tmp1[55]+z[82]*tmp1[56]+z[80]*tmp1[57]+z[78]*tmp1[58]+z[76]*tmp1[59]+z[74]*tmp1[60]
            +z[72]*tmp1[61]+z[70]*tmp1[62]+z[68]*tmp1[63]+z[66]*tmp1[64]+z[64]*tmp1[65]+z[62]*tmp1[66]+z[60]*tmp1[67]+z[58]*tmp1[68]+z[56]*tmp1[69]+z[54]*tmp1[70]
            +z[52]*tmp1[71]+z[50]*tmp1[72]+z[48]*tmp1[73]+z[46]*tmp1[74]+z[44]*tmp1[75]+z[42]*tmp1[76]+z[40]*tmp1[77]+z[38]*tmp1[78]+z[36]*tmp1[79]+z[34]*tmp1[80]
            +z[32]*tmp1[81]+z[30]*tmp1[82]+z[28]*tmp1[83]+z[26]*tmp1[84]+z[24]*tmp1[85]+z[22]*tmp1[86]+z[20]*tmp1[87]+z[18]*tmp1[88]+z[16]*tmp1[89]+z[14]*tmp1[90]
            +z[12]*tmp1[91]+z[10]*tmp1[92]+z[8]*tmp1[93]+z[6]*tmp1[94]+z[4]*tmp1[95]+z[2]*tmp1[96];
            tab[i+96*stg_first]=tab[i+96*stg_first]+z[0]*tmp1[0]
            +z[96]*tmp1[1]+z[95]*tmp1[2]+z[94]*tmp1[3]+z[93]*tmp1[4]+z[92]*tmp1[5]+z[91]*tmp1[6]+z[90]*tmp1[7]+z[89]*tmp1[8]+z[88]*tmp1[9]+z[87]*tmp1[10]
            +z[86]*tmp1[11]+z[85]*tmp1[12]+z[84]*tmp1[13]+z[83]*tmp1[14]+z[82]*tmp1[15]+z[81]*tmp1[16]+z[80]*tmp1[17]+z[79]*tmp1[18]+z[78]*tmp1[19]+z[77]*tmp1[20]
            +z[76]*tmp1[21]+z[75]*tmp1[22]+z[74]*tmp1[23]+z[73]*tmp1[24]+z[72]*tmp1[25]+z[71]*tmp1[26]+z[70]*tmp1[27]+z[69]*tmp1[28]+z[68]*tmp1[29]+z[67]*tmp1[30]
            +z[66]*tmp1[31]+z[65]*tmp1[32]+z[64]*tmp1[33]+z[63]*tmp1[34]+z[62]*tmp1[35]+z[61]*tmp1[36]+z[60]*tmp1[37]+z[59]*tmp1[38]+z[58]*tmp1[39]+z[57]*tmp1[40]
            +z[56]*tmp1[41]+z[55]*tmp1[42]+z[54]*tmp1[43]+z[53]*tmp1[44]+z[52]*tmp1[45]+z[51]*tmp1[46]+z[50]*tmp1[47]+z[49]*tmp1[48]+z[48]*tmp1[49]+z[47]*tmp1[50]
            +z[46]*tmp1[51]+z[45]*tmp1[52]+z[44]*tmp1[53]+z[43]*tmp1[54]+z[42]*tmp1[55]+z[41]*tmp1[56]+z[40]*tmp1[57]+z[39]*tmp1[58]+z[38]*tmp1[59]+z[37]*tmp1[60]
            +z[36]*tmp1[61]+z[35]*tmp1[62]+z[34]*tmp1[63]+z[33]*tmp1[64]+z[32]*tmp1[65]+z[31]*tmp1[66]+z[30]*tmp1[67]+z[29]*tmp1[68]+z[28]*tmp1[69]+z[27]*tmp1[70]
            +z[26]*tmp1[71]+z[25]*tmp1[72]+z[24]*tmp1[73]+z[23]*tmp1[74]+z[22]*tmp1[75]+z[21]*tmp1[76]+z[20]*tmp1[77]+z[19]*tmp1[78]+z[18]*tmp1[79]+z[17]*tmp1[80]
            +z[16]*tmp1[81]+z[15]*tmp1[82]+z[14]*tmp1[83]+z[13]*tmp1[84]+z[12]*tmp1[85]+z[11]*tmp1[86]+z[10]*tmp1[87]+z[9]*tmp1[88]+z[8]*tmp1[89]+z[7]*tmp1[90]
            +z[6]*tmp1[91]+z[5]*tmp1[92]+z[4]*tmp1[93]+z[3]*tmp1[94]+z[2]*tmp1[95]+z[1]*tmp1[96];
        }
    }
///////////////////////////////////////////////////





//this is new in that method and fi:

//when you want to have equal results that are in false modificator in normal FFT then change this:
/*
 fun_fourier_transform_FFT_radix_4_N_256_official
{
    for(int j=0;j<N;j++)
    {
      tab[j].real() =tab[j].real()*2/N;
      tab[j].imag() =tab[j].imag()*2/N;
    }
}
//and:

fun_inverse_fourier_transform_FFT_radix_4_N_256_official
{
    for(int j=0;j<N;j++)
    {
      tab[j].real() =tab[j].real()*0.5;
      tab[j].imag() =tab[j].imag()*0.5;
    }
}

//for official modificator that is only in inverse FFT:

 fun_fourier_transform_FFT_radix_4_N_256_official
{

}
fun_inverse_fourier_transform_FFT_radix_4_N_256_official
{
    for(int i=0;i<N;i++)
    {
        tablica1[0][i]=tablica1[0][i]*1/(float)N;
        tablica1[1][i]=tablica1[1][i]*1/(float)N;
    }
}

//haven't try it with other function that cos(x)+jsin(x)=sin(x+pi/2)+jsin(x) that is mirror inverse
*/

