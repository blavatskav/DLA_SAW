

#include <stdlib.h>
#include <time.h>
#include <stdio.h>
#include <math.h>

# define dioscuri 2
# define build 30

// number of SAW steps
# define Nmax 120
# define Nmin 10
// number of DLA steps
# define NNmax 1000
# define epsilon 0.01
# define maxrad 2000
# define pobud 30
# define sep 0


// DLA=1 - narostajut klastery 
double zvjazok[Nmax+2],bondfinal[Nmax+2][pobud+2],bond[Nmax+2][build+2],sup,veps,den,moment[10][Nmax+20],
weight[Nmax+20],weightfinal[Nmax+20],gir,girr,
  Rgfinal[NNmax+20], Rg[NNmax+20][pobud+2],dysrg[NNmax+20],
  periferi[Nmax+20][pobud+2],periferifinal[Nmax+20],dysper[Nmax+20],
  shuffle;
 
 int tak,visit[Nmax+2],Rrun,p1,p2,i1,i2,pp,dd,q,h,hh,zanjato[maxrad+10][maxrad+10],poslid[maxrad+10][maxrad+10],
 zanjatoSAW[maxrad+10][maxrad+10],x[Nmax+2],y[Nmax+2],
 adsorb[8],NN,dostup[8],seedx,fre,
 movex[8], movez[8],
 movey[8],seedy,ch,flight,coordx[NNmax+Nmax+2],coordy[NNmax+Nmax+2],
 zapuskx[maxrad+10],zapusky[maxrad+20],cluster,cl,c, gy,
 count,  i,nn,neww,test,enrich, zap,nomer,prob[Nmax+2],rr,vidvidano[Nmax+10],veb2,
 f,ff,stepx[9],stepy[9],stepz[9],func,kk[9],veb,dod,N,radx,rad,numeri,
 yes,dla,dlaSAW, xx,yy,krok,povt,cc,j,l,p,prun,g,bah,zag3,fff,
  gm, vklad[build+1];      
      
      
      
 #define IM1 2147483563
#define IM2 2147483399
#define AM (1.0/IM1)
#define IMM1 (IM1-1)
#define IA1 40014
#define IA2 40692
#define IQ1 53668
#define IQ2 52774
#define IR1 12211
#define IR2 3791
#define NTAB 32
#define NDIV (1+IMM1/NTAB)
#define EPS 1.2e-7

#define RNMX (1.0-EPS)
int M,jjj;
long k;
static long idum2=123456789;
static long iy=0;
static long iv[NTAB];
float temp;
static long idum;
      
 main ()  
{ 
       
// for random numbers generations


// start dlja vypadkovyh chysel
 idum = -100000;
 void rrr()
{
if (idum <= 0) 
{ 
if (-(idum) < 1) idum=1; 
else idum = -(idum);
idum2=(idum);
for (jjj=NTAB+7;jjj>=0;jjj--) {
k=(idum)/IQ1;
idum=IA1*(idum-k*IQ1)-k*IR1;
if (idum < 0) idum += IM1;
if (jjj < NTAB) iv[jjj] = idum;
}
iy=iv[0];
}
k=(idum)/IQ1; 
idum=IA1*(idum-k*IQ1)-k*IR1; 
if (idum < 0) idum += IM1; 
k=idum2/IQ2;
idum2=IA2*(idum2-k*IQ2)-k*IR2; 
if (idum2 < 0) idum2 += IM2;
jjj=iy/NDIV; 
iy=iv[jjj]-idum2; 
iv[jjj] = idum; 
if (iy < 1) iy += IMM1;
if ((temp=AM*iy) > RNMX) {
shuffle=RNMX;}
else {
shuffle=temp;}
}







int buildclaster;

buildclaster=build;

FILE  *data1, *data2,*data3,*data4,*data5, *data,*kroky, *percfile;
int tri,a,b,vybir,s;      

 srand((unsigned)time(NULL));
 a = 1;
 b = 4; 
 int x1,y1;
 
 char name[50]; char name1[50];char name2[50]; char name3[50]; char name5[50];char name4[50];
 char perc[50];
 char namekroky[50];  
 
 char str1[50];
 char str2 [50];  

if(dioscuri==1)
{sprintf (name2,"/mnt/dioscuri-nas/Viktoria/DLA/periferieps%1.2f.dat",epsilon);
data2 = fopen(name2,"w+");
 }
 
 
           
//else{sprintf(name2,"periferideps%1.2f.dat",epsilon);} 



 if(dioscuri==1)
{sprintf (name3,"/mnt/dioscuri-nas/Viktoria/DLA/momenteps%1.2f.dat",epsilon);
data3 = fopen(name3,"w+");
 } 
 
                
//else{sprintf(name3,"momenteps%1.2f.dat",epsilon);} 


rad=120;
radx=400;
// zapusk chastynok z poverhni sfery radiusa 1000 
nomer=0;
for(xx=0;xx<=rad;++xx)
{                      
yy=sqrt(rad*rad-xx*xx);
if(yy>=0 && yy<=rad)
{                                             
nomer=nomer+1;
zapuskx[nomer]=xx+radx;zapusky[nomer]=yy+radx;
//printf("xx %d yy %d zz %d nomer %d \n", xx,yy,zz,nomer);                                         
nomer=nomer+1;
zapuskx[nomer]=xx+radx;zapusky[nomer]=(-1)*yy+radx; 
//printf("xx %d yy %d nomer %d \n", xx,yy,nomer);                                         
nomer=nomer+1;
zapuskx[nomer]=(-1)*xx+radx;zapusky[nomer]=yy+radx; 
//printf("xx %d yy %d nomer %d \n", xx,yy,nomer);                                         
nomer=nomer+1;
zapuskx[nomer]=(-1)*xx+radx;zapusky[nomer]=(-1)*yy+radx; 
//printf("xx %d yy %d nomer %d \n", xx,yy,nomer);                                                                                          
}
}
       
                       
for(N=Nmin;N<=Nmax+1;N=N+10)
{

if(dioscuri==2)
{sprintf (name5,"/mnt/dioscuri-nas/Viktoria/DLA/bondN%deps%1.2f.dat",N,epsilon);
data5 = fopen(name5,"w+");
 }
/*if(dioscuri==1)
{sprintf (name1,"/mnt/dioscuri-nas/Viktoria/DLA/probDLAN%d.dat",N); }                
else{sprintf(name1,"weightDLAN%d.dat",N);} 
data1 = fopen(name1,"w+");
  */                          


               
//else{sprintf(name4,"RgDLAN%deps%1.2f.dat",N,epsilon);} 

                            
   zvjazok[N] =0;                        
  periferifinal[N]=0;
  dysper[N]=0;
  for(q=0;q<=10;++q)
  {moment[q][N]=0;}
  for(i=0;i<=NNmax;++i)
  {Rgfinal[i]=0;
  dysrg[i]=0;}
 
 

for(pp=1;pp<=pobud;++pp)
{
printf("SAW N %d  pobudova %d epsilon %f \n",N, pp, epsilon);                                          
                                          
 radx=400;                                          
 for(j=0;j<=N;++j)
 {weight[j]=0.0;
 }
 
 
 periferi[N][pp]=0;
 bondfinal[N][pp]=0;
  for(j=0;j<=NNmax;++j)
  {Rg[j][pp]=0.0;}
 
 
                                    
begin1: for (krok=0;krok<=N;++krok)
 { 
  x[krok]=0;y[krok]=0;
 } 
 
 x[0]= radx; y[0] = radx; 
 i  =  1; 
 
 while (i<=N)
 { 
//printf("i %d  N %d \n", i, N);

	nupa:  
          stepx[0] = x[i-1]+1; stepy[0] = y[i-1]; 
          stepx[1] = x[i-1]-1; stepy[1] = y[i-1]; 
          stepx[2] = x[i-1]; stepy[2] = y[i-1]+1;  
          stepx[3] = x[i-1]; stepy[3] = y[i-1]-1; 
          
        dd=0;    
	  prob[i] = 4;   
	  kk[0] = 1; kk[1] =1; kk[2] =1; kk[3] =1; kk[4] =1;kk[5] =1;
for (j=0; j<4; ++j)
         {   
	   for(krok=0;krok<i;++krok)
	   {
             
                      //    printf("dd %d j %d krok %d ff %d x %d y %d z %d \n",dd,j,krok,ff,x[krok][ff],y[krok][ff],z[krok][ff]);
       if (stepx[j] == x[krok] && stepy[j] == y[krok] )                         
                     {      
		     
	   
		     {
             kk[j]=0; 
             prob[i]=prob[i]-1;
             }
             
             }		       
                 }}		        
          
		   // printf("here \n");
		if(prob[i]<=0){goto begin1;}     
/* generacija vypadkovogo chysla */
vybirka1: vybir = rand()%(4-1+1)+1;
g = vybir -1;

 if (kk[g] == 1 )
     {
     x[i] =  stepx[g]; 
     y[i] =  stepy[g];
     
     }                        
 if (kk[g] == 0 )
    {
    goto vybirka1; 
    }
   // printf( "g %d %d %d  %d\n", g,x[i],y[i],z[i]);
i=i+1;
 }


numeri=1;

//if(dioscuri==1)
//{sprintf (name,"/mnt/dioscuri-nas/Viktoria/DLA/SAWN%deps%1.2f.dat",N,epsilon);}
//else{  sprintf(name,"SAWN%deps%1.2f.dat",N,epsilon);} 
//data = fopen(name,"w+");
for(j=0;j<=N;++j)
{              
coordx[numeri]=x[j];coordy[numeri]=y[j]; 
//printf("numer %d %d  %d %d\n", numeri, coordx[numeri],coordy[numeri],coordz[numeri]);
weight[numeri]=0;
//fprintf(data, "%d  %d \n", coordx[numeri],coordy[numeri]);
numeri=numeri+1;
}
//fclose(data);  

if(dioscuri==0)
{
sprintf(name,"SAWN%deps%1.2f.dat",N,epsilon);
data = fopen(name,"w+");
for(j=0;j<=N;++j)
{               
//printf("numer %d %d  %d %d\n", numeri, coordx[numeri],coordy[numeri],coordz[numeri]);
//weight[numeri]=0;
fprintf(data, "%d  %d \n", x[j],y[j]);
}
fclose(data);  
}

numeri=numeri-1;
// teper, na kozhnomu klasteri zapuskajemo 

printf("numeri %d \n",numeri);


/*sprintf(name,"Ring%d.dat",rad); 
data = fopen(name,"w+");
 for (krok=1; krok<=nomer;++krok)
 {
//printf("krok %d  fugacity %f popravka %f \n", krok,  fugacity[krok],sqrt(dysfug[krok]));
fprintf(data, "%d  %d \n", zapuskx[krok],zapusky[krok]);
 }
fclose(data);     
*/

cluster=1;
 while (cluster<=buildclaster) 
{
Rrun=0;
  // coordinate of seed
for(i=0;i<=maxrad;++i)
{ 
for(j=0;j<=maxrad;++j)
{ 
zanjato[i][j]=0;
zanjatoSAW[i][j]=0;
poslid[i][j]=0;
} 
}

for(i=0;i<=NNmax+N;++i)
{
coordx[i]=0;
coordy[i]=0;
}
numeri=1;
for(j=0;j<=N;++j)
{              
coordx[numeri]=x[j];coordy[numeri]=y[j]; 
//printf("num %d x %d y %d \n",numeri,coordx[numeri],coordy[numeri]);
numeri=numeri+1;
}
numeri=numeri-1;     
// zadajemo jak seed - trajektorija SAW
for(krok=0;krok<=N;++krok)
{    
zanjato[coordx[krok]][coordy[krok]]=1;
zanjatoSAW[coordx[krok]][coordy[krok]]=1;
poslid[coordx[krok]][coordy[krok]]=krok;
vidvidano[krok]=0;
bond[krok][cluster]=0;
}






for(j=0;j<=N;++j)
 {
 visit[j]=0;
}
printf("cluster %d N %d\n",cluster, N);
           
// nalitajuchi chastynky        
for(flight=1;flight<=NNmax;++flight)
 {     

//printf("flight %d radius %d nomer %d \n",flight, rad,nomer);






zapu: rrr();
//sup=shuffle;
//printf("sss %f nomer %d numeri %d \n",shuffle, nomer,flight+numeri);
zap=(shuffle)*10000+1;
//printf("zap %d  \n",zap);
if(zap>nomer){goto zapu;}
//printf(" !!! zap %d \n",zap);
coordx[flight+numeri]=zapuskx[zap];
coordy[flight+numeri]=zapusky[zap];
//printf("flight %d  x %d y %d \n",flight,zapuskx[zap],zapusky[zap]); 
yes=0; 
while(yes<1)
{   
    dla=0;
   dlaSAW=0;
   fre=0;
// chastynka zdijsnuje RW doky ne dosjagne seed
for(i=1;i<=6;++i){dostup[i]=0;}
        movex[1] = coordx[flight+numeri]+1; movey[1] =coordy[flight+numeri]; 
        movex[2] = coordx[flight+numeri]-1; movey[2] = coordy[flight+numeri]; 
        movex[3] = coordx[flight+numeri];   movey[3] =coordy[flight+numeri]+1;
        movex[4] = coordx[flight+numeri];   movey[4] = coordy[flight+numeri]-1; 
        // perevirka, chy ne vyhodymo za mezhi kola
for(i=1;i<=4;++i){ if (movex[i]<0 ^ movey[i]<0) 
 {dostup[i]=1;fre=fre+1;}}
 for(i=1;i<=4;++i){ if (sqrt((movex[i]-radx)*(movex[i]-radx)+(movey[i]-radx)*(movey[i]-radx))>rad) 
 {dostup[i]=1;fre=fre+1;}}
 for(i=1;i<=4;++i){ if (zanjato[movex[i]][movey[i]]==1)  {dla=1;dostup[i]=1;fre=fre+1;}}
 for(i=1;i<=4;++i){ if (zanjatoSAW[movex[i]][movey[i]]==1)  {dlaSAW=1;
//printf("krokSAW %d  coordx %d coordy %d coordz %d\n ",h,movex[i],movey[i],movez[i]);
 }}
 
 if(fre>=4)
 {goto zapu;}
 

 // nashtovhulysja na SAW
 if(dlaSAW==1)
 {
 checksaw: rrr();veb=shuffle*1000+1;
  //        printf("veb %d \n",veb);
          if(veb<=200){ch=1;}
          if(veb>200 && veb<=400){ch=2;}
          if(veb>400 && veb<=600){ch=3;}
          if(veb>600 && veb<=800){ch=4;}
          if(veb>800) {goto checksaw;}
          //if(dostup[ch]==1) {goto checksaw;}
if(zanjatoSAW[movex[ch]][movey[ch]]==1)
 {h=poslid[movex[ch]][movey[ch]];
 visit[h]=visit[h]+1.0;
// printf("h %d\n",h);
 //if(vidvidano[h]==0)
 {
 // number of particles attached to the SAW
 periferi[N][pp]=periferi[N][pp]+1.0; 
 }
 weight[h]=weight[h]+1.0;
 visit[h]=visit[h]+1;
 yes=1;
 zanjato[coordx[flight+numeri]][coordy[flight+numeri]]=1; 
// printf("zanjatoSAW  krok %d x %d y %d periferi %d weight %f \n", h,movex[ch],movey[ch],periferi[N],weight[h]);
 break;
 }
 else {goto checksaw;}
}
 
// nashtovhulysja na inshu chastynku 
 if(dla==1)
 {      
// z imovirnistju epsilon nalypajemo
rrr();veps=shuffle;
if(veps<=epsilon)
{                 
zanjato[coordx[flight+numeri]][coordy[flight+numeri]]=1;  
//printf("veps %f \n",veps);       
yes=1;
break;
}
else {checkeps: rrr();
          veb2=shuffle*1000+1;
          if(veb2<=200){ch=1;}
          if(veb2>200 && veb2<=400){ch=2;}
          if(veb2>400 && veb2<=600){ch=3;}
          if(veb2>600 && veb2<=800){ch=4;}
          if(veb2>800) {goto checkeps;}
           if(dostup[ch]==1){goto checkeps;}
 coordx[flight+numeri]=movex[ch];
 coordy[flight+numeri]=movey[ch];
  yes=0;}
}

 else
 {
 check2: rrr();
          veb2=shuffle*1000+1;
          if(veb2<=200){ch=1;}
          if(veb2>200 && veb2<=400){ch=2;}
          if(veb2>400 && veb2<=600){ch=3;}
          if(veb2>600 && veb2<=800){ch=4;}
          if(veb2>800) {goto check2;}
           if(dostup[ch]==1){goto check2;}
 coordx[flight+numeri]=movex[ch];
 coordy[flight+numeri]=movey[ch];
  yes=0;
}

// yes
}


   if(sep==1)
{
   girr=0.0;
	for (i1=numeri+1;i1<=(flight+numeri);++i1)
			  {
    //        printf("j %d coordx %d coordy %d \n",krok,coordx[krok],coordy[krok]);                     
     for (i2=i1+1;i2<=(flight+numeri);++i2)
			  {
		
            girr=girr+(coordx[i1]-coordx[i2])*(coordx[i1]-coordx[i2])+(coordy[i1]-coordy[i2])*(coordy[i1]-coordy[i2]);
		   }
           }  
           
           girr=girr/((flight)*(flight));  
		// rr=(gir/(2*(flight+numeri)*(flight+numeri)))+0.5;
        // probab[rr][N]=probab[rr][N]+1;
}
else
{
   girr=0.0;
	for (i1=1;i1<=(flight+numeri);++i1)
			  {
    //        printf("j %d coordx %d coordy %d \n",krok,coordx[krok],coordy[krok]);                     
     for (i2=i1+1;i2<=(flight+numeri);++i2)
			  {
		
            girr=girr+(coordx[i1]-coordx[i2])*(coordx[i1]-coordx[i2])+(coordy[i1]-coordy[i2])*(coordy[i1]-coordy[i2]);
		   }
           }  
           
           girr=girr/((flight+numeri)*(flight+numeri));  
		// rr=(gir/(2*(flight+numeri)*(flight+numeri)))+0.5;
        // probab[rr][N]=probab[rr][N]+1;
}


Rg[flight][pp]=Rg[flight][pp]+1.0*girr; 
Rrun=sqrt(1.0*girr); 
//printf("flight %d Rg %f  \n",flight,Rg[flight][pp]);
//printf("  Rrun %d\n",Rrun);
//flight


}

tak=0;
for(h=1;h<=N;++h)
 {
 if(visit[h]>0)
 {
 printf("krok %d visit %d \n",h,visit[h]);
 tak=tak+1;
 bond[N][cluster]=bond[N][cluster]+visit[h];
 }
}
bond[N][cluster]=1.0*bond[N][cluster]/tak;
printf("bond %f \n",bond[N][cluster]);
if(dioscuri==0)
{
 sprintf(name,"DLAN%deps%1.2fcluster%d.dat",N,epsilon,cluster);
data = fopen(name,"w+");
for(i=1;i<=NNmax;++i)
{
fprintf(data, "%d  %d  \n", coordx[i+numeri],coordy[i+numeri]);
}
fclose(data);   
}

bondfinal[N][pp]=bondfinal[N][pp]+bond[N][cluster];
cluster=cluster+1;
//cluster
}




for(i=1;i<=NNmax;++i)
{
Rg[i][pp]=Rg[i][pp]/build;
}
periferi[N][pp]=periferi[N][pp]/build; 
bondfinal[N][pp]=bondfinal[N][pp]/build;
printf("bondfinal %f \n",bondfinal[N][pp]);

printf("N %d \n", N);
den=0; 
for(i=0;i<=N;++i)
{
weight[i]=weight[i]/(build);
}

for(q=0;q<=9;q=q+1)
{
for(i=0;i<=N;++i)
{
moment[q][N]=moment[q][N]+pow(weight[i],q);
}
}


if(dioscuri==0)
{sprintf (name,"/mnt/dioscuri-nas/Viktoria/DLA/weightN%deps%1.2fpobud%d.dat",N,epsilon,pp); 
data = fopen(name,"w+");

for(i=0;i<=numeri;++i)
{
//if(weight[i]>0)
{                      
fprintf(data, " %f  \n", 1.0*weight[i]);
}
}
fclose(data);
}


//else{sprintf(name,"weightN%deps%1.2f.dat",N,epsilon);} 

// pobudova
}

for(krok=1;krok<=pobud;++krok)
{
periferifinal[N]=periferifinal[N]+periferi[N][krok];
zvjazok[N]=zvjazok[N]+bondfinal[N][krok];
  for(i=1;i<=NNmax;++i)
  {Rgfinal[i]=Rgfinal[i]+Rg[i][krok];}
}



 for(krok=1;krok<=pobud;++krok)
{ 
for(j=1;j<=NNmax;++j)
{                   
dysrg[j]=dysrg[j]+(Rgfinal[j]/pobud-Rg[j][krok])*(Rgfinal[j]/pobud-Rg[j][krok]);
}
dysper[N]=dysper[N]+(periferifinal[N]/pobud-periferi[N][krok])*(periferifinal[N]/pobud-periferi[N][krok]);
}

if(dioscuri==1)
{
fprintf(data3, "%d ", N);
}



if(dioscuri==1)
{
for(q=0;q<=9;++q)
{
fprintf(data3, "  %f ", moment[q][N]/pobud);
}
fprintf(data3, "\n");
}


if(dioscuri==1)
{                     
fprintf(data2, "%d  %f %f \n", N, 1.0*periferifinal[N]/(pobud), sqrt(dysper[N])/pobud);
}

printf( "N %d  zv %f  \n", N, 1.0*zvjazok[N]/(pobud));

if(dioscuri==2)
{                     
fprintf(data5, "%d  %f  \n", N, 1.0*zvjazok[N]/(pobud));
printf( "N %d  zv %f  \n", N, 1.0*zvjazok[N]/(pobud));
}

/*
if(sep==1)
{
sprintf (name4,"/mnt/dioscuri-nas/Viktoria/DLA/seplongRgDLAN%deps%1.2f.dat",N,epsilon); 
data4 = fopen(name4,"w+"); 
for(i=1;i<=NNmax;++i)
{
fprintf(data4, "%d  %f  %f \n", i, Rgfinal[i]/(pobud),sqrt(dysrg[i])/pobud);
printf( "N %d  %f   \n", i, Rgfinal[i]/(pobud));
}

fclose(data4);
}
else
{
sprintf (name4,"/mnt/dioscuri-nas/Viktoria/DLA/longRgDLAN%deps%1.2f.dat",N,epsilon); 
data4 = fopen(name4,"w+"); 
for(i=1;i<=NNmax;++i)
{
fprintf(data4, "%d  %f  %f \n", i, Rgfinal[i]/(pobud),sqrt(dysrg[i])/pobud);
printf( "N %d  %f   \n", i, Rgfinal[i]/(pobud));
}

fclose(data4);
}
*/

if(dioscuri==2)
{fclose(data5);
}

//N
}


if(dioscuri==1)
{
fclose(data3);
fclose(data2);
}




system ("pause"); 
return 0; 
}




