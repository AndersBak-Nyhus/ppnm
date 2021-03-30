#include<stdio.h>
#include<stdlib.h>
#include<math.h>
#include<pthread.h>

struct params {int n,count;unsigned int seed;};

void* throw_points(void* arg){
	struct params * p = (struct params*)arg;
	p->count=0;
	for(int i=0;i < p->n;i++){
		double x=(double)rand_r(&(p->seed))/RAND_MAX;
		double y=(double)rand_r(&(p->seed))/RAND_MAX;
		if(x*x+y*y<1)p->count++;
	}

	return NULL;
}

int main(int argc, char** argv){
	int N=(int)1e6;
	if(argc>1) N=(int)atof(argv[1]);
	struct params p1 = {.n=N/3, .count=0, .seed=1};
	struct params p2 = {.n=N/3, .count=0, .seed=13};
	struct params p3 = {.n=N/3, .count=0, .seed=42};
	pthread_t t1,t2,t3;
	pthread_create(&t1,NULL,throw_points,(void*)&p1);
	pthread_create(&t2,NULL,throw_points,(void*)&p2);
	pthread_create(&t3,NULL,throw_points,(void*)&p3);
	pthread_join(t1,NULL);
	pthread_join(t2,NULL);
	pthread_join(t3,NULL);
	int Nin = p1.count+p2.count+p3.count;
	int Ntot = p1.n+p2.n+p3.n;
	double pi=4*(double)Nin/Ntot;
	printf("%i %g %g\n",Ntot,pi,fabs(pi-M_PI));
return 0;
}
