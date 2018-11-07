#include <fcntl.h>
#include <ftw.h>
#include <errno.h>
#include <dirent.h>
#include <string.h>
#include <unistd.h>

#include "spawner.c"
#include <stdlib.h>     /* srand, rand */
#include <time.h>       /* time */

static pthread_mutex_t fork_mutex = PTHREAD_MUTEX_INITIALIZER;

/*** Launching simulations & LogLikelihood  ***/
/*--------------------------------------------*/
/*  0. initialize stuffs for torc
 1. create temporary folders
 2. fork process, where each fork does:
 2.1 enter tmp folder
 2.2 copy necessary stuff
 2.3 write input parametes to the simulation's input file
 2.4 write input parametes to the input file for likelihood evaluation
 2.5 run simulation script
 3. wait for simulation to finish
 4. read likelihood value from the file & store output in *res
 5. remove temporary direcotries
 --------------------------------------
 PARAMETERS :
 x[0] = Dw
 x[1] = rho
 x[2] = c ( Tend^2 = c / rho )
 x[3] = icx
 x[4] = icy
 x[5] = icz
 x[6] = PETsigma      
 x[7] = PETscale
 x[8] = ucT1
 x[9] = ucT2
 x[10] = Tisigma2
 */

#define USE_SCRATCH 0
#define REMOVE_DIRS 1
#define RESTRICT_WAITING_TIME 1 

void taskfun(double /*const*/ *x, int *pN, double *res, int winfo[4])
{
    // 0.
    int n  = *pN;
    int me = getpid();  /* spanwer_id : worker_id */
    char line[1024];
    char *largv[2];
    char taskname[256];
    double t0, t1; 
    int gen, chain, step, task;
    gen = winfo[0]; chain = winfo[1]; step = winfo[2]; task = winfo[3];
    

    // 1.
    //sprintf(taskname, "tmpdir.%d.%d", getpid(), torc_worker_id());
    char cwd[256], bindir[256];
    getcwd(cwd, 256);
    sprintf(bindir, "%s/TMPTMP", cwd);

#if USE_SCRATCH
    sprintf(taskname, "/scratch/tmpdir.%d.%d.%d.%d", gen, chain, step, task);
#else
    sprintf(taskname, "tmpdir.%d.%d.%d.%d", gen, chain, step, task);
#endif
    mkdir(taskname, S_IRWXU);

#if 1
    t0 = my_gettime();
    while (pthread_mutex_trylock(&fork_mutex) == EBUSY) usleep(500*1000);
#endif    
 
    // 2.
 #if 1 
    int rf = fork();
    if (rf < 0)
    {
        printf("spanwer(%d): fork failed!!!!\n", me);
        perror("fork");
        fflush(0);
    }
#else
    int tries = 0;
retry:
    rf = fork();
    if (rf < 0) {
        printf("spanwer(%d): fork failed!!!! [try = %d]\n", me, tries); fflush(0);
        tries++;
        if (tries < 10) {
            sleep(1);
            goto retry;
        }
    }
#endif
    
    
    
    if (rf == 0)
    {
        // 2.1 enter tmp folder
        chdir(taskname);
        
     // 2.2 copy necessary stuff
     //  int status = copy_from_dir("../TMPTMP");
     int status = copy_from_dir(bindir);	
     if (status!=0)
        {
            printf("Error in copy from dir\n");
            abort();
        }
        
        
        // Convert from log- to real space
        double param[n];
        param[0] = exp( x[0] );         // D
        param[1] = exp( x[1] );         // rho
        param[2] = exp( x[2] );         // c -> Tend^2 = c/rho
        param[3] = exp( x[3] );         // icx
        param[4] = exp( x[4] );         // icy
        param[5] = exp( x[5] );         // icz
        param[6] = exp( x[6] );         // PETsigma
        param[7] = exp( x[7] );         // PET scale
        param[8] = exp( x[8] );         // ucT1
        param[9] = exp( x[9] );         // ucT2
        param[10] = exp( x[10] );       // Ti_sigma2
        param[6] = param[6] * param[6];  // PETsigma2

	// account for parameter correlation
	param[2] = param[2] / param[1];  // Tend^2 = c / rho
        param[2] = sqrt(param[2]);


	// 1.3 Simulation input parametes
        int i;
	FILE *finp = fopen("InputParameters.txt", "w");
        for (i = 0; i < 3; i++) fprintf(finp, "%.16lf\n", param[i]);
        fclose(finp);
        
        // 1.4 Tumor initial position
        finp = fopen ("TumorIC.txt", "w");
	for (i = 3; i < 6; i++) fprintf(finp, "%.4f \n", param[i]);
        fclose (finp);
        
        // 1.5 write input parametes to the input file for likelihood evaluation
        // PETsigma2, PETscale, TIsigma2
        finp = fopen("LikelihoodInput.txt", "w");
        for (i = 6; i < n; i++) fprintf(finp, "%.16lf\n", param[i]);
        fclose(finp);
        
        /* 2. run simulation */
        sprintf(line, "./runAll.sh");
        parse(line, largv);
        
#if 1
        int fd = open("output", O_RDWR | O_CREAT, S_IRUSR | S_IWUSR);
        
        dup2(fd, 1);    // make stdout go to file
        dup2(fd, 2);    // make stderr go to file
        close(fd);	// fd no longer needed
#endif
        
        status = execvp(*largv, largv);
    }

    pthread_mutex_unlock(&fork_mutex);
    
    // 3.   
    if(RESTRICT_WAITING_TIME){
        if (rf > 0) {
            int status, wpid, ctr = 0;
            while ((wpid = waitpid(rf, &status, WNOHANG)) == 0) {
                 ctr++;
                 sleep(1);
                if (ctr == 1100) {
                    printf("XXXX: KILLING %d\n", rf); fflush(0);
                    kill(rf, SIGKILL);
                    break;
                }
	    	}
         }
	}
   else {
    int status;
    waitpid(rf, &status, 0);
   }
 
    // 4.
    double likelihood;
    
    FILE *pFile;
    // double gval = -1, gres;
    char fitname[256];
    
    //sprintf(fitname, "./%s/Likelihood.txt",taskname);
    sprintf(fitname, "%s/Likelihood.txt",taskname);  
    pFile = fopen(fitname, "r");
    
    if (pFile == NULL)
    {
        printf("spawner(%d): error with fopen for %s\n", me, fitname); fflush(0);
        //gres = -1;
        likelihood = 1e12; // penalty
    }
    else
    {
        while (!feof(pFile))
            fscanf(pFile, "%lf", &likelihood );
        
        fclose (pFile); 
    
    	if (isnan(likelihood) || isinf(likelihood))
    	likelihood = 1.0e+12; // penalty
    }
    
    likelihood = -likelihood;
    *res = likelihood;
    
    //inc_nfc();  /* increment function call counter*/
    
    // 5.
    if (REMOVE_DIRS) 
       rmrf(taskname);
    
    t1 = my_gettime();
    {
     	int i;
    	char msg[1024], buf[64];
    	sprintf(msg, "task(%d.%d.%d.%d):", gen, chain, step, task);
    	for (i = 0; i < n; i++)
    	{
     		sprintf(buf, "%.6lf ", x[i]);
        	strcat(msg, buf);
    	}
   
	 sprintf(buf, " = %.6lf", likelihood);
    	strcat(msg, buf);

    	sprintf(buf, " in %lf secs\n", t1-t0);
    	strcat(msg, buf);

    	printf("%s", msg); fflush(0);
    }

    return;
}

