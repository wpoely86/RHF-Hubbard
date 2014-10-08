/**
 * @mainpage 
 * This is a implementation of ThING using Libint for the electron repulsion integrals, and using faster
 * methods for the calculation of the S, T and U matrices.
 * @author Brecht Verstichel, Sebastian Wouters, Ward Poelmans
 * @date 08-06-2012
 */

#include <iostream>
#include <sstream>
#include <cmath>
#include <getopt.h>
#include <signal.h>

// if set, the signal has been given to stop the calculation and write current step to file
sig_atomic_t stopping = 0;

void stopcalcsignal(int sig);

//includes all important headers and defines which conditions are
//going to be used:
#include "include.h"

double Calculate_RHF(RSPM &);

/**
 * Here all the libraries come together to do the actual calculations.
 */
int main(int argc, char **argv)
{
   using std::cout;
   using std::endl;

   //initialize the random nr generator
   srand(time(NULL));

   cout.precision(10);


   int M = 10; //dim of spin restricted sp hilbert space
   int N = 10; //nr of particles
   double U = 10; // on-site interaction
   double step = 0.1;
   double shift = 0;
   bool UseMomSpace = true; // the default

   struct option long_options[] =
   {
      {"particles",  required_argument, 0, 'n'},
      {"sites",  required_argument, 0, 'm'},
      {"interaction", required_argument, 0, 'U'},
      {"step", required_argument, 0, 's'},
      {"site-space", no_argument, 0, 'b'},
      {"help",  no_argument, 0, 'h'},
      {0, 0, 0, 0}
   };

   int i,j;
   while( (j = getopt_long (argc, argv, "hn:m:U:s:b", long_options, &i)) != -1)
      switch(j)
      {
         case 'h':
         case '?':
            cout << "Usage: " << argv[0] << " [OPTIONS]\n"
               "\n"
               "    -n, --particles=particles    Set the number of particles\n"
               "    -m, --sites=sites            Set the number of sites\n"
               "    -U, --interaction=U          Set the interaction strength\n"
               "    -s, --step=s                 Set the step size. Will go from 0 to U in s steps\n"
               "    -b, --site-space             Use the site basis space instead of momentum space\n"
               "    -h, --help                   Display this help\n"
               "\n";
            return 0;
            break;
         case 'n':
            N = atoi(optarg);
            if( N <= 0)
            {
               std::cerr << "Invalid particle number!" << endl;
               return -1;
            }
            break;
         case 'm':
            M = atoi(optarg);
            if( M <= 0)
            {
               std::cerr << "Invalid particle number!" << endl;
               return -2;
            }
            break;
         case 'U':
            U = atof(optarg);
            break;
         case 's':
            step = atof(optarg);
            break;
         case 'b':
            UseMomSpace = false;
            break;
      }

   cout << "Starting with M=" << M << " N=" << N << " U=" << U << endl;

   Hubbard::init(U,M);

   if(UseMomSpace)
   {
      cout << "Working in momentum space" << endl;
      Hubbard::UseMomSpace();
   } else
      cout << "Working in site space" << endl;

   RSPM::init(M,N,shift);

   //initialize spm
   RSPM spm;
   spm.unit();

   // set up everything to handle SIGALRM
   struct sigaction act;
   act.sa_flags = 0;
   act.sa_handler = &stopcalcsignal;

   sigset_t blockset;
   sigemptyset(&blockset); // we don't block anything in the handler
   act.sa_mask = blockset;

   sigaction(SIGUSR1, &act, 0);

   double cur_U = 0;
   int runs = 0;
   while(cur_U < U)
   {
       cout << "Running with U = " << cur_U << endl;
       Hubbard::setU(cur_U);
       Calculate_RHF(spm);
       cout << endl;

       cur_U += step;
       RSPM::setshift(0);

       runs++;
       stopping = 0;
   }

   // actual calculation with correct U
   cout << "Final calculation" << endl;
   cout << "Running with U = " << U << endl;
   Hubbard::setU(U);
   Calculate_RHF(spm);

   std::stringstream h5_name;

   if(getenv("SAVE_H5_FILE"))
      h5_name << getenv("SAVE_H5_FILE");
   else
      h5_name << "result.h5";

   spm.SaveToFile(h5_name.str());

   cout << endl;
   cout << "single-particle spectrum is:" << endl;
   cout << endl;

   RSPM F;
   Vector v(M);

   RSPM::setshift(0);

   F.construct_Fock(spm);
   v.diagonalize(F);

   cout << v << endl;

   F.construct_sp_ham();

   double energy = 2.0*F.ddot(spm);

   for(int n = 0;n < N/2;++n)
      energy += 2.0*v[n];

   cout << "Ground-state energy:" << endl;
   cout << endl;

   cout << "E_0 = \t" << energy/2.0 << endl;

   cout << v << endl;

   v.diagonalize(spm);

   cout << endl;
   cout << v << endl;

   return 0;
}

double Calculate_RHF(RSPM &spm)
{
    const int M = spm.gM();
    const int N = spm.gN();
   //some variables needed in scf loop
   RSPM copy,F;
   Vector v(M);

   DIIS diis;

   RSPM commutator;

   double convergence = 1.0;

   int iter = 0;

   double Eprev = 0;

   //basic scf loop, easy peasy
   while(convergence > 1.0e-12)
   {
      ++iter;

      //construct the fock matrix
      F.construct_Fock(spm);

      if(!Hubbard::WhichSpace())
      {
         //add it to the diis object
         diis.push_F(F);

         //calculate the commutator between F and the spm
         commutator.commute(F,spm);

         //add it to diis object
         diis.push_comm(commutator);

         //construct the B matrix and solve the system: output = b coefficients
         diis.construct();

         //now construct the "relaxed" Fock matrix
         F.relax(diis);
      }

      v.diagonalize(F);

      copy = spm;

      spm.update(F);

      copy -= spm;

      convergence = std::sqrt(copy.ddot(copy));

      RSPM test;
      test.construct_Fock(spm);
      v.diagonalize(test);
      test.construct_sp_ham();

      double en= 2.0*test.ddot(spm);

      for(int n = 0;n < N/2;++n)
          en += 2.0*v[n];

      std::cout << Hubbard::getU() << "\t" << iter << "\t" << convergence << "\t" << en/2.0 << "\t" << Eprev-en << std::endl;
      Eprev = en;

      if(stopping && fabs(convergence) < 1e-3)
         break;
   }

   std::cout << "Current (convergened) point: " << Hubbard::getU() << "\t" << Eprev/2.0 << std::endl;

   return Eprev;
}

void stopcalcsignal(int sig)
{
   stopping = 1;
}

/* vim: set ts=3 sw=3 expandtab :*/
