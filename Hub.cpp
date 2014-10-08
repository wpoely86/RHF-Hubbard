#include <cmath>

#include "include.h"

double Hubbard::U = 0;
int Hubbard::M = 0;
bool Hubbard::momspace = false;

void Hubbard::init(double U, int M)
{
    Hubbard::U = U;
    Hubbard::M = M;
}

void Hubbard::clean()
{
}

double Hubbard::getU()
{
    return U;
}

void Hubbard::setU(double U)
{
    Hubbard::U = U;
}

/**
 * Return the matrix element for the one body term in the
 * 1D Hubbard hamiltonian: <a|T|b> where T the hopping term:
 * T = \sum_ij \sum_sigma a^+_i\sigma a_j\sigma , where i and j
 * are neighbouring sites and \sigma is the spin.
 * @param a first sp index
 * @param b second sp index
 * @return the matrix element
 */
double Hubbard::gT(int a, int b)
{
    if(momspace)
    {
        if(a == b)
            // in the momentum space, this term is diagonal
            return -2*cos(2.0*M_PI/M*a);
    }
    else
    {
        // the sp indices work as following:
        // a < M : spin up
        // a >= M : spin down
        // So a and a+M are the same sites but with
        // different spin

        // PBC: jump at the right end of the chain
        if(a == M-1 && b == 0)
            return -1;

        // PBC: jump at the left end of the chain
        if(b == M-1 && a == 0)
            return -1;

        // jump to right or the left
        if( (b < M-1 && a == b+1) || (b > 0 && a==b-1) )
            return -1;
    }

    return 0;
}

/**
 * The on-site interaction term of the 1D Hubbard model.
 * Gives: <ab|V|cd> with V = U *\sum_i n_i,up n_i,down
 * @param a the first sp index
 * @param a the second sp index
 * @param a the thirth sp index
 * @param a the fourth sp index
 * @return the matrix element <ab|V|cd>
 */
double Hubbard::gV(int a, int b, int c, int d)
{
    if(momspace)
    {
        // we need momentum convervation!
        // As we do RHF, all sites are doubled occupied.
        if( (a+b)%M == (c+d)%M )
            return Hubbard::U/M;
    } else
    {
        // the interaction is diagonal in the site basis
        // As we do RHF, all sites are doubled occupied.
        if( c==d && a==c && b==d )
            return Hubbard::U;
    }

    return 0;
}


void Hubbard::UseMomSpace()
{
    Hubbard::momspace = true;
}

void Hubbard::UseSitesSpace()
{
    Hubbard::momspace = false;
}

/**
 * @return true => momentum space, false => sites space
 */
bool Hubbard::WhichSpace()
{
    return Hubbard::momspace;
}
