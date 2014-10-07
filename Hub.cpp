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

double Hubbard::gT(int a, int b)
{
    if(momspace)
    {
        if(a == b)
            return -2*cos(2.0*M_PI/M*a);
    }
    else
    {
        if(a == M-1 && b == 0)
            return -1;

        if(b == M-1 && a == 0)
            return -1;

        if( (b < M-1 && a == b+1) || (b > 0 && a==b-1) )
            return -1;
    }

    return 0;
}

double Hubbard::gV(int a, int b, int c, int d)
{
    if(momspace)
    {
        if( (a+b)%M == (c+d)%M )
            return Hubbard::U/M;
    } else
    {
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
