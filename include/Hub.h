#ifndef HUB_H
#define HUB_H

class Hubbard
{
    public:
        static void init(double U,int M);
        static void clean();

        static void setU(double U);

        static double gT(int a, int b);
        static double gV(int a, int b, int c, int d);

        static double getU();

        static void UseMomSpace();

        static void UseSitesSpace();

        static bool WhichSpace();

    private:
        static double U;
        static int M;

        static bool momspace;

};

#endif /* HUB_H */
