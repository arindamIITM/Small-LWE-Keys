
////////////////////////////////////// Used functions, these are defined in the bottom of this code, after the main function.
long long int choose(long long int n, long long int k);

long long int projection(long long int a, long long int n, long long int w);

void projection_ternary(long long int a);

void brent(long long int x0, long long int (*f)(long long int), long long int &ptr_x1, long long int &ptr_x2,
           long long int &ptr_count);

void floyd(long long int x0, long long int (*f)(long long int), long long int &ptr_x1, long long int &ptr_x2,
           long long int &ptr_count);


long long int calc_number(long long int *out, long long int domain);

long long int flavour(long long int a);

long long int Previous_Prime(long long int n);

long long int collision_func(long long int a);


#include <iostream>
#include <cmath>
#include <functional>

using namespace std;

////////////////////////////////////////

const long long int n = 13;

const long long int w = 2;
const int w_s = 2, q = 16;


// Function multiplication is used to calculate the product of a matrix A together with a column vector s and output will be stored in array b_1. 

long long int **binom, A[n][n], s[n], b[n], b_1[n];

void multiplication(long long int *s);


long long int ell, bit_q, f1, f2, modul_p, domain;


//Inputs of merge are two integers & output is a vector of length n which is stored in array A1
long long int A1[n], A2[n - w];

void merge(long long int v1, long long int v2);

void calc_lwe_func(long long int (*A)[n], long long int *s, long long int flag);


//We used B1 in function projection_ternary
long long int B1 = choose(n, w);



/////////////////////////////////////////////////////////////////////////////////////// the main function starts here.


int main() {

    long long int i, j, k, x, y, lam, x0, flag, number_collision_per_lwe = 0;


    double average, Sample_Points[10000];

////////////////////////////////////// we are precalculating the binomial coefficients here.

    binom = new long long int *[n + 1];
    for (long long int i = 0; i < n + 1; ++i) {
        binom[i] = new long long int[w + 1];
        for (long long int j = 0; j < w + 1; ++j) {
            binom[i][j] = choose(i, j);
        }

    }

////////////////this is the domain size
    domain = binom[n][w] * binom[n - w][w];


/////////////////////////// this prime number is used for flavouring and is the just the prime smaller than domain size

    modul_p = Previous_Prime(domain);


    float ell_before = std::log(domain);
    ell_before = ell_before / std::log(q);

////////////////////////////// this is the ell that we used in the paper also where we know the error coordinates
    ell = round(ell_before);

    bit_q = round(std::log(q) / std::log(2));

//////////////////////////// for random number generator
    srand(time(0));


    long long int collision = 0;

/////////////////////////////////////// the for loop below gives random LWE instances, and for each instance we will try to find collisions.
    for (int num_exp = 0; num_exp < 1; num_exp++) {


        /////////////// here we are generating the random martix A of order nxn
        for (i = 0; i < n; i++)
            for (j = 0; j < n; j++)
                A[i][j] = rand() % q;


        //////////////// here we are generating the ternary secret s of length n containing w_s \pm 1, rest zero
        while (1) {

            for (i = 0; i < n; i++) {
                s[i] = (rand() % 3) - 1;

            }

            k = 0;
            for (j = 0; j < n; j++)
                if (s[j] == -1)
                    k++;
            if (k != w_s)
                continue;

            k = 0;
            for (j = 0; j < n; j++)
                if (s[j] == 1)
                    k++;
            if (k != w_s)
                continue;

            k = 0;
            for (j = 0; j < n; j++)
                if (s[j] == 0)
                    k++;
            if (k != (n - 2 * w_s))
                continue;

            break;
        }

        /*we are calculating the column vector b here. multiplication function is void function, it always stores the column vector in b_1[] array.*/

        multiplication(s);

        for (i = 0; i < n; i++)
            b[i] = b_1[i];

        long long int iteration, count = 0;

        long long int u1[n], found = 0;


        number_collision_per_lwe = 0;
        collision = 0;

        // here for each instance of lwe we are calculating how many collisions are needed to find solution

        for (iteration = 0; iteration < 1000000000; iteration++) { //

            // generating f1, f2 for flavours and x0 for collision search
            f1 = 1 + rand() % (modul_p - 2);
            f2 = 1 + rand() % (modul_p - 2);
            x0 = 1 + rand() % (domain - 2);


            // collision search
            floyd(x0, collision_func, x, y, lam);



            // checking if the collision is for function with same switch. this count should be half of the total iterations.

            if (flavour(x) % 2 == flavour(y) % 2) {
                count += 1;

            } else {
                // here we project the colliding points using ternary projection function and the add them coordinate wise

                projection_ternary(flavour(x));
                for (i = 0; i < n; i++)
                    u1[i] = A1[i];
                projection_ternary(flavour(y));


                for (i = 0; i < n; i++)
                    u1[i] = u1[i] + A1[i];

                // calculating A.u1 vector
                multiplication(u1);

                // matching A.u1 with the top 'ell' coordinates of b
                flag = 1;
                for (i = 0; i < ell; i++)
                    if (b_1[i] != b[i]) {
                        flag = 0;
                        break;

                    }

                if (flag == 1) {
                    Sample_Points[found] = lam;
                    found += 1;
                }


                if (found == 10000)
                    break;


            }


        }//iteration





    }
    cout << endl;
    cout << endl;

    cout << "All data points: " << endl;
    average = 0;
    for (i = 0; i < 10000; i++) {
        cout << Sample_Points[i] << ", ";
        average += Sample_Points[i];
    }

    cout << endl;
    cout << endl;

    cout << "Final average: " << (double) average / 10000 << endl;


}


//////////////////////////////////// the main function ends here.

/////////////////////////////////



void merge(long long int v1, long long int v2) {


    long long int i, j;


    for (i = 0; i < n; i++) {
        A1[i] = v1 & 1;
        v1 = v1 >> 1;

    }


    for (i = 0; i < n - w; i++) {
        A2[i] = v2 & 1;
        v2 = v2 >> 1;


    }


    j = 0;
    for (i = 0; i < n; i++)
        if (A1[i] == 0) {
            if (A2[j] == 1)
                A1[i] = -1;
            j++;
        }



//return A1;

}


//////////////////////////////////////////////

/////////////////////////////////////////

long long int choose(long long int n, long long int k) {
    if (k == 0)
        return 1;
    return (n * choose(n - 1, k - 1)) / k;
}

/////////////////////////////////////////

////////////////////////////////////////

long long int projection(long long int a, long long int n, long long int w) {
    long long int wn = n;
    long long int wk = w;
    long long int v = 0;
    long long int set = 0;
    while (wn != 0) {
        if (set == w)
            break;
        else if (wn + set == w) {
            v += (1ULL << (wn - 1));
            wn -= 1;
            set += 1;
        } else if (a < binom[wn - 1][wk])
            wn -= 1;
        else {
            a -= binom[wn - 1][wk];
            v += (1ULL << (wn - 1));
            wn -= 1;
            wk -= 1;
            set += 1;
        }
    }

    return v;
}

/////////////////////////////////////////

////////////////////////////////////////

void projection_ternary(long long int a) {


    long long int i1 = a % B1;
    long long int i2 = ((a - i1) / B1);
    long long int v1 = projection(i1, n, w);
    long long int v2 = projection(i2, n - w, w);
    merge(v1, v2);


}

/////////////////////////////////////////

////////////////////////////////////////


void
floyd(long long int x0, long long int (*collision_func)(long long int), long long int &ptr_x1, long long int &ptr_x2,
      long long int &ptr_count) {

    long long int tortoise = collision_func(x0);
    long long int hare = collision_func(collision_func(x0));
    long long int count = 0;
    while (tortoise != hare) {
        tortoise = collision_func(tortoise);
        hare = collision_func(collision_func(hare));
        count += 1;
    }
    long long int mu = 0;
    tortoise = x0;
    long long int x1 = tortoise;
    long long int x2 = hare;

    while (tortoise != hare) {
        x1 = tortoise;
        x2 = hare;
        tortoise = collision_func(tortoise);
        hare = collision_func(hare);
        mu += 1;
    }

    ptr_x1 = x1;
    ptr_x2 = x2;
    ptr_count = count;


}


/////////////////////////////////////////

////////////////////////////////////////



// Brent's algorithm for finding a collision in a  function
void
brent(long long int x0, long long int (*collision_func)(long long int), long long int &ptr_x1, long long int &ptr_x2,
      long long int &ptr_count) {
    // Set up the initial parameters
    long long int tortoise = x0;
    long long int hare = collision_func(x0);

    long long int power = 1;
    long long int lam = 1;
    long long int count = 1;

    // Search for a collision
    while (tortoise != hare) {
        if (power == lam) {
            tortoise = hare;
            power = power * 2;
            lam = 0;
        }
        hare = collision_func(hare);
        lam += 1;
        count += 1;


    }

    // Find the collision polong long int by moving the tortoise and hare at the same speed
    long long int mu = 0;
    tortoise = x0;
    hare = x0;
    long long int x1 = tortoise;
    long long int x2 = hare;
    for (long long int i = 0; i < lam; i++) {
        hare = collision_func(hare);
    }
    while (tortoise != hare) {
        x1 = tortoise;
        x2 = hare;
        tortoise = collision_func(tortoise);
        hare = collision_func(hare);
        mu += 1;
    }

    ptr_x1 = x1;
    ptr_x2 = x2;
    ptr_count = count;

}

////////////////////////////////////////

////////////////////////////////////////

void multiplication(long long int *s) {
    long long int i, j, k;
    for (i = 0; i < n; i++) {
        k = 0;
        for (j = 0; j < n; j++)
            k += s[j] * A[i][j];

        b_1[i] = (k + n * q) % q;
    }


}

////////////////////////////////////////

////////////////////////////////////////

void calc_lwe_func(long long int *s, long long int flag) {

    long long int i;
    //calculates A*s or b-A*s for given vector s depending on switch
    multiplication(s);
    if (flag == 1)
        for (i = 0; i < n; i++)
            b_1[i] = (b[i] - b_1[i] + q) % q;


}

////////////////////////////////////////

////////////////////////////////////////

long long int calc_number(long long int *out) {
    long long int a = out[0], i;
    for (i = 1; i < ell; i++)
        a += out[i] << (i * bit_q);


    return a % domain;
}



////////////////////////////////////////

////////////////////////////////////////

long long int flavour(long long int a) {


    return (f1 * a + f2 + f1 * modul_p) % modul_p;
}


////////////////////////////////////////

////////////////////////////////////////

long long int collision_func(long long int a) {


    a = flavour(a);


    projection_ternary(a);


    calc_lwe_func(A1, a % 2);
    a = calc_number(b_1);
    return a;
}


////////////////////////////////////////

////////////////////////////////////////


long long int Previous_Prime(long long int n) {

    int flag;
    while (1) {
        flag = 1;
        n = n - 1;
        // Check from 2 to n-1
        for (long long int i = 2; i < n; i++)
            if (n % i == 0) {
                flag = 0;
                break;
            }
        if (flag == 1)
            return (n);
    }


}
