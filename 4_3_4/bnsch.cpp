#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <time.h>
#include "util.h"
#include "bnsch.h"

#define C 30.0
#define N 2

int nx, ny, nz, n_level, c_relax;
float ***ct, ***cd, ***sc, ***smu, ***sor, ***mu, h, h2, h3, dt,
    xleft, xright, yleft, yright, zleft, zright, pi, gam, Cahn, sig, V0, V1, V2,
    betac, betad, M, M2, kappa, BBC, BBD, sum_c1, sum_c2, aphi, theta, s_har;

int main()
{
    extern int nx, ny, nz, n_level, c_relax;
    extern float ***ct, ***cd, ***sc, ***smu, ***sor, ***mu, h, h2, h3, dt,
        xright, yright, yleft, yright, zleft, zright, pi, gam, Cahn,
        betac, betad, M2, kappa, BBC, BBD, sum_c1, sum_c2, sig, V0, V1, V2;

    char buffer[20];
    int i, j, k, ns, it, max_it, count = 1;
    float ***oc, ***ooc, ***nc, value1_1, value2_1, sum_d1, diff, tol;
    float ***od, ***ood, ***nd, value1_2, value2_2, sum_d2;

    FILE *ddd, *fp;

    char buf_energy1[20], buf_area1[20], buf_volume1[20], buf_error1[20];
    char buf_energy2[20], buf_area2[20], buf_volume2[20], buf_error2[20];

    float surface_areac, energyc;
    float surface_aread, energyd;

    clock_t start;
    float elapsed;
    float max_error1, max_error2, temp;

    c_relax = 5;

    nx = gnx;
    ny = gny;
    nz = gnz;

    pi = 4.0 * atan(1.0);
    n_level = (int)(log(ny) / log(2) - 0.9);

    xleft = 0.0, xright = 4.0;
    yleft = 0.0, yright = xright * ny / (1.0 * nx);
    zleft = 0.0, zright = xright * nz / (1.0 * nx);

    max_it = 10000;

    ns = max_it / 1000;

    h = xright / (float)nx;
    h2 = pow(h, 2);
    h3 = pow(h, 3);

    gam = 3.0 * h / (2.0 * sqrt(2.0) * atanh(0.9));

    Cahn = pow(gam, 2);
    dt = 2.0e-5;

    tol = 1.0e-6;

    M = 300.0;
    kappa = 1.0;

    M2 = M * 3.0 / sqrt(2.0);

    sig = 1.0;

    printf("nx = %d , ny = %d , nz = %d\n", nx, ny, nz);
    printf("dt      = %f\n", dt);
    printf("max_it  = %d\n", max_it);
    printf("ns      = %d\n", ns);
    printf("n_level           = %d\n\n", n_level);

    ddd = fopen("data1/remarks.m", "w");
    fprintf(ddd, "nx = %d; ny = %d; nz = %d; \n", nx, ny, nz);
    fprintf(ddd, "xright  = %f; yright = %f; zright = %f; \n", xright, yright, zright);
    fprintf(ddd, "dt      = %f; \n", dt);
    fprintf(ddd, "compare_dt     = %f; \n", dt * kappa);
    fprintf(ddd, "gam     = %f; \n", gam);
    fprintf(ddd, "max_it  = %d; \n", max_it);
    fprintf(ddd, "ns      = %d; \n", ns);
    fprintf(ddd, "n_level = %d; \n\n", n_level);
    fprintf(ddd, "C = %f\n", float(C));
    fprintf(ddd, "data index: 1\n\n");
    fclose(ddd);

    ddd = fopen("data2/remarks.m", "w");
    fprintf(ddd, "nx = %d; ny = %d; nz = %d; \n", nx, ny, nz);
    fprintf(ddd, "xright  = %f; yright = %f; zright = %f; \n", xright, yright, zright);
    fprintf(ddd, "dt      = %f; \n", dt);
    fprintf(ddd, "compare_dt     = %f; \n", dt * kappa);
    fprintf(ddd, "gam     = %f; \n", gam);
    fprintf(ddd, "max_it  = %d; \n", max_it);
    fprintf(ddd, "ns      = %d; \n", ns);
    fprintf(ddd, "n_level = %d; \n\n", n_level);
    fprintf(ddd, "C = %f\n", float(C));
    fprintf(ddd, "data index: 2\n\n");
    fclose(ddd);

    oc = cube(0, nx + 1, 0, ny + 1, 0, nz + 1);
    ooc = cube(0, nx + 1, 0, ny + 1, 0, nz + 1);
    nc = cube(0, nx + 1, 0, ny + 1, 0, nz + 1);
    ct = cube(1, nx, 1, ny, 1, nz);

    od = cube(0, nx + 1, 0, ny + 1, 0, nz + 1);
    ood = cube(0, nx + 1, 0, ny + 1, 0, nz + 1);
    nd = cube(0, nx + 1, 0, ny + 1, 0, nz + 1);
    cd = cube(1, nx, 1, ny, 1, nz);

    mu = cube(0, nx + 1, 0, ny + 1, 0, nz + 1);

    sc = cube(1, nx, 1, ny, 1, nz);
    smu = cube(1, nx, 1, ny, 1, nz);

    sprintf(buf_energy1, "energy1.m");
    sprintf(buf_area1, "surface_areac1.m");
    sprintf(buf_volume1, "volume1.m");
    sprintf(buf_error1, "outerror1.m");

    sprintf(buf_energy2, "energy2.m");
    sprintf(buf_area2, "surface_areac2.m");
    sprintf(buf_volume2, "volume2.m");
    sprintf(buf_error2, "outerror2.m");

    fp = fopen(buf_energy1, "w");
    fclose(fp);
    fp = fopen(buf_area1, "w");
    fclose(fp);
    fp = fopen(buf_volume1, "w");
    fclose(fp);
    fp = fopen(buf_error1, "w");
    fclose(fp);

    fp = fopen(buf_energy2, "w");
    fclose(fp);
    fp = fopen(buf_area2, "w");
    fclose(fp);
    fp = fopen(buf_volume2, "w");
    fclose(fp);
    fp = fopen(buf_error2, "w");
    fclose(fp);

    zero_cube(mu, 1, nx, 1, ny, 1, nz);

    initialization(oc, 1);

    print_data(oc, count, 1);

    initialization(od, 2);
    print_data(od, count, 2);

    V1 = 0.0;
    V2 = 0.0;
    ijkloop
    {

        V1 = V1 + 3.0 * pow(0.5 * (oc[i][j][k] + 1.0), 2) - 2.0 * pow(0.5 * (oc[i][j][k] + 1.0), 3);
        V2 = V2 + 3.0 * pow(0.5 * (od[i][j][k] + 1.0), 2) - 2.0 * pow(0.5 * (od[i][j][k] + 1.0), 3);
    }

    V1 = V1 * h * h * h;
    V2 = V2 * h * h * h;

    float initial_massc = 0.0;
    ijkloop
    {
        initial_massc += oc[i][j][k];
    }

    float initial_massd = 0.0;
    ijkloop
    {
        initial_massd += od[i][j][k];
    }

    betac = Bphi(oc);
    betad = Bphi(od);

    start = clock();
    for (it = 1; it <= max_it; it++)
    {
        printf("iter = %d \n", it);

        cube_copy(ooc, oc, 1, nx, 1, ny, 1, nz);
        cube_copy(ood, od, 1, nx, 1, ny, 1, nz);

        surface_areac = Bphi(oc);
        surface_aread = Bphi(od);

        BBC = surface_areac - betac;
        BBD = surface_aread - betad;

        fp = fopen(buf_area1, "a");
        fprintf(fp, "%f \n", surface_areac);
        fclose(fp);

        fp = fopen(buf_area2, "a");
        fprintf(fp, "%f \n", surface_aread);
        fclose(fp);

        energyc = 0.0;
        energyd = 0.0;

        ijkloop
        {
            energyc += mu[i][j][k] * mu[i][j][k];
        }
        ijkloop
        {
            energyd += mu[i][j][k] * mu[i][j][k];
        }

        fp = fopen(buf_energy1, "a");

        energyc *= kappa * 3.0 / (2.0 * sqrt(2.0) * pow(gam, 3)) * h3;
        fprintf(fp, "%f ", energyc);

        energyc += M * BBC * BBC;
        fprintf(fp, "%f ", M * BBC * BBC);

        fprintf(fp, "%f \n", energyc);
        fclose(fp);

        fp = fopen(buf_energy2, "a");

        energyd *= kappa * 3.0 / (2.0 * sqrt(2.0) * pow(gam, 3)) * h3;
        fprintf(fp, "%f ", energyd);

        energyd += M * BBD * BBD;
        fprintf(fp, "%f ", M * BBD * BBD);

        fprintf(fp, "%f \n", energyd);
        fclose(fp);

        V0 = V1;

        cahn(oc, nc, oc, od, BBC, ct);

        V0 = V2;

        cahn(od, nd, od, oc, BBD, cd);

        /* sum_c1 = sum_d1 = 0.0;
         ijkloop {
             sum_c1 += nc[i][j][k];
             sum_d1 += pow(1-nc[i][j][k],2)*pow(1+nc[i][j][k],2);
         }

         sum_c2 = sum_d2 = 0.0;
         ijkloop {
             sum_c2 += nd[i][j][k];
             sum_d2 += pow(1-nd[i][j][k],2)*pow(1+nd[i][j][k],2);
         }

         sum_c3 = sum_d3 = 0.0;
         ijkloop {
             sum_c3 += ne[i][j][k];
             sum_d3 += pow(1-ne[i][j][k],2)*pow(1+ne[i][j][k],2);
         }

         value2_1 = (initial_massc - sum_c1)/(sum_d1);
         value2_2 = (initial_massd - sum_c2)/(sum_d2);
         value2_3 = (initial_masse - sum_c3)/(sum_d3);

          */

        ijkloop
        {
            oc[i][j][k] = nc[i][j][k];
        }

        ijkloop
        {
            od[i][j][k] = nd[i][j][k];
        }

        max_error1 = 0.0;
        ijkloop
        {
            temp = fabs(oc[i][j][k] - ooc[i][j][k]);
            if (max_error1 < temp)
                max_error1 = temp;
        }
        max_error1 = max_error1 / dt;

        fp = fopen(buf_error1, "a");
        fprintf(fp, "%f ", max_error1);
        fclose(fp);

        max_error1 = 0.0;
        ijkloop
        {
            temp = oc[i][j][k] - ooc[i][j][k];
            max_error1 += temp * temp;
        }
        max_error1 = max_error1 * h3 / (dt * dt);

        fp = fopen(buf_error1, "a");
        fprintf(fp, "%f \n", sqrt(max_error1));
        fclose(fp);

        max_error2 = 0.0;
        ijkloop
        {
            temp = fabs(od[i][j][k] - ood[i][j][k]);
            if (max_error2 < temp)
                max_error2 = temp;
        }
        max_error2 = max_error2 / dt;

        fp = fopen(buf_error2, "a");
        fprintf(fp, "%f ", max_error2);
        fclose(fp);

        max_error2 = 0.0;
        ijkloop
        {
            temp = od[i][j][k] - ood[i][j][k];
            max_error2 += temp * temp;
        }
        max_error2 = max_error2 * h3 / (dt * dt);

        fp = fopen(buf_error2, "a");
        fprintf(fp, "%f \n", sqrt(max_error2));
        fclose(fp);

        if (max_error1 < 1e-5 && max_error2 < 1e-5)
        {
            break;
        }

        sum_c1 = 0.0;
        ijkloop
        {
            sum_c1 += oc[i][j][k];
        }
        sum_c1 *= h3;

        fp = fopen(buf_volume1, "a");
        fprintf(fp, "%f \n", sum_c1);
        fclose(fp);

        sum_c2 = 0.0;
        ijkloop
        {
            sum_c2 += od[i][j][k];
        }
        sum_c2 *= h3;

        fp = fopen(buf_volume2, "a");
        fprintf(fp, "%f \n", sum_c2);
        fclose(fp);

        if (it % ns == 0)
        {
            count++;
            print_data(oc, count, 1);
            printf("print out counts %d , time %f, surface area deference=%f \n", count - 1, it * 1.0, BBC);
            print_data(od, count, 2);
            printf("print out counts %d , time %f, surface area deference=%f \n", count - 1, it * 1.0, BBD);
        }
    }

    count++;
    print_data(oc, count, 1);
    print_data(od, count, 2);

    printf("print out counts %d , time %f, surface area deference=%f \n", count - 1, it * 1.0, BBC);
    printf("print out counts %d , time %f, surface area deference=%f \n", count - 1, it * 1.0, BBD);

    elapsed = ((double)(clock() - start)) / CLOCKS_PER_SEC;
    printf("elapsed times = %f \n", elapsed);

    return 0;
}

void initialization(float ***phi, int sphere_id)
{
    extern float gam, xright, yright, zright;
    extern float ***mu;

    int i, j, k;
    float x, y, z, r, ff;
    float rra, rrb, rrc, rrd;

    float x_centera, y_centera, z_center;
    float x_centerb, y_centerb;
    float x_centerc, y_centerc;
    float x_centerd, y_centerd;

    pi = 2.0 * asin(1.0);

    x_centera = 2.3;
    y_centera = 2.3;

    x_centerb = 2.3;
    y_centerb = 1.7;

    x_centerc = 1.7;
    y_centerc = 2.3;

    x_centerd = 1.7;
    y_centerd = 1.7;

    if (sphere_id == 1)
        z_center = 1.7;

    else if (sphere_id == 2)
        z_center = 2.3;

    ijkloop
    {

        x = ((double)i - 0.5) * h;
        y = ((double)j - 0.5) * h;
        z = ((double)k - 0.5) * h;

        r = 0.3;
        rra = sqrt(pow(x - x_centera, 2) + pow(y - y_centera, 2) + pow(z - z_center, 2));
        rrb = sqrt(pow(x - x_centerb, 2) + pow(y - y_centerb, 2) + pow(z - z_center, 2));
        rrc = sqrt(pow(x - x_centerc, 2) + pow(y - y_centerc, 2) + pow(z - z_center, 2));
        rrd = sqrt(pow(x - x_centerd, 2) + pow(y - y_centerd, 2) + pow(z - z_center, 2));

        phi[i][j][k] = tanh((r - rra) / (sqrt(2.0) * gam)) + tanh((r - rrb) / (sqrt(2.0) * gam)) + tanh((r - rrc) / (sqrt(2.0) * gam)) + tanh((r - rrd) / (sqrt(2.0) * gam)) + 3;
    }

    float ***lap_phi = cube(1, nx, 1, ny, 1, nz);

    laplace_ch(phi, lap_phi, nx, ny, nz);

    ijkloop
    {
        mu[i][j][k] = pow(phi[i][j][k], 3) - phi[i][j][k] - Cahn * lap_phi[i][j][k];
    }
    free_cube(lap_phi, 1, nx, 1, ny, 1, nz);
}

void cahn(float ***c_old, float ***c_new, float ***o1, float ***o2, float BB, float ***ct)
{
    extern int nx, ny, nz;
    extern float ***sc, ***smu, ***mu;

    int i, j, k, max_it_CH = 100, it_mg = 1;
    float tol = 1.0e-6, resid = 1.0;

    source(c_old, o1, o2, sc, smu);

    cube_copy(ct, c_old, 1, nx, 1, ny, 1, nz);

    while (it_mg <= max_it_CH && resid > tol)
    {

        vcycle(c_new, c_old, mu, sc, smu, nx, ny, nz, 1, BB);

        resid = error(ct, c_new, nx, ny, nz);
        cube_copy(ct, c_new, 1, nx, 1, ny, 1, nz);

        it_mg++;
    }

    printf("error %15.12f  %d  \n", resid, it_mg - 1);
}

void source(float ***c_old, float ***o1, float ***o2, float ***src_c, float ***src_mu)
{
    extern float dt, sig, V0;
    int i, j, k;
    float C1, C2, V;
    C1 = C2 = C;

    float ***oo1 = cube(1, nx, 1, ny, 1, nz);

    float ***oo2 = cube(1, nx, 1, ny, 1, nz);

    ijkloop
    {

        if (o1[i][j][k] > 1)
            oo1[i][j][k] = 1;
        else if (o1[i][j][k] < -1)
            oo1[i][j][k] = -1;
        else
            oo1[i][j][k] = o1[i][j][k];

        if (o2[i][j][k] > 1)
            oo2[i][j][k] = 1;
        else if (o2[i][j][k] < -1)
            oo2[i][j][k] = -1;
        else
            oo2[i][j][k] = o2[i][j][k];

        V = 0.0;

        V = V + 3.0 * pow(0.5 * (oo1[i][j][k] + 1.0), 2) - 2.0 * pow(0.5 * (oo1[i][j][k] + 1.0), 3);
    }

    V = V * h * h * h;

    ijkloop
    {

        src_c[i][j][k] = (c_old[i][j][k] / dt) + (3.0 / (sqrt(2) * gam)) * C1 * oo1[i][j][k] * (oo2[i][j][k] * oo2[i][j][k] - 1) - sig * (6.0 * (0.5 * (oo1[i][j][k] + 1.0)) - 6 * pow(0.5 * (oo1[i][j][k] + 1.0), 2)) * (V - V0);

        src_mu[i][j][k] = -c_old[i][j][k];
    }

    free_cube(oo1, 1, nx, 1, ny, 1, nz);
    free_cube(oo2, 1, nx, 1, ny, 1, nz);
}

void vcycle(float ***uf_new, float ***uf_old, float ***wf_new,
            float ***su, float ***sw, int nxf, int nyf, int nzf, int ilevel, float BB)
{
    extern int n_level;

    relax(uf_new, uf_old, wf_new, su, sw, ilevel, nxf, nyf, nzf, BB);

    if (ilevel < n_level)
    {

        int nxc, nyc, nzc;
        float ***duc, ***dwc, ***uc_new, ***wc_new, ***uc_old,
            ***uc_def, ***wc_def, ***uf_def, ***wf_def;

        nxc = nxf / 2;
        nyc = nyf / 2;
        nzc = nzf / 2;

        duc = cube(1, nxc, 1, nyc, 1, nzc);
        dwc = cube(1, nxc, 1, nyc, 1, nzc);
        uc_new = cube(0, nxc + 1, 0, nyc + 1, 0, nzc + 1);
        uc_old = cube(1, nxc, 1, nyc, 1, nzc);
        wc_new = cube(0, nxc + 1, 0, nyc + 1, 0, nzc + 1);
        uf_def = cube(1, nxf, 1, nyf, 1, nzf);
        wf_def = cube(1, nxf, 1, nyf, 1, nzf);
        uc_def = cube(0, nxc + 1, 0, nyc + 1, 0, nzc + 1);
        wc_def = cube(0, nxc + 1, 0, nyc + 1, 0, nzc + 1);

        restrict3(uf_new, uc_new, uf_old, uc_old, wf_new, wc_new, nxc, nyc, nzc);

        defect(duc, dwc, uf_new, uf_old, wf_new, su, sw, nxf, nyf, nzf, uc_new, uc_old, wc_new, nxc, nyc, nzc, BB);

        cube_copy2(uc_def, uc_new, wc_def, wc_new, 1, nxc, 1, nyc, 1, nzc);

        vcycle(uc_def, uc_old, wc_def, duc, dwc, nxc, nyc, nzc, ilevel + 1, BB);

        cube_sub2(uc_def, uc_def, uc_new, wc_def, wc_def, wc_new, 1, nxc, 1, nyc, 1, nzc);

        prolong_ch(uc_def, uf_def, wc_def, wf_def, nxc, nyc, nzc);

        cube_add2(uf_new, uf_new, uf_def, wf_new, wf_new, wf_def, 1, nxf, 1, nyf, 1, nzf);

        relax(uf_new, uf_old, wf_new, su, sw, ilevel, nxf, nyf, nzf, BB);

        free_cube(duc, 1, nxc, 1, nyc, 1, nzc);
        free_cube(dwc, 1, nxc, 1, nyc, 1, nzc);
        free_cube(uc_new, 0, nxc + 1, 0, nyc + 1, 0, nzc + 1);
        free_cube(uc_old, 1, nxc, 1, nyc, 1, nzc);
        free_cube(wc_new, 0, nxc + 1, 0, nyc + 1, 0, nzc + 1);
        free_cube(uf_def, 1, nxf, 1, nyf, 1, nzf);
        free_cube(wf_def, 1, nxf, 1, nyf, 1, nzf);
        free_cube(uc_def, 0, nxc + 1, 0, nyc + 1, 0, nzc + 1);
        free_cube(wc_def, 0, nxc + 1, 0, nyc + 1, 0, nzc + 1);
    }
}

void relax(float ***c_new, float ***c_old, float ***mu_new, float ***su, float ***sw,
           int ilevel, int nxt, int nyt, int nzt, float BB)
{
    extern int c_relax;
    extern float dt, Cahn, xright, M2, gam, kappa;

    int i, j, k, iter;
    float ht2, a[4], f[2], det;
    float tmp_c, tmp_mu;
    float c[7];

    ht2 = pow(xright / (float)nxt, 2);

    c[0] = 1.0 / dt;
    c[1] = 3.0 * kappa / (sqrt(2.0) * pow(gam, 3));
    c[2] = 9.0 * sqrt(2.0) * kappa / (gam * ht2);

    c[3] = M2 / gam;
    c[4] = 3.0 * kappa / (sqrt(2.0) * gam * ht2);

    c[5] = 6.0 * pow(gam, 2) / ht2;
    c[6] = pow(gam, 2) / ht2;

    for (iter = 1; iter <= c_relax; iter++)
    {

        augmenc(c_new, nxt, nyt, nzt);
        augmenc(mu_new, nxt, nyt, nzt);

        ijkloopt
        {

            a[0] = c[0];
            a[1] = c[1] * (3.0 * pow(c_old[i][j][k], 2) - 1.0) + c[2] + c[3] * BB;
            a[2] = -3.0 * pow(c_new[i][j][k], 2) - c[5];
            a[3] = 1.0;

            tmp_c = c_new[i - 1][j][k] + c_new[i + 1][j][k] + c_new[i][j - 1][k] + c_new[i][j + 1][k] + c_new[i][j][k - 1] + c_new[i][j][k + 1];
            tmp_mu = mu_new[i - 1][j][k] + mu_new[i + 1][j][k] + mu_new[i][j - 1][k] + mu_new[i][j + 1][k] + mu_new[i][j][k - 1] + mu_new[i][j][k + 1];

            f[0] = su[i][j][k] + c[4] * tmp_mu;
            f[1] = sw[i][j][k] - 2.0 * pow(c_new[i][j][k], 3) - c[6] * tmp_c;

            det = a[0] * a[3] - a[1] * a[2];

            c_new[i][j][k] = (a[3] * f[0] - a[1] * f[1]) / det;
            mu_new[i][j][k] = (-a[2] * f[0] + a[0] * f[1]) / det;
        }
    }
}

void restrict2(float ***uf, float ***uc, float ***vf, float ***vc,
               int nxc, int nyc, int nzc)
{
    int i, j, k;

    for (i = 1; i <= nxc; i++)
        for (j = 1; j <= nyc; j++)
            for (k = 1; k <= nzc; k++)
            {

                uc[i][j][k] = 0.125 * (uf[2 * i][2 * j][2 * k] + uf[2 * i - 1][2 * j][2 * k] + uf[2 * i][2 * j - 1][2 * k] + uf[2 * i - 1][2 * j - 1][2 * k] + uf[2 * i][2 * j][2 * k - 1] + uf[2 * i - 1][2 * j][2 * k - 1] + uf[2 * i][2 * j - 1][2 * k - 1] + uf[2 * i - 1][2 * j - 1][2 * k - 1]);

                vc[i][j][k] = 0.125 * (vf[2 * i][2 * j][2 * k] + vf[2 * i - 1][2 * j][2 * k] + vf[2 * i][2 * j - 1][2 * k] + vf[2 * i - 1][2 * j - 1][2 * k] + vf[2 * i][2 * j][2 * k - 1] + vf[2 * i - 1][2 * j][2 * k - 1] + vf[2 * i][2 * j - 1][2 * k - 1] + vf[2 * i - 1][2 * j - 1][2 * k - 1]);
            }
}

void restrict3(float ***uf, float ***uc, float ***vf, float ***vc, float ***wf, float ***wc,
               int nxc, int nyc, int nzc)
{
    int i, j, k;

    for (i = 1; i <= nxc; i++)
        for (j = 1; j <= nyc; j++)
            for (k = 1; k <= nzc; k++)
            {

                uc[i][j][k] = 0.125 * (uf[2 * i][2 * j][2 * k] + uf[2 * i - 1][2 * j][2 * k] + uf[2 * i][2 * j - 1][2 * k] + uf[2 * i - 1][2 * j - 1][2 * k] + uf[2 * i][2 * j][2 * k - 1] + uf[2 * i - 1][2 * j][2 * k - 1] + uf[2 * i][2 * j - 1][2 * k - 1] + uf[2 * i - 1][2 * j - 1][2 * k - 1]);

                vc[i][j][k] = 0.125 * (vf[2 * i][2 * j][2 * k] + vf[2 * i - 1][2 * j][2 * k] + vf[2 * i][2 * j - 1][2 * k] + vf[2 * i - 1][2 * j - 1][2 * k] + vf[2 * i][2 * j][2 * k - 1] + vf[2 * i - 1][2 * j][2 * k - 1] + vf[2 * i][2 * j - 1][2 * k - 1] + vf[2 * i - 1][2 * j - 1][2 * k - 1]);
                wc[i][j][k] = 0.125 * (wf[2 * i][2 * j][2 * k] + wf[2 * i - 1][2 * j][2 * k] + wf[2 * i][2 * j - 1][2 * k] + wf[2 * i - 1][2 * j - 1][2 * k] + wf[2 * i][2 * j][2 * k - 1] + wf[2 * i - 1][2 * j][2 * k - 1] + wf[2 * i][2 * j - 1][2 * k - 1] + wf[2 * i - 1][2 * j - 1][2 * k - 1]);
            }
}

void prolong_ch(float ***uc, float ***uf, float ***vc, float ***vf,
                int nxc, int nyc, int nzc)
{
    int i, j, k;

    for (i = 1; i <= nxc; i++)
        for (j = 1; j <= nyc; j++)
            for (k = 1; k <= nzc; k++)
            {

                uf[2 * i][2 * j][2 * k] = uc[i][j][k];
                uf[2 * i - 1][2 * j][2 * k] = uc[i][j][k];
                uf[2 * i][2 * j - 1][2 * k] = uc[i][j][k];
                uf[2 * i - 1][2 * j - 1][2 * k] = uc[i][j][k];

                uf[2 * i][2 * j][2 * k - 1] = uc[i][j][k];
                uf[2 * i - 1][2 * j][2 * k - 1] = uc[i][j][k];
                uf[2 * i][2 * j - 1][2 * k - 1] = uc[i][j][k];
                uf[2 * i - 1][2 * j - 1][2 * k - 1] = uc[i][j][k];

                vf[2 * i][2 * j][2 * k] = vc[i][j][k];
                vf[2 * i - 1][2 * j][2 * k] = vc[i][j][k];
                vf[2 * i][2 * j - 1][2 * k] = vc[i][j][k];
                vf[2 * i - 1][2 * j - 1][2 * k] = vc[i][j][k];

                vf[2 * i][2 * j][2 * k - 1] = vc[i][j][k];
                vf[2 * i - 1][2 * j][2 * k - 1] = vc[i][j][k];
                vf[2 * i][2 * j - 1][2 * k - 1] = vc[i][j][k];
                vf[2 * i - 1][2 * j - 1][2 * k - 1] = vc[i][j][k];
            }
}

void defect(float ***duc, float ***dwc,
            float ***uf_new, float ***uf_old, float ***wf_new,
            float ***suf, float ***swf, int nxf, int nyf, int nzf,
            float ***uc_new, float ***uc_old, float ***wc_new,
            int nxc, int nyc, int nzc, float BB)
{

    float ***ruf, ***rwf, ***rruf, ***rrwf, ***ruc, ***rwc;

    ruc = cube(1, nxc, 1, nyc, 1, nzc);
    rwc = cube(1, nxc, 1, nyc, 1, nzc);
    ruf = cube(1, nxf, 1, nyf, 1, nzf);
    rwf = cube(1, nxf, 1, nyf, 1, nzf);
    rruf = cube(1, nxc, 1, nyc, 1, nzc);
    rrwf = cube(1, nxc, 1, nyc, 1, nzc);

    nonL(ruc, rwc, uc_new, uc_old, wc_new, nxc, nyc, nzc, BB);
    nonL(ruf, rwf, uf_new, uf_old, wf_new, nxf, nyf, nzf, BB);

    cube_sub2(ruf, suf, ruf, rwf, swf, rwf, 1, nxf, 1, nyf, 1, nzf);

    restrict2(ruf, rruf, rwf, rrwf, nxc, nyc, nzc);

    cube_add2(duc, ruc, rruf, dwc, rwc, rrwf, 1, nxc, 1, nyc, 1, nzc);

    free_cube(ruc, 1, nxc, 1, nyc, 1, nzc);
    free_cube(rwc, 1, nxc, 1, nyc, 1, nzc);
    free_cube(ruf, 1, nxf, 1, nyf, 1, nzf);
    free_cube(rwf, 1, nxf, 1, nyf, 1, nzf);
    free_cube(rruf, 1, nxc, 1, nyc, 1, nzc);
    free_cube(rrwf, 1, nxc, 1, nyc, 1, nzc);
}

void nonL(float ***ru, float ***rw, float ***c_new, float ***c_old, float ***mu_new,
          int nxt, int nyt, int nzt, float BB)
{
    extern float dt, Cahn, M2, kappa, gam;

    int i, j, k;
    float ***lap_mu, ***lap_c;
    float c[4];

    lap_mu = cube(1, nxt, 1, nyt, 1, nzt);
    lap_c = cube(1, nxt, 1, nyt, 1, nzt);

    laplace_ch(c_new, lap_c, nxt, nyt, nzt);
    laplace_ch(mu_new, lap_mu, nxt, nyt, nzt);

    c[0] = 1.0 / dt;
    c[1] = 3.0 * kappa / (sqrt(2.0) * pow(gam, 3));
    c[2] = 3.0 * kappa / (sqrt(2.0) * gam);
    c[3] = M2 / gam;

    ijkloopt
    {

        ru[i][j][k] = c[0] * c_new[i][j][k] + c[1] * (3.0 * pow(c_old[i][j][k], 2.0) - 1.0) * mu_new[i][j][k] - c[2] * lap_mu[i][j][k] + c[3] * BB * mu_new[i][j][k];
        rw[i][j][k] = mu_new[i][j][k] - pow(c_new[i][j][k], 3.0) + Cahn * lap_c[i][j][k];
    }

    free_cube(lap_mu, 1, nxt, 1, nyt, 1, nzt);
    free_cube(lap_c, 1, nxt, 1, nyt, 1, nzt);
}

void laplace_ch(float ***a, float ***lap_a, int nxt, int nyt, int nzt)
{
    extern float xright;

    int i, j, k;
    float ht2, dadx_L, dadx_R, dady_B, dady_T, dadz_D, dadz_U;

    ht2 = pow(xright / (float)nxt, 2);

    augmenc(a, nxt, nyt, nzt);

    ijkloopt
    {

        dadx_L = a[i][j][k] - a[i - 1][j][k];
        dadx_R = a[i + 1][j][k] - a[i][j][k];
        dady_B = a[i][j][k] - a[i][j - 1][k];
        dady_T = a[i][j + 1][k] - a[i][j][k];
        dadz_D = a[i][j][k] - a[i][j][k - 1];
        dadz_U = a[i][j][k + 1] - a[i][j][k];

        lap_a[i][j][k] = (dadx_R - dadx_L + dady_T - dady_B + dadz_U - dadz_D) / ht2;
    }
}

float error(float ***c_old, float ***c_new, int nxt, int nyt, int nzt)
{
    float ***r, res;

    r = cube(1, nxt, 1, nyt, 1, nzt);

    cube_sub(r, c_new, c_old, 1, nxt, 1, nyt, 1, nzt);
    res = cube_max(r, 1, nxt, 1, nyt, 1, nzt);

    free_cube(r, 1, nxt, 1, nyt, 1, nzt);

    return res;
}

void augmenc(float ***c, int nxt, int nyt, int nzt)
{
    int i, j, k;

    for (j = 1; j <= nyt; j++)
    {
        for (k = 1; k <= nzt; k++)
        {
            c[0][j][k] = c[1][j][k];
            c[nxt + 1][j][k] = c[nxt][j][k];
        }
    }

    for (i = 0; i <= nxt + 1; i++)
    {
        for (k = 1; k <= nzt; k++)
        {
            c[i][0][k] = c[i][1][k];
            c[i][nyt + 1][k] = c[i][nyt][k];
        }
    }

    for (i = 0; i <= nxt + 1; i++)
    {
        for (j = 0; j <= nyt + 1; j++)
        {
            c[i][j][0] = c[i][j][1];
            c[i][j][nzt + 1] = c[i][j][nzt];
        }
    }
}

float Bphi(float ***oc)
{
    int i, j, k;
    float value1, value2;

    value1 = value2 = 0.0;

    augmenc(oc, nx, ny, nz);

    ijkloop
    {
        value1 += (pow(oc[i + 1][j][k] - oc[i][j][k], 2) + pow(oc[i][j + 1][k] - oc[i][j][k], 2) + pow(oc[i][j][k + 1] - oc[i][j][k], 2)) / (h * h);
        value2 += pow(pow(oc[i][j][k], 2) - 1.0, 2);
    }

    return (0.5 * gam * value1 + 0.25 * value2 / gam) * h3 * 3.0 / (2.0 * sqrt(2.0));
}
