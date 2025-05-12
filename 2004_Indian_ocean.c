#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <string.h>
#include <sys/stat.h>
#include <assert.h>

#ifdef _WIN32
#include <direct.h>
#define mkdir(dir, mode) _mkdir(dir)
#endif

// ======================
// REALISTIC PARAMETERS
// ======================
#define NX 200       // 200 points (100 km at Dx=500m)
#define NY 200       // 200 points (100 km at Dx=500m)
#define NT 6000      // 3000 sec at Dt=0.5s

// Physical parameters - 2004 Indian Ocean Event
const double g = 9.81;       // Gravity [m/s²] (universal)
const double ho = 4000.0;    // Initial depth at Sunda Trench [m] (Satake et al., 2006)
const double fbo = 0.002;    // Bed friction (silty seabed near Sumatra) (Jaffe et al., 2006)
const double zbmax = 10.0;   // Max seabed uplift [m] (Lay et al., 2005 - peak displacement)
const double nd = 50.0;      // Rupture duration [s] (8-10 min total rupture time)

// Numerical parameters
const double Dx = 1000.0;    // Spatial step [m] (1km resolution for regional scale)
const double Dt = 1.0;       // Time step [s] (CFL-stable for Dx=1000m, ho=4000m)
// ======================
// HELPER FUNCTIONS
// ======================
void prepare_output_dir(const char* dirname) {
    struct stat st = {0};
    if (stat(dirname, &st) == -1) {
        #ifdef _WIN32
        _mkdir(dirname);
        #else
        mkdir(dirname, 0700);
        #endif
    }
    char command[256];
    #ifdef _WIN32
    snprintf(command, sizeof(command), "del /Q \"%s\\*\"", dirname);
    #else
    snprintf(command, sizeof(command), "rm -f %s/*", dirname);
    #endif
    system(command);
}

void write_csv(const char* filename, double* data, int nx, int ny) {
    FILE* fp = fopen(filename, "w");
    if (!fp) {
        perror("Error opening file");
        exit(1);
    }
    for (int i = 0; i < nx; i++) {
        for (int j = 0; j < ny; j++) {
            fprintf(fp, "%.4f", data[i*ny + j]);
            if (j < ny - 1) fprintf(fp, ",");
        }
        fprintf(fp, "\n");
    }
    fclose(fp);
}

// ======================
// MAIN SIMULATION
// ======================
int main() {
    // Derived parameters
    const double r = g * pow(Dt/Dx, 2) / 2.0;
    const double c = sqrt(g * ho);       //wave speed ~200m/s
    const double alpha = c * Dt / Dx;    // Courant number
    assert(alpha < 1.0 && "CFL condition violated!");

    const int cx = NX/2;                 // Epicenter X
    const int cy = NY/2;                 // Epicenter Y
    const int radius = 30;               // ~30 km radius

    prepare_output_dir("wave_data_india");

    // Allocate arrays
    double* h = malloc(NX*NY*sizeof(double));
    double* fb = malloc(NX*NY*sizeof(double));
    double* z = malloc(NX*NY*sizeof(double));
    double* zo = malloc(NX*NY*sizeof(double));
    double* zn = malloc(NX*NY*sizeof(double));
    double* bed = malloc(NX*NY*sizeof(double));

    // Initialize
    for (int i = 0; i < NX*NY; i++) {
        h[i] = ho;
        fb[i] = fbo;
        z[i] = (i == cx*NY + cy) ? 0.1 * zbmax : 0.0;
        zo[i] = 0.0;
        bed[i] = 0.0;
    }

    double max_crest = -INFINITY;
    double max_trough = INFINITY;
    double zb = 0.0, zbo = 0.0;

    // Time loop
    for (int k = 0; k < NT; k++) {
        double t = (k+1)*Dt;

        // Seabed motion (linear ramp)
        double zbn = (t <= nd) ? zbmax*t/nd : zbmax;

        // Apply seabed uplift (Gaussian distribution)
        for (int i = cx-radius; i <= cx+radius; i++) {
            for (int j = cy-radius; j <= cy+radius; j++) {
                if (i >= 0 && i < NX && j >= 0 && j < NY) {
                    double dist = sqrt(pow(i-cx,2) + pow(j-cy,2));
                    if (dist <= radius) {
                        bed[i*NY+j] = (zbn - 2*zb + zbo) * exp(-pow(dist/radius,2));
                    }
                }
            }
        }

        // Wave update (2D wave equation)
        for (int i = 1; i < NX-1; i++) {
            for (int j = 1; j < NY-1; j++) {
                double h_avg_x = (h[(i+1)*NY+j] + 2*h[i*NY+j] + h[(i-1)*NY+j])/4.0;
                double h_avg_y = (h[i*NY+j+1] + 2*h[i*NY+j] + h[i*NY+j-1])/4.0;

                double dzdx = (z[(i+1)*NY+j] - z[(i-1)*NY+j])/(2.0*Dx);
                double dzdy = (z[i*NY+j+1] - z[i*NY+j-1])/(2.0*Dx);

                double d2zdx2 = (z[(i+1)*NY+j] - 2*z[i*NY+j] + z[(i-1)*NY+j])/pow(Dx,2);
                double d2zdy2 = (z[i*NY+j+1] - 2*z[i*NY+j] + z[i*NY+j-1])/pow(Dx,2);

                zn[i*NY+j] = 2*z[i*NY+j] - zo[i*NY+j]
                           + g*h[i*NY+j]*Dt*Dt*(d2zdx2 + d2zdy2)
                           - fbo*Dt*(z[i*NY+j] - zo[i*NY+j])
                           + bed[i*NY+j];

                // Track extremes
                max_crest = fmax(max_crest, zn[i*NY+j]);
                max_trough = fmin(max_trough, zn[i*NY+j]);
            }
        }

        // Sommerfeld boundary conditions
        for (int j = 0; j < NY; j++) {
            zn[0*NY+j] = zo[0*NY+j] + alpha * (zn[1*NY+j] - zo[0*NY+j]); // Left
            zn[(NX-1)*NY+j] = zo[(NX-1)*NY+j] + alpha * (zn[(NX-2)*NY+j] - zo[(NX-1)*NY+j]); // Right
        }
        for (int i = 0; i < NX; i++) {
            zn[i*NY+0] = zo[i*NY+0] + alpha * (zn[i*NY+1] - zo[i*NY+0]); // Bottom
            zn[i*NY+NY-1] = zo[i*NY+NY-1] + alpha * (zn[i*NY+NY-2] - zo[i*NY+NY-1]); // Top
        }

        // Update states
        memcpy(zo, z, NX*NY*sizeof(double));
        memcpy(z, zn, NX*NY*sizeof(double));
        zbo = zb;
        zb = zbn;



        char filename[128];
        snprintf(filename, sizeof(filename), "wave_data_india/frame_%05d.csv", k);
        write_csv(filename, z, NX, NY);

    }

    printf("Simulation complete. Domain: %dx%d km\n", (int)(NX*Dx/1000), (int)(NY*Dx/1000));
    printf("Maximum wave crest: %.2f m \n", max_crest);
    printf("Maximum wave trough: %.2f m\n", max_trough);

    free(h); free(fb); free(z); free(zo); free(zn); free(bed);
    return 0;
}
