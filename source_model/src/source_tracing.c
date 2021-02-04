/*
 * source_tracing.c
 *
 * Part of the SHeM ray tracing simulation.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <sys/time.h>
#include "mtwister.h"
#include "atom_ray_tracing3D.h"
#include "math.h"

typedef struct _param {
	double skimmer_d;
	double pinhole_d;
	double skimmer_pinhole_dist;

} SimParam;

void write_beam(Rays3D const * const beam, char const * const fname, double ds, double dp, double x);

void gen_pos_circle(double r0, double centre[3], double normal[3], double pos[3],
		MTRand * myrng);

void get_hit_rays(int const * const hits, Rays3D const * const all_rays, Rays3D * const beam);

int main(int argc, char * argv []) {
	// Input arguments
	int n_rays;
	double x1;
	if (argc == 1)
		printf("Defaulting to 100 rays and 21.5cm source distance.\n");
	if (argc == 2)
		printf("Defaulting to 21.5cm source distance\n");
	else if (argc != 3)
		printf("Unexpected number of arguments, ignoring after the first\n");

	n_rays = argc > 1 ? atoi(argv[1]) : 100;
	x1 = argc > 2 ? atof(argv[2]) : 21.5e3;

	// Other variables
	int i;
	unsigned long t;
	struct timeval tv;
	MTRand myrng;

	// Overall parameters
	// TODO: move these to command line arguments
	double d_s = 100e-3;
	double d_p = 0.5e-3;
	double skim_centre[3] = {0, 0, 0};

	double theta_max = atan(d_s/x1);
	Circle pinhole;

	// Seed random number generator
    gettimeofday(&tv, 0);
    t = (unsigned long)tv.tv_sec + (unsigned long)tv.tv_usec;
	seedRand(t, &myrng);

	// Parameters for generating rays
	double new_dir[3];
	double normal[3] = {0, 0, 1};

	// Create a Circle struct for the pinhole
	double pinhole_c[3] = {0, 0, x1};
	double pinhole_normal[3] = {1/sqrt(2), 0, -1/sqrt(2)};
	Material M;
    M.name = "diffuse";
    M.func_name = "cosine";
    M.params = 0;
    M.n_params = 0;
    M.func = distribution_by_name("cosine");
	set_up_circle(1, pinhole_c, d_p/2, pinhole_normal, M, 0, &pinhole);

	// Which rays go through the pinhole
	int * hits;
	hits = calloc(n_rays, sizeof(int));

	// Generate rays
	Rays3D all_rays;
	all_rays.nrays = n_rays;
	all_rays.rays = calloc(all_rays.nrays, sizeof(Ray3D));
	for (i = 0; i < all_rays.nrays; i++) {
		double init_dir[3];
		double pos[3];
		Ray3D the_ray;

		uniform_scatter(normal, init_dir, new_dir, &theta_max, &myrng);
		gen_pos_circle(d_s/2, skim_centre, normal, pos, &myrng);
		new_Ray(&the_ray, pos, new_dir);
		all_rays.rays[i] = the_ray;
	}

	printf("Constructed a beam at the skimmer\n");

	for (i = 0; i < all_rays.nrays; i++) {
		bool meets_pinhole;
		int which_surface, tri_hit;
		double pinhole_inter[3];
		double pinhole_n[3];
		double min_dist = 1e20;

		// Test if the rays pass through the pinhole
		meets_pinhole = false; // By default do not meet the pinhole
		which_surface = -1;    // Start out not on any surface
		tri_hit = -1;          // Always not on a triangle because there aren't any

		scatterCircle(&all_rays.rays[i], pinhole, &min_dist, pinhole_inter, pinhole_n,
				&tri_hit, &which_surface, &meets_pinhole);
		hits[i] = meets_pinhole ? 1 : 0;
		if (meets_pinhole) {
			update_ray_position(&all_rays.rays[i], pinhole_inter);
		}
	}

	printf("Finished tracing the beam to the pinhole, saving beam state\n");

	// Resulting beam
	Rays3D beam;
	get_hit_rays(hits, &all_rays, &beam);
	free(hits);
	clean_up_rays(all_rays);

	printf("Here");
	// Write out the beam state to file
	char fname[50] = "results/beam.csv";
	write_beam(&beam, fname, d_s, d_p, x1);

	clean_up_rays(beam);
}

void write_beam(Rays3D const * const beam, char const * const fname, double ds, double dp, double x) {
	FILE * fid;
	int i;

	fid = fopen(fname, "w");
	fprintf(fid, "Source ray simulation with 45deg pinhole.\n");
	fprintf(fid, "Skimmer diameter = %f mm.\n", ds);
	fprintf(fid, "Pinhole diameter = %f mm.\n", dp);
	fprintf(fid, "Skimmer - pinhole distance = %f mm.\n", x);
	fprintf(fid, "Resulting beam in the pinhole:\n");
	fprintf(fid, "x,y,z,vx,vy,vz\n");
	for (i = 0; i < beam->nrays; i++) {
		double * e = beam->rays[i].position;
	    double * d = beam->rays[i].direction;
		fprintf(fid, "%f,%f,%f,%f,%f,%f\n", e[0], e[1], e[2], d[0], d[1], d[2]);
	}
	fclose(fid);
}

void get_hit_rays(int const * const hits, Rays3D const * const all_rays, Rays3D * const beam) {
	int i, j, n_beam;

	// Number of rays in the beam at the pinhole
	sum_array_int(hits, all_rays->nrays, &n_beam);
	printf("We collected %i rays in the pinhole.\n", n_beam);
	beam->rays = calloc(n_beam, sizeof(Ray3D));
	beam->nrays = n_beam;
	j = 0;
	for (i = 0; i < all_rays->nrays; i++) {
		if (hits[i]) {
			beam->rays[j] = all_rays->rays[i];
			j++;
		}
	}
}

void gen_pos_circle(double r0, double centre[3], double normal[3], double pos[3],
		MTRand * myrng) {
	double r, phi, uni_rand, c_theta;
	int i;
	double u[3];
	double n0[3] = {0, 0, 1};
	double pos2[3];

	// Create random position in a circle
	genRand(myrng, &uni_rand);
	genRand(myrng, &phi);
	phi = phi*2*M_PI;
	r = cbrt(r0*r0*r0 * uni_rand);

	pos[0] = r*cos(phi);
	pos[1] = r*sin(phi);
	pos[2] = 0;

	// Rotate to match the normal of the circle
	cross(normal, n0, u);
	dot(normal, n0, &c_theta);
	general_rotation(pos, u, pos2, c_theta);

	// Translate to the centre of the circle
	for (i = 0; i < 3; i++)
		pos[i] = pos2[i] + centre[i];
}

