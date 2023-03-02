/*
 * Copyright (c) 2018-20, Sam Lambrick.
 * All rights reserved.
 * This file is part of the Sub-beam Ray Tracing simulation, subject to the
 * GNU/GPL-3.0-or-later.
 *
 * Intersect the ray with a particular surface. combinations of surfaces are
 * used to create a single interaction of the ray path.
 */
#include "intersect_detection3D.h"
#include "ray_tracing_core3D.h"
#include "distributions3D.h"
#include "mtwister.h"
#include <math.h>
#include <stdbool.h>

static void intersectPlane(Ray3D * const the_ray, Triangle element, double new_loc[3],
		bool * const hit, bool * const within);

/*
 * Finds the distance to, the normal to, and the position of a rays intersection
 * with an analytically defined sphere. Returns 0 if the ray does not intersect
 * the sphere. The inclusion of the analytic sphere is not very general and can
 * only really be used for the specific probable it was written for without
 * alteration--a single sphere places just touching a flat surface in the centre
 * of the scan region.
 *
 * INPUTS:
 *  the_ray    - pointer to a ray3D struct, contains information on the ray
 *  normal     - double array, a variable to store the normal to the point
 *               of intersection in
 *  t          - double, a variable to store the distance to the
 *               intersection point in
 *  the_sphere - AnalytSphere, contains all the important information about the
 *               analytic sphere
 *
 * OUTPUTS:
 *  intersect - int, 1 or 0 depending on if the ray intersects the sphere.
 *
 */
void scatterSphere(Ray3D * the_ray, AnalytSphere * the_sphere, double * const min_dist,
        double nearest_inter[3], double nearest_n[3], double nearest_b[6], int * const tri_hit,
        int * const which_surface, bool * const meets_sphere) {
    double a,b,c;
    double beta, gamma;
    double distance;
    double *e;
    double *d;

    /* Get the present position and direction of the ray */
    e = the_ray->position;
    d = the_ray->direction;

    /* Centre of the sphere */
    a = the_sphere->sphere_c[0];
    b = the_sphere->sphere_c[1];
    c = the_sphere->sphere_c[2];

    /* Coefficients of the quadratic equation */
    beta = 2*(d[0]*(e[0] - a) + d[1]*(e[1] - b) + d[2]*(e[2] - c));
    gamma = -the_sphere->sphere_r*the_sphere->sphere_r - 2*e[0]*a - 2*e[1]*b -
        2*e[2]*c + e[0]*e[0] + e[1]*e[1] + e[2]*e[2] + a*a + b*b + c*c;

    /* Do we hit the sphere */
    if (beta*beta - 4*gamma < 0) {
        *meets_sphere = false;
        return;
    }

    /* Solve the quadratic equation. Take the smaller root */
    distance = (-beta - sqrt(beta*beta - 4*gamma))/2;

    if ((distance*distance < *min_dist) && (distance > 0)) {
        /* Shortest intersection yet found */

        /* Normal to the sphere at that point */
        nearest_n[0] = (e[0] + distance*d[0] - a)/the_sphere->sphere_r;
        nearest_n[1] = (e[1] + distance*d[1] - b)/the_sphere->sphere_r;
        nearest_n[2] = (e[2] + distance*d[2] - c)/the_sphere->sphere_r;

        /* Intersection with the sphere */
        nearest_inter[0] = e[0] + d[0]*distance;
        nearest_inter[1] = e[1] + d[1]*distance;
        nearest_inter[2] = e[2] + d[2]*distance;

        *min_dist = distance*distance;

        normalise(nearest_n);

        /* We are not on a triangle */
        *tri_hit = -1;

        /* We are now on the sphere */
        *which_surface = the_sphere->surf_index;

        /* We do hit the sphere and it is the closest interesction*/
        *meets_sphere = true;
        return;
    }

    /* The sphere is not the closest intersection */
    *meets_sphere = false;
    return;
}

void scatterCircle(Ray3D * the_ray, Circle * the_circle, double * const min_dist,
        double nearest_inter[3], double nearest_n[3], double nearest_b[6], int * const tri_hit,
        int * const which_surface, bool * const meets_circle) {
    int i;
	double new_loc[3];
    double distance = 0;
    double * e;
    double * d;
    bool hit, within;
    Triangle element;

    /* Get the present position and direction of the ray */
    e = the_ray->position;
    d = the_ray->direction;

    get_elementCircle(the_circle, &element);
	intersectPlane(the_ray, element, new_loc, &hit, &within);

	if (!hit) {
		*meets_circle = false;
		return;
	}

    // TODO: generalise for a circle not in the xz plane
    double dist_x = fabs(the_circle->centre[0] - new_loc[0]);
    double dist_z = fabs(the_circle->centre[2] - new_loc[2]);
    if (dist_x*dist_x + dist_z*dist_z > the_circle->r*the_circle->r) {
        *meets_circle = false;
        return;
    }

	//
	for (i = 0; i < 3; i++)
		distance += (new_loc[i] - e[i])*(new_loc[i] - e[i]);

    if ((distance*distance < *min_dist) && (distance > 0)) {
        /* Shortest intersection yet found */

        nearest_n[0] = the_circle->normal[0];
        nearest_n[1] = the_circle->normal[1];
        nearest_n[2] = the_circle->normal[2];

        nearest_inter[0] = e[0] + d[0]*distance;
        nearest_inter[1] = e[1] + d[1]*distance;
        nearest_inter[2] = e[2] + d[2]*distance;

        *min_dist = distance*distance;

        /* We are not on a triangle */
        *tri_hit = -1;

        /* We are now on the sphere */
        *which_surface = the_circle->surf_index;

        /* We do hit the circle and it is the closest interesction*/
        *meets_circle = true;
        return;
    }

    /* The circle is not the closest intersection */
    *meets_circle = false;
    return;
}

/*
 * Find the point of intersection between a ray and a plane defined by three
 * points and a normal vector. Also determines if the intersection point is within
 * the triangle bounded by the three points.
 *
 * INPUTS:
 *  the_ray - Information on the current position of the ray
 *  a       - first of the 3 points
 *  b       - second of the 3 points
 *  c       - third of the 3 points
 *  normal  - normal vector
 *  new_loc - variable to store intersection point (if one exists)
 *  hit     - is the plane intersected
 *  within  - is the intersection within the triangle bounded by the 3 points
 */
static void intersectPlane(Ray3D * const the_ray, Triangle element, double new_loc[3],
		bool * const hit, bool * const within) {
	double AA[3][3];
	double v[3], u[3];
	double epsilon;
	double * e, * d;

	e = the_ray->position;
	d = the_ray->direction;

	/* If the triangle is 'back-facing' then the ray cannot hit it */
	double test;
	dot(element.normal, the_ray->direction, &test);
	if (test > 0) {
		*hit = false;
		return;
	}

    /*
     * Construct the linear equation
     * AA u = v, where u contains (alpha, beta, t) for the propagation
     * equation:
     * e + td = a + beta(b - a) + gamma(c - a)
     */
    //propagate3D(a, e, -1, v); // <- simpler to write, heavier computation
    v[0] = element.v1[0] - e[0];
    v[1] = element.v1[1] - e[1];
    v[2] = element.v1[2] - e[2];

    /* This could be pre-calculated and stored, however it would involve an
     * array of matrices
     */
    AA[0][0] = element.v1[0] - element.v2[0];
    AA[0][1] = element.v1[0] - element.v3[0];
    AA[0][2] = d[0];
    AA[1][0] = element.v1[1] - element.v2[1];
    AA[1][1] = element.v1[1] - element.v3[1];
    AA[1][2] = d[1];
    AA[2][0] = element.v1[2] - element.v2[2];
    AA[2][1] = element.v1[2] - element.v3[2];
    AA[2][2] = d[2];

    /*
     * Tests to see if this triangle is parallel to the ray, if it is the
     * determinant of the matrix AA will be zero, we must set a tolerance for
     * size of determinant we will allow.
     */
    epsilon = 0.0000000001;
    int success = 0; // Default to no success, this was causing problems some how...
    solve3x3(AA, u, v, epsilon, &success); // <- NOTE: this is the biggest computation
    if (!success) {
    	*hit = false;
    	return;
    }

    new_loc[0] = e[0] + (u[2]*d[0]);
    new_loc[1] = e[1] + (u[2]*d[1]);
    new_loc[2] = e[2] + (u[2]*d[2]);
    *hit = true;

    // Is the intersection point within the triangle defined by the 3 points
    *within = (u[0] >= 0) && (u[1] >= 0) && ((u[0] + u[1]) <= 1) && (u[2] > 0);

}

void scatterPlane(Ray3D * the_ray, Plane plane, double * const min_dist,
		double nearest_inter[3], double nearest_n[3], double nearest_b[6], int * const meets,
		int * const which_surface) {
	double *e;
	bool hit, within;
	double new_loc[3];
	Triangle element;

	e = the_ray->position;
	*meets = 0;

	get_elementPlane(&plane, &element);
	intersectPlane(the_ray, element, new_loc, &hit, &within);

	if (hit) {
		int dist;
		double movment[3];

		/* Movement is the vector from the current location to the possible
		 * new location */
		movment[0] = new_loc[0] - e[0];
		movment[1] = new_loc[1] - e[1];
		movment[2] = new_loc[2] - e[2];

		/* NOTE: we are comparing the square of the distance */
		dist = movment[0]*movment[0] + movment[1]*movment[1] +
				movment[2]*movment[2];
        if (dist < *min_dist) {
        	*meets = 1;
            /* This is the smallest intersection found so far */
            *min_dist = dist;

            nearest_n[0] = element.normal[0];
            nearest_n[1] = element.normal[1];
            nearest_n[2] = element.normal[2];
            nearest_inter[0] = new_loc[0];
            nearest_inter[1] = new_loc[1];
            nearest_inter[2] = new_loc[2];

            *which_surface = plane.surf_index;
        }
	}
}

/*
 * Finds the distance to, the normal to, and the position of a ray's intersection
 * with an triangulated surface.
 *
 * INPUTS:
 *
 *  current_tri     - int, the index of the triangle the current ray is on, -1
 *                    indicates that the ray is not on any triangle
 *  current_surface - int, an index stating which surface the ray is on,
 *                    0=sample, 1=pinhole plate, -1=none, -2=sphere
 *  min_dist        - double pointer, to store the minimum distance to a surface
 *                    in
 *  nearest_inter   - double array, to store the location of the nearest
 *                    intersection in
 *  nearest_n       - double array, to store the normal to the surface at the
 *                    nearest intersection
 *  meets           - int, 1 or 0 to store if the ray hits any surfaces
 *  tri_hit         - int pointer, to store the index of the triangle that the
 *                    ray hits (if it hits any)
 *  which_surface   - int pointer, which surface does the ray intersect (if any)
 *
 * NOTE: this function is messy as attempts (mostly successful) have been made
 *       to improve the speed of the simulation as this is the section of code
 *       called the highest number of times, hence the rather low level looking
 *       code.
 */
void scatterTriag(Ray3D * the_ray, Surface3D * sample, double * const min_dist,
        double nearest_inter[3], double nearest_n[3], double nearest_b[6], int * const meets, int * const tri_hit,
        int * const which_surface) {
    int j;
    double *e;

    /* Position and direction of the ray */
    e = the_ray->position;

    /* Loop through all triangles in the surface */
    for (j = 0; j < sample->n_faces; j++) {
    	double new_loc[3];
        bool hit, within;
        Triangle element;

        /* Skip this triangle if the ray is already on it */
        if ((the_ray->on_element == j) && (the_ray->on_surface == sample->surf_index)) {
            continue;
        }

        /*
         * Specify which triangle and get its normal.
         */
        get_element3D(sample, j, &element);
        intersectPlane(the_ray, element, new_loc, &hit, &within);

        if (!hit)
        	continue;

        /* Find if the point of intersection is inside the triangle */
        /* Must also find if the ray is propagating forwards */
        if (within) {
            double movment[3];
            double dist;

            /* We have hit a triangle */
            *meets = 1; // <- I think I've found the problem....

            /* Movement is the vector from the current location to the possible
             * new location */
            movment[0] = new_loc[0] - e[0];
            movment[1] = new_loc[1] - e[1];
            movment[2] = new_loc[2] - e[2];

            /* NOTE: we are comparing the square of the distance */
            dist = movment[0]*movment[0] + movment[1]*movment[1] +
                movment[2]*movment[2];

            if (dist < *min_dist) {
                int k;
                /* This is the smallest intersection found so far */
                *min_dist = dist;

                *tri_hit = j;
                for (k = 0; k < 3; k++) {
                    nearest_n[k] = element.normal[k];
                    nearest_b[k] = element.lattice[k];
                    nearest_inter[k] = new_loc[k];
                }
                for (k = 3; k < 6; k++)
                    nearest_b[k] = element.lattice[k];

                *which_surface = sample->surf_index;
            }
        }
    }
}

/*
 * Attempt to scatter off the total sample specification, that is a triangulated surface
 * an analytical sphere, and a circular surface.
 */
void scatterSample(Ray3D * the_ray, Sample overall_sample, double * min_dist, double * nearest_inter,
                   double * nearest_n, double * nearest_b, int * meets, int * tri_hit, int * which_surface) {
    bool meets_sphere = 0;
    bool meets_circle = 0;

    scatterTriag(the_ray, overall_sample.triag_sample, min_dist, nearest_inter, nearest_n,
        nearest_b, meets, tri_hit, which_surface);

    for (int i = 0; i < overall_sample.n_sphere; i ++) {
        bool meets_this_sphere = 0;
        if (overall_sample.the_sphere[i].make_sphere) {
            if (the_ray->on_surface != overall_sample.the_sphere[i].surf_index) {
                scatterSphere(the_ray, &overall_sample.the_sphere[i], min_dist, nearest_inter,
                    nearest_n, nearest_b, tri_hit, which_surface, &meets_this_sphere);
                meets_sphere = meets_sphere || meets_this_sphere;
            }
        }
    }

    if (overall_sample.the_circle->make_circle) {
        if (the_ray->on_surface != overall_sample.the_circle->surf_index) {
            scatterCircle(the_ray, overall_sample.the_circle, min_dist, nearest_inter,
            	nearest_n, nearest_b, tri_hit, which_surface, &meets_circle);
        }
    }

    *meets = *meets || meets_circle || meets_sphere;
}

/*
 * Scatter off a back wall with n detectors. Returns 0 if the ray is not
 * detected or the index of the detector is has entered.
 */
void multiBackWall(Ray3D * the_ray, NBackWall wallPlate, double * const min_dist,
        double nearest_inter[3], double nearest_n[3], double nearest_b[6], int * meets,
        int * const tri_hit, int * const which_surface, int * const which_aperture) {
    double *d;
    d = the_ray->direction;

    /*
     * If the ray is travelling in the positive y-direction then it may hit
     * the 'back wall'.
     */
    *meets = 0;
    if (d[1] > 0) {
        double alpha, x;
        double wall_hit[3];
        double backNormal[3];
        double test;
        double r, x_disp, y_disp;
        double *e;
        *which_aperture = 0;
        int i;

        e = the_ray->position;

        backNormal[0] = 0;
        backNormal[1] = -1;
        backNormal[2] = 0;

        /*
         * Find where the ray hits the back wall, the back wall is defined to be
         * in the plane y = 0
         */
        alpha = (0 - e[1])/d[1];

        /*
         * If the distance to the back wall is longer than a previous
         * intersection then the ray does not hit the back wall.
         */
        if (alpha*alpha > *min_dist) {
            *which_aperture = 0;
            return;
        }

        /* propagate the ray to that position */
        propagate(e, d, alpha, wall_hit);

        test = 10; /* Number greater than 1 */
        for (i = 0; i < wallPlate.n_detect; i++) {
            BackWall plate;
            double tmp;
            get_nth_aperture(i, &wallPlate, &plate);
            /*
             * See if it goes into the aperture. The aperture is defined by a centre,
             * along with two axes. Uses the formula for an ellipse:
             *     (x-h)^2/a^2 + (z-k)^2/b^2 = 1
             * where h and k are the x and z positions of the centre of the ellipse.
             */
			// TODO: now we have a rotation of the elliptical aperture...
			// how to do that, some maths is in order

			// Move to centre
            x_disp = wall_hit[0] - plate.aperture_c[0];
            y_disp = wall_hit[2] - plate.aperture_c[1];

			// Rotate the point in the opposite direction to the direction of
			// rotation of the ellipse
			tmp = x_disp*cos(- plate.aperture_rotate*M_PI/180) - y_disp*sin(- plate.aperture_rotate*M_PI/180);
			y_disp = x_disp*sin(- plate.aperture_rotate*M_PI/180) + y_disp*cos(- plate.aperture_rotate*M_PI/180);
			x_disp = tmp;

            test = x_disp*x_disp/
                (0.25*plate.aperture_axes[0]*plate.aperture_axes[0]) +
                y_disp*y_disp/
                (0.25*plate.aperture_axes[1]*plate.aperture_axes[1]);

            if (test < 1) {
                *which_aperture = i + 1;
                break;
            }
        }

        if (test < 1) {
            /* Goes into detector aperture */
            /* Check that the distance to the aperture is the smallest so far */
            if (alpha*alpha < *min_dist) {
                *min_dist = alpha*alpha;

                /* -1 is not being on a triangle */
                *tri_hit = -1;
                nearest_n[0] = backNormal[0];
                nearest_n[1] = backNormal[1];
                nearest_n[2] = backNormal[2];
                nearest_inter[0] = wall_hit[0];
                nearest_inter[1] = wall_hit[1];
                nearest_inter[2] = wall_hit[2];

                *which_surface = wallPlate.surf_index;

                /*
                 * If we have gone into the aperture then the ray is both dead and
                 * should be counted.
                 */
                return;
            }
        }

        /* If the position is not then we can scatter of the back wall if that is
         * what is wanted */
        x = wall_hit[0]*wall_hit[0] + wall_hit[2]*wall_hit[2];
        r = wallPlate.circle_plate_r;
        if ((x <=  r*r) && wallPlate.plate_represent) {
            /* We have met a surface */
            *meets = 1;

            if (alpha*alpha < *min_dist) {
                *min_dist = alpha*alpha;

                /* -1 is not being on a triangle */
                *tri_hit = -1;
                nearest_n[0] = backNormal[0];
                nearest_n[1] = backNormal[1];
                nearest_n[2] = backNormal[2];
                nearest_inter[0] = wall_hit[0];
                nearest_inter[1] = wall_hit[1];
                nearest_inter[2] = wall_hit[2];

                *which_surface = wallPlate.surf_index;
            }
        }
    }

    /* The ray has not been detected. */
    *which_aperture = 0;
    return;
}

/* Scatter off an abstract hemisphere */
void abstractScatter(Ray3D * the_ray, AbstractHemi const * plate, double * const min_dist,
        int * const tri_hit, int * const which_surface, int * const which_aperture) {
    double angle_sep, tmp;

    // Get the angle between the ray and the
    dot(the_ray->direction, plate->det_dir, &tmp);
    angle_sep = acos(tmp);

    if (angle_sep > plate->half_cone_angle) {
        return;
    } else if (*min_dist > 1e6) { // 1e6 = 1km
        *which_aperture = 1;
        *which_surface = plate->surf_index;
        *tri_hit = -1;
    }
}


/*
 * Scatters the ray off a corrugated surface.
 *
 * TODO: make this
 */
