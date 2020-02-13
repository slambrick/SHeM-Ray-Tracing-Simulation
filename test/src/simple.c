/*
 * Copyright (c) 2020, Dan Seremet.
 * All rights reserved.
 * This file is part of the SHeM ray tracing simulation, subject to the
 * GNU/GPL-3.0-or-later.
 *
 * Test for the simple 3D functions, such as vector normalisation,
 * reflection and finding a parametrization of a perpendicular plane.
 */

#include <stdio.h>
#include <math.h>
#include <assert.h>

#include "small_functions3D.h"

const double epsilon = 1e-9;


int main() {
    double normal[] = {0, 0, 1};
    double direction[] = {1, 0, -1};
    double new_dir[3];

    double t1[3], t2[3];

    // test normalisation
    normalise(direction);
    assert(fabs(direction[0] - 1/sqrt(2)) < epsilon);
    assert(fabs(direction[1]) < epsilon);
    assert(fabs(direction[2] + 1/sqrt(2)) < epsilon);

    // test reflection
    reflect3D(normal, direction, new_dir);
    assert(fabs(new_dir[0] - 1/sqrt(2)) < epsilon);
    assert(fabs(new_dir[1]) < epsilon);
    assert(fabs(new_dir[2] - 1/sqrt(2)) < epsilon);

    // test that perpendicular plane produces perpendicular planes
    perpendicular_plane(new_dir, t1, t2);
    assert(fabs(dot(new_dir, t1)) < epsilon);
    assert(fabs(dot(new_dir, t2)) < epsilon);

    // test that these directions are normalised
    assert(fabs(norm2(t1) - 1) < epsilon);
    assert(fabs(norm2(t2) - 1) < epsilon);

    printf("\nAll tests passed.\n");
    return 0;
}
